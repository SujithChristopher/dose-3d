# %% [markdown]
"""
# Quick EPID Dose Reconstruction Verification

Fast EPID dose reconstruction for verification purposes.
"""

# %% 
import math
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from pydicom import dcmread
from tqdm import tqdm
from scipy.fft import fft, ifft
from concurrent.futures import ThreadPoolExecutor
import warnings
warnings.filterwarnings('ignore')

PI = math.pi

# Try numba for acceleration
try:
    from numba import jit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def jit(func):
        return func
    prange = range

# %% [markdown]
"""
## Fast Data Loading and Preprocessing
"""

# %%
def load_epid_subset(data_path, max_files=None, every_nth=1, max_workers=8):
    dicom_files = [f for f in os.listdir(data_path) if f.endswith('.dcm')]
    dicom_files.sort()
    
    if every_nth > 1:
        dicom_files = dicom_files[::every_nth]
    if max_files:
        dicom_files = dicom_files[:max_files]
    
    file_paths = [os.path.join(data_path, f) for f in dicom_files]
    
    def read_epid_fast(file_path):
        try:
            dcm = dcmread(file_path)
            dose_data = dcm.pixel_array.astype(np.float32)
            angle = float(getattr(dcm, 'GantryAngle', 0.0))
            return dose_data, angle
        except:
            return None, None
    
    dose_images, angles = [], []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        with tqdm(total=len(file_paths), desc="Loading") as pbar:
            futures = {executor.submit(read_epid_fast, path): path for path in file_paths}
            for future in futures:
                dose_data, angle = future.result()
                if dose_data is not None:
                    dose_images.append(dose_data)
                    angles.append(angle)
                pbar.update(1)
    
    return np.array(dose_images), np.array(angles)

def extract_dose_increments(dose_images, angles):
    n_images = len(dose_images)
    dose_increments = np.zeros_like(dose_images)
    
    dose_increments[0] = dose_images[0].copy()
    for i in range(1, n_images):
        increment = dose_images[i] - dose_images[i-1]
        dose_increments[i] = np.maximum(increment, 0)
    
    sort_indices = np.argsort(angles)
    return dose_increments[sort_indices], angles[sort_indices]

# %% [markdown]
"""
## Optimized reconstruction functions with numba acceleration
"""

# %%
if NUMBA_AVAILABLE:
    @jit(nopython=True)
    def fast_dose_weighting(dose_increment, sid, pixel_size):
        nrows, ncols = dose_increment.shape
        col_coords = pixel_size * (np.arange(ncols) - ncols/2 + 0.5)
        row_coords = pixel_size * (np.arange(nrows) - nrows/2 + 0.5)
        weighted_dose = np.zeros_like(dose_increment)
        
        for i in prange(nrows):
            for j in prange(ncols):
                r_sq = row_coords[i]**2 + col_coords[j]**2
                distance = np.sqrt(sid**2 + r_sq)
                weight = (sid / distance)**2
                weighted_dose[i, j] = dose_increment[i, j] * weight
        return weighted_dose
    
    @jit(nopython=True)
    def fast_backproject_dose(weighted_dose, angle_rad, sad, sid, pixel_size, vol_size, vol_extent):
        nrows, ncols = weighted_dose.shape
        dose_volume = np.zeros((vol_size, vol_size, vol_size))
        coords = np.linspace(-vol_extent, vol_extent, vol_size)
        cos_angle, sin_angle = np.cos(angle_rad), np.sin(angle_rad)
        
        for i in prange(vol_size):
            for j in prange(vol_size):
                for k in prange(vol_size):
                    x, y, z = coords[i], coords[j], coords[k]
                    source_to_voxel = sad + x * sin_angle - y * cos_angle
                    
                    if source_to_voxel > 0.1:
                        det_x = (x * cos_angle + y * sin_angle) * (sid / source_to_voxel)
                        det_y = z * (sid / source_to_voxel)
                        px = det_x / pixel_size + ncols / 2
                        py = det_y / pixel_size + nrows / 2
                        
                        if 0 <= px < ncols-1 and 0 <= py < nrows-1:
                            px0, py0 = int(px), int(py)
                            wx, wy = px - px0, py - py0
                            dose_val = ((1-wx)*(1-wy)*weighted_dose[py0, px0] +
                                       wx*(1-wy)*weighted_dose[py0, px0+1] +
                                       (1-wx)*wy*weighted_dose[py0+1, px0] +
                                       wx*wy*weighted_dose[py0+1, px0+1])
                            distance_weight = (sad / source_to_voxel)**2
                            dose_volume[i, j, k] += dose_val * distance_weight
        return dose_volume
else:
    def fast_dose_weighting(dose_increment, sid, pixel_size):
        nrows, ncols = dose_increment.shape
        col_coords = pixel_size * (np.arange(ncols) - ncols/2 + 0.5)
        row_coords = pixel_size * (np.arange(nrows) - nrows/2 + 0.5)
        row_grid, col_grid = np.meshgrid(row_coords, col_coords, indexing='ij')
        distance = np.sqrt(sid**2 + row_grid**2 + col_grid**2)
        return dose_increment * (sid / distance)**2
    
    def fast_backproject_dose(weighted_dose, angle_rad, sad, sid, pixel_size, vol_size, vol_extent):
        return np.zeros((vol_size, vol_size, vol_size))

print(f"Reconstruction functions loaded (Numba: {NUMBA_AVAILABLE})")

# %%
class QuickEPIDReconstructor:
    def __init__(self, sad=1000.0, sid=1500.0, pixel_size=0.172, vol_size=25):
        self.sad = sad
        self.sid = sid
        self.pixel_size = pixel_size
        self.vol_size = vol_size
        self.vol_extent = 50.0
    
    def reconstruct(self, dose_increments, angles, chunk_size=50):
        n_projections = len(dose_increments)
        angles_rad = angles * np.pi / 180.0
        dose_volume = np.zeros((self.vol_size, self.vol_size, self.vol_size), dtype=np.float32)
        
        chunks = [list(range(i, min(i + chunk_size, n_projections))) 
                 for i in range(0, n_projections, chunk_size)]
        
        for chunk_indices in tqdm(chunks, desc="Reconstructing"):
            chunk_contribution = np.zeros_like(dose_volume)
            
            for proj_idx in chunk_indices:
                dose_increment = dose_increments[proj_idx]
                if np.max(dose_increment) <= 0:
                    continue
                
                weighted_dose = fast_dose_weighting(dose_increment, self.sid, self.pixel_size)
                
                if NUMBA_AVAILABLE:
                    dose_contrib = fast_backproject_dose(
                        weighted_dose, angles_rad[proj_idx], self.sad, self.sid, 
                        self.pixel_size, self.vol_size, self.vol_extent
                    )
                    chunk_contribution += dose_contrib
            
            dose_volume += chunk_contribution
        
        return dose_volume

# %% [markdown] 
"""
## Configuration and Data Loading
"""

# %%
# Configuration
QUICK_TEST = True
MAX_FILES = 100 if QUICK_TEST else None
EVERY_NTH = 4 if QUICK_TEST else 1
VOL_SIZE = 25 if QUICK_TEST else 50

data_path = r"dataset/VMAT 2025 - 6. SIB COMPLEX TARGET/T1/873251691"

if os.path.exists(data_path):
    # Load data
    dose_images, angles = load_epid_subset(data_path, MAX_FILES, EVERY_NTH)
    dose_increments, sorted_angles = extract_dose_increments(dose_images, angles)
    
    # Get geometry parameters
    first_file = os.path.join(data_path, [f for f in os.listdir(data_path) if f.endswith('.dcm')][0])
    sample_dcm = dcmread(first_file)
    sad = float(getattr(sample_dcm, 'RadiationMachineSAD', 1000.0))
    sid = float(getattr(sample_dcm, 'RTImageSID', 1500.0))
    pixel_size = 0.172
    
    # Reconstruct
    reconstructor = QuickEPIDReconstructor(sad, sid, pixel_size, VOL_SIZE)
    
    start_time = time.time()
    dose_volume = reconstructor.reconstruct(dose_increments, sorted_angles)
    recon_time = time.time() - start_time
    
    print(f"Reconstruction: {recon_time:.2f}s, {len(dose_increments)/recon_time:.1f} proj/s")
else:
    print("Data path not found")

# %% [markdown]
"""
## Visualization
"""

# %%
# Visualization
center = VOL_SIZE // 2
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Cross sections
im1 = axes[0,0].imshow(dose_volume[:, :, center].T, cmap='hot', origin='lower')
axes[0,0].set_title('Axial')
plt.colorbar(im1, ax=axes[0,0])

im2 = axes[0,1].imshow(dose_volume[center, :, :].T, cmap='hot', origin='lower')
axes[0,1].set_title('Sagittal')
plt.colorbar(im2, ax=axes[0,1])

# Profiles
coords = np.arange(VOL_SIZE) - center
axes[1,0].plot(coords, dose_volume[:, center, center], 'b-', label='X')
axes[1,0].plot(coords, dose_volume[center, :, center], 'r-', label='Y')
axes[1,0].plot(coords, dose_volume[center, center, :], 'g-', label='Z')
axes[1,0].set_title('Profiles')
axes[1,0].legend()
axes[1,0].grid(True)

# Statistics
axes[1,1].axis('off')
stats = f"""Volume: {VOL_SIZE}³
Projections: {len(dose_increments)}
Time: {recon_time:.2f}s
Rate: {len(dose_increments)/recon_time:.1f} proj/s

Dose Range:
Min: {dose_volume.min():.2f}
Max: {dose_volume.max():.2f}
Mean: {dose_volume.mean():.2f}

Coverage: {100*np.count_nonzero(dose_volume)/dose_volume.size:.1f}%
Numba: {'Yes' if NUMBA_AVAILABLE else 'No'}"""

axes[1,1].text(0.1, 0.9, stats, transform=axes[1,1].transAxes, 
               fontsize=10, verticalalignment='top', fontfamily='monospace')

plt.tight_layout()
plt.show()
