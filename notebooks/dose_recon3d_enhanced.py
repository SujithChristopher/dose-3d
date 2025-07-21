# %% [markdown]
"""
# Enhanced EPID Dose Reconstruction with Detector Translation Correction

Enhanced EPID dose reconstruction that accounts for X-Ray Image Receptor Translation 
and detector positioning corrections for improved geometric accuracy.
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
## Enhanced Data Loading with Detector Translation Extraction
"""

# %%
def load_epid_with_geometry(data_path, max_files=None, every_nth=1, max_workers=8):
    """Load EPID data with enhanced geometry information including detector translations."""
    dicom_files = [f for f in os.listdir(data_path) if f.endswith('.dcm')]
    dicom_files.sort()
    
    if every_nth > 1:
        dicom_files = dicom_files[::every_nth]
    if max_files:
        dicom_files = dicom_files[:max_files]
    
    file_paths = [os.path.join(data_path, f) for f in dicom_files]
    
    def read_epid_enhanced(file_path):
        try:
            dcm = dcmread(file_path)
            dose_data = dcm.pixel_array.astype(np.float32)
            angle = float(getattr(dcm, 'GantryAngle', 0.0))
            
            # Extract detector translation if available
            detector_translation = None
            if hasattr(dcm, 'XRayImageReceptorTranslation'):
                # DICOM tag (3002,000D)
                translation = dcm.XRayImageReceptorTranslation
                if len(translation) >= 3:
                    detector_translation = np.array([float(x) for x in translation[:3]])
            
            # Extract additional geometry parameters
            sad = float(getattr(dcm, 'RadiationMachineSAD', 1000.0))
            sid = float(getattr(dcm, 'RTImageSID', 1500.0))
            
            # Get pixel spacing
            pixel_spacing = getattr(dcm, 'ImagePlanePixelSpacing', [0.172, 0.172])
            pixel_size = float(pixel_spacing[0]) if len(pixel_spacing) > 0 else 0.172
            
            return {
                'dose_data': dose_data,
                'angle': angle,
                'detector_translation': detector_translation,
                'sad': sad,
                'sid': sid,
                'pixel_size': pixel_size
            }
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            return None
    
    results = []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        with tqdm(total=len(file_paths), desc="Loading") as pbar:
            futures = {executor.submit(read_epid_enhanced, path): path for path in file_paths}
            for future in futures:
                result = future.result()
                if result is not None:
                    results.append(result)
                pbar.update(1)
    
    return results

def extract_dose_increments_enhanced(epid_data):
    """Extract dose increments with geometry tracking."""
    n_images = len(epid_data)
    
    # Sort by gantry angle
    epid_data.sort(key=lambda x: x['angle'])
    
    dose_images = np.array([data['dose_data'] for data in epid_data])
    angles = np.array([data['angle'] for data in epid_data])
    
    # Extract dose increments
    dose_increments = np.zeros_like(dose_images)
    dose_increments[0] = dose_images[0].copy()
    
    for i in range(1, n_images):
        increment = dose_images[i] - dose_images[i-1]
        dose_increments[i] = np.maximum(increment, 0)
    
    # Extract geometry data
    detector_translations = []
    sads = []
    sids = []
    pixel_sizes = []
    
    for data in epid_data:
        detector_translations.append(data['detector_translation'])
        sads.append(data['sad'])
        sids.append(data['sid'])
        pixel_sizes.append(data['pixel_size'])
    
    return {
        'dose_increments': dose_increments,
        'angles': angles,
        'detector_translations': detector_translations,
        'sads': np.array(sads),
        'sids': np.array(sids),
        'pixel_sizes': np.array(pixel_sizes)
    }

# %% [markdown]
"""
## Enhanced Reconstruction with Detector Translation Correction
"""

# %%
def apply_detector_translation_correction(pixel_coords, detector_translation):
    """Apply detector translation correction to pixel coordinates."""
    if detector_translation is None:
        return pixel_coords
    
    # Apply translation in detector coordinate system
    # detector_translation is [lateral, longitudinal, normal] in mm
    corrected_coords = pixel_coords.copy()
    corrected_coords[0] += detector_translation[0]  # lateral shift
    corrected_coords[1] += detector_translation[1]  # longitudinal shift
    
    return corrected_coords

if NUMBA_AVAILABLE:
    @jit(nopython=True)
    def enhanced_dose_weighting(dose_increment, sid, pixel_size, translation_x=0.0, translation_y=0.0):
        """Enhanced dose weighting with detector translation correction."""
        nrows, ncols = dose_increment.shape
        col_coords = pixel_size * (np.arange(ncols) - ncols/2 + 0.5) + translation_x
        row_coords = pixel_size * (np.arange(nrows) - nrows/2 + 0.5) + translation_y
        weighted_dose = np.zeros_like(dose_increment)
        
        for i in prange(nrows):
            for j in prange(ncols):
                r_sq = row_coords[i]**2 + col_coords[j]**2
                distance = np.sqrt(sid**2 + r_sq)
                weight = (sid / distance)**2
                weighted_dose[i, j] = dose_increment[i, j] * weight
        return weighted_dose
    
    @jit(nopython=True)
    def enhanced_backproject_dose(weighted_dose, angle_rad, sad, sid, pixel_size, 
                                 vol_size, vol_extent, translation_x=0.0, translation_y=0.0):
        """Enhanced backprojection with detector translation correction."""
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
                        # Project to detector with translation correction
                        det_x = (x * cos_angle + y * sin_angle) * (sid / source_to_voxel)
                        det_y = z * (sid / source_to_voxel)
                        
                        # Apply detector translation correction
                        det_x_corrected = det_x - translation_x
                        det_y_corrected = det_y - translation_y
                        
                        px = det_x_corrected / pixel_size + ncols / 2
                        py = det_y_corrected / pixel_size + nrows / 2
                        
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
    def enhanced_dose_weighting(dose_increment, sid, pixel_size, translation_x=0.0, translation_y=0.0):
        nrows, ncols = dose_increment.shape
        col_coords = pixel_size * (np.arange(ncols) - ncols/2 + 0.5) + translation_x
        row_coords = pixel_size * (np.arange(nrows) - nrows/2 + 0.5) + translation_y
        row_grid, col_grid = np.meshgrid(row_coords, col_coords, indexing='ij')
        distance = np.sqrt(sid**2 + row_grid**2 + col_grid**2)
        return dose_increment * (sid / distance)**2
    
    def enhanced_backproject_dose(weighted_dose, angle_rad, sad, sid, pixel_size, 
                                 vol_size, vol_extent, translation_x=0.0, translation_y=0.0):
        return np.zeros((vol_size, vol_size, vol_size))

print(f"Enhanced reconstruction functions loaded (Numba: {NUMBA_AVAILABLE})")

# %%
class EnhancedEPIDReconstructor:
    """Enhanced EPID reconstructor with detector translation correction."""
    
    def __init__(self, vol_size=100, vol_extent=50.0):
        self.vol_size = vol_size
        self.vol_extent = vol_extent
        self.translation_stats = {}
    
    def analyze_detector_translations(self, detector_translations):
        """Analyze detector translation patterns."""
        valid_translations = [t for t in detector_translations if t is not None]
        
        if not valid_translations:
            self.translation_stats = {
                'has_translations': False,
                'mean_translation': np.zeros(3),
                'max_translation': np.zeros(3),
                'std_translation': np.zeros(3)
            }
            return
        
        translations_array = np.array(valid_translations)
        
        self.translation_stats = {
            'has_translations': True,
            'count': len(valid_translations),
            'mean_translation': np.mean(translations_array, axis=0),
            'max_translation': np.max(np.abs(translations_array), axis=0),
            'std_translation': np.std(translations_array, axis=0),
            'range_translation': np.ptp(translations_array, axis=0)
        }
        
        print(f"Detector Translation Analysis:")
        print(f"  Valid translations: {len(valid_translations)}")
        print(f"  Mean translation: [{self.translation_stats['mean_translation'][0]:.4f}, "
              f"{self.translation_stats['mean_translation'][1]:.4f}, "
              f"{self.translation_stats['mean_translation'][2]:.4f}] mm")
        print(f"  Max |translation|: [{self.translation_stats['max_translation'][0]:.4f}, "
              f"{self.translation_stats['max_translation'][1]:.4f}, "
              f"{self.translation_stats['max_translation'][2]:.4f}] mm")
        print(f"  Translation std: [{self.translation_stats['std_translation'][0]:.4f}, "
              f"{self.translation_stats['std_translation'][1]:.4f}, "
              f"{self.translation_stats['std_translation'][2]:.4f}] mm")
    
    def reconstruct(self, reconstruction_data, chunk_size=10):
        """Enhanced reconstruction with detector translation correction."""
        dose_increments = reconstruction_data['dose_increments']
        angles = reconstruction_data['angles']
        detector_translations = reconstruction_data['detector_translations']
        sads = reconstruction_data['sads']
        sids = reconstruction_data['sids']
        pixel_sizes = reconstruction_data['pixel_sizes']
        
        # Analyze detector translations
        self.analyze_detector_translations(detector_translations)
        
        n_projections = len(dose_increments)
        angles_rad = angles * np.pi / 180.0
        dose_volume = np.zeros((self.vol_size, self.vol_size, self.vol_size), dtype=np.float32)
        
        # Use mean geometry parameters
        mean_sad = np.mean(sads)
        mean_sid = np.mean(sids)
        mean_pixel_size = np.mean(pixel_sizes)
        
        print(f"Using geometry: SAD={mean_sad:.1f}mm, SID={mean_sid:.1f}mm, PixelSize={mean_pixel_size:.3f}mm")
        print(f"Processing {n_projections} projections for complete reconstruction")
        
        chunks = [list(range(i, min(i + chunk_size, n_projections))) 
                 for i in range(0, n_projections, chunk_size)]
        
        for chunk_indices in tqdm(chunks, desc="Enhanced Reconstruction"):
            chunk_contribution = np.zeros_like(dose_volume)
            
            for proj_idx in chunk_indices:
                dose_increment = dose_increments[proj_idx]
                if np.max(dose_increment) <= 0:
                    continue
                
                # Get detector translation for this projection
                translation = detector_translations[proj_idx]
                translation_x = translation[0] if translation is not None else 0.0
                translation_y = translation[1] if translation is not None else 0.0
                
                # Apply enhanced weighting with translation correction
                weighted_dose = enhanced_dose_weighting(
                    dose_increment, sids[proj_idx], pixel_sizes[proj_idx],
                    translation_x, translation_y
                )
                
                if NUMBA_AVAILABLE:
                    dose_contrib = enhanced_backproject_dose(
                        weighted_dose, angles_rad[proj_idx], 
                        sads[proj_idx], sids[proj_idx], pixel_sizes[proj_idx],
                        self.vol_size, self.vol_extent,
                        translation_x, translation_y
                    )
                    chunk_contribution += dose_contrib
            
            dose_volume += chunk_contribution
        
        return dose_volume

# %% [markdown] 
"""
## Configuration and Enhanced Reconstruction
"""

# %%
# Configuration - Use all images for complete reconstruction
QUICK_TEST = False  # Set to False to include all images
MAX_FILES = None    # Include all available files
EVERY_NTH = 1       # Use every image, no skipping
VOL_SIZE = 100      # High resolution reconstruction

data_path = r"dataset/VMAT 2025 - 6. SIB COMPLEX TARGET/T1/873251691"

if os.path.exists(data_path):
    # Load enhanced data - all images
    print("Loading ALL EPID data with enhanced geometry extraction...")
    epid_data = load_epid_with_geometry(data_path, MAX_FILES, EVERY_NTH)
    print(f"Loaded {len(epid_data)} DICOM files for complete reconstruction")
    
    if epid_data:
        reconstruction_data = extract_dose_increments_enhanced(epid_data)
        
        # Enhanced reconstruction
        reconstructor = EnhancedEPIDReconstructor(vol_size=VOL_SIZE)
        
        start_time = time.time()
        dose_volume = reconstructor.reconstruct(reconstruction_data)
        recon_time = time.time() - start_time
        
        n_projections = len(reconstruction_data['dose_increments'])
        print(f"Complete enhanced reconstruction: {recon_time:.2f}s, {n_projections/recon_time:.1f} proj/s")
        print(f"Total projections processed: {n_projections}")
        print(f"Reconstruction volume: {VOL_SIZE}³ = {VOL_SIZE**3:,} voxels")
    else:
        print("No valid EPID data loaded")
else:
    print("Data path not found")

# %% [markdown]
"""
## Enhanced Visualization with Translation Analysis
"""

# %%
if 'dose_volume' in locals():
    # Enhanced visualization
    center = VOL_SIZE // 2
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # Cross sections
    im1 = axes[0,0].imshow(dose_volume[:, :, center].T, cmap='hot', origin='lower')
    axes[0,0].set_title('Axial (Enhanced)')
    plt.colorbar(im1, ax=axes[0,0])
    
    im2 = axes[0,1].imshow(dose_volume[center, :, :].T, cmap='hot', origin='lower')
    axes[0,1].set_title('Sagittal (Enhanced)')
    plt.colorbar(im2, ax=axes[0,1])
    
    im3 = axes[0,2].imshow(dose_volume[:, center, :].T, cmap='hot', origin='lower')
    axes[0,2].set_title('Coronal (Enhanced)')
    plt.colorbar(im3, ax=axes[0,2])
    
    # Profiles
    coords = np.arange(VOL_SIZE) - center
    axes[1,0].plot(coords, dose_volume[:, center, center], 'b-', label='X')
    axes[1,0].plot(coords, dose_volume[center, :, center], 'r-', label='Y')
    axes[1,0].plot(coords, dose_volume[center, center, :], 'g-', label='Z')
    axes[1,0].set_title('Dose Profiles')
    axes[1,0].legend()
    axes[1,0].grid(True)
    
    # Statistics
    axes[1,1].axis('off')
    stats = f"""Complete Enhanced Reconstruction:
Volume: {VOL_SIZE}³ ({VOL_SIZE**3:,} voxels)
Projections: {n_projections} (ALL)
Time: {recon_time:.2f}s
Rate: {n_projections/recon_time:.1f} proj/s

Dose Range:
Min: {dose_volume.min():.3f}
Max: {dose_volume.max():.3f}
Mean: {dose_volume.mean():.3f}
Std: {dose_volume.std():.3f}

Coverage: {100*np.count_nonzero(dose_volume)/dose_volume.size:.1f}%
Numba: {'Yes' if NUMBA_AVAILABLE else 'No'}
Translation Corrected: {'Yes' if reconstructor.translation_stats['has_translations'] else 'No'}"""
    
    axes[1,1].text(0.1, 0.9, stats, transform=axes[1,1].transAxes, 
                   fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    # Translation analysis
    axes[1,2].axis('off')
    if reconstructor.translation_stats['has_translations']:
        trans_stats = f"""Detector Translation Analysis:
Valid: {reconstructor.translation_stats['count']} projections

Mean Translation (mm):
X: {reconstructor.translation_stats['mean_translation'][0]:.4f}
Y: {reconstructor.translation_stats['mean_translation'][1]:.4f}
Z: {reconstructor.translation_stats['mean_translation'][2]:.4f}

Max |Translation| (mm):
X: {reconstructor.translation_stats['max_translation'][0]:.4f}
Y: {reconstructor.translation_stats['max_translation'][1]:.4f}
Z: {reconstructor.translation_stats['max_translation'][2]:.4f}

Translation Range (mm):
X: {reconstructor.translation_stats['range_translation'][0]:.4f}
Y: {reconstructor.translation_stats['range_translation'][1]:.4f}
Z: {reconstructor.translation_stats['range_translation'][2]:.4f}"""
    else:
        trans_stats = """Detector Translation Analysis:
No detector translations found
Using standard geometry"""
    
    axes[1,2].text(0.1, 0.9, trans_stats, transform=axes[1,2].transAxes, 
                   fontsize=9, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.show()

# %%
# Optional: Save enhanced reconstruction results
if 'dose_volume' in locals():
    output_dir = "results"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    timestamp = int(time.time())
    
    # Save volume
    np.save(f"{output_dir}/enhanced_dose_volume_{timestamp}.npy", dose_volume)
    
    # Save metadata
    metadata = {
        'vol_size': VOL_SIZE,
        'n_projections': n_projections,
        'reconstruction_time': recon_time,
        'translation_stats': reconstructor.translation_stats,
        'geometry': {
            'mean_sad': float(np.mean(reconstruction_data['sads'])),
            'mean_sid': float(np.mean(reconstruction_data['sids'])),
            'mean_pixel_size': float(np.mean(reconstruction_data['pixel_sizes']))
        }
    }
    
    import json
    with open(f"{output_dir}/enhanced_metadata_{timestamp}.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"Enhanced reconstruction saved to {output_dir}/")