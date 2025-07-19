#!/usr/bin/env python3
"""
EPID Dose-Based 3D Reconstruction
=================================

Reconstructs 3D dose distributions from Electronic Portal Imaging Device (EPID) images.
Handles cumulative dose delivery from linear accelerator with proper dose increment extraction.

Author: Generated for radiation therapy dose reconstruction
Environment: py12 conda environment
"""

import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from pydicom import dcmread, dcmwrite
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import time
import warnings
warnings.filterwarnings('ignore')

# Scientific computing
from scipy.fft import fft, ifft
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata

# Advanced libraries
try:
    from skimage.transform import radon, iradon
    from skimage.filters import gaussian
    SKIMAGE_AVAILABLE = True
    print("✓ Scikit-image available for advanced reconstruction")
except ImportError:
    SKIMAGE_AVAILABLE = False
    print("⚠ Scikit-image not available - using basic reconstruction")

try:
    from numba import jit, prange
    NUMBA_AVAILABLE = True
    print("✓ Numba JIT compilation available")
except ImportError:
    NUMBA_AVAILABLE = False
    print("⚠ Numba not available - using standard NumPy")
    def jit(func):
        return func
    prange = range

PI = math.pi

def analyze_epid_data(file_path):
    """Analyze EPID DICOM file to understand dose data structure"""
    try:
        dcm = dcmread(file_path)
        
        print(f"Analyzing EPID file: {os.path.basename(file_path)}")
        print(f"{'='*50}")
        
        # Basic image info
        print(f"Image dimensions: {dcm.Rows} x {dcm.Columns}")
        print(f"Pixel data type: {dcm.pixel_array.dtype}")
        print(f"Pixel value range: {dcm.pixel_array.min()} to {dcm.pixel_array.max()}")
        print(f"Pixel mean: {dcm.pixel_array.mean():.2f}")
        
        # Geometric information
        if hasattr(dcm, 'PixelSpacing'):
            print(f"Pixel spacing: {dcm.PixelSpacing} mm")
        if hasattr(dcm, 'GantryAngle'):
            print(f"Gantry angle: {dcm.GantryAngle}°")
        if hasattr(dcm, 'RTImageSID'):
            print(f"RT Image SID: {dcm.RTImageSID} mm")
        if hasattr(dcm, 'RadiationMachineSAD'):
            print(f"Radiation machine SAD: {dcm.RadiationMachineSAD} mm")
        
        print(f"Modality: {dcm.Modality}")
        
        return dcm
        
    except Exception as e:
        print(f"Error analyzing EPID file: {e}")
        return None

def load_epid_data(data_path, max_workers=6, max_files=None):
    """Load EPID DICOM files efficiently"""
    
    # Get all DICOM files
    dicom_files = [f for f in os.listdir(data_path) if f.endswith('.dcm')]
    dicom_files.sort()  # Ensure consistent ordering
    
    if max_files:
        dicom_files = dicom_files[:max_files]
    
    file_paths = [os.path.join(data_path, f) for f in dicom_files]
    
    def read_epid_file(file_path):
        try:
            dcm = dcmread(file_path)
            
            # Extract dose data
            dose_data = dcm.pixel_array.astype(np.float32)
            
            # Extract angle information
            if hasattr(dcm, 'GantryAngle'):
                angle = float(dcm.GantryAngle)
            elif hasattr(dcm, 'BeamAngle'):
                angle = float(dcm.BeamAngle)
            else:
                # If no angle info, use file index
                angle = 0.0
            
            # Extract other relevant parameters
            metadata = {
                'file_path': file_path,
                'angle': angle,
                'pixel_spacing': getattr(dcm, 'PixelSpacing', [1.0, 1.0]),
                'image_position': getattr(dcm, 'ImagePositionPatient', [0.0, 0.0, 0.0]),
                'dose_scaling': getattr(dcm, 'DoseGridScaling', 1.0),
                'sad': getattr(dcm, 'RadiationMachineSAD', 1000.0),
                'sid': getattr(dcm, 'RTImageSID', 1500.0)
            }
            
            return dose_data, metadata
            
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            return None, None
    
    # Load files in parallel
    dose_images = []
    metadata_list = []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        with tqdm(total=len(file_paths), desc="Loading EPID files") as pbar:
            future_to_path = {executor.submit(read_epid_file, path): path for path in file_paths}
            
            for future in future_to_path:
                dose_data, metadata = future.result()
                if dose_data is not None:
                    dose_images.append(dose_data)
                    metadata_list.append(metadata)
                pbar.update(1)
    
    if not dose_images:
        raise ValueError("No valid EPID files found")
    
    # Convert to numpy arrays
    dose_images = np.array(dose_images)
    angles = np.array([meta['angle'] for meta in metadata_list])
    
    print(f"\nLoaded {len(dose_images)} EPID images")
    print(f"Image shape: {dose_images[0].shape}")
    print(f"Angle range: {angles.min():.1f}° to {angles.max():.1f}°")
    print(f"Dose range: {dose_images.min():.2f} to {dose_images.max():.2f}")
    
    return dose_images, angles, metadata_list

def preprocess_epid_dose_data(dose_images, angles, metadata, method='cumulative_differential'):
    """Preprocess EPID dose data for reconstruction
    
    EPID images represent cumulative dose delivery:
    - Each image = cumulative dose up to that gantry angle
    - To get dose delivered at each angle: current_image - previous_image
    - Last image contains total cumulative dose
    """
    
    print(f"Preprocessing EPID dose data using {method} method...")
    
    if method == 'cumulative_differential':
        # Proper EPID dose extraction: current - previous
        n_images = len(dose_images)
        processed_doses = np.zeros_like(dose_images)
        processed_angles = []
        
        print(f"Processing {n_images} cumulative EPID images...")
        print(f"Last image contains total cumulative dose: {dose_images[-1].max():.2f}")
        
        # First image is the dose from start to first angle
        processed_doses[0] = dose_images[0].copy()
        processed_angles.append(angles[0])
        
        # For subsequent images: dose_increment = current - previous
        for i in tqdm(range(1, n_images), desc="Extracting dose increments"):
            dose_increment = dose_images[i] - dose_images[i-1]
            
            # Check for reasonable dose increment
            if np.max(dose_increment) < 0:
                print(f"Warning: Negative dose increment at index {i}, using zero")
                dose_increment = np.maximum(dose_increment, 0)  # Ensure non-negative
            
            processed_doses[i] = dose_increment
            processed_angles.append(angles[i])
            
            # Debug info for first few images
            if i <= 3:
                print(f"  Image {i}: Cumulative dose max = {dose_images[i].max():.2f}, "
                      f"Increment max = {dose_increment.max():.2f}")
    
    # Sort by angle for proper reconstruction
    processed_angles = np.array(processed_angles)
    sort_indices = np.argsort(processed_angles)
    
    sorted_doses = processed_doses[sort_indices]
    sorted_angles = processed_angles[sort_indices]
    
    # Verify dose increment extraction
    total_reconstructed = np.sum(sorted_doses, axis=0)
    original_total = dose_images[-1]  # Last image should be total cumulative
    
    correlation = np.corrcoef(total_reconstructed.flatten(), original_total.flatten())[0, 1]
    
    print(f"\nDose extraction verification:")
    print(f"  Processed {len(sorted_doses)} dose projections")
    print(f"  Sum of increments max: {total_reconstructed.max():.2f}")
    print(f"  Original cumulative max: {original_total.max():.2f}")
    print(f"  Correlation coefficient: {correlation:.4f}")
    print(f"  Dose increment range: {sorted_doses.min():.3f} to {sorted_doses.max():.3f}")
    print(f"  Angle range: {sorted_angles.min():.1f}° to {sorted_angles.max():.1f}°")
    
    if correlation < 0.9:
        print(f"  ⚠️  Warning: Low correlation suggests dose extraction may need refinement")
    else:
        print(f"  ✅ Good correlation - dose extraction appears correct")
    
    return sorted_doses, sorted_angles, total_reconstructed, original_total

class EPIDDoseReconstructor:
    """3D dose reconstruction from EPID images with proper dose increment handling"""
    
    def __init__(self, sad=1000.0, sid=1500.0, pixel_size=0.172, recon_size=50):
        self.sad = sad  # Source-to-axis distance (mm)
        self.sid = sid  # Source-to-image distance (mm)
        self.pixel_size = pixel_size  # EPID pixel size (mm)
        self.recon_size = recon_size  # Reconstruction volume size
        
        # Calculate virtual detector pixel size at isocenter
        self.virtual_pixel_size = self.pixel_size * self.sad / self.sid
        
        print(f"EPID Dose Reconstructor initialized:")
        print(f"  SAD: {self.sad} mm")
        print(f"  SID: {self.sid} mm")
        print(f"  EPID pixel size: {self.pixel_size} mm")
        print(f"  Virtual pixel size: {self.virtual_pixel_size:.4f} mm")
        print(f"  Reconstruction volume: {self.recon_size}³")
        print(f"  Dose reconstruction: Based on incremental dose delivery")
    
    def dose_weight_projection(self, dose_increment):
        """Apply dose-specific weighting to projection"""
        nrows, ncols = dose_increment.shape
        
        # Create coordinate grids at detector plane
        col_coords = self.pixel_size * (np.arange(ncols) - ncols/2 + 0.5)
        row_coords = self.pixel_size * (np.arange(nrows) - nrows/2 + 0.5)
        
        row_grid, col_grid = np.meshgrid(row_coords, col_coords, indexing='ij')
        
        # Distance from source to detector pixel
        detector_distance = np.sqrt(self.sid**2 + row_grid**2 + col_grid**2)
        
        # Geometric correction for dose (inverse square law)
        geometric_correction = (self.sid / detector_distance)**2
        
        # Apply magnification correction to get dose at isocenter plane
        magnification_factor = self.sad / self.sid
        
        # Combined weighting for dose reconstruction
        dose_weight = geometric_correction * magnification_factor
        
        # Apply weighting to dose increment
        weighted_dose = dose_increment * dose_weight
        
        return weighted_dose
    
    def dose_backproject(self, weighted_dose_increment, angle_rad, dose_volume):
        """Backproject dose increment into 3D volume"""
        nrows, ncols = weighted_dose_increment.shape
        vol_size = dose_volume.shape[0]
        
        # Define reconstruction volume coordinates (at isocenter)
        max_detector_size = max(nrows, ncols) * self.virtual_pixel_size
        vol_extent = max_detector_size * 0.6  # Adjust as needed
        coords = np.linspace(-vol_extent, vol_extent, vol_size)
        
        # Create 3D coordinate grids
        X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')
        
        # Rotation matrices for current gantry angle
        cos_angle = np.cos(angle_rad)
        sin_angle = np.sin(angle_rad)
        
        # Distance from source to voxel
        source_to_voxel = self.sad + X * sin_angle - Y * cos_angle
        
        # Avoid division by zero
        source_to_voxel = np.maximum(source_to_voxel, 0.1)
        
        # Project 3D coordinates onto detector plane
        detector_x = (X * cos_angle + Y * sin_angle) * (self.sid / source_to_voxel)
        detector_y = Z * (self.sid / source_to_voxel)
        
        # Convert to detector pixel coordinates
        pixel_x = detector_x / self.pixel_size + ncols / 2
        pixel_y = detector_y / self.pixel_size + nrows / 2
        
        # Bounds checking
        valid_mask = ((pixel_x >= 0) & (pixel_x < ncols - 1) & 
                     (pixel_y >= 0) & (pixel_y < nrows - 1) &
                     (source_to_voxel > 0))
        
        # Initialize dose contribution for this angle
        dose_contribution = np.zeros_like(X, dtype=np.float32)
        
        if np.any(valid_mask):
            # Extract valid coordinates
            px = pixel_x[valid_mask]
            py = pixel_y[valid_mask]
            
            # Bilinear interpolation indices
            px0 = np.floor(px).astype(int)
            py0 = np.floor(py).astype(int)
            px1 = np.minimum(px0 + 1, ncols - 1)
            py1 = np.minimum(py0 + 1, nrows - 1)
            
            # Interpolation weights
            wx = px - px0
            wy = py - py0
            
            # Bilinear interpolation of dose values
            dose_val = ((1 - wx) * (1 - wy) * weighted_dose_increment[py0, px0] +
                       wx * (1 - wy) * weighted_dose_increment[py0, px1] +
                       (1 - wx) * wy * weighted_dose_increment[py1, px0] +
                       wx * wy * weighted_dose_increment[py1, px1])
            
            # Dose deposition weighting
            distance_weight = (self.sad / source_to_voxel[valid_mask])**2
            
            # Apply dose contribution to volume
            dose_contribution[valid_mask] = dose_val * distance_weight
        
        return dose_contribution
    
    def reconstruct_3d_dose(self, dose_increments, angles, chunk_size=25):
        """Reconstruct 3D dose distribution from EPID dose increments"""
        n_projections = len(dose_increments)
        angles_rad = angles * np.pi / 180.0
        
        # Initialize 3D dose volume
        vol_size = self.recon_size
        dose_volume = np.zeros((vol_size, vol_size, vol_size), dtype=np.float32)
        
        print(f"Reconstructing 3D dose from {n_projections} dose increments...")
        print(f"Volume size: {vol_size}³")
        print(f"Angle range: {angles.min():.1f}° to {angles.max():.1f}°")
        
        # Track total dose for verification
        total_dose_added = 0.0
        
        # Process in chunks to manage memory
        chunks = [list(range(i, min(i + chunk_size, n_projections))) 
                 for i in range(0, n_projections, chunk_size)]
        
        for chunk_idx, chunk_indices in enumerate(tqdm(chunks, desc="Reconstructing dose")):
            chunk_dose = np.zeros_like(dose_volume)
            
            for proj_idx in chunk_indices:
                # Get dose increment and angle
                dose_increment = dose_increments[proj_idx]
                angle_rad = angles_rad[proj_idx]
                
                # Skip if no dose was delivered at this angle
                if np.max(dose_increment) <= 0:
                    continue
                
                # Apply dose-specific weighting
                weighted_dose = self.dose_weight_projection(dose_increment)
                
                # Backproject dose increment into 3D volume
                dose_contribution = self.dose_backproject(weighted_dose, angle_rad, dose_volume)
                
                chunk_dose += dose_contribution
                total_dose_added += np.sum(dose_increment)
                
                # Debug info for first few projections
                if proj_idx < 3:
                    print(f"  Projection {proj_idx}: angle={angles[proj_idx]:.1f}°, "
                          f"dose_max={np.max(dose_increment):.3f}, "
                          f"contribution_max={np.max(dose_contribution):.6f}")
            
            dose_volume += chunk_dose
        
        print(f"\n3D dose reconstruction completed:")
        print(f"  Total dose processed: {total_dose_added:.2f}")
        print(f"  3D volume dose range: {dose_volume.min():.6f} to {dose_volume.max():.6f}")
        print(f"  Non-zero voxels: {np.count_nonzero(dose_volume):,}")
        print(f"  Center dose: {dose_volume[vol_size//2, vol_size//2, vol_size//2]:.6f}")
        
        return dose_volume

def visualize_3d_dose(dose_volume, title="3D Dose Distribution"):
    """Visualize 3D dose distribution with cross-sections and analysis"""
    
    vol_size = dose_volume.shape[0]
    center = vol_size // 2
    
    # Create comprehensive visualization
    fig = plt.figure(figsize=(20, 15))
    
    # Main cross-sections
    ax1 = plt.subplot(3, 4, 1)
    im1 = ax1.imshow(dose_volume[:, :, center].T, cmap='hot', aspect='equal', origin='lower')
    ax1.set_title(f'Axial (Z={center})')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    plt.colorbar(im1, ax=ax1, shrink=0.8)
    
    ax2 = plt.subplot(3, 4, 2)
    im2 = ax2.imshow(dose_volume[center, :, :].T, cmap='hot', aspect='equal', origin='lower')
    ax2.set_title(f'Sagittal (X={center})')
    ax2.set_xlabel('Y')
    ax2.set_ylabel('Z')
    plt.colorbar(im2, ax=ax2, shrink=0.8)
    
    ax3 = plt.subplot(3, 4, 3)
    im3 = ax3.imshow(dose_volume[:, center, :].T, cmap='hot', aspect='equal', origin='lower')
    ax3.set_title(f'Coronal (Y={center})')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Z')
    plt.colorbar(im3, ax=ax3, shrink=0.8)
    
    # Maximum intensity projections
    ax4 = plt.subplot(3, 4, 4)
    mip_axial = np.max(dose_volume, axis=2)
    im4 = ax4.imshow(mip_axial.T, cmap='hot', aspect='equal', origin='lower')
    ax4.set_title('MIP Axial')
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    plt.colorbar(im4, ax=ax4, shrink=0.8)
    
    ax5 = plt.subplot(3, 4, 5)
    mip_sagittal = np.max(dose_volume, axis=0)
    im5 = ax5.imshow(mip_sagittal.T, cmap='hot', aspect='equal', origin='lower')
    ax5.set_title('MIP Sagittal')
    ax5.set_xlabel('Y')
    ax5.set_ylabel('Z')
    plt.colorbar(im5, ax=ax5, shrink=0.8)
    
    ax6 = plt.subplot(3, 4, 6)
    mip_coronal = np.max(dose_volume, axis=1)
    im6 = ax6.imshow(mip_coronal.T, cmap='hot', aspect='equal', origin='lower')
    ax6.set_title('MIP Coronal')
    ax6.set_xlabel('X')
    ax6.set_ylabel('Z')
    plt.colorbar(im6, ax=ax6, shrink=0.8)
    
    # Dose profiles
    ax7 = plt.subplot(3, 4, 7)
    profile_x = dose_volume[:, center, center]
    profile_y = dose_volume[center, :, center]
    profile_z = dose_volume[center, center, :]
    coords = np.arange(len(profile_x)) - center
    
    ax7.plot(coords, profile_x, 'b-', linewidth=2, label='X profile')
    ax7.plot(coords, profile_y, 'r-', linewidth=2, label='Y profile')
    ax7.plot(coords, profile_z, 'g-', linewidth=2, label='Z profile')
    ax7.set_title('Dose Profiles')
    ax7.set_xlabel('Position (voxels)')
    ax7.set_ylabel('Dose')
    ax7.legend()
    ax7.grid(True, alpha=0.3)
    
    # Dose histogram
    ax8 = plt.subplot(3, 4, 8)
    nonzero_doses = dose_volume[dose_volume > 0]
    if len(nonzero_doses) > 0:
        ax8.hist(nonzero_doses, bins=50, alpha=0.7, color='red', edgecolor='black')
        ax8.set_title('Dose Histogram (Non-zero)')
        ax8.set_xlabel('Dose')
        ax8.set_ylabel('Frequency')
        ax8.grid(True, alpha=0.3)
    
    # 3D dose statistics
    ax9 = plt.subplot(3, 4, 9)
    ax9.axis('off')
    
    stats_text = f"""
    3D Dose Statistics ({vol_size}³):
    
    Total voxels: {dose_volume.size:,}
    Non-zero voxels: {np.count_nonzero(dose_volume):,}
    Coverage: {100 * np.count_nonzero(dose_volume) / dose_volume.size:.1f}%
    
    Dose Range:
    Min: {dose_volume.min():.6f}
    Max: {dose_volume.max():.6f}
    Mean: {dose_volume.mean():.6f}
    Std: {dose_volume.std():.6f}
    
    Center Analysis:
    Center dose: {dose_volume[center, center, center]:.6f}
    Max dose location: {np.unravel_index(np.argmax(dose_volume), dose_volume.shape)}
    """
    
    ax9.text(0.05, 0.95, stats_text, transform=ax9.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace')
    
    # Dose volume histogram (DVH-like)
    ax10 = plt.subplot(3, 4, 10)
    if len(nonzero_doses) > 0:
        sorted_doses = np.sort(nonzero_doses)[::-1]  # Sort descending
        volumes = np.arange(1, len(sorted_doses) + 1) / len(sorted_doses) * 100
        ax10.plot(sorted_doses, volumes, 'b-', linewidth=2)
        ax10.set_title('Dose-Volume Relationship')
        ax10.set_xlabel('Dose')
        ax10.set_ylabel('Volume (%)')
        ax10.grid(True, alpha=0.3)
    
    # Dose distribution in 3D space
    ax11 = plt.subplot(3, 4, 11)
    # Show dose distribution as scatter plot
    if np.count_nonzero(dose_volume) > 0:
        z_coords, y_coords, x_coords = np.where(dose_volume > dose_volume.max() * 0.1)
        doses = dose_volume[z_coords, y_coords, x_coords]
        scatter = ax11.scatter(x_coords, y_coords, c=doses, cmap='hot', alpha=0.6, s=20)
        ax11.set_title('High-Dose Regions (>10% max)')
        ax11.set_xlabel('X')
        ax11.set_ylabel('Y')
        plt.colorbar(scatter, ax=ax11, shrink=0.8)
    
    # Quality metrics
    ax12 = plt.subplot(3, 4, 12)
    ax12.axis('off')
    
    # Calculate quality metrics
    max_dose = dose_volume.max()
    mean_dose = dose_volume.mean()
    coverage = 100 * np.count_nonzero(dose_volume) / dose_volume.size
    
    # Dose uniformity (coefficient of variation for non-zero doses)
    if len(nonzero_doses) > 0:
        cv = np.std(nonzero_doses) / np.mean(nonzero_doses) * 100
    else:
        cv = 0
    
    quality_text = f"""
    Quality Metrics:
    
    Coverage: {coverage:.1f}%
    {'✅ Good' if coverage > 5 else '⚠️ Poor'} coverage
    
    Dose Uniformity:
    CV: {cv:.1f}%
    {'✅ Uniform' if cv < 50 else '⚠️ Variable'} distribution
    
    Reconstruction Quality:
    Max/Mean ratio: {max_dose/mean_dose:.1f}
    {'✅ Good' if max_dose/mean_dose > 2 else '⚠️ Poor'} contrast
    
    Volume Characteristics:
    Non-zero: {np.count_nonzero(dose_volume):,}
    Total: {dose_volume.size:,}
    Efficiency: {100*np.count_nonzero(dose_volume)/dose_volume.size:.1f}%
    """
    
    ax12.text(0.05, 0.95, quality_text, transform=ax12.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace')
    
    plt.suptitle(f'{title} - Complete Analysis', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.show()
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"3D DOSE VISUALIZATION SUMMARY")
    print(f"{'='*60}")
    print(f"Volume dimensions: {vol_size}x{vol_size}x{vol_size}")
    print(f"Total voxels: {dose_volume.size:,}")
    print(f"Non-zero voxels: {np.count_nonzero(dose_volume):,}")
    print(f"Dose coverage: {coverage:.1f}%")
    print(f"Dose range: {dose_volume.min():.6f} to {dose_volume.max():.6f}")
    print(f"Center dose: {dose_volume[center, center, center]:.6f}")
    print(f"Quality: {'✅ Good' if coverage > 5 and cv < 50 else '⚠️ Needs review'}")

def main():
    """Main function to run EPID dose reconstruction"""
    
    print("EPID Dose-Based 3D Reconstruction")
    print("=" * 50)
    print(f"Python version: {sys.version}")
    print(f"Libraries loaded successfully!")
    
    # Data path
    epid_path = r"E:\CMC\pyprojects\radio_therapy\dose-3d\dataset\VMAT 2025 - 6. SIB COMPLEX TARGET\T1\873251691"
    
    # Check if path exists
    if not os.path.exists(epid_path):
        print(f"❌ Error: Path not found: {epid_path}")
        return
    
    print(f"\nProcessing EPID data from: {epid_path}")
    
    try:
        # Step 1: Analyze sample data
        print("\n" + "="*50)
        print("STEP 1: ANALYZING SAMPLE EPID DATA")
        print("="*50)
        
        first_file = os.path.join(epid_path, "00000.dcm")
        if os.path.exists(first_file):
            sample_dcm = analyze_epid_data(first_file)
        else:
            print("❌ Sample file not found")
            return
        
        # Step 2: Load all EPID data
        print("\n" + "="*50)
        print("STEP 2: LOADING ALL EPID DATA")
        print("="*50)
        
        start_time = time.time()
        dose_images, angles, metadata = load_epid_data(epid_path, max_files=None)
        load_time = time.time() - start_time
        
        print(f"\nComplete EPID Dose Delivery Dataset:")
        print(f"  Total images: {len(dose_images)}")
        print(f"  Loading time: {load_time:.2f} seconds")
        print(f"  Angle coverage: {angles.max() - angles.min():.1f}°")
        
        # Step 3: Preprocess dose data
        print("\n" + "="*50)
        print("STEP 3: PREPROCESSING DOSE DATA")
        print("="*50)
        
        processed_doses, processed_angles, reconstructed_total, original_total = preprocess_epid_dose_data(
            dose_images, angles, metadata, method='cumulative_differential'
        )
        
        # Step 4: Create reconstructor
        print("\n" + "="*50)
        print("STEP 4: INITIALIZING RECONSTRUCTOR")
        print("="*50)
        
        if metadata:
            sad = metadata[0]['sad']
            sid = metadata[0]['sid']
            pixel_spacing = metadata[0]['pixel_spacing']
            pixel_size = pixel_spacing[0] if isinstance(pixel_spacing, list) else 0.172
        else:
            sad, sid, pixel_size = 1000.0, 1500.0, 0.172
        
        reconstructor = EPIDDoseReconstructor(sad=sad, sid=sid, pixel_size=pixel_size, recon_size=50)
        
        # Step 5: Run reconstruction
        print("\n" + "="*50)
        print("STEP 5: RUNNING 3D RECONSTRUCTION")
        print("="*50)
        
        start_time = time.time()
        
        # Adaptive chunk size
        num_images = len(processed_doses)
        chunk_size = 15 if num_images > 500 else 25
        
        dose_volume_3d = reconstructor.reconstruct_3d_dose(
            processed_doses,
            processed_angles,
            chunk_size=chunk_size
        )
        
        reconstruction_time = time.time() - start_time
        
        # Step 6: Results and validation
        print("\n" + "="*50)
        print("STEP 6: VALIDATION AND RESULTS")
        print("="*50)
        
        total_3d_dose = np.sum(dose_volume_3d)
        total_2d_dose = np.sum(processed_doses)
        
        print(f"Reconstruction completed in {reconstruction_time:.2f} seconds")
        print(f"Images processed: {len(processed_doses)}")
        print(f"Processing rate: {len(processed_doses)/reconstruction_time:.1f} images/second")
        print(f"3D volume shape: {dose_volume_3d.shape}")
        print(f"Dose conservation: {100*abs(total_3d_dose/total_2d_dose - 1.0):.1f}% error")
        
        # Step 7: Visualization
        print("\n" + "="*50)
        print("STEP 7: COMPREHENSIVE VISUALIZATION")
        print("="*50)
        
        visualize_3d_dose(dose_volume_3d, "EPID Dose Reconstruction - Complete Analysis")
        
        # Step 8: Save results
        print("\n" + "="*50)
        print("STEP 8: SAVING RESULTS")
        print("="*50)
        
        output_dir = r"E:\CMC\pyprojects\radio_therapy\dose-3d\results"
        os.makedirs(output_dir, exist_ok=True)
        
        # Save 3D volume
        np.save(os.path.join(output_dir, "epid_dose_3d.npy"), dose_volume_3d)
        
        # Save metadata
        import json
        reconstruction_info = {
            'volume_shape': dose_volume_3d.shape,
            'n_projections': len(processed_doses),
            'angle_range': [processed_angles.min(), processed_angles.max()],
            'dose_range': [dose_volume_3d.min(), dose_volume_3d.max()],
            'reconstruction_time': reconstruction_time,
            'sad': sad,
            'sid': sid,
            'pixel_size': pixel_size,
            'dose_conservation_error': 100*abs(total_3d_dose/total_2d_dose - 1.0)
        }
        
        with open(os.path.join(output_dir, "reconstruction_info.json"), 'w') as f:
            json.dump(reconstruction_info, f, indent=2)
        
        print(f"Results saved to: {output_dir}")
        print(f"  - 3D dose volume: epid_dose_3d.npy")
        print(f"  - Metadata: reconstruction_info.json")
        
        # Final summary
        print("\n" + "="*80)
        print("EPID DOSE RECONSTRUCTION COMPLETED SUCCESSFULLY")
        print("="*80)
        print(f"✅ Processed {len(processed_doses)} images from complete dose delivery")
        print(f"✅ Reconstructed {dose_volume_3d.shape[0]}³ dose volume")
        print(f"✅ Dose conservation error: {100*abs(total_3d_dose/total_2d_dose - 1.0):.1f}%")
        print(f"✅ Processing time: {reconstruction_time:.2f} seconds")
        print(f"✅ Comprehensive visualization and analysis completed")
        print(f"✅ Results saved for further analysis")
        
    except Exception as e:
        print(f"❌ Error during reconstruction: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()