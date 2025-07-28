#!/usr/bin/env python3
"""
Filtered Backprojection EPID Dose Reconstruction

This program implements true filtered backprojection for EPID dose reconstruction
with configurable filter types and sizes.

Features:
- True filtered backprojection (backprojects filtered projection)
- Configurable filter types: Shepp-Logan, Ram-Lak, Cosine, Hamming
- Adjustable filter sizes (0.1 to 1.0 factor)
- Global mask for consistent processing
- Filter visualization capabilities
- DICOM I/O support

Author: Generated with Claude Code
Date: 2025
"""

import math
import os
import time
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from pydicom import dcmread, dcmwrite
from scipy.fft import fft, ifft
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import cv2

warnings.filterwarnings('ignore')

# Constants
PI = math.pi

# Global mask storage
GLOBAL_MASK = None

class FilteredBackprojectionReconstructor:
    """
    True Filtered Backprojection reconstructor with configurable filters
    
    This class implements filtered backprojection where the filtered projection
    (not the weighted projection) is backprojected into the 3D volume.
    """
    
    def __init__(self, SOD=1000.0, delta_dd=0.1075, Nimage=100):
        """
        Initialize the reconstructor
        
        Args:
            SOD: Source-to-Object Distance (mm)
            delta_dd: Detector pixel spacing (mm)
            Nimage: Reconstruction image size (voxels)
        """
        self.SOD = SOD
        self.delta_dd = delta_dd
        self.Nimage = Nimage
        self.global_mask_initialized = False
        
    def filter_SL_configurable(self, N, d, filter_type='shepp_logan', filter_size_factor=1.0):
        """
        Generate configurable reconstruction filters
        
        Args:
            N: Filter length
            d: Pixel spacing
            filter_type: 'shepp_logan', 'ram_lak', 'cosine', 'hamming'
            filter_size_factor: 0.1 to 1.0 (controls filter extent)
            
        Returns:
            fh: Filter array
        """
        fh = np.zeros(N)
        
        # Determine effective filter size
        effective_N = int(N * filter_size_factor)
        center = N // 2
        start_idx = max(0, center - effective_N // 2)
        end_idx = min(N, center + effective_N // 2)
        
        if filter_type == 'shepp_logan':
            for k1 in range(start_idx, end_idx):
                denominator = 4*(k1-N/2.0)**2-1
                if abs(denominator) > 1e-10:
                    fh[k1] = -2.0/(PI*PI*d*d*denominator)
        
        elif filter_type == 'ram_lak':
            # Ramp filter (Ram-Lak)
            for k1 in range(start_idx, end_idx):
                if k1 == center:
                    fh[k1] = 1.0 / (4.0 * d * d)
                else:
                    fh[k1] = -1.0 / (PI * PI * d * d * (k1 - center) ** 2)
        
        elif filter_type == 'cosine':
            # Cosine filter
            for k1 in range(start_idx, end_idx):
                if k1 == center:
                    fh[k1] = 1.0 / (4.0 * d * d)
                else:
                    freq = (k1 - center) / (N * d)
                    fh[k1] = -1.0 / (PI * PI * d * d * (k1 - center) ** 2) * np.cos(PI * freq * d)
        
        elif filter_type == 'hamming':
            # Hamming windowed filter
            for k1 in range(start_idx, end_idx):
                if k1 == center:
                    fh[k1] = 1.0 / (4.0 * d * d)
                else:
                    base_filter = -1.0 / (PI * PI * d * d * (k1 - center) ** 2)
                    window_val = 0.54 + 0.46 * np.cos(2 * PI * (k1 - center) / effective_N)
                    fh[k1] = base_filter * window_val
        
        return fh
    
    def nearestPowerOf2(self, N):
        """Find nearest power of 2 >= N"""
        a = int(math.log2(N))
        if 2**a == N:
            return N
        return 2**(a + 1)
    
    def weight_projection(self, projection_beta, SOD, delta_dd):
        """Apply geometric weighting to projection"""
        Nrows, Ncolumns = projection_beta.shape
        dd_column = delta_dd*np.arange(-Ncolumns/2+0.5, (Ncolumns/2+1)-0.5, 1.0)
        dd_row = delta_dd*np.arange(-Nrows/2+0.5, (Nrows/2+1)-0.5, 1.0)
        dd_row2D, dd_column2D = np.meshgrid(dd_row, dd_column, indexing='ij')
        weighted_projection = projection_beta*SOD/np.sqrt(SOD*SOD+np.power(dd_row2D, 2.0)+np.power(dd_column2D, 2.0))
        return weighted_projection
    
    def convolve_with_filter(self, projection, filter_kernel):
        """Convolve projection with filter kernel"""
        Nrows, Ncolumns = projection.shape
        Nfft = self.nearestPowerOf2(2 * Ncolumns - 1)
        
        # Ensure filter is the right size
        if len(filter_kernel) != Nfft:
            fh_padded = np.zeros(Nfft)
            fh_padded[:min(len(filter_kernel), Nfft)] = filter_kernel[:min(len(filter_kernel), Nfft)]
            filter_kernel = fh_padded
        
        fh_fft = fft(filter_kernel / 2.0)
        
        projection_padded = np.zeros((Nrows, Nfft))
        projection_padded[:, :Ncolumns] = projection

        projection_fft = fft(projection_padded, axis=1)
        convoluted_freq = projection_fft * fh_fft
        convoluted_time = ifft(convoluted_freq, axis=1).real
        filtered_projection = convoluted_time[:, :Ncolumns]
        
        return filtered_projection
    
    def initialize_global_mask(self, projections):
        """Initialize the global mask once for all projections"""
        global GLOBAL_MASK
        
        Nrows, Ncolumns = projections.shape[1], projections.shape[2]
        MX, MZ = self.Nimage, int(self.Nimage * Nrows / Ncolumns)
        
        # Create coordinate grids
        roi = self.delta_dd*np.array([-Ncolumns/2.0+0.5, Ncolumns/2.0-0.5, -Nrows/2.0+0.5, Nrows/2.0-0.5])
        hx = (roi[1]-roi[0])/(MX-1)
        xrange = roi[0]+hx*np.arange(0, MX)
        hy = (roi[3]-roi[2])/(MZ-1)
        yrange = roi[2]+hy*np.arange(0, MZ)
        XX, YY, ZZ = np.meshgrid(xrange, xrange, yrange, indexing='ij')
        
        # Compute conservative mask
        U = (self.SOD + XX*np.sin(0) - YY*np.cos(0))/self.SOD
        a = (XX*np.cos(0) + YY*np.sin(0))/U
        xx = np.int32(np.floor(a/self.delta_dd))
        b = ZZ/U
        yy = np.int32(np.floor(b/self.delta_dd))
        xx = xx + int(Ncolumns/2)
        yy = yy + int(Nrows/2)
        
        # Create mask with safety margin
        margin = 5
        GLOBAL_MASK = np.where((xx >= margin) & (xx < Ncolumns-1-margin) & 
                              (yy >= margin) & (yy < Nrows-1-margin))
        
        self.global_mask_initialized = True
        print(f"Global mask initialized with {len(GLOBAL_MASK[0])} valid voxels")
        return GLOBAL_MASK
    
    def backproject_with_global_mask(self, filtered_projection, SOD, beta_num, beta_m, delta_dd, Nimage):
        """Backproject using global mask"""
        global GLOBAL_MASK
        
        Nrows, Ncolumns = filtered_projection.shape
        MX, MZ = Nimage, int(Nimage*Nrows/Ncolumns)
        
        roi = delta_dd*np.array([-Ncolumns/2.0+0.5, Ncolumns/2.0-0.5, -Nrows/2.0+0.5, Nrows/2.0-0.5])
        hx = (roi[1]-roi[0])/(MX-1)
        xrange = roi[0]+hx*np.arange(0, MX)
        hy = (roi[3]-roi[2])/(MZ-1)
        yrange = roi[2]+hy*np.arange(0, MZ)
        XX, YY, ZZ = np.meshgrid(xrange, xrange, yrange, indexing='ij')
        temp_rec = np.zeros((MX, MX, MZ))
        
        U = (SOD+XX*np.sin(beta_m)-YY*np.cos(beta_m))/SOD
        a = (XX*np.cos(beta_m)+YY*np.sin(beta_m))/U
        xx = np.int32(np.floor(a/delta_dd))
        u1 = a/delta_dd-xx
        b = ZZ/U
        yy = np.int32(np.floor(b/delta_dd))
        u2 = b/delta_dd-yy
        xx = xx+int(Ncolumns/2)
        yy = yy+int(Nrows/2)

        if GLOBAL_MASK is not None:
            mask = GLOBAL_MASK
        else:
            mask = np.where((xx >= 0) & (xx < Ncolumns-1) & (yy >= 0) & (yy < Nrows-1))
        
        xx_masked = xx[mask]
        yy_masked = yy[mask]
        u1_masked = u1[mask]
        u2_masked = u2[mask]
        
        # Ensure indices are within bounds
        valid_indices = (xx_masked >= 0) & (xx_masked < Ncolumns-1) & (yy_masked >= 0) & (yy_masked < Nrows-1)
        xx_masked = xx_masked[valid_indices]
        yy_masked = yy_masked[valid_indices]
        u1_masked = u1_masked[valid_indices]
        u2_masked = u2_masked[valid_indices]
        
        if len(xx_masked) > 0:
            temp = ((1-u1_masked)*(1-u2_masked)*filtered_projection[yy_masked, xx_masked]+
                    (1-u1_masked)*u2_masked*filtered_projection[yy_masked+1, xx_masked]+
                    (1-u2_masked)*u1_masked*filtered_projection[yy_masked, xx_masked+1]+
                    u1_masked*u2_masked*filtered_projection[yy_masked+1, xx_masked+1])
            
            final_mask = tuple(arr[valid_indices] for arr in mask)
            temp_rec[final_mask] = temp_rec[final_mask] + temp/(np.power(U[final_mask], 2))*2*PI/beta_num
        
        return temp_rec
    
    def reconstruct_filtered_backprojection(self, projections, angles, 
                                          filter_type='shepp_logan',
                                          filter_size_factor=0.5,
                                          use_weighting=True,
                                          chunk_size=25, 
                                          use_global_mask=True,
                                          verbose=True):
        """
        Perform filtered backprojection reconstruction
        
        Args:
            projections: 3D array of projection images
            angles: Array of projection angles (degrees) 
            filter_type: 'shepp_logan', 'ram_lak', 'cosine', 'hamming'
            filter_size_factor: 0.1 to 1.0 (controls filter extent)
            use_weighting: Apply geometric weighting before filtering
            chunk_size: Number of projections to process at once
            use_global_mask: Use consistent mask across projections
            verbose: Print progress information
            
        Returns:
            rec_image: 3D reconstructed volume
        """
        n_angles = len(angles)
        beta_rad = angles * PI / 180.0
        
        # Initialize global mask if requested
        if use_global_mask and not self.global_mask_initialized:
            self.initialize_global_mask(projections)
        
        # Prepare filter
        Ncolumns = projections.shape[2]
        Nrows = projections.shape[1]
        Nfft = self.nearestPowerOf2(2 * Ncolumns - 1)
        
        filter_kernel = self.filter_SL_configurable(
            Nfft, self.delta_dd, 
            filter_type=filter_type, 
            filter_size_factor=filter_size_factor
        )
        
        if verbose:
            print(f"FILTERED BACKPROJECTION RECONSTRUCTION:")
            print(f"  Projections: {n_angles}")
            print(f"  Filter type: {filter_type}")
            print(f"  Filter size factor: {filter_size_factor}")
            print(f"  Filter non-zero elements: {np.count_nonzero(filter_kernel)}")
            print(f"  Filter range: {filter_kernel.min():.6f} to {filter_kernel.max():.6f}")
            print(f"  Use weighting: {use_weighting}")
            print(f"  Use global mask: {use_global_mask}")
        
        # Initialize result
        MX = self.Nimage
        MZ = int(self.Nimage * Nrows / Ncolumns)
        rec_image = np.zeros((MX, MX, MZ))
        
        # Process in chunks
        chunks = [list(range(i, min(i + chunk_size, n_angles))) 
                 for i in range(0, n_angles, chunk_size)]
        
        if verbose:
            print(f"  Processing in {len(chunks)} chunks...")
        
        for chunk_idx, chunk_indices in enumerate(tqdm(chunks, desc="Filtered backprojection", disable=not verbose)):
            chunk_result = np.zeros((MX, MX, MZ))
            
            for angle_idx in chunk_indices:
                projection_beta = projections[angle_idx, :, :]
                
                # Step 1: Optional geometric weighting
                if use_weighting:
                    input_projection = self.weight_projection(projection_beta, self.SOD, self.delta_dd)
                else:
                    input_projection = projection_beta
                
                # Step 2: Filter the projection
                filtered_projection = self.convolve_with_filter(input_projection, filter_kernel)
                
                # Step 3: FILTERED BACKPROJECTION - backproject the filtered data
                temp_rec = self.backproject_with_global_mask(
                    filtered_projection,  # Use filtered projection!
                    self.SOD, n_angles, beta_rad[angle_idx], 
                    self.delta_dd, self.Nimage
                )
                
                chunk_result += temp_rec
            
            rec_image += chunk_result
        
        return rec_image

def load_dicom_fast(path_list, max_workers=6):
    """Fast DICOM loading with threading"""
    def read_dicom_data(fname):
        try:
            dcm = dcmread(fname)
            return dcm.pixel_array.astype(np.uint16), float(dcm.GantryAngle)
        except Exception as e:
            print(f"Error reading {fname}: {e}")
            return None, None
    
    images, angles = [], []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        with tqdm(total=len(path_list), desc="Loading DICOM files") as pbar:
            future_to_path = {executor.submit(read_dicom_data, path): path for path in path_list}
            
            for future in future_to_path:
                img, angle = future.result()
                if img is not None:
                    images.append(img)
                    angles.append(angle)
                pbar.update(1)
    
    return np.array(images), np.array(angles)

def process_differential_unrotated(images, angles, threshold=10000):
    """Process differential images without rotation correction"""
    n_images = len(images)
    shape = images[0].shape
    processed_images = np.zeros((n_images, shape[0], shape[0]), dtype=np.uint16)
    processed_angles = []
    
    prev = np.zeros((shape[0], shape[0]), dtype=np.uint16)
    
    for idx in tqdm(range(n_images), desc="Processing differences"):
        curr = images[idx]
        _m = curr - prev
        
        if np.max(_m) > threshold:
            # Use previous processed image when threshold exceeded
            if idx > 0:
                processed_images[idx, :, :] = processed_images[idx-1, :, :]
                processed_angles.append(processed_angles[idx-1])
            else:
                processed_images[idx, :, :] = curr - prev
                processed_angles.append(angles[idx])
        else:
            # Normal differential processing
            processed_images[idx, :, :] = curr - prev
            processed_angles.append(angles[idx])
        
        prev = curr
    
    return processed_images, np.array(processed_angles)

def visualize_filters(delta_dd=0.1, N=512, save_path=None):
    """
    Visualize different filter types and sizes
    
    Args:
        delta_dd: Pixel spacing for filter calculation
        N: Filter length
        save_path: Optional path to save the plot
    """
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    fig.suptitle('Reconstruction Filter Comparison', fontsize=16)
    
    x_axis = np.arange(N) - N//2
    
    # Filter configurations to test
    filter_configs = [
        ('shepp_logan', 1.0), ('shepp_logan', 0.5), ('shepp_logan', 0.25),
        ('ram_lak', 0.5), ('cosine', 0.5), ('hamming', 0.5),
        ('shepp_logan', 0.1), ('ram_lak', 0.25)
    ]
    
    reconstructor = FilteredBackprojectionReconstructor(delta_dd=delta_dd)
    
    for i, (filter_type, size_factor) in enumerate(filter_configs):
        row = i // 4
        col = i % 4
        
        if row < 2 and col < 4:
            # Generate filter
            test_filter = reconstructor.filter_SL_configurable(
                N, delta_dd, filter_type=filter_type, filter_size_factor=size_factor
            )
            
            # Plot filter
            axes[row, col].plot(x_axis, test_filter, 'b-', linewidth=1.5)
            axes[row, col].set_title(f'{filter_type}\nSize: {size_factor}', fontsize=10)
            axes[row, col].grid(True, alpha=0.3)
            axes[row, col].set_xlim(-N//4, N//4)
            
            # Set y-axis limits
            if np.any(test_filter != 0):
                y_max = np.max(np.abs(test_filter[test_filter != 0]))
                axes[row, col].set_ylim(-y_max*1.1, y_max*1.1)
            
            # Add statistics
            non_zero = np.count_nonzero(test_filter)
            axes[row, col].text(0.02, 0.98, f'Non-zero: {non_zero}', 
                              transform=axes[row, col].transAxes, 
                              verticalalignment='top', fontsize=8,
                              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Remove empty subplots
    for i in range(len(filter_configs), 8):
        row = i // 4
        col = i % 4
        if row < 2 and col < 4:
            fig.delaxes(axes[row, col])
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Filter visualization saved to: {save_path}")
    
    plt.show()

def display_reconstruction_results(rec_image, title="Filtered Backprojection Results", save_path=None):
    """
    Display reconstruction results in 3 cross-sectional views
    
    Args:
        rec_image: 3D reconstructed volume
        title: Plot title
        save_path: Optional path to save the plot
    """
    X_c = rec_image.shape[0] // 2
    Y_c = rec_image.shape[1] // 2
    Z_c = rec_image.shape[2] // 2
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(title, fontsize=14)
    
    # Cross sections
    im1 = axes[0].imshow(rec_image[30, :, :].T, cmap='hot', origin='lower')
    axes[0].set_title('YZ cross section (X=30)')
    axes[0].axis('off')
    plt.colorbar(im1, ax=axes[0], shrink=0.8)
    
    im2 = axes[1].imshow(rec_image[:, Y_c, :].T, cmap='hot', origin='lower')
    axes[1].set_title(f'XZ cross section (Y={Y_c})')
    axes[1].axis('off')
    plt.colorbar(im2, ax=axes[1], shrink=0.8)
    
    im3 = axes[2].imshow(rec_image[:, :, Z_c].T, cmap='hot', origin='lower')
    axes[2].set_title(f'XY cross section (Z={Z_c})')
    axes[2].axis('off')
    plt.colorbar(im3, ax=axes[2], shrink=0.8)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Reconstruction results saved to: {save_path}")
    
    plt.show()
    
    # Print statistics
    print(f"\nReconstruction Statistics:")
    print(f"  Shape: {rec_image.shape}")
    print(f"  Min: {rec_image.min():.6f}")
    print(f"  Max: {rec_image.max():.6f}")
    print(f"  Mean: {rec_image.mean():.6f}")
    print(f"  Std: {rec_image.std():.6f}")
    print(f"  Non-zero voxels: {np.count_nonzero(rec_image)} / {rec_image.size}")
    print(f"  Coverage: {np.count_nonzero(rec_image)/rec_image.size*100:.1f}%")

def main():
    """
    Main function demonstrating filtered backprojection reconstruction
    """
    print("="*70)
    print("FILTERED BACKPROJECTION EPID DOSE RECONSTRUCTION")
    print("="*70)
    print(f"Started at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Example data path - modify as needed
    data_path = r"E:\CMC\pyprojects\radio_therapy\dose-3d\dataset\VMAT 2025 - 6. SIB COMPLEX TARGET\T1\873251691"
    
    if not os.path.exists(data_path):
        print(f"❌ Data path not found: {data_path}")
        print("Please modify the data_path variable in main() function")
        return
    
    # Load DICOM files
    print("\n📁 STEP 1: Loading DICOM data...")
    print(f"   Scanning directory: {data_path}")
    
    files = [f for f in os.listdir(data_path) if f.endswith('.dcm')]
    files.sort()
    paths = [os.path.join(data_path, f) for f in files]
    
    print(f"   Found {len(files)} DICOM files")
    print(f"   First file: {files[0] if files else 'None'}")
    print(f"   Last file: {files[-1] if files else 'None'}")
    
    # Load images
    print("   Loading DICOM pixel data and angles...")
    start_load = time.time()
    raw_images, raw_angles = load_dicom_fast(paths, max_workers=6)
    load_time = time.time() - start_load
    
    print(f"   ✅ Loaded {len(raw_images)} images in {load_time:.1f}s")
    print(f"   Image shape: {raw_images[0].shape}")
    print(f"   Data type: {raw_images.dtype}")
    print(f"   Angle range: {raw_angles.min():.1f}° to {raw_angles.max():.1f}°")
    
    # Process differential images
    print("\n🔄 STEP 2: Processing differential images...")
    print("   Computing frame-to-frame differences...")
    
    start_diff = time.time()
    processed_images, processed_angles = process_differential_unrotated(raw_images, raw_angles)
    diff_time = time.time() - start_diff
    
    print(f"   ✅ Processed differences in {diff_time:.1f}s")
    
    # Sort by angle
    print("   Sorting projections by gantry angle...")
    sorted_indices = np.argsort(processed_angles)
    sorted_images = processed_images[sorted_indices]
    sorted_angles = processed_angles[sorted_indices]
    
    print(f"   ✅ Sorted {len(sorted_images)} projections")
    print(f"   Final angle range: {sorted_angles.min():.1f}° to {sorted_angles.max():.1f}°")
    print(f"   Angular spacing: mean={np.mean(np.diff(sorted_angles)):.1f}°, std={np.std(np.diff(sorted_angles)):.1f}°")
    
    # Get reconstruction parameters from DICOM
    print("\n⚙️  STEP 3: Setting up reconstruction parameters...")
    print("   Reading DICOM metadata...")
    
    dcm = dcmread(paths[0])
    SID = dcm.RTImageSID
    SAD = dcm.RadiationMachineSAD
    SOD = SAD
    width = 0.172  # mm
    delta_dd = width * SOD / SID
    Nimage = 100
    
    print(f"   Source-to-Axis Distance (SAD): {SOD} mm")
    print(f"   Source-to-Image Distance (SID): {SID} mm")
    print(f"   Effective pixel size (delta_dd): {delta_dd:.6f} mm")
    print(f"   Reconstruction grid size: {Nimage}³ voxels")
    print(f"   Expected output volume: {Nimage*Nimage*int(Nimage*sorted_images.shape[1]/sorted_images.shape[2])} voxels")
    
    # Create reconstructor
    print("   Initializing reconstructor...")
    reconstructor = FilteredBackprojectionReconstructor(SOD, delta_dd, Nimage)
    print("   ✅ Reconstructor ready")
    
    # Test different filter configurations
    print("\n🔧 STEP 4: Running reconstructions with different filters...")
    
    test_configs = [
        {'filter_type': 'shepp_logan', 'filter_size_factor': 0.5, 'name': 'Shepp-Logan (medium)'},
        {'filter_type': 'hamming', 'filter_size_factor': 0.5, 'name': 'Hamming (smooth)'},
        {'filter_type': 'shepp_logan', 'filter_size_factor': 0.25, 'name': 'Shepp-Logan (small)'},
    ]
    
    results = {}
    total_start = time.time()
    
    for i, config in enumerate(test_configs):
        print(f"\n   🎯 Test {i+1}/3: {config['name']}")
        print(f"      Filter: {config['filter_type']}")
        print(f"      Size factor: {config['filter_size_factor']}")
        print("      Starting reconstruction...")
        
        start_time = time.time()
        rec_image = reconstructor.reconstruct_filtered_backprojection(
            sorted_images, sorted_angles,
            filter_type=config['filter_type'],
            filter_size_factor=config['filter_size_factor'],
            use_weighting=True,
            chunk_size=25,
            use_global_mask=True,
            verbose=True  # Enable verbose output for reconstruction details
        )
        
        reconstruction_time = time.time() - start_time
        coverage = np.count_nonzero(rec_image) / rec_image.size * 100
        
        config_key = f"{config['filter_type']}_{config['filter_size_factor']}"
        results[config_key] = {
            'image': rec_image,
            'time': reconstruction_time,
            'coverage': coverage,
            'max_val': rec_image.max(),
            'mean_val': rec_image.mean(),
            'name': config['name']
        }
        
        print(f"      ✅ Completed in {reconstruction_time:.1f}s")
        print(f"      Coverage: {coverage:.1f}% ({np.count_nonzero(rec_image):,} non-zero voxels)")
        print(f"      Value range: {rec_image.min():.6f} to {rec_image.max():.6f}")
        print(f"      Mean value: {rec_image.mean():.6f}")
    
    total_reconstruction_time = time.time() - total_start
    print(f"\n   ✅ All reconstructions completed in {total_reconstruction_time:.1f}s")
    
    # Display results
    print("\n📊 STEP 5: Analyzing and displaying results...")
    
    # Show comparison table
    print("\n   RECONSTRUCTION COMPARISON:")
    print("   " + "="*80)
    print(f"   {'Configuration':<25} {'Time(s)':<8} {'Coverage(%)':<12} {'Max Value':<12} {'Mean Value':<12}")
    print("   " + "-"*80)
    
    for key, result in results.items():
        print(f"   {result['name']:<25} {result['time']:<8.1f} {result['coverage']:<12.1f} {result['max_val']:<12.6f} {result['mean_val']:<12.6f}")
    
    # Find best result
    best_config = max(results.keys(), key=lambda k: results[k]['coverage'])
    best_result = results[best_config]
    
    print(f"\n   🏆 Best configuration: {best_result['name']}")
    print(f"      Coverage: {best_result['coverage']:.1f}%")
    print(f"      Reconstruction time: {best_result['time']:.1f}s")
    print(f"      Max value: {best_result['max_val']:.6f}")
    
    # Show filter visualization
    print("\n   Generating filter visualization...")
    visualize_filters(delta_dd=delta_dd)
    
    # Show best reconstruction
    print("   Displaying reconstruction results...")
    best_image = best_result['image']
    
    display_reconstruction_results(
        best_image, 
        title=f"Best Result: {best_result['name']}"
    )
    
    print(f"\n📋 STEP 6: Final Summary")
    print("   " + "="*50)
    print(f"   Total processing time: {time.time() - total_start + load_time + diff_time:.1f}s")
    print(f"     - Data loading: {load_time:.1f}s")
    print(f"     - Differential processing: {diff_time:.1f}s") 
    print(f"     - Reconstructions: {total_reconstruction_time:.1f}s")
    print(f"   Best filter configuration: {best_result['name']}")
    print(f"   Final reconstruction coverage: {best_result['coverage']:.1f}%")
    print(f"   Completed at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    print(f"\n✨ FILTERED BACKPROJECTION RECONSTRUCTION COMPLETE! ✨")
    print("   You can experiment with different filter parameters:")
    print("   - Adjust filter_size_factor (0.1 to 1.0)")
    print("   - Try different filter_type ('shepp_logan', 'ram_lak', 'cosine', 'hamming')")
    print("   - Modify reconstruction size (Nimage parameter)")
    
    return best_image, results

if __name__ == "__main__":
    # Run the main reconstruction
    reconstructed_volume, all_results = main()