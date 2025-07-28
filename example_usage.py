#!/usr/bin/env python3
"""
Example Usage of Filtered Backprojection Reconstruction

This script demonstrates how to use the filtered backprojection
reconstruction program with your own data.

Author: Generated with Claude Code
"""

from filtered_backprojection_reconstruction import *
import os

def simple_reconstruction_example():
    """Simple example of running a filtered backprojection reconstruction"""
    
    print("Simple Filtered Backprojection Example")
    print("="*50)
    
    # 1. Set your data path
    data_path = r"E:\CMC\pyprojects\radio_therapy\dose-3d\dataset\VMAT 2025 - 6. SIB COMPLEX TARGET\T1\873251691"
    
    # Check if path exists
    if not os.path.exists(data_path):
        print(f"⚠️  Data path not found: {data_path}")
        print("Please update the data_path variable with your DICOM folder")
        return None
    
    # 2. Load your DICOM data
    print("Loading DICOM files...")
    files = [f for f in os.listdir(data_path) if f.endswith('.dcm')]
    files.sort()
    paths = [os.path.join(data_path, f) for f in files]
    
    raw_images, raw_angles = load_dicom_fast(paths[:50])  # Use first 50 for speed
    
    # 3. Process differential images
    print("Processing differential images...")
    processed_images, processed_angles = process_differential_unrotated(raw_images, raw_angles)
    
    # Sort by angle
    sorted_indices = np.argsort(processed_angles)
    sorted_images = processed_images[sorted_indices]
    sorted_angles = processed_angles[sorted_indices]
    
    # 4. Get reconstruction parameters
    dcm = dcmread(paths[0])
    SOD = dcm.RadiationMachineSAD
    SID = dcm.RTImageSID
    delta_dd = 0.172 * SOD / SID
    
    # 5. Create reconstructor and run reconstruction
    print("Running filtered backprojection...")
    reconstructor = FilteredBackprojectionReconstructor(SOD, delta_dd, Nimage=50)  # Smaller for speed
    
    result = reconstructor.reconstruct_filtered_backprojection(
        sorted_images, sorted_angles,
        filter_type='shepp_logan',      # Try: 'shepp_logan', 'ram_lak', 'cosine', 'hamming'
        filter_size_factor=0.5,         # Try: 0.1 to 1.0
        use_weighting=True,
        chunk_size=10,
        use_global_mask=True
    )
    
    # 6. Display results
    display_reconstruction_results(result, "Simple Example Results")
    
    return result

def compare_filters_example():
    """Example comparing different filter types"""
    
    print("Filter Comparison Example")
    print("="*30)
    
    # Just visualize filters without running full reconstruction
    print("Visualizing different filter types...")
    visualize_filters(delta_dd=0.1, N=512)
    
    print("\nFilter Recommendations:")
    print("• shepp_logan (0.5): Good general purpose")
    print("• hamming (0.3): Smoothest, best for noisy data")  
    print("• ram_lak (0.5): Sharpest, best for high detail")
    print("• cosine (0.4): Good compromise between sharp and smooth")

def custom_reconstruction_example():
    """Example showing how to customize reconstruction parameters"""
    
    print("Custom Reconstruction Example")  
    print("="*35)
    
    # Example parameters - modify these for your data
    data_path = "your/dicom/folder/path"
    
    reconstruction_params = {
        'filter_type': 'hamming',        # Filter type
        'filter_size_factor': 0.3,       # Filter size (smaller = smoother)
        'use_weighting': True,           # Apply geometric weighting
        'chunk_size': 25,                # Memory management
        'use_global_mask': True,         # Consistent masking
        'Nimage': 100                    # Reconstruction grid size
    }
    
    print("Reconstruction parameters:")
    for key, value in reconstruction_params.items():
        print(f"  {key}: {value}")
    
    print("\nTo run with your data:")
    print("1. Update 'data_path' variable")
    print("2. Adjust reconstruction_params as needed")
    print("3. Load data and run reconstruction")
    
    # Example code structure:
    example_code = '''
# Load your data
files = [f for f in os.listdir(data_path) if f.endswith('.dcm')]
paths = [os.path.join(data_path, f) for f in files]
raw_images, raw_angles = load_dicom_fast(paths)

# Process data  
processed_images, processed_angles = process_differential_unrotated(raw_images, raw_angles)
sorted_indices = np.argsort(processed_angles)
sorted_images = processed_images[sorted_indices]
sorted_angles = processed_angles[sorted_indices]

# Get parameters from DICOM
dcm = dcmread(paths[0])
SOD = dcm.RadiationMachineSAD
SID = dcm.RTImageSID
delta_dd = 0.172 * SOD / SID

# Run reconstruction
reconstructor = FilteredBackprojectionReconstructor(SOD, delta_dd, reconstruction_params['Nimage'])
result = reconstructor.reconstruct_filtered_backprojection(
    sorted_images, sorted_angles, **reconstruction_params
)

# Display results
display_reconstruction_results(result)
'''
    
    print("Example code:")
    print(example_code)

if __name__ == "__main__":
    print("Filtered Backprojection Examples")
    print("="*40)
    print()
    print("Choose an example to run:")
    print("1. Simple reconstruction (if data path is correct)")
    print("2. Compare different filters (visualization only)")
    print("3. Show custom reconstruction template")
    print()
    
    choice = input("Enter choice (1, 2, or 3): ").strip()
    
    if choice == "1":
        simple_reconstruction_example()
    elif choice == "2":
        compare_filters_example()
    elif choice == "3":
        custom_reconstruction_example()
    else:
        print("Invalid choice. Running filter comparison...")
        compare_filters_example()