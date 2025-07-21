"""
TPS Reference Data
Hardcoded reference dose distribution from EPID_12_t0.dcm for dose verification.
This avoids repeatedly reading the large DICOM file.
"""

import numpy as np

# TPS Reference metadata
TPS_METADATA = {
    'shape': (100, 100, 100),
    'data_type': 'uint32',
    'min_value': 5536,
    'max_value': 1038022,
    'source_file': 'EPID_12_t0.dcm',
    'description': 'Reference dose distribution for EPID reconstruction verification'
}

def load_tps_dicom_template():
    """
    Load the TPS DICOM file to use as template for exporting reconstructed data.
    The file should be in the same directory as this Python script.
    
    Returns:
        pydicom.Dataset: TPS DICOM dataset for header copying, or None if unavailable
    """
    try:
        from pydicom import dcmread
        import os
        from pathlib import Path
        
        # Get the directory where this Python file is located
        script_dir = Path(__file__).parent
        
        # TPS reference file should be in the same directory
        tps_path = script_dir / "EPID_12_t0.dcm"
        
        if tps_path.exists():
            dcm = dcmread(str(tps_path))
            return dcm
        else:
            print(f"TPS template file not found: {tps_path}")
            print("Please place EPID_12_t0.dcm in the same directory as the application")
            return None
            
    except Exception as e:
        print(f"Error loading TPS template: {e}")
        return None

def export_reconstructed_dicom(reconstructed_volume, output_path, scale_factor=1.0):
    """
    Export reconstructed volume as DICOM using TPS template headers.
    
    Args:
        reconstructed_volume (numpy.ndarray): Reconstructed dose volume
        output_path (str): Path to save the exported DICOM file
        scale_factor (float): Scale factor for intensity values
        
    Returns:
        bool: True if export successful, False otherwise
    """
    try:
        from pydicom import dcmwrite
        import datetime
        
        # Load TPS template
        tps_template = load_tps_dicom_template()
        if tps_template is None:
            return False
        
        # Create copy of template
        export_dcm = tps_template.copy()
        
        # Scale and convert reconstructed data
        scaled_volume = (reconstructed_volume * scale_factor).astype(np.uint32)
        
        # Update DICOM fields with new data
        export_dcm.PixelData = scaled_volume.tobytes()
        export_dcm.Rows = reconstructed_volume.shape[0]
        export_dcm.Columns = reconstructed_volume.shape[1]
        export_dcm.NumberOfFrames = str(reconstructed_volume.shape[2])
        
        # Update metadata
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        export_dcm.SeriesDescription = f"EPID_Reconstruction_{timestamp}"
        export_dcm.StudyDescription = "EPID Dose Reconstruction"
        
        # Save DICOM file
        dcmwrite(output_path, export_dcm)
        return True
        
    except Exception as e:
        print(f"Export error: {e}")
        return False

def calculate_volume_statistics(volume):
    """
    Calculate comprehensive statistics for the reconstructed volume.
    
    Args:
        volume (numpy.ndarray): 3D dose volume
        
    Returns:
        dict: Dictionary containing volume statistics
    """
    if volume is None:
        return {}
    
    # Basic statistics
    stats = {
        'shape': volume.shape,
        'total_voxels': volume.size,
        'min_value': float(volume.min()),
        'max_value': float(volume.max()),
        'mean_value': float(volume.mean()),
        'std_value': float(volume.std()),
        'median_value': float(np.median(volume)),
        'nonzero_voxels': int(np.count_nonzero(volume)),
        'zero_voxels': int(volume.size - np.count_nonzero(volume)),
        'coverage_percent': float(100 * np.count_nonzero(volume) / volume.size)
    }
    
    # Percentiles
    stats['percentile_25'] = float(np.percentile(volume, 25))
    stats['percentile_75'] = float(np.percentile(volume, 75))
    stats['percentile_90'] = float(np.percentile(volume, 90))
    stats['percentile_95'] = float(np.percentile(volume, 95))
    
    # Additional metrics
    stats['range'] = stats['max_value'] - stats['min_value']
    stats['coefficient_of_variation'] = stats['std_value'] / stats['mean_value'] if stats['mean_value'] > 0 else 0
    
    return stats