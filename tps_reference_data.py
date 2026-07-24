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

# Stored-integer full scale for the exported dose grid. Eclipse uses ~1e6 for the
# same purpose: the absolute dose lives entirely in DoseGridScaling, and the
# integers only need enough range to carry the relative distribution.
_STORED_MAX = 1_000_000


def export_reconstructed_dicom(reconstructed_volume, output_path, scale_factor=1.0,
                               voxel_size_mm=None, origin_mm=None,
                               dose_units='GY', frame_of_reference_uid=None,
                               patient_id=None, patient_name=None,
                               series_description=None):
    """
    Export a reconstructed volume as an RTDOSE using the TPS file for boilerplate only.

    The template supplies patient/study identification and the general RTDOSE
    skeleton; everything that a gamma-analysis tool (VeriSoft et al.) actually
    reads is recomputed from this volume, because inheriting it produces a file
    that describes the template's dose grid rather than this one.

    Args:
        reconstructed_volume (numpy.ndarray): 3D volume indexed [frame, row, column]
            i.e. (z, y, x), matching the DICOM multi-frame byte order.
        output_path (str): Path to save the exported DICOM file
        scale_factor (float): Multiplier applied to the volume before encoding.
            Use it to convert into the unit named by `dose_units` (e.g. 0.01 for
            a cGy volume exported as GY).
        voxel_size_mm (tuple): (dz, dy, dx) grid spacing. Defaults to the template's.
        origin_mm (tuple): (z0, y0, x0) position of voxel [0,0,0]. Defaults to the
            template's ImagePositionPatient.
        dose_units (str): 'GY' if the scaled volume is absolute dose in gray,
            otherwise 'RELATIVE'.
        frame_of_reference_uid (str): Set this to the reference dose's
            FrameOfReferenceUID so the two datasets register against each other.
        patient_id (str), patient_name (str): Override the template's identifiers.
        series_description (str): Defaults to a timestamped EPID reconstruction tag.

    Returns:
        bool: True if export successful, False otherwise
    """
    try:
        from pydicom import dcmwrite
        from pydicom.uid import generate_uid
        import datetime

        tps_template = load_tps_dicom_template()
        if tps_template is None:
            return False

        volume = np.asarray(reconstructed_volume, dtype=np.float64)
        if volume.ndim != 3:
            print(f"Export error: expected a 3D volume, got shape {volume.shape}")
            return False

        export_dcm = tps_template.copy()

        # --- pixel data -------------------------------------------------------
        # RTDOSE stored values are unsigned, but filtered backprojection produces
        # negative lobes; casting those straight to uint32 wraps them to ~4e9.
        scaled = volume * float(scale_factor)
        n_negative = int((scaled < 0).sum())
        if n_negative:
            print(f"Export note: clipped {n_negative} negative voxels "
                  f"({100.0 * n_negative / scaled.size:.2f}%) to zero")
        scaled = np.clip(scaled, 0.0, None)

        peak = float(scaled.max())
        dose_grid_scaling = peak / _STORED_MAX if peak > 0 else 1.0
        stored = np.rint(scaled / dose_grid_scaling).astype(np.uint32)

        n_frames, n_rows, n_cols = stored.shape
        export_dcm.PixelData = stored.tobytes()
        export_dcm.NumberOfFrames = str(n_frames)
        export_dcm.Rows = n_rows
        export_dcm.Columns = n_cols
        export_dcm.BitsAllocated = 32
        export_dcm.BitsStored = 32
        export_dcm.HighBit = 31
        export_dcm.PixelRepresentation = 0
        export_dcm.SamplesPerPixel = 1
        export_dcm.PhotometricInterpretation = 'MONOCHROME2'

        # --- dose scaling: the only tags that turn integers back into dose ----
        export_dcm.DoseGridScaling = f"{dose_grid_scaling:.10g}"
        export_dcm.DoseUnits = dose_units
        export_dcm.DoseType = 'PHYSICAL'
        export_dcm.DoseSummationType = 'PLAN'

        # --- geometry ---------------------------------------------------------
        if voxel_size_mm is not None:
            dz, dy, dx = (float(v) for v in voxel_size_mm)
        else:
            dy, dx = (float(v) for v in export_dcm.PixelSpacing)
            offsets = export_dcm.GridFrameOffsetVector
            dz = float(offsets[1] - offsets[0]) if len(offsets) > 1 else dy

        export_dcm.PixelSpacing = [dy, dx]
        export_dcm.SliceThickness = dz
        # Must be exactly NumberOfFrames long - the template's 120-entry vector on
        # a 100-frame volume is what makes importers reject or mis-stack the grid.
        export_dcm.GridFrameOffsetVector = [i * dz for i in range(n_frames)]
        export_dcm.FrameIncrementPointer = 0x3004000C

        if origin_mm is not None:
            z0, y0, x0 = (float(v) for v in origin_mm)
            export_dcm.ImagePositionPatient = [x0, y0, z0]
        export_dcm.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]

        if frame_of_reference_uid:
            export_dcm.FrameOfReferenceUID = frame_of_reference_uid

        # --- identity ---------------------------------------------------------
        # Reusing the template's UIDs made every export collide as one instance.
        sop_uid = generate_uid()
        export_dcm.SOPInstanceUID = sop_uid
        export_dcm.SeriesInstanceUID = generate_uid()
        export_dcm.file_meta.MediaStorageSOPInstanceUID = sop_uid
        export_dcm.InstanceNumber = 1

        if patient_id is not None:
            export_dcm.PatientID = patient_id
        if patient_name is not None:
            export_dcm.PatientName = patient_name

        now = datetime.datetime.now()
        timestamp = now.strftime("%Y%m%d_%H%M%S")
        export_dcm.SeriesDescription = series_description or f"EPID_Reconstruction_{timestamp}"
        export_dcm.StudyDescription = "EPID Dose Reconstruction"
        export_dcm.ContentDate = now.strftime("%Y%m%d")
        export_dcm.ContentTime = now.strftime("%H%M%S")
        export_dcm.SeriesDate = export_dcm.ContentDate
        export_dcm.SeriesTime = export_dcm.ContentTime

        # Inherited from the template and describing a different dose distribution.
        for tag in ('DVHSequence', 'ReferencedRTPlanSequence',
                    'ReferencedStructureSetSequence', 'TissueHeterogeneityCorrection'):
            if tag in export_dcm:
                delattr(export_dcm, tag)

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