import pydicom
import numpy as np
import os

def analyze_epid_dicom(file_path):
    """
    Analyze EPID DICOM file and extract relevant information for 3D dose reconstruction
    """
    print(f"Analyzing DICOM file: {file_path}")
    print("="*60)
    
    # Read the DICOM file
    ds = pydicom.dcmread(file_path)
    
    # 1. DICOM Headers - Basic Information
    print("1. BASIC DICOM INFORMATION:")
    print(f"   Modality: {getattr(ds, 'Modality', 'Not found')}")
    print(f"   SOP Class UID: {getattr(ds, 'SOPClassUID', 'Not found')}")
    print(f"   Study Date: {getattr(ds, 'StudyDate', 'Not found')}")
    print(f"   Study Time: {getattr(ds, 'StudyTime', 'Not found')}")
    print(f"   Manufacturer: {getattr(ds, 'Manufacturer', 'Not found')}")
    print(f"   Station Name: {getattr(ds, 'StationName', 'Not found')}")
    print()
    
    # 2. Gantry Angle and Treatment Parameters
    print("2. GANTRY AND TREATMENT PARAMETERS:")
    print(f"   Gantry Angle: {getattr(ds, 'GantryAngle', 'Not found')}")
    print(f"   Patient Position: {getattr(ds, 'PatientPosition', 'Not found')}")
    print(f"   Beam Limiting Device Angle: {getattr(ds, 'BeamLimitingDeviceAngle', 'Not found')}")
    print(f"   Patient Support Angle: {getattr(ds, 'PatientSupportAngle', 'Not found')}")
    print(f"   Table Top Eccentric Angle: {getattr(ds, 'TableTopEccentricAngle', 'Not found')}")
    print(f"   Table Top Pitch Angle: {getattr(ds, 'TableTopPitchAngle', 'Not found')}")
    print(f"   Table Top Roll Angle: {getattr(ds, 'TableTopRollAngle', 'Not found')}")
    print()
    
    # 3. Image Dimensions and Pixel Information
    print("3. IMAGE DIMENSIONS AND PIXEL INFORMATION:")
    print(f"   Rows: {getattr(ds, 'Rows', 'Not found')}")
    print(f"   Columns: {getattr(ds, 'Columns', 'Not found')}")
    print(f"   Bits Allocated: {getattr(ds, 'BitsAllocated', 'Not found')}")
    print(f"   Bits Stored: {getattr(ds, 'BitsStored', 'Not found')}")
    print(f"   High Bit: {getattr(ds, 'HighBit', 'Not found')}")
    print(f"   Pixel Representation: {getattr(ds, 'PixelRepresentation', 'Not found')}")
    print(f"   Samples Per Pixel: {getattr(ds, 'SamplesPerPixel', 'Not found')}")
    print(f"   Photometric Interpretation: {getattr(ds, 'PhotometricInterpretation', 'Not found')}")
    print()
    
    # 4. Pixel Spacing and Geometry
    print("4. PIXEL SPACING AND GEOMETRY:")
    if hasattr(ds, 'PixelSpacing'):
        print(f"   Pixel Spacing: {ds.PixelSpacing} mm")
    else:
        print("   Pixel Spacing: Not found")
    
    if hasattr(ds, 'ImagerPixelSpacing'):
        print(f"   Imager Pixel Spacing: {ds.ImagerPixelSpacing} mm")
    else:
        print("   Imager Pixel Spacing: Not found")
    
    print(f"   Slice Thickness: {getattr(ds, 'SliceThickness', 'Not found')}")
    print(f"   Slice Location: {getattr(ds, 'SliceLocation', 'Not found')}")
    
    # Image Position and Orientation
    if hasattr(ds, 'ImagePositionPatient'):
        print(f"   Image Position Patient: {ds.ImagePositionPatient}")
    else:
        print("   Image Position Patient: Not found")
    
    if hasattr(ds, 'ImageOrientationPatient'):
        print(f"   Image Orientation Patient: {ds.ImageOrientationPatient}")
    else:
        print("   Image Orientation Patient: Not found")
    print()
    
    # 5. RT-specific Information
    print("5. RT-SPECIFIC INFORMATION:")
    print(f"   RT Image Label: {getattr(ds, 'RTImageLabel', 'Not found')}")
    print(f"   RT Image Name: {getattr(ds, 'RTImageName', 'Not found')}")
    print(f"   RT Image Description: {getattr(ds, 'RTImageDescription', 'Not found')}")
    print(f"   RT Image Plane: {getattr(ds, 'RTImagePlane', 'Not found')}")
    print(f"   X-Ray Image Receptor Translation: {getattr(ds, 'XRayImageReceptorTranslation', 'Not found')}")
    print(f"   X-Ray Image Receptor Angle: {getattr(ds, 'XRayImageReceptorAngle', 'Not found')}")
    print(f"   RT Image Position: {getattr(ds, 'RTImagePosition', 'Not found')}")
    print(f"   Radiation Machine Name: {getattr(ds, 'RadiationMachineName', 'Not found')}")
    print(f"   Radiation Machine SAD: {getattr(ds, 'RadiationMachineSAD', 'Not found')}")
    print(f"   RT Image SID: {getattr(ds, 'RTImageSID', 'Not found')}")
    print()
    
    # 6. Dose-related Information
    print("6. DOSE-RELATED INFORMATION:")
    print(f"   Primary Dosimeter Unit: {getattr(ds, 'PrimaryDosimeterUnit', 'Not found')}")
    print(f"   Exposure: {getattr(ds, 'Exposure', 'Not found')}")
    print(f"   Exposure Time: {getattr(ds, 'ExposureTime', 'Not found')}")
    print(f"   Rescale Intercept: {getattr(ds, 'RescaleIntercept', 'Not found')}")
    print(f"   Rescale Slope: {getattr(ds, 'RescaleSlope', 'Not found')}")
    print(f"   Rescale Type: {getattr(ds, 'RescaleType', 'Not found')}")
    print()
    
    # 7. Pixel Data Analysis
    print("7. PIXEL DATA ANALYSIS:")
    if hasattr(ds, 'pixel_array'):
        pixel_array = ds.pixel_array
        print(f"   Pixel Array Shape: {pixel_array.shape}")
        print(f"   Data Type: {pixel_array.dtype}")
        print(f"   Min Value: {np.min(pixel_array)}")
        print(f"   Max Value: {np.max(pixel_array)}")
        print(f"   Mean Value: {np.mean(pixel_array):.2f}")
        print(f"   Std Dev: {np.std(pixel_array):.2f}")
        
        # If rescale slope/intercept exist, show calibrated values
        if hasattr(ds, 'RescaleSlope') and hasattr(ds, 'RescaleIntercept'):
            calibrated = pixel_array * ds.RescaleSlope + ds.RescaleIntercept
            print(f"   Calibrated Min: {np.min(calibrated):.2f}")
            print(f"   Calibrated Max: {np.max(calibrated):.2f}")
            print(f"   Calibrated Mean: {np.mean(calibrated):.2f}")
    else:
        print("   Pixel Array: Not accessible")
    print()
    
    # 8. Beam and Collimator Information
    print("8. BEAM AND COLLIMATOR INFORMATION:")
    if hasattr(ds, 'BeamLimitingDeviceSequence'):
        for i, device in enumerate(ds.BeamLimitingDeviceSequence):
            print(f"   Device {i+1}:")
            print(f"     RT Beam Limiting Device Type: {getattr(device, 'RTBeamLimitingDeviceType', 'Not found')}")
            if hasattr(device, 'LeafJawPositions'):
                print(f"     Leaf/Jaw Positions: {device.LeafJawPositions}")
    else:
        print("   Beam Limiting Device Sequence: Not found")
    print()
    
    # 9. Additional RT Image Information
    print("9. ADDITIONAL RT IMAGE INFORMATION:")
    print(f"   Fluence Map Sequence: {hasattr(ds, 'FluenceMapSequence')}")
    print(f"   Fluence Data Source: {getattr(ds, 'FluenceDataSource', 'Not found')}")
    print(f"   Fluence Data Scale: {getattr(ds, 'FluenceDataScale', 'Not found')}")
    
    # Check for any private tags that might contain dose information
    print("\n10. PRIVATE TAGS (first 10):")
    private_tags = [elem for elem in ds if elem.tag.is_private]
    for i, elem in enumerate(private_tags[:10]):
        print(f"   {elem.tag}: {elem.keyword} = {elem.value}")
    
    if len(private_tags) > 10:
        print(f"   ... and {len(private_tags) - 10} more private tags")
    
    print("\n" + "="*60)
    print("Analysis complete!")

if __name__ == "__main__":
    file_path = r"E:\CMC\pyprojects\radio_therapy\dose-3d\dataset\10x10_TRIAL 1\873251696\00000.dcm"
    
    if os.path.exists(file_path):
        analyze_epid_dicom(file_path)
    else:
        print(f"Error: File not found at {file_path}")