#!/usr/bin/env python3
"""
Build script for EPID Reconstruction GUI using PyInstaller
Usage: python build_gui.py
"""

import os
import subprocess
import shutil
from pathlib import Path

def build_epid_gui():
    """Build the EPID GUI executable using PyInstaller"""
    
    print("Building EPID Reconstruction GUI...")
    
    # Check required files exist
    required_files = ["epid_reconstruction_gui.py", "tps_reference_data.py", "EPID_12_t0.dcm"]
    for file in required_files:
        if not Path(file).exists():
            print(f"Error: Required file '{file}' not found!")
            return False
    
    # Clean previous build
    if Path("dist").exists():
        print("Cleaning previous build...")
        shutil.rmtree("dist")
    
    # PyInstaller command
    cmd = [
        "pyinstaller",
        "--windowed",                           # No console window
        "--name=EPID_Reconstruction_GUI",       # Output name
        "--add-data", "EPID_12_t0.dcm;.",      # Include DICOM template
        "--add-data", "tps_reference_data.py;.", # Include TPS module
        "epid_reconstruction_gui.py"           # Main script
    ]
    
    print("Running PyInstaller...")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Build completed successfully!")
        
        # Show build info
        dist_path = Path("dist/EPID_Reconstruction_GUI")
        if dist_path.exists():
            size_mb = sum(f.stat().st_size for f in dist_path.rglob('*') if f.is_file()) / (1024*1024)
            print(f"Output directory: {dist_path}")
            print(f"Total size: {size_mb:.0f} MB")
            print(f"Executable: {dist_path}/EPID_Reconstruction_GUI.exe")
            print(f"DICOM template included: EPID_12_t0.dcm")
            print(f"TPS reference included: tps_reference_data.py")
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"Build failed: {e}")
        if e.stderr:
            print(f"Error output: {e.stderr}")
        return False

if __name__ == "__main__":
    success = build_epid_gui()
    if success:
        print("\nBuild completed! You can now distribute the 'dist/EPID_Reconstruction_GUI' folder.")
    else:
        print("\nBuild failed. Please check the errors above.")