"""
Create a portable application package without compilation
This approach creates a portable package with Python and dependencies
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path
import zipfile

def create_portable_package():
    """Create a portable application package"""
    print("Creating portable EPID Reconstruction GUI package...")
    
    # Create package directory
    package_dir = Path("portable_epid_gui")
    if package_dir.exists():
        shutil.rmtree(package_dir)
    package_dir.mkdir()
    
    print("Setting up portable package structure...")
    
    # Copy main application files
    app_files = [
        "epid_reconstruction_gui.py",
        "tps_reference_data.py",
        "test_gui_components.py",
        "GUI_README.md",
        "EPID_12_t0.dcm"  # Include TPS DICOM template
    ]
    
    for file in app_files:
        if Path(file).exists():
            shutil.copy2(file, package_dir)
            print(f"  Copied: {file}")
    
    # Create launcher script
    launcher_content = '''@echo off
echo Starting EPID Reconstruction GUI...
echo.
echo Checking Python environment...

python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found!
    echo Please install Python 3.8+ and ensure it's in your PATH
    echo.
    pause
    exit /b 1
)

echo Python found. Checking dependencies...

python -c "import numpy, scipy, matplotlib, pydicom, PySide6, cv2, tqdm" >nul 2>&1
if errorlevel 1 (
    echo ERROR: Required dependencies not found!
    echo Please install the required packages:
    echo   pip install numpy scipy matplotlib pydicom PySide6 opencv-python tqdm numba
    echo.
    pause
    exit /b 1
)

echo All dependencies found. Starting application...
echo.
python epid_reconstruction_gui.py

if errorlevel 1 (
    echo.
    echo Application exited with error. Press any key to close...
    pause >nul
)
'''
    
    with open(package_dir / "run_epid_gui.bat", "w") as f:
        f.write(launcher_content)
    
    # Create Linux launcher
    linux_launcher = '''#!/bin/bash
echo "Starting EPID Reconstruction GUI..."
echo

echo "Checking Python environment..."
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Python3 not found!"
    echo "Please install Python 3.8+ and ensure it's in your PATH"
    read -p "Press Enter to continue..."
    exit 1
fi

echo "Python found. Checking dependencies..."
if ! python3 -c "import numpy, scipy, matplotlib, pydicom, PySide6, cv2, tqdm" &> /dev/null; then
    echo "ERROR: Required dependencies not found!"
    echo "Please install the required packages:"
    echo "  pip install numpy scipy matplotlib pydicom PySide6 opencv-python tqdm numba"
    echo
    read -p "Press Enter to continue..."
    exit 1
fi

echo "All dependencies found. Starting application..."
echo
python3 epid_reconstruction_gui.py
'''
    
    with open(package_dir / "run_epid_gui.sh", "w") as f:
        f.write(linux_launcher)
    
    # Make Linux script executable
    os.chmod(package_dir / "run_epid_gui.sh", 0o755)
    
    # Create requirements file
    requirements = '''# EPID Reconstruction GUI Requirements
numpy>=1.21.0
scipy>=1.7.0
matplotlib>=3.5.0
pydicom>=2.0.0
PySide6>=6.0.0
opencv-python>=4.5.0
tqdm>=4.60.0
numba>=0.55.0
cx_Freeze>=6.0.0
'''
    
    with open(package_dir / "requirements.txt", "w") as f:
        f.write(requirements)
    
    # Create installation instructions
    install_instructions = '''EPID Reconstruction GUI - Portable Package
==========================================

This is a portable version of the EPID Reconstruction GUI that requires
Python to be installed on the target system.

## Quick Start

### Windows:
1. Double-click "run_epid_gui.bat"
2. Follow any installation prompts if dependencies are missing

### Linux/Mac:
1. Open terminal in this folder
2. Run: ./run_epid_gui.sh
3. Follow any installation prompts if dependencies are missing

## Manual Installation

If the launcher scripts don't work, install dependencies manually:

1. Ensure Python 3.8+ is installed
2. Install required packages:
   pip install -r requirements.txt
3. Run the application:
   python epid_reconstruction_gui.py

## System Requirements
- Python 3.8 or higher
- Windows 10/11, Linux, or macOS
- Minimum 4GB RAM
- 2GB free disk space

## Usage
1. Select EPID data folder containing DICOM files
2. Configure reconstruction parameters
3. Start reconstruction and view results
4. Export results as DICOM with TPS headers

## Files Included
- epid_reconstruction_gui.py: Main application
- tps_reference_data.py: TPS utilities and export functions
- test_gui_components.py: Component testing script
- GUI_README.md: Detailed documentation
- requirements.txt: Python package requirements
- run_epid_gui.bat: Windows launcher
- run_epid_gui.sh: Linux/Mac launcher

For technical support, refer to GUI_README.md
'''
    
    with open(package_dir / "INSTALL.txt", "w") as f:
        f.write(install_instructions)
    
    print(f"\nPASS - Portable package created: {package_dir}")
    print(f"Package size: {get_dir_size(package_dir):.1f} MB")
    
    return package_dir

def get_dir_size(path):
    """Calculate directory size in MB"""
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            total_size += os.path.getsize(file_path)
    return total_size / (1024 * 1024)

def create_zip_distribution(package_dir):
    """Create ZIP file for distribution"""
    zip_name = "EPID_Reconstruction_GUI_Portable.zip"
    
    print(f"Creating ZIP distribution: {zip_name}")
    
    with zipfile.ZipFile(zip_name, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(package_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, package_dir.parent)
                zipf.write(file_path, arcname)
    
    zip_size = os.path.getsize(zip_name) / (1024 * 1024)
    print(f"PASS - ZIP created: {zip_name} ({zip_size:.1f} MB)")
    
    return zip_name

def main():
    """Main packaging process"""
    print("EPID Reconstruction GUI - Portable Package Creator")
    print("=" * 60)
    
    # Create portable package
    package_dir = create_portable_package()
    
    # Create ZIP distribution
    zip_file = create_zip_distribution(package_dir)
    
    print("\n" + "=" * 60)
    print("Portable package creation completed!")
    print("\nDistribution options:")
    print(f"1. Folder: {package_dir}")
    print(f"2. ZIP file: {zip_file}")
    print("\nUsers need Python 3.8+ and will auto-install dependencies")
    print("This approach avoids cx_Freeze compilation issues")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())