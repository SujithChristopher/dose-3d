"""
Build script for creating standalone executable
Handles the build process and provides user feedback
"""

import subprocess
import sys
import os
import shutil
from pathlib import Path

def check_requirements():
    """Check if cx_Freeze is installed"""
    try:
        import cx_Freeze
        print(f"PASS - cx_Freeze version {cx_Freeze.__version__} found")
        return True
    except ImportError:
        print("FAIL - cx_Freeze not found. Install with: pip install cx_Freeze")
        return False

def clean_build_directory():
    """Clean previous build artifacts"""
    build_dir = Path("build")
    if build_dir.exists():
        print("Cleaning previous build...")
        shutil.rmtree(build_dir)
        print("PASS - Build directory cleaned")

def build_executable():
    """Build the standalone executable"""
    print("Building standalone executable...")
    print("This may take several minutes...")
    
    try:
        # Run cx_Freeze build
        result = subprocess.run([
            sys.executable, "setup_freeze.py", "build"
        ], capture_output=True, text=True, cwd=os.getcwd())
        
        if result.returncode == 0:
            print("PASS - Build completed successfully!")
            
            # Check if executable was created
            exe_path = Path("build/epid_reconstruction_gui/EPID_Reconstruction_GUI.exe")
            if exe_path.exists():
                print(f"PASS - Executable created: {exe_path}")
                print(f"PASS - Executable size: {exe_path.stat().st_size / (1024*1024):.1f} MB")
                return True
            else:
                print("FAIL - Executable not found after build")
                return False
                
        else:
            print("FAIL - Build failed!")
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            return False
            
    except Exception as e:
        print(f"FAIL - Build error: {e}")
        return False

def create_installer_info():
    """Create information file for distribution"""
    info_content = """
EPID Reconstruction GUI - Standalone Application
================================================

Installation:
1. Extract all files to a folder (e.g., C:\\EPID_Reconstruction)
2. Run EPID_Reconstruction_GUI.exe
3. No Python installation required!

System Requirements:
- Windows 10/11 (64-bit)
- Minimum 4GB RAM
- 2GB free disk space

Usage:
1. Select EPID data folder containing DICOM files
2. Configure reconstruction parameters
3. Start reconstruction and view results
4. Export results as DICOM with TPS headers

For support or issues, contact the development team.

Built with cx_Freeze and Python {python_version}
""".format(python_version=sys.version.split()[0])
    
    with open("build/epid_reconstruction_gui/README.txt", "w") as f:
        f.write(info_content)
    
    print("PASS - Installation info created")

def main():
    """Main build process"""
    print("EPID Reconstruction GUI - Standalone Build")
    print("=" * 50)
    
    # Check requirements
    if not check_requirements():
        return 1
    
    # Clean previous build
    clean_build_directory()
    
    # Build executable
    if not build_executable():
        return 1
    
    # Create installer info
    create_installer_info()
    
    print("\n" + "=" * 50)
    print("Build completed successfully!")
    print("\nStandalone application location:")
    print("  build/epid_reconstruction_gui/")
    print("\nTo distribute:")
    print("  1. Zip the entire 'epid_reconstruction_gui' folder")
    print("  2. Users extract and run EPID_Reconstruction_GUI.exe")
    print("  3. No Python installation required on target machines")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())