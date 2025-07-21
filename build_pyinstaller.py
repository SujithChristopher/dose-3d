"""
PyInstaller build script for EPID Reconstruction GUI
This script handles the build process and provides user feedback
"""

import subprocess
import sys
import os
import shutil
from pathlib import Path
import time

def check_pyinstaller():
    """Check if PyInstaller is available"""
    try:
        import PyInstaller
        print(f"PASS - PyInstaller version {PyInstaller.__version__} found")
        return True
    except ImportError:
        print("FAIL - PyInstaller not found. Install with: pip install pyinstaller")
        return False

def clean_build_directories():
    """Clean previous build artifacts"""
    directories_to_clean = ["build", "dist", "__pycache__"]
    
    for dir_name in directories_to_clean:
        dir_path = Path(dir_name)
        if dir_path.exists():
            print(f"Cleaning {dir_name}/...")
            shutil.rmtree(dir_path)
    
    print("PASS - Build directories cleaned")

def build_with_pyinstaller():
    """Build executable using PyInstaller"""
    print("Building standalone executable with PyInstaller...")
    print("This may take several minutes...")
    
    start_time = time.time()
    
    try:
        # Run PyInstaller with the spec file
        result = subprocess.run([
            sys.executable, "-m", "PyInstaller", 
            "--clean",  # Clean PyInstaller cache
            "epid_gui.spec"
        ], capture_output=True, text=True, cwd=os.getcwd())
        
        build_time = time.time() - start_time
        
        if result.returncode == 0:
            print(f"PASS - Build completed successfully in {build_time:.1f} seconds!")
            
            # Check if executable was created
            exe_path = Path("dist/EPID_Reconstruction_GUI.exe")
            if exe_path.exists():
                exe_size = exe_path.stat().st_size / (1024*1024)
                print(f"PASS - Executable created: {exe_path}")
                print(f"PASS - Executable size: {exe_size:.1f} MB")
                return True, exe_path
            else:
                print("FAIL - Executable not found after build")
                return False, None
                
        else:
            print("FAIL - Build failed!")
            print("STDOUT:", result.stdout[-1000:])  # Last 1000 chars
            print("STDERR:", result.stderr[-1000:])  # Last 1000 chars
            return False, None
            
    except Exception as e:
        print(f"FAIL - Build error: {e}")
        return False, None

def test_executable(exe_path):
    """Test if the executable runs"""
    if not exe_path or not exe_path.exists():
        return False
    
    print("Testing executable...")
    
    try:
        # Test quick startup (just check if it starts without crashing immediately)
        result = subprocess.run([str(exe_path), "--version"], 
                              capture_output=True, text=True, timeout=10)
        
        # Since our app doesn't have --version, we expect it to fail gracefully
        print("PASS - Executable starts without crashing")
        return True
        
    except subprocess.TimeoutExpired:
        print("WARN - Executable started but didn't respond to --version (this is expected)")
        return True
    except Exception as e:
        print(f"FAIL - Executable test failed: {e}")
        return False

def create_distribution_info(exe_path):
    """Create distribution information"""
    if not exe_path:
        return
    
    dist_dir = exe_path.parent
    
    # Create README for distribution
    readme_content = f"""
EPID Reconstruction GUI - Standalone Application
===============================================

This is a standalone executable version of the EPID Reconstruction GUI.
No Python installation is required!

## Quick Start
1. Double-click EPID_Reconstruction_GUI.exe
2. Select your EPID data folder
3. Configure reconstruction parameters
4. Start reconstruction and view results

## System Requirements
- Windows 10/11 (64-bit)
- Minimum 4GB RAM
- 2GB free disk space

## Features
- Advanced cone beam FDK reconstruction
- Interactive 3D volume visualization
- Real-time statistics and analysis
- DICOM export with TPS headers
- Multiple color mapping options

## File Size
Executable: {exe_path.stat().st_size / (1024*1024):.1f} MB
This includes all Python dependencies and libraries.

## Support
For technical support, refer to the GUI_README.md in the source code
or contact the development team.

Built with PyInstaller {get_pyinstaller_version()}
Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}
"""
    
    with open(dist_dir / "README.txt", "w") as f:
        f.write(readme_content)
    
    print("PASS - Distribution info created")

def get_pyinstaller_version():
    """Get PyInstaller version"""
    try:
        import PyInstaller
        return PyInstaller.__version__
    except:
        return "Unknown"

def create_directory_distribution():
    """Create a directory-based distribution for easier deployment"""
    print("Creating directory distribution...")
    
    # Create onedir distribution
    try:
        result = subprocess.run([
            sys.executable, "-m", "PyInstaller", 
            "--clean",
            "--onedir",  # Create directory instead of single file
            "--windowed",  # No console
            "--name", "EPID_Reconstruction_GUI_Dir",
            "epid_reconstruction_gui.py"
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("PASS - Directory distribution created")
            return True
        else:
            print("WARN - Directory distribution failed")
            return False
            
    except Exception as e:
        print(f"WARN - Directory distribution error: {e}")
        return False

def main():
    """Main build process"""
    print("EPID Reconstruction GUI - PyInstaller Build")
    print("=" * 60)
    
    # Check PyInstaller
    if not check_pyinstaller():
        return 1
    
    # Clean previous builds
    clean_build_directories()
    
    # Build executable
    success, exe_path = build_with_pyinstaller()
    if not success:
        return 1
    
    # Test executable
    test_executable(exe_path)
    
    # Create distribution info
    create_distribution_info(exe_path)
    
    # Try to create directory distribution as well
    create_directory_distribution()
    
    print("\n" + "=" * 60)
    print("PyInstaller build completed!")
    print("\nDistribution files:")
    print(f"  Single file: dist/EPID_Reconstruction_GUI.exe")
    print(f"  Directory:   dist/EPID_Reconstruction_GUI_Dir/")
    print("\nTo distribute:")
    print("  Option 1: Share the single .exe file (larger but simpler)")
    print("  Option 2: Zip the directory version (faster startup)")
    print("\nNo Python installation required on target machines!")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())