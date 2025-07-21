"""
Simple PyInstaller build script for EPID Reconstruction GUI
Using direct command approach instead of spec file
"""

import subprocess
import sys
import os
import shutil
from pathlib import Path
import time

def clean_build():
    """Clean previous build artifacts"""
    dirs_to_clean = ["build", "dist", "__pycache__"]
    for dir_name in dirs_to_clean:
        if Path(dir_name).exists():
            shutil.rmtree(dir_name)
    print("PASS - Build directories cleaned")

def build_onefile():
    """Build single file executable"""
    print("Building single file executable...")
    
    cmd = [
        sys.executable, "-m", "PyInstaller",
        "--onefile",
        "--windowed",
        "--name", "EPID_Reconstruction_GUI",
        "--add-data", "tps_reference_data.py;.",
        "--add-data", "EPID_12_t0.dcm;.",
        "--hidden-import", "numpy",
        "--hidden-import", "scipy.fft",
        "--hidden-import", "matplotlib.backends.backend_qtagg",
        "--hidden-import", "PySide6",
        "--hidden-import", "cv2",
        "--hidden-import", "pydicom",
        "--hidden-import", "tqdm",
        "--hidden-import", "numba",
        "--exclude-module", "tkinter",
        "--exclude-module", "unittest",
        "--exclude-module", "test",
        "epid_reconstruction_gui.py"
    ]
    
    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    build_time = time.time() - start_time
    
    if result.returncode == 0:
        print(f"PASS - Single file build completed in {build_time:.1f} seconds")
        return True
    else:
        print("FAIL - Single file build failed")
        print("STDERR:", result.stderr[-1000:])
        return False

def build_onedir():
    """Build directory distribution"""
    print("Building directory distribution...")
    
    cmd = [
        sys.executable, "-m", "PyInstaller",
        "--onedir",
        "--windowed", 
        "--name", "EPID_Reconstruction_GUI_Dir",
        "--add-data", "tps_reference_data.py;.",
        "--add-data", "EPID_12_t0.dcm;.",
        "--hidden-import", "numpy",
        "--hidden-import", "scipy.fft",
        "--hidden-import", "matplotlib.backends.backend_qtagg",
        "--hidden-import", "PySide6",
        "--hidden-import", "cv2",
        "--hidden-import", "pydicom",
        "--hidden-import", "tqdm",
        "--hidden-import", "numba",
        "--exclude-module", "tkinter",
        "--exclude-module", "unittest",
        "--exclude-module", "test",
        "epid_reconstruction_gui.py"
    ]
    
    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    build_time = time.time() - start_time
    
    if result.returncode == 0:
        print(f"PASS - Directory build completed in {build_time:.1f} seconds")
        return True
    else:
        print("WARN - Directory build failed")
        print("STDERR:", result.stderr[-500:])
        return False

def check_results():
    """Check build results"""
    results = []
    
    # Check single file
    onefile_path = Path("dist/EPID_Reconstruction_GUI.exe")
    if onefile_path.exists():
        size_mb = onefile_path.stat().st_size / (1024*1024)
        print(f"PASS - Single file: {onefile_path} ({size_mb:.1f} MB)")
        results.append(("Single file", str(onefile_path), f"{size_mb:.1f} MB"))
    
    # Check directory
    onedir_path = Path("dist/EPID_Reconstruction_GUI_Dir")
    if onedir_path.exists():
        exe_path = onedir_path / "EPID_Reconstruction_GUI_Dir.exe"
        if exe_path.exists():
            size_mb = exe_path.stat().st_size / (1024*1024)
            print(f"PASS - Directory: {onedir_path} (main exe: {size_mb:.1f} MB)")
            results.append(("Directory", str(onedir_path), f"{size_mb:.1f} MB"))
    
    return results

def main():
    """Main build process"""
    print("EPID Reconstruction GUI - Simple PyInstaller Build")
    print("=" * 60)
    
    # Check if required files exist
    required_files = ["epid_reconstruction_gui.py", "tps_reference_data.py", "EPID_12_t0.dcm"]
    for file in required_files:
        if not Path(file).exists():
            print(f"FAIL - Required file missing: {file}")
            return 1
    
    print("PASS - All required files found")
    
    # Clean build
    clean_build()
    
    # Try single file build
    onefile_success = build_onefile()
    
    # Try directory build
    onedir_success = build_onedir()
    
    # Check results
    results = check_results()
    
    print("\n" + "=" * 60)
    if results:
        print("PyInstaller build completed!")
        print("\nBuilt distributions:")
        for name, path, size in results:
            print(f"  {name}: {path} ({size})")
        
        print("\nDistribution notes:")
        print("- Single file: Slower startup but easier to distribute")
        print("- Directory: Faster startup but more files to manage")
        print("- TPS template (EPID_12_t0.dcm) included in both")
        print("- No Python installation required on target machines")
        
        return 0
    else:
        print("Build failed - no executables created")
        return 1

if __name__ == "__main__":
    sys.exit(main())