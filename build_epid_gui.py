"""
PyInstaller build script for EPID Reconstruction GUI
Directory distribution only (no --onefile for smaller size)
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
            print(f"Cleaned: {dir_name}")
    print("✓ Build directories cleaned")

def build_gui():
    """Build EPID GUI as directory distribution"""
    print("Building EPID Reconstruction GUI...")
    
    cmd = [
        sys.executable, "-m", "PyInstaller",
        "--onedir",  # Directory distribution only
        "--windowed",  # No console window
        "--name", "EPID_Reconstruction_GUI",
        "--add-data", "tps_reference_data.py;.",
        "--hidden-import", "numpy",
        "--hidden-import", "scipy.fft", 
        "--hidden-import", "matplotlib.backends.backend_qtagg",
        "--hidden-import", "matplotlib.backends.backend_qt5agg",
        "--hidden-import", "PySide6.QtWidgets",
        "--hidden-import", "PySide6.QtCore", 
        "--hidden-import", "PySide6.QtGui",
        "--hidden-import", "cv2",
        "--hidden-import", "pydicom",
        "--hidden-import", "tqdm",
        "--hidden-import", "numba",
        "--hidden-import", "concurrent.futures",
        "--exclude-module", "tkinter",
        "--exclude-module", "unittest",
        "--exclude-module", "test",
        "--exclude-module", "jupyter",
        "--exclude-module", "IPython",
        "--exclude-module", "notebook",
        "--exclude-module", "torch",
        "--exclude-module", "torchvision",
        "--exclude-module", "tensorflow",
        "--exclude-module", "keras",
        "--exclude-module", "sklearn",
        "--exclude-module", "pandas",
        "--exclude-module", "seaborn",
        "--exclude-module", "plotly",
        "epid_reconstruction_gui.py"
    ]
    
    # Add TPS template if it exists
    if Path("EPID_12_t0.dcm").exists():
        cmd.extend(["--add-data", "EPID_12_t0.dcm;."])
        print("✓ Including TPS template: EPID_12_t0.dcm")
    
    start_time = time.time()
    print("Running PyInstaller...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    build_time = time.time() - start_time
    
    if result.returncode == 0:
        print(f"✓ Build completed in {build_time:.1f} seconds")
        return True
    else:
        print("✗ Build failed")
        print("Error output:")
        print(result.stderr[-1000:])
        return False

def check_build_result():
    """Check and report build results"""
    dist_dir = Path("dist/EPID_Reconstruction_GUI")
    exe_path = dist_dir / "EPID_Reconstruction_GUI.exe"
    
    if dist_dir.exists() and exe_path.exists():
        # Calculate directory size
        total_size = sum(f.stat().st_size for f in dist_dir.rglob('*') if f.is_file())
        size_mb = total_size / (1024*1024)
        
        exe_size_mb = exe_path.stat().st_size / (1024*1024)
        file_count = len(list(dist_dir.rglob('*')))
        
        print(f"✓ Build successful!")
        print(f"  Location: {dist_dir}")
        print(f"  Total size: {size_mb:.1f} MB")
        print(f"  Main executable: {exe_size_mb:.1f} MB")
        print(f"  Files: {file_count}")
        
        return True
    else:
        print("✗ Build failed - executable not found")
        return False

def create_launcher_script():
    """Create a simple launcher script"""
    launcher_content = '''@echo off
cd /d "%~dp0"
cd dist\\EPID_Reconstruction_GUI
EPID_Reconstruction_GUI.exe
pause
'''
    
    with open("launch_epid_gui.bat", "w") as f:
        f.write(launcher_content)
    print("✓ Created launcher script: launch_epid_gui.bat")

def main():
    """Main build process"""
    print("EPID Reconstruction GUI - PyInstaller Build")
    print("=" * 50)
    
    # Check if main file exists
    if not Path("epid_reconstruction_gui.py").exists():
        print("✗ epid_reconstruction_gui.py not found")
        return 1
    
    print("✓ Main GUI file found")
    
    # Check optional dependencies
    if Path("tps_reference_data.py").exists():
        print("✓ TPS reference module found")
    else:
        print("⚠ TPS reference module not found (optional)")
    
    # Clean previous builds
    clean_build()
    
    # Build GUI
    if not build_gui():
        return 1
    
    # Check results
    if not check_build_result():
        return 1
    
    # Create launcher
    create_launcher_script()
    
    print("\n" + "=" * 50)
    print("Build complete! Directory distribution created.")
    print("\nTo run:")
    print("  1. Use: launch_epid_gui.bat")
    print("  2. Or navigate to: dist/EPID_Reconstruction_GUI/")
    print("  3. Run: EPID_Reconstruction_GUI.exe")
    print("\nDirectory benefits:")
    print("  - Smaller total size than --onefile")
    print("  - Faster startup time")
    print("  - Easier to debug if needed")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())