"""
Minimal cx_Freeze setup for EPID GUI - avoiding problematic packages
"""

import sys
from cx_Freeze import setup, Executable

# Minimal build options
build_exe_options = {
    "packages": ["numpy", "matplotlib", "PySide6"],
    "excludes": [
        "tkinter", "unittest", "test", "distutils", "setuptools", 
        "pip", "wheel", "email", "html", "http", "urllib", "xml"
    ],
    "include_files": [("tps_reference_data.py", "tps_reference_data.py")],
    "build_exe": "build_minimal",
    "silent": True
}

# GUI base for Windows
base = "Win32GUI" if sys.platform == "win32" else None

setup(
    name="EPID_GUI_Minimal",
    version="1.0",
    description="EPID Reconstruction GUI - Minimal Build",
    options={"build_exe": build_exe_options},
    executables=[Executable("epid_reconstruction_gui.py", base=base, target_name="epid_gui.exe")]
)