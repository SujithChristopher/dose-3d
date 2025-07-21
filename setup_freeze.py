"""
cx_Freeze setup script for EPID Reconstruction GUI
Creates a standalone executable from the Python application
"""

import sys
from cx_Freeze import setup, Executable
import os

# Simplified build options to avoid circular dependencies
build_exe_options = {
    "packages": [
        "numpy", "scipy", "matplotlib", "pydicom", "PySide6", 
        "cv2", "tqdm", "concurrent.futures"
    ],
    "excludes": [
        "tkinter", "unittest", "email", "html", "http", "urllib", 
        "xml", "pydoc", "doctest", "argparse", "pickle", "test",
        "distutils", "setuptools", "pip", "wheel"
    ],
    "include_files": [
        # Include the TPS reference data module
        ("tps_reference_data.py"),
    ],
    "build_exe": "build/epid_reconstruction_gui",
    "optimize": 2,
    "zip_include_packages": ["*"],
    "zip_exclude_packages": []
}

# GUI applications require a different base on Windows
base = None
if sys.platform == "win32":
    base = "Win32GUI"  # Use Win32GUI for windowed applications

# Main executable
executables = [
    Executable(
        "epid_reconstruction_gui.py",
        base=base,
        target_name="EPID_Reconstruction_GUI.exe",
        icon=None,  # Add icon file path if you have one
        shortcut_name="EPID Reconstruction GUI",
        shortcut_dir="DesktopFolder"
    )
]

setup(
    name="EPID Reconstruction GUI",
    version="1.0.0",
    description="Advanced EPID Dose Reconstruction with Interactive Visualization",
    author="EPID Reconstruction Team",
    options={"build_exe": build_exe_options},
    executables=executables
)