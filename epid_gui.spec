# -*- mode: python ; coding: utf-8 -*-

# PyInstaller spec file for EPID Reconstruction GUI
# This file configures the build process for creating a standalone executable

import sys
from pathlib import Path

# Analysis phase - what files and modules to include
a = Analysis(
    ['epid_reconstruction_gui.py'],  # Main script
    pathex=[],  # Additional paths to search for modules
    binaries=[],  # Binary files to include
    datas=[
        # Include the TPS reference data module
        ('tps_reference_data.py', '.'),
        # Include the TPS DICOM template file
        ('EPID_12_t0.dcm', '.'),
    ],
    hiddenimports=[
        # Explicitly include modules that might be missed
        'numpy',
        'scipy',
        'scipy.fft',
        'matplotlib',
        'matplotlib.backends.backend_qt5agg',
        'matplotlib.backends.backend_qtagg',
        'pydicom',
        'pydicom.encoders',
        'pydicom.decoders',
        'PySide6',
        'PySide6.QtCore',
        'PySide6.QtGui', 
        'PySide6.QtWidgets',
        'cv2',
        'tqdm',
        'concurrent.futures',
        'numba',
        'pathlib',
        'datetime',
        'warnings'
    ],
    hookspath=[],  # Additional hook directories
    hooksconfig={},  # Hook configuration
    runtime_hooks=[],  # Runtime hooks
    excludes=[
        # Exclude unnecessary modules to reduce size
        'tkinter',
        'unittest',
        'test',
        'email',
        'html',
        'http',
        'urllib',
        'xml',
        'pydoc',
        'doctest',
        'argparse',
        'pickle',
        'sqlite3'
    ],
    noarchive=False,
    optimize=0,
)

# PYZ phase - create Python archive
pyz = PYZ(a.pure, a.zipped_data, cipher=None)

# EXE phase - create executable
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='EPID_Reconstruction_GUI',
    debug=False,  # Set to True for debugging
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,  # Compress executable with UPX if available
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,  # Set to True if you want console window
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,  # Add icon file path if you have one (e.g., 'icon.ico')
)

# Optional: Create a directory distribution instead of single file
# COLLECT creates a directory with all files
# Uncomment the following lines if you prefer directory distribution:

# coll = COLLECT(
#     exe,
#     a.binaries,
#     a.datas,
#     strip=False,
#     upx=True,
#     upx_exclude=[],
#     name='EPID_Reconstruction_GUI'
# )