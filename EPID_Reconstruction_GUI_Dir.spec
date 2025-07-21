# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['epid_reconstruction_gui.py'],
    pathex=[],
    binaries=[],
    datas=[('tps_reference_data.py', '.'), ('EPID_12_t0.dcm', '.')],
    hiddenimports=['numpy', 'scipy.fft', 'matplotlib.backends.backend_qtagg', 'PySide6', 'cv2', 'pydicom', 'tqdm', 'numba'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['tkinter', 'unittest', 'test'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='EPID_Reconstruction_GUI_Dir',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='EPID_Reconstruction_GUI_Dir',
)
