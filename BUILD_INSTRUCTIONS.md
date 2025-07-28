# EPID Reconstruction GUI - Build Instructions

## Quick Build

Simply run the build script:

```bash
python build_gui.py
```

## Requirements

Before building, ensure you have:

1. **Python environment** with all dependencies installed:
   - PySide6
   - numpy, scipy, matplotlib
   - pydicom, opencv-python, tqdm
   - PyInstaller

2. **Required files** in the project directory:
   - `epid_reconstruction_gui.py` (main application)
   - `tps_reference_data.py` (TPS reference module)
   - `EPID_12_t0.dcm` (DICOM template file)

## Manual Build Command

If you prefer to run PyInstaller directly:

```bash
pyinstaller --windowed --name="EPID_Reconstruction_GUI" --add-data "EPID_12_t0.dcm;." --add-data "tps_reference_data.py;." epid_reconstruction_gui.py
```

## Output

The build creates:
- **Directory**: `dist/EPID_Reconstruction_GUI/` (~434 MB)
- **Executable**: `EPID_Reconstruction_GUI.exe`
- **Dependencies**: All required libraries and files in `_internal/` folder

## Distribution

To distribute the application:
1. Copy the entire `dist/EPID_Reconstruction_GUI/` folder
2. Share with end users
3. Users run `EPID_Reconstruction_GUI.exe` (no Python installation required)

## Notes

- **Directory mode** (not onefile) is used for smaller size and better performance
- **DICOM template** (`EPID_12_t0.dcm`) is automatically included
- **TPS reference data** is embedded for dose verification
- **Windowed mode** means no console window appears