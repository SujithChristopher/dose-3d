# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an EPID (Electronic Portal Imaging Device) dose reconstruction project that reconstructs 3D dose distributions from linear accelerator imaging data using cone beam CT reconstruction algorithms. The project is primarily developed using Jupyter notebooks for iterative research and analysis.

## Core Architecture

### Main Components

- **epid_dose_reconstruction.py**: Complete 3D dose reconstruction pipeline with the `EPIDDoseReconstructor` class
- **main.py**: Simple GUI entry point using PySide6
- **analyze_epid_dicom.py**: DICOM analysis utilities
- **notebooks/**: Research notebooks for algorithm development and testing
- **dataset/**: DICOM data from linear accelerator EPID acquisitions
- **results/**: Output files (3D dose volumes, reconstruction metadata)

### Key Algorithms

The reconstruction process handles cumulative dose delivery from EPID images:
1. **Dose Increment Extraction**: Converts cumulative EPID images to dose increments per gantry angle
2. **Geometric Correction**: Applies inverse square law and magnification corrections
3. **3D Backprojection**: Projects dose increments into 3D volume using cone beam geometry
4. **Dose Conservation**: Validates total dose preservation throughout reconstruction

### Data Structure

- EPID images are sequential DICOM files representing cumulative dose at different gantry angles
- Each dataset contains hundreds to thousands of projection images
- Reconstruction produces 3D dose volumes (typically 50³ voxels)
- Results saved as NumPy arrays with JSON metadata

## Development Environment

### Python Dependencies
- Scientific: numpy, scipy, matplotlib
- Medical imaging: pydicom
- Optional performance: numba (JIT compilation), scikit-image (advanced reconstruction)
- GUI: PySide6
- Parallel processing: concurrent.futures
- Progress tracking: tqdm

### Common Commands

**Run main reconstruction pipeline:**
```bash
python epid_dose_reconstruction.py
```

**Launch GUI interface:**
```bash
python main.py
```

**Start Jupyter for notebook development:**
```bash
jupyter notebook
# or
jupyter lab
```

**Process specific dataset (modify path in script):**
Edit the `epid_path` variable in `epid_dose_reconstruction.py` main() function

### Key Notebooks

- `EPID_Dose_Reconstruction_3D.ipynb`: Main reconstruction development
- `ConeBeam_FDK_*.ipynb`: FDK algorithm variants and optimizations
- `dose_check.ipynb`: Dose validation and analysis
- `reading_the_files.ipynb`: DICOM file handling

## Working with EPID Data

### Dataset Organization
```
dataset/
├── VMAT 2025 - 1. CENTRAL TARGET/
├── VMAT 2025 - 6. SIB COMPLEX TARGET/
├── epid-1-arc-vmat/
└── 3DDose/
```

### Reconstruction Parameters
- SAD (Source-to-Axis Distance): typically 1000mm
- SID (Source-to-Image Distance): typically 1500mm  
- EPID pixel size: typically 0.172mm
- Reconstruction volume: configurable (default 50³)

### Performance Considerations
- Large datasets (1000+ images) require chunked processing
- Parallel loading with ThreadPoolExecutor (default 6 workers)
- Memory management for 3D volume reconstruction
- Optional Numba JIT compilation for speed

## Research Workflow

1. **Data Exploration**: Use `reading_the_files.ipynb` to examine new DICOM datasets
2. **Algorithm Development**: Iterate in relevant ConeBeam notebooks
3. **Reconstruction Testing**: Run `epid_dose_reconstruction.py` with new parameters
4. **Validation**: Use `dose_check.ipynb` for quality assessment
5. **Results Analysis**: Examine outputs in `results/` directory

The project emphasizes iterative development through notebooks while maintaining production-ready reconstruction pipelines in the main Python scripts.