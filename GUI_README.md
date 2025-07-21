# EPID Dose Reconstruction GUI

Advanced PySide6 GUI application for EPID dose reconstruction using cone beam FDK algorithm with real-time dose verification.

## Features

### 🔬 **Advanced Reconstruction**
- **Cone Beam FDK Algorithm**: Based on corrected ConeBeam_FDK_corrected2.ipynb
- **Quick vs Full Mode**: Quick mode (every 4th slice) for fast testing, Full mode for complete reconstruction
- **Rotation Correction Toggle**: Enable/disable rotation correction when threshold exceeded
- **Optimized Performance**: Multi-threaded DICOM loading and chunked reconstruction  
- **Memory Efficient**: Processes large datasets without memory overflow
- **Real-time Progress**: Live progress tracking and status updates

### 🖼️ **Interactive Volume Viewer**
- **3D Slice Navigation**: View XY (axial), XZ (coronal), YZ (sagittal) planes
- **Mouse Cursor Tracking**: Real-time voxel position and dose value display
- **Scroll Wheel Navigation**: Navigate slices using mouse scroll in any view
- **Crosshair Display**: Visual guides showing slice intersections
- **Multiple Color Maps**: Hot, Jet, Viridis, Plasma, Coolwarm, CMRmap_r, Gray
- **Consistent Scaling**: Fixed colormap range across all slices to prevent false bright regions

### 📊 **Volume Statistics**
- **Real-time Statistics**: Min, max, mean, standard deviation
- **Coverage Analysis**: Non-zero voxel count and coverage percentage
- **Comprehensive Metrics**: Percentiles, coefficient of variation
- **Live Updates**: Statistics update automatically with new reconstructions

### 📁 **DICOM Export**
- **TPS Header Integration**: Uses TPS template DICOM headers for compatibility
- **Scale Factor Control**: Adjustable intensity scaling for export
- **Automatic Metadata**: Timestamps and series descriptions
- **One-click Export**: Simple export workflow with file browser

### 🎛️ **User Interface**
- **Clean Single-Panel Design**: All controls in one streamlined interface
- **Intuitive Controls**: Easy-to-use sliders and dropdown menus
- **Compact Layout**: Efficient use of screen space
- **Status Monitoring**: Real-time reconstruction status and error handling
- **Parameter Control**: Adjustable reconstruction settings

## Files Structure

```
dose-3d/
├── epid_reconstruction_gui.py      # Main GUI application
├── tps_reference_data.py           # TPS reference data and verification functions
├── test_gui_components.py          # Component testing script
├── launch_gui.py                   # Simple launcher script
├── GUI_README.md                   # This documentation
└── notebooks/
    ├── ConeBeam_FDK_corrected2.ipynb  # Source reconstruction algorithm
    └── dose_recon3d_quick.py           # Quick verification logic
```

## Requirements

### Python Environment
- **Python 3.8+**
- **Conda environment**: py12 (as specified in your setup)

### Required Packages
```bash
pip install PySide6 matplotlib numpy scipy pydicom numba opencv-python tqdm
```

### Optional Packages
- **numba**: For JIT compilation acceleration (recommended)
- **concurrent.futures**: For parallel processing (usually included)

## Usage

### Quick Start
1. **Activate Environment**:
   ```bash
   conda activate py12
   ```

2. **Run Component Tests**:
   ```bash
   python test_gui_components.py
   ```

3. **Launch GUI**:
   ```bash
   python epid_reconstruction_gui.py
   # or
   python launch_gui.py
   ```

### Reconstruction Workflow

1. **Select Data Folder**:
   - Click "Browse Data Folder"
   - Select folder containing EPID DICOM files
   - Example: `dataset/VMAT 2025 - 6. SIB COMPLEX TARGET/T1/873251691`

2. **Configure Parameters**:
   - **Mode**: Quick (every 4th slice, 50³ volume) or Full (all slices, 100³ volume)
   - **Image Size**: Reconstruction volume size (auto-adjusted by mode)
   - **Chunk Size**: Memory management parameter (auto-adjusted by mode)
   - **Rotation Correction**: Enable/disable rotation correction when threshold exceeded

3. **Start Reconstruction**:
   - Click "Start Reconstruction"
   - Monitor progress bar and status messages
   - Reconstruction runs in background thread

4. **View Results**:
   - Switch to "Volume Viewer" tab
   - Use sliders to navigate slices
   - Mouse over images to see voxel values
   - Scroll wheel to change slices

5. **Review Statistics**:
   - View real-time statistics in the Volume Statistics panel
   - Monitor min/max values, coverage percentage
   - Check non-zero voxel distribution

6. **Export Results**:
   - Set desired scale factor for intensity values
   - Click "Browse..." to select output path
   - Click "Export DICOM" to save with TPS headers

### Interactive Controls

#### Volume Viewer
- **Slice Sliders**: X, Y, Z slice selection (0-99)
- **Colormap Dropdown**: Select visualization color scheme
- **Crosshairs Checkbox**: Toggle intersection guides
- **Mouse Cursor**: Shows real-time position and dose value
- **Scroll Wheel**: Navigate slices in focused view

#### Statistics & Export
- **Statistics Panel**: Real-time volume statistics display
- **Scale Factor**: Adjust intensity scaling for DICOM export
- **Export Path**: File browser for output location selection
- **Export Button**: One-click DICOM export with TPS headers

## Technical Details

### Reconstruction Algorithm
Based on the corrected FDK (Feldkamp-Davis-Kress) cone beam reconstruction:

1. **DICOM Loading**: Multi-threaded loading with progress tracking
2. **Differential Processing**: Extract dose increments from cumulative images  
3. **Geometric Correction**: Apply inverse square law and magnification
4. **Filtering**: Shepp-Logan filter in frequency domain
5. **Backprojection**: 3D volume reconstruction with bilinear interpolation

### DICOM Export System
- **TPS Template**: Uses EPID_12_t0.dcm as header template
- **Header Copying**: Preserves DICOM metadata for compatibility
- **Automatic Updates**: Timestamps and series descriptions
- **Scale Control**: Configurable intensity scaling for different systems

### Performance Optimizations
- **Threading**: Parallel DICOM loading
- **Chunking**: Memory-efficient reconstruction
- **Numba JIT**: Optional acceleration for compute-intensive functions
- **Caching**: TPS reference data cached in memory

## Troubleshooting

### Common Issues

**Import Errors**:
```bash
# Test components first
python test_gui_components.py

# If matplotlib/Qt issues:
pip install --upgrade matplotlib PySide6
```

**Memory Issues**:
- Reduce Image Size (50 instead of 100)
- Increase Chunk Size for fewer chunks
- Close other applications

**Performance Issues**:
- Install numba: `pip install numba`
- Use Quick mode for faster testing
- Enable hardware acceleration if available

**Colormap Issues**:
- Fixed: Bright regions in empty areas now correctly scaled
- All slices use consistent colormap range (vmin/vmax from entire volume)
- Toggle different colormaps to see dose distribution clearly

**Rotation Correction**:
- Enable: Applies cv2.flip(-abs(curr-prev), 1) when threshold exceeded
- Disable: Uses previous image data when threshold exceeded  
- Only affects processing when np.max(difference) > threshold

**GUI Not Displaying**:
- Check if running in headless environment
- Verify PySide6 installation
- Try different matplotlib backend

### Export Troubleshooting

**"TPS Template Not Available"**:
- Check EPID_12_t0.dcm exists in dataset/3DDose/
- Verify file permissions and accessibility
- Update path in tps_reference_data.py if needed

**Export Fails**:
- Check output directory permissions
- Ensure adequate disk space
- Verify scale factor is reasonable (0.1-10.0)

**DICOM Compatibility**:
- Exported DICOM uses original TPS headers
- Compatible with most medical imaging software
- Scale factor affects intensity values only

## Development

### Adding New Features

**New Color Maps**:
```python
# In InteractiveVolumeViewer.setup_ui()
self.colormap_combo.addItems(['hot', 'jet', 'new_colormap'])
```

**Custom Verification Metrics**:
```python
# In tps_reference_data.py, modify calculate_dose_metrics()
def calculate_dose_metrics(reconstructed, reference):
    # Add new metric calculations
    custom_metric = your_calculation(reconstructed, reference)
    metrics['custom_metric'] = custom_metric
```

**Additional Views**:
```python
# In InteractiveVolumeViewer, add new subplot
ax4 = self.figure.add_subplot(224)
# Implement new visualization
```

### Code Structure

- **`EPIDReconstructionGUI`**: Main window class
- **`ReconstructionWorker`**: Background reconstruction thread  
- **`InteractiveVolumeViewer`**: 3D volume visualization widget
- **`DoseVerificationWidget`**: TPS comparison interface
- **`OptimizedFDKReconstructor`**: Core reconstruction algorithm

## Support

For issues or questions:
1. Check component tests: `python test_gui_components.py`
2. Review error messages in GUI status panel
3. Verify all dependencies are installed
4. Check DICOM data integrity

## Citation

Based on cone beam FDK reconstruction algorithm with optimizations for EPID dose reconstruction applications.