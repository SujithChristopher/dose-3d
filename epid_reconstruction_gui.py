"""
EPID Dose Reconstruction GUI
Advanced PySide6 GUI for cone beam FDK reconstruction with dose verification.

Features:
- Interactive 3D volume viewer with slice navigation
- Multiple color mapping options
- Real-time dose verification against TPS reference
- Mouse cursor position tracking
- Scroll wheel slice navigation
- Optimized reconstruction pipeline
"""

import sys
import os
import time
import math
import warnings
from pathlib import Path

import numpy as np
import matplotlib
# Try to use appropriate Qt backend for PySide6
try:
    matplotlib.use('QtAgg')  # Modern backend that should work with PySide6
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
except ImportError:
    try:
        matplotlib.use('Qt5Agg')  # Fallback to Qt5Agg
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    except ImportError:
        matplotlib.use('Agg')  # Final fallback
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import matplotlib.patches as patches

from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                              QHBoxLayout, QGridLayout, QPushButton, QLabel, 
                              QSlider, QComboBox, QProgressBar, QTextEdit,
                              QFileDialog, QGroupBox, QSpinBox, QDoubleSpinBox,
                              QCheckBox, QSplitter, QTabWidget, QTableWidget,
                              QTableWidgetItem, QStatusBar, QMessageBox)
from PySide6.QtCore import Qt, QThread, QTimer, Signal, QObject
from PySide6.QtGui import QFont, QPixmap, QCursor

# Import reconstruction modules
try:
    from pydicom import dcmread
    from scipy.fft import fft, ifft
    from concurrent.futures import ThreadPoolExecutor
    from tqdm import tqdm
    import cv2
    PYDICOM_AVAILABLE = True
except ImportError as e:
    print(f"Import warning: {e}")
    PYDICOM_AVAILABLE = False

# Import TPS utilities
try:
    from tps_reference_data import export_reconstructed_dicom, calculate_volume_statistics, load_tps_dicom_template
    TPS_AVAILABLE = True
except ImportError:
    print("TPS utilities not available")
    TPS_AVAILABLE = False

warnings.filterwarnings('ignore')

# Try numba for performance
try:
    from numba import jit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def jit(func):
        return func
    prange = range

PI = math.pi

class OptimizedFDKReconstructor:
    """Optimized FDK reconstruction from ConeBeam_FDK_corrected2.ipynb"""
    
    def __init__(self, SOD=1000.0, delta_dd=0.1075, Nimage=100):
        self.SOD = SOD
        self.delta_dd = delta_dd
        self.Nimage = Nimage
        
    def filter_SL(self, N, d):
        """Shepp-Logan filter"""
        fh_SL = np.zeros(N)
        for k1 in range(0, N, 1):
            fh_SL[k1] = -2.0/(PI*PI*d*d*(4*(k1-N/2.0)**2-1))
        return fh_SL

    def nearestPowerOf2(self, N):
        """Find nearest power of 2"""
        a = int(math.log2(N))
        if 2**a == N:
            return N
        return 2**(a + 1)

    def Fun_Weigth_Projection(self, projection_beta, SOD, delta_dd):
        """Weight projection"""
        Nrows, Ncolumns = projection_beta.shape
        dd_column = delta_dd*np.arange(-Ncolumns/2+0.5, (Ncolumns/2+1)-0.5, 1.0)
        dd_row = delta_dd*np.arange(-Nrows/2+0.5, (Nrows/2+1)-0.5, 1.0)
        dd_row2D, dd_column2D = np.meshgrid(dd_row, dd_column, indexing='ij')
        weighted_projection = projection_beta*SOD/np.sqrt(SOD*SOD+np.power(dd_row2D, 2.0)+np.power(dd_column2D, 2.0))
        return weighted_projection

    def optimize_convolution(self, weighted_projection, fh_RL):
        """Optimized convolution"""
        Nrows, Ncolumns = weighted_projection.shape
        Nfft = self.nearestPowerOf2(2 * Ncolumns - 1)
        fh_RL_padded = np.zeros(Nfft)
        fh_RL_padded[:len(fh_RL)] = fh_RL / 2.0
        
        fh_RL_fft = fft(fh_RL_padded)
        
        projection_padded = np.zeros((Nrows, Nfft))
        projection_padded[:, :Ncolumns] = weighted_projection

        projection_fft = fft(projection_padded, axis=1)
        convoluted_freq = projection_fft * fh_RL_fft
        convoluted_time = ifft(convoluted_freq, axis=1).real
        filtered_projection = convoluted_time[:, :Ncolumns]
        
        return filtered_projection

    def Fun_BackProjection(self, filtered_projection, SOD, beta_num, beta_m, delta_dd, Nimage):
        """Backprojection"""
        Nrows, Ncolumns = filtered_projection.shape
        MX, MZ = Nimage, int(Nimage*Nrows/Ncolumns)
        
        roi = delta_dd*np.array([-Ncolumns/2.0+0.5, Ncolumns/2.0-0.5, -Nrows/2.0+0.5, Nrows/2.0-0.5])
        hx = (roi[1]-roi[0])/(MX-1)
        xrange = roi[0]+hx*np.arange(0, MX)
        hy = (roi[3]-roi[2])/(MZ-1)
        yrange = roi[2]+hy*np.arange(0, MZ)
        XX, YY, ZZ = np.meshgrid(xrange, xrange, yrange, indexing='ij')
        temp_rec = np.zeros((MX, MX, MZ))
        U = (SOD+XX*np.sin(beta_m)-YY*np.cos(beta_m))/SOD
        a = (XX*np.cos(beta_m)+YY*np.sin(beta_m))/U
        xx = np.int32(np.floor(a/delta_dd))
        u1 = a/delta_dd-xx
        b = ZZ/U
        yy = np.int32(np.floor(b/delta_dd))
        u2 = b/delta_dd-yy
        xx = xx+int(Ncolumns/2)
        yy = yy+int(Nrows/2)

        mask = np.where((xx >= 0) & (xx < Ncolumns-1) & (yy >= 0) & (yy < Nrows-1))
        xx = xx[mask]
        yy = yy[mask]
        u1 = u1[mask]
        u2 = u2[mask]
        
        temp = ((1-u1)*(1-u2)*filtered_projection[yy, xx]+(1-u1)*u2*filtered_projection[yy+1, xx]+
                (1-u2)*u1*filtered_projection[yy, xx+1]+u1*u2*filtered_projection[yy+1, xx+1])
        temp_rec[mask] = temp_rec[mask]+temp/(np.power(U[mask], 2))*2*PI/beta_num
        
        return temp_rec

class ReconstructionWorker(QThread):
    """Worker thread for reconstruction to keep GUI responsive"""
    progress_update = Signal(int)
    status_update = Signal(str)
    result_ready = Signal(object)
    error_occurred = Signal(str)
    
    def __init__(self, data_path, reconstruction_params):
        super().__init__()
        self.data_path = data_path
        self.params = reconstruction_params
        
    def load_dicom_fast(self, path_list, max_workers=6):
        """Fast DICOM loading"""
        def read_dicom_data(fname):
            try:
                dcm = dcmread(fname)
                return dcm.pixel_array.astype(np.uint16), float(dcm.GantryAngle)
            except Exception as e:
                return None, None
        
        images, angles = [], []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_path = {executor.submit(read_dicom_data, path): path for path in path_list}
            
            for i, future in enumerate(future_to_path):
                img, angle = future.result()
                if img is not None:
                    images.append(img)
                    angles.append(angle)
                self.progress_update.emit(int(20 * i / len(path_list)))
        
        return np.array(images), np.array(angles)
    
    def process_differential_original(self, images, angles, threshold=10000, enable_rotation=True):
        """Process differential images with optional rotation correction"""
        n_images = len(images)
        shape = images[0].shape
        processed_images = np.zeros((n_images, shape[0], shape[0]), dtype=np.uint16)
        processed_angles = []
        
        prev = np.zeros((shape[0], shape[0]), dtype=np.uint16)
        
        for idx in range(n_images):
            curr = images[idx]
            _m = curr - prev
            
            if np.max(_m) > threshold:
                if enable_rotation:
                    # Apply rotation correction (flip and abs) when threshold exceeded
                    processed_images[idx, :, :] = cv2.flip(-abs(curr - prev), 1)
                    processed_angles.append(angles[idx])
                else:
                    # No rotation - use previous image data when threshold exceeded
                    if idx > 0:
                        processed_images[idx, :, :] = processed_images[idx-1, :, :]
                        processed_angles.append(processed_angles[idx-1])
                    else:
                        processed_images[idx, :, :] = curr - prev
                        processed_angles.append(angles[idx])
            else:
                # Normal differential processing (no rotation applied regardless of flag)
                processed_images[idx, :, :] = curr - prev
                prev = curr
                processed_angles.append(angles[idx])
            
            prev = curr
            self.progress_update.emit(20 + int(30 * idx / n_images))
        
        return processed_images, np.array(processed_angles)
    
    def run(self):
        try:
            self.status_update.emit("Loading DICOM files...")
            
            # Get DICOM files
            files = [f for f in os.listdir(self.data_path) if f.endswith('.dcm')]
            files.sort()  # Ensure consistent ordering
            
            # Apply every_nth sampling for quick mode
            every_nth = self.params.get('every_nth', 1)
            if every_nth > 1:
                files = files[::every_nth]
                self.status_update.emit(f"Quick mode: Using every {every_nth} files ({len(files)} total)")
            
            paths = [os.path.join(self.data_path, f) for f in files]
            
            # Load images
            raw_images, raw_angles = self.load_dicom_fast(paths)
            
            self.status_update.emit("Processing differential images...")
            processed_images, processed_angles = self.process_differential_original(
                raw_images, raw_angles, enable_rotation=self.params.get('enable_rotation', True))
            
            # Get reconstruction parameters from first DICOM
            dcm = dcmread(paths[0])
            SID = dcm.RTImageSID
            SAD = dcm.RadiationMachineSAD
            SOD = SAD
            SDD = SID
            width = 0.172
            delta_dd = width * SOD / SDD
            
            # Create reconstructor
            reconstructor = OptimizedFDKReconstructor(SOD, delta_dd, self.params['image_size'])
            
            self.status_update.emit("Running FDK reconstruction...")
            self.progress_update.emit(50)
            
            # Reconstruct
            rec_image = self.reconstruct_chunked(reconstructor, processed_images, processed_angles)
            
            self.progress_update.emit(100)
            self.status_update.emit("Reconstruction complete!")
            
            result = {
                'reconstructed_volume': rec_image,
                'parameters': {
                    'SOD': SOD,
                    'SDD': SDD,
                    'delta_dd': delta_dd,
                    'n_projections': len(processed_angles),
                    'image_size': self.params['image_size']
                }
            }
            
            self.result_ready.emit(result)
            
        except Exception as e:
            self.error_occurred.emit(str(e))
    
    def reconstruct_chunked(self, reconstructor, projections, angles, chunk_size=25):
        """Chunked reconstruction"""
        n_angles = len(angles)
        beta_rad = angles * PI / 180.0
        
        # Prepare filter
        Ncolumns = projections.shape[2]
        Nrows = projections.shape[1]
        Nfft = reconstructor.nearestPowerOf2(2 * Ncolumns - 1)
        fh_RL = reconstructor.filter_SL(Nfft, reconstructor.delta_dd)
        
        # Initialize result
        MX = reconstructor.Nimage
        MZ = int(reconstructor.Nimage * Nrows / Ncolumns)
        rec_image = np.zeros((MX, MX, MZ))
        
        # Process in chunks
        chunks = [list(range(i, min(i + chunk_size, n_angles))) 
                 for i in range(0, n_angles, chunk_size)]
        
        for chunk_idx, chunk_indices in enumerate(chunks):
            chunk_result = np.zeros((MX, MX, MZ))
            
            for angle_idx in chunk_indices:
                projection_beta = projections[angle_idx, :, :]
                
                weighted_projection = reconstructor.Fun_Weigth_Projection(
                    projection_beta, reconstructor.SOD, reconstructor.delta_dd)
                
                filtered_projection = reconstructor.optimize_convolution(weighted_projection, fh_RL)
                
                temp_rec = reconstructor.Fun_BackProjection(
                    weighted_projection, reconstructor.SOD, n_angles, 
                    beta_rad[angle_idx], reconstructor.delta_dd, reconstructor.Nimage)
                
                chunk_result += temp_rec
            
            rec_image += chunk_result
            
            # Update progress
            progress = 50 + int(40 * (chunk_idx + 1) / len(chunks))
            self.progress_update.emit(progress)
        
        return rec_image

class InteractiveVolumeViewer(QWidget):
    """Interactive 3D volume viewer with slice navigation, colormaps, and statistics"""
    
    def __init__(self):
        super().__init__()
        self.volume = None
        self.current_slice_x = 50
        self.current_slice_y = 50  
        self.current_slice_z = 50
        self.current_colormap = 'hot'
        self.show_crosshairs = True
        self.cursor_position = (0, 0)
        self.volume_stats = {}
        
        self.setup_ui()
        
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Controls
        controls_layout = QHBoxLayout()
        
        # Colormap selection
        colormap_group = QGroupBox("Colormap")
        colormap_layout = QHBoxLayout(colormap_group)
        self.colormap_combo = QComboBox()
        self.colormap_combo.addItems(['hot', 'jet', 'viridis', 'plasma', 'coolwarm', 'CMRmap_r', 'gray'])
        self.colormap_combo.currentTextChanged.connect(self.update_colormap)
        colormap_layout.addWidget(self.colormap_combo)
        controls_layout.addWidget(colormap_group)
        
        # View controls
        view_group = QGroupBox("View Controls")
        view_layout = QHBoxLayout(view_group)
        self.crosshairs_check = QCheckBox("Show Crosshairs")
        self.crosshairs_check.setChecked(True)
        self.crosshairs_check.stateChanged.connect(self.toggle_crosshairs)
        view_layout.addWidget(self.crosshairs_check)
        controls_layout.addWidget(view_group)
        
        layout.addLayout(controls_layout)
        
        # Slice controls
        slice_controls = QHBoxLayout()
        
        # X slice
        x_group = QGroupBox("X Slice")
        x_layout = QVBoxLayout(x_group)
        self.x_slider = QSlider(Qt.Horizontal)
        self.x_slider.setRange(0, 99)
        self.x_slider.setValue(50)
        self.x_slider.valueChanged.connect(self.update_x_slice)
        self.x_label = QLabel("X: 50")
        x_layout.addWidget(self.x_label)
        x_layout.addWidget(self.x_slider)
        slice_controls.addWidget(x_group)
        
        # Y slice  
        y_group = QGroupBox("Y Slice")
        y_layout = QVBoxLayout(y_group)
        self.y_slider = QSlider(Qt.Horizontal)
        self.y_slider.setRange(0, 99)
        self.y_slider.setValue(50)
        self.y_slider.valueChanged.connect(self.update_y_slice)
        self.y_label = QLabel("Y: 50")
        y_layout.addWidget(self.y_label)
        y_layout.addWidget(self.y_slider)
        slice_controls.addWidget(y_group)
        
        # Z slice
        z_group = QGroupBox("Z Slice")
        z_layout = QVBoxLayout(z_group)
        self.z_slider = QSlider(Qt.Horizontal)
        self.z_slider.setRange(0, 99)
        self.z_slider.setValue(50)
        self.z_slider.valueChanged.connect(self.update_z_slice)
        self.z_label = QLabel("Z: 50")
        z_layout.addWidget(self.z_label)
        z_layout.addWidget(self.z_slider)
        slice_controls.addWidget(z_group)
        
        layout.addLayout(slice_controls)
        
        # Matplotlib figure
        self.figure = Figure(figsize=(12, 4))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        layout.addWidget(self.canvas)
        
        # Cursor position label
        self.cursor_label = QLabel("Cursor: (0, 0) | Value: 0.00")
        layout.addWidget(self.cursor_label)
        
        # Statistics panel
        stats_group = QGroupBox("Volume Statistics")
        stats_layout = QGridLayout(stats_group)
        
        # Basic stats
        self.min_label = QLabel("Min: --")
        self.max_label = QLabel("Max: --")
        self.mean_label = QLabel("Mean: --")
        self.std_label = QLabel("Std: --")
        
        stats_layout.addWidget(self.min_label, 0, 0)
        stats_layout.addWidget(self.max_label, 0, 1)
        stats_layout.addWidget(self.mean_label, 1, 0)
        stats_layout.addWidget(self.std_label, 1, 1)
        
        # Coverage stats
        self.coverage_label = QLabel("Coverage: --%")
        self.nonzero_label = QLabel("Non-zero: --")
        stats_layout.addWidget(self.coverage_label, 2, 0)
        stats_layout.addWidget(self.nonzero_label, 2, 1)
        
        layout.addWidget(stats_group)
        
    def set_volume(self, volume):
        """Set the volume to display"""
        self.volume = volume
        if volume is not None:
            # Set slider ranges based on actual volume dimensions
            max_x = volume.shape[0] - 1
            max_y = volume.shape[1] - 1
            max_z = volume.shape[2] - 1
            
            self.x_slider.setRange(0, max_x)
            self.y_slider.setRange(0, max_y)
            self.z_slider.setRange(0, max_z)
            
            # Set current slices to center or clamp to valid range
            self.current_slice_x = min(max_x // 2, max_x)
            self.current_slice_y = min(max_y // 2, max_y)
            self.current_slice_z = min(max_z // 2, max_z)
            
            # Update slider positions
            self.x_slider.setValue(self.current_slice_x)
            self.y_slider.setValue(self.current_slice_y)
            self.z_slider.setValue(self.current_slice_z)
            
            # Update labels
            self.x_label.setText(f"X: {self.current_slice_x}")
            self.y_label.setText(f"Y: {self.current_slice_y}")
            self.z_label.setText(f"Z: {self.current_slice_z}")
            
            # Calculate and update statistics
            self.update_statistics()
            self.update_display()
    
    def update_statistics(self):
        """Update volume statistics display"""
        if self.volume is not None and TPS_AVAILABLE:
            self.volume_stats = calculate_volume_statistics(self.volume)
            
            # Update basic stats
            self.min_label.setText(f"Min: {self.volume_stats['min_value']:.2f}")
            self.max_label.setText(f"Max: {self.volume_stats['max_value']:.2f}")
            self.mean_label.setText(f"Mean: {self.volume_stats['mean_value']:.2f}")
            self.std_label.setText(f"Std: {self.volume_stats['std_value']:.2f}")
            
            # Update coverage stats
            self.coverage_label.setText(f"Coverage: {self.volume_stats['coverage_percent']:.1f}%")
            self.nonzero_label.setText(f"Non-zero: {self.volume_stats['nonzero_voxels']}")
        else:
            # Reset labels
            self.min_label.setText("Min: --")
            self.max_label.setText("Max: --")
            self.mean_label.setText("Mean: --")
            self.std_label.setText("Std: --")
            self.coverage_label.setText("Coverage: --%")
            self.nonzero_label.setText("Non-zero: --")
    
    def set_reference_volume(self, reference):
        """Set reference volume for comparison"""
        self.reference_volume = reference
        self.update_display()
        
    def update_colormap(self, colormap):
        """Update colormap"""
        self.current_colormap = colormap
        self.update_display()
        
    def toggle_crosshairs(self, state):
        """Toggle crosshairs display"""
        self.show_crosshairs = state == Qt.Checked
        self.update_display()
        
    def update_x_slice(self, value):
        """Update X slice"""
        self.current_slice_x = value
        self.x_label.setText(f"X: {value}")
        self.update_display()
        
    def update_y_slice(self, value):
        """Update Y slice"""
        self.current_slice_y = value
        self.y_label.setText(f"Y: {value}")
        self.update_display()
        
    def update_z_slice(self, value):
        """Update Z slice"""
        self.current_slice_z = value
        self.z_label.setText(f"Z: {value}")
        self.update_display()
    
    def update_display(self):
        """Update the display with consistent colormap scaling"""
        if self.volume is None:
            return
            
        self.figure.clear()
        
        # Ensure slice indices are within bounds
        self.current_slice_x = max(0, min(self.current_slice_x, self.volume.shape[0] - 1))
        self.current_slice_y = max(0, min(self.current_slice_y, self.volume.shape[1] - 1))
        self.current_slice_z = max(0, min(self.current_slice_z, self.volume.shape[2] - 1))
        
        # Get global min/max for consistent colormap scaling across all slices
        vmin = self.volume.min()
        vmax = self.volume.max()
        
        # Create subplots
        ax1 = self.figure.add_subplot(131)
        ax2 = self.figure.add_subplot(132)
        ax3 = self.figure.add_subplot(133)
        
        # YZ plane (sagittal)
        slice_yz = self.volume[self.current_slice_x, :, :].T
        im1 = ax1.imshow(slice_yz, cmap=self.current_colormap, origin='lower', vmin=vmin, vmax=vmax)
        ax1.set_title(f'YZ Plane (X={self.current_slice_x})')
        ax1.set_xlabel('Y')
        ax1.set_ylabel('Z')
        
        # XZ plane (coronal)  
        slice_xz = self.volume[:, self.current_slice_y, :].T
        im2 = ax2.imshow(slice_xz, cmap=self.current_colormap, origin='lower', vmin=vmin, vmax=vmax)
        ax2.set_title(f'XZ Plane (Y={self.current_slice_y})')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Z')
        
        # XY plane (axial)
        slice_xy = self.volume[:, :, self.current_slice_z].T
        im3 = ax3.imshow(slice_xy, cmap=self.current_colormap, origin='lower', vmin=vmin, vmax=vmax)
        ax3.set_title(f'XY Plane (Z={self.current_slice_z})')
        ax3.set_xlabel('X')
        ax3.set_ylabel('Y')
        
        # Add crosshairs
        if self.show_crosshairs:
            # YZ plane crosshairs
            ax1.axhline(self.current_slice_z, color='white', linestyle='--', alpha=0.7)
            ax1.axvline(self.current_slice_y, color='white', linestyle='--', alpha=0.7)
            
            # XZ plane crosshairs
            ax2.axhline(self.current_slice_z, color='white', linestyle='--', alpha=0.7)
            ax2.axvline(self.current_slice_x, color='white', linestyle='--', alpha=0.7)
            
            # XY plane crosshairs
            ax3.axhline(self.current_slice_y, color='white', linestyle='--', alpha=0.7)
            ax3.axvline(self.current_slice_x, color='white', linestyle='--', alpha=0.7)
        
        # Add colorbars with consistent scaling
        self.figure.colorbar(im1, ax=ax1, shrink=0.8)
        self.figure.colorbar(im2, ax=ax2, shrink=0.8)
        self.figure.colorbar(im3, ax=ax3, shrink=0.8)
        
        self.figure.tight_layout()
        self.canvas.draw()
    
    def on_mouse_move(self, event):
        """Handle mouse movement for cursor tracking"""
        if event.inaxes is not None and self.volume is not None:
            x, y = int(event.xdata or 0), int(event.ydata or 0)
            
            # Determine which subplot and get value
            if event.inaxes == self.figure.axes[0]:  # YZ plane
                if 0 <= y < self.volume.shape[2] and 0 <= x < self.volume.shape[1]:
                    value = self.volume[self.current_slice_x, x, y]
                    self.cursor_label.setText(f"YZ Plane | Cursor: ({x}, {y}) | Value: {value:.2f}")
            elif event.inaxes == self.figure.axes[1]:  # XZ plane
                if 0 <= y < self.volume.shape[2] and 0 <= x < self.volume.shape[0]:
                    value = self.volume[x, self.current_slice_y, y]
                    self.cursor_label.setText(f"XZ Plane | Cursor: ({x}, {y}) | Value: {value:.2f}")
            elif event.inaxes == self.figure.axes[2]:  # XY plane
                if 0 <= y < self.volume.shape[1] and 0 <= x < self.volume.shape[0]:
                    value = self.volume[x, y, self.current_slice_z]
                    self.cursor_label.setText(f"XY Plane | Cursor: ({x}, {y}) | Value: {value:.2f}")
    
    def on_scroll(self, event):
        """Handle scroll wheel for slice navigation"""
        if event.inaxes is not None and self.volume is not None:
            direction = 1 if event.step > 0 else -1
            
            # Determine which subplot and update appropriate slice with bounds checking
            if len(self.figure.axes) >= 3:
                if event.inaxes == self.figure.axes[0]:  # YZ plane - change X
                    new_x = max(0, min(self.volume.shape[0]-1, self.current_slice_x + direction))
                    if new_x != self.current_slice_x:
                        self.x_slider.setValue(new_x)
                elif event.inaxes == self.figure.axes[1]:  # XZ plane - change Y
                    new_y = max(0, min(self.volume.shape[1]-1, self.current_slice_y + direction))
                    if new_y != self.current_slice_y:
                        self.y_slider.setValue(new_y)
                elif event.inaxes == self.figure.axes[2]:  # XY plane - change Z
                    new_z = max(0, min(self.volume.shape[2]-1, self.current_slice_z + direction))
                    if new_z != self.current_slice_z:
                        self.z_slider.setValue(new_z)

class ExportWidget(QWidget):
    """Export reconstructed volume as DICOM"""
    
    def __init__(self):
        super().__init__()
        self.reconstructed_volume = None
        self.setup_ui()
        
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Export controls
        export_group = QGroupBox("Export to DICOM")
        export_layout = QGridLayout(export_group)
        
        # Scale factor
        export_layout.addWidget(QLabel("Scale Factor:"), 0, 0)
        self.scale_spin = QDoubleSpinBox()
        self.scale_spin.setRange(0.1, 10.0)
        self.scale_spin.setValue(1.0)
        self.scale_spin.setDecimals(2)
        self.scale_spin.setToolTip("Scale factor for intensity values")
        export_layout.addWidget(self.scale_spin, 0, 1)
        
        # Output path
        export_layout.addWidget(QLabel("Output Path:"), 1, 0)
        self.path_edit = QLabel("No path selected")
        self.path_edit.setStyleSheet("border: 1px solid gray; padding: 2px;")
        export_layout.addWidget(self.path_edit, 1, 1)
        
        self.browse_path_button = QPushButton("Browse...")
        self.browse_path_button.clicked.connect(self.browse_output_path)
        export_layout.addWidget(self.browse_path_button, 1, 2)
        
        # Export button
        self.export_button = QPushButton("Export DICOM")
        self.export_button.clicked.connect(self.export_dicom)
        self.export_button.setEnabled(False)
        export_layout.addWidget(self.export_button, 2, 0, 1, 3)
        
        layout.addWidget(export_group)
        
        # Status and info
        self.export_status = QLabel("Ready to export")
        layout.addWidget(self.export_status)
        
        # TPS template info
        if TPS_AVAILABLE:
            tps_info = QLabel("TPS Template: Available (EPID_12_t0.dcm)")
            tps_info.setStyleSheet("color: green;")
        else:
            tps_info = QLabel("TPS Template: Not Available")
            tps_info.setStyleSheet("color: red;")
        layout.addWidget(tps_info)
        
        layout.addStretch()
        
    def set_reconstructed_volume(self, volume):
        """Set reconstructed volume for export"""
        self.reconstructed_volume = volume
        self.update_export_status()
        
    def browse_output_path(self):
        """Browse for output file path"""
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save DICOM File", 
            "epid_reconstruction.dcm",
            "DICOM Files (*.dcm);;All Files (*)"
        )
        if file_path:
            self.path_edit.setText(file_path)
            self.output_path = file_path
            self.update_export_status()
    
    def update_export_status(self):
        """Update export button status"""
        can_export = (self.reconstructed_volume is not None and 
                     hasattr(self, 'output_path') and 
                     TPS_AVAILABLE)
        self.export_button.setEnabled(can_export)
        
        if can_export:
            self.export_status.setText("Ready to export")
        elif self.reconstructed_volume is None:
            self.export_status.setText("No reconstructed volume available")
        elif not hasattr(self, 'output_path'):
            self.export_status.setText("Please select output path")
        elif not TPS_AVAILABLE:
            self.export_status.setText("TPS template not available")
    
    def export_dicom(self):
        """Export reconstructed volume as DICOM"""
        if not hasattr(self, 'output_path') or self.reconstructed_volume is None:
            return
            
        self.export_status.setText("Exporting...")
        self.export_button.setEnabled(False)
        
        try:
            scale_factor = self.scale_spin.value()
            success = export_reconstructed_dicom(
                self.reconstructed_volume, 
                self.output_path, 
                scale_factor
            )
            
            if success:
                self.export_status.setText(f"Export successful: {Path(self.output_path).name}")
                QMessageBox.information(self, "Export Complete", 
                                      f"DICOM file exported successfully!\n\n{self.output_path}")
            else:
                self.export_status.setText("Export failed")
                QMessageBox.critical(self, "Export Error", "Failed to export DICOM file.")
                
        except Exception as e:
            self.export_status.setText(f"Export error: {str(e)}")
            QMessageBox.critical(self, "Export Error", f"An error occurred:\n{str(e)}")
        
        finally:
            self.export_button.setEnabled(True)

class EPIDReconstructionGUI(QMainWindow):
    """Main GUI application"""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("EPID Dose Reconstruction - Advanced GUI")
        self.setGeometry(100, 100, 1400, 900)
        
        self.current_volume = None
        self.worker = None
        
        self.setup_ui()
        self.setup_menu()
        self.setup_status_bar()
        
    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QHBoxLayout(central_widget)
        
        # Left panel - controls
        left_panel = QWidget()
        left_panel.setMaximumWidth(350)
        left_panel.setMinimumWidth(300)
        left_layout = QVBoxLayout(left_panel)
        
        # Data selection
        data_group = QGroupBox("Data Selection")
        data_layout = QVBoxLayout(data_group)
        
        self.path_label = QLabel("No data selected")
        data_layout.addWidget(self.path_label)
        
        self.browse_button = QPushButton("Browse Data Folder")
        self.browse_button.clicked.connect(self.browse_data_folder)
        data_layout.addWidget(self.browse_button)
        
        left_layout.addWidget(data_group)
        
        # Reconstruction parameters
        params_group = QGroupBox("Reconstruction Parameters")
        params_layout = QGridLayout(params_group)
        
        params_layout.addWidget(QLabel("Mode:"), 0, 0)
        self.recon_mode_combo = QComboBox()
        self.recon_mode_combo.addItems(['Quick (Every 4th)', 'Full (All Slices)'])
        self.recon_mode_combo.currentTextChanged.connect(self.update_mode_params)
        params_layout.addWidget(self.recon_mode_combo, 0, 1)
        
        params_layout.addWidget(QLabel("Image Size:"), 1, 0)
        self.image_size_spin = QSpinBox()
        self.image_size_spin.setRange(25, 200)
        self.image_size_spin.setValue(100)  # Default to 100 as standard size
        params_layout.addWidget(self.image_size_spin, 1, 1)
        
        params_layout.addWidget(QLabel("Chunk Size:"), 2, 0)
        self.chunk_size_spin = QSpinBox()
        self.chunk_size_spin.setRange(10, 100)
        self.chunk_size_spin.setValue(25)
        params_layout.addWidget(self.chunk_size_spin, 2, 1)
        
        # Rotation correction toggle
        params_layout.addWidget(QLabel("Rotation Correction:"), 3, 0)
        self.rotation_check = QCheckBox("Enable")
        self.rotation_check.setChecked(False)  # Unchecked by default
        self.rotation_check.setToolTip("Apply rotation correction when threshold exceeded")
        params_layout.addWidget(self.rotation_check, 3, 1)
        
        left_layout.addWidget(params_group)
        
        # Reconstruction controls
        recon_group = QGroupBox("Reconstruction")
        recon_layout = QVBoxLayout(recon_group)
        
        self.reconstruct_button = QPushButton("Start Reconstruction")
        self.reconstruct_button.clicked.connect(self.start_reconstruction)
        self.reconstruct_button.setEnabled(False)
        recon_layout.addWidget(self.reconstruct_button)
        
        self.progress_bar = QProgressBar()
        recon_layout.addWidget(self.progress_bar)
        
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(80)
        self.status_text.setReadOnly(True)
        recon_layout.addWidget(self.status_text)
        
        left_layout.addWidget(recon_group)
        
        # Add export widget to left panel
        self.export_widget = ExportWidget()
        left_layout.addWidget(self.export_widget)
        
        left_layout.addStretch()
        main_layout.addWidget(left_panel)
        
        # Right panel - volume viewer only
        self.volume_viewer = InteractiveVolumeViewer()
        main_layout.addWidget(self.volume_viewer)
        
    def setup_menu(self):
        """Setup menu bar"""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu('File')
        
        # Add actions
        file_menu.addAction('Open Data Folder', self.browse_data_folder)
        file_menu.addSeparator()
        file_menu.addAction('Exit', self.close)
        
        # View menu
        view_menu = menubar.addMenu('View')
        view_menu.addAction('Reset View', self.reset_view)
        
        # Help menu
        help_menu = menubar.addMenu('Help')
        help_menu.addAction('About', self.show_about)
        
    def setup_status_bar(self):
        """Setup status bar"""
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")
        
    def update_mode_params(self, mode_text):
        """Update parameters based on reconstruction mode"""
        if "Quick" in mode_text:
            # Quick mode defaults - keep 100 as standard size
            self.image_size_spin.setValue(100)
            self.chunk_size_spin.setValue(25)
        else:
            # Full mode defaults
            self.image_size_spin.setValue(100)
            self.chunk_size_spin.setValue(15)  # Smaller chunks for larger volumes
    
    def browse_data_folder(self):
        """Browse for data folder"""
        folder = QFileDialog.getExistingDirectory(self, "Select EPID Data Folder")
        if folder:
            self.path_label.setText(f"Selected: {Path(folder).name}")
            self.data_path = folder
            self.reconstruct_button.setEnabled(True)
            self.status_bar.showMessage(f"Data folder selected: {folder}")
            
    def start_reconstruction(self):
        """Start reconstruction in worker thread"""
        if not hasattr(self, 'data_path'):
            return
            
        # Disable controls
        self.reconstruct_button.setEnabled(False)
        self.progress_bar.setValue(0)
        self.status_text.clear()
        
        # Setup worker
        is_quick_mode = "Quick" in self.recon_mode_combo.currentText()
        params = {
            'image_size': self.image_size_spin.value(),
            'chunk_size': self.chunk_size_spin.value(),
            'enable_rotation': self.rotation_check.isChecked(),
            'quick_mode': is_quick_mode,
            'every_nth': 4 if is_quick_mode else 1
        }
        
        self.worker = ReconstructionWorker(self.data_path, params)
        self.worker.progress_update.connect(self.progress_bar.setValue)
        self.worker.status_update.connect(self.update_status)
        self.worker.result_ready.connect(self.reconstruction_complete)
        self.worker.error_occurred.connect(self.reconstruction_error)
        self.worker.start()
        
    def update_status(self, message):
        """Update status text"""
        self.status_text.append(message)
        self.status_bar.showMessage(message)
        
    def reconstruction_complete(self, result):
        """Handle reconstruction completion"""
        self.current_volume = result['reconstructed_volume']
        
        # Update viewers
        self.volume_viewer.set_volume(self.current_volume)
        self.export_widget.set_reconstructed_volume(self.current_volume)
        
        # Show summary
        params = result['parameters']
        summary = f"""
Reconstruction Complete!
Volume shape: {self.current_volume.shape}
Value range: {self.current_volume.min():.2f} to {self.current_volume.max():.2f}
Projections: {params['n_projections']}
SOD: {params['SOD']:.1f} mm
SDD: {params['SDD']:.1f} mm
"""
        self.update_status(summary)
        
        # Re-enable controls
        self.reconstruct_button.setEnabled(True)
        self.status_bar.showMessage("Reconstruction complete")
        
    def reconstruction_error(self, error_message):
        """Handle reconstruction error"""
        self.update_status(f"Error: {error_message}")
        self.reconstruct_button.setEnabled(True)
        self.status_bar.showMessage("Reconstruction failed")
        
        QMessageBox.critical(self, "Reconstruction Error", f"An error occurred:\n{error_message}")
        
    def reset_view(self):
        """Reset volume viewer"""
        if self.current_volume is not None:
            self.volume_viewer.set_volume(self.current_volume)
            
    def show_about(self):
        """Show about dialog"""
        about_text = """
EPID Dose Reconstruction GUI

Advanced PySide6 application for EPID dose reconstruction using cone beam FDK algorithm.

Features:
• Interactive 3D volume visualization
• Multiple color mapping options  
• Real-time dose verification
• Optimized reconstruction pipeline
• Mouse cursor tracking and scroll navigation

Based on corrected FDK reconstruction algorithm.
        """
        QMessageBox.about(self, "About", about_text)

def main():
    """Main application entry point"""
    app = QApplication(sys.argv)
    
    # Set application properties
    app.setApplicationName("EPID Reconstruction GUI")
    app.setApplicationVersion("1.0")
    
    # Check dependencies
    missing_deps = []
    if not PYDICOM_AVAILABLE:
        missing_deps.append("pydicom")
    
    if missing_deps:
        QMessageBox.critical(None, "Missing Dependencies", 
                           f"Missing required packages: {', '.join(missing_deps)}\n\n"
                           "Please install with: pip install " + " ".join(missing_deps))
        return 1
    
    # Create and show main window
    window = EPIDReconstructionGUI()
    window.show()
    
    return app.exec()

if __name__ == "__main__":
    sys.exit(main())