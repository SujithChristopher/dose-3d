"""Load two DICOM volumes, sum them, view all three, and export the sum."""
from pathlib import Path

import numpy as np
from PySide6.QtWidgets import (QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QLabel,
                                QGroupBox, QComboBox, QCheckBox, QSlider, QPushButton,
                                QFileDialog, QMessageBox, QDoubleSpinBox)
from PySide6.QtCore import Qt

from matplotlib.figure import Figure

from ..optional_deps import FigureCanvas, dcmread, TPS_AVAILABLE, export_reconstructed_dicom


class DicomAdditionViewer(QWidget):
    """Viewer for showing two DICOM volumes and their sum"""

    def __init__(self):
        super().__init__()
        self.volume1 = None
        self.volume2 = None
        self.summed_volume = None
        self.current_slice_x = 50
        self.current_slice_y = 50
        self.current_slice_z = 50
        self.current_colormap = 'hot'
        self.show_crosshairs = True
        self.transpose_mode = False

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

        # Matplotlib figure - optimized for 3x3 grid with larger subplots
        self.figure = Figure(figsize=(15, 12), dpi=80)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        layout.addWidget(self.canvas)

        # Cursor position label
        self.cursor_label = QLabel("Cursor: (0, 0) | Value: 0.00")
        layout.addWidget(self.cursor_label)

    def set_volumes(self, volume1, volume2):
        """Set the two volumes to display and compute sum"""
        self.volume1 = volume1
        self.volume2 = volume2

        if volume1 is not None and volume2 is not None:
            # Ensure volumes are same shape
            if volume1.shape != volume2.shape:
                print(f"Warning: Volume shapes don't match: {volume1.shape} vs {volume2.shape}")
                # Resize to minimum common shape
                min_shape = tuple(min(s1, s2) for s1, s2 in zip(volume1.shape, volume2.shape))
                self.volume1 = volume1[:min_shape[0], :min_shape[1], :min_shape[2]]
                self.volume2 = volume2[:min_shape[0], :min_shape[1], :min_shape[2]]

            # Compute sum
            self.summed_volume = self.volume1 + self.volume2

            # Set slider ranges based on actual volume dimensions
            max_x = self.volume1.shape[0] - 1
            max_y = self.volume1.shape[1] - 1
            max_z = self.volume1.shape[2] - 1

            self.x_slider.setRange(0, max_x)
            self.y_slider.setRange(0, max_y)
            self.z_slider.setRange(0, max_z)

            # Set current slices to center
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

            self.update_display()

    def set_transpose_mode(self, transpose_mode):
        """Set transpose mode and update display"""
        self.transpose_mode = transpose_mode
        self.update_display()

    def get_display_volumes(self):
        """Get volumes for display, applying transpose if enabled"""
        if self.volume1 is None or self.volume2 is None or self.summed_volume is None:
            return None, None, None

        if self.transpose_mode:
            vol1 = np.transpose(self.volume1, (2, 1, 0))
            vol2 = np.transpose(self.volume2, (2, 1, 0))
            vol_sum = np.transpose(self.summed_volume, (2, 1, 0))
        else:
            vol1 = self.volume1
            vol2 = self.volume2
            vol_sum = self.summed_volume

        return vol1, vol2, vol_sum

    def get_summed_volume(self, transpose=False):
        """Get the summed volume, optionally transposed"""
        if self.summed_volume is None:
            return None

        if transpose:
            return np.transpose(self.summed_volume, (2, 1, 0))
        else:
            return self.summed_volume.copy()

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
        """Update the display showing all three volumes"""
        if self.volume1 is None or self.volume2 is None or self.summed_volume is None:
            return

        self.figure.clear()

        # Get display volumes (potentially transposed)
        vol1, vol2, vol_sum = self.get_display_volumes()
        if vol1 is None:
            return

        # Ensure slice indices are within bounds
        self.current_slice_x = max(0, min(self.current_slice_x, vol1.shape[0] - 1))
        self.current_slice_y = max(0, min(self.current_slice_y, vol1.shape[1] - 1))
        self.current_slice_z = max(0, min(self.current_slice_z, vol1.shape[2] - 1))

        # Get global min/max for consistent colormap scaling
        all_volumes = [vol1, vol2, vol_sum]
        vmin = min(vol.min() for vol in all_volumes)
        vmax = max(vol.max() for vol in all_volumes)

        # Create 3x3 subplot grid with improved spacing
        for row, (volume, title) in enumerate([(vol1, "DICOM 1"),
                                              (vol2, "DICOM 2"),
                                              (vol_sum, "Sum")]):

            # YZ plane (sagittal)
            ax1 = self.figure.add_subplot(3, 3, row*3 + 1)
            slice_yz = volume[self.current_slice_x, :, :].T
            im1 = ax1.imshow(slice_yz, cmap=self.current_colormap, origin='lower', vmin=vmin, vmax=vmax)
            ax1.set_title(f'{title} - YZ (X={self.current_slice_x})', fontsize=10, pad=8)
            ax1.set_xlabel('Y', fontsize=9)
            ax1.set_ylabel('Z', fontsize=9)
            ax1.tick_params(axis='both', which='major', labelsize=8)

            # XZ plane (coronal)
            ax2 = self.figure.add_subplot(3, 3, row*3 + 2)
            slice_xz = volume[:, self.current_slice_y, :].T
            im2 = ax2.imshow(slice_xz, cmap=self.current_colormap, origin='lower', vmin=vmin, vmax=vmax)
            ax2.set_title(f'{title} - XZ (Y={self.current_slice_y})', fontsize=10, pad=8)
            ax2.set_xlabel('X', fontsize=9)
            ax2.set_ylabel('Z', fontsize=9)
            ax2.tick_params(axis='both', which='major', labelsize=8)

            # XY plane (axial)
            ax3 = self.figure.add_subplot(3, 3, row*3 + 3)
            slice_xy = volume[:, :, self.current_slice_z].T
            im3 = ax3.imshow(slice_xy, cmap=self.current_colormap, origin='lower', vmin=vmin, vmax=vmax)
            ax3.set_title(f'{title} - XY (Z={self.current_slice_z})', fontsize=10, pad=8)
            ax3.set_xlabel('X', fontsize=9)
            ax3.set_ylabel('Y', fontsize=9)
            ax3.tick_params(axis='both', which='major', labelsize=8)

            # Add crosshairs
            if self.show_crosshairs:
                # YZ plane crosshairs
                ax1.axhline(self.current_slice_z, color='white', linestyle='--', alpha=0.7, linewidth=1)
                ax1.axvline(self.current_slice_y, color='white', linestyle='--', alpha=0.7, linewidth=1)

                # XZ plane crosshairs
                ax2.axhline(self.current_slice_z, color='white', linestyle='--', alpha=0.7, linewidth=1)
                ax2.axvline(self.current_slice_x, color='white', linestyle='--', alpha=0.7, linewidth=1)

                # XY plane crosshairs
                ax3.axhline(self.current_slice_y, color='white', linestyle='--', alpha=0.7, linewidth=1)
                ax3.axvline(self.current_slice_x, color='white', linestyle='--', alpha=0.7, linewidth=1)

        # Improved layout with better spacing
        self.figure.subplots_adjust(left=0.05, bottom=0.08, right=0.95, top=0.92,
                                   wspace=0.25, hspace=0.35)
        self.canvas.draw()

    def on_mouse_move(self, event):
        """Handle mouse movement for cursor tracking"""
        if event.inaxes is not None and self.summed_volume is not None:
            x, y = int(event.xdata or 0), int(event.ydata or 0)

            # Determine which volume and plane
            subplot_idx = None
            for i, ax in enumerate(self.figure.axes):
                if event.inaxes == ax:
                    subplot_idx = i
                    break

            if subplot_idx is not None:
                row = subplot_idx // 3  # 0=vol1, 1=vol2, 2=sum
                col = subplot_idx % 3   # 0=YZ, 1=XZ, 2=XY

                # Get display volumes (potentially transposed)
                vol1, vol2, vol_sum = self.get_display_volumes()
                if vol1 is None:
                    return

                volume = [vol1, vol2, vol_sum][row]
                volume_name = ["DICOM 1", "DICOM 2", "Sum"][row]
                plane_name = ["YZ", "XZ", "XY"][col]

                # Get value based on plane
                if col == 0:  # YZ plane
                    if 0 <= y < volume.shape[2] and 0 <= x < volume.shape[1]:
                        value = volume[self.current_slice_x, x, y]
                elif col == 1:  # XZ plane
                    if 0 <= y < volume.shape[2] and 0 <= x < volume.shape[0]:
                        value = volume[x, self.current_slice_y, y]
                elif col == 2:  # XY plane
                    if 0 <= y < volume.shape[1] and 0 <= x < volume.shape[0]:
                        value = volume[x, y, self.current_slice_z]
                else:
                    value = 0

                self.cursor_label.setText(f"{volume_name} {plane_name} | Cursor: ({x}, {y}) | Value: {value:.2f}")

    def on_scroll(self, event):
        """Handle scroll wheel for slice navigation"""
        if event.inaxes is not None and self.summed_volume is not None:
            direction = 1 if event.step > 0 else -1

            # Get display volumes to check bounds
            vol1, vol2, vol_sum = self.get_display_volumes()
            if vol1 is None:
                return

            # Determine which subplot
            subplot_idx = None
            for i, ax in enumerate(self.figure.axes):
                if event.inaxes == ax:
                    subplot_idx = i
                    break

            if subplot_idx is not None:
                col = subplot_idx % 3  # 0=YZ, 1=XZ, 2=XY

                if col == 0:  # YZ plane - change X
                    new_x = max(0, min(vol1.shape[0]-1, self.current_slice_x + direction))
                    if new_x != self.current_slice_x:
                        self.x_slider.setValue(new_x)
                elif col == 1:  # XZ plane - change Y
                    new_y = max(0, min(vol1.shape[1]-1, self.current_slice_y + direction))
                    if new_y != self.current_slice_y:
                        self.y_slider.setValue(new_y)
                elif col == 2:  # XY plane - change Z
                    new_z = max(0, min(vol1.shape[2]-1, self.current_slice_z + direction))
                    if new_z != self.current_slice_z:
                        self.z_slider.setValue(new_z)


class DicomAdditionWidget(QWidget):
    """Widget for DICOM addition functionality"""

    def __init__(self):
        super().__init__()
        self.dicom1_data = None
        self.dicom2_data = None
        self.dicom1_path = None
        self.dicom2_path = None

        self.setup_ui()

    def setup_ui(self):
        main_layout = QHBoxLayout(self)

        # Left panel - controls
        left_panel = QWidget()
        left_panel.setMaximumWidth(350)
        left_panel.setMinimumWidth(300)
        left_layout = QVBoxLayout(left_panel)

        # DICOM Import section
        import_group = QGroupBox("DICOM Import")
        import_layout = QGridLayout(import_group)

        # DICOM 1
        import_layout.addWidget(QLabel("DICOM 1:"), 0, 0)
        self.dicom1_label = QLabel("No file selected")
        self.dicom1_label.setStyleSheet("border: 1px solid gray; padding: 2px;")
        import_layout.addWidget(self.dicom1_label, 0, 1)

        self.browse1_button = QPushButton("Browse...")
        self.browse1_button.clicked.connect(self.browse_dicom1)
        import_layout.addWidget(self.browse1_button, 0, 2)

        # DICOM 2
        import_layout.addWidget(QLabel("DICOM 2:"), 1, 0)
        self.dicom2_label = QLabel("No file selected")
        self.dicom2_label.setStyleSheet("border: 1px solid gray; padding: 2px;")
        import_layout.addWidget(self.dicom2_label, 1, 1)

        self.browse2_button = QPushButton("Browse...")
        self.browse2_button.clicked.connect(self.browse_dicom2)
        import_layout.addWidget(self.browse2_button, 1, 2)

        left_layout.addWidget(import_group)

        # Processing section
        process_group = QGroupBox("Processing Options")
        process_layout = QVBoxLayout(process_group)

        # Transpose option
        self.transpose_check = QCheckBox("Apply Transpose (2,1,0)")
        self.transpose_check.setToolTip("Apply np.transpose(volume, (2,1,0)) for display and export")
        self.transpose_check.stateChanged.connect(self.on_transpose_changed)
        process_layout.addWidget(self.transpose_check)

        # Add volumes button
        self.add_button = QPushButton("Add Volumes")
        self.add_button.clicked.connect(self.add_volumes)
        self.add_button.setEnabled(False)
        process_layout.addWidget(self.add_button)

        left_layout.addWidget(process_group)

        # Export section
        export_group = QGroupBox("Export")
        export_layout = QGridLayout(export_group)

        # Scale factor
        export_layout.addWidget(QLabel("Scale Factor:"), 0, 0)
        self.scale_spin = QDoubleSpinBox()
        self.scale_spin.setRange(0.1, 10.0)
        self.scale_spin.setValue(1.0)
        self.scale_spin.setDecimals(2)
        export_layout.addWidget(self.scale_spin, 0, 1)

        # Output path
        export_layout.addWidget(QLabel("Output Path:"), 1, 0)
        self.output_path_label = QLabel("No path selected")
        self.output_path_label.setStyleSheet("border: 1px solid gray; padding: 2px;")
        export_layout.addWidget(self.output_path_label, 1, 1)

        self.browse_output_button = QPushButton("Browse...")
        self.browse_output_button.clicked.connect(self.browse_output_path)
        export_layout.addWidget(self.browse_output_button, 1, 2)

        # Export button
        self.export_button = QPushButton("Export Summed DICOM")
        self.export_button.clicked.connect(self.export_summed_dicom)
        self.export_button.setEnabled(False)
        export_layout.addWidget(self.export_button, 2, 0, 1, 3)

        left_layout.addWidget(export_group)

        # Status
        self.status_label = QLabel("Import two DICOM files to begin")
        left_layout.addWidget(self.status_label)

        left_layout.addStretch()
        main_layout.addWidget(left_panel)

        # Right panel - viewer
        self.viewer = DicomAdditionViewer()
        main_layout.addWidget(self.viewer)

    def browse_dicom1(self):
        """Browse for first DICOM file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select DICOM 1",
            "",
            "DICOM Files (*.dcm);;All Files (*)"
        )
        if file_path:
            self.dicom1_path = file_path
            self.dicom1_label.setText(Path(file_path).name)
            self.load_dicom1()

    def browse_dicom2(self):
        """Browse for second DICOM file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select DICOM 2",
            "",
            "DICOM Files (*.dcm);;All Files (*)"
        )
        if file_path:
            self.dicom2_path = file_path
            self.dicom2_label.setText(Path(file_path).name)
            self.load_dicom2()

    def load_dicom1(self):
        """Load first DICOM file"""
        try:
            if self.dicom1_path:
                dcm = dcmread(self.dicom1_path)
                self.dicom1_data = dcm.pixel_array.astype(np.float32)
                self.status_label.setText(f"DICOM 1 loaded: {self.dicom1_data.shape}")
                self.check_ready_to_add()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load DICOM 1: {str(e)}")

    def load_dicom2(self):
        """Load second DICOM file"""
        try:
            if self.dicom2_path:
                dcm = dcmread(self.dicom2_path)
                self.dicom2_data = dcm.pixel_array.astype(np.float32)
                self.status_label.setText(f"DICOM 2 loaded: {self.dicom2_data.shape}")
                self.check_ready_to_add()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load DICOM 2: {str(e)}")

    def check_ready_to_add(self):
        """Check if both DICOMs are loaded and enable add button"""
        if self.dicom1_data is not None and self.dicom2_data is not None:
            self.add_button.setEnabled(True)
            self.status_label.setText("Ready to add volumes")

    def add_volumes(self):
        """Add the two volumes and update viewer"""
        if self.dicom1_data is not None and self.dicom2_data is not None:
            self.viewer.set_volumes(self.dicom1_data, self.dicom2_data)
            self.viewer.set_transpose_mode(self.transpose_check.isChecked())
            self.export_button.setEnabled(hasattr(self, 'output_path'))
            self.status_label.setText("Volumes added successfully")

    def on_transpose_changed(self, state):
        """Handle transpose checkbox change"""
        if self.viewer.summed_volume is not None:
            self.viewer.set_transpose_mode(self.transpose_check.isChecked())
            self.status_label.setText(f"Transpose {'enabled' if self.transpose_check.isChecked() else 'disabled'}")

    def browse_output_path(self):
        """Browse for output file path"""
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Summed DICOM",
            "summed_volume.dcm",
            "DICOM Files (*.dcm);;All Files (*)"
        )
        if file_path:
            self.output_path = file_path
            self.output_path_label.setText(Path(file_path).name)
            self.export_button.setEnabled(self.viewer.summed_volume is not None)

    def export_summed_dicom(self):
        """Export the summed volume as DICOM"""
        if not hasattr(self, 'output_path') or self.viewer.summed_volume is None:
            return

        try:
            # Get summed volume with optional transpose
            summed_volume = self.viewer.get_summed_volume(transpose=self.transpose_check.isChecked())
            scale_factor = self.scale_spin.value()

            if TPS_AVAILABLE:
                success = export_reconstructed_dicom(summed_volume, self.output_path, scale_factor)
                if success:
                    self.status_label.setText(f"Export successful: {Path(self.output_path).name}")
                    QMessageBox.information(self, "Export Complete",
                                          f"Summed DICOM exported successfully!\n\n{self.output_path}")
                else:
                    self.status_label.setText("Export failed")
                    QMessageBox.critical(self, "Export Error", "Failed to export DICOM file.")
            else:
                QMessageBox.critical(self, "Export Error", "TPS utilities not available for export.")

        except Exception as e:
            self.status_label.setText(f"Export error: {str(e)}")
            QMessageBox.critical(self, "Export Error", f"An error occurred:\n{str(e)}")
