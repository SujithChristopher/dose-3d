"""Interactive 3-plane slice viewer for a reconstructed 3D volume."""
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
                                QComboBox, QCheckBox, QSlider, QLabel, QGroupBox)
from PySide6.QtCore import Qt

from matplotlib.figure import Figure

from ..optional_deps import FigureCanvas, TPS_AVAILABLE, calculate_volume_statistics


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
