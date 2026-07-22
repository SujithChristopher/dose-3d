"""Panel for exporting a reconstructed volume as a TPS-reference-templated DICOM."""
from pathlib import Path

from PySide6.QtWidgets import (QWidget, QVBoxLayout, QGridLayout, QLabel, QGroupBox,
                                QDoubleSpinBox, QPushButton, QFileDialog, QMessageBox)

from ..optional_deps import TPS_AVAILABLE, export_reconstructed_dicom


class ExportWidget(QWidget):
    """Export reconstructed volume as DICOM"""

    def __init__(self):
        super().__init__()
        self.reconstructed_volume = None
        self.geometry = None
        self.dose_units = 'RELATIVE'
        self.unit_scale = 1.0
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

    def set_reconstructed_volume(self, volume, geometry=None,
                                 dose_units='RELATIVE', unit_scale=1.0):
        """Set reconstructed volume for export.

        geometry is the reconstruction 'parameters' dict (for voxel_size_mm /
        origin_mm); dose_units is 'GY' once a calibration has been applied, and
        unit_scale converts the volume into that unit (0.01 for cGy -> Gy).
        """
        self.reconstructed_volume = volume
        self.geometry = geometry or {}
        self.dose_units = dose_units
        self.unit_scale = unit_scale
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
            scale_factor = self.scale_spin.value() * self.unit_scale
            success = export_reconstructed_dicom(
                self.reconstructed_volume,
                self.output_path,
                scale_factor,
                voxel_size_mm=self.geometry.get('voxel_size_mm'),
                origin_mm=self.geometry.get('origin_mm'),
                dose_units=self.dose_units
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
