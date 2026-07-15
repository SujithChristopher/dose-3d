"""Batch-reconstruct dose-labeled EPID acquisitions and fit pixel<->dose linear calibration."""
import os
import re
import json
from pathlib import Path
from datetime import datetime

import numpy as np
from scipy.stats import linregress
from PySide6.QtWidgets import (QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QLabel,
                                QGroupBox, QPushButton, QSpinBox, QCheckBox, QFileDialog,
                                QTableWidget, QTableWidgetItem, QProgressBar, QTextEdit,
                                QMessageBox)
from PySide6.QtCore import Qt

from matplotlib.figure import Figure

from ..optional_deps import FigureCanvas
from ..workers.reconstruction_worker import ReconstructionWorker
from ..workers.calibration_worker import BatchCalibrationWorker


def find_dicom_directory(root):
    """Return (dir_path, file_count) for the first directory at or under root containing .dcm files."""
    if not os.path.isdir(root):
        return root, 0
    direct = [f for f in os.listdir(root) if f.lower().endswith('.dcm')]
    if direct:
        return root, len(direct)
    for dirpath, _dirnames, filenames in os.walk(root):
        dcm_files = [f for f in filenames if f.lower().endswith('.dcm')]
        if dcm_files:
            return dirpath, len(dcm_files)
    return root, 0


class DoseCalibrationWidget(QWidget):
    """Batch-reconstructs multiple dose-labeled EPID acquisitions and fits a
    linear pixel-value <-> dose relationship from their center voxels."""

    DOSE_PATTERN = re.compile(r'(\d+(?:\.\d+)?)\s*cGy', re.IGNORECASE)

    def __init__(self):
        super().__init__()
        self.dataset_root = None
        self.batch_worker = None
        self.results = []
        self.fit_result = None  # (slope, intercept, r_value)
        self.last_params = {}
        self.setup_ui()

    def setup_ui(self):
        main_layout = QHBoxLayout(self)

        # Left panel - controls
        left_panel = QWidget()
        left_panel.setMaximumWidth(400)
        left_panel.setMinimumWidth(350)
        left_layout = QVBoxLayout(left_panel)

        # Dataset selection
        data_group = QGroupBox("Dataset Selection")
        data_layout = QVBoxLayout(data_group)

        self.dataset_label = QLabel("No dataset folder selected")
        self.dataset_label.setWordWrap(True)
        data_layout.addWidget(self.dataset_label)

        browse_row = QHBoxLayout()
        self.browse_dataset_button = QPushButton("Browse Dataset Folder...")
        self.browse_dataset_button.clicked.connect(self.browse_dataset_folder)
        browse_row.addWidget(self.browse_dataset_button)

        self.scan_button = QPushButton("Scan Folder")
        self.scan_button.clicked.connect(self.scan_folder)
        self.scan_button.setEnabled(False)
        browse_row.addWidget(self.scan_button)
        data_layout.addLayout(browse_row)

        left_layout.addWidget(data_group)

        # Acquisitions table
        table_group = QGroupBox("Detected Acquisitions")
        table_layout = QVBoxLayout(table_group)
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(
            ["Include", "Folder", "Dose (cGy)", "# Files", "Status", "Pixel Value"])
        self.table.horizontalHeader().setStretchLastSection(True)
        table_layout.addWidget(self.table)
        left_layout.addWidget(table_group)

        # Reconstruction parameters
        params_group = QGroupBox("Reconstruction Parameters")
        params_layout = QGridLayout(params_group)

        params_layout.addWidget(QLabel("Image Size:"), 0, 0)
        self.image_size_spin = QSpinBox()
        self.image_size_spin.setRange(25, 200)
        self.image_size_spin.setValue(100)
        params_layout.addWidget(self.image_size_spin, 0, 1)

        params_layout.addWidget(QLabel("Every Nth Frame:"), 1, 0)
        self.every_nth_spin = QSpinBox()
        self.every_nth_spin.setRange(1, 20)
        self.every_nth_spin.setValue(4)
        params_layout.addWidget(self.every_nth_spin, 1, 1)

        params_layout.addWidget(QLabel("Rotation Correction:"), 2, 0)
        self.rotation_check = QCheckBox("Enable")
        self.rotation_check.setChecked(False)
        params_layout.addWidget(self.rotation_check, 2, 1)

        left_layout.addWidget(params_group)

        # Run / export controls
        run_group = QGroupBox("Batch Calibration")
        run_layout = QVBoxLayout(run_group)

        self.run_button = QPushButton("Run Batch Reconstruction")
        self.run_button.clicked.connect(self.run_batch_reconstruction)
        self.run_button.setEnabled(False)
        run_layout.addWidget(self.run_button)

        self.progress_bar = QProgressBar()
        run_layout.addWidget(self.progress_bar)

        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(80)
        self.status_text.setReadOnly(True)
        run_layout.addWidget(self.status_text)

        self.export_button = QPushButton("Export Calibration (JSON)")
        self.export_button.clicked.connect(self.export_calibration)
        self.export_button.setEnabled(False)
        run_layout.addWidget(self.export_button)

        left_layout.addWidget(run_group)

        self.status_label = QLabel("Select a dataset folder to begin")
        left_layout.addWidget(self.status_label)

        left_layout.addStretch()
        main_layout.addWidget(left_panel)

        # Right panel - fit plot
        self.figure = Figure(figsize=(5, 4))
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self._reset_plot()
        main_layout.addWidget(self.canvas)

    def _reset_plot(self):
        self.ax.clear()
        self.ax.set_xlabel("Dose (cGy)")
        self.ax.set_ylabel("Pixel Value (center 3x3x3 mean)")
        self.ax.set_title("Dose Calibration")
        self.ax.grid(True, alpha=0.3)
        self.canvas.draw()

    def browse_dataset_folder(self):
        default_dir = ""
        candidate = Path(__file__).resolve().parent.parent.parent / "dataset" / "10x10 FS"
        if candidate.exists():
            default_dir = str(candidate)
        folder = QFileDialog.getExistingDirectory(self, "Select Dataset Folder", default_dir)
        if folder:
            self.dataset_root = folder
            self.dataset_label.setText(folder)
            self.scan_button.setEnabled(True)
            self.status_label.setText("Dataset folder selected - click Scan Folder")

    def scan_folder(self):
        if not self.dataset_root:
            return

        self.table.setRowCount(0)
        self.results = []
        self.fit_result = None
        self.export_button.setEnabled(False)
        self._reset_plot()

        subfolders = sorted(
            d for d in os.listdir(self.dataset_root)
            if os.path.isdir(os.path.join(self.dataset_root, d))
        )

        for name in subfolders:
            full_path = os.path.join(self.dataset_root, name)
            dcm_dir, n_files = find_dicom_directory(full_path)
            match = self.DOSE_PATTERN.search(name)
            dose_value = float(match.group(1)) if match else None

            row = self.table.rowCount()
            self.table.insertRow(row)

            include_item = QTableWidgetItem()
            if n_files > 0:
                include_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
                include_item.setCheckState(Qt.Checked)
            else:
                include_item.setFlags(Qt.ItemIsEnabled)
                include_item.setCheckState(Qt.Unchecked)
            self.table.setItem(row, 0, include_item)

            folder_item = QTableWidgetItem(name)
            folder_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            folder_item.setData(Qt.UserRole, dcm_dir)
            self.table.setItem(row, 1, folder_item)

            dose_item = QTableWidgetItem(f"{dose_value:.1f}" if dose_value is not None else "")
            dose_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable)
            self.table.setItem(row, 2, dose_item)

            files_item = QTableWidgetItem(str(n_files))
            files_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            self.table.setItem(row, 3, files_item)

            status_item = QTableWidgetItem("Ready" if n_files > 0 else "No DICOM files")
            status_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            self.table.setItem(row, 4, status_item)

            pixel_item = QTableWidgetItem("")
            pixel_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            self.table.setItem(row, 5, pixel_item)

        self.table.resizeColumnsToContents()
        self.run_button.setEnabled(self.table.rowCount() > 0)
        self.status_label.setText(f"Found {self.table.rowCount()} acquisition folder(s)")

    def run_batch_reconstruction(self):
        acquisitions = []
        row_by_label = {}
        for row in range(self.table.rowCount()):
            include_item = self.table.item(row, 0)
            if include_item.checkState() != Qt.Checked:
                continue
            folder_item = self.table.item(row, 1)
            dose_item = self.table.item(row, 2)
            files_item = self.table.item(row, 3)
            try:
                dose = float(dose_item.text())
            except ValueError:
                self.table.item(row, 4).setText("Invalid dose value")
                continue
            if int(files_item.text()) == 0:
                continue
            label = folder_item.text()
            dcm_dir = folder_item.data(Qt.UserRole)
            acquisitions.append((dcm_dir, label, dose))
            row_by_label[label] = row

        if not acquisitions:
            QMessageBox.warning(self, "No Acquisitions",
                                 "No valid acquisitions selected. Check dose values and file counts.")
            return

        self._row_by_label = row_by_label
        self.results = []
        self.fit_result = None
        self.run_button.setEnabled(False)
        self.export_button.setEnabled(False)
        self.progress_bar.setValue(0)
        self.status_text.clear()
        self._reset_plot()

        self.last_params = {
            'image_size': self.image_size_spin.value(),
            'chunk_size': 25,
            'enable_rotation': self.rotation_check.isChecked(),
            'every_nth': self.every_nth_spin.value(),
        }

        self.batch_worker = BatchCalibrationWorker(acquisitions, self.last_params, self.dataset_root)
        self.batch_worker.progress_update.connect(self.progress_bar.setValue)
        self.batch_worker.status_update.connect(self.append_status)
        self.batch_worker.acquisition_done.connect(self.on_acquisition_done)
        self.batch_worker.acquisition_failed.connect(self.on_acquisition_failed)
        self.batch_worker.finished_all.connect(self.on_batch_finished)
        self.batch_worker.start()

    def append_status(self, message):
        self.status_text.append(message)

    def on_acquisition_done(self, label, dose, pixel_value):
        row = self._row_by_label.get(label)
        if row is not None:
            self.table.item(row, 4).setText("Done")
            self.table.item(row, 5).setText(f"{pixel_value:.2f}")
        self.ax.scatter([dose], [pixel_value], color='tab:blue')
        self.canvas.draw()

    def on_acquisition_failed(self, label, error):
        row = self._row_by_label.get(label)
        if row is not None:
            self.table.item(row, 4).setText(f"Error: {error}")
        self.append_status(f"Failed: {label} - {error}")

    def on_batch_finished(self, results):
        self.results = results
        self.run_button.setEnabled(True)

        if len(results) < 2:
            self.status_label.setText("Need at least 2 successful acquisitions to fit a calibration.")
            return

        pixel_values = np.array([r['pixel_value'] for r in results])
        doses = np.array([r['dose_cGy'] for r in results])
        fit = linregress(doses, pixel_values)
        self.fit_result = (fit.slope, fit.intercept, fit.rvalue)

        x_line = np.linspace(doses.min(), doses.max(), 50)
        y_line = fit.slope * x_line + fit.intercept
        self.ax.plot(x_line, y_line, color='tab:red',
                     label=f"pixel = {fit.slope:.4g}*dose + {fit.intercept:.4g}\nR^2 = {fit.rvalue**2:.4f}")
        self.ax.legend(loc='best', fontsize=8)
        self.canvas.draw()

        self.export_button.setEnabled(True)
        self.status_label.setText(
            f"Fit complete: slope={fit.slope:.4g}, intercept={fit.intercept:.4g}, R^2={fit.rvalue**2:.4f}")

        if fit.slope <= 0 or fit.rvalue ** 2 < 0.8:
            QMessageBox.warning(
                self, "Calibration Warning",
                "Pixel value does not increase cleanly with dose (weak or negative correlation).\n"
                "Review the acquisitions and reconstruction parameters before trusting this calibration.")

    def export_calibration(self):
        if not self.results or self.fit_result is None:
            return

        default_name = f"calibration_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Calibration", default_name, "JSON Files (*.json);;All Files (*)")
        if not file_path:
            return

        slope, intercept, r_value = self.fit_result
        data = {
            'created': datetime.now().isoformat(),
            'dataset_folder': self.dataset_root,
            'reconstruction_params': self.last_params,
            'voxel_sampling': '3x3x3_mean',
            'points': self.results,
            'fit': {
                'slope': slope,
                'intercept': intercept,
                'r_squared': r_value ** 2,
                'formula': 'pixel_value = slope * dose_cGy + intercept'
            }
        }

        try:
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2)
            self.status_label.setText(f"Exported: {Path(file_path).name}")
            QMessageBox.information(self, "Export Complete", f"Calibration exported:\n{file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Failed to export calibration:\n{str(e)}")
