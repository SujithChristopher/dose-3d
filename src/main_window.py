"""Main application window: tabs for reconstruction, DICOM addition, and dose calibration."""
from pathlib import Path

from PySide6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
                                QTabWidget, QGroupBox, QLabel, QPushButton, QComboBox,
                                QSpinBox, QCheckBox, QProgressBar, QTextEdit, QFileDialog,
                                QStatusBar, QMessageBox)

from .workers.reconstruction_worker import ReconstructionWorker
from .calibration import load_calibration, apply_calibration
from .updater import check_for_update_manual
from .widgets.volume_viewer import InteractiveVolumeViewer
from .widgets.export_widget import ExportWidget
from .widgets.dicom_addition import DicomAdditionWidget
from .widgets.dose_calibration import DoseCalibrationWidget


class EPIDReconstructionGUI(QMainWindow):
    """Main GUI application"""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("EPID Dose Reconstruction - Advanced GUI")
        self.setGeometry(100, 100, 1400, 900)

        self.current_volume = None
        self.worker = None
        self.calibration = None

        self.setup_ui()
        self.setup_menu()
        self.setup_status_bar()

    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Main layout with tab widget
        main_layout = QVBoxLayout(central_widget)

        # Create tab widget
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)

        # Tab 1: EPID Reconstruction
        self.setup_reconstruction_tab()

        # Tab 2: DICOM Addition
        self.setup_dicom_addition_tab()

        # Tab 3: Dose Calibration
        self.setup_dose_calibration_tab()

    def setup_reconstruction_tab(self):
        """Setup the EPID reconstruction tab"""
        recon_tab = QWidget()
        self.tab_widget.addTab(recon_tab, "EPID Reconstruction")

        # Main layout for reconstruction tab
        main_layout = QHBoxLayout(recon_tab)

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

        # Dose calibration
        calib_group = QGroupBox("Dose Calibration (optional)")
        calib_layout = QVBoxLayout(calib_group)

        self.calibration_label = QLabel("No calibration loaded")
        self.calibration_label.setWordWrap(True)
        calib_layout.addWidget(self.calibration_label)

        self.load_calibration_button = QPushButton("Load Calibration (JSON)...")
        self.load_calibration_button.clicked.connect(self.load_calibration_file)
        calib_layout.addWidget(self.load_calibration_button)

        self.apply_calibration_check = QCheckBox("Apply Calibration to Reconstruction")
        self.apply_calibration_check.setEnabled(False)
        calib_layout.addWidget(self.apply_calibration_check)

        self.include_intercept_check = QCheckBox("Include Intercept Offset")
        self.include_intercept_check.setToolTip(
            "Off (default) = slope-only scaling: dose = pixel / slope.\n"
            "On = dose = (pixel - intercept) / slope.")
        self.include_intercept_check.setEnabled(False)
        calib_layout.addWidget(self.include_intercept_check)

        left_layout.addWidget(calib_group)

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

        # Right panel - volume viewer
        self.volume_viewer = InteractiveVolumeViewer()
        main_layout.addWidget(self.volume_viewer)

    def setup_dicom_addition_tab(self):
        """Setup the DICOM addition tab"""
        addition_tab = QWidget()
        self.tab_widget.addTab(addition_tab, "DICOM Addition")

        # Add the DICOM addition widget to the tab
        layout = QVBoxLayout(addition_tab)
        self.dicom_addition_widget = DicomAdditionWidget()
        layout.addWidget(self.dicom_addition_widget)

    def setup_dose_calibration_tab(self):
        """Setup the dose calibration tab"""
        calibration_tab = QWidget()
        self.tab_widget.addTab(calibration_tab, "Dose Calibration")

        layout = QVBoxLayout(calibration_tab)
        self.dose_calibration_widget = DoseCalibrationWidget()
        layout.addWidget(self.dose_calibration_widget)

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
        help_menu.addAction('Check for Updates...', self.check_for_updates)
        help_menu.addSeparator()
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

    def load_calibration_file(self):
        """Load a pixel<->dose calibration JSON exported from the Dose Calibration tab"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Calibration", "", "JSON Files (*.json);;All Files (*)")
        if not file_path:
            return

        try:
            self.calibration = load_calibration(file_path)
            fit = self.calibration['fit']
            self.calibration_label.setText(
                f"Loaded: {Path(file_path).name}\n"
                f"slope={fit['slope']:.4g}, intercept={fit['intercept']:.4g}, R^2={fit['r_squared']:.4f}")
            self.apply_calibration_check.setEnabled(True)
            self.apply_calibration_check.setChecked(True)
            self.include_intercept_check.setEnabled(True)
            self.status_bar.showMessage(f"Calibration loaded: {Path(file_path).name}")
        except Exception as e:
            self.calibration = None
            QMessageBox.critical(self, "Calibration Error", f"Failed to load calibration:\n{str(e)}")

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
        volume = result['reconstructed_volume']
        units = "raw pixel value"

        if self.calibration is not None and self.apply_calibration_check.isChecked():
            include_intercept = self.include_intercept_check.isChecked()
            volume = apply_calibration(volume, self.calibration, include_intercept=include_intercept)
            units = "cGy, calibrated" if include_intercept else "cGy, slope-only calibrated"

        self.current_volume = volume

        # Update viewers
        self.volume_viewer.set_volume(self.current_volume)
        self.export_widget.set_reconstructed_volume(self.current_volume)

        # Show summary
        params = result['parameters']
        summary = f"""
                    Reconstruction Complete! ({units})
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

    def check_for_updates(self):
        """Manually check for a new release (Help > Check for Updates...)"""
        check_for_update_manual(self)

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
