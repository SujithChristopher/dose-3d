"""Background worker that batch-reconstructs dose-labeled EPID acquisitions for calibration."""
import numpy as np
from PySide6.QtCore import QThread, Signal

from .reconstruction_worker import ReconstructionWorker


class BatchCalibrationWorker(QThread):
    """Worker thread that reconstructs a batch of dose-labeled EPID acquisitions
    and reports the center-voxel pixel value for each, for pixel<->dose calibration."""
    progress_update = Signal(int)      # overall 0-100 across all acquisitions
    status_update = Signal(str)
    acquisition_done = Signal(str, float, float)   # folder_label, dose_cGy, pixel_value
    acquisition_failed = Signal(str, str)          # folder_label, error message
    finished_all = Signal(list)        # list of {'folder', 'dose_cGy', 'pixel_value'}

    def __init__(self, acquisitions, reconstruction_params):
        """acquisitions: list of (dcm_folder_path, folder_label, dose_cGy)"""
        super().__init__()
        self.acquisitions = acquisitions
        self.params = reconstruction_params

    def run(self):
        results = []
        n = len(self.acquisitions)
        for i, (path, label, dose) in enumerate(self.acquisitions):
            try:
                self.status_update.emit(f"Reconstructing {label} ({dose} cGy)...")

                worker = ReconstructionWorker(path, self.params)
                worker.progress_update.connect(
                    lambda p, i=i: self.progress_update.emit(int(100 * (i + p / 100) / n)))
                captured = {}
                worker.result_ready.connect(lambda r: captured.update(result=r))
                worker.error_occurred.connect(lambda e: captured.update(error=e))
                worker.run()  # synchronous call - reuses full FDK pipeline, no extra thread

                if 'error' in captured:
                    raise RuntimeError(captured['error'])

                vol = captured['result']['reconstructed_volume']
                cz, cy, cx = (s // 2 for s in vol.shape)
                patch = vol[cz - 1:cz + 2, cy - 1:cy + 2, cx - 1:cx + 2]
                pixel_value = float(np.mean(patch))

                results.append({'folder': label, 'dose_cGy': dose, 'pixel_value': pixel_value})
                self.acquisition_done.emit(label, dose, pixel_value)

            except Exception as e:
                self.acquisition_failed.emit(label, str(e))

            self.progress_update.emit(int(100 * (i + 1) / n))

        self.finished_all.emit(results)
