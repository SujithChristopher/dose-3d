"""Background worker that reconstructs a single EPID acquisition folder into a 3D dose volume."""
import os

import numpy as np
from PySide6.QtCore import QThread, Signal

from ..algorithms.fdk import OptimizedFDKReconstructor
from ..optional_deps import dcmread, cv2, ThreadPoolExecutor, PI


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
            rec_image = np.transpose(rec_image, (2,1,0))

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
