"""Background worker that reconstructs a single EPID acquisition folder into a 3D dose volume."""
import os

import numpy as np
from PySide6.QtCore import QThread, Signal

from ..algorithms.fdk import OptimizedFDKReconstructor
from ..optional_deps import dcmread, cv2, ThreadPoolExecutor, PI

# Backproject the ramp-filtered projections instead of the weighted ones.
#
# Deliberately off, and deliberately not exposed in the GUI. The filter itself is
# correct (it reproduces the analytic ramp response on a disc projection to within
# 3%), but it is the wrong operator for this data: an EPID differential frame is a
# broad fluence map whose content sits at ~DC, while the ramp has zero gain at DC
# and rises linearly to the 0.336 mm-pitch Nyquist, where only detector noise
# lives. Measured on real frames it takes projection SNR from ~86-180 down to
# ~4-6, and the reconstruction correlates only 0.47 with the unfiltered one.
#
# Unfiltered backprojection is a 1/r-blurred superposition of the delivered
# fluence, which is what the dose calibration was fitted against. Flip this to
# True only alongside detector binning and filter apodisation - see the v0.1.8
# regression analysis.
USE_RAMP_FILTER = False


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
        """Fast DICOM loading, rescaled to a single physical count scale.

        The linac rescales EPID gain as the cumulative signal grows (RescaleSlope
        steps 1 -> 2 -> 4 part-way through an acquisition to keep the running total
        inside 16 bits). Raw stored values from different frames are therefore NOT
        on the same scale, and differencing them across a gain step produces a
        full-frame negative jump. Applying RescaleSlope/Intercept here makes every
        frame comparable, which is what the differential and the dose calibration
        both assume.
        """
        def read_dicom_data(fname):
            try:
                dcm = dcmread(fname)
                slope = float(getattr(dcm, 'RescaleSlope', 1.0) or 1.0)
                intercept = float(getattr(dcm, 'RescaleIntercept', 0.0) or 0.0)
                counts = dcm.pixel_array.astype(np.float32) * slope + intercept
                return counts, float(dcm.GantryAngle)
            except Exception as e:
                return None, None

        # Frames go straight into one preallocated float32 stack, filled in
        # bounded blocks. Submitting every frame at once would keep the whole
        # acquisition alive twice over (in the futures and in the stack), which
        # is ~4.5 GB for a 400-frame 1190x1190 run.
        n_paths = len(path_list)
        block = max(max_workers * 4, 1)
        stack, angles, count = None, [], 0

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for start in range(0, n_paths, block):
                futures = [executor.submit(read_dicom_data, path)
                           for path in path_list[start:start + block]]
                for offset, future in enumerate(futures):
                    img, angle = future.result()
                    if img is not None:
                        if stack is None:
                            stack = np.empty((n_paths,) + img.shape, dtype=np.float32)
                        stack[count] = img
                        angles.append(angle)
                        count += 1
                    self.progress_update.emit(int(20 * (start + offset) / n_paths))

        if stack is None:
            raise RuntimeError(f"No readable DICOM frames in {n_paths} files")

        return stack[:count], np.array(angles)

    def process_differential_original(self, images, angles, threshold=10000, enable_rotation=True):
        """Process differential images with optional rotation correction.

        Operates on rescaled float counts from load_dicom_fast, so `threshold` now
        only guards against genuinely implausible frame-to-frame jumps (beam
        interruption / restart). It used to fire on every gain step, because
        differencing raw uint16 across a step wrapped around to ~2^15.
        """
        # Differenced in place: rescaling to float32 doubled the size of the frame
        # stack, so holding a separate output array as well would double peak memory
        # again on a full 1190x1190 acquisition. Only one frame is copied at a time.
        n_images = len(images)
        n_rows, n_cols = images[0].shape
        processed_angles = []

        prev = np.zeros((n_rows, n_cols), dtype=np.float32)

        for idx in range(n_images):
            curr = images[idx].copy()
            diff = curr - prev

            if np.max(diff) > threshold:
                if enable_rotation:
                    # Apply rotation correction (flip and abs) when threshold exceeded
                    images[idx] = cv2.flip(-np.abs(diff), 1)
                    processed_angles.append(angles[idx])
                else:
                    # No rotation - use previous image data when threshold exceeded
                    if idx > 0:
                        images[idx] = images[idx-1]
                        processed_angles.append(processed_angles[-1])
                    else:
                        images[idx] = diff
                        processed_angles.append(angles[idx])
            else:
                # Normal differential processing (no rotation applied regardless of flag)
                images[idx] = diff
                processed_angles.append(angles[idx])

            prev = curr
            self.progress_update.emit(20 + int(30 * idx / n_images))

        processed_images = images

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
            SID = float(dcm.RTImageSID)
            SAD = float(dcm.RadiationMachineSAD)
            SOD = SAD
            SDD = SID
            # Detector pixel pitch must come from the panel, not a hardcoded guess -
            # this dataset reports 0.336 mm at SID 1600, so a fixed 0.172 mm shrinks
            # the reconstructed FOV by ~2x and invalidates any distance-to-agreement.
            width = self.params.get('detector_pixel_mm') or float(dcm.ImagePlanePixelSpacing[0])
            delta_dd = width * SOD / SDD   # detector pitch projected back to isocentre

            # Create reconstructor
            reconstructor = OptimizedFDKReconstructor(SOD, delta_dd, self.params['image_size'])

            self.status_update.emit("Running FDK reconstruction...")
            self.progress_update.emit(50)

            # Reconstruct
            rec_image = self.reconstruct_chunked(
                reconstructor, processed_images, processed_angles,
                chunk_size=self.params.get('chunk_size', 25))
            rec_image = np.transpose(rec_image, (2,1,0))   # (x,y,z) -> (z,y,x)

            self.progress_update.emit(100)
            self.status_update.emit("Reconstruction complete!")

            # Physical extent of the reconstruction grid, mirroring the roi/hx/hy
            # arithmetic inside Fun_BackProjection. Reported in volume-axis order
            # (z, y, x) so the DICOM exporter can write real geometry instead of
            # inheriting the template's.
            n_rows, n_cols = processed_images.shape[1], processed_images.shape[2]
            MX = self.params['image_size']
            MZ = int(MX * n_rows / n_cols)
            dx = delta_dd * (n_cols - 1) / (MX - 1)
            dz = delta_dd * (n_rows - 1) / (MZ - 1)
            x0 = delta_dd * (-n_cols / 2.0 + 0.5)
            z0 = delta_dd * (-n_rows / 2.0 + 0.5)

            result = {
                'reconstructed_volume': rec_image,
                'parameters': {
                    'SOD': SOD,
                    'SDD': SDD,
                    'delta_dd': delta_dd,
                    'detector_pixel_mm': width,
                    'ramp_filter': USE_RAMP_FILTER,
                    'n_projections': len(processed_angles),
                    'image_size': self.params['image_size'],
                    'voxel_size_mm': (dz, dx, dx),
                    'origin_mm': (z0, x0, x0),
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
        fh_RL = None
        if USE_RAMP_FILTER:
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

                if USE_RAMP_FILTER:
                    projection_to_backproject = reconstructor.optimize_convolution(
                        weighted_projection, fh_RL)
                else:
                    projection_to_backproject = weighted_projection

                temp_rec = reconstructor.Fun_BackProjection(
                    projection_to_backproject, reconstructor.SOD, n_angles,
                    beta_rad[angle_idx], reconstructor.delta_dd, reconstructor.Nimage)

                chunk_result += temp_rec

            rec_image += chunk_result

            # Update progress
            progress = 50 + int(40 * (chunk_idx + 1) / len(chunks))
            self.progress_update.emit(progress)

        return rec_image
