"""Load and apply pixel-value <-> dose linear calibrations produced by the Dose Calibration tab."""
import json


def load_calibration(path):
    """Load a calibration JSON (as exported by DoseCalibrationWidget) and return its dict."""
    with open(path, 'r') as f:
        return json.load(f)


def apply_calibration(volume, calibration, include_intercept=False):
    """Scale a reconstructed pixel-value volume into dose (cGy) using a loaded calibration.

    dose = slope * pixel_value [+ intercept]

    intercept is a per-ROI background offset from the calibration fit, not a
    per-voxel physical quantity, so applying it everywhere is opt-in and off
    by default (slope-only scaling).
    """
    slope = calibration['fit']['slope']
    intercept = calibration['fit']['intercept'] if include_intercept else 0.0
    return slope * volume + intercept
