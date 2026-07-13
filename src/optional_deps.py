"""Centralized third-party / optional dependency imports and availability flags.

Every module in this package imports from here instead of repeating
try/except boilerplate, so all of them see the same import-success state.
Fallback names are always bound (to None where there's no safe default) so
that importing a module never fails just because an optional dependency is
missing - failure is deferred to the point where the feature is actually used,
matching the original single-file script's behavior.
"""
import math
import warnings

warnings.filterwarnings('ignore')

PI = math.pi

# Matplotlib Qt backend selection
import matplotlib
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

# DICOM loading / image processing
try:
    from pydicom import dcmread
    from concurrent.futures import ThreadPoolExecutor
    from tqdm import tqdm
    import cv2
    PYDICOM_AVAILABLE = True
except ImportError as e:
    print(f"Import warning: {e}")
    PYDICOM_AVAILABLE = False
    dcmread = None
    ThreadPoolExecutor = None
    tqdm = None
    cv2 = None

# TPS reference utilities (repo-root module, not part of this package)
try:
    from tps_reference_data import export_reconstructed_dicom, calculate_volume_statistics, load_tps_dicom_template
    TPS_AVAILABLE = True
except ImportError:
    print("TPS utilities not available")
    TPS_AVAILABLE = False
    export_reconstructed_dicom = None
    calculate_volume_statistics = None
    load_tps_dicom_template = None

# Numba JIT (optional performance boost, unused directly by the GUI today)
try:
    from numba import jit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def jit(func):
        return func
    prange = range
