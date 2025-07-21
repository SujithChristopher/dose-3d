"""
Test script to verify GUI components work correctly
"""

import sys
import numpy as np

# Test imports
print("Testing imports...")

try:
    from PySide6.QtWidgets import QApplication
    print("PASS - PySide6 import successful")
except ImportError as e:
    print(f"FAIL - PySide6 import failed: {e}")

try:
    import matplotlib
    matplotlib.use('QtAgg')
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    print("PASS - Matplotlib QtAgg backend successful")
except ImportError:
    try:
        matplotlib.use('Qt5Agg')
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        print("PASS - Matplotlib Qt5Agg backend successful")
    except ImportError as e:
        print(f"FAIL - Matplotlib Qt backend failed: {e}")

try:
    from pydicom import dcmread
    print("PASS - PyDICOM import successful")
except ImportError as e:
    print(f"FAIL - PyDICOM import failed: {e}")

try:
    from tps_reference_data import export_reconstructed_dicom, calculate_volume_statistics, load_tps_dicom_template
    print("PASS - TPS utilities import successful")
except ImportError as e:
    print(f"FAIL - TPS utilities import failed: {e}")

try:
    from scipy.fft import fft, ifft
    print("PASS - SciPy FFT import successful")
except ImportError as e:
    print(f"FAIL - SciPy FFT import failed: {e}")

try:
    from concurrent.futures import ThreadPoolExecutor
    print("PASS - Threading support available")
except ImportError as e:
    print(f"FAIL - Threading import failed: {e}")

try:
    from numba import jit
    print("PASS - Numba JIT compilation available")
except ImportError:
    print("WARN - Numba not available (optional)")

# Test basic functionality
print("\nTesting basic functionality...")

# Test volume creation and viewing
test_volume = np.random.rand(50, 50, 50) * 1000
print(f"PASS - Test volume created: {test_volume.shape}")

# Test TPS utilities
try:
    from tps_reference_data import calculate_volume_statistics, load_tps_dicom_template
    stats = calculate_volume_statistics(test_volume)
    template = load_tps_dicom_template()
    if stats and 'min_value' in stats:
        print("PASS - Volume statistics working")
        print(f"  Stats: min={stats['min_value']:.2f}, max={stats['max_value']:.2f}")
    if template is not None:
        print("PASS - TPS template loaded")
    else:
        print("WARN - TPS template not available")
except Exception as e:
    print(f"FAIL - TPS utilities failed: {e}")

print("\nAll component tests completed!")
print("\nTo run the GUI:")
print("python epid_reconstruction_gui.py")