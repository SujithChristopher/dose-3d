"""
EPID Dose Reconstruction GUI
Advanced PySide6 GUI for cone beam FDK reconstruction with dose verification.

Features:
- Interactive 3D volume viewer with slice navigation
- Multiple color mapping options
- Real-time dose verification against TPS reference
- Mouse cursor position tracking
- Scroll wheel slice navigation
- Optimized reconstruction pipeline
- Batch dose calibration (pixel value <-> dose linear fit)
"""
import sys

from PySide6.QtWidgets import QApplication, QMessageBox

from src.optional_deps import PYDICOM_AVAILABLE
from src.main_window import EPIDReconstructionGUI
from src._version import __version__
from src.updater import check_for_update_async


def main():
    """Main application entry point"""
    app = QApplication(sys.argv)

    # Set application properties
    app.setApplicationName("EPID Reconstruction GUI")
    app.setApplicationVersion(__version__)

    # Check dependencies
    missing_deps = []
    if not PYDICOM_AVAILABLE:
        missing_deps.append("pydicom")

    if missing_deps:
        QMessageBox.critical(None, "Missing Dependencies",
                           f"Missing required packages: {', '.join(missing_deps)}\n\n"
                           "Please install with: pip install " + " ".join(missing_deps))
        return 1

    # Create and show main window
    window = EPIDReconstructionGUI()
    window.show()

    check_for_update_async(window)

    return app.exec()

if __name__ == "__main__":
    sys.exit(main())
