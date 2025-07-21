#!/usr/bin/env python3
"""
Simple launcher for EPID Reconstruction GUI
This script handles environment setup and launches the main GUI application.
"""

import os
import sys
import subprocess

def main():
    """Launch the EPID Reconstruction GUI"""
    
    print("EPID Dose Reconstruction GUI Launcher")
    print("=" * 40)
    
    # Check if we're in the correct directory
    if not os.path.exists("epid_reconstruction_gui.py"):
        print("Error: epid_reconstruction_gui.py not found in current directory")
        print("Please run this script from the dose-3d directory")
        return 1
    
    # Check if TPS reference data exists
    if not os.path.exists("tps_reference_data.py"):
        print("Warning: TPS reference data not found")
        print("Dose verification will be disabled")
    
    # Check test components first
    print("Running component tests...")
    try:
        result = subprocess.run([sys.executable, "test_gui_components.py"], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            print("Component tests failed:")
            print(result.stderr)
            return 1
        print("Component tests passed!")
    except Exception as e:
        print(f"Could not run component tests: {e}")
    
    # Launch GUI
    print("\nLaunching EPID Reconstruction GUI...")
    print("Close the GUI window to exit.")
    
    try:
        subprocess.run([sys.executable, "epid_reconstruction_gui.py"])
    except KeyboardInterrupt:
        print("\nGUI launch cancelled")
    except Exception as e:
        print(f"Error launching GUI: {e}")
        return 1
    
    print("GUI closed successfully")
    return 0

if __name__ == "__main__":
    sys.exit(main())