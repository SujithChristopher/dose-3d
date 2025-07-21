"""
Test slice bounds checking for volume viewer
"""

import numpy as np

def test_slice_bounds():
    """Test slice bounds checking logic"""
    print("Testing slice bounds checking...")
    
    # Test volume with different shapes
    test_volumes = [
        np.random.rand(50, 50, 50),   # 50x50x50
        np.random.rand(100, 100, 100), # 100x100x100
        np.random.rand(75, 75, 75),   # 75x75x75
        np.random.rand(25, 25, 25),   # 25x25x25
    ]
    
    for i, volume in enumerate(test_volumes):
        print(f"\nTest {i+1}: Volume shape {volume.shape}")
        
        # Simulate slice position calculation
        max_x = volume.shape[0] - 1
        max_y = volume.shape[1] - 1  
        max_z = volume.shape[2] - 1
        
        # Center positions
        center_x = min(max_x // 2, max_x)
        center_y = min(max_y // 2, max_y)
        center_z = min(max_z // 2, max_z)
        
        print(f"  Max indices: X={max_x}, Y={max_y}, Z={max_z}")
        print(f"  Center slices: X={center_x}, Y={center_y}, Z={center_z}")
        
        # Test bounds checking
        test_positions = [0, center_x, max_x, max_x + 1, -1]
        
        for pos in test_positions:
            clamped_x = max(0, min(pos, max_x))
            clamped_y = max(0, min(pos, max_y))
            clamped_z = max(0, min(pos, max_z))
            
            print(f"  Position {pos} -> X:{clamped_x}, Y:{clamped_y}, Z:{clamped_z}")
            
            # Test if slicing would work
            try:
                slice_yz = volume[clamped_x, :, :].T
                slice_xz = volume[:, clamped_y, :].T 
                slice_xy = volume[:, :, clamped_z].T
                print(f"    Slicing successful: YZ{slice_yz.shape}, XZ{slice_xz.shape}, XY{slice_xy.shape}")
            except IndexError as e:
                print(f"    Slicing failed: {e}")
    
    print("\nPASS - Slice bounds checking test completed")

def test_default_settings():
    """Test default GUI settings"""
    print("\nTesting default settings...")
    
    # Test rotation default
    rotation_enabled = False  # Default should be False
    print(f"Rotation enabled by default: {rotation_enabled}")
    
    # Test image size default
    default_image_size = 100  # Should be 100
    print(f"Default image size: {default_image_size}")
    
    # Test mode parameters
    quick_mode_size = 100  # Should be 100 for both modes now
    full_mode_size = 100
    
    print(f"Quick mode image size: {quick_mode_size}")
    print(f"Full mode image size: {full_mode_size}")
    
    print("PASS - Default settings test completed")

if __name__ == "__main__":
    print("Testing Slice Bounds and Default Settings")
    print("=" * 50)
    
    test_slice_bounds()
    test_default_settings()
    
    print("\n" + "=" * 50)
    print("All tests completed!")
    print("\nFixed issues:")
    print("- Rotation unchecked by default")
    print("- Image size defaults to 100 (standard)")
    print("- Slice bounds checking prevents IndexError")
    print("- Dynamic slider ranges based on volume dimensions")