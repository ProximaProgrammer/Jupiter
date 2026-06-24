import numpy as np

# Exact parameters matching the NASA JIRAM LBL configuration
DAT_FILE = 'JIR_SPE_RESPONSIVITY_V02.DAT'
BANDS = 336
SAMPLES = 256

print("--- RE-EVALUATING JIRAM SPE BINARY MATRIX ---")

# 1. Read the raw array directly using Big-Endian float ('>f4')
raw_data = np.fromfile(DAT_FILE, dtype='>f4')
print(f"Total values extracted from file: {len(raw_data)}")

# 2. FUNDAMENTAL FIX: Reshape using (SAMPLES, BANDS) based on AXIS_ITEMS ordering
try:
    # 256 spatial rows, 336 spectral columns
    responsivity_matrix = raw_data.reshape(SAMPLES, BANDS)
    print(f"Successfully loaded true Responsivity Matrix: {responsivity_matrix.shape}")
    
    # 3. Save the clean, non-corrupted responsivity table grid
    headers = "Sample_Row\t" + "\t".join([f"Bin_{i+1}" for i in range(BANDS)])
    row_indices = np.arange(1, SAMPLES + 1).reshape(-1, 1)
    final_grid = np.hstack((row_indices, responsivity_matrix))
    
    # We save using '%e' because calibration units exist at scientific scales
    np.savetxt("jiram4.txt", final_grid, fmt=["%d"] + ["%e"]*BANDS, delimiter="\t", header=headers)
    print(" -> File written: 'jiram4.txt'")
    
    # Check a localized sample - real values should track closely with clean scientific boundaries
    print(f"\nVerification check (Row 1, Bin 1-3): {responsivity_matrix[0, 0:3]}")

except Exception as e:
    print(f"Error structuring matrix: {e}")
