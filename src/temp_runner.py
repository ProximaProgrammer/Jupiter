import numpy as np
import pvl

# 1. Load the PDS structural instructions from the label file
lbl_path = 'JIR_IMG_M_RESPONSIVITY_V03.LBL'
dat_path = 'JIR_IMG_M_RESPONSIVITY_V03.DAT'

label = pvl.load(lbl_path)

# 2. Find the structural layout object inside the PDS tracking dict
# JIRAM calibration files pack matrix variables under a 'TABLE' or 'COMPRESS_TABLE' object block
table_meta = label.get('TABLE') or label.get('MATRIX')
if not table_meta:
    # Fallback to general file characteristics if root-level definitions are missing
    columns = 432  # Standard JIRAM detector width axis mapping
    rows = 128     # Standard JIRAM L-band imaging window block
else:
    columns = table_meta.get('COLUMNS', 432)
    rows = table_meta.get('ROWS', 128)

print(f"Extracting JIRAM matrix layout based on LBL specifications: {rows} rows x {columns} columns")

# 3. Read the Big Endian data stream
raw_data = np.fromfile(dat_path, dtype='>f4')

# 4. Extract just the active table metrics, ignoring trailing frame layout padding
total_elements = rows * columns
active_matrix = raw_data[:total_elements].reshape(rows, columns)

# 5. Clean mask check: Identify and print structural data bounds 
null_mask_value = -1.213305e-38
valid_points = ~np.isclose(active_matrix, null_mask_value, rtol=1e-5)

print(f"Total entries in active slice: {active_matrix.size}")
print(f"Valid operational coefficients found: {np.sum(valid_points)}")

# 6. Save as an explicitly aligned ASCII text file mapping the true calibration coefficients
np.savetxt("GOOGOOGAAGAA.txt", active_matrix, fmt="%e", delimiter="\t")
print("Successfully generated 'GOOGOOGAAGAA.txt'")
