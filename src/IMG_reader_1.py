# print("First a preview test")

# try:
#     # JIRAM data is usually 32-bit Big-Endian floats ('>f4')
#     # If the image looks like 'noise', try '>u2' (16-bit Unsigned Integer)
#     data = np.fromfile('data/JIR_IMG_RDR_2023250T010133_V01.IMG', dtype='<f4')
#     image = data.reshape((256, 432))

#     # Calculate the 2nd and 98th percentile to ignore extreme white dots
#     vmin, vmax = np.percentile(image, [2, 98])

#     plt.figure(figsize=(10, 6))
# # Apply the vmin/vmax to fix the contrast
#     plt.imshow(image, cmap='magma', vmin=vmin, vmax=vmax)
#     plt.colorbar(label='Intensity (Clipped)')
#     plt.title("Juno JIRAM: Contrast Enhanced")
#     from scipy.ndimage import median_filter
#     denoised_image = median_filter(image, size=3)
#     plt.imshow(denoised_image, cmap='magma', vmin=np.percentile(denoised_image, 2), vmax=np.percentile(denoised_image, 98))
#     plt.show()
    
# except Exception as e:
#     print(f"Error: {e}")

print("Running...")

'''
:::Agenda:::
add all IMG files (first section) of https://atmos.nmsu.edu/PDS/data/PDS4/juno_jiram_bundle/data_calibrated/orbit10/
'''

# Import libraries
import numpy as np
import glob, os
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, find_peaks, peak_widths
import pywt
from sklearn.decomposition import PCA

from astropy.io import fits
from planetaryimage import PDS3Image

# Load all .IMG files safely
data_folder = "data"
lines = 256
samples = 432
dtype = ">f4"

img_files = sorted(glob.glob(os.path.join(data_folder, "*.IMG")))
spectra = []

def resample(arr, n=1000):
    x_old = np.linspace(0, 1, len(arr))
    x_new = np.linspace(0, 1, n)
    return np.interp(x_new, x_old, arr)

for count,f in enumerate(img_files):
    print(count+1)
    try:
        data = np.fromfile(f, dtype=dtype)

        if data.size != 256*432:
            print(f"Skipping {f}: wrong size")
            continue

        data = data.reshape((256, 432))
        intensity = data.flatten()

        # Remove NaN/Inf
        intensity = intensity[np.isfinite(intensity)]

        # Remove extreme outliers (top 1%)
        cutoff = np.percentile(intensity, 99)
        intensity = intensity[intensity < cutoff]

        if len(intensity) < 100:
            print(f"Skipping {f}: too little valid data")
            continue

        # Normalize
        max_val = np.max(intensity)
        if max_val != 0:
            intensity = intensity / max_val

        # Resample
        intensity = resample(intensity)

        spectra.append(intensity)

    except Exception as e:
        print("Skipped:", f, e)

spectra = np.array(spectra)
print("Loaded spectra shape:", spectra.shape)

# ANALYSIS & PLOTS (Step 3 goes here) --> is there stuff missing??

# Use first spectrum for example --> I thought we were using all the data...
intensity = spectra[0]
wavelength = np.arange(intensity.size)

# Smooth and get residual
smooth = savgol_filter(intensity, 11, 3)
residual = intensity - smooth

peaks, _ = find_peaks(residual, prominence=0.003)
dips, _ = find_peaks(-residual, prominence=0.003)
widths = peak_widths(residual, peaks, rel_height=0.5)

# Wavelet transform
scales = np.arange(1, 80)
coeffs, freqs = pywt.cwt(residual, scales, 'morl')

pca = PCA(n_components=2)
components = pca.fit_transform(spectra)


# Plots
fig, axs = plt.subplots(4, 1, figsize=(10, 8))

axs[0].plot(wavelength, intensity)
axs[0].set_title("Normalized Spectrum")

axs[1].plot(wavelength, residual, label="Residual")
axs[1].scatter(wavelength[peaks], residual[peaks], color='r', label="Peaks")
axs[1].scatter(wavelength[dips], residual[dips], color='g', label="Dips")
axs[1].legend()
axs[1].set_title("Spectral Features")

axs[2].imshow(coeffs, aspect='auto', extent=[wavelength.min(), wavelength.max(), scales.max(), scales.min()])
axs[2].set_title("Wavelet Transform")

axs[3].scatter(components[:,0], components[:,1])
axs[3].set_title("PCA of Spectra")
axs[3].set_xlabel("PC1")
axs[3].set_ylabel("PC2")

plt.tight_layout()
plt.show()


'''
# =========================
# SETTINGS
# =========================
data_folder = "data"  # folder containing .IMG/.LBL files
img_files = sorted(glob.glob(os.path.join(data_folder, "*.IMG")))

if len(img_files) == 0:
    raise Exception("No .IMG files found in folder!")

# =========================
# LOAD MULTIPLE SPECTRA
# =========================
spectra = []
wavelength = None

for f in img_files:
    try:
        # PDS3 .IMG often works with astropy.io.fits for direct array reading
        # If fails, np.fromfile may be needed
        data = fits.getdata(f)
        # Assuming data has 2 columns: wavelength, intensity
        if data.ndim == 2 and data.shape[1] >= 2:
            wl = data[:,0]
            intensity = data[:,1]
        else:
            # fallback for single column -> just intensity
            intensity = data.flatten()
            wl = np.arange(len(intensity))

        # Remove NaNs
        mask = np.isfinite(wl) & np.isfinite(intensity)
        wl = wl[mask]
        intensity = intensity[mask]

        # Normalize
        intensity = intensity / np.max(intensity)

        spectra.append(intensity)
        if wavelength is None:
            wavelength = wl

    except Exception as e:
        print("Skipped file:", f, "| Error:", e)

spectra = np.array(spectra)
print("Loaded spectra:", spectra.shape)

# =========================
# USE FIRST SPECTRUM FOR FEATURE ANALYSIS
# =========================
intensity = spectra[0]

# Smooth + residual
smooth = savgol_filter(intensity, 11, 3)
residual = intensity - smooth

# Peaks and dips
peaks, _ = find_peaks(residual, prominence=0.003)
dips, _ = find_peaks(-residual, prominence=0.003)

# Peak widths
widths = peak_widths(residual, peaks, rel_height=0.5)

# Wavelet
scales = np.arange(1, 80)
coeffs, freqs = pywt.cwt(residual, scales, 'morl')

# PCA across spectra
pca = PCA(n_components=2)
components = pca.fit_transform(spectra)

# Autocorrelation
corr = np.correlate(residual, residual, mode='full')

# =========================
# PLOTS
# =========================

# Raw normalized spectrum
plt.figure()
plt.plot(wavelength, intensity)
plt.title("Normalized Spectrum")
plt.xlabel("Wavelength")
plt.ylabel("Normalized Intensity")

# Residual + peaks/dips
plt.figure()
plt.plot(wavelength, residual, label="Residual")
plt.scatter(wavelength[peaks], residual[peaks], label="Peaks", color='r')
plt.scatter(wavelength[dips], residual[dips], label="Dips", color='g')
plt.title("Spectral Features")
plt.xlabel("Wavelength")
plt.ylabel("Residual Intensity")
plt.legend()

# Wavelet transform
plt.figure()
plt.imshow(coeffs, aspect='auto', extent=[wavelength.min(), wavelength.max(), scales.max(), scales.min()])
plt.title("Wavelet Transform")
plt.xlabel("Wavelength")
plt.ylabel("Scale")

# PCA
plt.figure()
plt.scatter(components[:,0], components[:,1])
plt.title("PCA of Spectra")
plt.xlabel("PC1")
plt.ylabel("PC2")

# Autocorrelation
plt.figure()
plt.plot(corr)
plt.title("Autocorrelation of Residual")

plt.show()

# =========================
# OUTPUT
# =========================
print("Number of spectra loaded:", len(spectra))
print("Number of peaks in first spectrum:", len(peaks))
print("Average peak width:", np.mean(widths[0]))
'''