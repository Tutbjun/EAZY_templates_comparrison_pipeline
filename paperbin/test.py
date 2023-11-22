import numpy as np
from scipy.integrate import simps

# Sample spectrum and filter response data (replace with your actual data)
wavelengths_spectrum = np.linspace(300, 900, 100)
flux_spectrum = np.random.rand(100)
wavelengths_filter = np.linspace(400, 700, 100)
filter_response = np.random.rand(100)

# Interpolate the filter response to match the wavelength grid of the spectrum
filter_response_interp = np.interp(wavelengths_spectrum, wavelengths_filter, filter_response)

# Calculate the product of the spectrum and the interpolated filter response
product = flux_spectrum * filter_response_interp

# Integrate the product using the Simpson's rule
bandpass_value = simps(product, wavelengths_spectrum)

print(f"Bandpass value: {bandpass_value:.3e}")