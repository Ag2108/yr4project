import numpy as np
import matplotlib.pyplot as plt

# 1. Configuration
data_file = 'tbtest.dat'
output_image = 'band_structure.png'

# 2. Load the data
# np.loadtxt automatically handles varying whitespace between columns
try:
    data = np.loadtxt(data_file)
except FileNotFoundError:
    print(f"Error: Could not find '{data_file}'. Please ensure it is in the same directory.")
    exit()

# 3. Extract columns
# Column index 0 is the flux (\phi)
phi = data[:, 0]

# Columns from index 1 to the end are the eigenvalues (Energy bands)
energies = data[:, 1:]

# 4. Set up the plot
plt.figure(figsize=(8, 5))

# Plot all energy bands at once
# plt.plot gracefully handles 2D arrays by plotting each column as a line
plt.plot(phi, energies, color='forestgreen', linewidth=2)

# 5. Styling
plt.title('Energy Bands vs. Magnetic Flux, $N=7$', fontsize=14)
plt.xlabel(r'$\phi$ ($\phi_0$)', fontsize=12)
plt.ylabel('Energy (eV)', fontsize=12)

# Create a dashed grid similar to your reference images
plt.grid(True, linestyle='--', color='gray', alpha=0.7)

# Adjust axes to fit tightly around the data
plt.xlim([phi.min(), phi.max()])

# 6. Display and save
plt.tight_layout()
plt.savefig(output_image, dpi=300)
print(f"Plot successfully saved to {output_image}")

# Uncomment the line below if you want a pop-up window to appear when running the script
# plt.show()
