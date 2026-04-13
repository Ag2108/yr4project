import numpy as np
import matplotlib.pyplot as plt
import os

# 1. Parse the read.in file to determine the system size
sys_size = 1
cell_size = 1
rng = 1

try:
    with open('read.in', 'r') as f:
        for line in f:
            # Skip comments and empty lines
            if line.strip().startswith('#') or '=' not in line:
                continue
            
            # Split the line by the equals sign
            key, val = line.split('=', 1)
            key = key.strip()
            # Extract the first value after the equals sign (ignoring trailing comments)
            val = val.split()[0].strip()
            
            if key == 'SIZE':
                sys_size = int(val)
            elif key == 'CELLSIZE':
                cell_size = int(val)
            elif key == 'RANGE':
                rng = int(val)
except FileNotFoundError:
    print("Warning: 'read.in' not found in the directory. Defaulting N to 8.")
    sys_size = 1
    cell_size = 8

# Calculate total number of sites
N = sys_size * cell_size

# 2. Configuration
data_file = 'tbtest.dat'
output_image = f'band_structure_N{N}_R{rng}SSH.png' # Dynamically named output file

# 3. Load the data
try:
    data = np.loadtxt(data_file)
except FileNotFoundError:
    print(f"Error: Could not find '{data_file}'. Please ensure it is in the same directory.")
    exit()

# 4. Extract columns
phi = data[:, 0]
energies = data[:, 1:]

# 5. Set up the plot
plt.figure(figsize=(8, 5))

# Plot all energy bands at once
plt.plot(phi, energies, color='forestgreen', linewidth=2)

# 6. Styling
# Dynamically update the title with the calculated N
plt.title(f'Energy Bands vs. Magnetic Flux, $N={N}$, $NN={rng}$', fontsize=14)
plt.xlabel(r'$\phi$ ($\phi_0$)', fontsize=12)
plt.ylabel('Energy (eV)', fontsize=12)

# Create a dashed grid
plt.grid(True, linestyle='--', color='gray', alpha=0.7)

# Adjust axes to fit tightly around the data
plt.xlim([phi.min(), phi.max()])

# 7. Display and save
plt.tight_layout()
plt.savefig(output_image, dpi=300)
print(f"Plot successfully saved to {output_image}")

# plt.show()