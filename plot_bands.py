import numpy as np
import matplotlib.pyplot as plt
import os

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
    print("Warning: 'read.in' not found in the directory. Defaulting to N=8.")
    sys_size = 1
    cell_size = 8

N = sys_size * cell_size

data_file = 'tbtest.dat'
output_image = f'band_structure_N{N}_R{rng}SSH.png' # Dynamically named output file

try:
    data = np.loadtxt(data_file)
except FileNotFoundError:
    print(f"Error: Could not find '{data_file}'. Please ensure it is in the same directory.")
    exit()


phi = data[:, 0]
energies = data[:, 1:]

plt.figure(figsize=(8, 5))

plt.plot(phi, energies, color='forestgreen', linewidth=2)

plt.title(f'Energy Bands vs. Magnetic Flux, $N={N}$, $NN={rng}$', fontsize=14)
plt.xlabel(r'$\phi$ ($\phi_0$)', fontsize=12)
plt.ylabel('Energy (eV)', fontsize=12)

plt.grid(True, linestyle='--', color='gray', alpha=0.7)

plt.xlim([phi.min(), phi.max()])

plt.tight_layout()
plt.savefig(output_image, dpi=300)
print(f"Plot successfully saved to {output_image}")

# plt.show()