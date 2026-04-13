import numpy as np
import matplotlib.pyplot as plt
import os

# 1. Parse the read.in file to determine the system parameters
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
            # Extract the first value after the equals sign
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
energy_file = 'tbtest.dat'
# Automatically detect if you used currents.dat or currenttest.dat in Main.f90
current_file = 'currents.dat' if os.path.exists('currents.dat') else 'currenttest.dat'
output_image = f'band_structure_current_N{N}_R{rng}.png' # Dynamically named output file

# 3. Load the data
try:
    energy_data = np.loadtxt(energy_file)
except FileNotFoundError:
    print(f"Error: Could not find '{energy_file}'. Please ensure it is in the same directory.")
    exit()

try:
    current_data = np.loadtxt(current_file)
except FileNotFoundError:
    print(f"Error: Could not find '{current_file}'. Did the Fortran script generate it?")
    exit()

# 4. Extract columns
phi_e = energy_data[:, 0]
energies = energy_data[:, 1:]

phi_c = current_data[:, 0]
currents = current_data[:, 1:]

# 5. Set up the plot with 2 stacked subplots sharing the X-axis
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

# --- TOP PLOT: Energy Bands ---
ax1.plot(phi_e, energies, color='forestgreen', linewidth=2)
ax1.set_title(f'Energy Bands and Persistent Current, $N={N}$, $NN={rng}$', fontsize=14)
ax1.set_ylabel('Energy $-E_F$ (eV)', fontsize=12)
ax1.grid(True, linestyle='--', color='gray', alpha=0.7)
ax1.set_xlim([phi_e.min(), phi_e.max()])

# --- BOTTOM PLOT: Persistent Currents ---
# Using a contrasting color (firebrick red) for the currents
ax2.plot(phi_c, currents, color='firebrick', linewidth=1.5)
ax2.set_xlabel(r'$\phi$ ($\phi_0$)', fontsize=12)
ax2.set_ylabel(r'Current ($I/I_0$)', fontsize=12)
ax2.grid(True, linestyle='--', color='gray', alpha=0.7)
ax2.set_xlim([phi_c.min(), phi_c.max()])

# 6. Display and save
plt.tight_layout()
plt.savefig(output_image, dpi=300)
print(f"Plot successfully saved to {output_image}")

# Uncomment below to show the pop-up window
# plt.show()