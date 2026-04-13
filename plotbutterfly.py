import matplotlib.pyplot as plt

# Lists to store our coordinates for the scatter plot
phi_vals = []
energy_vals = []

# Path to your master data file
file_path = 'master_1.571.dat'

try:
    with open(file_path, 'r') as f:
        for line in f:
            # Clean up the line and skip empty lines or comments
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) > 1:
                try:
                    # The first column is our Magnetic Flux Ratio (phi)
                    phi = float(parts[0])
                    
                    # The rest of the columns are the Energy eigenvalues
                    # We iterate through them and pair each one with the phi value
                    energies = [float(e) for e in parts[1:]]
                    for e in energies:
                        phi_vals.append(phi)
                        energy_vals.append(e)
                except ValueError:
                    # Skips lines that don't contain valid numbers
                    continue

    print(f"Successfully loaded {len(phi_vals)} data points.")
    
    # Create the figure
    plt.figure(figsize=(10, 8))
    
    # Adaptive point size: if we have hundreds of thousands of points, 
    # we make them tiny so the fractal gaps remain sharp.
    point_size = 0.2 if len(phi_vals) > 100000 else 1.0
    
    # Plotting the data
    plt.scatter(phi_vals, energy_vals, s=point_size, color='teal', alpha=0.6)
    
    # Formatting the plot with LaTeX notation
    plt.title("Hofstadter Butterfly Cross-Section (AB Phase = 1.571)", fontsize=16)
    plt.xlabel(r"AAH Ratio ($\alpha = p/q$)", fontsize=14)
    plt.ylabel(r"Energy ($E$)", fontsize=14)
    
    # Visual aids
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.xlim(0, 1)
    
    plt.tight_layout()
    
    # Save the output
    plt.savefig('hofstadter_plot4.png', dpi=300)
    print("Plot saved as 'hofstadter_plot.png'")

except Exception as e:
    print(f"An error occurred: {e}")