import time

# System parameters
L = 3 # Width and height of the grid (100x100 = 10,000 sites)
      # L^2 x L^2 Hamiltonian; BEWARE
N = L * L

print(f"Generating read.in for a {L}x{L} grid ({N} sites)...")
start_time = time.time()

with open('read.in', 'w') as f:
    f.write(f"#READ.IN file for {L}x{L} 2D Square Lattice (Landau Gauge)\n")
    f.write("SIZE    =1\n")
    f.write("RANGE   =1\n")
    f.write(f"CELLSIZE={N}\n\n")
    
    # 1. Write the Sites
    f.write("# --- SITES (site_num, epsilon, t_vals(intra), t_vals(inter)) ---\n")
    for i in range(1, N + 1):
        f.write(f"SITE    ={i} 0.0 1.0 1.0\n")
        
    f.write("\n")
    
    # 2. Write the Connections
    f.write("# --- CONNECTIONS (horizontal and vertical) ---\n")
    for i in range(1, N + 1):
        # Calculate grid coordinates (0-indexed) to figure out neighbors
        row = (i - 1) // L
        col = (i - 1) % L
        
        # Connect to the Right (Horizontal)
        if col < L - 1:
            f.write(f"CONNECT ={i} {i + 1}\n")
            
        # Connect to the Top (Vertical)
        if row < L - 1:
            f.write(f"CONNECT ={i} {i + L}\n")

    for i in range(1, N + 1):
        # Calculate grid coordinates (0-indexed) to figure out neighbors
        row = (i - 1) // L
        col = (i - 1) % L

        if col == L - 1:
            f.write(f"CONNECTN ={i} {i - L + 1}\n")
            
    # CRITICAL: Always end with EXIT
    f.write("\nEXIT\n")

print(f"Done in {time.time() - start_time:.3f} seconds! 'read.in' is ready.")