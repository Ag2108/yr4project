
import numpy as np
from math import gcd
from scipy.stats import gaussian_kde
import os
from collections import defaultdict
import subprocess

# System parameters
max_q=100

if not isinstance(max_q, int) or max_q <= 1:
        raise ValueError("max_q must be a positive integer, greater than 1")


for q in range(2, max_q+1):

    for p in range(1, q+1):

        if gcd(p, q) == 1:
            a = p / q

            N = q
            L = q

            print(f"Generating read.in for {N} sites)...")

            with open('read.in', 'w') as f:
                f.write("PHIMAX  =2\n")
                f.write("SIZE    =1\n")
                f.write("RANGE   =2\n")
                f.write(f"CELLSIZE={N}\n\n")
                
                # 1. Write the Sites
                f.write("# --- SITES (site_num, epsilon, t_vals(intra), t_vals(inter)) ---\n")
                for i in range(1, N + 1):
                    f.write(f"SITE    ={i} 2.0 1.0 0.1 1.0 0.1\n")
                    
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

                for i in range(1, N + 1):
                    # Calculate grid coordinates (0-indexed) to figure out neighbors
                    row = (i - 1) // L
                    col = (i - 1) % L

                    if col == L - 1:
                        f.write(f"CONNECTN ={i} {i - L + 1}\n")

                f.write(f"AAH = {a} 1")
                        
                # CRITICAL: Always end with EXIT
                f.write("\nEXIT\n")

            print("Running Fortran", q)
            try:
                subprocess.run(["make", "run"], check=True)
            except subprocess.CalledProcessError:
                print("Error: Fortran program crashed!")
                exit()

            if os.path.exists("tbtest.dat"):
                with open("tbtest.dat", "r") as f:
                    # Read the file line-by-line instead of all at once
                    lines = f.readlines()
                
                for line in lines:
                    line = line.strip()
                    if not line:
                        continue  # Skip any blank lines
                        
                    parts = line.split()
                    
                    # Extract the AB phase (first column) for this specific row
                    raw_ab_phase = parts[0]
                    ab_phase = f"{float(raw_ab_phase):.3f}"
                    
                    # Format the phi_key (assuming 'a' is defined in your outer loop)
                    phi_key = f"{a:.6f}"
                    
                    # Gather the energy values for this row
                    energies = "\t".join(parts[1:])
                    
                    # Route to the correct file based on the AB phase
                    target_filename = f"master_{ab_phase}.dat"
                    
                    # Append to that specific file
                    with open(target_filename, "a") as target_f:
                        target_f.write(f"{phi_key}\t{energies}\n")