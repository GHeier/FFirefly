import numpy as np
import tbmodels

def Hr_to_Ek(filename, num_pts, mu, window=1.0):

# Load tight-binding model
    model = tbmodels.Model.from_wannier_files(hr_file=filename)

# Create reduced k-point grid
    grid = np.linspace(0, 1, n_k, endpoint=False)
    kx, ky, kz = np.meshgrid(grid, grid, grid, indexing='ij')
    kx_flat = kx.flatten()
    ky_flat = ky.flatten()
    kz_flat = kz.flatten()

# Determine number of bands
    sample_eigvals = model.eigenval([0.0, 0.0, 0.0])
    num_bands = len(sample_eigvals)

# Pass 1: Detect active bands
    active_bands = np.zeros(num_bands, dtype=bool)

    num_kpts = len(kx_flat)
    eigvals_all = np.empty((num_kpts, num_bands), dtype=np.float64)

# Fast loop to fill eigenvalues
    for i in range(num_kpts):
        kpt = [kx_flat[i], ky_flat[i], kz_flat[i]]
        eigvals_all[i] = model.eigenval(kpt)
        perc = i / num_kpts
        print(f"  {perc * 50:.1f}% ", end='\r') 

# Find active bands
    within_window = np.abs(eigvals_all - mu) <= window
    active_bands = np.any(within_window, axis=0)

# Map original band indices to filtered output indices starting from 1
    active_band_indices = [i for i, active in enumerate(active_bands) if active]
    band_index_map = {orig_n: new_n+1 for new_n, orig_n in enumerate(active_band_indices)}

# Map original band index -> output band index (1-based)
    active_band_indices = np.where(active_bands)[0]
    band_index_map = {orig: i + 1 for i, orig in enumerate(active_band_indices)}

# Prepare output
    output_lines = ["# kx\tky\tkz\tn\tf"]

# Format helper
    def fmt(x): return f"{x: .6f}"

# Vectorized output generation
    for i in range(num_kpts):
        kx_val = kx_flat[i]
        ky_val = ky_flat[i]
        kz_val = kz_flat[i]
        for orig_n in active_band_indices:
            new_n = band_index_map[orig_n]
            e = eigvals_all[i, orig_n]
            line = f"{fmt(kx_val)}\t{fmt(ky_val)}\t{fmt(kz_val)}\t{new_n}\t{fmt(e)}"
            output_lines.append(line)
        perc = i / num_kpts
        print(f"  {50 + perc * 50:.1f}% ", end='\r') 

# Write to file
    with open("band_kgrid_filtered.dat", "w") as f:
        f.write("\n".join(output_lines))

    print()
    print(f"Saved {len(active_band_indices)} active band(s) to 'band_kgrid_filtered.dat'")

#num_pts = 25  
#mu = 11.2
#window = 1.0  # +/- 1 eV around mu
#Hr_to_Ek('Pb_hr.dat', num_pts, mu, window=1.0)
