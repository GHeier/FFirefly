### ** 🔹 Current Calculation Categories**
- cateogry: `test`(default) - Runs tests to confirm all categories are functioning correctly
- cateogry: `superconductor` - Calculating gap size and symmetry
   - method: `bcs`, `eliashberg`
   - calculation: `diagonalization` of exact gap functions or `projection` of gap functions onto a basis set
   - FS_only: `true` to calculate the gap only on the Fermi surface or `false` for the entire Brillouin zone
- cateogry: `response` - Calculating bare susceptibility
   - method: `sparse_ir`(sparse matsubaras from convolution theorem) operates at finite T and requires a dense k-grid. `libtetrabz`(analytic tetrahedra integration) operates at 0 T and allows for a sparser grid. If it throws a "STOP NESTING" then alter the grid size
   - dynamic: `true` to save data at finite $\omega$ or `false` for $\omega=0$
   - wpts: number of Matsubara frequencies to use. wpts=1 calculates iv=0
- cateogry: `vertex` - Calculating a given 2-particle vertex
   - interaction: `FLEX` calculates the FLEX vertex. Takes in bare susceptibility (_chi.dat) as input.
   - scf: `true` to find self-consistent vertex (and output self-energy) or `false` for 0th-order
- cateogry: `DOS` - Calculating Density of States
   - method: `libtetrabz`(analytic tetrahedra integration) or `surface_sum` both give similar results.
- cateogry: `fermi_surface` - Calculates Fermi Surface


