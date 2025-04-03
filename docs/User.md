### ** 🔹 Current Calculation Categories**
- Superconducting Gap calculations, both BCS and Eliashberg
   - cateogry name: `superconductor`
   - methods: `bcs`, `eliashberg`
   - calculation: `diagonalization` of exact gap functions or `projection` of gap functions onto a basis set
   - FS_only: `true` or `false` to calculate the gap only on the Fermi surface or the entire Brillouin zone
- Compute **self-energy (Σ)** and **2PI vertices (Γ) (with possible corrections)**  
- Response functions (polarization, etc)
    - category name: `response`
    - methods: `libtetrabz`(analytic tetrahedra integration), `sparse_ir` (sparse matsubaras with convolution theorem)
    - wpts: number of Matsubara frequencies to use. wpts=1 calculates iv=0
- Band structure/density of states/Fermi surface calculations and displays


