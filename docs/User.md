### 🔹 **Input Variables**

---
#### 🔸 CONTROL
    category - options listed below
    method - options listed below
    outdir - where data is saved
    indir - where data is read
    prefix - beginning of all relevant input/output files
    verbosity - low, high(default)

#### 🔸 SYSTEM
    interaction - vertex function. Options: (FLEX, EPC)
    Temperature - units of eV
    cell_type - Simple Cubic, Body Centered Cubic, Face Centered Cubic, Orthorhombic, Tetragonal, Hexagonal
    U0 - Hubbard Repulsion value
    nbnd - Number of bands
    dimension - Dimension of system

#### 🔸 MESH
    k_mesh - Number of points in kx, ky, kz. Used as input for calculations on meshes (ie $E(k)$)
    q_mesh - Number of points in qx, qy, qz. Used for mesh output (ie $\chi(q)$)
    w_pts - Number of points in frequency. Used for Matsubara frequencies as well

#### 🔸 CELL
    3x3 lattice vectors (implicitly defined if given BRILLOUIN_ZONE)

#### 🔸 BRILLOUIN_ZONE
    3x3 reciprocal lattice vectors (implicitly defined if given CELL)

#### 🔸 BANDS

    FORMAT:

    nbnd bands are listed in the format of:
        band1 = 'band_name'
            var = value
        band2 = 'band_name'
            var = value
        ...
    
    VARIABLES:

    band_names - 'fermi_gas', 'tight_binding'
        fermi_gas vars
            eff_mass
            shift
        tight_binding vars
            t0
            t1
            ...
            t10

### 🔹 **Current Calculation Categories**

---

#### 🔸 `test` *(default)*
- **Purpose**: Runs basic tests to confirm that all categories are functioning correctly.

---







#### 🔸 `hamiltonian`

- **DOS**
  - [gaussian](src/hamiltonian/DOS/gaussian/README.md) - Computes the Density of States with gaussian spreading for smooth results.
  - [tetrahedra](src/hamiltonian/DOS/tetrahedra/README.md) - Computes the Density of States using surface construction at discrete w-points.
- **FS**
  - [tetrahedra](src/hamiltonian/FS/tetrahedra/README.md)
- **generate**
  - [band_structure](src/hamiltonian/generate/band_structure/README.md) - Computes Hamiltonian, H(k) based on an H(r) tight-binding construction
  - [hk_from_hr](src/hamiltonian/generate/hk_from_hr/README.md) - Computes Hamiltonian, H(k) based on an H(r) tight-binding construction

---

#### 🔸 `many_body`

- **many_body**
  - [sparse_ir](src/many_body/many_body/sparse_ir/README.md) - Performs FLEX calculations using DLR sparse_ir package
  - [triqs](src/many_body/many_body/triqs/README.md) - Solves FLEX or FLEX+DMFT with DLR calculations at finite Temperature
- **renormalization**
  - [FS_approx](src/many_body/renormalization/FS_approx/README.md) - Computes quasiparticle weight Z(k) across the Fermi Surface using V(w)=V(0) approximation.
  - [analytic](src/many_body/renormalization/analytic/README.md) - Calculates quasiparticle weight Z approximating V(w)=V(0)
  - [from_sigma](src/many_body/renormalization/from_sigma/README.md) - Calculates Z(k) based on the slope of a given Sigma(iω,k) at ω→0
- **response**
  - [sparse_ir](src/many_body/response/sparse_ir/README.md) - Computes the non-interacting response function by via Green's function convolution
  - [tetrahedra](src/many_body/response/tetrahedra/README.md) - Calculates non-interacting response function chi0(w,q) using recursive tetrahedron integration
- **self_energy**
  - [sparse_ir](src/many_body/self_energy/sparse_ir/README.md) - Calculates Self-Energy from given Vertex and non-interacting green's function
  - [triqs](src/many_body/self_energy/triqs/README.md) - Calculates the self-energy using Iterated Perturbation Theory (IPT) on the imaginary axis.
- **vertex**
  - [from_susceptibility](src/many_body/vertex/from_susceptibility/README.md) - Analytically calculates FLEX vertex from chi(q,w)

---

#### 🔸 `superconductor`

- **bcs**
  - [convolution](src/superconductor/bcs/convolution/README.md) - Solves the linearized BCS gap equation across the Brillouin Zone
  - [matrix](src/superconductor/bcs/matrix/README.md) - Solves the linearized BCS gap equation on the Fermi Surface using standard matrix diagonalization.
- **bcs_w**
  - [hmatrix](src/superconductor/bcs_w/hmatrix/README.md) - Solves the linearized BCS gap equation across the Fermi Surface using compressed Hierarchical Matrices and a lanczos matrix solver.
- **eliashberg**
  - [convolution](src/superconductor/eliashberg/convolution/README.md) - Uses Lanczos eigensolver to find the leading gap solutions in full Brillouin Zone
  - [hmatrix](src/superconductor/eliashberg/hmatrix/README.md) - Solves Eliashberg equation on real axis using HMatrix compression and Lanczos solver.
  - [sparse_ir](src/superconductor/eliashberg/sparse_ir/README.md) - Solves the linearized Eliashberg equation on imaginary axis using convolution and the power iteration / Krylov projection approach.

---
