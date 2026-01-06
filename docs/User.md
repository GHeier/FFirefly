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

#### 🔸 `superconductor`
- **Purpose**: Calculates superconducting gap size and symmetry.
- **Options**:
  - `method`:  
    - `bcs` – Basic BCS approximation  
        - Outputs `_gap.dat` file with gap values along surface
    - `eliashberg` – Eliashberg theory
        - Outputs `_gap.dat` file with gap values across brillouin zone
  - `calculation`:  
    - `diagonalization` – Direct diagonalization of gap functions  
    - `projection` – Projection onto a predefined basis set
  - `FS_only`:  
    - `true` – Only compute gap on the Fermi surface  
    - `false` – Compute over the entire Brillouin zone

---

#### 🔸 `response`
! **Purpose**: Computes bare susceptibility.
! **Options**:
  ! `method`:  
    ! `sparse_ir` – Uses sparse Matsubara frequencies (finite T; requires dense k-grid)  
        - Outputs `_chi.dat` file with chi values across brillouin zone
    ! `libtetrabz` – Uses analytic tetrahedra integration (0 K; allows sparse grid)  
        - Outputs `_chi.dat` file with chi values across brillouin zone
        ! ⚠️ *If you encounter "STOP NESTING", adjust the k-grid size.*
  ! `dynamic`:  
    ! `true` – Calculate and save data for finite $\omega$  
    ! `false` – Only compute at $\omega = 0$
  ! `wpts`:  
    ! Number of Matsubara frequencies.  
    ! `wpts = 1` calculates only at $i\omega = 0$.

---

#### 🔸 `vertex`
- **Purpose**: Calculates a specific two-particle vertex function.
- **Options**:
  - `interaction`:  
    - `FLEX` – Computes the FLEX vertex; requires bare susceptibility file (`_chi.dat`)
- Outputs `_vertex.dat` file with vertex values across brillouin zone

---

#### 🔸 `DOS`
- **Purpose**: Computes the Density of States.
- **Options**:
  - `method`:  
    - `libtetrabz` – Analytic tetrahedra integration  
    - `surface_sum` – Alternative method with comparable results
  - w_pts:
    - Number of energy points across bandwidth
  - k_mesh
    - Density of points when defining tetrahedra/surface points
- Outputs `_DOS.dat` file with Density of States vs Energy
