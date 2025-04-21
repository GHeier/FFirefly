### 🔹 **Module Objects**

---

#### 🔸 `Bands`

- **Purpose**: Provides band structure data from `_bands.dat`, or falls back to the `.cfg` definition (`epsilon` function) if the file is not present.
- **Behavior**:
  - Automatically reads `_bands.dat` upon initialization.
  - Defaults to `.cfg` if no file is found.
  - Internally wraps a `Field_R` object.

- **Usage**:
  - **C++**:
    ```cpp
    Bands band;
    band(n, k);
    ```
  - **Python**:
    ```python
    band = Bands()
    band(n, k)
    ```
  - **Parameters**:
    - `n`: *integer* – Band index (starting at **1**)
    - `k`: *vector<float>* – 1D, 2D, or 3D momentum vector

---

#### 🔸 `Vertex`

- **Purpose**: Returns the two-particle vertex $\Gamma$ using data from `_2PI.dat` or the fallback `interaction` defined in `.cfg`.
- **Behavior**:
  - Automatically loads `_2PI.dat` if available.
  - Falls back to `.cfg` definition if `_2PI.dat` is missing and the interaction is analytically defined.
  - Accepts optional label arguments (e.g., spin, valley).

- **Usage**:
  - **C++**:
    ```cpp
    Vertex V;
    V(q, w, label1, label2);
    ```
  - **Python**:
    ```python
    V = Vertex()
    V(q, w, label1, label2)
    ```
  - **Parameters**:
    - `q`: *vector<float>* – Momentum transfer vector  
    - `w`: *float* – Energy transfer (pass as imaginary value for Matsubara)  
    - `label1`, `label2`: *optional* – Labels for degrees of freedom (e.g., spin, orbital)

---

#### 🔸 `Field_R`

- **Purpose**: Wraps a `CMField` and returns **real-valued** scalar field data.

- **Usage**:
  ```cpp
  Field_R field;
  field(x, y, z);

#### 🔸 `Field_C`

- **Purpose**: Wraps a `CMField` and returns **complex-valued** scalar field data.

- **Usage**:
  ```cpp
  Field_C field;
  field(x, y, z);

#### 🔸 `Field_VR`

- **Purpose**: Wraps a `CMField` and returns **real-valued** vector field data.

- **Usage**:
  ```cpp
  Field_VR field;
  field(x, y, z);

#### 🔸 `Field_VC`

- **Purpose**: Wraps a `CMField` and returns **complex-valued** vector field data.

- **Usage**:
  ```cpp
  Field_VC field;
  field(x, y, z);

#### 🔸 `CMData`

- **Purpose**: Manages structured datasets and handles reading/writing from disk with a defined column format.
- **Behavior**:
  - Accepts a set of data points (typically from a file).
  - Data format:  
    ```
    x, y, z, w, n, f
    ```
    where the header defines each column.
  - `w` and `n` columns must be specified during initialization.
  - If initialized with a file, it will automatically:
    - Read the header
    - Create and map internal variables accordingly

- **Usage**:
  - **C++**:
    ```cpp
    CMData data("input.dat");
    CMData data(x, y, z, f);
    data.save(filename);
    ```
  - **Python** (if applicable):
    ```python
    data = CMData("input.dat")
    data = CMData(x, y, z, f)
    data.save(filename);

- Save Format:
    x            y         n         f    
-1.000000    -1.000000     1     -2.000000      
    ...
    ```

---

#### 🔸 `CMField`

- **Purpose**: Performs structured field operations such as interpolation, based on data loaded via `CMData`.
- **Behavior**:
  - Uses `x`, `y`, `z` columns from `CMData` to build a regular mesh grid.
  - Interpolates `w` values using binary search over the grid.
  - Uses `n` values (if available) as supplementary input.
  - Automatically assumes **periodic boundary conditions**: queries outside the defined grid wrap around rather than erroring.
  - Designed for efficient and safe high-dimensional lookups.

- **Usage**:
  - **C++**:
    ```cpp
    CMField field("file.dat");
    auto value = field(x, y, z);
    ```
  - **Python** (if supported):
    ```python
    field = CMField("file.dat")
    value = field(n, x, y, z, w)
    ```

  - **Parameters**:
    - `x`, `y`, `z`: *float* – Spatial coordinates used for mesh-based interpolation
    - `n`, `w`: *int, float* – Index and 4th dimension coordinate used for binary-search interpolation

---
