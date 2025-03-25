### ** 🔹 Module Objects **
- CMData
    - Accepts a set of points and values for reading and writing files
    - Data saved in the format x, y, z, w, n, f (header listed)
    - w and n columns are specified on declaration
    - Can be called with a file, where it will read the header and create the necessary variables automatically

- CMField
    - Takes CMData as a component for data reading/writing
    - Reads x, y, z as a mesh grid used for interpolation
    - w values are interpolated via binary search
    - n values are given as input if the n column exists
