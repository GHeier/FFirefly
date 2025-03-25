# Module Objects
## CMData
### Overview
Class for data input, output, and storage. Data is stored in columns with automatically generated headers of the form # x y z w n f. x y z are the grid points, w the 4th dimensional point, n the index indicating which function it is, and f the result of the function. f may also be Re(f) Im(f) for complex values, fx fy fz for vector values, and Re(fx) Im(fx) ... for complex vector values.

---

### Public Member Variables

    - points (vector<Vec>)
    - w_points (vector<float>)
    - 
- Accepts a set of points and values for reading and writing files
- Data saved in the format x, y, z, w, n, f (header listed)
- w and n columns are specified on declaration
- Can be called with a file, where it will read the header and create the necessary variables automatically

- CMField
    - Takes CMData as a component for data reading/writing
    - Reads x, y, z as a mesh grid used for interpolation
    - w values are interpolated via binary search
    - n values are given as input if the n column exists
