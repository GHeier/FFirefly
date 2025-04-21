### ** 🔹 Module Objects **
- Bands
    - Reads in "_bands.dat" automatically, and upon calling the object returns the energy at index n and vector k. If no "_bands.dat" is found, it defaults to the bands defined in .cfg (aka the epsilon function)
    - Wraps Field_R
    - use: Bands band; band(n, k) 
        - n:integer is the band index. Index starts at 1
        - k:vector{float} is the 1, 2, or 3D momentum vector

- Vertex
    - Reads in "_2PI.dat" automatically, and upon calling the object returns $\Gamma$ at vector k. If no "_2PI.dat" is found, it defaults to the vertex defined as `interaction` in .cfg, if no heavy calculation is required.
    - Accepts label arguments as well. (spin, valley, etc)
    - use: Vertex V; V(q, w, label1, label2)  
        - q:vector{float} momentum transfer vector
        - w:float energy transfer value. If working with matsubara, pass in the imaginary value here

- Field_R
    - Wraps CMField, returning real values

- Field_C
    - Wraps CMField, returning complex values

- Field_VR
    - Wraps CMField, returning real vector values

- Field_VC
    - Wraps CMField, returning complex vector values

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
    - Assumes periodicity, so there are no "out of bounds" for interpolation
