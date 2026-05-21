# many_body / self_energy / sparse_ir

## Overview

Calculates Self-Energy from Vertex and non-interacting green's function. Takes the Vertex as input, all calculations are on imaginary axis

## Quick Description

Calculates Self-Energy from given Vertex and non-interacting green's function

## Dependencies
- sparse_ir
- FFTW
- PyCall
- Printf
- Interpolations

## Install Instructions

```bash
```
Install relevant packages using julia repl

### Parameters

Configuration parameters from `input.cfg`:

- `prefix` : Prefix for input/output files
- `outdir` : Directory for input/output files
- `T` : Temperature
- `dim` : Dimensionality of the system (1, 2, or 3)
- `q_mesh` : Number of q-points in each dimension
- `k_mesh` : Number of k-points in each dimension
- `nbnd` : Number of bands
- `fermi_energy` : Fermi energy
- `brillouin_zone` : BZ matrix


## Results Saved

- `{outdir}_{prefix}_.self_energy.{filetype}` - self-energy on iw-k grid

## Testing

None

## Calculation Details

Performs FFTs to compute self-energy from vertex function and non-interacting Green's function on imaginary frequency axis. Matsubara frequencies are handled using sparse IR basis for efficiency, where imaginary frequencies are transformed to Discrete Lehman Representation and then to imaginary time for computing fourier transforms.

### Algorithm

1. Find G from non-interacting Hamiltonian, V from input
2. Compute G(r, tau) and V(r, tau) using sparse IR basis and FFTs
3. Calculate self-energy in (r, tau) space
4. Transform self-energy back to (k, iw) space using inverse FFTs and sparse IR basis

### Implementation Notes

- Single band 

## References

1. https://spm-lab.github.io/sparse-ir-tutorial/src/FLEX_jl.html
