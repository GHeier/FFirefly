# many_body / vertex / from_susceptibility

## Overview

Analytically calculates FLEX vertex from chi(q,w) using U and the w-q points in chi(w,q)

## Quick Description

Analytically calculates FLEX vertex from chi(q,w)

## Dependencies
- None

## Install Instructions

```bash
No special installation required - built automatically by fly-build.sh
```

### Parameters

outdir
prefix
filetype
U0 - U strength

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}_{prefix}_vertex.{ext}` - Saved vertex file in same w-k format as chi file

## Testing

None

## Calculation Details

### Algorithm


### Implementation Notes

- Only saves as hdf5 file at the moment

## References

1. Scalapino et al., "d-wave pairing near a spin-density-wave instability", Phys. Rev. B. 34, (1986). 
