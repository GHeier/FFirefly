# many_body / renormalization / from_sigma

## Overview

Calculates Z(k) based on the slope of Sigma(iω,k) at ω→0, using data from a previous self-energy calculation. Saves Z(k) across the points given in Sigma(iw,k)

## Quick Description

Calculates Z(k) based on the slope of Sigma(iω,k) at ω→0, using data from a previous self-energy calculation.

## Dependencies
- None

## Install Instructions

```bash
None
```

### Parameters

`k_mesh` - Number of k-points in each dimension (int or list of int)
`brillouin_zone` - Matrix defining the Brillouin zone 

## Results Saved

- `{outdir}_{prefix}_renormalization.{ext}` - Description of what this file contains

## Testing

None

## Calculation Details

### Algorithm

Description of the algorithm:
1. Step 1
2. Step 2
3. etc.

### Implementation Notes

- Any important implementation details
- Performance considerations
- Known limitations

## References

1. Author et al., "Paper Title", Journal Volume, Pages (Year). DOI/arXiv
2. Additional references as needed
