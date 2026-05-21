# many_body / renormalization / analytic

## Overview

Calculates renormalization constant Z analytically based on the FLEX formula for self-energy. Takes in chi(w,q) as input and calculates either on q-mesh or on stored points, whichever is given by the chi(w,q) used for calculations

## Quick Description

Calculates quasiparticle weight Z approximating V(w)=V(0)

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
brillouin_zone

## Results Saved

- `{outdir}_{prefix}_renormalization.{ext}` - Saves renormalization as a function of k, on the same k-points as chi(w,q) input

## Testing

None

## Calculation Details

### Algorithm


### Implementation Notes

- Known to disagree with computational FLEX self-energy, perhaps Scalapino made a typo

## References

1. Scalapino et al., "d-wave pairing near a spin-density-wave instability", Phys. Rev. B. 34, (1986). 
