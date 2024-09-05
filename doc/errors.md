
# Common mistakes

### Common Errors

```
(plo + ihi*CellSize(idim)) < (plo + (ihi + 1)*CellSize(idim))'
```

This (may) happens when xmin = xmax in a particular dimensions, even if that dimension is not used
(for example a 2D run but using 3D compile options).
Check the lines in inputs

```
geometry.prob_lo     = 0.0    0.0    0.0
geometry.prob_hi     = 1.0    1.0    1.0
```