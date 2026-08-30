# Euler post-processing

This folder contains the Python post-processing used for the WENO5-HLLC Euler reference case.

Run from the repository root:

```bash
python postprocessing/euler/postprocess_weno_hllc.py     --dir results/euler/reference_raw     --out results/euler/reference/post
```

Add `--animate` only when snapshot CSV files are available.

The scripts produce field contours, numerical schlieren, shock-angle checks,
surface pressure comparisons, force and residual histories, analytical reference
plots and compact CSV extracts.
