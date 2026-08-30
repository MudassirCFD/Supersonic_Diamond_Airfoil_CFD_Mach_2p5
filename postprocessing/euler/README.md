# Euler post-processing

This folder contains the Python post-processing used for the WENO5-HLLC Euler reference case.

Run from the repository root:

```bash
python postprocessing/euler/postprocess_weno_hllc.py \
    --dir results/euler/reference_raw \
    --out results/euler/reference/post
