# PocketFilter: Step 2 software reference

Software baseline for later AUP-ZU3 / Vitis HLS verification. Adapted from
[TimD1/SquiggleFilter](https://github.com/TimD1/SquiggleFilter),
`sdtw_analysis.ipynb`. Credit: Dunn et al., *SquiggleFilter: An Accelerator for
Portable Virus Detection*, MICRO-54, 2021. This is an adaptation, not a new DTW
algorithm or a completed FPGA accelerator. No new license is asserted over
upstream code; preserve attribution and verify upstream redistribution terms.

## Useful code

| File | Purpose |
|---|---|
| golden_model.py | Callable sDTW cost and threshold classification; Numba acceleration |
| smoke_test.py | 30 random comparisons against an independently written dense recurrence |
| run_real_data.py | Reference generation, FAST5 preprocessing, scoring and vector export |
| upstream_sdtw_analysis.ipynb | Original source, retained for provenance and preprocessing helpers |
| reported_initial_scores.csv | The 20 scores supplied by the student from the Windows run |

`run_real_data.py` packages the previously executed PowerShell/Python workflow;
it is newly organized here, not a claim of a separate real-data rerun. It adds
read-ID deduplication. The original `original_core.py` excerpt is unnecessary
for this runner. Do not run the entire upstream notebook for this workflow.

## Changes relative to upstream

- Explicit reference/query/length arguments replace notebook-global algorithm inputs.
- Preserves the original bonus=10, run-count cap, recurrence and initialization;
  this is not an interchangeable textbook DTW implementation.
- Corrects output `min_cost[l]` to `min_cost[l-1]` for an L-sample prefix.
  Original full-length indexing can go out of bounds; results are not claimed
  bit-identical to the unmodified notebook.
- Adds input validation, zero-deviation handling, a keep/eject interface and exports.
- Real-data runner preserves original stall trimming, skips the following 1000
  samples, and uses cumulative-prefix normalization at 500-sample boundaries.
- Preserves original reference-generation behavior, including trailing padding.
- Uses integer-valued normalized signals, but floating-point NumPy cost arrays.
  Hardware bit widths and bit-accurate integer behavior remain Step 3 work.
- Original preprocessing is loaded from the included notebook; this is not yet
  a fully notebook-independent end-to-end implementation.

## Run on Windows

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pip check
.\.venv\Scripts\python.exe smoke_test.py
.\.venv\Scripts\python.exe run_real_data.py
```

Place the original reference and data at:

- data/dna_kmer_model.txt
- data/covid/reference.fasta
- data/covid/fast5/*.fast5
- data/human/fast5/*.fast5

Reference/model files are under the upstream repository's data directory.
Raw signal sources are documented in upstream setup.sh. They are not bundled.
The student's downloaded samples can be reused without downloading again.
Copy this package's files into a separate folder or back up existing results
before rerunning: the runner overwrites reference.npy, initial_scores.csv and
numbered vector files. Set MAX_READS to expand the run, with disjoint read IDs
for calibration and testing. Dependencies are not a frozen reproduction of
the student's environment; save `pip freeze` from the working Windows venv.

## Observed results (student-reported)

- Smoke test: PASS for all 30 random comparisons, including full-length boundary.
- Reference: 29,903 bases, 4,096 k-mers, 59,806 normalized samples, range 0–216.
- Query prefix: 1,000 samples; 10 COVID-dataset and 10 human-dataset reads.

| Dataset label | n | Min cost | Median cost | Max cost |
|---|---:|---:|---:|---:|
| COVID | 10 | 2941 | 3812.5 | 8830 |
| Human | 10 | 5532 | 6230 | 7984 |

Scores and IDs in reported_initial_scores.csv are transcribed from the student's
console output, not independently recomputed here. Original FAST5 files and
NPZ test vectors remain on the student's machine. Labels indicate dataset
origin; per-read biological identity has not been independently established.

The lower COVID median is encouraging, but the distributions overlap. No final
threshold, held-out accuracy, real-time throughput or power claim is made.
Synthetic recurrence agreement tests implementation consistency, not biological
accuracy or independent correctness of the upstream algorithm.

## Step 3

Translate the chosen recurrence into C++ for Vitis HLS; compare its cost against
saved vectors; select signed integer widths with overflow checks; run C/RTL
co-simulation and synthesis; then integrate with PYNQ/DMA. Measure resources
and throughput before claiming the 2M-samples/s project goal. Keep algorithm
and preprocessing versions fixed across software and hardware comparisons.

## Upload

Create an empty GitHub repository and upload this folder's files, preserving
attribution. Do not upload `.venv`, raw FAST5 datasets or downloaded archives.
Include selected local NPZ vectors and reference.npy separately if desired;
their hashes and provenance should accompany them. Replace the email's URL
placeholder with the actual repository link after uploading.
