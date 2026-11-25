# Estimating silent logical failure in small-set-flip decoding

Reproduction code and numerical certificates for the paper by Vicenzo Scavino Alfaro. The complete archive, including the 1,412 binary strip-catalog fragments, is available on Zenodo.

[Paper (PDF)](manuscrito/Silent_logical_failure_in_small_set_flip_decoding.pdf) · [Reproduction guide](REPRODUCTION.md) · [Archived release on Zenodo](https://doi.org/10.5281/zenodo.24521353)

## Run the quick check

Use 64-bit Python 3.12 with assertions enabled. The quick check uses the standard library and runs from a fresh clone.

```sh
git clone --depth 1 https://github.com/vicenzoscavino1999/ssf-qec-reproduction.git
cd ssf-qec-reproduction
python reproducir.py rapido
```

On Windows, `py -3.12` also selects Python 3.12. The command verifies the SHA-256 hashes of the 1,383 files in the checkout. It then recomputes the Bernoulli intervals from the archived counts, regenerates 128 strip-table rows and checks 62 numerical values in Figures 4–6. The final line reports `OK` and the path to `resultado.json`.

For Table 2 at `p = 10^-6` and `q = 10^-4`, the deterministic bounds on the silent-failure probability, rounded outward, are

| Qubits | Lower bound | Upper bound |
| ---: | ---: | ---: |
| 225 | 3.24009 × 10⁻¹⁴ | 3.27319 × 10⁻¹⁴ |
| 900 | 8.63039 × 10⁻¹⁴ | 9.69644 × 10⁻¹⁴ |

## Rebuild the two censuses

Run `python reproducir.py sectores` to compile the C++17 decoder, build the two codes from their seeds and repeat the counts used in Table 2. The silent counts in the sector with two data faults and one readout fault are 324 and 864. The reproduction guide lists the compiler and CPU requirements.

## Reproduce the full strip runs

Download the archived ZIP and extract it into an empty directory. From the extracted directory, `python reproducir.py rapido` verifies all 2,796 files in the archive. The guide gives the commands, time and storage requirements for `tiras-1`, `tiras-2` and `catalogos`. The computational source files in this checkout match the archived copies byte for byte. Appendix F in the repository manuscript links to the archive and the source repository. Its SHA-256 manifest covers the files distributed here. [Download the full ZIP from Zenodo](https://zenodo.org/records/24521353/files/REPRODUCCION_SSF_QEC_20261124.zip?download=1) · [Full reproduction instructions](REPRODUCTION.md)

The SHA-256 of `REPRODUCCION_SSF_QEC_20261124.zip` is `4c4bf3b511714db87015952f4edb0b41100985aee681a513d5ed4412682c174c`.

## Citation and license

Cite the Zenodo release using DOI **10.5281/zenodo.24521353**. The author's original material is distributed under Creative Commons Attribution 4.0 International. [Citation metadata](CITATION.cff) · [License](LICENSE)
