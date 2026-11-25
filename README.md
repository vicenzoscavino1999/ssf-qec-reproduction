# Certified estimates of silent logical failure in small-set-flip decoding with noisy readout

Reproduction code and numerical certificates for the paper by Vicenzo Scavino Alfaro. The complete archive, including the 1,412 binary strip-catalog fragments, is available on Zenodo.

[Paper (PDF)](manuscrito/Silent_logical_failure_in_small_set_flip_decoding.pdf) · [Reproduction guide](REPRODUCTION.md) · [Archived release on Zenodo](https://doi.org/10.5281/zenodo.24521353)

## Run the quick check

Use 64-bit Python 3.12 with assertions enabled. The quick check uses the standard library and runs from a fresh clone.

```sh
git clone --branch circuit-postselection-v2.1 --depth 1 https://github.com/vicenzoscavino1999/ssf-qec-reproduction.git
cd ssf-qec-reproduction
python reproducir.py rapido
```

On Windows, `py -3.12` also selects Python 3.12. The command verifies the SHA-256 hashes of the files listed in the checkout manifest. It then recomputes the Bernoulli intervals from the archived counts, replays the two SSF trajectories in Figure 1, regenerates 128 strip-table rows and checks 62 numerical values in Figures 5–7. The final line reports `OK` and the path to `resultado.json`.

For Table 2 at `p = 10^-6` and `q = 10^-4`, the deterministic bounds on the silent-failure probability, rounded outward, are

| Qubits | Lower bound | Upper bound |
| ---: | ---: | ---: |
| 225 | 3.24009 × 10⁻¹⁴ | 3.27319 × 10⁻¹⁴ |
| 900 | 8.63039 × 10⁻¹⁴ | 9.69644 × 10⁻¹⁴ |

## Estimate at your own Bernoulli rates

From this checkout, choose the 225- or 900-qubit code and enter independent data and readout rates. The command evaluates the interval from Corollary 6 and reports a point estimate with a certified relative-error bound.

```sh
python reproducir.py regional --qubits 225 --p 1e-6 --q 1e-4
```

The region, boundary cases and output fields are described in the [regional calculation](REPRODUCTION.md#evaluate-the-regional-interval).

## Rebuild the two censuses

Run `python reproducir.py sectores` to compile the C++17 decoder, build the two codes from their seeds and repeat the counts used in Table 2. The silent counts in the sector with two data faults and one readout fault are 324 and 864. The reproduction guide lists the compiler and CPU requirements.

Run `python reproducir.py marca-tasas` to rebuild the one-round abort-flag bounds on 225 and 900 qubits and the finite Bernoulli interval at `p=1e-4, q=1e-3` on 225 qubits. The quick check evaluates the same bounds from the stored counts. [Flag and high-rate reproduction](REPRODUCTION.md#reproduce-the-abort-flag-and-higher-bernoulli-rates)

Run `python reproducir.py marca-circuito` for the eight deterministic circuit postselection intervals on 225 qubits. At the VV-first point with operation and readout rates of `2e-7` and `1e-5`, fewer than 48 discards per million rounds give at least a 257-fold reduction in silent-failure risk per retained round. [Circuit postselection](REPRODUCTION.md#circuit-postselection-on-225-qubits)

## Reproduce the full strip runs

Download the archived ZIP and extract it into an empty directory. From the extracted directory, `python reproducir.py rapido` verifies all 2,796 files in the archive. The guide gives the commands, time and storage requirements for `tiras-1`, `tiras-2` and `catalogos`. The archived release contains the original manuscript and its quick check. This checkout includes the opening-figure replay, the regional calculator, the postselection and high-rate experiments, and the revised manuscript. Appendix F in the repository manuscript links to the archive and the source repository. Its SHA-256 manifest covers the files distributed here. [Download the full ZIP from Zenodo](https://zenodo.org/records/24521353/files/REPRODUCCION_SSF_QEC_20261124.zip?download=1) · [Full reproduction instructions](REPRODUCTION.md)

The SHA-256 of `REPRODUCCION_SSF_QEC_20261124.zip` is `4c4bf3b511714db87015952f4edb0b41100985aee681a513d5ed4412682c174c`.

## Citation and license

Cite the Zenodo release using DOI **10.5281/zenodo.24521353**. The author's original material is distributed under Creative Commons Attribution 4.0 International. [Citation metadata](CITATION.cff) · [License](LICENSE)
