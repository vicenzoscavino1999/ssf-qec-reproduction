# Reproduce the silent logical failure estimates

## 1. Reproduce the intervals in Table 2

Use 64-bit Python 3.12 with assertions enabled (on Windows, you can use `py -3.12` for `python`). If you have the full `REPRODUCCION_SSF_QEC_20261125.zip`, extract it into an empty directory and leave its nested ZIPs intact. Allow 15 GB for download and extraction together. From the extracted folder or the compact source copy, open a terminal beside `reproducir.py` and run the quick check. It uses only the Python standard library.

```sh
python reproducir.py rapido
```

At independent data and readout fault rates `p = 1e-6` and `q = 1e-4`, the command reconstructs Table 2 from the stored counts. The expected bounds on silent logical failure are rounded outward.

| Code | Lower bound | Upper bound |
| --- | ---: | ---: |
| 225 qubits | 3.24009 × 10⁻¹⁴ | 3.27319 × 10⁻¹⁴ |
| 900 qubits | 8.63039 × 10⁻¹⁴ | 9.69644 × 10⁻¹⁴ |

The final line prints `OK` and the path to the run summary, `resultado.json`. In that directory, `01-calculos_del_articulo.log` gives both intervals at full output precision and `marca-tasas.json` gives the Bernoulli flag bounds from the stored counts. The quick check also builds the strip tables, replays both SSF trajectories in Figure 1, verifies the negative zero in Appendix B.3 and checks the values plotted in Figures 5–7.

Each run creates a directory under `resultados_ssf/` beside the source directory. To choose where it goes, pass `--salida` followed by a new directory outside the source tree.

## 2. Choose a result to recompute

Run `python reproducir.py ACTION` with the action in the middle column. Output paths below are relative to the run directory.

| Result in the manuscript | Action | Expected output |
| --- | --- | --- |
| Family-wide noise rates, Table 1 and Appendices B–C | `tasas` | The three (3,4) rate rectangles, the transfer from data/readout to circuit noise, and the (17,18) certificate. The twelve ratios for the second rectangle are in `06-verificar_H5A_3_4_esquina_2.log`. [Derivation](certificados/fibras/proofs/H5A_3_4_esquina_2.txt) |
| Distance-six codes, Table 2 and Figure 5(b) | `sectores` | `censo/resultado_225.json` and `censo/resultado_900.json` give the counts from the two code seeds. Two data faults and one wrong readout give 324 silent inputs on 225 qubits and 864 on 900. |
| Bernoulli flag, Table 3 and Appendix E.1.3. Higher rates, Appendix E.2.1 | `marca-tasas` | Repeats the census and orbit counts. `marca-tasas.json` gives the discard and retained-risk bounds. At `(p,q) = (1e-4,1e-3)` the 225-qubit estimate has relative error below 5.396%. [Output fields](sectores/flag_rates/README.txt) |
| Circuit flag, Table 3 and Appendix E.1.4 | `marca-circuito` | `marca-circuito.json` gives all eight circuit rows. At `(r_c,r_m) = (2e-7,1e-5)` with VV first, fewer than 48 rounds per million are discarded and the silent-failure risk per retained round falls by a factor of at least 257. |
| Strip cutoff and resource comparisons, Tables 4, 7 and 8 | `rapido` | In `tablas_tiras/`, `circuit_rows.tex` and `phenomenological_rows.tex` give the table entries. `resources.json` gives the recorded times and memory. |
| Affine F17 instance, Table 5 and Appendix D | `f17` | Logical masks, implicit comparisons, 1,214 circuit traces in `certificate.jsonl`, and sufficient truncation orders 253,000 and 126,000. [Finite-sector inputs](certificados/sector_f17/README.txt) |
| 225-qubit pilot, Table 6 and Figure 5(a) | `piloto` | `piloto/piloto_HGP_225_esquinas_20261124.json` contains three deterministic Bernoulli intervals and eight stratified circuit intervals. The latter have simultaneous coverage of at least `1 − 16 exp(−7)` under independent sampling within each stratum. |
| PATH25, Appendix E.3 | `path25` | Two decoder implementations agree on 4,096 syndromes. On the test circuit, enumeration covers 524,288 assignments and the component cutoff is complete at eight noisy locations. [Estimator outputs](path25/estimadores/README.txt) |
| Rare-event comparisons, Tables 9–10 and Appendix E.4 | `comparaciones` | Recomputes the estimates and mixing diagnostics from stored samples. The 61-qubit reference is deterministic. The 225-qubit circuit comparison includes both deterministic and stratified reference intervals. |
| Manuscript and Figures 1–7 | `pdf` | Compiles the PDF with three pdfLaTeX passes. The figures and bibliography are in the TeX source. |

**Software.** The reference runs used Windows x64 and Python 3.12.10. Use Windows for `path25`, the 61-qubit chain checks in `comparaciones`, both splitting runs and strip production. Four actions compile C++17 code. On Windows, `sectores`, `marca-tasas`, `dss-225` and `splitting-225` use Visual Studio x64 Build Tools on a CPU supporting AVX2 and BMI2. On x86-64, `sectores`, `marca-tasas` and `dss-225` also support GNU C++. The PDF requires pdfLaTeX, TikZ, PGFPlots and the other packages in the [preamble](manuscrito/Silent_logical_failure_in_small_set_flip_decoding.tex).

## 3. Evaluate another pair of Bernoulli rates

```sh
python reproducir.py regional --qubits 900 --p 4e-7 --q 5e-5
```

Choose `--qubits 225` or `900`. For positive rates, Corollary 6 requires `q <= 1e-4` and `z_p/z_q <= 9999/999999`, where `z_x = x/(1-x)`. You can enter a decimal, fraction or scientific notation. When `p = 0` and `0 <= q <= 1e-4`, the probability is zero.

For the example above, the interval is approximately `[6.7622e-15, 7.9308e-15]` (rounded outward). In `resultado.json`, `regional.lower` is the two-data, one-readout mass rounded down. The upper endpoint allows a remainder of up to 1.80% of that mass on 225 qubits or 17.28% on 900. `regional.estimate` gives the harmonic estimate and `regional.max_relative_percent` bounds its relative error (`null` at exact zero). At the rates of Table 2, the table uses additional sectors, so its endpoints can differ from the regional calculation.

## 4. Repeat the sampling and strip runs

`dss-225` repeats one million draws with replacement from the two-data, one-readout stratum, then enumerates all 2,721,600 inputs in that stratum. In `dss225.json` the reference sample has 119 silent failures and the census has 324. The sampling commands use the recorded seeds. To regenerate the splitting chains, run `splitting-225` or `splitting-61`. The recorded runs took 161.91 s on 225 qubits (28 rate levels) and 2.49 h on 61 qubits.

For strip production, use the full ZIP with its 1,412 catalog fragments on Windows x64. The runs allow 24 GiB of aggregate memory. A 32 GiB machine is recommended. To verify the recorded runs, use `python tiras/reproduce.py verify-recorded`. This command includes the catalog checks. To recompute the intervals, run the two commands below in order. The second reads the output directory of the first.

```powershell
python reproducir.py tiras-1 --salida ../ssf-results/strip-1
python reproducir.py tiras-2 --previa ../ssf-results/strip-1 --salida ../ssf-results/strip-2
```

The first run took 3.32 h on the reference machine, with a 6.5 h budget. The second adds the 523-qubit strip and has a 9 h budget. Each run allows up to 60 GB of output. To regenerate the catalogs themselves, use `python reproducir.py catalogos` (8 h budget, 24 GiB aggregate RAM, 60 GB output). [Strip files](tiras/README.txt)

## Supplements cited in Appendix F

| Supplement | Location in this distribution |
| --- | --- |
| `SOURCE_CUMULANT_CERTIFICATE.zip` | `certificados/f17/` |
| `FIBRE_CERTIFICATES.zip` | [Circuit fibers and first contact](certificados/fibras/README.txt) in `certificados/fibras/` |
| `REPRODUCCION_PAPER_20261117_v3.zip` | `certificados/supplement/h1/` and `certificados/supplement/h3/` |
| `F17_CIRCUIT_SECTOR_20261107.zip` | `certificados/sector_f17/` |
| `SSF_SECTORES_225_900_20261124.zip` | `sectores/` |
| `HGP225_PILOT.zip` | `comparaciones/` |
| `QEC_RARE_EVENT_COMPARISONS_20261124.zip` | `comparaciones/` and the 61-qubit run records in `tiras/workspace/` |
| `PATH25_VALIDATION_20261106.zip` and `PATH25_CERTIFIED_ESTIMATORS_20261112.zip` | `path25/validacion/` and `path25/estimadores/`, respectively |
| `STRIP_COMPONENT_CERTIFICATES_20261116.zip` and `T2_1_REPRODUCCION_COMPLETA_20261116.zip` | `tiras/tables/` and `tiras/`, respectively |

Version 3 of the full package contains the manuscript, the launcher and the strip catalogs. The compact source is tagged `circuit-postselection-v3`. [Zenodo archive](https://doi.org/10.5281/zenodo.24542801) · [GitHub release](https://github.com/vicenzoscavino1999/ssf-qec-reproduction/releases/tag/circuit-postselection-v3)
