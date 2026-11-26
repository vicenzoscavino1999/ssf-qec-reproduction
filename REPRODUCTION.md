# Reproduce the silent logical failure estimates

## 1. Check the archived results

Use software release **4** (Git tag `circuit-postselection-v4`), cited in
Appendix F of the [manuscript](manuscrito/Silent_logical_failure_in_small_set_flip_decoding.pdf).
Download `ssf-qec-reproduction-v4.zip` from [Zenodo release 4](https://doi.org/10.5281/zenodo.24548710) or use the tagged Git checkout. Both include the 225-qubit postselection comparisons and their reference data. Section 4 covers the large catalogues for the full strip runs.

Open a terminal beside `reproducir.py` and run these two commands with
64-bit Python 3.12 and assertions enabled. Both use the standard library.
On Windows, `py -3.12` can replace `python`.

```sh
python reproducir.py rapido
python reproducir.py revision-v2
```

`rapido` reconstructs Table 2 from the stored counts at independent data and readout fault rates `p = 1e-6` and `q = 1e-4`. The expected bounds on silent logical failure are rounded outward.

| Code | Lower bound | Upper bound |
| --- | ---: | ---: |
| 225 qubits | 3.24009 × 10⁻¹⁴ | 3.27319 × 10⁻¹⁴ |
| 900 qubits | 8.63039 × 10⁻¹⁴ | 9.69644 × 10⁻¹⁴ |

Each command prints `OK` and the path to its own `resultado.json`. In the `rapido` directory, `01-calculos_del_articulo.log` gives both intervals at full output precision and `marca-tasas.json` gives the Bernoulli flag bounds. This check also builds the strip tables, replays the two SSF trajectories in Figure 1, verifies the negative zero in Appendix B.3 and checks the values plotted in Figures 5–7. In the `revision-v2` directory, `04-comparison_table.log` gives the nine rows of Table 6 assembled from the archived numerical results. Section 2 gives the commands to recompute them.

Each run creates a directory under `resultados_ssf/` beside the source directory. To choose where it goes, pass `--salida` followed by a new directory outside the source tree.

## 2. Choose a result to recompute

Run `python reproducir.py ACTION` with the action in the middle column. Output paths below are relative to the run directory.

| Result in the manuscript | Action | Expected output |
| --- | --- | --- |
| Family-wide noise rates, Table 1 and Appendices B–C | `tasas` | The three (3,4) rate rectangles, the transfer from data/readout to circuit noise, and the (17,18) certificate. The twelve ratios for the second rectangle are in `06-verificar_H5A_3_4_esquina_2.log`. [Derivation](certificados/fibras/proofs/H5A_3_4_esquina_2.txt) |
| Distance-six codes, Table 2 and Figure 5(b) | `sectores` | `censo/resultado_225.json` and `censo/resultado_900.json` give the counts from the two code seeds. Two data faults and one wrong readout give 324 silent inputs on 225 qubits and 864 on 900. |
| Bernoulli flag, Table 3 and Appendix E.1.3. Higher rates, Appendix E.2.1 | `marca-tasas` | Repeats the census and orbit counts. `marca-tasas.json` gives the discard and retained-risk bounds. At `(p,q) = (1e-4,1e-3)` the 225-qubit estimate has relative error below 5.396%. [Output fields](sectores/flag_rates/README.txt) |
| Circuit flag, Table 3 and Appendix E.1.4 | `marca-circuito` | `marca-circuito.json` gives all eight circuit rows. At `(r_c,r_m) = (2e-7,1e-5)` with VV first, fewer than 48 rounds per million are discarded and the silent-failure risk per retained round falls by a factor of at least 257. |
| Strip cutoff and resource comparisons, Tables 4, 9 and 10 | `rapido` | In `tablas_tiras/`, `circuit_rows.tex` and `phenomenological_rows.tex` give the table entries. `resources.json` gives the recorded times and memory. |
| Affine F17 instance, Table 5 and Appendix D | `f17` | Logical masks, implicit comparisons, 1,214 circuit traces in `certificate.jsonl`, and sufficient truncation orders 253,000 and 126,000. [Finite-sector inputs](certificados/sector_f17/README.txt) |
| 225-qubit pilot, Table 7 | `piloto` | `piloto/piloto_HGP_225_esquinas_20261124.json` contains three deterministic Bernoulli intervals and eight stratified circuit intervals. The latter have simultaneous coverage of at least `1 − 16 exp(−7)` under independent sampling within each stratum. |
| Twelve family-wide row inequalities, the Bernoulli density rule and gate order at matched rates (Table 8 and Figure 5(a)) | `revision-v2` | The first three logs give the twelve row ratios, the 225-qubit density bounds and the four pairs of matched schedules. The fourth checks Table 6 against the manuscript. |
| 225-qubit postselection, Table 6 | `postseleccion-v2` | Rebuilds `comparison.json` with the discard bounds for one Bernoulli and eight circuit points, and the circuit risk bounds. The next command computes the tighter Bernoulli risks used in the table. |
| Lee and the HGP adaptation of Meister at equal expected discard, Appendix E.1.6 | `meister-225` | Rebuilds `meister-225.json` from the affine seed at `p=1e-6, q=1e-4`, including both risk intervals and the seven-check refinement. Uses the same optional environment as the preceding command. |
| PATH25, Appendix E.3 | `path25` | Two decoder implementations agree on 4,096 syndromes. On the test circuit, enumeration covers 524,288 assignments and the component cutoff is complete at eight noisy locations. [Estimator outputs](path25/estimadores/README.txt) |
| Rare-event comparisons, Tables 11–12 and Appendix E.4 | `comparaciones` | Recomputes the estimates and mixing diagnostics from stored samples. The 61-qubit reference is deterministic. The 225-qubit circuit comparison includes both deterministic and stratified reference intervals. |
| Manuscript and Figures 1–7 | `pdf` | Compiles the PDF with three pdfLaTeX passes. The figures and bibliography are in the TeX source. |

**Software.** The reference runs used Windows x64 and Python 3.12.10. Use Windows for `path25`, the 61-qubit chain checks in `comparaciones`, both splitting runs and strip production. Four actions compile C++17 code. On Windows, `sectores`, `marca-tasas`, `dss-225` and `splitting-225` use Visual Studio x64 Build Tools on a CPU supporting AVX2 and BMI2. On x86-64, `sectores`, `marca-tasas` and `dss-225` also support GNU C++. The PDF requires pdfLaTeX, TikZ, PGFPlots and the other packages in the [preamble](manuscrito/Silent_logical_failure_in_small_set_flip_decoding.tex).

### Recompute the postselection comparisons in Table 6

With Python 3.12 and Git on your path, create a virtual environment
from the repository root.

~~~sh
python -m venv ../ssf-env
~~~

Activate it with `..\ssf-env\Scripts\Activate.ps1` in PowerShell or
`. ../ssf-env/bin/activate` in a POSIX shell. In PowerShell you can
also use `..\ssf-env\Scripts\python.exe` for `python` in the commands below.

~~~sh
python -m pip install -r comparaciones/postselection_v2/requirements.txt
python reproducir.py postseleccion-v2
python reproducir.py meister-225
~~~

The pinned Lee–English–Bartlett decoder scores the recorded input.
SSF supplies the correction and the terminal test in all nine rows.
Graph columns are sorted by the integer encoding the recorded bits
in positions 0–107 and the final data bits in positions 108–332.
BP uses the parallel product-sum schedule, one thread per worker and
at most 30 iterations. Order-zero LSD runs even when BP converges.
The threshold is calibrated from the input distribution before classifying
silent histories. Each row records the cutoff in millionths, discards
per million rounds, and the silent risk per retained round. Circuit
`risk_stratified` intervals reuse the sixteen bounds from the original
pilot, with joint coverage above 98.54%. `risk_deterministic` uses the
full deterministic remainder. Both coincide under Bernoulli noise.
The nine-point run took about seven minutes with four workers on the
reference Windows machine. The Bernoulli refinement took about two
minutes with one BP+LSD thread. Both commands can run directly after
installation and check their scientific outputs against the archived values.
The risk bounds in the JSON files are probabilities per retained round.
Table 6 displays those probabilities in units of `1e-18`.

- **Bernoulli, `p=1e-6, q=1e-4`.** The upper endpoint of
  `risk_per_retained_round.lee` in `meister-225.json` is
  `1.68889439986e-25`, which gives the abstract’s bound of `1.69e-25`.
  The same file gives the risk of the HGP adaptation of Meister’s score
  in `risk_per_retained_round.cluster_geometry`. Both rules discard
  2.58762–2.60134 rounds per million in expectation. [Calculation and reference output](comparaciones/postselection_v2/meister_225/README.md)
- **Circuit, VV first, `r_c=2e-7, r_m=1e-5`.** Select the matching
  schedule and rates in both files. In `comparison.json`, the upper
  endpoint of `cluster.risk_stratified` is below `5.70e-19`. In
  `marca-circuito.json` (from `python reproducir.py marca-circuito`), the SSF lower endpoint
  of `risk_per_retained_round` exceeds `1.17e-18`, and `gain_lower`
  exceeds 257 relative to keeping all rounds. Lee discards fewer than
  30.43 rounds per million and SSF discards more than 47.75. The risk
  comparison uses the joint stratified coverage above 98.54%.

To recompute just this circuit point with the cluster scorer, run
`python comparaciones/postselection_v2/compare_rules.py --point 6 --output ../vv-first-225.json`.

## 3. Evaluate another pair of Bernoulli rates

```sh
python reproducir.py regional --qubits 900 --p 4e-7 --q 5e-5
```

Choose `--qubits 225` or `900`. For positive rates, Corollary 6 requires `q <= 1e-4` and `z_p/z_q <= 9999/999999`, where `z_x = x/(1-x)`. You can enter a decimal, fraction or scientific notation. When `p = 0` and `0 <= q <= 1e-4`, the probability is zero.

For the example above, the interval is approximately `[6.7622e-15, 7.9308e-15]` (rounded outward). In `resultado.json`, `regional.lower` is the two-data, one-readout mass rounded down. The upper endpoint allows a remainder of up to 1.80% of that mass on 225 qubits or 17.28% on 900. `regional.estimate` gives the harmonic estimate and `regional.max_relative_percent` bounds its relative error (`null` at exact zero). At the rates of Table 2, the table uses additional sectors, so its endpoints can differ from the regional calculation.

## 4. Repeat the sampling and strip runs

`dss-225` repeats one million draws with replacement from the two-data, one-readout stratum, then enumerates all 2,721,600 inputs in that stratum. In `dss225.json` the reference sample has 119 silent failures and the census has 324. The sampling commands use the recorded seeds. To regenerate the splitting chains, run `splitting-225` or `splitting-61`. The recorded runs took 161.91 s on 225 qubits (28 rate levels) and 2.49 h on 61 qubits.

For strip production on Windows x64, extract `REPRODUCCION_SSF_QEC_20261125.zip` from the same Zenodo record into a separate directory (15 GB for download and extraction together). Run the strip commands in that directory. Its 1,412 catalog fragments and strip programs are the ones used for this revision. The runs allow 24 GiB of aggregate memory. A 32 GiB machine is recommended. To verify the recorded runs, use `python tiras/reproduce.py verify-recorded`. This command includes the catalog checks. To recompute the intervals, run the two commands below in order. The second reads the output directory of the first.

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

The 1,412 full strip catalog fragments are in `REPRODUCCION_SSF_QEC_20261125.zip` in release 4. The same file is also available in the [version 3.1 Zenodo archive](https://doi.org/10.5281/zenodo.24542801). The code and reference output for the postselection comparisons are in `comparaciones/postselection_v2/` in tag `circuit-postselection-v4`.
