# Reproduction guide

This guide accompanies *Estimating silent logical failure in small-set-flip decoding* by Vicenzo Scavino Alfaro. It covers the 34-page manuscript, the numerical certificates and the programs that reproduce the tables and figures. The GitHub checkout contains the code and the inputs for the quick check and the two finite-code censuses. The complete Zenodo archive also contains all 1,412 binary strip-catalog fragments used by the production runs. [Source repository](https://github.com/vicenzoscavino1999/ssf-qec-reproduction) · [Complete archive](https://doi.org/10.5281/zenodo.24521353)

## Run the quick check

Use 64-bit Python 3.12 with assertions enabled. The quick check uses the Python standard library. From the directory containing `reproducir.py`, run

```sh
python reproducir.py rapido
```

If you downloaded the ZIP, allow 15 GB for the archive and its extracted contents. These PowerShell commands extract it into an empty directory with a short path and run the same check. On Windows, `py -3.12` selects Python 3.12 explicitly.

```powershell
py -3.12 -m zipfile -e "REPRODUCCION_SSF_QEC_20261124.zip" "D:\SSF"
Set-Location -LiteralPath "D:\SSF"
py -3.12 reproducir.py rapido
```

A ZIP64-compatible archive utility can also extract the package. Leave the smaller ZIP libraries inside the modules in place. Python reads them directly.

The final line reports `OK` and the path to `resultado.json`. The quick check verifies the SHA-256 manifest, recomputes the Bernoulli intervals from the archived counts, and regenerates 128 strip-table rows in `tablas_tiras/`. It also checks the F17 and (3,4) truncation orders, the negative zero, the circuit-fiber example, the chain diagnostics and 62 numerical values in Figures 4–6.

For Table 2 at `p = 10^-6` and `q = 10^-4`, the deterministic bounds on the silent-failure probability, rounded outward, are

| Qubits | Lower bound | Upper bound |
| ---: | ---: | ---: |
| 225 | 3.24009 × 10⁻¹⁴ | 3.27319 × 10⁻¹⁴ |
| 900 | 8.63039 × 10⁻¹⁴ | 9.69644 × 10⁻¹⁴ |

The log `01-calculos_del_articulo.log` retains more digits. To run the decoder and rebuild these two censuses from their seeds, use `python reproducir.py sectores`.

## Find the calculation behind a result

Use `python reproducir.py ACTION` with an action from the table. Enter the command and option names exactly as printed. By default, each run creates a fresh directory under `resultados_ssf/`, next to the checkout or extracted archive. To choose a location, pass `--salida` followed by a new directory outside the source tree. Each run writes its own logs and `resultado.json`. If a step fails, the program identifies the relevant log.

| Paper location | Action | What it reproduces |
| --- | --- | --- |
| Section 5, Table 1, Appendices B and C | `tasas` | Verifies the three (3,4) rate corners using rational arithmetic, the two circuit-order columns, the direct transfer through fibers and the Bernoulli and circuit certificates for (17,18). Includes the local first-contact checks. |
| Section 5, negative zero and truncation orders | `rapido` | Counts the degree-four polynomial in two ways and brackets its negative zero in `(-0.097080, -0.097079)` using exact rational signs. The weight-six witness gives sufficient orders 6,482 and 6,575 for 225 and 900 qubits. The two-data, one-readout sector gives 3,028 and 3,057, respectively. |
| Section 6, Appendix D and Table 4 | `f17` | Constructs the H1 masks, repeats the H3 implicit comparison and generates and checks traces for 1,214 affine-sector configurations. Verifies the sufficient orders 253,000 and 126,000 and a circuit-sector contribution above `1.02 × 10^-499` at the rates in the paper. |
| Section 7.1, Table 2, Figure 4(b) and Appendix E.1 | `sectores` | Compiles the C++17 decoder, builds both seeds, performs the orbit-based censuses and compares the logical-class histograms and intervals with the archived direct censuses. The two-data, one-readout silent counts are 324 and 864. |
| Section 7.1, Table 5, Figure 4(a) and Appendix E.2 | `piloto` | Recomputes the three Bernoulli rows and eight circuit rows for 225 qubits. A second implementation checks the intervals, gate schedules and counts for the 16 samples. Their joint error probability is bounded by `16 exp(-7)`. |
| Section 7.2, Tables 3, 6 and 7, Figures 5(a) and 6 | `rapido` | Regenerates the strip error bounds from the archived interval endpoints and retrieves the resource measurements. At 523 qubits, equal rates of `10^-5` and cutoff four, the rounded bounds are 1.221% for the component cutoff and 99.513% for the global cutoff. |
| Section 7, Tables 8 and 9, Appendix E.4 | `comparaciones` | Reanalyzes the archived runs. For 225 qubits, it reconstructs the Bennett estimate and its chain diagnostics. For 61 qubits, the original verifier checks the full run record. It also checks the exact frequency and the interval for the subset-sampling stratum. |
| Appendix E.3, PATH25 | `path25` | Compares two decoders on 4,096 syndromes and repeats the calculation over 524,288 assignments. The circuit estimator has zero remainder at cutoff eight. The verifier compares the deterministic outputs with the reference files. |
| Figures 1–6 and the manuscript | `pdf` | Runs pdfLaTeX three times in the output directory. TikZ draws the diagrams and the TeX source contains the plotted values and bibliography. |

The uniform (17,18) certificates apply to the entire family. The circuit point `r_c = 5 × 10^-33, r_m = 2 × 10^-11` and the sector contribution in Table 4 use the specific affine F17 instance. Both finite codes in Table 2 have distance six. The distance-ten product in the discussion is a comparison drawn from the literature.

## Repeat the sampling or rebuild the strip catalogs

| Action | Calculation and resources |
| --- | --- |
| `splitting-225` | Compiles the sampler, builds the VV-first circuit model, repeats the anchor estimate and the 28 levels down to `r_c = 2 × 10^-7, r_m = 10^-5`, reanalyzes the chains and times the reference estimator. The published splitting run took 161.91 s. |
| `dss-225` | Rebuilds the model tables and compiles the sampler. It makes one million draws with the recorded seed and counts all 2,721,600 configurations in the (2,1) stratum. The published counts are 119 sample hits and 324 silent configurations in the census. |
| `splitting-61` | Runs the original seeds, anchor and rate levels for the 61-qubit strip. The published run took 2.49 h. To analyze its archived records, use `comparaciones`. |
| `tiras-1` and `tiras-2` | Recompute the strip production results from the complete catalogs. The original first run took 3.32 h under a 6.5 h budget. The second run includes 523 qubits and has a 9 h budget. |
| `catalogos` | Regenerates the 1,412 strip-catalog fragments and checks their hashes against the archived catalog. The original budget was 8 h, 24 GiB of aggregate memory and up to 60 GB of output. |

Run the two strip calculations from the extracted Zenodo archive. Choose new output directories.

```powershell
py -3.12 reproducir.py tiras-1 --salida "D:\ssf_results\strip-run-1"
py -3.12 reproducir.py tiras-2 --previa "D:\ssf_results\strip-run-1" --salida "D:\ssf_results\strip-run-2"
```

The second command reads the first run's result and completion record from the directory supplied to `--previa`. The production programs retain partial outputs when a resource limit is reached. Allow at least 32 GiB of RAM and additional disk space for new outputs. The published timings refer to the machine used for the paper. The `tiras/` directory contains the individual production specifications and resource budgets.

To check the archived strip production separately, run `python tiras/reproduce.py verify-recorded` from the extracted archive on Windows. It verifies the catalog hashes and the 128 result rows.

## Software and interpretation

- The full reproduction was tested on Windows x64 with Python 3.12.10. Python assertions must be enabled for the certificate checks. The quick check and rational-arithmetic calculations use the standard library.
- The C++17 actions are `sectores`, `splitting-225` and `dss-225`. The Windows builds use Visual Studio Build Tools with the x64 C++ component and a processor supporting AVX2/BMI2. The package also provides GNU C++17 build commands on x86-64 for `sectores` and `dss-225`. The `path25`, `splitting-61`, `splitting-225`, `tiras-1`, `tiras-2` and `catalogos` actions use the original Windows setup.
- The `pdf` action requires pdfLaTeX and the packages in the TeX preamble, including Latin Modern, TikZ, PGFPlots, microtype and flushend. A different TeX installation or build timestamp can change the rebuilt PDF's byte-level hash.
- The verifiers check integer counts and rational inequalities. Figure checks use the precision of the printed numerical labels. The hybrid intervals and Markov-chain experiments include their seeds, orders, counts and diagnostics. An `OK` for a splitting run means that the published estimate and its mixing diagnostics were reproduced. The 61-qubit reference interval is deterministic. The 225-qubit circuit reference interval uses the pilot's joint coverage bound. Actions that report runtime or memory measure them again on the current machine.

## Validation of the archived calculations

A fresh Windows extraction reproduced both orbit-based censuses, the 11 pilot intervals, the subset-sampling and splitting runs for 225 qubits, and the PATH25 checks. The rebuilt PDF had 34 pages and its extracted text matched the supplied PDF page by page. For the multi-hour strip calculations and the 61-qubit splitting run, validation used the archived outputs. The production actions above rerun those calculations from their inputs.

## Locate the supplements cited in Appendix F

| Supplement cited in the paper | Location in the complete package |
| --- | --- |
| `SOURCE_CUMULANT_CERTIFICATE.zip` | `certificados/f17/` |
| `FIBRE_CERTIFICATES.zip` | `certificados/fibras/` |
| `REPRODUCCION_PAPER_20261117_v3.zip` | The two cited calculations are in `certificados/supplement/h1/` and `certificados/supplement/h3/`. |
| `F17_CIRCUIT_SECTOR_20261107.zip` | `certificados/sector_f17/` contains the geometry and traces. `calculos_del_articulo.py` evaluates the rates in the paper. |
| `SSF_SECTORES_225_900_20261124.zip` | `sectores/` |
| `HGP225_PILOT.zip` | `comparaciones/` contains the generator and independent verifier for the 11 intervals. |
| `QEC_RARE_EVENT_COMPARISONS_20261124.zip` | `comparaciones/`, with the 61-qubit programs and per-chain records under `tiras/workspace/T2_2_*/`. |
| `PATH25_VALIDATION_20261106.zip` | `path25/validacion/` |
| `PATH25_CERTIFIED_ESTIMATORS_20261112.zip` | `path25/estimadores/` |
| `STRIP_COMPONENT_CERTIFICATES_20261116.zip` | `tiras/tables/` |
| `T2_1_REPRODUCCION_COMPLETA_20261116.zip` | `tiras/`, including the inputs for both production runs. |

The archive retains the original filenames and internal paths. Its root manifest, `SHA256.json`, covers every distributed file except the manifest itself.
