# Reproduction guide

This guide accompanies *Certified estimates of silent logical failure in small-set-flip decoding with noisy readout* by Vicenzo Scavino Alfaro. It covers the manuscript in this checkout, the numerical certificates and the programs that reproduce the tables and figures. The GitHub checkout contains the code and the inputs for the quick check, the finite-code censuses and the flag and high-rate experiments. The complete Zenodo archive also contains all 1,412 binary strip-catalog fragments used by the production runs. [Source for this manuscript, tag `circuit-postselection-v2.1`](https://github.com/vicenzoscavino1999/ssf-qec-reproduction/tree/circuit-postselection-v2.1) · [Complete archive](https://doi.org/10.5281/zenodo.24521353)

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

The final line reports `OK` and the path to `resultado.json`. The quick check verifies the SHA-256 manifest, recomputes the Bernoulli intervals from the archived counts, and regenerates 128 strip-table rows in `tablas_tiras/`. It also checks the F17 and (3,4) truncation orders, the negative zero, the circuit-fiber example, the chain diagnostics and 62 numerical values in Figures 5–7. It also replays both decoder trajectories in Figure 1.

For Table 2 at `p = 10^-6` and `q = 10^-4`, the deterministic bounds on the silent-failure probability, rounded outward, are

| Qubits | Lower bound | Upper bound |
| ---: | ---: | ---: |
| 225 | 3.24009 × 10⁻¹⁴ | 3.27319 × 10⁻¹⁴ |
| 900 | 8.63039 × 10⁻¹⁴ | 9.69644 × 10⁻¹⁴ |

The log `01-calculos_del_articulo.log` retains more digits. To run the decoder and rebuild these two censuses from their seeds, use `python reproducir.py sectores`.

## Evaluate the regional interval

Use this checkout to evaluate Corollary 6 at your own independent Bernoulli rates. Select one of the two distance-six codes. Both examples use the decoder and tie order from Section 6.1.

```sh
python reproducir.py regional --qubits 225 --p 1e-6 --q 1e-4
python reproducir.py regional --qubits 900 --p 4e-7 --q 5e-5
```

For positive rates, require `q <= 1e-4` and `z_p/z_q <= 9999/999999`, where `z_x = x/(1-x)`. Supply rates as decimals, scientific notation or fractions. The command checks these inequalities exactly. For example, keeping `p = 1e-6` and lowering `q` to `5e-5` exceeds the odds-ratio bound. At `p = 0` and `0 <= q <= 1e-4`, the silent-failure probability is exactly zero.

The lower endpoint is the mass of the sector with exactly two data faults and one wrong readout. The upper endpoint allows a remainder of 1.80% of that mass on 225 qubits or 17.28% on 900. The program reads the archived counts, evaluates both endpoints with rational arithmetic and displays them rounded outward. It prints the result and records it under `regional` in `resultado.json`. The `estimate` field is the rounded harmonic estimate and `max_relative_percent` bounds the relative error of the displayed value. For an exact zero this field is `null`. The Table 2 intervals at the corner use additional sector counts.

## Reproduce the abort flag and higher Bernoulli rates

From this checkout, run

```sh
python reproducir.py marca-tasas
```

This command builds the 225- and 900-qubit models from their seeds. It compiles the C++17 decoder, counts all four-data inputs on 225 qubits, checks the low-weight abort rule with an independent Python decoder and recomputes the orbit bounds. The same independent decoder also checks all 7,920 reduced inputs in the 225-qubit row lemma with the published physical tie order and repeats the 292-representative census giving the three-data/one-readout ceiling. The command compares the complete recomputed result with the reference and saves it as `marca-tasas.json` in the output directory. `python reproducir.py rapido` checks the same rational intervals from the stored tallies.

At `p=1e-6, q=1e-4` the gain is the unfiltered silent-failure probability divided by the probability of silent failure conditional on keeping the round. Its lower bounds are 909.5 and 608.9 on the two codes. On 225 qubits at `p=1e-4, q=1e-3` the relative error bound is 5.396%. The latter calculation uses the affine row lemma, the exact four-data count of 14,904, the three-data/one-readout ceiling of 101,088 and upper counts on the two mixed sectors (3,2) and (4,1). [Proof outline, commands and output keys](sectores/flag_rates/README.txt)

## Circuit postselection on 225 qubits

Run the following command from the checkout to regenerate all eight circuit rows of Table 3.

```sh
python reproducir.py marca-circuito
```

The command uses Python's standard library. It builds both gate schedules, compares two SSF implementations on the observed syndromes and evaluates both implementations of the joint wire law. The discard census includes all one- and two-fiber responses and all 204,156 triples of check fibers. The 36 silent data triples at zero effective offset determine the leading retained risk. The program uses the deterministic lower endpoints archived for Table 6, which the separate `piloto` action rebuilds. It checks the output against the reference and writes `marca-circuito.json`. For each schedule and rate pair, the file reports discards per million attempted rounds, silent-failure risk per retained round and a lower bound on the reduction factor. All endpoints are rounded outward. [Circuit calculation](sectores/flag_rates/circuit_flag.py) · [Reference values](sectores/flag_rates/circuit-flag.json)

## Find the calculation behind a result

Use `python reproducir.py ACTION` with an action from the table. Enter the command and option names exactly as printed. By default, each run creates a fresh directory under `resultados_ssf/`, next to the checkout or extracted archive. To choose a location, pass `--salida` followed by a new directory outside the source tree. Each run writes its own logs and `resultado.json`. If a step fails, the program identifies the relevant log.

| Paper location | Action | What it reproduces |
| --- | --- | --- |
| Figure 1 | `rapido` | Replays the two data faults on sites 3 and 4, first with correct readouts and then with check 2 inverted. The log gives the ordered flips, the residual and the terminal recorded checks for each input. |
| Section 5, Table 1, Appendices B and C | `tasas` | Verifies the three (3,4) rate corners using rational arithmetic, the two circuit-order columns, the direct transfer through fibers and the Bernoulli and circuit certificates for (17,18). Includes the local first-contact checks. The [Bernoulli row derivation](certificados/fibras/proofs/H5A_3_4_esquina_2.txt) gives the twelve-ratio test in Appendix B.3. The run reports the maximum in `06-verificar_H5A_3_4_esquina_2.log`. |
| Section 5, negative zero and truncation orders | `rapido` | Counts the degree-four polynomial in two ways and brackets its negative zero in `(-0.097080, -0.097079)` using exact rational signs. The weight-six witness gives sufficient orders 6,482 and 6,575 for 225 and 900 qubits. The two-data, one-readout sector gives 3,028 and 3,057, respectively. |
| Section 5.1, Appendix D and Table 5 | `f17` | Constructs the H1 masks, repeats the H3 implicit comparison and generates and checks traces for 1,214 affine-sector configurations. Verifies the sufficient orders 253,000 and 126,000 and a circuit-sector contribution above `1.02 × 10^-499` at the rates in the paper. |
| Section 6.1, Table 2, Figure 5(b) and Appendix E.1 | `sectores` | Compiles the C++17 decoder, builds both seeds, performs the orbit-based censuses and compares the logical-class histograms and intervals with the archived direct censuses. The two-data, one-readout silent counts are 324 and 864. |
| Section 6.1, Table 6, Figure 5(a) and Appendix E.2 | `piloto` | Recomputes the three Bernoulli rows and eight circuit rows for 225 qubits. A second implementation checks the intervals, gate schedules and counts for the 16 samples. Their joint error probability is bounded by `16 exp(-7)`. |
| Sections 6.1–6.2, Table 3 and Appendices E.1–E.2 | `marca-tasas` | Rebuilds the flag, the 225-qubit affine row patterns and the finite high-rate interval from the seeds. |
| Section 6.2, Table 3 and Appendix E.1.4 | `marca-circuito` | Rebuilds the deterministic discard and retained-risk bounds on 225 qubits for all eight circuit points of Table 6. |
| Section 6.3, Tables 4, 7 and 8, Figures 6(a) and 7 | `rapido` | Regenerates the strip error bounds from the archived interval endpoints and retrieves the resource measurements. At 523 qubits, equal rates of `10^-5` and cutoff four, the rounded bounds are 1.221% for the component cutoff and 99.513% for the global cutoff. |
| Section 6.4, Tables 9 and 10, Appendix E.4 | `comparaciones` | Reanalyzes the archived runs. For 225 qubits, it reconstructs the Bennett estimate and its chain diagnostics. For 61 qubits, the original verifier checks the full run record. It also checks the exact frequency and the interval for the subset-sampling stratum. |
| Appendix E.3, PATH25 | `path25` | Compares two decoders on 4,096 syndromes and repeats the calculation over 524,288 assignments. The circuit estimator has zero remainder at cutoff eight. The verifier compares the deterministic outputs with the reference files. |
| Figures 1–7 and the manuscript | `pdf` | Runs pdfLaTeX three times in the output directory. TikZ draws the diagrams and the TeX source contains the plotted values and bibliography. |

The uniform (17,18) certificates apply to the entire family. The circuit point `r_c = 5 × 10^-33, r_m = 2 × 10^-11` and the sector contribution in Table 5 use the specific affine F17 instance. Both finite codes in Table 2 have distance six. The distance-ten product in the discussion is a comparison drawn from the literature.

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
- The C++17 actions are `sectores`, `marca-tasas`, `splitting-225` and `dss-225`. The Windows builds use Visual Studio Build Tools with the x64 C++ component and a processor supporting AVX2/BMI2. The package also provides GNU C++17 build commands on x86-64 for `sectores`, `marca-tasas` and `dss-225`. The `path25`, `splitting-61`, `splitting-225`, `tiras-1`, `tiras-2` and `catalogos` actions use the original Windows setup.
- The `pdf` action requires pdfLaTeX and the packages in the TeX preamble, including Latin Modern, TikZ, PGFPlots, microtype and flushend. A different TeX installation or build timestamp can change the rebuilt PDF's byte-level hash.
- The verifiers check integer counts and rational inequalities. Figure checks use the precision of the printed numerical labels. The hybrid intervals and Markov-chain experiments include their seeds, orders, counts and diagnostics. An `OK` for a splitting run means that the published estimate and its mixing diagnostics were reproduced. The 61-qubit reference interval is deterministic. The 225-qubit circuit reference interval uses the pilot's joint coverage bound. Actions that report runtime or memory measure them again on the current machine.

## Validation of the archived calculations

A fresh Windows extraction reproduced both orbit-based censuses, the 11 pilot intervals, the subset-sampling and splitting runs for 225 qubits, and the PATH25 checks. The archived 34-page PDF was rebuilt and its extracted text matched the archived reference page by page. Section and figure numbers in the table above refer to the revised manuscript in this checkout. For the multi-hour strip calculations and the 61-qubit splitting run, validation used the archived outputs. The production actions above rerun those calculations from their inputs.

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
