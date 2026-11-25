T2-1 -- COMPLETE REPRODUCTION MATERIALS
16 November 2026

START HERE / LEER PRIMERO

This release adds the production dependencies and all 1,412 catalogue
shards to the compact strip-table release. It covers T2-1 (both circuit
runs, m = 6, 12, 24, 48), with its accepted T2-0/T2-1M/B4/SOLO_U inputs.
T2-2 is excluded. No experimental run was repeated to prepare this package.

The scientific sources, sealed plan, inputs, results and catalogue files
are byte-identical copies. The new reproduce.py only selects commands,
sets the working directory, checks integrity and prevents overwriting.
It does not change the decoder, ordering, noise laws, precision or budgets.

REQUIREMENTS

- Windows 64-bit; Python 3.12 (the production version), standard library.
- No Codex, Claude, network, original D: drive layout, pip or LaTeX needed.
- Extract the ZIP completely, to any short path (for example C:\T21).
  Use a ZIP64-capable extractor. Python's zipfile supports this archive.
- Catalogue inputs occupy 6,830,509,776 bytes. Allow about 7 GB extracted
  plus the ZIP. Future production outputs require additional disk space.
- Production retains up to 14 workers, 24 GiB aggregate memory and the
  original per-process guards. A 32 GiB machine or larger is recommended.
- Do not use python -O. Run from a terminal; quote paths containing spaces.

READ-ONLY CHECKS AND TABLE REPRODUCTION (NO EXPERIMENTS)

From the extracted release directory:

  python -I -B reproduce.py integrity

Checks every supplied file and all 1,412 catalogue hashes and byte counts.
This establishes file integrity, not a mathematical proof or a fresh census.

  python -I -B reproduce.py tables --out C:\T21_tables_new

Reconstructs the paper's tables and scalar intervals from the supplied exact
endpoints. Does not run the decoder, gas DP, catalogues or Monte Carlo.
It uses the unchanged compact postprocessor in tables/.

  python -I -B reproduce.py verify-recorded

Runs the original SOLO_U verifier and the original integral verifiers of
the recorded runs 1 and 2 (48, 96 and 128 rows respectively), without DP,
decoder or sampling. It reconstructs the finite interval arithmetic and
decisions; it does not independently recompute catalogue completeness or
the Monte Carlo binomial interval. Outputs go to the terminal. Historical
log files and recorded resources are evidence, not newly measured results.

FUTURE SCIENTIFIC REPRODUCTION (NOT RUN WHEN THIS RELEASE WAS PREPARED)

These are expensive commands. They never start through integrity, tables,
verify-recorded or --help. Each requires --execute-production. Use new
output folders OUTSIDE the release. On a failed control or resource limit,
the original program stops and retains evidence. There are no retries or
automatic budget changes. Do not interpret an incomplete output as success.

1. Run 1, using the included verified catalogues:

  python -I -B reproduce.py run-1 --out C:\T21_run1_new --execute-production

This reruns N1, the gas DP, interval evaluation and the fixed direct-sampling
control (7,000,000,000 histories; at most 200,000,000 decoder calls, including
the 279,049 shortcut checks). SOLO_U and the calibrations remain accepted
inputs, as in the published protocol. Run 1 has a 6.5-hour budget; the
Monte Carlo work budget is 1.5 hours within that total. Historical run 1
took about 3.32 hours on the author's machine, not a runtime guarantee.

2. Run 2, only after run 1 has passed its fixed continuation rule:

  python -I -B reproduce.py run-2 --out C:\T21_run2_new --previous C:\T21_run1_new --execute-production

Run 2 adds m = 48 and consolidates all four strips. Its original budget is
9 hours. The conductor rechecks run 1's result, seal and continuation rule.
Each run preserves the original 60 GB output cap, 384-bit arithmetic and
D = 10^16 for the two new k = 4 sum rules. The inherited D = 10^14 remains
unchanged. All precise controls are in the supplied PLAN_T2_1.txt.

3. Optional regeneration of the catalogues and their original engine:

  python -I -B reproduce.py rebuild-engine --out C:\T21_engine_new --execute-production

This invokes the unchanged T2-1M producer: catalogues plus its own DP and
checks. It is NOT a cheap hash check. Its original budget is 8 hours,
7.5 billion decoder calls, 24 GiB aggregate memory and 60 GB output.
It writes a new trabajo/ and compares its manifest with the supplied one.
It never replaces the frozen catalogues automatically. A resource stop or
different manifest requires inspection. This optional command is provided
but was not executed to build this release.

VERIFY NEW RUNS

The original verifier can also read future results. From workspace/:

  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_intervalos.py --verificar C:\T21_run1_new\RESULTADO_T2_1_C1.json
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_intervalos.py --verificar C:\T21_run2_new\RESULTADO_T2_1.json

The original entry points for SOLO_U, T2-0, B4 and the decoder dependencies
are preserved for inspection. They are not automatically rerun: their
accepted calibrations and reference results are explicit inputs to T2-1.
This release does not claim to replay every historical predecessor project.

CONTENTS AND INTERPRETATION

workspace/    Original relative layout: frozen code, dependencies, accepted
              calibration inputs, catalogues, both production results/logs,
              the sealed plan and the final interpretation addendum.
tables/       Unchanged compact strip package: table script, endpoints and
              phenomenological B4/B4R reference data. Its README describes
              that compact subpackage, not the full release around it.
PROVENANCE.json  Copy origins, source hashes, roles and catalogue inventory.
MANIFEST.json    SHA-256 and size of every release file except itself.
VALIDATION.json  Checks actually performed during packaging, with the
                 distinction between packaging checks and scientific reruns.

Read workspace/T2_1_CIERRE_DEC058_20261115/ADENDA_PRECISIONES_T2_1.txt for
the published interpretation: A2 is bounded growth, not constant error;
nine cells certify 10% with components but not the specified rival at the
same cutoff; equal-time comparisons are projections. Some cells do not
meet the accuracy target. Direct sampling is a compatibility check, with
nominal binomial coverage, not a proof of the certified upper endpoint.

The 128 circuit rows and phenomenological reference tables can be checked
without regenerating inputs. Full production commands regenerate circuit
T2-1 from the explicit accepted inputs; the optional engine command also
regenerates its catalogue inputs. Resource times and memory depend on the
machine. Result JSONs include timing/path provenance, so future full JSON
files need not be byte-identical even when scientific quantities agree.

Historical absolute paths occur in reports, metadata and an unused helper
in rect.py. The supported entry points use the bundled relative paths.
Do not run archived administrative scripts; none are included as entry
points. No protected manuscript, result or project-control file is changed.

RESUMEN EN ESPANOL

El paquete incluye los catalogos completos, las dependencias y las dos
corridas guardadas. "integrity", "tables" y "verify-recorded" no repiten
experimentos. Los comandos de produccion quedan preparados para otro
momento y requieren --execute-production. T2-2 queda fuera. La comprobacion
de portabilidad de esta entrega no equivale a una nueva corrida cientifica.
