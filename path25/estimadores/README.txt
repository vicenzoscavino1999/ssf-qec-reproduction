PATH25_CERTIFIED_ESTIMATORS_20261112
Certified finite estimators on connected PATH25

Supplementary material for the section "Certified silent-failure estimates
on connected PATH25" and the appendix "Certified finite estimators on
PATH25" of "Exact canonical recomposition and silent-failure estimation for
small-set-flip decoding" (draft of 12 November 2026).

1. CONTENTS
reproduce_h6.py          launcher (see section 3)
MANIFEST.sha256.txt      SHA-256 of every other file in this archive
qec/SSF_QEM_extension_v2.zip
                         unchanged copy of the encoded-observable package;
                         read by contrastar_h6_u4.py (pinned by hash)
qec_revision_matematica/05_COMPROBACIONES/
  PILOTO_CONECTADO_PATH25_20261028/   three files of the earlier PATH25
                         calculations (PREFLIGHT.json, REFERENCE_TABLE.json,
                         REFERENCE.json), pinned by hash in the programs
  H6_U2_CIRCUITO_PATH25_MULTICOMPONENTE_20261112/
                         one-round circuit on PATH25 with eight noisy
                         locations: verificar_h6_u2.py and RESULTADO_H6_U2.json
                         (all 524288 assignments, oracle, gas identities)
  H6_U3_ESTIMADOR_CERTIFICADO_CIRCUITO_PATH25_20261112/
                         circuit estimator (estimador_h6_u3.py), its sealed
                         estimate (ESTIMACION_H6_U3.json, SELLO_ESTIMACION.json)
                         and the later comparison (contrastar_h6_u3.py,
                         RESULTADO_H6_U3.json)
  H6_U4_ESTIMADOR_CERTIFICADO_FENOMENOLOGICO_PATH25_20261112/
                         the same for mixed phenomenological noise
  H6_U5_RECURSOS_G4_Y_MATRIZ_FINAL_20261112/
                         resource measurements (medir_h6_u5.py,
                         MEDICIONES_H6_U5.json, mediciones_registradas/) and
                         the weight-four cancellation control
                         (control_v14_fenomenologico.py,
                         CONTROL_V14_FENOMENOLOGICO.json)
Each unit folder also holds PLAN_UNIDAD.txt, fixed before its programs were
written and run, and INFORME_*.txt, the working report. Both are internal
documents in Spanish; section 5 lists how some of their statements must be
read. Exact fractions are in the JSON files.

2. REQUIREMENTS
Python 3 with its standard library only (the runs were made with Python
3.12 on Windows). Run every program with "python -I -B". No program writes
over an existing file.

3. ONE-COMMAND REPRODUCTION
    python -I -B reproduce_h6.py --out <new directory outside this archive>
The launcher checks the manifest, then runs the stored programs without
modifying them and compares byte by byte:
  RESULTADO_H6_U2.json, ESTIMACION_H6_U3.json, RESULTADO_H6_U3.json,
  ESTIMACION_H6_U4.json, RESULTADO_H6_U4.json, CONTROL_V14_FENOMENOLOGICO.json.
The comparison programs (contrastar_*) are run on the stored sealed
estimates, as in the original order: the estimate was written and hashed
first, and the reference was computed afterwards. Seal files, resource
records and times are not deterministic and are not compared. The launcher
writes REPRODUCTION_SUMMARY.json and exits with status 0 only if every
comparison holds. Expected run time: about half a minute.

4. INDIVIDUAL COMMANDS (from the folder of each program)
    python -I -B verificar_h6_u2.py --output <new.json> [--tiempos <new.txt>]
    python -I -B estimador_h6_u3.py --output <new.json> --sello <new.json> --recursos <new.json> [--kmax K]
    python -I -B contrastar_h6_u3.py --estimacion <json> --sello <json> --recursos-estimador <json> --output <new.json> --recursos <new.json>
    python -I -B estimador_h6_u4.py --output <new.json> --sello <new.json> --recursos <new.json> [--kmax K]
    python -I -B contrastar_h6_u4.py --estimacion <json> --sello <json> --recursos-estimador <json> --output <new.json> --recursos <new.json>
    python -I -B control_v14_fenomenologico.py --output <new.json>
    python -I -B medir_h6_u5.py --output <new.json>
medir_h6_u5.py repeats 28 runs of the stored programs and writes them into
a new folder "mediciones" next to itself; the folder of the recorded runs
is shipped as "mediciones_registradas" so that this command can run.
Times and memory are machine dependent.

5. HOW TO READ THE INTERNAL REPORTS AND RESOURCE RECORDS
These precisions follow an independent review and take precedence over the
wording of INFORME_*.txt; the stored programs and outputs are unchanged.
a) Monte Carlo. The sample counts quoted for direct Monte Carlo use the
   normal approximation 384(1-P)/P with the exact P known for this
   benchmark. They are an orientation only, not a finite-confidence
   guarantee, and no sampling was run.
b) Rounding. Certified errors must be rounded upward: 9.66 % (exactly
   9.6522...%) rather than "<= 9.65 %" for the circuit at e=1/100, and
   0.282 % (exactly 0.28111...%) rather than "<= 0.28 %" for the
   phenomenological estimator at k=5. The 10 % target is met in all
   reported cases.
c) Resources. The estimator reaches the 10 % target in less time in some
   finite cases; at the principal phenomenological pair the exact reference
   is cheaper. These tests establish no general or scaling advantage. The
   circuit row k=8 is the same estimator with complete cutoff, not a
   separate enumeration program. The rows k=8 (circuit) and k=5
   (phenomenological) are single runs; other times are medians of five
   whole-process runs. Reference times are those of the comparison
   processes, including their controls and input/output; internal times
   are in the resource records. Each process evaluates a fixed batch of
   rates, not one isolated query.
d) Contract, stopping and reproduction. The circuit uses e=1/100 and
   e=1/1000 on eight selected locations of one extraction round; the
   phrase "tasas de lectura del 1-2 %" in the internal reports mixes the
   two contracts and should be disregarded. The phenomenological pairs are
   (p,q) = (1/50,1/100), (1/100,1/200) and (1/100,0). The estimators sweep
   the cutoff up to kmax and then select the first k meeting the 10 %
   criterion; the selection uses no reference value but is not an
   automatic stop. In the 28 runs of medir_h6_u5.py, 13 complete outputs
   are byte-identical to the stored files and 15 runs with a partial
   cutoff reproduce exactly the corresponding rows.

6. SCOPE
One nonregular instance (PATH25: 25 data qubits, 12 checks, one logical
qubit). The circuit has eight noisy locations out of 217 and one round.
The certificate is a finite component certificate, not the truncated
partition algorithm of the main theorems. Nothing here covers F17, a
growing regular family, noise at every circuit location or several rounds.
