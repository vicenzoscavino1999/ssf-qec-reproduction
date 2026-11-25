STRIP COMPONENT CERTIFICATES -- 16 November 2026

Companion data and code snapshot for the growing-strip sections of
Exact canonical recomposition and silent-failure estimation for
small-set-flip decoding.

PORTABLE POSTPROCESSING (Python 3.11 or newer; standard library only)
  python reproduce_tables.py --output /absolute/path/to/new-output
Use a new output directory. This checks the package hashes, reconstructs
the 128 circuit intervals from exact B and fixed-point U/E/O endpoints,
checks rational phenomenological interval identities, and regenerates
the two displayed error tables and the resource data. It also checks the
stored direct-sampling totals and compatibility. It does NOT rerun the
decoder, the catalogue producers, the gas DP, or Monte Carlo. It does not
independently recompute the binomial interval or catalogue completeness.

DATA
  results/circuit.json: final combined 32-cell, 4-cutoff circuit output.
  results/phenomenological.json: corrected phenomenological output.
  results/completeness_phenomenological.json: auxiliary-weight sum rules.
  results/calibration.json: circuit fibres and small-cutoff calibration.
  results/engine.json: catalogue and engine output.
  results/direct_sampling.json: the T2-1 direct-sampling control, with
    integer counts and the nominal 99.9% binomial interval.
  results/circuit_resources.json: measured production resources.
  results/catalogue_manifest.json: identities of the large catalogue shards.
Exact fractions and fixed-point endpoints are authoritative; display
decimals are rounded summaries. Resource measurements are historical
observations and are not reproduced by the postprocessing command.

PRODUCTION SNAPSHOT
  source_snapshot/ preserves selected original producer programs and
    their original paths, hashes and import conventions.
  methods/ contains the frozen protocol and derivations.
This is a compact evidence release, not a self-contained full-production
launcher. The production code uses additional original workspace modules
and about 6.8 GB of catalogue shards, which are NOT included. The original
catalogues and their dependencies remain in the research workspace. A
fresh production reproduction must restore those dependencies and follow
the original resource limits. Do not run these scripts as a side effect
of table postprocessing. No new scientific computation was performed
while preparing this release.

SCOPE
Four fixed-width strips with one logical qubit, beta=1/2, cutoffs 1..4;
phenomenological noise and one specified noisy extraction round. This
release supports finite cutoff comparisons, not a global runtime
speedup or new F17 noise or asymptotic guarantees. The independent T2-2
transport/splitting pilot and subsequent H5 exploration are outside this
release and are not used in the manuscript's new claims.
