Strip table inputs and outputs
=============================

The quick check in ../../REPRODUCTION.md builds the tables from the stored
endpoints. To run just the table calculation, use Python 3.11 or later from
this directory and choose a new output directory.

    python reproduce_tables.py --output ../../../strip-tables

circuit_rows.tex and phenomenological_rows.tex contain the table entries.
resources.json contains the recorded timing and memory data.
POSTPROCESSING.json reports how many circuit rows were checked (128).
The endpoints are exact fractions or fixed-point bounds. The tables round
the displayed error bounds upward.

Reference inputs under results/
-------------------------------
circuit.json gives the 32 circuit cells at four cutoffs. phenomenological.json
and completeness_phenomenological.json give the independent-noise endpoints
and their auxiliary weight sums. calibration.json and engine.json contain the
small-cutoff checks and catalog counts. circuit_resources.json records the
production costs. direct_sampling.json contains the sampling counts and their
nominal 99.9% binomial interval. catalogue_manifest.json identifies the binary
catalog fragments.

methods/ gives the derivations and source_snapshot/ contains the corresponding
producer code. To repeat the original strip runs, use the full distribution
and section 4 of the main guide.
