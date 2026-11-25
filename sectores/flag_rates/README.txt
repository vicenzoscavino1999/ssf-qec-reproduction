Postselection and higher Bernoulli rates on the finite codes
==========================================================

The one-round rule flags a run when SSF makes a move and stops with a
nonzero observed syndrome. Table 3 reports discarded rounds per million, the silent-failure risk
conditional on keeping the round, and the gain. The gain divides the
unfiltered silent-failure probability by that conditional risk. The separate
comparison with input-syndrome postselection uses the two discard rates.

Bernoulli outputs
-----------------
The main guide (../../REPRODUCTION.md) gives the commands for evaluating
the stored tallies and for rebuilding the census. Both produce
marca-tasas.json, whose reference is pilot-results.json. To evaluate only
the Bernoulli flag and high-rate estimates, run this command from this
directory with a fresh output filename.

    python -B run_pilot.py ../.. --output ../../../ssf-bernoulli-flag.json

The following entries locate the reported results in that JSON.

    flag225_direct                   225-qubit flag at p=1e-6, q=1e-4
    orbit_bounds.900_at_reference    900-qubit flag at the same rates
    flag_regional.225                225-qubit gain over the odds region
    orbit_bounds.900_gain_on_cone    900-qubit gain over the odds region
    orbit_bounds.225_at_high_rates   225-qubit interval at p=1e-4, q=1e-3
    finite_region_225                225-qubit small-rate rectangle

Appendix E.1.3 derives the discard and conditional-risk bounds. Counts
through two faults set the discard interval. In the zero-readout sector
the weight-three count is exact on both codes. The weight-four count is
exact on 225 qubits and bounded on 900, where a weight-five bound also
enters. feasibility.py combines those counts and the binomial tails.

For the high-rate point in Appendix E.2.1, the inputs are the exact count
C40 = 14,904, the affine row bound for two data faults, the upper count
C31 <= 101,088, and orbit bounds for the (3,2) and (4,1) sectors.
orbit_tables.py sets up the orbits. bounded_orbits.cpp allows tied first
moves and feasibility.py charges genuine removals to smaller sectors.
Appendix E.1.1 derives these bounds. rectangle225() evaluates the
small-rate rectangle with the certified harmonic-estimate error.

Circuit outputs
---------------
marca-circuito.json gives the flag bounds for four rate pairs with VV first and
four with CC first on 225 qubits. Its reference is circuit-flag.json.
Each row contains discards_per_million, unflagged_joint,
risk_per_retained_round and gain_lower. circuit_flag.py computes the
bounds from the wire and check probabilities in Appendix E.1.4.
