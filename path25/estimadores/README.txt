Exact component estimates on PATH25
==================================

The path25 action in ../../REPRODUCTION.md runs the syndrome census and
the component estimators. This directory contains the estimators for the
nonregular code with 25 data qubits, 12 checks and one logical qubit.
To run only the estimators, use 64-bit Python 3.12 from this directory
and choose a new output directory.

    python -I -B reproduce_h6.py --out ../../../path25-results

Circuit noise
-------------
Eight of the 217 locations in the extraction round are noisy. The two
rates are 1/100 and 1/1000. The exact reference enumerates 524,288 fault
assignments. At rate 1/100 the first cutoff meeting the 10% target has
certified error at most 9.66%. At k = 8 the remainder is zero.

Independent data and readout noise
----------------------------------
The three rate pairs (p,q) are (1/50,1/100), (1/100,1/200) and (1/100,0).
The first cutoff meeting the 10% relative-error target is k = 4 for
(1/50,1/100) and k = 3 for the other two pairs. At k = 5 the largest
certified relative error across the three pairs is 0.282%, rounded upward.

REPRODUCTION_SUMMARY.json lists the six deterministic JSON outputs and
whether each matches its reference byte for byte. The outputs contain the
probabilities and interval endpoints. The separate timings cover fixed
batches of rates. Direct-sampling workload estimates use the normal
approximation 384(1-P)/P at the known reference probability.

The archived U2 report gives the implementation credits for the circuit verifier.
https://github.com/vicenzoscavino1999/ssf-qec-reproduction/blob/7d53e0764ea732743df65ee940b4023992d87d1f/path25/estimadores/qec_revision_matematica/05_COMPROBACIONES/H6_U2_CIRCUITO_PATH25_MULTICOMPONENTE_20261112/INFORME_U2.txt
