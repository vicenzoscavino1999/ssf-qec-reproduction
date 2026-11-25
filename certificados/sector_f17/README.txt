F17: a finite weighted circuit-failure sector
===========================================

This supplement accompanies Proposition "A finite circuit contribution"
and Appendices "A finite weighted circuit sector" / "A sufficient sampling
budget for the primitive circuit" of the manuscript dated 7 November 2026.

Requirements: Python 3.10 or later, standard library. No operating-system
specific modules, circuit simulator, external decoder, or absolute paths.

From this extracted directory:

    python finite_sector.py verify certificate.jsonl
    python verify_budget.py

To reproduce the finite certificate in a NEW file:

    python finite_sector.py produce reproduced.jsonl
    python finite_sector.py verify reproduced.jsonl

The producer refuses to overwrite an existing file. Each finite operation
has a 60-second bound, checked at each background. Slower hardware can
raise a timeout; this does not change the mathematical domain. The recorded
certificate contains 1214 backgrounds and 28829 states, including stops.
The producer uses bit-mask intersections; the verifier independently builds
point-to-line incidences and accumulates overlaps. Every maximum, tie, chosen
line, stop, correction XOR and physical residual is verified. The domains
are specified by two separate predicates and compared in their full order.
The outcome is 810 silent nonstabilizer residuals and 404 detectable ones.
All silent residuals are the 17 horizontal plus 17 vertical lines.

The manuscript proves the reduction from the full FGL candidate family to
the 306 line comparisons used here. It also proves the primitive-noise
mixture, the nonstabilizer criterion and the factor 306 across VV rows.
Those mathematical proofs are separate from these finite program checks.

The retained domain includes the empty background, the two stated families
of weights at most five, and the specified remaining singletons, pairs and
triples. No labels are inferred for the 3730 omitted backgrounds. Their
exact weight polynomial yields the local interval in the manuscript.
Other triples V, larger backgrounds and the full circuit failure probability
are outside this finite certificate.

rational_inputs.json contains the exact operation/readout rates at the
manuscript's rational corner, reference interval endpoints and the sufficient
attempt count. verify_budget.py reconstructs theta and the mixture, uses
90-digit outward rational brackets and Bernoulli lower bounds for the two
large powers, then recomputes the global proposal normalizer and the ceiling
2000*Zbar/L_lower. The displayed decimal values are approximations. The
integer ceiling is sufficient for 10% relative error at 95% confidence for
the specified proposal. No random samples or global decoder runs are used
to obtain that ceiling; it is neither a measured runtime nor a necessary
sample complexity.

The files are independently portable. SHA256SUMS.txt identifies their exact
contents and is an integrity check, separate from the mathematical proofs.
