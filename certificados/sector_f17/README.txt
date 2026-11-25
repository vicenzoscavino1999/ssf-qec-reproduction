Finite circuit sector on the affine F17 instance
===============================================

Appendix D.3 uses a finite classification at r_c = 5e-33 and r_m = 2e-11.
The f17 action in ../../REPRODUCTION.md rebuilds its traces and evaluates
the finite-sector lower bound at those rates.

To inspect the supplied traces, run the following command from this directory
with Python 3.10 or later and assertions enabled.

    python finite_sector.py verify certificate.jsonl

The 1,214 retained backgrounds contain 28,829 states including the stopping
states. Of these backgrounds, 810 are silent and 404 are detectable. Each
silent history leaves the same 34-line residual. Table 5 breaks down the
536 backgrounds outside the empty case and the first two retained families.

To generate a new trace file, run the producer and the verifier in order.

    python finite_sector.py produce reproduced.jsonl
    python finite_sector.py verify reproduced.jsonl

Choose a filename that does not already exist. The producer uses bit masks.
The verifier builds the point–line incidence sets independently and checks
each maximum, tie, stopping decision and residual. Both operations have a
60-second limit.

Earlier sampling-budget calculation
-----------------------------------
The optional command below evaluates the sampling budget stored in
rational_inputs.json at that file's rates (r_c about 9.66e-82 and r_m about
1.47e-20). Its output is a sufficient number of attempts for 10% relative
error at 95% confidence under the specified proposal.

    python verify_budget.py
