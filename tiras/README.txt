Strip calculations on 61, 127, 259 and 523 qubits
=================================================

Section 4 of ../REPRODUCTION.md gives the commands, run order and resource
budgets for fresh calculations. The compact source copy is sufficient for
rebuilding the published tables from their endpoints.

Inputs and results
------------------
tables/results/ contains the exact endpoints and the recorded time and memory
measurements. tables/README.txt describes the table files. workspace/ contains
the decoder and interval programs, the two recorded runs and their input
inventories. In the full ZIP it also contains all 1,412 catalog fragments.

The recorded-run verifier checks each catalog hash and then the interval
arithmetic for the 48-row calibration and the 96- and 128-row run outputs.
For a file-integrity check alone, run this command from this directory.

    python -I -B reproduce.py integrity

Production uses 384-bit arithmetic and up to 14 workers. The numerical
interpretation of the strip comparisons is in
workspace/T2_1_CIERRE_DEC058_20261115/ADENDA_PRECISIONES_T2_1.txt.
