PATH25 syndrome census
======================

The path25 action in ../../REPRODUCTION.md runs the validation on a fresh copy
and compares its results with expected/. The literal and optimized decoders
agree on all 4,096 syndromes and all 8,192 quotient-class comparisons.

To run just the syndrome census on Windows x64, copy run/ into a fresh
working directory. Keep worker.py and inputs/ beside supervisor.py. From the
folder containing supervisor.py, run the following command with assertions
enabled.

    python -B -X utf8 supervisor.py

The supervisor writes attempt_01/RESULTS.json. The single worker has limits of
60 seconds, 256 MiB of memory and 32 MiB of output. A completed reference run
has status PASS_BOUNDED_PATH25.

Reference values
----------------
The census makes 12,438 oracle queries on 1,664 supports, with 376 decorations.
Of 565 recomposition sources, 123 have multiple components and 111 of those
contain decoder actions. The first three total-failure coefficients by weight
are [0, 4, 115], and the silent-failure coefficients are [0, 0, 14].

expected/ORACLE.json and expected/RECOMPOSITION.json contain the detailed
counts. For cross-machine comparisons use the JSON values. Raw binary table
hashes depend on native byte order, and timings depend on the machine.
