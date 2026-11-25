PATH25 PILOT REPRODUCTION

From the root of the reproduction package, run

    python reproducir.py path25

Use 64-bit Python 3.12 on Windows. The driver copies the validation inputs to a fresh output directory, runs the supervisor and compares the deterministic scientific results with expected/. It also reproduces the finite-estimator calculation described in the main guide.

To run the validation supervisor on its own, copy run/ to a fresh directory and run the following command from that copy.

    python -B -X utf8 supervisor.py

Keep worker.py and the complete inputs/ directory beside supervisor.py. The supervisor creates attempt_01/ and requires that directory to be new. It imports the frozen P1 functions and runs a single worker. The original 60 s limit can stop a run on a slower machine. If the time limit is reached, retain the partial outputs for inspection.

Expected exact results

    4,096 syndromes and 8,192 quotient-class comparisons
    12,438 local oracle queries
    376 decorations on 1,664 queried supports
    565 recomposition sources, including 123 multicomponent cases and 111 with actions
    103 candidates and 43 distinct footprints
    eta = 1 and refined density 1/3 after fixed-order dominance
    First total-failure coefficients by weight 0, 1, 2: [0, 4, 115]
    First silent-failure coefficients at the same weights: [0, 0, 14]

The original reference tables, spectra and traces are in expected/. The files INFORME_PATH25.txt and LEMA_DENSIDAD_ORDEN.txt contain the historical analysis and the proof of the uniform fixed-order bound. REVISION_RESULTADOS.json records the subsequent mathematical review. That review made no new decoder calls. Binary uint64 table hashes depend on native byte order. The integer spectra and JSON tables support comparisons across machines.
