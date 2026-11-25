PATH25 CONNECTED-HOST VALIDATION
Archived results: 28 October 2026. Package: 6 November 2026.

This package copies the existing scientific code, frozen inputs and outputs
byte for byte. No decoder or oracle was executed during packaging.
Requirements: Windows and Python 3.10+; standard library only.

From a fresh extraction run:
  python -B -X utf8 run/supervisor.py
Never use python -O (assertions are scientific checks).
The supervisor creates run/attempt_01 and refuses an existing directory.
It uses a hidden single worker with 60 s, 256 MiB and 32 MiB stop ceilings.
A slow computer can reach the stop limit; partial results are not success.
The immutable reference outputs are in expected/, separate from new output.
RESULTS.json records 4096 matched syndromes, 12438 oracle queries, 1664
supports and 565 recomposition sources, including 123 multi-component
cases and 111 with actions. ORACLE.json and RECOMPOSITION.json give scope.
Timing, process memory, paths and native-endian binary hashes are provenance
fields; compare integer/rational scientific results directly. The archived
revision records three stationary cancellation pairs (empty traces).

The catalogue covers the supports specified in the manuscript; it is not
a complete instance catalogue or a global probability estimator. The code
imports pinned P1 functions, without invoking that worker's historical main.
Original scientific files retain historical names and metadata.
MANIFEST_SHA256.json covers every supplied file except itself.
