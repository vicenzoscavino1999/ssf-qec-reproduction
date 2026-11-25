Uniform source-cumulant certificate for small-set-flip decoding

Extract the entire archive and run from the extracted directory:

    python verify.py

Requires Python 3.10 or later and its standard library. Use the default
interpreter mode (assertions enabled). Keep the nested ZIP intact in
certificate/. Python imports supporting modules directly from it.

For the joint rectangle 0 <= p <= 5e-33 and 0 <= q <= 2e-11, the command
checks the deterministic profile at degree (17,18), the local charges and
group recurrences, the infinite tail, fifteen rooted-tree inequalities,
and the exterior loads. The three source rows are printed in VV, CC and
measurement order. The last two numbers printed are bounds on the
higher-order difference from a reference in the working notes. Those
two numbers are distinct from the full physical loads in the paper.
Numerical assertions use integers or rational arithmetic. Output decimals
give the sizes of the bounds. No decoder data or simulation is required.

Proof guide (in proofs/). The general result is in sections 1 and 2 of
H5A_fusiones_sin_repetir.txt. The later sections discuss a separate
exploratory rate. For the counting proof use H5A_p_general.txt for the
chronological dual and support size, H5A_interacciones_por_fuentes.txt
and H5A_primera_fusion.txt for removing isolated pieces, and
H5A_consumo_y_alcance.txt for the first-crossing inequalities. The
causal encoding and the fifteen inequalities are in sections 1 to 5 of
H5A_contactos_entre_ejecuciones.txt. The marked series and its tail are
in sections 1 to 3 of H5A_primer_apagado_por_grupos.txt. The normalized
partition and the effective approximation procedure are in the first
two sections of H5A_cumulantes_y_componentes.txt and in the paper.
The general result relies on the uniform decoder and physical-filter
contract in the manuscript. Times for a probability estimator are not
part of the certificate.

expected.txt gives the reference console output. SHA256SUMS identifies
all included files other than the checksum list itself.
