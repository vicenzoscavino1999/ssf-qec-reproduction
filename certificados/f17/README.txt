Uniform (17,18) source-cumulant certificate
=========================================

The Bernoulli rectangle for the whole family is p <= 5e-33 and q <= 2e-11.
From this directory, run the following command with Python 3.10 or later
and assertions enabled. The verifier imports modules from the ZIP inside
certificate/ and uses only the standard library.

    python verify.py

expected.txt gives the reference output. The three source rows are
printed in VV, CC and measurement order. The last two values bound
the higher-order remainders for the data and measurement rows at the
upper corner. The verifier uses integer and rational arithmetic.

Derivations in proofs/
----------------------
H5A_fusiones_sin_repetir.txt, sections 1–2, states the joint bound.
H5A_p_general.txt defines the chronological construction and support-size
bounds. H5A_interacciones_por_fuentes.txt and H5A_primera_fusion.txt
treat isolated pieces. H5A_consumo_y_alcance.txt gives the first-crossing
inequalities. H5A_contactos_entre_ejecuciones.txt, sections 1–5, derives
the causal encoding and the fifteen tree inequalities.
H5A_primer_apagado_por_grupos.txt, sections 1–3, bounds the marked
series and its tail. H5A_cumulantes_y_componentes.txt, sections 1–2,
gives the normalized partition used in the approximation.
