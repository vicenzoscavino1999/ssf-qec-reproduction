Circuit fibers and first contact
================================

The tasas action in ../../REPRODUCTION.md runs the checks listed below.
To inspect a single result, run the corresponding script from this
directory with Python and assertions enabled. For example, the following
command reproduces the one-round fiber calculation.

    python -B verify_fibres.py

The script compares the joint circuit response and signed filtered sums
on all 512 subhosts of a two-wire circuit. Proposition 24 in Appendix C.5
gives the factorization for the stated one-round binary CSS model.

Uniform rates in Table 1
------------------------
verificar_FASE5_circuito_propagacion.py checks the direct transfer from
p = 3.2e-7 and q = 1e-3 to r_c = 4.2e-8 and r_m = 1e-3 in the (3,4)
family. To reproduce the sharper schedule-dependent rows of Table 1, use
the following scripts for the indicated upper bound on q or r_m.

    1e-3    verificar_H5A_3_4_primer_apagado_de_ternas.py --global
    1e-4    verificar_H5A_3_4_esquina_2.py
    1e-5    verificar_H5A_3_4_esquina_3.py
            verificar_H5A_3_4_esquina_3_circuito.py

At 1e-5 the first script gives the Bernoulli row. The second evaluates
the VV-first and admitted-order circuit bounds at both degree tilts.
The derivations in proofs/ share the script suffixes. The direct-transfer
argument is in sections 2–5 of proofs/FASE5_CIRCUITO_3_4_PROPAGACION.txt.

For the (17,18) family, verificar_FASE6_F17_circuito.py checks the
uniform transfer and the stronger bound for the finite affine F17
instance. Its Bernoulli input is SOURCE_CUMULANT_CERTIFICATE.zip.
The finite and uniform arguments are in sections 2 and 3–6, respectively,
of proofs/FASE6_F17_TRANSFERENCIA_PQ_A_RUIDO_CIRCUITAL.txt.

First contact
-------------
Lemma 19 and its examples are in Appendix A.3. The original verifier calls
this result Lemma 14. verify_first_contact.py
checks the local inequalities on 815 degree pairs, the boundary cases
(2,3), (2,4) and (3,6), and the seed matrices and full decoder histories
for the (7,10) and (8,11) examples. It computes maxima by sorted row
gains and by enumeration of the eligible XOR footprints.
