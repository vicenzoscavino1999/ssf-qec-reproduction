Exact data and check fibres in one noisy CSS extraction round
Companion to Appendix C.6, 23 November 2026

Extract this archive and run the commands below from its root with Python
3.10 or later and its standard library. Use the default assertion settings.
The tests use fractions. Decimal output is rounded for display.

General proposition
    python -I -B -X utf8 verify_fibres.py

This checks the projected CNOT law by comparing all primitive Pauli
categories with the independent-fibre construction on a two-wire
circuit. It compares joint laws and signed, physically filtered sums for
every primitive/action subhost. It also evaluates the activity test at the
directly transferred (3,4) point and at the third measurement corner.
The size-independent proposition and its proof are in the manuscript.

Uniform (3,4) self-products

Transfer from p0=3.2e-7, q0=1e-3 to rc=4.2e-8, rm=1e-3
    python -I -B -X utf8 verificar_FASE5_circuito_propagacion.py

Table 5, measurement endpoint 1e-3
    python -I -B -X utf8 verificar_H5A_3_4_primer_apagado_de_ternas.py --global

Table 5, measurement endpoint 1e-4
    python -I -B -X utf8 verificar_H5A_3_4_esquina_2.py

Table 5, measurement endpoint 1e-5
    python -I -B -X utf8 verificar_H5A_3_4_esquina_3.py
    python -I -B -X utf8 verificar_H5A_3_4_esquina_3_circuito.py

The last two commands cover Bernoulli and both circuit orders. The circuit
run checks both degree tilts. An open numerical witness is recorded with
its largest failing ratio. The smaller rectangles listed in Table 5 have
closed witnesses. The relevant Spanish derivations in proofs/ have the
same suffixes as the verifiers. The direct transfer proof is in
FASE5_CIRCUITO_3_4_PROPAGACION.txt, sections 2 through 5.

Uniform (17,18) self-products and the single affine F17 instance
    python -I -B -X utf8 verificar_FASE6_F17_circuito.py

Section 2 of proofs/FASE6_F17_TRANSFERENCIA_PQ_A_RUIDO_CIRCUITAL.txt
derives the finite-volume certificate. Sections 3 through 6 derive the
uniform certificate. The verifier uses the original
SOURCE_CUMULANT_CERTIFICATE.zip from this directory and checks its hash.
That embedded archive also supplies the underlying Bernoulli source proof
and its own verify.py entry point.

First contact, Lemma 14 (section 5.3 and appendices C.6 and N.1)
    python -I -B -X utf8 verify_first_contact.py
The manuscript gives the degree-uniform proof. The verifier checks the local
inequalities on 815 degree pairs and the full histories of the (2,3), (2,4)
and (3,6) boundary examples. For the interior (7,10) and (8,11) examples it
constructs the seed matrices, verifies regularity, distinct rows and columns,
and isolation in the physical host. It reconstructs the observed syndromes
and checks the fixed order and the stable endpoint. It computes every maximum
by sorting row gains and independently by enumerating eligible XOR footprints.
