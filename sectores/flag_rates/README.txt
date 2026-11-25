Finite-code pilot for a one-round SSF abort flag and higher Bernoulli rates
=======================================================================

Inputs and scope
----------------
Use the paper's ssf-qec-reproduction checkout at tag circuit-postselection-v2.
Python 3.10+ is sufficient
to check the arithmetic from the archived tallies. To rebuild them, use
a C++17 compiler, as for the paper's sector census.

From the integrated repository, python reproducir.py marca-tasas runs the
full Bernoulli census. The quick check python reproducir.py rapido recomputes the
same intervals from the stored tallies. Each action writes marca-tasas.json
to its output directory. When using this standalone folder, run either of
the following commands from the folder containing run_pilot.py.

    python -B run_pilot.py path/to/ssf-qec-reproduction
    python -B run_pilot.py path/to/ssf-qec-reproduction --recount

The second command regenerates the two seeds, physical decoder tables, all
new direct and orbit tallies, and checks the low-weight 225-qubit flag census
with the second Python implementation. It also verifies the 324 affine row
patterns among all 7,920 reduced inputs against that decoder and reruns
the 292-representative single-readout census giving C31 <= 101088. Both modes recompute pilot-results.json
in a temporary directory and compare all reported values. The present bounds
apply to the two specified distance-six codes with the prescribed SSF/FGL
tie order, independent data and readout bits, one round, and the X sector.

A. Abort flag
-------------
Abort if SSF made at least one move and the observed stopping syndrome is
nonzero. Strict descent lets one implement the rule as 0 < |s_stop| < |s_in|.
For silent failure the physical residual is a nontrivial logical, hence
s_stop = D. With nonempty D and no move the input data error must itself be
logical and have at least six sites. The included audit exhibits data sites
0..5 and the single flipped readout 0 on the 225-qubit code.

The exact weight-at-most-two flag counts plus the whole independent-binomial
tail bound Pr(abort). For D=0 use exact C30, exact C40 on 225 and upper C40 on
900. Bound higher data weights by the binomial tail. At 900, C50 also has a
finite upper count. The positive-D, empty-trace contribution is at most
Pr(Bin(n,p)>=6) [1-(1-q)^m]. Divide the joint unflagged upper endpoint by
1-Pr(abort)_upper to bound the conditional risk per kept round. Input-zero
post-selection rejects whenever the initial observed syndrome is nonempty.
The discard comparison measures yield and does not equate residual risks.

For an odds-region bound, divide the same unflagged estimate by T30.
Its upper coefficients increase with x=p/(1-p) and y=q/(1-q). Use the
positive lower bound T30+T21 for the total failure probability. This gives
a uniform lower gain on 0<y<=y0, 0<x/y<=x0/y0, anchored at (1e-6,1e-4).

B. Probability interval at higher rates
--------------------------------------
All arithmetic uses fractions. For 225 combine the necessary nine-readout
row pattern (comparaciones/H5A_piloto_HGP_225.txt, sections 2 to 4)
with sectores/intervalo.py. The total exactly-two-data upper bound is
324 p^2(1-p)^223 q(1-q)^8. Subtract the known (2,1) and (2,2) masses
before bounding the (2,>=3) remainder. The (3,1) upper count is 101088 in comparaciones/piloto_HGP_225_esquinas_20261124.json.
The direct 103962600-configuration data-only census makes C40=14904
an exact contribution to both ends of the interval.

The small-rate rectangle uses the lower sum T21+T30. After factoring out
(1-p)^225 (1-q)^108, the uncertain readout terms are bounded by T21 times
a polynomial with nonnegative coefficients in the odds. For the pair tail
the polynomial is (1+y)^99-1-99y. Charge the >=5-data union bound to T30.
Both normalized majorants increase up to the upper corner. The harmonic
relative error for [L,U] is (U-L)/(U+L). If p=0 the probability is zero.
See rectangle225() for the individual charges.

Orbit bounds used in the high-rate 225-qubit interval and the 900-qubit flag
--------------------------------------------------------------------------
Generate ordered coordinate tuples separately on the two HGP axes under the
verified seed automorphisms. Quotient by permutations within VV sites, CC
sites and readouts. The product of axis orbit sizes divided by the internal
permutation stabilizer gives each unordered orbit's population. Total masses
must equal binom(n,a) binom(m,b) for the sector.

Allow every first move tied after score and size maximization. If its
syndrome and logical tag equal a nonempty subset of the initial faults,
the prescribed continuation is that of the remaining r data faults with
the same readouts. Charge these histories by binom(n-r,a-r) C_rb, summed
over r<a. At distance six, r=0 and r=1 contribute zero here. For all
other first moves explore tied descendants. A terminal silent state must
have observed syndrome D and a nonzero residual logical tag. D must be
stable. Below |D| no later state can be silent by strict descent. An
unfinished branch is counted as potentially silent. An orbit is excluded
only when every exceptional path fails. Zero/nonzero logical class, syndrome,
score and action size are equivariant for this relaxed search. The physical
mask tie used by the real decoder need not be invariant. The procedure
reproduces the two previously published exceptional three-data/one-readout
counts, 93312 and 1260288.

For five data and D=0, first charge genuine moves to the known three-
and four-data failures. If no original fault is fully lit, a distinguished
original overlaps another, with probability at most 4*24/(n-1) in this
sector. Otherwise a first move outside all original subsets must have a
fully lit footprint. For each generator, enumerate its local masks. Each
of the twelve cells has five external data neighbours, each incident to
one cell at most, verified directly by check_geometry(). For a footprint
of d cells, u of which require odd external parity, count the missing
external faults with the coefficient of degree 5-|I| in

 (5t+10t^3+t^5)^u (1+10t^2+5t^4)^(d-u) (1+t)^(n-7-5d).

I is the set of faulty internal data qubits. Omit patterns containing a
representative of the move or its complement in the generator. Summing the
remaining patterns over every generator and footprint gives an upper count
without dividing out duplicate generators. This also checks against the
released local footprint majorant. Charge all six-or-more data faults in
the zero-readout class by the exact binomial tail.

Circuit postselection on 225 qubits
----------------------------------
From the complete checkout, python reproducir.py marca-circuito rebuilds the
eight one-round circuit intervals with the two gate orders of manuscript
Table 6. The action checks its result against circuit-flag.json. The full
reproduction guide describes the output fields and the deterministic inputs.

Main entries in pilot-results.json
----------------------------------
flag225_direct and orbit_bounds.900_at_reference give the abort, unflagged
and conditional intervals at p=1e-6, q=1e-4. orbit_bounds.225_at_high_rates
gives the refined interval at p=1e-4, q=1e-3. finite_region_225 gives the
small-rate rectangle. For the odds-cone abort result use flag_regional.225
and orbit_bounds.900_gain_on_cone. Earlier comparison bounds are retained
under separate keys. Values are rounded outwards where used as bounds.
