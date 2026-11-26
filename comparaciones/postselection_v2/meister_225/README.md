# A cluster-distance pilot on the affine 225-qubit code

At independent Bernoulli rates \(p=10^{-6}\) and \(q=10^{-4}\), Lee's published cluster score gives lower silent risk per retained round than the geometry score defined below. Both rules see the same noisy syndrome and use the same BP+LSD clusters. The physical SSF decoder performs every correction under the manuscript's fixed tie order.

| Rule at matched discard | Deterministic interval for silent risk per retained round |
| --- | --- |
| Lee cluster score at the calibrated cutoff | \([4.27266\times10^{-28},\,1.68890\times10^{-25}]\) |
| HGP adaptation of Meister's cluster distance | \([3.20481\times10^{-20},\,1.65255\times10^{-18}]\) |

The common discard probability lies between 2.58762 and 2.60134 per million attempts. Dividing the unrounded lower endpoint for the geometry score by the upper endpoint for Lee gives a risk factor greater than 189,750. These are bounds on the true probabilities of the specified rules on this finite code.

## What the geometry score measures

For each recorded syndrome \(s\), collect the data mechanisms in the final active BP+LSD clusters into a set \(C(s)\). Every data fault has the same log-odds cost here. In units of that cost the score is

\[
g(s)=\min_{L\in\ker H_X\setminus S_X}|\operatorname{supp}L\setminus C(s)|.
\]

Small values indicate that the clustered data nearly cover a nontrivial physical \(X\) logical. Readout faults influence cluster growth. They contribute zero physical data support to a logical operator. This choice adapts the zero cost of sites inside clusters in Meister's surface-code construction to the silent-failure observable of a one-round HGP. [Meister, Pattison and Preskill, Section III.1](https://arxiv.org/html/2405.07433v2#S3.SS1).

The calculation includes all logical classes and stabilizer representatives. Write a zero-syndrome error on the product as \((A,B)\), with \(A\) the VV block and \(B\) the CC block. Then \(AH_0^T+H_0^TB=0\). The nine logical tags imply that for a nontrivial class there is a seed kernel word \(b\) with \(b^TA\ne0\). That contraction is a nonzero seed word. Its occupied columns give a lower bound on the number of paid VV sites. Single-row logicals give the upper bound. When \(B=0\), each nonzero row of \(A\) is a seed word. When \(B\ne0\), there are at least six VV sites and one CC site. These bounds agree for all 80,194 distinct data-cluster supports encountered. A separate binary elimination and search over up to four paid columns checked 87 representative cases against the geometric computation.

## Matching the native Lee discard probability

Let \(a_L\) be the discard probability of the published Lee score at the cutoff calibrated for this code, \(3368/10^6\). Discard every \(g<5\), retain every \(g>5\), and at \(g=5\) discard with an independent coin of probability

\[
\alpha=\frac{a_L-\Pr(g<5)}{\Pr(g=5)}.
\]

This gives equality of the two discard probabilities in expectation over the input and the coin. The calibration uses the input distribution alone. Exact enumeration by counts of data and readout faults covers \((0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(2,0)\). Assigning the remaining input mass in either direction gives \(0.0178794\le\alpha\le0.0180717\). The risk bounds cover the entire interval. The decisive examples have \(g=6\), so how the coin falls leaves them all retained.

## What separates the rules

Both rules discard every silent history in sectors \((2,1),(2,2),(3,0)\), with 324, 32,076 and 36 histories respectively. For \((3,1)\), expand the existing affine orbits that allow silent failure under some maximum-score tie. The expansion yields exactly 101,088 distinct configurations. Two implementations of physical SSF agreed on the failure decision for every candidate. Of 59,184 silent histories, Lee discards all and the geometry rule retains 324. Alone these give a joint silent mass of \(324p^3(1-p)^{222}q(1-q)^{107}\). For example, data sites 0, 3, 4 and wrong readout 3 give lit checks 1, 4, 5, 8. BP+LSD places only readout mechanisms in the active clusters for this input. Its geometry score is six while its Lee score triggers discard. Site indices start at zero.

The tighter Lee upper bound follows from its cutoff. Let \(w_p=\log 999999\), \(w_q=\log 9999\). Their ratio lies between \(3/2\) and \(8/5\), so the total weight \(225w_p+108w_q<468w_q\). If the clusters contain three reads, two data sites, or one of each, the cluster weight 2-norm is at least \(\sqrt{3}w_q\), even when they are in separate clusters. Its normalized value exceeds the cutoff. BP+LSD assembles a syndrome-consistent solution on the union of its active clusters. Thus every accepted syndrome is zero, one or two individual readout bits, or one data column. Both SSF implementations verified on all 6,112 such syndrome candidates that the correction is stabilizer-equivalent to at most one physical data site. [Lee, English and Bartlett, Sections 2.2.2 and 4.2](https://arxiv.org/html/2510.05795v4).

A retained silent failure under Lee therefore needs at least five erroneous data sites since \(d_X=6\). Paying the entire five-data tail already makes the geometry risk at least 6.97 times the Lee risk. For the tighter bound, at data weight five the single-site correction must complete a weight-six logical. The seed projection above forces every weight-six logical onto a single VV row, with six seed words for each of the twelve rows. Removing one site gives 432 distinct possible five-data supports. All 432 give retained silent examples at zero readout faults. If \(A\sim\mathrm{Bin}(225,p)\), the joint retained silent mass under Lee is consequently between \(432p^5(1-p)^{220}(1-q)^{108}\) and \(432p^5(1-p)^{220}+\Pr(A\ge6)\). Divide by \(1-a_L\) for retained-round risk.

An additional paired count checks the direction of the result using the archived sector tails. In sector \((2,3)\) physical SSF finds 1,300,536 silent histories among the 1,571,724 candidates allowed by the row lemma. Charge all of their mass in the tail. Since \(g<5\) requires at least two clustered data sites, Lee always discards in that case. On the unclassified tail the geometry rule can gain over Lee only through the boundary coin. Even allowing it that entire potential gain leaves a strictly positive difference in retained risk of at least \(2.7629\times10^{-21}\) in Lee's favor.

The same census tightens the Bernoulli row for the seven-check rule. It keeps 324 silent \((3,1)\) inputs with four lit checks and 864 with six. Replacing the previous upper count of 101,088 for this sector by the actual 1,188 kept gives a risk interval of \([1.1750\times10^{-19},4.4631\times10^{-19}]\), with 0.0133897–0.0133917 discards per million.

## Reproduce

From tag `circuit-postselection-v4`, install the pinned Python 3.12 dependencies described in the [reproduction guide](../../../REPRODUCTION.md). From the repository root run

~~~sh
python reproducir.py meister-225
~~~

The launcher builds the input catalogue and scores 347,548 distinct syndromes using one BP+LSD thread. It writes `meister-225.json` in the run directory and compares the scientific output with `results.json` next to this guide. The original isolated run took 131 seconds on the Windows CPU. Stage timings and the processor are recorded separately for each run.
