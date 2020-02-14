# Subtype Descriptions

Analysis of the subtypes produced by hierarchical clustering on PC2:20, using cosine distance and dynamic tree cutting.
The deepSplit=0 clustering is used as the primary dataset, with additional subtypes created by deepSplit=1 noted.

## A - Alanine

A1 - Selection against anything larger than alanine. Buried

A2 - Selection against strong charge (and proline). Buried

A3 - Weak selection against aromatic and (iso)leucine - large hydrophobic AAs? More surface accessible

A4 - Anything but proline

A5 - Highly conserved, generally buried

A6 - Appears to improve with mutation, very surface accessible

A7 - Strong selection against tyrosine, weak selection against arginine, histidine, lysine, cysteine & glutamate.

A8 - Not methionine, otherwise weak selection against glycine & phenylalanine.

A4 and A6 are both in the first subcluster of subtypes, which is likely other permissive/anything but proline subtypes.

## C - Cysteine

C0 - Outliers, generally permissive with some selection against aromatics and proline.
Somewhat more likely to occur in helices compared to other subtypes?

C1 - Highly selective, FoldX/chemical environment suggests disulphide bonds generally go here, although not exclusively.
Possibly also contains other positions where only C functions.
Least selected against group is aromatics - maybe because the aromatic/sulphur interaction partially makes up for the disulphide bond?

C2 - Selective for hydrophobic amino acids, and generally not surface accessible. These positions likely rely on cysteine's hydrophobicity and possibly contains the group of sulphur/aromatic interactions, as their characteristic pattern occurs in the mean profile and somewhat in the chemical environment

## D - Aspartate

D1 - Highly selective, most permissive to E and N (both somewhat similar to D)

D2 - Selects against basic groups and larger side chains

D3 - Fully permissive, very accessible

D4 - Anything but proline

D5 - Appears to improve with substitution

D4 and D5 are both in the overall permissive subcluster, which seems to be subdivided into "anything but proline" and very permissive/improving subclusters.

## E - Glutamate

E1 - Most selective, most tolerated substitution is D

E2 - Generally permissive

E3 - Anything but proline

E4 - Everything is tolerated apart from large aromatics

E5 - Weak selection against some hydrophobic and basic residues, strong selection against proline and tyrosine.
Tolerates smaller hydrophobic residues (G, A, etc.) and polar side-chains as well as D

E6 - Mutations appear to improve

E7 - Generally Permissive, but most ER scores slightly positive

E2 and E7 appear similar (generally permissive) but since E7 is slightly positive it is grouped with E6 and E2 within the main cluster.

## F - Phenylalanine

F0 -  Small group of outliers, seems to want to mutate to aliphatic hydrophobic residues

F1 - Selective against charged/polar residues, tolerates aromatic and hydrophobic

F2 - Most selective, tolerates tyrosine. Least tolerant to charged and proline

F3 - Broadly tolerant, somewhat selective against cysteine and valine?

F1 and F2 are very intermixed in the true dendrogram, although F2's mean profile is more similar to F3's

## G - Glycine

G0 - Outliers, on average non-selective

G1 - Seems selective against everything but the smallest AAs.

G2 - Highly selective against everything, potentially needs to be small

G3 - Nonselective

G4 - Selective against larger hydrophobic residues and aromatics, otherwise nonselective

G5 - Strong selection against proline, weak selection against aromatics and isoleucine

All glycine residues appear enriched in turn and bend secondary structures

## H - Histidine

H1 - Selective, strong against proline and glycine, tolerant of Q, N and Y (all have polar groups)

H2 - Mostly nonselective, with weak selection against A, M, N and W (not a clear trend in these)

H3 - Other than proline, most substitutions improve fitness

Minimal clustering dendrogram and profile dendrogram are the same

## I - Isoleucine

I1 - Intolerant towards charged and polar residues, and proline. Tolerates hydrophobic substitutions

I2 - Highly selective, most tolerant to valine, leucine and methionine

I3 - Nonselective

I4 - Selective against negative charge and proline only

I5 - Selective against charged residues (both positive and negative), but not polar.

## K - Lysine

K1 - Weakly selective against everything except arginine, strongly selective against proline

K2 - Nonselective

K3 - Everything but proline

K4 - Selective against threonine and glutamine only

K5 - Substitutions generally improve fitness

## L - Leucine

L1 - Only tolerates hydrophobic residues and phenylalanine

L2 - Selective against everything but isoleucine, methionine and valine (inc. alanine, oddly)

L3 - Nonselective

L4 - Selective against charged, weakly selective against polar

L5 - Weak selection against isoleucine (?) and arginine

L6 - Strong selection against aspartate only (only weak against glutamate)

L7 - Tolerates aromatic and larger hydrophobic residues

L8 - Weak selection against aromatics only

L9 - Anything but proline

L10 - Strong selection against proline and arginine. Weak selection against negative charge, lysine and Y/W

## M - Methionine

M0 - Outliers, nonselective

M1 - Tolerates leucine, isoleucine, valine and threonine (longer aliphatic chains?).
Strongly intolerant to negative charge and proline

M2 - Moderately selection against everything, strong selection against cysteine

M3 - Strong selection against proline, weak selection against negative charge

## N - Asparagine

N0 - Outliers, on average seems to prefer being hydrophobic

N1 - Moderately selective against everything

N2 - Nonselective

N3 - Anything but proline

## P - Proline

P0 - Substitutions generally improve fitness

P1 - Selective against everything, strongest against aromatics

P2 - Mostly nonselective, some selection against arginine and glutamine

## Q - Glutamine

Q1 - Nonselective

Q2 - Weakly selective, tolerates basic and polar.
Particularly selective against aromatics

Q3 - Anything but proline

Q4 - Strong selection, only really tolerates H, K & E.
Strongest selection against proline and glycine

Q5 - Weak effects, mostly small improvement on substitution

Q6 - Strong selection against negative charge

## R - Arginine

R1 - Most selective, only tolerates lysine

R2 - Strongly selective against proline, weakly selective against negative charge

R3 - Nonselective

R4 - Intolerant to negative charge, proline, and hydrophobic residues

## S - Serine

S1 - Strongly selective, tolerates threonine, alanine and glycine.
Least tolerant to aromatics and hydrophobic residues

S2 - Nonselective

S3 - Selective against charged (both polarities) and proline

S4 - Anything by alanine (weak selection)

S5 - Anything but proline (strong selection)

S6 - Selective against negative charge and glutamine

S7 - Selective against aromatics and proline

## T - Threonine

T1 - Nonselective

T2 - Strongly selection, tolerates serine and alanine

T3 - Selective against negative charge and aromatics

T4 - Anything but proline

T5 - Selective against negative charge, proline and glycine

T6 - Tolerates hydrophobic residues

T7 - Strong selection against lysine and arginine

T8 - Fitness increases with other polar/charged residues

T9 - Strong selection against asparagine, weak selection against charged and aromatic

## V - Valine

V1 - Tolerates hydrophobic residues and threonine

V2 - Similar to V1 but doesn't tolerate threonine and weaker tolerance of non-valine hydrophobic residues

V3 - Nonselective

V4 - Intolerant of positive charge, proline and polar residues with NH2 groups (Q and N, which possibly tend to have positive partial charges on N?

V5 - Intolerant of aromatics and charged residues

V6 - Anything but proline

## W - Tryptophan

W1 - Tolerates other aromatics

W2 - Tolerates cysteine, glycine, leucine, arginine and serine (no clear pattern?)

## Y - Tyrosine

Y0 - Outliers, nonselective

Y1 - Generally selective, tolerates phenylalanine and histidine

Y2 - Tolerates aromatics and hydrophobic aliphatics

Subtypes are fairly intermixed in the original hierarchical clustering dendrogram.

