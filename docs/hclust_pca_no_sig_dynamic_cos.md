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

C2 - Selective for hydrophobic amino acids, and generally not surface accessible. These positions likely rely on cysteine's hydrophobicity and possibly contains the group of sulphur/aromatic interactions, as their characteristic pattern occurs in the mean profile and somewhat in the chemical environment.

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
Tolerates smaller hydrophobic residues (G, A, etc.) and polar side-chains as well as D.

E6 - Mutations appear to improve

E7 - Generally Permissive, but most ER scores slightly positive

E2 and E7 appear similar (generally permissive) but since E7 is slightly positive it is grouped with E6 and E2 within the main cluster.

## F - Phenylalanine

## G - Glycine

## H - Histidine

## I - Isoleucine

## K - Lysine

## L - Leucine

## M - Methionine

## N - Asparagine

## P - Proline

## Q - Glutamine

## R - Arginine

## S - Serine

## T - Threonine

## V - Valine

## W - Tryptophan

## Y - Tyrosine
