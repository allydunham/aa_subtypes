# Subtype Descriptions

Analysis of the subtypes produced by hierarchical clustering on PC2:20, using cosine distance and dynamic tree cutting.
The deepSplit=0 clustering is used as the primary dataset, with additional subtypes created by deepSplit=1 noted.

## A - Alanine

Summary:

* A1 - Nothing larger than alanine, buried
* A2 - Selection against strong charge (and proline), buried
* A3 - Weak selection against aromatic and (iso)leucine - large hydrophobic AAs? More surface accessible
* A4 - Anything but proline
* A5 - Highly conserved, generally buried
* A6 - Appears to improve with mutation, very surface accessible
* A7 - Strong selection against tyrosine, weak selection against arginine, histidine, lysine, cysteine & glutamate.
* A8 - Not methionine, otherwise weak selection against glycine & phenylalanine.

### A1

Appears to select against any residues more bulky than alanine, tolerating cysteine,
glycine etc.
It generally has a low surface accessibility.

* ADRB2 - 181 is inwards facing towards a phenylalanine (potential clash),and 76 & 78 are in the transmembrane helix bundle, although not obviously lacking space.
* amiE - 35, 39, 90, 94 are on the interface between two helices, all others are in a buried cluster of sheets and helices.
* APH3II - all internal to the globule, appear reasonably cramped
* CALM1 - 74 & 129 are internal, 47 & 104 are not, but maybe in the protein complex that CALM1 is part of.
* CBS - generally internal, often in secondary structure interfaces but not obviously cramped.
* CCR5 - All in helix interfaces
* CP - 31, 108 are near larger residues, all are internal and somewhat cramped based on spheres representation.
* CXCR4 - 162 is in a cramped helix interface, but others actually seem to be external, although as CXCR4 is a 7TM receptor they are membrane facing.
* HA - 85 is near a bulky aromatic, some others are buried and in helix interfaces although not all
* HSP90 - 41 is buried and cramped, but other positions less so
* MAPK1 - 171 is buried and cramped, 219 is in a helix adjacent to a tryptophan
* NP - Often between alpha helices, although one (218) is external to the structure (but maybe part of an interface in the full NP super-structure). 411 is very buried in the middle of a lobe, reasonably packed with some bulky residues nearby (several aromatics).
The lobe is actually in the middle of a different monomer, so this residue is likely part of the interaction stabilising the nucleoprotein multimer.
* PAB1 - Only position appears external and non-cramped
* Ras - 18 & 155 are buried in cramped helices, 66 & 122 appear exposed and unrestricted
* Src - Both positions, in particular 433, are cramped and buried.
* TEM1 - 133 is buried and cramped
* TP53 - 276 interfaces with DNA, 159 & 161 are in a cramped, buried strand
* TPK1 - 145 appears to be in a kink in a helix, potentially requiring small size. 39 is in a cramped area, 45 is buried
* TPMT - neither are that buried or cramped
* UBE2I - 133 is buried near a bulky aromatic, 106 appears exposed.

## C - Cysteine

Summary:

* C0 - Outliers, generally permissive with some selection against aromatics and proline.
Somewhat more likely to occur in helices compared to other subtypes?
* C1 - Highly selective, cysteine only functions such as disulphide bonds
* C2 - Selective for hydrophobic amino acids, generally buried, includes cysteine/aromatic interactions.

### C1

FoldX and chemical environment results suggest disulphide bonds generally go here, although not exclusively.
The subtype likely also contains other positions where only C functions.
The subtype is least selected against aromatics - maybe because the aromatic/sulphur interaction partially makes up for the disulphide bond?

### C2

This subtypes appears to use cysteine's hydrophobicity as well as the sulphur/aromatic interaction that allows cysteine to create stabilising interactions with aromatics (Orabi & English 2016)
The characteristic aromatic interaction pattern occurs in the subtypes mean profile and somewhat in the chemical environment, although the enrichment for aromatic neighbours is weaker than when it was an exclusive group (likely because it is swamped by residues without the interaction).
In some sense this subtype could be manually split into two, although that is generally not identified algorithmically (apart from a single lucky run of kmeans)

The identified aromatic interactions are:

* Src 501 - In a helix bundle surrounded by hydrophobic residues, with an aromatic nearby but not adjacent
* TEM1 75 - Also in a bundle and near an aromatic ring, but also one half of a disulphide bond (the other half is not in C4)
* MAPK1 216 - Buried and near an aromatic
* MAPK1 161 - Exterior but near an aromatic
* MAPK1 166 - In an apparent binding cleft, appears to conjugate the bound molecule in this structure (which is aromatic), but is also a modified residue in the crystal structure
* ADRB2 327 - Semi-buried residue near hydrophobic residues and an aromatic
* TPMT 133 - Residue is on the exterior but the side chain faces inwards into a hydrophobic pocket
* UBE2I 43/75 - Near each other in 3D space, both side chains face into the protein and there are hydrophobic and aromatic residues nearby, although not directly conjugating like in some examples.
* APH3II 209 - near to a phenylalanine
* APH3II 209 - In a cleft with a NA ion and an aromatic drug in the crystal structure
* CCR5 213 - directly conjugates a tyrosine ring
* TPK1 104 - Directly conjugates a phenylalanine ring
* TPK1 88 - Near but not directly adjacent to a tyrosine ring

## D - Aspartate

Summary:

* D1 - Highly selective, most permissive to E and N (both somewhat similar to D)
* D2 - Selects against basic groups and larger side chains
* D3 - Fully permissive, very accessible
* D4 - Anything but proline
* D5 - Appears to improve with substitution

D4 and D5 are both in the overall permissive subcluster, which seems to be subdivided into "anything but proline" and very permissive/improving subclusters.

## E - Glutamate

Summary:

* E1 - Most selective, most tolerated substitution is D
* E2 - Generally permissive
* E3 - Anything but proline
* E4 - Everything is tolerated apart from large aromatics
* E5 - Weak selection against some hydrophobic and basic residues, strong selection against proline and tyrosine.
Tolerates smaller hydrophobic residues (G, A, etc.) and polar side-chains as well as D
* E6 - Mutations appear to improve
* E7 - Generally Permissive, but most ER scores slightly positive

E2 and E7 appear similar (generally permissive) but since E7 is slightly positive it is grouped with E6 and E2 within the main cluster.

## F - Phenylalanine

Summary:

* F0 -  Small group of outliers, seems to want to mutate to aliphatic hydrophobic residues
* F1 - Selective against charged/polar residues, tolerates aromatic and hydrophobic
* F2 - Most selective, tolerates tyrosine. Least tolerant to charged and proline
* F3 - Broadly tolerant, somewhat selective against cysteine and valine?

F1 and F2 are very intermixed in the true dendrogram, although F2's mean profile is more similar to F3's

## G - Glycine

Summary:

* G0 - Outliers, on average non-selective
* G1 - Seems selective against everything but the smallest AAs.
* G2 - Highly selective against everything, potentially needs to be small
* G3 - Nonselective
* G4 - Selective against larger hydrophobic residues and aromatics, otherwise nonselective
* G5 - Strong selection against proline, weak selection against aromatics and isoleucine

All glycine residues appear enriched in turn and bend secondary structures

## H - Histidine

Summary:

* H1 - Selective, strong against proline and glycine, tolerant of Q, N and Y (all have polar groups)
* H2 - Mostly nonselective, with weak selection against A, M, N and W (not a clear trend in these)
* H3 - Other than proline, most substitutions improve fitness

Minimal clustering dendrogram and profile dendrogram are the same

## I - Isoleucine

Summary:

* I1 - Intolerant towards charged and polar residues, and proline. Tolerates hydrophobic substitutions
* I2 - Highly selective, most tolerant to valine, leucine and methionine
* I3 - Nonselective
* I4 - Selective against negative charge and proline only
* I5 - Selective against charged residues (both positive and negative), but not polar.

## K - Lysine

Summary:

* K1 - Weakly selective against everything except arginine, strongly selective against proline
* K2 - Nonselective
* K3 - Everything but proline
* K4 - Selective against threonine and glutamine only
* K5 - Substitutions generally improve fitness

## L - Leucine

Summary:

* L1 - Only tolerates hydrophobic residues and phenylalanine
* L2 - Selective against everything but isoleucine, methionine and valine (inc. alanine, oddly)
* L3 - Nonselective
* L4 - Selective against charged, weakly selective against polar
* L5 - Weak selection against isoleucine (?) and arginine
* L6 - Strong selection against aspartate only (only weak against glutamate)
* L7 - Tolerates aromatic and larger hydrophobic residues
* L8 - Weak selection against aromatics only
* L9 - Anything but proline
* L10 - Strong selection against proline and arginine. Weak selection against negative charge, lysine and Y/W

## M - Methionine

Summary:

* M0 - Outliers, nonselective
* M1 - Tolerates leucine, isoleucine, valine and threonine (longer aliphatic chains?).
Strongly intolerant to negative charge and proline
* M2 - Moderately selection against everything, strong selection against cysteine
* M3 - Strong selection against proline, weak selection against negative charge

## N - Asparagine

Summary:

* N0 - Outliers, on average seems to prefer being hydrophobic
* N1 - Moderately selective against everything
* N2 - Nonselective
* N3 - Anything but proline

## P - Proline

Summary:

* P0 - Substitutions generally improve fitness
* P1 - Selective against everything, strongest against aromatics
* P2 - Mostly nonselective, some selection against arginine and glutamine

## Q - Glutamine

Summary:

* Q1 - Nonselective
* Q2 - Weakly selective, tolerates basic and polar.
Particularly selective against aromatics
* Q3 - Anything but proline
* Q4 - Strong selection, only really tolerates H, K & E.
Strongest selection against proline and glycine
* Q5 - Weak effects, mostly small improvement on substitutionQ6 - Strong selection against negative charge

## R - Arginine

Summary:

* R1 - Most selective, only tolerates lysine
* R2 - Strongly selective against proline, weakly selective against negative charge
* R3 - Nonselective
* R4 - Intolerant to negative charge, proline, and hydrophobic residues

## S - Serine

Summary:

* S1 - Strongly selective, tolerates threonine, alanine and glycine.
Least tolerant to aromatics and hydrophobic residues
* S2 - Nonselective
* S3 - Selective against charged (both polarities) and proline
* S4 - Anything by alanine (weak selection)
* S5 - Anything but proline (strong selection)
* S6 - Selective against negative charge and glutamine
* S7 - Selective against aromatics and proline

## T - Threonine

Summary:

* T1 - Nonselective
* T2 - Strongly selection, tolerates serine and alanine
* T3 - Selective against negative charge and aromatics
* T4 - Anything but proline
* T5 - Selective against negative charge, proline and glycine
* T6 - Tolerates hydrophobic residues
* T7 - Strong selection against lysine and arginine
* T8 - Fitness increases with other polar/charged residues
* T9 - Strong selection against asparagine, weak selection against charged and aromatic

## V - Valine

Summary:

* V1 - Tolerates hydrophobic residues and threonine
* V2 - Similar to V1 but doesn't tolerate threonine and weaker tolerance of non-valine hydrophobic residues
* V3 - Nonselective
* V4 - Intolerant of positive charge, proline and polar residues with NH2 groups (Q and N, which possibly tend to have positive partial charges on N?
* V5 - Intolerant of aromatics and charged residues
* V6 - Anything but proline

## W - Tryptophan

Summary:

* W1 - Tolerates other aromatics
* W2 - Tolerates cysteine, glycine, leucine, arginine and serine (no clear pattern?)

## Y - Tyrosine

Summary:

* Y0 - Outliers, nonselective
* Y1 - Generally selective, tolerates phenylalanine and histidine
* Y2 - Tolerates aromatics and hydrophobic aliphatics

Subtypes are fairly intermixed in the original hierarchical clustering dendrogram.
