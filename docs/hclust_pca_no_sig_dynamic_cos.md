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

### A4

Tolerates any substitution apart from proline.

Examples:

* MAPK1 - 3 positions (143, 309 & 352) are at the C terminus of alpha helices, where the proline kink likely destabilises the helix. 281 is in the turn between two helices. 174 is in an apparently free coil.
* TP53 - 3 position (119, 129 & 138) are in the turns of a beta-sheet. 347 is in the middle of an alpha-helix

## C - Cysteine

Summary:

* C0 - Outliers, generally permissive with some selection against aromatics and proline.
Somewhat more likely to occur in helices compared to other subtypes?
* C1 - Highly selective, cysteine only functions such as disulphide bonds
* C2 - Selective for hydrophobic amino acids, generally buried, includes cysteine/aromatic interactions.

### C1

FoldX and chemical environment results suggest disulphide bonds generally go here, although not exclusively (24/73 positions vs 3/37 in C2).
The subtype likely also contains other positions where only C functions, such as ion conjugation, and potentially includes similar aromatic interactions to C2 where the interaction is highly important.
The subtype is least selected against aromatics - maybe because the aromatic/sulphur interaction partially makes up for the disulphide bond, or also fulfils some of the other functions.

Examples without disulphide bond:

* NP - 44 & 279 are adjacent to each other with a tyrosine ring between them - potentially some interaction through the aromatic? 223 is near a tyrosine, potentially having an aromatic interaction.
* TP53 - 176, 238 & 242 all appear to be conjugating a Z ion (Pace and Weerapana, 2014). Others unclear
* CCR5 - 224 is on the receptor surface, potentially with an active role
* DRB2 - Several membrane facing cysteine residues on the transmembrane helices?
* GAL4 - 4 residues (14, 21, 28, 38) conjugating two zinc ions
* APH3II - 131 appears to be interacting with a tryptophan residue
* BRCA1 - 39 & 64 and 27 & 47 each conjugate a zinc ion
* CBS - 52 is bonded to the Fe of a haem group in the crystal structure (which is a known binding in the native protein)
* SUMO1 - 52 is directed at a phenylalanine

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

### D1

Often appear to have a specific role (see examples), although not always.
In other cases it could simply be we cannot easily see the role.

Examples:

* HA - 93 is directly adjacent to an arginine (likely ionic interaction?)
* CBS - 221 is directly adjacent to an arginine. 49, 309 & 198 are all grouped around a phenylalanine aromatic group; possibly they induce a dipole? 388 is near an arginine
* TEM1 - 189, 208 neighbour a K ion
* HSP90 - All four (40, 79, 113, 143) all border an apparent ATP binding pocket, as do 3 D3 positions
* PAB1 - 184 possibly bonding to lysine 180, 138 has a similar arrangement with lysine 140

### D4

Tolerates any substitution other than proline.

Examples:

* APH3II - 118 & 238 are in helices, 23 & 52 at the start and end of beta strand turns
* PAB1 - 144 is in a helix, 136 & 160 are in coils around this helix
* TPMT - 31 & 162 are both in helices

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

### E1

As with D1, E1 positions often seem to have a recognisable role utilising the charge properties of glutamate, although in many cases it is less obvious too.

Examples:

* ADRB2 - Not clear what role they play but 188 and 306 are both near the extracellular binding site
* amiE - 105 & 108 both extend towards a similar pair of lysine and arginines. 248, 249 & 250 are arrayed at the interface between subunits, next to the equivalent residues on the other subunit, creating a very negative pocket. Are many other positions where role isn't clear
* APH3II - 262 is in the cleft where kanamycin is bound (APH3II is a kanamycin resistance protein). 182 is at the negative end of an alpha helix, possible stabilising its dipole (Sali et al. 1988).
* UBE4B - 1083 and 1084 are on a very polar exterior helix face (2 other E, 2 K and 1 R)

## F - Phenylalanine

Summary:

* F0 - Small group of outliers, seems to want to mutate to aliphatic hydrophobic residues
* F1 - Selective against charged/polar residues, tolerates aromatic and hydrophobic
* F2 - Most selective, tolerates tyrosine. Least tolerant to charged and proline
* F3 - Broadly tolerant, somewhat selective against cysteine and valine?

F1 and F2 are very intermixed in the true dendrogram, although F2's mean profile is more similar to F3's

### F2

Seems to utilise specific aromatic properties as well as hydrophobicity.
They occur relatively frequently in groups.

Examples:

* BRCA1 - 43 surrounded by hydrophobic residues and filling a ring shaped space
* CALM1 - 17 & 69 are potentially pi-stacking
* CBS - 332, 385 & 396 are positioned in a cleft together, along with several other aromatic residues; more potential pi-stacking? 197 & 310 are also adjacent to each other
* TP53 - 270 & 113 are adjacent and face one another. 134 is near them but in the middle of an apparently polar pocket?

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

### M1

Examples:

* HSP90 - 84 is in a binding pocket for ATP, 116 fills a hydrophobic cleft, 16 looks to interact with a phenylalanine ring, 589 is near to a tryptophan residue (interestingly the 589 position is classified M1 from the profile in Hietpas et al. 2011 and M3 when from Jiang et al. 2013).
* MAPK1 - 221 & 293 appear to be filling hydrophobic space
* NP - series of positions in a helix bundle
* TP53 - 243 faces towards bound DNA, 133 is in a space near two aromatics, 160 seems to fill a large hole
* Src - both fill holes near aromatics

### M2

Examples:

* Src - 317 in a hydrophobic pocket near two aromatics
* Ras - 72 in a pocket near a phenylalanine
* NP - Clustered in a different helix bundle to M1 positions, but not obvious what the difference is other than potentially more aromatics
* HA - 460 in a fairly tight looking pocket near a tyrosine.
* CBS - 464 is in a tight pocket, 458 in a similar pocket but with two aromatics, 337 is near two aromatics in an exposed cleft
* CALM1 - 72 & 73 are in a helix helix interface adjacent to two M1 positions

The physical difference between M1 and M2 is not completely clear.
They have similar distributions of distance to the nearest aromatic.
To some extent M1 positions appear to physically group together in some proteins and M2 positions group separately, potentially due to physical characteristics or the importance of the two domains.

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
