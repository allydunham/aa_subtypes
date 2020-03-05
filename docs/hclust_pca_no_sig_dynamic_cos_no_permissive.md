# Subtype Descriptions

Analysis of the subtypes produced by hierarchical clustering on PC2:20, using cosine distance, dynamic tree cutting and with permissive positions (all |ER| < 0.4) excluded into their own subtype before the main clustering step.
The deepSplit=0 clustering is used as the primary dataset, with additional subtypes created by deepSplit=1 noted (WIP).

## A - Alanine

Summary:

* A1 -
* A2 -
* A3 -
* A4 -
* A5 -
* AP - Permissive positions
* AO - Outliers

### A1

### A2

### A3

### A4

### A5

### AP

### AO

## C - Cysteine

Summary:

* C1 - Highly selective, cysteine only functions such as disulphide bonds
* C2 - Selective against hydrophilic amino acids, much more buried, includes cysteine/aromatic interactions.
* CP - Permissive positions
* CO - Outliers

### C1

<!-- FoldX and chemical environment results suggest disulphide bonds generally go here, although not exclusively (24/73 positions vs 3/37 in C2).
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
* SUMO1 - 52 is directed at a phenylalanine -->

### C2

<!-- This subtypes appears to use cysteine's hydrophobicity as well as the sulphur/aromatic interaction that allows cysteine to create stabilising interactions with aromatics (Orabi & English 2016)
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
* TPK1 88 - Near but not directly adjacent to a tyrosine ring -->

### CP

### CO

<!-- Two general categories; some pattern of negative substitutions and permissive positions that tend to be somewhat improved on mutation.
Positions with positive ER largely come from MAPK1.
Among positions with selective substitutions there doesn't appear to be any particular pattern.
Outliers are somewhat more likely to occur in helices compared to other subtype. -->

## D - Aspartate

Summary:

* D1 -
* D2 -
* DP - Fully permissive, very accessible
* DO - Outliers

D4 and D5 are both in the overall permissive subcluster, which seems to be subdivided into "anything but proline" and very permissive/improving subclusters.

### D1

<!-- Often appear to have a specific role (see examples), although not always.
In other cases it could simply be we cannot easily see the role.

Examples:

* HA - 93 is directly adjacent to an arginine (likely ionic interaction?)
* CBS - 221 is directly adjacent to an arginine. 49, 309 & 198 are all grouped around a phenylalanine aromatic group; possibly they induce a dipole? 388 is near an arginine
* TEM1 - 189, 208 neighbour a K ion
* HSP90 - All four (40, 79, 113, 143) all border an apparent ATP binding pocket, as do 3 D3 positions
* PAB1 - 184 possibly bonding to lysine 180, 138 has a similar arrangement with lysine 140 -->

### D2

<!-- Tolerates any substitution other than proline.

Examples:

* APH3II - 118 & 238 are in helices, 23 & 52 at the start and end of beta strand turns
* PAB1 - 144 is in a helix, 136 & 160 are in coils around this helix
* TPMT - 31 & 162 are both in helices -->

### DP

### D0

<!-- A single position in ADRB2 (29), which appears to prefer being a serine.
This is potentially because there's a repulsive charge interaction with a nearby glutamate, although that wouldn't explain why serine specifically and not any other residue. -->

## E - Glutamate

Summary:

* E1 -
* E2 -
* E3 -
* EP -
* EO - Outliers

### E1

<!-- As with D1, E1 positions often seem to have a recognisable role utilising the charge properties of glutamate, although in many cases it is less obvious too.

Examples:

* ADRB2 - Not clear what role they play but 188 and 306 are both near the extracellular binding site
* amiE - 105 & 108 both extend towards a similar pair of lysine and arginines. 248, 249 & 250 are arrayed at the interface between subunits, next to the equivalent residues on the other subunit, creating a very negative pocket. Are many other positions where role isn't clear
* APH3II - 262 is in the cleft where kanamycin is bound (APH3II is a kanamycin resistance protein). 182 is at the negative end of an alpha helix, possible stabilising its dipole (Sali et al. 1988).
* UBE4B - 1083 and 1084 are on a very polar exterior helix face (2 other E, 2 K and 1 R) -->

### E2

### E3

<!-- Tolerates any substitution apart from proline.

Examples:

* ADRB2 - 180 & 225 are both at the end of helices (quite common in general for E though) and 237 is in an intracellular loop, although this part of the protein seems to vary between structures.
* CALM1 - All three positions are in helices
* CBS - 399 is on the end of a helix and 436 in the middle of a short helix
* MAPK1 - 33 is on the end of a beta sheet, 250 near the end of a short helix and 186 in a coil
* Ras - 91 & 98 are in a single helix, 162 in a helix, 76 in the coil at the end of a sheet and 62 & 63 in a fairly tight turn in a beta sheet structure
* TPMT - 205 is in a curved beta sheet (saddle point like curve) and 203 in a small alpha helix -->

### EP

### EO

## F - Phenylalanine

Summary:

* F1 -
* F2 -
* FP - Permissive positions
* FO - Outliers

F1 and F2 are very intermixed in the true dendrogram, although F2's mean profile is more similar to F3's

### F1

### F2

<!-- Seems to utilise specific aromatic properties as well as hydrophobicity.
They occur relatively frequently in groups.

Examples:

* BRCA1 - 43 surrounded by hydrophobic residues and filling a ring shaped space
* CALM1 - 17 & 69 are potentially pi-stacking
* CBS - 332, 385 & 396 are positioned in a cleft together, along with several other aromatic residues; more potential pi-stacking? 197 & 310 are also adjacent to each other
* TP53 - 270 & 113 are adjacent and face one another. 134 is near them but in the middle of an apparently polar pocket? -->

### FP

### FO

<!-- Generally have a mixture of positive substitutions (again largely in MAPK1) and permissive.
Two positions have mostly permissive substitutions with a few somewhat selected against. -->

## G - Glycine

Summary:

* G1 -
* G2 -
* G3 -
* G4 -
* G5 -
* GP - Permissive positions
* GO - Outliers

All glycine residues appear enriched in turn and bend secondary structures

### G1

### G2

### G3

### G4

### G5

### GP

### GO

<!-- Outliers split into 3 broad types: positive ER (often in MAPK1 again, but not exclusively); permissive positions; and positions where a few substitutions are strongly selected against.
Four positions have a single very strongly deleterious substitution. -->

## H - Histidine

Summary:

* H1 -
* HP - Permissive positions
* HO - Strong improvement when replaced with methionine or alanine, otherwise mix of weak  selection/improvement

### H1

### HP

### HO

## I - Isoleucine

Summary:

* I1 -
* I2 -
* I3 -
* IP - Permissive positions
* IO - Outliers

### I1

### I2

### I3

### IP

### IO

## K - Lysine

Summary:

* K1 -
* K2 -
* K3 -
* KP - Permissive positions
* KO - Outliers

## L - Leucine

Summary:

* L1 -
* L2 -
* L3 -
* L4 -
* L5 -
* L6 -
* L7 -
* L8 -
* LP - Permissive positions
* LO - Outliers

### L1

### L2

### L3

### L4

### L5

### L6

### L7

### L8

### LP

### LO

## M - Methionine

Summary:

* M1 - Tolerates leucine, isoleucine, valine and threonine (longer aliphatic chains?). Strongly intolerant to negative charge and proline
* M2 - Strong selection against proline, weak selection against negative charge
* MP - Permissive positions
* MO - Outliers

### M1

<!-- Examples:

* HSP90 - 84 is in a binding pocket for ATP, 116 fills a hydrophobic cleft, 16 looks to interact with a phenylalanine ring, 589 is near to a tryptophan residue (interestingly the 589 position is classified M1 from the profile in Hietpas et al. 2011 and M3 when from Jiang et al. 2013).
* MAPK1 - 221 & 293 appear to be filling hydrophobic space
* NP - series of positions in a helix bundle
* TP53 - 243 faces towards bound DNA, 133 is in a space near two aromatics, 160 seems to fill a large hole
* Src - both fill holes near aromatics -->

### M2

<!-- Examples:

* Src - 317 in a hydrophobic pocket near two aromatics
* Ras - 72 in a pocket near a phenylalanine
* NP - Clustered in a different helix bundle to M1 positions, but not obvious what the difference is other than potentially more aromatics
* HA - 460 in a fairly tight looking pocket near a tyrosine.
* CBS - 464 is in a tight pocket, 458 in a similar pocket but with two aromatics, 337 is near two aromatics in an exposed cleft
* CALM1 - 72 & 73 are in a helix helix interface adjacent to two M1 positions

The physical difference between M1 and M2 is not completely clear.
They have similar distributions of distance to the nearest aromatic.
To some extent M1 positions appear to physically group together in some proteins and M2 positions group separately, potentially due to physical characteristics or the importance of the two domains. -->

### MP

### MO

Outliers are generally permissive with a few exceptions, for example ADRB2 96 is very selective against asparagine, Src 286 prefers glycine & UBE2I 1 is somewhat selective.

## N - Asparagine

Summary:

* N1 - Moderately selective against everything, weakest against serine and threonine
* N2 - Nonselective
* N3 - Anything but proline
* N4 - Selective against positive charge and proline

### NO

Outliers are all positive ER positions from MAPK1, with a particular preference for hydrophobic residues.
This suggests addition of hydrophobic residues here creates a GoF phenotype in MAPK1, but no other asparagine's have this property.

## P - Proline

Summary:

* P1 - Moderately selective against everything
* P2 - Selective against aromatics
* P3 - Tolerates alanine, leucine, glutamine, serine and threonine, strong selection against negative charge and methionine
* P4 - Selects against negative charge and asparagine
* P5 - Selects against glutamine and (less so) threonine
* PO - Outliers, on average weakly improved by substitution

### PO

<!-- Generally nonselective outliers as well as some positive ER MAPK1 positions -->

## Q - Glutamine

Summary:

* Q1 - Moderately selective for basic and polar
* Q2 - Anything but proline
* Q3 - Nonselective
* Q4 - Weak improvement with threonine/serine and hydrophobic substitutions
* QO - Outliers, weak selection against negative charge, strong selection against lysine and histidine

## R - Arginine

Summary:

* R1 - Most selective, only tolerates lysine
* R2 - Strongly selective against proline, weakly selective against negative charge
* R3 - Nonselective
* R4 - Intolerant to negative charge, proline, and aromatic residues
* RO - Outliers, nonselective

### RO

<!-- Four outliers; 2 nonselective and 2 in ADRB2 with some selection, in particular 228 which strongly selects against lysine and proline. -->

## S - Serine

Summary:

* S1 - Strongly selective, tolerates threonine and somewhat alanine and glycine.
* S2 - Nonselective
* S3 - Selective against negative charge and aromatic
* S4 - Anything by proline
* S5 - Selection against positive charge and proline, and less so negative charge
* S6 - Anything but alanine
* S7 - Selective against tryptophan and proline
* SO - Outliers

## T - Threonine

Summary:

* T1 - Tolerates serine only, strongest selection against aromatic and proline
* T2 - Anything but proline
* T3 - Nonselective
* T4 - Selects against negative charge and aromatics
* T5 - Selective against negative charge, proline and glycine
* T6 - Mostly nonselective, weak selection against tyrosine, tryptophan and asparagine
* T7 - Intolerant to charge

## V - Valine

Summary:

* V1 - Tolerates hydrophobic residues and threonine
* V2 - Similar to V1 but doesn't tolerate threonine and weaker tolerance of non-valine hydrophobic residues
* V3 - Intolerant of positive charge, proline and polar residues with NH2 groups (Q and N, which possibly tend to have positive partial charges on N?
* V4 - Nonselective
* V5 - Weaker selection against aromatics and lysine
* V6 - Anything but proline
* VO - Outliers

### VO

<!-- Four MAPK1 outliers with somewhat positive ER scores and one nonselective CP position. -->

## W - Tryptophan

Summary:

* W1 - Tolerates phenylalanine and tyrosine, and is less intolerant other large residues.
* W2 - Tolerates cysteine, glycine, leucine, arginine and serine (no clear pattern?). Less intolerant of other aromatics
* WO - Outliers

### WO

<!-- Two positions in Src, one of which has most substitutions strongly favoured particularly for hydrophobic residues and one position that is nonselective apart from preferring to be lysine. -->

### W1

<!-- Tolerates phenylalanine and tyrosine as well as being less intolerant of methionine, leucine, isoleucine and histidine (all also reasonably large)

Examples:

* ADRB2 - 4 positions face into the membrane (32, 105, 158 & 173), 99 & 109 look to potentially be interacting via pi orbitals, 313 fills a large-ish space in the helix bundle and 286 is in the bundle at the base of a ligand binding site.
* amiE - 209 is buried deep in the core, potentially interacting with a histidine nearby (although it is positioned perpendicularly)
* CBS - 7 positions, generally in positions where they fill a large space. 408 and 410 are part of a tryptophan triplet, with a W2 position in between, which seems to fill a big space and potentially create pi orbital interactions.
* CXCR4 - 161 is on the exterior of a helix bundle. 252 faces into the helix bundle into a region with many aromatics. 195 is at the edge of the bundle in the interface with a second subunit (in the homodimer structure), close to two other aromatics.
* HA - 196 & 250 are in the internal pocket of the binding domain, various other positions are buried in the binding domain with several buried in the transmembrane helix domains.
* MAPK1 - Two positions in generically hydrophobic positions.
* TEM1 - 226 & 287 are in a pair facing each other with their sides on the protein surface, 163 is directly facing the surface and 208 is buried near a disulphide bond, several charged residues, a leucine and an isoleucine. -->

### W2

<!-- Tolerates aromatics and some other (rather discordant) substitutions.

Examples:

* amiE - 5 positions, generally buried and surrounded by a range of different residue types.
* APH3II - 69 faces into the protein towards various hydrophobic residues.
* CCR5 - 5 positions face into the core of a helix bundle, with a number of other aromatic residues. 190 is on the outside of the bundle, facing into the membrane.
* TPMT - 29 and 33 face the binding pocket of SAH (in this structure), 230 conjugates pi orbitals with a histidine, 78 is on the protein surface
* CXCR4 - 94 is near various other aromatics, 102 is in the helix bundle and 283 faces out from it (into the membrane) -->

## Y - Tyrosine

Summary:

* Y1 - Generally selective, tolerates phenylalanine and histidine
* Y2 - Intolerant of charged and polar residues, most tolerates phenylalanine, histidine and hydrophobic residues
* Y3 - anything but proline
* YO - Outliers

Subtypes are fairly intermixed in the original hierarchical clustering dendrogram, and there is not a large, obvious difference when looking at Y1 and Y2 positions.

### YO

<!-- Most positions are nonselective or only weakly selective, with MAPK1 positions having weakly positive ER scores.
ADRB2 141 strongly selects against aspartate.
CXCR4 121 is somewhat selective against polar residues. -->

### Y1

<!-- Tolerates phenylalanine and histidine, the other two amino acids with a single aromatic ring.
The intolerance of tryptophan potentially means they are in spaces where it is too bulky.
The requirement for an aromatic ring suggest pi interactions or similar.

Examples:

* amiE - Many positions in a cluster, large number of aromatics maybe means tryptophan is too bulky?
* CBS - Again a group appears in a sort of cluster (4 along a line), with another two independent: 233 is part of a haem binding pocket and 163 appears alone in a reasonably large pocket (potentially needing the size of the aromatic ring)
* CCR5 - A number of positions on the inner faces of the helix bundle (not adjacent to each other however), then two pairs of nearby residues (184/187 & 14/15) plus 10 is sulphated and near to an N-Acetyl-2-Deoxy-2-Amino-Galactose ligand (at least in the cryoEM structure), so potentially has an active role.
* CP - 3 positions in the beta sheet. 86 facing externally and the other two buried. As this is a bacteriophage coat protein it is likely that 86 actually faces other complex subunits.
* HA - many positions, somewhat clustered? 502 & 507 seem to be potentially stacking?
* MAPK1 - not a clear pattern, some internal in large-ish pockets, a few on the surface.
* TP53 - 327 externally facing, other positions clumped in a lobe but not immediately next to one another. 126 is near two other aromatics, potentially interacting, others appear to be in simple hydrophobic pockets.
* UBE2I - 137 & 144 are adjacent, potentially interacting. 68 is near a phenylalanine, with another potential interaction. -->

### Y2

<!-- Tolerant of all aromatic and hydrophobic substitutions

Examples:

* ADRB2 - 209 & 132 face outwards in the transmembrane domain
* amiE - 255 is buried on the face of a beta sheet, but doesn't seem to be hugely hydrophobic
* APH3II - both are internal, in hydrophobic pockets but also are predicted to make h-bonds to the backbone
* CBS - 484 is buried but also has a predicted h-bond
* HSP90 - 125 is internal, 47 externally facing, although with the aliphatic ring surrounded
* PTEN - 6 positions are near the surface but not really accessible. Generally seem to be in reasonably hydrophobic positions, but not entirely.
* Src - 329 is buried although not in a hydrophobic region (makes polar contacts with 2 glutamine residues), 439 & 482 are on the surface.
* TPMT - 131 is somewhat exposed and contacts a ligand in the crystal structure. 180 & 240 are on a beta sheet, in apparently hydrophobic conditions. -->
