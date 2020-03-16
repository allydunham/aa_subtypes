# Subtype Descriptions

Analysis of the subtypes produced by hierarchical clustering on PC2:20, using cosine distance, dynamic tree cutting and with permissive positions (all |ER| < 0.4) excluded into their own subtype before the main clustering step.
The deepSplit parameter for dynamic tree cutting is selected for each AA qualitatively, looking at the clusterings produced by both 0 and 1 and the mean rank correlation of position profiles within clusters.

## A - Alanine

deepSplit = 0

Summary:

* A1 - Tolerates smaller residues
* A2 - Selects against (most) polar, charged & aromatic
* A3 - (Weak) anything but proline
* A4 - Selects against aromatic and larger hydrophobic
* A5 - Strong selection, tolerates glycine & serine
* AP - Permissive positions
* AO - Outliers

### A1

These positions tolerate smaller residues, generally being buried in the protein and tending to be near other hydrophobic residues.

Examples:

* ADRB2 - Several positions in the core of the transmembrane helices (76, 78 & 118), two residues in the cytoplasm side (13 & 19), plus various others.
* amiE - cluster of residues in the core of a fold and some between two alpha helices (35, 39, 90 & 94)
* TP53 - 276 faces towards the DNA, 159 & 161 face into the protein core from a beta sheet and 189 is unclear
* PTEN - 72 is in a fairly tight turn
* Src - 433 is in a turn between two short alpha helices, 374 is in the core of a helix bundle

### A2

These positions select for hydrophobic residues, as well as tolerating serine, threonine and cysteine, all of which are small and likely don't disrupt the structure with polarity (potentially h-bonding to the backbone to resolve dipoles).
They are also usually buried.

Examples:

* APH3II - 4 positions in the core of a helix bundle
* CALM1 - Most positions are core facing, with one (104) facing outside
* CXCR4 - generally found within the helix bundles, difference with A1 is not clear
* Ras - 18 & 146 face each other, several other positions on the interior of a fold,

### A3

Spread throughout proteins these positions are only selective against proline.
They are less conserved than other subtypes.

Examples:

* ADRB2 - All four positions are at the start and end or the transmembrane helices
* APH3II - many positions, many at start/end of secondary structures with some in the turns between them
* GAL4 - Both positions are at the start of helices
* TEM1 - 4 positions in helices, one in a loop in between a helix and sheet

### A4

These positions select against aromatics and larger hydrophobic residues.
Often surface accessible they likely interact badly with the solvent/prefer to be buried and disrupt the structure when a larger hydrophobic group is added.

Examples:

* CBS - Many positions, mixture between surface residues and buried ones, with some seemingly cramped
* TPK1 - 4 & 84 face outwards, 43 & 129 are internal and could potentially clash if they were larger
* UBE2I - 2 positions, both external
* NP - Many positions, almost all exposed

### A5

Highly selective, these positions only tolerate glycine and serine.
They are often buried, although not as strongly as A1/2

Examples:

* NP - many positions, generally buried. Most don't look unusually cramped
* TPK1 - Buried, but doesn't appear cramped
* UBE2I - 152 is buried and near to two large aromatics, 44 is buried but not cramped, 129 is on the surface
* Src - 314 is buried with a phenylalanine nearby

### AP

Permissive positions, very depleted in beta strands compared to other subtypes.
They are very surface accessible

Examples:

* amiE - all on the protein surface
* CXCR4 - three positions on the cytoplasmic end of two alpha helices and the loop between them
* TEM1 - many positions all surface accessible
* HSP90 - five positions, all near each other on the surface

### AO

Two outlier positions, TP53 83/84, both of which are generally permissive (presumably one or two substitution ERs just over the threshold of being a permissive position).

## C - Cysteine

deepSplit = 0

Summary:

* C1 - Highly selective, cysteine only functions such as disulphide bonds
* C2 - Selective against hydrophilic amino acids, much more buried, includes cysteine/aromatic interactions.
* CP - Permissive positions
* CO - Outliers

### C1

FoldX and chemical environment results suggest disulphide bonds generally go here, although not exclusively (24/73 positions vs 3/37 in C2).
The subtype likely also contains other positions where only C functions, such as ion conjugation, and potentially includes similar aromatic interactions to C2 where the interaction is highly important.
The subtype is least selected against aromatics - maybe because the aromatic/sulphur interaction partially makes up for the disulphide bond, or also fulfils some of the other functions.

Examples without disulphide bond:

* NP - 44 & 279 are adjacent to each other with a tyrosine ring between them - potentially some interaction through the aromatic? 223 is near a tyrosine, potentially having an aromatic interaction.
* TP53 - 176, 238 & 242 all appear to be conjugating a Z ion (Pace and Weerapana, 2014). Others unclear
* CCR5 - 224 is on the receptor surface, potentially with an active role
* ADRB2 - Several membrane facing cysteine residues on the transmembrane helices?
* GAL4 - 4 residues (14, 21, 28, 38) conjugating two zinc ions
* APH3II - 131 appears to be interacting with a tryptophan residue
* BRCA1 - 39 & 64 and 27 & 47 each conjugate a zinc ion
* CBS - 52 is bonded to the Fe of a haem group in the crystal structure (which is a known binding in the native protein)
* SUMO1 - 52 is directed at a phenylalanine

### C2

This subtypes appears to use cysteine's hydrophobicity as well as the sulphur/aromatic interaction that allows cysteine to create stabilising interactions with aromatics (Orabi & English 2016)
The characteristic aromatic interaction pattern occurs in the subtypes mean profile and somewhat in the chemical environment, although the enrichment for aromatic neighbours is weaker than when it was an exclusive group (likely because it is swamped by residues without the interaction).
In some sense this subtype could be manually split into two, although that is generally not identified algorithmically (apart from a single lucky run of kmeans)
These positions tend to be more strongly buried than C1 positions, although both generally are.

The identified aromatic interactions are:

* Src 501 - In a helix bundle surrounded by hydrophobic residues, with an aromatic nearby but not adjacent
* TEM1 75 - Also in a bundle and near an aromatic ring, but also one half of a disulphide bond (the other half is not in C2)
* ADRB2 327 - Semi-buried residue near hydrophobic residues and an aromatic
* TPMT 133 - Residue is on the exterior but the side chain faces inwards into a hydrophobic pocket
* UBE2I 43/75 - Near each other in 3D space, both side chains face into the protein and there are hydrophobic and aromatic residues nearby, although not directly conjugating like in some examples.
* APH3II 209 - near to a phenylalanine
* APH3II 192 - In a cleft with a NA ion and an aromatic drug in the crystal structure
* CCR5 213 - directly conjugates a tyrosine ring
* TPK1 104 - Directly conjugates a phenylalanine ring
* TPK1 88 - Near but not directly adjacent to a tyrosine ring

### CP

Three permissive positions, which are weakly conserved and with low structural impact.

### CO

Outliers

## D - Aspartate

deepSplit = 1

Summary:

* D1 - Selects for negative charge
* D2 - Selects for negative charge or polarity
* D3 - Anything but proline
* DP - Fully permissive, very accessible
* DO - Outliers, often improves

### D1

Often appear to have a specific role (see examples), although not always.
In other cases it could simply be we cannot easily see the role.

Examples:

* HA - 93 is directly adjacent to an arginine (likely ionic interaction?)
* CBS - 221 is directly adjacent to an arginine. 49, 309 & 198 are all grouped around a phenylalanine aromatic group; possibly they induce a dipole? 388 is near an arginine
* TEM1 - 212, 231 neighbour a K ion
* HSP90 - 40 & 143 all border an apparent ATP binding pocket, as do 3 DP positions and a D2 position

### D2

Generally tolerates some polarity as well as negative charge.

Examples:

* HSP90 - 79 is near the ATP binding pocket

### D3

Tolerates any substitution other than proline.

Examples:

* APH3II - 118 & 238 are in helices, 23 & 52 at the start and end of beta strand turns
* PAB1 - 144 is in a helix, 136 & 160 are in coils around this helix
* TPMT - 31 & 162 are both in helices

### DP

### D0

## E - Glutamate

deepSplit = 0

Summary:

* E1 - (Weak) selection for negative charge
* E2 - Anything but proline
* E3 - Selects against aromatics and larger hydrophobic residues
* EP - Permissive positions
* EO - Outliers

### E1

As with D1, E1 positions often seem to have a recognisable role utilising the charge properties of glutamate, although in many cases it is less obvious too.

Examples:

* ADRB2 - Not clear what role they play but 188 and 306 are both near the extracellular binding site
* amiE - 105 & 108 both extend towards a similar pair of lysine and arginines. 248, 249 & 250 are arrayed at the interface between subunits, next to the equivalent residues on the other subunit, creating a very negative pocket. Are many other positions where role isn't clear
* APH3II - 262 is in the cleft where kanamycin is bound (APH3II is a kanamycin resistance protein). .
* UBE4B - 1083 and 1084 are on a very polar exterior helix face (2 other E, 2 K and 1 R)

### E2

Tolerates any substitution apart from proline.

Examples:

* Ras - 91 & 98 are in a single helix, 162 in a helix and 76 in the coil at the end of a sheet
* TPMT - 203 in a small alpha helix

### E3

Examples:

* APH3II - 182 is at the negative end of an alpha helix, possible stabilising its dipole (Sali et al. 1988)

### EP

### EO

## F - Phenylalanine

deepSplit = 1

Summary:

* F1 - Strong selection against charge & polarity
* F2 - Strong selection against non-aromatics
* F3 - Tolerates Cysteine, isoleucine, leucine, serine valine and tyrosine?
* FP - Permissive positions
* FO - Outliers

F1 and F2 are very intermixed in the true dendrogram, although F2's mean profile is more similar to F3's

### F1

### F2

Seems to utilise specific aromatic properties as well as hydrophobicity.
They occur relatively frequently in groups.

Examples:

* BRCA1 - 43 surrounded by hydrophobic residues and filling a ring shaped space
* CBS - 332, 385 & 396 are positioned in a cleft together, along with several other aromatic residues; more potential pi-stacking? 197 & 310 are also adjacent to each other
* TP53 - 270 & 113 are adjacent and face one another. 134 is near them but in the middle of an apparently polar pocket?

### F3

Examples:

* CALM1 - 17 & 69 are potentially pi-stacking

### FP

### FO

## G - Glycine

deepSplit = 1

Summary:

* G1 - Strong selection against all
* G2 - Selects for smaller hydrophobic residues
* G3 - Selects against larger hydrophobic/aromatics
* G4 - Anything but proline or isoleucine
* G5 - Selects against polarity, but sometimes only weakly
* G6 - Largely permissive, on average
* G7 - Strong selection against aspartate and lysine, with weaker selection against other charge/polarity
* GP - Permissive positions
* GO - Outliers

All glycine residues appear enriched in turn and bend secondary structures

### G1

### G2

### G3

### G4

### G5

### G6

### G7

### GP

### GO

## H - Histidine

deepSplit = 0

Summary:

* H1 - Overall strong selection, most tolerant of tyrosine, glutamine and asparagine
* HP - Permissive positions
* HO - Strong improvement when replaced with methionine or alanine, otherwise mix of weak positive and negative effects

### H1

### HP

### HO

## I - Isoleucine

deepSplit = 0

Summary:

* I1 - Selects for hydrophobicity
* I2 - Selects for larger hydrophobic residues
* I3 - Strong selection against proline, weaker selection against leucine and lysine (?)
* IP - Permissive positions
* IO - Outliers

### I1

### I2

### I3

### IP

### IO

## K - Lysine

deepSplit = 1

Summary:

* K1 - Weak selection for polarity
* K2 - Anything but proline
* K3 - Selects against threonine, glutamine & glutamate
* KP - Permissive positions
* KO - Outliers

### K1

### K2

### K3

### KP

### KO

## L - Leucine

deepSplit = 0

Summary:

* L1 - Selects against charge and polarity
* L2 - Only tolerates methionine and isoleucine
* L3 - Selects for hydrophobic residues
* L4 - Strong selection against lysine, weaker against other charge
* L5 - Selection against negative charge (much weaker against positive)
* L6 - Anything but proline
* L7 - Weak selection against isoleucine, proline, tyrosine & tryptophan
* L8 - Strong selection against arginine
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

deepSplit = 0

Summary:

* M1 - Tolerates leucine, isoleucine, valine and threonine (longer aliphatic chains?). Strongly intolerant to negative charge and proline
* M2 - Strong selection against proline, weak selection against negative charge
* MP - Permissive positions
* MO - Outliers

### M1

Examples:

* HSP90 - 16 looks to interact with a phenylalanine ring, 589 is near to a tryptophan residue (interestingly the 589 position is classified M1 from the profile in Hietpas et al. 2011 and M2 when from Jiang et al. 2013).
* TP53 - 243 faces towards bound DNA, 133 is in a space near two aromatics, 160 seems to fill a large hole
* NP - Many positions throughout two helix bundles, often group
* Ras - 72 in a pocket near a phenylalanine
* Src - 286 on surface near an aromatic, 484 and 317 in a hydrophobic pocket near two aromatics
* CBS - 464 is in a tight pocket, 337 is near two aromatics in an exposed cleft
* CALM1 - 37, 52, 72, 73 & 77 are adjacent to each other in a helix helix interface
* HA - 460 in a fairly tight looking pocket near a tyrosine.

### M2

Examples:

* HSP90 84 is in a binding pocket for ATP, 116 fills a hydrophobic cleft,

### MP

### MO

Outliers are generally permissive with a few exceptions, for example ADRB2 96 is very selective against asparagine, Src 286 prefers glycine & UBE2I 1 is somewhat selective.

## N - Asparagine

deepSplit = 0

Summary:

* N1 - Moderately selective against everything, weakest against serine and threonine
* N2 - Anything but proline
* NP - Permissive positions
* NO - Outliers

### N1

### N2

### NP

### NO

Outliers are all positive ER positions from MAPK1, with a particular preference for hydrophobic residues.
This suggests addition of hydrophobic residues here creates a GoF phenotype in MAPK1, but no other asparagine's have this property.

## P - Proline

deepSplit = 1

Summary:

* P1 - Selects moderately against everything
* P2 - Tolerates serine, threonine, glutamine, leucine and alanine
* P3 - Strong selection against aromatic, weaker against other large hydrophobic residues
* P4 - Selects against charge and polarity
* PP - Permissive positions
* PO - Outliers, on average weakly improved by substitution

### P1

### P2

### P3

### P4

### PP

### PO

## Q - Glutamine

deepSplit = 1

Summary:

* Q1 - Generally selective, most tolerant to polarity or charge, strongest against aromatics
* Q2 - Anything but proline
* Q3 - Weak selection against hydrophobic residues
* Q4 - Strong selection against negative charge
* Q5 - Improves on average, apart from lysine
* QP - Permissive positions

### Q1

### Q2

### Q3

### Q4

### Q5

### QP

## R - Arginine

deepSplit = 1

Summary:

* R1 - intolerant of everything apart from lysine
* R2 - Anything but proline
* R3 - Selects against negative charge and proline
* R4 - Weak selection
* R5 - Selects against negative charge, aromatics and proline
* RP - Nonselective
* RO - Outliers, nonselective

### R1

### R2

### R3

### R4

### R5

### RP

### RO

## S - Serine

deepSplit = 0

Summary:

* S1 - Only tolerates threonine
* S2 - Anything but proline
* S3 - Selects against negative charge
* S4 - Weak selection against polar
* S5 - Anything but alanine
* SP - Permissive positions
* SO - Outliers

### S1

### S2

### S3

### S4

### S5

### SP

### SO

## T - Threonine

deepSplit = 1

Summary:

* T1 - Only tolerates serine
* T2 - Anything but proline
* T3 - Selects against aromatics and (less so) negative charge
* T4 - Selects against negative charge, glycine and proline
* T5 - Weak selection against cysteine only
* T6 - Strong selection against asparagine, weaker selection against other positive charge, negative charge, aromatics and cysteine
* TP - Permissive positions
* TO - Outliers0

### T1

### T2

### T3

### T4

### T5

### T6

### TP

### TO

## V - Valine

deepSplit = 0

Summary:

* V1 - Selects for hydrophobic residues and threonine
* V2 - Only tolerates isoleucine
* V3 - Selects against charged (apart from glutamate) and some polar, particularly histidine
* V4 - Similar to V3 but selects against glutamate instead of aspartate and particularly disfavours glutamine
* V5 - Anything but proline
* VP - Permissive positions
* VO - Outliers

### V1

### V2

### V3

### V4

### V5

### VP

### VO

## W - Tryptophan

deepSplit = 0

Summary:

* W1 - Tolerates phenylalanine and tyrosine, and is less intolerant other large residues.
* W2 - Tolerates cysteine, glycine, leucine, arginine and serine (no clear pattern?). Less intolerant of other aromatics
* WP - Permissive positions
* WO - Outliers

### W1

Tolerates phenylalanine and tyrosine as well as being less intolerant of methionine, leucine, isoleucine and histidine (all also reasonably large)

Examples:

* ADRB2 - 4 positions face into the membrane (32, 105, 158 & 173), 99 & 109 look to potentially be interacting via pi orbitals, 313 fills a large-ish space in the helix bundle and 286 is in the bundle at the base of a ligand binding site.
* amiE - 209 is buried deep in the core, potentially interacting with a histidine nearby (although it is positioned perpendicularly)
* CBS - 43 & 208 are adjacent to each other, 390 is near a phenylalanine and 323 fills a large space
* CXCR4 - 161 is on the exterior of a helix bundle. 252 faces into the helix bundle into a region with many aromatics. 195 is at the edge of the bundle in the interface with a second subunit (in the homodimer structure), close to two other aromatics.
* HA - 196 & 250 are in the internal pocket of the binding domain, various other positions are buried in the binding domain, 437 is buried in the transmembrane helix domain and 359 & 366 are near to each (but not adjacent) at the domain on the other end of the helix bundle.
* TEM1 - 227 & 286 are in a pair facing each other with their sides on the protein surface and 208 is buried near a disulphide bond, several charged residues, a leucine and an isoleucine.

### W2

Tolerates cysteine, glycine, leucine, arginine and serine (no clear pattern?). Less intolerant of other aromatics.

Examples:

* amiE - 5 positions, generally buried and surrounded by a range of different residue types.
* APH3II - 69 faces into the protein towards various hydrophobic residues.
* CCR5 - 4 positions face into the core of a helix bundle, with a number of other aromatic residues. 190 is on the outside of the bundle, facing into the membrane.
* TPMT - 29 and 33 face the binding pocket of SAH (in this structure), 230 conjugates pi orbitals with a histidine, 78 is on the protein surface
* CXCR4 - 94 is near various other aromatics, 102 is in the helix bundle and 283 faces out from it (into the membrane)

### WP

### WO

## Y - Tyrosine

deepSplit = 1

Summary:

* Y1 - Tolerates aromatics, particularly intolerant of charge (other than aspartate) and glycine (?)
* Y2 - Intolerant of charged and polar residues, most tolerates phenylalanine, histidine and hydrophobic residues
* Y3 - Generally selective, somewhat tolerates phenylalanine and histidine
* Y4 - Anything but proline
* YP - Permissive positions

Subtypes are fairly intermixed in the original hierarchical clustering dendrogram, and there is not a large, obvious difference when looking at Y1 and Y2 positions.

### Y1

Tolerates phenylalanine and histidine, the other two amino acids with a single aromatic ring.
The intolerance of tryptophan potentially means they are in spaces where it is too bulky.
The requirement for an aromatic ring suggest pi interactions or similar.

Examples:

* amiE - Many positions in a cluster, large number of aromatics maybe means tryptophan is too bulky?
* CBS - 233 is part of a haem binding pocket
* CCR5 - A number of positions on the inner faces of the helix bundle (not adjacent to each other however), then two pairs of nearby residues (184/187 & 14/15) plus 10 is sulphated and near to an N-Acetyl-2-Deoxy-2-Amino-Galactose ligand (at least in the cryoEM structure), so potentially has an active role.
* HA - 502 & 507 seem to be potentially stacking?
* TP53 - 327 externally facing, 126 is near two other aromatics, potentially interacting
* UBE2I - 134 is adjacent to a Y3 position (144), potentially interacting

### Y2

Tolerant of all aromatic and hydrophobic substitutions

Examples:

* ADRB2 - 209 & 132 face outwards in the transmembrane domain
* amiE - 255 is buried on the face of a beta sheet, but doesn't seem to be hugely hydrophobic
* APH3II - 218 & 244 are internal, in hydrophobic pockets but also are predicted to make h-bonds to the backbone. 22 is near the surface although not truly accessible, and again makes h-bonds
* HSP90 - 125 & 24 are internal and near to one another
* PTEN - 6 positions are near the surface but not really accessible. Generally seem to be in reasonably hydrophobic positions, but not entirely.
* Src - 482 is on the surface.
* TPMT - 131 is somewhat exposed and contacts a ligand in the crystal structure. 180 & 240 are on a beta sheet, in apparently hydrophobic conditions.

### Y3

Examples:

* CBS - Again a group appears in a sort of cluster (223, 308, 301) along a line, with a Y1 position (381) as well. 163 appears alone in a reasonably large pocket (potentially needing the size of the aromatic ring) and 484 fills a pocket near a pair of aromatics (phenylalanine/histidine)
* HA - many positions, somewhat clustered at either in the binding domain lobe.

### Y4

### YO
