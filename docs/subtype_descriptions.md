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

Nine outlier positions, most of which have a different selection of weakly selective substitutions

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

Generally tolerates some polarity as well as negative charge. Often surface accessible

Examples:

* HSP90 - 79 is near the ATP binding pocket
* amiE - Many positions, generally near a polar residue. E.g. 167 faces a glutamine and 239 is near a glutamine and arginine
* CALM1 - 94 is near a serine, mostly surface accessible
* TEM1 - three surface accessible positions
* TPMT - two surface accessible positions

### D3

Tolerates any substitution other than proline.

Examples:

* APH3II - 118 & 238 are in helices, 23 & 52 at the start and end of beta strand turns
* PAB1 - 144 is in a helix, 136 & 160 are in coils around this helix
* TPMT - 31 & 162 are both in helices

### DP

Reasonably large number of positions (51), mostly surface accessible.

Examples:

* HSP90 - in one lobe of the protein, mostly in loops between secondary structures and all accessible
* PAB1 - Again clustered together, with 3 of 4 in loops and 1 in a helix. All accessible
* UBI - All aspartate positions in UBI are found to be permissive, which is likely incorrect

### D0

Ten positions, generally largely permissive

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
* UBI - 34 is at the end of an alpha helix
* TEM1 - Two positions in an alpha helix (270 & 277) and 210 is at the start of another, 56 is in a beta sheet and 61 & 87 are in loops
* GAL4 - three positions (56, 62 & 65) are in an alpha helix involved in DNA binding

### E3

Selects against aromatics and larger hydrophobic residues

Examples:

* APH3II - 182 is at the negative end of an alpha helix, possible stabilising its dipole (Sali et al. 1988)
* ADRB2 - various positions in on the intracellular domain, exposed to the cytoplasm. Some are potentially involved with binding to the coupled G-protein.
* NP - four positions all surface accessible.
* TEM1 - two surface accessible residues (119 & 195) and one semi buried (35) in a cramped position, likely interacting with a nearby arginine but also surrounded by hydrophobic residues.
* TPMT - 205 is surface accessible in a helix, 237 is in a cleft with the side chain partially exposed

### EP

Permissive positions, generally very surface accessible

Examples:

* HSP90 - many very accessible positions in one lobe of the protein. Mixture of loops and secondary structure elements.
* TPMT - Two positions, both accessible
* TP53 - Two positions, both accessible and far from the DNA binding site

### EO

Outlier positions, generally largely permissive, but a few positions that don't tolerate aspartate substitutions, all in ADRB2.
One is on the extracellular side of ADRB2 (17), although not particularly near the binding site.
A second is on the other end of the same helix on the intracellular side, although again not near the G-protein binding region.
THe other two are not in the structure, but likely in the G-protein binding domain.

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

Strong selection against charge and polarity, generally buried.

Examples:

* ADRB2 - 290 & 193 are in the extracellular binding pocket. A number of positions are in the transmembrane helices, most facing out into the membrane as well as some facing into the core of the helix bundle. 10 is on the intracellular domain, likely in the G-protein binding domain.
* amiE - many positions, all buried and in largely in hydrophobic pockets
* CBS - 443 is partially exposed, 112 and 334 are both near to another phenylalanine (both F2). 112 is very buried but 334 is partially exposed.

### F2

Seems to utilise specific aromatic properties as well as hydrophobicity.
They occur relatively frequently in groups.

Examples:

* BRCA1 - 43 surrounded by hydrophobic residues and filling a ring shaped space
* CBS - 332, 385 & 396 are positioned in a cleft together, along with several other aromatic residues; more potential pi-stacking? 197 & 310 are also adjacent to each other
* TP53 - 270 & 113 are adjacent and face one another. 134 is near them but in the middle of an apparently polar pocket?

### F3

Tolerates Cysteine, isoleucine, leucine, serine valine and tyrosine

Examples:

* ADRB2 - Two positions, one in the intracellular end of the transmembrane helix (223) and one in the G-protein binding domain (240). Both are accessible, although 223 might be exposed to the membrane.
* CALM1 - 17 & 69 are potentially pi-stacking
* CCR5 - many positions in the transmembrane helices, largely facing into the membrane
* TPK1 - A number of positions, all buried in generally hydrophobic pockets

### FP

Only ten phenylalanine positions are permissive.

### FO

Two outliers: APH3II 20 is generally improved by substitution (although weakly) and TEM1 17 is broadly permissive.

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

Strong overall selection, often near to other glycine residues and often appears at the ends of helices/sheets or in tight turns, like an extreme version of the 'anything but proline' positions.

Examples:

* ADRB2 - 315 is in a transmembrane helix, 83 is at the kink in another helix, 255 is in a tight turn
* CALM1 - 99 is in a tight turn
* CBS - 139, 162, 185, 345 are all in tight turns at the end of helices, various other positions in looser loops and in helices/strands
* HA - Many positions, lots at the start/end of helices/sheets or in tight loops near these ends, plus various other positions in looser loops or the middle of secondary structural features
* HSP90 - 94 is in a tight turn at the end of a helix, 118, 121 & 123 all near each other in the loop at the end of an alpha helix and 81 is at the apex of the loop in a fold
* Protein G - 9 & 41 are both shortly after beta strands and 14 is at the start of a strand.
* TP53 - 226 is in a tight loop, 266 at the start of a strand and 108 just before a hairpin loop.
* UBI - 75 & 76 are on the C-terminus, which are used to link into other ubiquitin molecules in polyubiquitin chains. These are therefore very functionally important

### G2

Selects for smaller hydrophobic residues

Examples:

* CP - in three loops and an alpha helix, near the edges of the protein, but likely not exposed as this is where the virus capsid monomers join.
* HSP90 - 83 is in the base of the ATP binding pocket, which also contains many other glycine residues. The side where an R group would protrude is cramped
* Ras - 138 is in a tight turn after a helix, 75 is in a cramped space, 115 is at the end of a beta sheet in a reasonably cramped space (several aromatics nearby)
* TEM1 - A number of positions on a beta sheet sandwiched between helices, potentially where there is not space between the secondary structure for larger residues
* UBI - 10 is next to a lysine used in polyubiquitination, which may be blocked by bulky residues. 35 is in a loop after a helix

### G3

Selects against larger hydrophobic/aromatics

Examples:

* CBS - a collection of buried positions
* HSP90 - 100 is at the mouth of the ATP binding pocket, potentially blocking it if it were larger. 153 is on the interface between HSP90 monomers, potentially interfering with this binding if mutated to a bulky residue.
* CBS - Various buried positions and some on the surface
* TPK1 - Mixture of surface and buried positions, not clear pattern
* Ubi - 47 is next to lysine 48, where binding can occur in ubiquitination. This position is mostly open, but it likely blocked by bulky aromatics

### G4

Anything but proline or isoleucine

Examples:

* APH3II - many positions in loops
* HSP90 - 170 is at the end of a beta strand
* Src - 392 is at the end of a beta strand, 440 is in the loop at the end of a helix and 347 in a free looking loop

### G5

Selects against polarity, but sometimes only weakly

Examples:

* ADRB2 - 320 is on the interface of two transmembrane helices, 257 is in the intracellular domain where the G-protein binds and 16 is in the extracellular domain
* CALM1 - 24, 26 & 62 are in loops on one end of the protein, 41 is in a loop between helices and 97 is between a helix and strand
* TPK1 - 15 is in an exterior loop and 47 at the end of a helix

### G6

Largely permissive, on average

Examples:

* ADRB2 - 35 & 276 are in transmembrane helices facing the membrane, several other positions are in the extra/intracellular domain regions not in the structure.
* APH3II - 10 & 21 are in surface loops, 234 is on the end of an alpha helix
* TPMT - 161 is at the end of a helix and 45 is in a loop

### G7

Strong selection against aspartate and lysine, with weaker selection against other charge/polarity

Examples:

* infA - three positions in a lobe near to, but not adjacent to the bound DNA
* TEM1 - two positions on helices facing into the core, but into fairly large spaces
* CXCR4 - Two positions in transmembrane helices, one facing the membrane (57) and one into the helix bundle (258). 273 is in an extracellular loop, away from the binding site.

### GP

38 positions, spread over many proteins.

### GO

One, largely permissive, outlier position (CXCR4 306)

## H - Histidine

deepSplit = 0

Summary:

* H1 - Overall strong selection, most tolerant of tyrosine, glutamine and asparagine
* HP - Permissive positions
* HO - Strong improvement when replaced with methionine or alanine, otherwise mix of weak positive and negative effects

### H1

Overall strong selection, most tolerant of tyrosine, glutamine and asparagine. The cluster appears somewhat under-split, looking at the actual position profiles

Examples:

* ADRB2 - Positions in the intracellular and extracellular domains, but none in the transmembrane helices. Looking at the actual profile all these positions have strong selection pressures beyond the subtype pattern. Four positions (92, 172, 178 & 296) are around the binding pocket
* amiE - 3, 232 & 275 are deeply buried, 281 is on the surface adjacent to a phenylalanine, 107 is just below the surface next to a tyrosine, 26 is on the surface
* BRCA1 - 41 is next to a Zn+ ion, which is also coordinated by 3 cysteines. This is a known interaction - Zhou et al. 2013
* HSP90 - 197 is in a very positive pocket near two arginines
* UBE2I - 20 is in an open cleft in the protein, 83 is buried near various, with two tyrosines somewhat nearby.

### HP

16 permissive positions across a range of proteins, generally much more accessible than H1 positions.

Examples:

* Ras - 166 is partially exposed at the end of an alpha helix
* SUMO1 - Two positions (43 & 98) both project out of the protein
* TEM1 - 24 & 156 are near to H1 positions, potentially able to interact. 94 projects out of the protein, but near a K+ ion.
* TPMT - 227 points into the core, but into a reasonably large void.

### HO

Most outliers are broadly permissive, two positions in CXCR4 (228 & 232) are intolerant to asparagine, proline, glutamine and threonine.

## I - Isoleucine

deepSplit = 0

Summary:

* I1 - Selects for hydrophobicity
* I2 - Selects for larger hydrophobic residues
* I3 - Strong selection against proline, weaker selection against leucine and lysine (?)
* IP - Permissive positions
* IO - Outliers

### I1

Selects for hydrophobicity

Examples:

* ADRB2 - Many positions in the transmembrane helices, mostly facing out into the membrane but a few also face into the helix bundle
* BRCA1 - 15 & 90 face into a helix bundle, 42 & 68 are general hydrophobic positions
* CBS - A range of internal hydrophobic core positions
* CP - 105 & 61 are in interfaces between secondary structural features, 34 is on the edge of the monomer, likely interfacing with other capsid subunits
* HSP90 - Many positions facing into the core

### I2

Selects for larger hydrophobic residues

Examples:

* ADRB2 - again a collection of positions in the transmembrane helices, but this time all facing into the core of the bundle
* CBS - Many internal positions, potentially in larger spaces than I1 positions. For example, 143 & 166 are in a beta sheet, projecting onto one side while their neighbours 142 & 167 (I1 positions) project onto the other side. The I2 side appears to be more open.
* CXCR4 - Many positions in the transmembrane helices, not clearly distributed differently to I1 positions however.
* HA - many core residues, especially in the lower lobe
* TPK1 - 153-156 are a run of four I2 positions in a row on a beta sheet, seeming to have a lot of space above and below the strand. There are various other internal positions too

### I3

Strong selection against proline, weaker selection against leucine and lysine, which is odd.

Examples:

* CXCR4 - many transmembrane helix positions, generally facing out into the membrane
* HSP90 - three positions near the start of an alpha helix (20) and two beta strands (147 & 172)
* Protein G - near the end of a beta strand (n.b. is a valine in the structure)
* ADRB2 - 58 is at the end of a transmembrane helix, 205 is at the kink in a transmembrane helix and 154 is just in the middle of a helix

### IP

22 permissive positions, less buried than other isoleucines

Examples:

* HSP90 - A number of positions, without a particular pattern
* CCR5 - one transmembrane helix position (54), fairly near the extracellular side
* BRCA1 - three positions, no strong pattern

### IO

9 outlier positions, some trend for selecting against methionine, particularly TPK1 122, but not a strong pattern of profiles

## K - Lysine

deepSplit = 1

Summary:

* K1 - Weak selection for polarity, most selective for positive charge
* K2 - Strong selection against aspartate, slightly weaker selection against glutamate, proline, glycine and aromatics
* K3 - Anything but proline
* K4 - Strongly selects against threonine & glutamine
* K5 - Strong selection against glutamate, weak selection against aspartate
* KP - Permissive positions
* KO - Outliers

### K1

Weak selection for polarity

Examples:

* BRCA1 - 45 makes a polar contact with a nearby glutamate, 75 faces towards an asparagine, but is not predicted to make a polar contact in the side chain conformations found in the structure (which may not be correct)
* CCR5 - Many positions in the intra-/extracellular domains, a couple make explicit polar contacts, most just face the solvent
* GAL4 - all three positions are adjacent to the bound DNA, with 10 & 23 predicted to make hydrogen bonds to it, and 20 would be able to with a small change in side chain configuration
* infA - two positions, facing the solvent
* Protein G - 31 makes a polar contact with a glutamate in the next turn of the alpha-helix, as well as being solvent accessible

### K2

Strong selection against aspartate, weak selection against glutamate, proline, glycine and aromatics. Some trend towards nucleic acid binding, so potentially these are substitutions that specifically disrupt this and other similar roles

Examples:

* CBS - many positions at the ends of helices or in the loops just after
* GAL4 - 18 also faces the DNA, and is between two K1 residues seemingly contacting the DNA. Potentially it isn't as important for binding but proline substitutions disturb the binding residue positions too much. 45 is also next to a arginine that looks to contact the bound DNA. Other positions are in helices
* PAB1 - three positions (131, 156 & 166) in a series of loops adjacent to the bound RNA, 166 is predicted to h-bond to the RNA but the others are not
* Src - 301 is in the loop after a beta strand

### K3

Anything but proline

Examples:

* CXCR4 - two transmembrane helix positions
* HSP90 - two nearby positions, 98 & 102, in a loop between helices
* TEM1 - five positions in secondary structures, with three being towards the ends
* UBE2I - two secondary structure positions and one in a long loop between beta sheet strands

### K4

Selects against threonine & glutamine, majority of these positions are in CXCR4

Examples:

* CXCR4 - many positions at the base of the transmembrane helices
* CBS - three surface positions

### K5

Strong selection against glutamate, weak selection against aspartate

Examples:

* ADRB2 - three positions at intracellular end of a transmembrane helix (267, 270 & 273) with another nearby (227), potentially a region needing to be positively charged, or at least not negatively charged. Two other positions further away on the intracellular domain (60 & 147)
* Src - 346 is near to a glutamate residue, which would repel if another negative residue were substituted
* Ras - 16 is near an Mg ion and the triphosphate end of GCP (in the structure), which is potentially doesn't interact well with negative charge. This position also tolerates other polar residues and aspartate.
* PAB1 - 140 near to an aspartate
* CXCR4 - 67 is somewhat near an aspartate

### KP

78 permissive positions, very surface accessible (although all K positions are)

Examples:

* BRCA1 - positions on the surface of the domain
* CBS - Surface positions, generally projecting outwards away from other side chains
* HSP90 - Mostly surface positions, although not all

### KO

Six positions, mostly broadly permissive, apart from ADRB2 372 which is strongly selective against arginine.

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

Selects against charge and polarity

Examples:

* ADRB2 - transmembrane helix positions face the membrane, plus a position in the extracellular domain
* BRCA1 - 82 faces into a helix bundle, other positions face into the core
* CBS - Many core positions, 230 faces a bound haem group
* GAL4 - 61 & 64 are on the face of the alpha helix interface between the two monomers, forming a leucine zipper type structure
* PTEN - 320 extends between two beta sheets, 23, 25 & 57 project into an alpha/beta fold

### L2

Only tolerates methionine and isoleucine (larger non-aromatic hydrophobics), but less conserved on average than L1 or L3.

Examples:

* CCR5 - three transmembrane helix positions (255, 257 & 285), 255 faces into the bundle, the others into the membrane
* HA - Many positions, all facing into the core, potentially tend to fill larger spaces
* SUMO1 - Three internally facing positions
* TP53 - Again a number of internally facing positions

### L3

Selects for hydrophobic residues

Examples:

* APH3II - Many hydrophobic internal positions
* HSP90 - Three internal hydrophobic positions, two on an alpha helix (15 & 18) and one between beta sheets (175)
* Ras - Two internal positions
* TEM1 - Large proportion of positions in the core are L3 positions
* UBI - Three core positions (50, 56 & 67) and one surface position (15)

### L4

Strong selection against lysine, weaker against other charge. More often surface accessible than L1/2/3

Examples:

* ADRB2 - 311 faces the membrane, 53 interfaces between the transmembrane helix and intracellular domain, 230 and 342 protrude into the cytoplasm
* CBS - Three standard surface positions (419, 456 & 468), plus one position (402) on the inner face of a loop that protrudes out from the protein, creating an void between loop and the main protein globule.
* CCR5 - Main positions in the transmembrane helices, most of which face the membrane, but also several internal positions
* CXCR4 - Same general pattern as CCR5
* TPK1 - Three internal positions, all near phenylalanine residues (24 & 105 around one, 28 near another

### L5

Selection against negative charge (much weaker against positive)

Examples:

* ADRB2 - two membrane facing transmembrane helix positions (45 & 84), two intracellular domain positions (145 & 266) and one extracellular position (11)
* APH3II - 103 is at the mouth of a binding cleft, with a Na & Mg ion plus a KAN molecule in the structure, 94 is at the base if this cleft but away from any of the ligands in the structure, 184 is at the dimer interface, near an arginine - potentially too strong a binding here makes dimerisation too favoured
* CALM1 - 2 surface positions and 1 buried. All look to be broadly in hydrophobic surroundings
* CXCR4 - two positions inside the transmembrane helix bundle and two on the dimerisation interface of the transmembrane helices

### L6

Anything but proline

Examples:

* GAL4 - 16 & 32 are at the end of alpha helices coordinating two zinc ions and binding DNA
* HSP90 - Several general secondary structure positions
* infA - Two loop positions
* UBI - two positions in the final beta sheet, just before the polyubiquitination point

### L7

Weak selection against isoleucine, proline, tyrosine & tryptophan (larger/branched residues)

Examples:

* TP53 - a collection of surface accessible positions
* TPK1 - Two surface accessible positions, and various generic hydrophobic internal positions
* CBS - 77 is on the face of a cleft in the protein and 397 & 423 is on the protein surface
* CCR5 - 137 is on the extracellular surface, and 174 on the intracellular domain surface

### L8

Strong selection against arginine

Examples:

* CXCR4 - 4 membrane facing positions, two facing the binding site (120 & 466), 317 is on the surface of the intracellular domain
* TPMT - 235 is in a mostly hydrophobic pocket, apart from vaguely nearby a glutamate that's already facing an arginine
* HSP90 - Two positions in the interface between helices and a beta sheet
* CALM1 - Two symmetrical positions (49 & 106) either side of the fold, pointing into the core of helix turns

### LP

46 permissive positions, often surface accessible ones

Examples:

* HSP90 - Many positions in the lower lobe, both internal and external. Potentially a region not strongly selected in the experiment
* BRCA1 - 3 surface positions and one internal

### LO

32 outlier positions, some with meaningful selection. Six positions are intolerant to isoleucine and two to valine, plus a number of weaker modes.
These are more likely to be surface positions than most leucine positions.

## M - Methionine

deepSplit = 0

Summary:

* M1 - Tolerates leucine, isoleucine, valine and threonine (longer aliphatic chains?). Strongly intolerant to negative charge and proline
* M2 - Strong selection against proline, weak selection against negative charge
* MP - Permissive positions
* MO - Outliers

### M1

Tolerates leucine, isoleucine, valine and threonine (longer aliphatic chains?). Strongly intolerant to negative charge and proline

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

Strong selection against proline, weak selection against negative charge

Examples:

* ADRB2 - 4 transmembrane helix position, two facing the membrane and two facing into the bundle
* APH3II - 148 is on the end of an alpha helix next to a glutamate
* HSP90 - 84 is in a binding pocket for ATP, 116 fills a hydrophobic cleft, 589 is in an alpha helix adjacent to a glutamate (would repel with a negative charge)
* NP - 105 is on a surface beta sheet, 222 faces inwards on a helix, 238 faces inwards on a helix next to an aspartate, 456 and 481 are on the surface
* TEM1 - 127 & 209 are on the ends of alpha helices, 67 & 180 are buried in turns at the end of helices

### MP

10 permissive positions

Examples:

* TPK1 - 1 is found to be permissive in the study, but this is potentially an artefact of the experiment or alternate start positions might be available. 193 is in a surface loop before a sheet
* SUMO1/CALM1 - both also have position 1 listed as permissive, again potentially an artefact
* HSP90 - 105 is in an open surface pocket

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

Moderately selective against everything, weakest against serine and threonine, range of surface accessibility

Examples:

* ADRB2 - 301 & 196 are in the extracellular domain facing the surface, 293 is in the extracellular binding site facing the ketone group of the P0G bound in the structure, 59, 61 & 322 are in the core of the helix bundle, all predicted to make hydrogen bonds
* APH3II - 195 looks to be coordinating a Mg ion, 58 is near to an aspartate
* CALM1 - 54 & 112 are on the surface
* CCR5 - various transmembrane positions, several internal making h-bonds and positions in the intra/extracellular domains, some on the surface and others forming internal h-bonds
* CP - 56 & 88 face into the viral capsid; 25 & 37 are on the edge of the protein, at the interface with other subunits; 126 is on another interface, also making internal h-bonds; 4 is shown covalently bonded to proline 118; 13 & 117 are on the external surface
* HSP90 - 37 & 92 face an ATP binding pocket
* PAB1 - 127 is in the RNA binding pocket, but not directly facing the RNA, 139 is a buried residue making internal h-bonds

### N2

Anything but proline, much more surface accessible

Examples:

* ADRB2 - 148 is at the top of a transmembrane helix, 183 is in a short helix at the top of the binding site
* GAL4 - 27 is in a helix near the DNA binding domain
* infA - 43 is in a short surface helix
* Protein G - 8 is at the end of a beta sheet and 37 at the end of an alpha helix
* TP53 - 263 is in a loop between beta sheet strands, 345 is in a helix, 235 & 268 are in beta strands, 239 & 246 are in a long loop between strands
* UBE2I - 37 is in a loop between strands and 124 is in a loop between helices

### NP

31 permissive positions, vey surface accessible

Examples:

* CCR5 - Two positions in the extracellular domain
* HSP90 - Various surface accessible positions
* TEM1 - Four surface positions
* UBI - 25 in on the surface in an alpha helix and 60 on a surface loop - potentially these are artefacts and are important in some contexts as all of UBI is highly conserved

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

Selects moderately against everything, somewhat buried but not all

Examples:

* CBS - 170 is in a tight turn after a strand, 59 & 64 form two tight turns creating the shape of the haem binding pocket
* GAL4 - 42 & 48 are along a loop running along the length of the bound DNA
* Src - 465 & 467 form a tight hairpin in the middle of a loop
* TPK1 - 167 & 234 are in the turns between strands of two sandwiched beta sheets

### P2

Tolerates serine, threonine, glutamine, leucine and alanine. More surface accessible

Examples:

* ADRB2 - Selection of membrane facing and extracellular domain positions
* CBS - Many positions, some surface others buried. No standout positions
* TPMT - 193 & 196 are in a loop on the edge of the binding pocket, 139 is on a surface helix
* SUMO1 - 58 is facing inwards on a surface loop, with the side chain in a hydrophobic pocket
* CP - 117 is shown covalently bonded to asparagine 3

### P3

Strong selection against aromatic, weaker against other large hydrophobic residues, surface accessible.
They tend to be in a position that could clash or on the protein surface.
Sometimes they seem to avoid a clash with a sharp turn.

Examples:

* amiE - 80 is on the surface, 333 is on a loop close to a tyrosine, which would clash
* APH3II - 76 is near a phenylalanine, 45 & 98 are on the surface, 109 is on the surface near a tryptophan, 194 is adjacent to a bound NA ion
* PTEN - 30 & 89 are on the surface
* TP53 - 142 is in a surface pocket with potential clashes with overhanging protein, 223 faces another backbone loop, which would clash, 250 is on the surface that contacts DNA
* UBE2I - 28 faces several other sidechains close by, 72 faces another backbone region, 128 faces a tyrosine, 88 is on the surface

### P4

Selects against charge and polarity, mostly buried

Examples:

* CXCR4 - Many transmembrane helix positions, generally facing into the bundle but some also face the membrane/the other monomer
* HA - 231 & 237 are buried in hydrophobic pockets
* TPMT - 68 is deeply buried, 189 is on the surface
* UBE2I - 73 is buried surrounded by hydrophobic and aromatic residues, 105 is on the surface

### PP

56 permissive positions, very surface accessible

Examples:

* HSP90 - Surface positions
* TEM1 - Surface positions
* UBI - 37 & 38 are in an external helix, but facing inwards
* CBS - Surface positions (but so are many other P positions)

### PO

Three positions are broadly weakly selective, particularly against histidine, Src 506 appears to be improved be becoming positively charged (K/R), the other 8 are generally permissive.

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

Generally selective, most tolerant to polarity or charge, strongest against aromatics

Examples:

* ADRB2 - 4 positions on the intracellular domain, 142, 247 & 337 are on the surface, 229 is slightly buried making an internal h-bond
* amiE - 8 buried positions making internal h-bond
* HSP90 - 14 is on the surface of a cleft into the protein
* PTEN - 149 is on the surface
* Src - 278 & 394 are surface positions
* SUMO1 - 4 surface positions, 55 & 92 are predicted to make h-bonds as well whereas 69 & 94 are not
* TPK1 - 134 is in a surface binding cleft, where PYI is bound in the structure. 109 is predicted to make h-bonds to E113 up the alpha helix

### Q2

Anything but proline

Examples:

* ADRB2 - 224 is at the end of a transmembrane helix, 231 is in the intracellular domain alpha helix
* HSP90 - 9 is in a chain between secondary structure elements, 119 is in a long loop, 135 is on the surface in the middle of a beta strand
* Ras - 25 is at a kinked end of a helix, 99 & 129 are near the end of regular helices
* TP53 - 104 & 165 are in hairpin turns, 136 is just after a strand and 331 at the end of one
* UBI - 31 & 41 are before/after a tight turn with a short helix,

### Q3

Weak selection against hydrophobic residues

Examples:

* CALM1 - 6 surface positions
* CCR5 - One extracellular surface accessible position (313), 7 intracellular surface positions, 261 is a buried position making an h-bond
* Ras - Two surface positions
* UBE4B - Two surface positions

### Q4

Strong selection against negative charge

Examples:

* ADRB2 - 170 & 179 are surface positions on the extracellular domain, 65 is directed at the base of a transmembrane helix, 250 is internal in the intracellular domain and 243 is on its surface
* GAL4 - 9 is near the backbone of the bound DNA
* SUMO1 - 4 & 29 are surface positions
* TPK1 - 9 is in the PYI binding pocket and near an Mg ion, 166 is accessible but directed inwards at a hydrophobic helix (it h-bonds to a threonine)

### Q5

Improves on average, apart from lysine

Examples:

* YAP1 - 186 is on the surface but next to a tryptophan
* CXCR4 - 210 is in the extracellular binding domain, 212 is on the membrane surface, in the dimerisation interface
* TPK1 - 26 and 157 are adjacent to each other near the surface, 144 & 190 face phenylalanine residues, 194 is very surface accessible

### QP

54 permissive positions.

Examples:

* amiE - surface positions
* CCR5 - 194 is inside the bottom of the bottom of the transmembrane domain, mostly near hydrophobic positions
* HSP90 - partially exposed positions, generally have hydrophobic and hydrophilic  residues nearby
* TEM1 - 5 surface positions

## R - Arginine

deepSplit = 1

Summary:

* R1 - Intolerant of everything apart from lysine
* R2 - Anything but proline
* R3 - Selects against negative charge and proline
* R4 - Weak selection
* R5 - Selects against negative charge, aromatics and proline
* RP - Nonselective
* RO - Outliers, nonselective

### R1

Intolerant of everything apart from lysine

Examples:

* APH3II - 211 & 226 are in the binding pocket with KAN, 66 & 183 face aspartate residues
* CALM1 - 75 & 87 face glutamate residues, 91 is on the surface
* CBS - 51, 224 & 266 are in the haem binding pocket, 182 is near an aspartate, 161 is near a glutamate, other residues on the surface
* GAL4 - 15, 46 & 51 are near the DNA backbone
* NP - many surface positions

### R2

Anything but proline

Examples:

* APH3II - Four positions in the middle of helices, 122 is at the beginning of a helix and 202 at the beginning of a beta sheet, 44 is in a tight turn between strands in a sheet
* HA - 472 is in a tight turn between a helix a strand, 224 is in a turn between strands
* Src - 480 is at the end of a helix with a tight turn after it, 419 is in a turn between helices, 463 is in a generic loop
* TPK1 - 33 & 51 are in alpha helices, 171 is in a surface beta sheet

### R3

Selects against negative charge and proline

Examples:

* CCR5 - Two surface positions at the bottom of the transmembrane domain
* Ras - 41 faces an aspartate in another bound Ras monomer
* TEM1 - 255 is next to a glutamate, 271 is somewhat near the bound CB4
* TPMT - 140 & 215 are on the surface, 163 is somewhat near a glutamate, 152 is at the base of a binding pocket

### R4

Weak selection

Examples:

* ADRB2 - 304 is exposed in the extracellular domain, 151 is membrane facing at the bottom of the transmembrane domain, 221 is near a glutamate in the intracellular domain, 260 is on the surface of the intracellular domain.
* TPMT - 226 is partially exposed with both hydrophobic and hydrophilic residues
* SUMO1 - 70 is on the surface, near an R2
* TP53 - 213 is buried with mainly backbone nearby

### R5

Selects against negative charge, aromatics and proline

Examples:

* ADRB2 - 328 is exposed in the intracellular domain, potentially near where the G-protein binds
* APH3II - 253 is near two aspartate residues, 258 is buried in a reasonably cramped position
* TEM1 - 241 is in a tight interface between a sheet and helices, 162 is on the surface near another arginine and a glutamate, 202 is mostly buried with the amine group exposed
* UBE2I - 104 is very exposed on the surface

### RP

43 permissive positions

Examples:

* TP53 - Three very exposed residues
* HSP90 - Many surface positions
* TEM1 - Many exposed positions

### RO

Three outlier positions, all slightly improved by mutation.

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

Only tolerates threonine, generally the most buried subtype with the strongest decrease in hydrophobic solvation when mutated

Examples:

* ADRB2 - 111, 161, 165 & 220 are in transmembrane helices, predicted to form polar contacts with the backbone, 246 forms an h-bond linking the intracellular domain to the transmembrane helices, 329 is at the base of the short perpendicular intracellular helix, 203 is in the binding site forming an h-bond to the bound P0G in the crystal structure
* amiE - A large number of buried positions forming polar contacts
* CBS - 63 is in the haem binding pocket, many other internal positions making internal h-bonds
* CCR5 -  75, 114 & 160 form internal h-bonds in the transmembrane domains, 6, 17 & 270 are surface positions in the intracellular domain
* GAL4 - 22 & 47 form h-bonds on the protein surface
* HSP90 - 4 internal positions making h-bonds
* Ras - 17 coordinates a Mg ion, 89 makes internal h-bonds
* TEM1 - 68 binds the bound CB4 ligand, 51 is on the surface, all other positions form internal h-bonds
* TPK1 - 216 & 218 are in a binding pocket, making polar contacts with the bound PYI ligand in the crystal structure, 62 & 243 are on the surface and 74 is exposed but predicted to form an internal h-bond

### S2

Anything but proline, generally more surface accessible

Examples:

* ADRB2 - 207 is in a kink in a transmembrane helix, 137 in a tight turn between helices, 120 in the middle of a transmembrane helix, 261 & 262 are in a coil between transmembrane helix and intracellular domain
* APH3II - 68 & 114 are in alpha helices, 40 at the end of a strand, 204 in the middle of a strand, 32 & 104 sit in coils
* HSP90 - 109 is in a turn before a helix and 198 at the end of a helix
* TP53 - 240 is just before a S1 position contacting bound DNA, 261 in a turn between secondary structures and 215 is in the middle of a sheet

### S3

Selects against negative charge, wide range of surface accessibility

Examples:

* CALM1 - 82 is between a D and E and 39 is at the base of a helix
* CCR5 - 38 faces the membrane, 180 is buried
* HSP90 - 36 & 39 are adjacent on one side of a helix, near to several negatively charged residues
* PAB1 - 155 is buried
* TEM1 - 233 is buried, near to an K ion and various other charged residues

### S4

Weak selection against polar, very surface accessible

Examples:

* amiE - 233 is deeply buried surrounded by hydrophobic residues, 272 is buried near an asparagine and aspartate and 290 is on the surface but facing inwards towards hydrophobic residues
* CALM1 - 18 & 102 are exposed on the surface
* CXCR4 - 131 is buried within the transmembrane helices, 217 is exposed to the membrane and on the dimerisation interface, 285 is in the binding pocket, 178 is exposed on the extracellular surface
* TPK1 - 92 is on the surface, 138 is buried surrounded by hydrophobic residues
* TPMT - 229 is exposed on the surface

### S5

Anything but alanine, very surface accessible

Examples:

* ADRB2 - 204 is int he wall of the binding site, 236 is in the intracellular domain
* Ras - 106 is on the protein surface
* TP53 - 4 surface accessible positions
* UBE2I - 7 is on the protein surface

### SP

60 permissive positions

Examples:

* CCR5 - 63 is in the loop between transmembrane helices
* HSP90 - 17, 25 & 155 are on the dimerisation interface, other positions are in a wide range of conditions
* TEM1 - Four surface positions
* UBI - Five externally facing positions

### SO

14 outlier positions, mostly permissive apart from two positions in ADRB2 - 74 & 41 select against asparagine with 74 also selecting against proline.

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
* TO - Outliers

### T1

Only tolerates serine, generally buried

Examples:

* ADRB2 - 25 is an exposed extracellular position, 68 & 136 form h-bonds in the transmembrane helices
* BRCA1 - 37 makes internal h-bonds, 69 & 77 are on the surface of the domain in the crystal structure, but it's not known where they are in the overall structure
* CALM1 - 35 is fully exposed, 27 is on the surface but making h-bonds to other residues
* CBS - 296, 383, 428 & 494 are surface accessible, the other positions are buried making internal h-bonds, 150,  257 & 260 interact with a bound PLP
* CCR5 - 5 internal positions in the intracellular domain, various other positions make internal h-bonds between transmembrane helices
* HSP90 - 171 is in the ATP binding pocket, 101 makes internal h-bonds near the bound ATP, 22 is buried

### T2

Anything but proline

Examples:

* ADRB2 - Four transmembrane helix positions, one of which (110) is in the binding domain wall
* APH3II - 73 is at the end of an alpha helix, 85 is at the end of a beta sheet, 130 is in a short alpha helix turn, 28 is in a beta turn and 146 is at the end of a helix
* HSP90 - 157 is in a strand on the surface, 13 is at the start of an alpha helix and 95 in the hairpin turn after a helix
* UBE2I - 135 is in a short surface helix, 29 at the end of a strand

### T3

Selects against aromatics and (less so) negative charge

Examples:

* CALM1 - 69 & 77 are both on the surface of the protein and have nearby residues that could clash
* SUMO1 - 41 & 42 are on the surface, near the dimerisation domain in the crystal structure, 95 is exposed but with residues that can potentially clash
* HA - 147 and 151 are exposed, 228 is on the surface of a large cleft, partly facing another lobe of the protein
* TEM1 - 27 is on the interface between two helices, 139 is on the surface but facing several large residues, 147 is on the surface but facing parallel to it towards other residues

### T4

Selects against negative charge, glycine and proline

Examples:

* ADRB2 - 283 faces the membrane, 164 faces another transmembrane helix and 189 faces the extracellular surface
* CALM1 - 45 is on the surface of the cleft between lobes
* NP - 232 and 296 face each other on the protein surface, 23 & 130 are on the surface
* Src - 292 is at the end of a strand on the surface, but facing another protein surface, 443 is on the protein surface
* TP53 - 284 directly faces the DNA backbone, 25 is at the surface near the DNA backbone, 230 is at the end of a strand near various cramped residues, 256 is on the surface in a sheet

### T5

Weak selection against cysteine only, very surface accessible

Examples:

* ADRB2 - 100 & 177 are on the extracellular domain surface, 66 & 146 are on the surface of the intracellular domain
* CXCR4 - 73 is at the interface between transmembrane helices, 142 is in a turn between helices at the intracellular end, facing several charged residues in the intracellular domain
* Ras - Two exposed surface positions
* Src - 293 is partially exposed near charged and polar residues and 304 is on the protein surface

### T6

Strong selection against asparagine, weaker selection against other positive charge, negative charge, aromatics and cysteine

Examples:

* ADRB2 - 281 faces the membrane and 195 faces into the binding domain
* CALM1 - Five surface accessible positions
* CBS - 333 is very accessible, 460 is buried with two negatively charged residues close by as well as a tyrosine
* TPK1 - 5 surface accessible positions, 141 is in a narrow interface between helices, 149 is partially buried forming h-bonds to a nearby glutamine

### TP

58 permissive positions

Examples:

* CBS - 53 is on the surface just below the haem binding pocket, 193 & 318 are on the surface, 141 is in a mostly hydrophobic pocket between a helix and sheet
* HSP90 - 5 positions, mostly exposed
* Protein G - Mostly surface positions, apart from 16 which is semi exposed but also partly buried between the sheet and helix
* TPMT - 190 is exposed on the surface, 38 is buried with aromatic and hydrophobic residues nearby

### TO

Seven permissive positions

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

Selects for hydrophobic residues and threonine, strongly buried

Examples:

* ADRB2 - Many transmembrane helix domains, some facing the membrane and some in the helix interfaces
* APH3II - Many mostly buried domains, near to other hydrophobics
* CBS - Many buried hydrophobic positions
* HSP90 - Buried hydrophobic core positions
* PTEN - 343 is on the hydrophobic interface between beta sheets, 53 & 54 face inwards from a surface alpha helix
* Src - One buried core position, three surface positions

### V2

Only tolerates isoleucine, buried

Examples:

* ADRB2 - similar distribution to V1 residues, not clear what the difference is
* CP - 11, 19, 21 & 45 are in the interface between beta sheets, 29 is on the internal capsid face and 74 is on the interface to other capsid monomers
* HA - Appear to be general buried positions, similar to V1 again
* PTEN - 217 is in the interface between sheets
* Src - 284 lines a cleft, 397 and 405 are buried in hydrophobic pockets
* TP53 - 122 directly faces the bound DNA, 143, 197 & 216 are in the interface between sheets, 274 is in a hydrophobic pocket

### V3

Selects against charged (apart from glutamate) and some polar, particularly histidine, less buried

Examples:

* amiE - buried hydrophobic positions, again look visually similar to V1/2
* CALM1 - 36 & 143 are on interfaces between helices
* CCR5 - Many transmembrane domain positions, most of which face the membrane
* TPK1 - 208 & 226 are on the surface, 139 & 188 are deeply buried in the hydrophobic core

### V4

Similar to V3 but selects against glutamate instead of aspartate and particularly disfavours glutamine, also more surface accessible

Examples:

* ADRB2 - seven membrane facing positions in the transmembrane helices
* amiE - Many buried positions interfacing between helices/sheets
* CALM1 - 109 & 122 are both partially accessible but interfacing with other hydrophobic residues
* TEM1 - 29 & 101 are both on the protein surface
* TPMT - 23 & 34 are surface accessible and 156 is buried in a hydrophobic pocket

### V5

Anything but proline, somewhat buried

Examples:

* APH3II - 242 is in a surface helix, 83, 84 & 198 are in strands, 46 is just before a strand and 75 is in the turn after a helix
* CXCR4 - Two transmembrane helix positions and 177 is in a strand forming part of the binding domain
* HSP90 - 74, 136 & 158 are at the ends of strands in a beta sheet, 201 is in the turn between a helix and strand
* TEM1 - 82 & 106 are both at the ends of helices

### VP

25 permissive positions

Examples:

* HSP90 - 6 surface positions
* Protein G - 21 is exposed on one end
* CBS - 252 is in the interface between a sheet and a helix, 414 is in a buried pocket that's mostly backbone

### VO

29 outlier positions, mostly nonselective apart from a few positions with a strong selection pressure - for example ADRB2 114 strongly selects against leucine only and Src 326 is improved by substituting arginine

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

Three permissive positions:

* TEM1 163 is flat to the surface
* CBS 54 is buried but its backbone next to the haem binding pocket
* CBS 410 is partially exposed flat to a an alpha helix

### WO

Two outliers, both in Src. 285 is on the underside of a beta sheet, but seems to almost clash with the backbone of the neighbouring strand and is strongly improved when mutated to aliphatic hydrophobic residues. 289 is at the end of the same strand and is broadly permissive apart from lysine mutations which are measured to slightly improve fitness, potentially because of a nearby glutamate.

## Y - Tyrosine

deepSplit = 1

Summary:

* Y1 - Tolerates aromatics, particularly intolerant of charge (other than aspartate) and glycine (?)
* Y2 - Intolerant of charged and polar residues, most tolerates phenylalanine, histidine and hydrophobic residues
* Y3 - Generally selective, somewhat tolerates phenylalanine and histidine
* Y4 - Anything but proline
* YP - Permissive positions
* YO - Outliers

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

Generally selective, somewhat tolerates phenylalanine and histidine

Examples:

* ADRB2 - 308 is in the extracellular binding site, 219 is in the transmembrane bundle, h-bonding to a fairly distant arginine, 70 is at the intracellular end of the transmembrane helices facing outward
* CBS - Again a group appears in a sort of cluster (223, 308, 301) along a line, with a Y1 position (381) as well. 163 appears alone in a reasonably large pocket (potentially needing the size of the aromatic ring) and 484 fills a pocket near a pair of aromatics (phenylalanine/histidine)
* HA - many positions, somewhat clustered at either in the binding domain lobe.
* TPK1 - 53 is in a buried cleft with the alcohol group exposed and the aromatic group facing hydrophobic residues, 163 & 221 are buried near other aromatic/hydrophobic residues
* TP53 - 163 is facing two glutamate residues, 205 is buried in a hydrophobic pocket

### Y4

Anything but proline

Examples:

* CCR5 - 214 is buried in the transmembrane helix, at a kink in one of the helices, 307 is in the extracellular alpha helix
* CXCR4 - 4 transmembrane helix positions, some buried, some exposed
* PAB1 - 197 is in a buried loop, 143 is in a short helix, both near the bound RNA
* Src - 439 is at the end of an exposed helix, 329 is at the end of a buried strand with a tight turn after it

### YP

Seven permissive positions:

* Two infA positions (44 & 60) that are different amino acids in the homology model.
* HSP90 - 47 is in a wide cleft between an alpha helix and the rest of the protein, 146 is on the surface, 184 is partially buried with a glutamate nearby
* SUMO1 - 21 is full exposed
* amiE - 339 is mostly exposed

### YO

Two largely permissive outlier positions
