# Combined Mutational Scans Dataset EV2

Full deep landscape dataset, containing all positions normalised ER scores and additional annotations

## Columns

study: Study the position comes from
gene: Gene
position: Position in gene (relative to cannonical Uniprot sequence)
wt: Wild-type amino acid
A-Y: Normalised ER score when the position is mutated to each amino acid
log10_sift_A-Y: Log10 SIFT4G score for each substitution at this position
mean_score: Mean normalised ER score across all substitutions
mean_sift: Mean log10 SIFT4G score across all substitutions
total_energy: Mean FoldX ddG prediction (kJ/mol) across all substitutions
backbone_hbond-energy_ionisation: Mean prediction for each energy term in FoldX's force field model
phi: Phi backbone angle
psi: Psi backbone angle
all_atom_abs: All atom absolute surface accessibility as calculated by Naccess 
all_atom_rel: All atom relative surface accessibility as calculated by Naccess 
side_chain_abs: Side chain absolute surface accessibility as calculated by Naccess 
side_chain_rel: Side chain relative surface accessibility as calculated by Naccess 
backbone_abs: Backbone atom absolute surface accessibility as calculated by Naccess 
backbone_rel: Backbone atom relative surface accessibility as calculated by Naccess 
non_polar_abs: Non-polar absolute surface accessibility as calculated by Naccess 
non_polar_rel: Non-polar relative surface accessibility as calculated by Naccess 
polar_abs: Polar absolute surface accessibility as calculated by Naccess 
polar_rel: Polar relative  surface accessibility as calculated by Naccess 
within_10_0_A-Y: Count of each amino acid within a 10 angstrom sphere centered at the position
angstroms_to_A-Y: Distance to each amino acid, in angstroms
ss_g-t: Porter5 predictions for each DSSP secondary structure class
hydrophobicity: Hydrophobicity score from Bandyopadhyay adn Mehler (2008)
PC1-20: Principal component projection of the normalised ER score profile
tSNE1-2: tSNE coordinates of the normalised ER score profile
umap1-2: UMAP projection of the normalised ER score profile
