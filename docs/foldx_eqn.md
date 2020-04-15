# Force Field
## Overall equation:
ΔG = Wvdw \* ΔGvdw + WsolvH \* ΔGsolvH +  WsolvP \* ΔGsolvP + ΔGwb + ΔGhbond + ΔGel + ΔGKon + Wmc \* T \* ΔSmc + Wsc \* T \* ΔSsc

Basically:
ΔG = ΔH - TΔS

## Terms:
ΔGvdw = Sum of van der waals of all atoms wrt same interactions with solvent 
ΔGsolvH = Difference in solvation energy for apolar groups 
ΔGsolvP = Difference in solvation energy for polar groups
ΔGwbi = Extra stabalising energy for H2O making 2 h-bonds to protein 
ΔGhbond = Energy change from forming h-bonds within the protein rather than with solvent 
ΔGel = Electrostatic contribution of charged groups, including helix dipole
ΔGKon = Electrostatic effect on association constant k[on]
T = Temperature (K) 
ΔSmc = Entropy of fixing main chain in conformation
ΔSsc = Entropy of fixing side chains in conformation
Wxxx = Weighting terms, all 1 apart from Wvdv = 0.33

# BuildModel output terms
From: http://foldxsuite.crg.eu/command/BuildModel

Total Energy - predicted overall stability
Backbone Hbond - contribution of backbone Hbonds
Sidechain Hbond - contribution of sidechain-sidechain and sidechain-backbone Hbonds
Van der Waals - contribution of the VanderWaals
Electrostatics - electrostatic interactions
Solvation Polar - penalization for burying polar groups
Solvation Hydrophobic - contribution of hydrophobic groups
Van der Waals clashes - energy penalization due to interresidue VanderWaals’ clashes
Entropy Side Chain - entropy cost of fixing the side chain
Entropy Main Chain - entropy cost of fixing the main chain
Sloop Entropy - ONLY FOR ADVANCED USERS
Mloop Entropy - ONLY FOR ADVANCED USERS
Cis Bond - cost of having a cis peptide bond
Torsional Clash - intraresidue VanderWaals’ torsional clashes
Backbone Clash - backbone-backbone VanderWaals. *not considered in the total*
Helix Dipole - electrostatic contribution of helix dipole
Water Bridge - water bridges
Disulfide - disulfide bonds
Electrostatic Kon - electrostatic interaction between molecules in the precomplex
Partial Covalent Bonds - interactions with bound metals
Energy Ionisation - ionisation energy
Entropy Complex - entropy cost of forming a complex
