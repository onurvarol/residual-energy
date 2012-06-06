Residual Energies Contributions
Copyright (c) 2012, Onur Varol 
http://onurvarol.com

Usage: residual_energy [OPTIONS]
  -nc   Nresidue: number of residue
  -f    Nrame   : number of frame
  -cmap CMAP    : name of contact map file which has 0-1 values in every column and rows (default "cmap_temp.pr")
  -dr   Dr file : name of the displacement value file which has nc row, f column (default "dr_temp.pr")
  -out  Out file: name of the output file of the algoritm (default "U_all.pr")

Algorithm: Algorithm generates contribution of each residues to internal energy
using displacement of residues during MD simulations. Pairwise correlations are
generated using GNM [1].

[1] Burak Erman. 2011. Relationship between ligand binding sites, protein architecture
and correlated paths of energy and conformational fluctuations.

Please see the file LICENSE for terms of use.