load structures/2w22_c64s_c295s_context.pdb, BTL2_2W22_cysfree
hide everything
show cartoon, BTL2_2W22_cysfree
color gray80, BTL2_2W22_cysfree

select godoy_btl2_sites, chain A and resi 236+333+342+39+93+187+195
show sticks, godoy_btl2_sites
color magenta, godoy_btl2_sites

select btl2_top_final_sites, chain A and resi 181+182+200+222+218+180+278+226+194+198
show spheres, btl2_top_final_sites
color cyan, btl2_top_final_sites

select cysfree_background_sites, chain A and resi 65+296
show sticks, cysfree_background_sites
color orange, cysfree_background_sites

select btl2_lysines, resn LYS
show sticks, btl2_lysines
color yellow, btl2_lysines

set sphere_scale, 0.35
set cartoon_transparency, 0.15
zoom godoy_btl2_sites, 12

# Magenta: Godoy experimental Cys sites.
# Cyan: representative top-ranked CysMutML final-score candidates.
# Orange: PDB-numbered Cys-free background substitutions corresponding to experimental C64S/C295S.
# Yellow: Lys residues.
