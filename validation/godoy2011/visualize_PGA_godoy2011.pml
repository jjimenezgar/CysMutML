load structures/1k5q.pdb, PGA_1K5Q
hide everything
show cartoon, PGA_1K5Q
color gray80, PGA_1K5Q

select godoy_pga_sites, (chain A and resi 86) or (chain B and resi 9+201+112+361+380)
show sticks, godoy_pga_sites
color magenta, godoy_pga_sites

select pga_top_final_sites, chain A and resi 208+128+204+3+112+168+201+109+113
show spheres, pga_top_final_sites
color cyan, pga_top_final_sites

select pga_lysines, resn LYS
show sticks, pga_lysines
color yellow, pga_lysines

set sphere_scale, 0.35
set cartoon_transparency, 0.15
zoom godoy_pga_sites, 12

# Magenta: Godoy experimental Cys sites.
# Cyan: representative top-ranked CysMutML final-score candidates.
# Yellow: Lys residues.
