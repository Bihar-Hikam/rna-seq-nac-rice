# =========================================
# PyMOL Script - Protein Structure Comparison
# LOC_Os03g60080 (Nipponbare vs Pokkali vs N22)
# =========================================

# ===== LOAD STRUCTURES =====
load Nipponbare_Os03g60080.pdb, Nipponbare
load Pokkali_Os03g60080.pdb, Pokkali
load N22_Os03g60080.pdb, N22

# ===== BASIC REPRESENTATION =====
hide everything
show cartoon

# ===== COLORING =====
color cyan, Nipponbare
color magenta, Pokkali
color gold, N22

# ===== ALIGNMENT =====
align Pokkali, Nipponbare
align N22, Nipponbare

# ===== IMPROVE VISUALIZATION =====
bg_color white
set cartoon_transparency, 0.2

# ===== OPTIONAL: SHOW DIFFERENCES =====
# highlight regions (optional if you know mutation position)
# select mut_site, resi 100
# show sticks, mut_site
# color red, mut_site

# ===== RENDER HIGH QUALITY IMAGE =====
ray 1200,1200

# ===== SAVE OUTPUT =====
png LOC_Os03g60080_structure_comparison.png, dpi=300

# ===== SAVE SESSION =====
save LOC_Os03g60080_session.pse
