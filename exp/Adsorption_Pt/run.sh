#step 1 get site list from surface slab
./../../bin/Surfkit --site3D Pt211.POSCAR 3D.POSCAR 1 9
#step 2 adsorb H atom to surface slab
./../../bin/Surfkit --adsorb -ele Pt211.POSCAR 3D.POSCAR Slab.POSCAR 3 1 H 
