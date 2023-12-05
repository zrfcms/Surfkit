#step 1 get site list from surface slab
./../../bin/Surfkit --site2D Graphene.POSCAR 2D.POSCAR 1 3
#step 2 adsorb H2O to surface slab
./../../bin/Surfkit --adsorb -file Graphene.POSCAR 2D.POSCAR Slab.POSCAR 1 5 H20.poscar H 1 1 0 0 0 1 0 0 0 1
