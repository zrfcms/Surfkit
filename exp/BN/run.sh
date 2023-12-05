#step 1 convert BN.POSCAR to conventioal cell
./../../bin/Surfkit --conv BN.POSCAR BN_conv.POSCAR 
#step 1.5 convert BN.POSCAR to 2D primitive cell
./../../bin/Surfkit --prim2D BN.POSCAR BN_prim2D.POSCAR 
#step 2 convert conventioal cell to supercell
./../../bin/Surfkit --proj -mat BN_conv.POSCAR BN_slab0.POSCAR 1 -1 0 1 0 -1 1 1 1 
#step 3 scan atomic layer of BN_slab0.POSCAR
./../../bin/Surfkit  --scan BN_slab0.POSCAR -5 5
#step 4 redefine slab thickness
./../../bin/Surfkit  --cleav BN_slab0.POSCAR BN_slab1.POSCAR -2 3
#step 5 redefine Vacuum thickness
./../../bin/Surfkit  --vacuum BN_slab1.POSCAR BN_slab2.POSCAR 10
