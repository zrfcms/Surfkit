#step 1 align Al2O3.POSCAR to desired condition
./../../bin/Surfkit --align Al2O3.POSCAR step1.POSCAR a x b y
#step 2 align Al2O3.POSCAR to desired condition
./../../bin/Surfkit --prim2D step1.POSCAR step2.POSCAR
#step 3 scan Al2O3.POSCAR for atomic layer
./../../bin/Surfkit --scan step2.POSCAR -7 7 
#step 4 we found a non-polar center at 1 so we shift model by -1
./../../bin/Surfkit --shift step2.POSCAR step3.POSCAR -1
#step 5 scan modified Al2O3.POSCAR for atomic layer
./../../bin/Surfkit --scan step3.POSCAR -3 3 
#step 6 cleave model to Al2O3 slab
./../../bin/Surfkit --cleav step3.POSCAR step4.POSCAR -2 2
#step 7 adding vacuum layer for Al2O3 
./../../bin/Surfkit --vacuum step4.POSCAR step5.POSCAR 5


