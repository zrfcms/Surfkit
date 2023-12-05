#add adsorbtion on FeS2
./../../bin/Surfkit --proj -mat FeS2.POSCAR step1.POSCAR 1 0 0 0 1 0 0 0 1
./../../bin/Surfkit --align step1.POSCAR step1.POSCAR a x b y 
./../../bin/Surfkit --cleav step1.POSCAR step2.POSCAR 0 15
./../../bin/Surfkit --vacuum step2.POSCAR step3.POSCAR 15
./../../bin/Surfkit --site3D step3.POSCAR 3D.POSCAR 3 9
./../../bin/Surfkit --adsorb -ele step3.POSCAR 3D.POSCAR Slab.POSCAR 1 1 H 
rm step*.POSCAR
