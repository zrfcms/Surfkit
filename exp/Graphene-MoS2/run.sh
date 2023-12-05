#step1 search for twist angle
#miln is suggested when searching large scale system
./../../bin/Surfkit --match  GRAPHENE.POSCAR MoS2.POSCAR mlen 50 mbox 50 mtwi 10 mvec 0.02 magm 125 migm 115 mang 0.1

#step2 open excel select one twist angle you want to build
#this step can also controled by python or any other program which support CSV format
#sort searched result is permitted however changing data is not allowed

#step3 pick twisted slab 
./../../bin/Surfkit --twist GRAPHENE.POSCAR MoS2.POSCAR S37.POSCAR 1.871770

#step4 adjust spacing between layers
./../../bin/Surfkit --space S37.POSCAR S37.POSCAR 3 1000 auto
