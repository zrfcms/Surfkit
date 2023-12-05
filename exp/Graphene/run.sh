#step1 search for twist angle
#miln is suggested when searching large scale system
./../../bin/Surfkit --match GRAPHENE.POSCAR GRAPHENE.POSCAR maln 230 miln 220 mbox 230 mtwi 6 mvec 0.00001 magm 122 migm 118 mang 0.01 

#step2 open excel select one twist angle you want to build
#this step can also controled by python or any other program which support CSV format
#sort searched result is permitted however changing data is not allowed

#step3 pick twisted slab 
./../../bin/Surfkit --twist GRAPHENE.POSCAR GRAPHENE.POSCAR S8011.POSCAR 1.1
