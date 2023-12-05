QB = ./QB/QB_dump.cpp ./QB/QB_function.cpp ./QB/QB_neighbor.cpp ./QB/QB_network.cpp ./QB/QB_pbc.cpp ./QB/QB_Readlmc.cpp ./QB/QB_Readlmp.cpp ./QB/QB_Readvasp.cpp ./QB/QB_create.cpp 
SPGLIB =./QSPG/spglib/arithmetic.cpp ./QSPG/spglib/cell.cpp ./QSPG/spglib/debug.cpp ./QSPG/spglib/delaunay.cpp ./QSPG/spglib/hall_symbol.cpp \
 ./QSPG/spglib/kgrid.cpp ./QSPG/spglib/kpoint.cpp ./QSPG/spglib/mathfunc.cpp ./QSPG/spglib/niggli.cpp ./QSPG/spglib/pointgroup.cpp \
./QSPG/spglib/primitive.cpp ./QSPG/spglib/refinement.cpp ./QSPG/spglib/site_symmetry.cpp ./QSPG/spglib/sitesym_database.cpp ./QSPG/spglib/spacegroup.cpp \
./QSPG/spglib/spg_database.cpp ./QSPG/spglib/spglib.cpp ./QSPG/spglib/spin.cpp ./QSPG/spglib/symmetry.cpp ./QSPG/spglib/determination.cpp ./QSPG/spglib/overlap.cpp 

QSPG = ./QSPG/QSPG_link.cpp ./QSPG/QSPG_minimize.cpp ./QSPG/QSPG_symmetry.cpp ./QSPG/QSPG_database.cpp ./QSPG/QSPG_setbox.cpp ./QSPG/QSPG_refine.cpp ./QSPG/QSPG_wyckoff.cpp

SK = ./SK/SK_orientation.cpp ./SK/SK_minimize.cpp ./SK/SK_cellscan.cpp ./SK/SK_math.cpp ./SK/SK_layermatch.cpp ./SK/SK_cylindermatch.cpp ./SK/SK_match.cpp ./SK/SK_check.cpp ./SK/SK_rigid.cpp ./SK/SK_Guassian.cpp

SRC = main_layer.cpp main_slab.cpp main_scan.cpp main.cpp $(QB) $(QSPG) $(SPGLIB) $(SK)

OBJ = $(patsubst  %.cpp, %.o, $(SRC))

Surfkit: $(OBJ)
	g++ $(OBJ) -g -o $@
	cp ./Surfkit ./../bin

%.o: %.cpp
	g++ -c -MMD $< -g -o  $@ 

	
clean:
	rm $(OBJ)