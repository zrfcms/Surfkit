#include"./spglib/spglib.h"
typedef struct
{
	double (*pos)[3];
	int *type;
	int num;
	double mat[3][3];
	int ele_n;
	char (*elename)[32];
}SPG_tools;
typedef struct
{
	int  number;
	char international_short[11];
	char choice[6];
	char international_full[20];
	char pointgroup_schoenflies[4];
	int  holohedry;
	char pointgroup_symbol[6];
}QSPG_data;
#define QSPGDataset SpglibDataset
#define QSPGDataset_free_dataset spg_free_dataset
extern void SPG_free(SPG_tools* SPG);
extern void SPG2QB(SPG_tools* SPG,QB_tools *QB);
extern void QB2SPG(QB_tools *QB,SPG_tools* SPG);

extern void QSPG_minimize(QB_tools* output,QB_tools* input,double symprec);
extern void QSPG_refined(QB_tools* output,QB_tools* input,double symprec);
extern void QSPG_primitive(QB_tools* output,QB_tools* input,double symprec);

extern void QSPG_get_symmetry(QB_tools*QB,QSPGDataset**data,double symprec);

extern void QSPG_add_symmetry_atom(QB_tools*input,QB_tools*output,int hall_number,int type,double x,double y,double z,double symprec);
extern void QSPG_add_operation_atom(QB_tools*input,QB_tools*output,int symmetry_num,int (*rot)[3][3],double (*tra)[3],int type,double x,double y,double z,double symprec);

extern void QSPG_exchange_symmetry_atom(QB_tools*input,QB_tools*output,int hall_number,int type,double x,double y,double z,double symprec);
extern void QSPG_exchange_operation_atom(QB_tools*input,QB_tools*output,int symmetry_num,int (*rot)[3][3],double (*tra)[3],int type,double x,double y,double z,double symprec);

extern void QSPG_simplify_symmetry_atom(QB_tools*input,QB_tools*output,int hall_number,double symprec);
extern void QSPG_simplify_operation_atom(QB_tools*input,QB_tools*output,int symmetry_num,int (*rot)[3][3],double (*tra)[3],double symprec);

extern void QB_f2c(QB_tools *QB,double input[3],double output[3]);
extern void QB_c2f(QB_tools *QB,double input[3],double output[3]);

extern void SPG_pick_site_symmetry( int* multiplicity,
									int* wyckoff_letter,
									char site_sym_symbol[7],
									const double position[3],
									SPGCONST double bravais_lattice[3][3],
									const int hall_number,
									const double symprec);
//set box
extern void QSPG_getbox(int mode,int R_flag,double lattice[3][3],double ia,double ib,double ic,double ialpha,double ibeta,double igamma);
extern void QSPG_setbox(QB_tools* output,QB_tools* input,double lattice[3][3],int hall_number,double symprec);

//redifine orientation
extern void QSPG_redefine(QB_tools* output,QB_tools* input,int rot[3][3],double symprec);
extern void SPG_wrap_box(SPG_tools* spg,double symprec);
extern void QSPG_reshapebox(QB_tools* output,QB_tools* input,int a1,int x1,int a2,int x2);
extern QSPG_data QSPG_database[531];
extern int QSPG_spacegroup2hallnum[231];
extern int QSPG_spacegroup_choice[231];