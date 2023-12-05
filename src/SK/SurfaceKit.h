//get orientation according to miller indices
extern void SK_auto_Frot(int c[3],int F_rot[3][3]);
//redefine current box
extern void SK_redefine(SPG_tools* spg,int F_rot[3][3],double symprec);
//get minimum area of selected plane
//defined by in-plane vector
//0 -> a
//1 -> b
//2 -> c
extern void SK_SPG_minimize_face(SPG_tools* spg,int v1,int v2,double symprec);
//get minimum height of selected cylinder
//defined by in-plane vector
//0 -> a
//1 -> b
//2 -> c
extern void SK_SPG_minimize_cylinder(SPG_tools* spg,int v1,int v2,double symprec);
//get 3D repetition structure of selected cylinder
//defined by in-plane vector
//0 -> a
//1 -> b
//2 -> c
extern int  SK_SPG_pick_3Dlattice(SPG_tools* spg,int v1,int v2,double symprec);
//rot matrix is dump out to get new orientation
extern int  SPG_delaunay_reduce_surface(double red_lattice[3][3],
									    int    F_rot[3][3],
									    double lattice[3][3],
									    double symprec);

//check whether particle exist in selected position
//1=yes 0=no
extern int SK_SPG_overloop_vec(SPG_tools* spg,int type,double v_x,double v_y,double v_z,double symprec);

//save analyzed equivalent surface
#define SK_ATOMIC_LAYER 0
#define SK_POLAR_CENTER 1
typedef struct
{
	char name[1024];
	char integration[1024];
	int type;
	double height;
}SK_layer_list;										 
//save analyzed high symmetric points
typedef struct 
{
	int num;
	double (*pos)[3];
	int *wyckoff;
	int *equivalent_id;
	int *equivalent_num;
}SK_site_list;

//save 
typedef struct
{
   double dis;
   double theta;
   double dtheta;
   int id_guest;
   int vec[2];
}SK_crystal_face;

//save data of in plane vector pairs
typedef struct
{
   double theta1;
   double theta2;
   double sigma1;
   double sigma2;
   int trans1[3][3];
   int trans2[3][3];
   double mismatch[5];
   int cylinder_flag;
}SK_scaned_data;

//get possible solute site according to wyckoff symbol									   
extern void SK_get_sitelist_3D(SK_site_list *list,
							   const int sym_size,
							   int (*sym_rot)[3][3],
							   const double (*sym_trans)[3],
							   SPGCONST double bravais_lattice[3][3],
							   const int hall_number,
							   double symprec);
							   
//get possible solute site according to wyckoff symbol							   
extern void SK_get_sitelist_2D(SK_site_list *list,
							   int v1,int v2,
							   const int sym_size,
							   int (*sym_rot)[3][3],
							   const double (*sym_trans)[3],
							   SPGCONST double bravais_lattice[3][3],
							   const int hall_number,
							   double symprec);


extern void SK_cross_product(double out[3],double in1[3],double in2[3]);
extern double SK_min_height(double C_mat[3][3]);
extern void SPG_displace(SPG_tools*spg,double dsp[3]);
extern void SK_rotate_matrix(double F_mat[3][3],double vv[3], double angle);
extern void SPG_copy(SPG_tools*input,SPG_tools*output);
extern void SK_copy_matrix(double output[3][3],double input[3][3]);
extern int  SK_getrp(int a, int b);

extern void SK_mismatch_set(double max_d,double min_d,double max_box,double limitdt,double min_gamma,
							double max_gamma,double max_mismatch,double max_dangle);
extern int  SK_mismatch_scan(double C_mat1[3][3],double C_mat2[3][3],SK_scaned_data** rotate_list);
extern int  SK_mismatch_scan_half(double C_mat1[3][3],double C_mat2[3][3],double A_vec1[2],double A_vec2[2],SK_scaned_data** rotate_list);
extern void SK_cylinder_scan(double C_mat1[3][3],double C_mat2[3][3],SK_scaned_data** rotate_list,
							 int selected,double thickness1,double thickness2);
							 
extern void SK_layerdata_save(char filename[1024],int  list_n,SK_scaned_data*  list);		
extern void SK_layerdata_load(char filename[1024],int *list_n,SK_scaned_data** list);		

extern void SK_transposition_matrix(double input[3][3],double output[3][3]);
extern void SK_matrix_multiplication(double input1[3][3],double input2[3][3],double output[3][3]);
extern void SK_inverse_matrix(double input[3][3],double output[3][3]);

extern void SK_matrix_left_multiplication(double input1[3][3],double input2[3],double output[3]);

extern void SK_polar_center(int *list_n,SK_layer_list **list,double min_height,double max_height,SPG_tools* spg,double symprec);
extern void SK_atomic_layers(int *list_n,SK_layer_list **list,double min_height,double max_height,QB_tools *QB,double symprec);
extern void SK_layer_dump(char output_name[1024],int list_n,SK_layer_list *list);
extern void SK_layer_sort(int list_n,SK_layer_list*list);

extern void SK_matrix_print(double input[3][3]);

extern double SK_rigid_average_r(QB_tools *QB,double r,double min,double max);
extern void SK_rigid_alias(QB_tools *QB,double h,double r,int max_step,double symprec);

extern void Gaussian_lattice_reduction(double vec1[3],double vec2[3]);
extern void Right_hand_reduction(double mat[3][3]);