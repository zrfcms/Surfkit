/*===============================================================
*   Copyright[c] 2015-2016, Z. R. Liu and R. F. Zhang
*
*   This file is part of the
*   AACSD - Atomic Analyzer for Crystal Structure and Defect
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*================================================================*/
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<vector>
#ifdef SPaMD_viewer
#include <FL/Fl.H>
#include "FL/Fl_Value_Input.H"
#include "FL/Fl_Choice.H"
#define QB_fprintf fprintf
#define QB_fopen fl_fopen
#define QB_fscanf fscanf
#define QB_printf printf
#define QB_fclose fclose
#define QB_fgets fgets
#else
#define QB_fprintf fprintf
#define QB_fopen fopen
#define QB_fscanf fscanf
#define QB_printf printf
#define QB_fclose fclose
#define QB_fgets fgets
/*
#define QB_fprintf SPaMD_fprintf
#define QB_fopen  SPaMD_fopen

#define QB_fclose SPaMD_fclose
#define QB_fgets  SPaMD_fgets */
#endif
typedef struct{
    //atom (outside program) ID
    int     id;
    //atom type
    int     type;
    //position
    double  x;
    double  y;
    double  z;
}atomdata;//save atom data

typedef struct{
    double  x;
    double  y;
    double  z;
}QB_vector;//save atom data

typedef struct
{
	char name[32];
}VASP_ELENAME;
typedef struct
{
    char        **exch;
    double      **exdat;
    int         exchannel;
    atomdata    *atom;
    double      boundx,boundy,boundz;
	double      mat[3][3];
	double      zerox,zeroy,zeroz;
	int 		box;//whether use xy
    double      startx,starty,startz;
    int         TypeNumber;
    int	        TotalNumber;
    double	  **par_exdouble;
    int 	  **par_exint;
    int	        num_exint;
    char	  **name_exint;
    int 	    num_exdouble;
    char	  **name_exdouble;
    int         px,py,pz;//1 = periodic boundary condition
    ///------------------------------------------------------
    ///boost net system
    ///------------------------------------------------------
    //lenth of block in every direction
    double      tbx,tby,tbz;
    //number of block in every direction
    int         xlot,ylot,zlot;
    //ID of atom in block
    int         ****network_ID;
    //number of atom in each block
    int         ***network_num;
    int         network_flag;
    int     	list_flag;
    int     	MCN;
    int         extra_last;    
	struct
    {
        //how many type name in the list
        int     n;
        VASP_ELENAME* t;//type name
    }ele;//save neighbor list
    struct
    {
        //how many neighbor in the list
        int     num;
        //position ID (containing PBC) of neighbor
        int     *id;
        double  *dis;
    }neb;//save neighbor list
}QB_tools;

extern void	    	QB_init(QB_tools *QB); 
//read LAMMPS input or dump data
extern int          QB_read_lmp(QB_tools *QB,const char name[]);
extern int          QB_read_lmc(QB_tools *QB,const char name[]);
extern int          QB_read_cif(QB_tools *QB,const char name[]);
extern int          QB_read_vasp(QB_tools *QB,const char name[]);
extern void 		QB_free_exchannel(QB_tools* QB);
#define QB_F_UNKNOWN 0
#define QB_S_LMP 1
#define QB_F_LMPC 2
#define QB_S_LMC 3
#define QB_F_LMCC 4
#define QB_S_VASP 5
#define QB_F_VASPC 6
#define QB_S_CIF 7
extern int          QB_read_file(QB_tools *QB,const char name[50]);
//create single crystal
extern double 	    QB_create_fcc(QB_tools *QB,double a,double M11,double M12,double M13,double M21,double M22,double M23,double M31,double M32,double M33);
extern double       QB_create_bcc(QB_tools *QB,double a,double M11,double M12,double M13,double M21,double M22,double M23,double M31,double M32,double M33);
extern double       QB_create_hcp(QB_tools *QB,double a,double c,double M11,double M12,double M13,double M21,double M22,double M23,double M31,double M32,double M33);
extern double       QB_create_any(QB_tools *QB,double a,double M11,double M12,double M13,double M21,double M22,double M23,double M31,double M32,double M33,int n,double single_list[][4]);
//free extra channel
extern void         QB_free_exchannel(QB_tools *QB);
extern double 	    QB_get_data(QB_tools *QB,int i,const char *name);
extern void 	    QB_save_data(QB_tools *QB,int i,const char *name,double x);
///------------------------------------------------------
///data control system
///------------------------------------------------------
extern void         QB_duplicate(QB_tools *QB,int nx,int ny,int nz);
extern void 	    QB_recover_atom(QB_tools *QB,int i);
extern void 	    QB_delete_atom(QB_tools *QB,int i);
extern void 		QB_affine_delete(QB_tools *QB);
extern void	        QB_create_atom(QB_tools *QB,int type,double x,double y,double z);
extern void         QB_move_box(QB_tools *QB,double movement_x,double movement_y,double movement_z);
extern void         QB_move_atom(QB_tools *QB,int i,double movement_x,double movement_y,double movement_z);
extern void 	    QB_rotate(QB_tools *QB,int i,QB_vector center,QB_vector vv, float angle);
extern void 	    QB_mirror(QB_tools *QB,int i,QB_vector p0, double dis);
extern void 	    QB_swichxz(QB_tools *QB);
extern void 	    QB_swichxy(QB_tools *QB);
extern void         QB_swichyz(QB_tools *QB);
extern void 	    QB_wrap(QB_tools *QB,int i);//wrap an atom back to box;
extern void 		QB_wrap_vec(QB_tools *QB,double *v_x,double *v_y,double *v_z);//wrap any point back to box;
extern int 	        QB_inbox(QB_tools *QB,int i);
extern void 	    QB_swap(QB_tools *QB,int i,int j);
extern void         QB_sort_atom(QB_tools *QB,int deal(QB_tools *QB,int a,int b));
extern void 		QB_combine(QB_tools *QB1,QB_tools *QB2);

///------------------------------------------------------
///data saving system
///------------------------------------------------------

extern void   	    QB_free_atom(QB_tools *QB);
extern void 		QB_free_all(QB_tools *QB);
///------------------------------------------------------
///adding data system
///------------------------------------------------------
extern void	        QB_add_exint(QB_tools *QB,const char *name);	
extern void 	    QB_add_exdouble(QB_tools *QB,const char *name);	

///------------------------------------------------------
///data slot(faster data) system
///------------------------------------------------------
extern int          QB_slot_sum(QB_tools *QB);
extern int          QB_slot_get(QB_tools *QB,const char *name);
extern int          QB_slot_name(QB_tools *QB,int i,char *name);
extern double       QB_slot_data(QB_tools *QB,int slot,int i);
extern void 	    QB_slot_save(QB_tools *QB,int slot,int i,double x);
///update and free
extern void 	    QB_update_extra(QB_tools *QB);//update extra data so that we can keep extra data's number equal with atom 
extern void         QB_free_extra(QB_tools *QB);//free extra data 
//lmc dump system
extern void	        QB_dump_lmc(QB_tools *QB,const char *name);
extern void	        QB_dump_lmp(QB_tools *QB,const char *name);
extern void	        QB_dump_cif(QB_tools *QB,const char *name);
#ifdef SPaMD_model
extern int  	    QB_SPaMD_upload(QB_tools *QB);
#endif

///------------------------------------------------------
///predioc boundary condition modifer
///------------------------------------------------------
extern void        QB_pbc(QB_tools *QB,double rmax); 
extern void        QB_pbc_clean(QB_tools *QB);
//init boost network
extern void        QB_network_init(QB_tools *QB,double cutoff);
extern void 	   QB_network_update(QB_tools *QB,double cutoff,int i);
//free networf
extern void        QB_network_free(QB_tools *QB);

///------------------------------------------------------
///chain table system
///------------------------------------------------------
extern void        QB_list_build(QB_tools *QB,int i,double cutoff);//get chain table based on specific atom
extern void        QB_nearest_atom(QB_tools *QB,double r_x,double r_y,double r_z,double cutoff);//get chain table based on specific position
extern void        QB_pbclist_build(QB_tools *QB,int i,double cutoff);
extern void        QB_list_clear(QB_tools *QB);
extern void 	   QB_copy_system(QB_tools *output,QB_tools *input);
extern void 	   QB_copy_box(QB_tools *output,QB_tools *input);
extern void 	   QB_combine_system(QB_tools *output,QB_tools *input);
///math functions
//return x*x+y*y+z*z
extern double  	   QB_squaredlenths(double x,double y,double z);
///check data change charactor in to float number
extern double      QB_checkdat(char c[50]);
extern int 	       QB_checkdat_err;//if there is error when QB_checkdat
extern int         QB_checkchar(char c[50]);
//min and max function
extern int         QB_Min(int x,int y);
extern int         QB_Max(int x,int y);

extern void 	   QB_save_elename(QB_tools *QB,int i,char* name);
extern void 	   QB_load_elename(QB_tools *QB,int i,char* name);
extern int 		   QB_get_elename(QB_tools *QB,char*name);
extern void 	   QB_free_elename(QB_tools *QB);

extern int 		   QB_read_vasp(QB_tools *QB,const char name[50]);
extern void		   QB_dump_vasp_Cartesian(QB_tools *QB,const char name[50]);
extern void		   QB_dump_vasp_Direct(QB_tools *QB,const char name[50]);

extern void 	   QB_mat_update(QB_tools *QB);
extern void 	   QB_bound_update(QB_tools *QB);

#ifdef SPaMD_viewer
//region system
extern int         DEAL_ALL(int id,atomdata pt);
extern int         DEAL_NONE(int id,atomdata pt);
extern double      VOL_BOX();
extern double      VOL_BLANK();
extern int 		   QB_region_add(int deal(int id,atomdata pt),double vol);
extern int         QB_region_extract(int atom_id,atomdata pt,int region_id);
extern int         QB_region_and(int a,int b,double vol);
extern int         QB_region_or(int a,int b,double vol);
extern int         QB_region_comple(int a,int b,double vol);
extern int         QB_region_neither(int a,int b,double vol);
extern int         QB_region_not(int a,double vol);
extern int         QB_region_cuboid(double minx,double miny,double minz,double maxx,double maxy,double maxz);
extern int         QB_region_sphere(double cenx,double ceny,double cenz,double r);
extern int         QB_region_cylinder(double cenx,double ceny,double cenz,double r,double h,int axis);
extern int         QB_region_plane(double x,double y,double z,double px,double py,double pz,double thickness);
extern int         QB_region_cone(double x,double y,double z,double r,double h,int dir);
extern int 		   QB_region_rotate(int region_id,double cx,double cy,double cz,double lx,double ly,double lz,double angle);
extern int         QB_region_move(int region_id,int unit,double x,double y,double z);
extern int 		   QB_region_polyhedron(int mode,int num,
						double x,double y,double z,
						double dx1,double dy1,double dz1,
						double dx2,double dy2,double dz2,
						double dx3,double dy3,double dz3,
						double distance1,double distance2,double distance3);
extern int         QB_region_prism(double x,double y,double z,double distance,double h,int n,int axis);
extern int         QB_region_periodic(int region_id);
extern int         QB_region_expression(const char slot[1024],double min,double max);
extern void        QB_region_set_invest_system(QB_tools*QB);
extern int		   QB_region_deform(int region_id,double x,double y,double z,double xx,double xy,double xz,double yx,double yy,double yz,double zx,double zy,double zz);
extern void		   QB_region_free();
extern int 	       QB_region_inbox();

//type def for region system
typedef struct{
    int (*deal)(int id,atomdata pt);
    double vol;
}REGION_CTRL;
extern REGION_CTRL* QB_region_control;
#define QB_REGION_AND     1
#define QB_REGION_OR      2
#define QB_REGION_COMPLE  3
#define QB_REGION_NEITHER 4
#define QB_REGION_NOT     5
typedef struct{
    int id[2];
    int mode;
}REGION_MATH;
extern REGION_MATH* QB_region_math;
typedef struct 
{
    double minx;
    double miny;
    double minz;
    double maxx;
    double maxy;
    double maxz;
}REGION_CUBOID;
extern REGION_CUBOID* QB_region_cuboiddata;
typedef struct 
{
    double R;
    double X;
    double Y;
    double Z;
}REGION_SPHERE;
extern REGION_SPHERE* QB_region_spheredata;
typedef struct 
{
    double R;
    double X;
    double Y;
    double Z;
    double H;
    int axis;
}REGION_CYLINDER;
extern REGION_CYLINDER* QB_region_cylinderdata;
typedef struct 
{
    double X;
    double Y;
    double Z;
    double B;
	double T;
}REGION_PLANE;
extern REGION_PLANE* QB_region_planedata;
typedef struct 
{
    double X;
    double Y;
    double Z;
    double H;
    double R;
    int  DIR;
}REGION_CONE;
extern REGION_CONE* QB_region_conedata;
typedef struct 
{
    double C_X;
    double C_Y;
    double C_Z;
    double L_X;
    double L_Y;
    double L_Z;
    double cos_angle;
    double sin_angle;
    int    id;
}REGION_ROTATE;
extern REGION_ROTATE* QB_region_rotatedata;
typedef struct 
{
    double X;
    double Y;
    double Z;
    int    id;
	int    unit;
}REGION_MOVE;
extern REGION_MOVE* QB_region_movedata;
typedef struct 
{
    double X;
    double Y;
    double Z;
    double mat[3][3];
    int    id;
}REGION_DEFORM;
extern REGION_DEFORM* QB_region_deformdata;
typedef struct 
{
    int    id;
}REGION_PERIODIC;
extern REGION_PERIODIC* QB_region_periodicdata;

typedef struct 
{
    char   slot[1024];
	double max;
	double min;
}REGION_EXPRESSION;
extern REGION_EXPRESSION* QB_region_expressiondata;

typedef struct 
{
    int    mode;
	int    num;
	double x,y,z;
	double dir_x[3],dir_y[3],dir_z[3];
	double distance[3];
	Fl_Value_Input *dx_ui[3];
	Fl_Value_Input *dy_ui[3];
	Fl_Value_Input *dz_ui[3];
	Fl_Value_Input *d_ui[3];
	Fl_Choice *num_ui;
}REGION_POLYHEDRON;
extern REGION_POLYHEDRON* QB_region_polyhedrondata;

typedef struct 
{
    int    mode;//X Y Z
	double x,y,z;
	double distance;
	double h;
	int n;
}REGION_PRISM;
extern REGION_PRISM* QB_region_prismdata;
#endif