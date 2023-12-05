#include"./../QB/QB.h"
#include"QSPG.h"
#include "spglib/cell.h"
#include "spglib/mathfunc.h"
#include "spglib/symmetry.h"
#include "spglib/site_symmetry.h"
#include "spglib/sitesym_database.h"
extern void SPG_add_atom(SPG_tools* spg,int type,double x,double y,double z);
static int SPG_get_Wyckoff_notation(char site_sym_symbol[7],
									const double position[3],
									const int sym_size,
									int (*sym_rot)[3][3],
									const double (*sym_trans)[3],
									const int ref_multiplicity,
									SPGCONST double bravais_lattice[3][3],
									const int hall_number,
									const double symprec);
static int SPG_get_equivalent_atoms(int symmetry_num,int (*rot)[3][3],double (*tra)[3],double x,double y,double z,double symprec);							
void SPG_pick_site_symmetry( int* multiplicity,
							 int* wyckoff_letter,
							 char site_sym_symbol[7],
							 const double position[3],
							 SPGCONST double bravais_lattice[3][3],
							 const int hall_number,
							 const double symprec)
{
	
	int rotations[192][3][3];
    double translations[192][3];
	int size=spg_get_symmetry_from_database(rotations,translations,hall_number);
	*multiplicity=SPG_get_equivalent_atoms(size,rotations,translations,position[0],position[1],position[2],symprec);
	*wyckoff_letter=SPG_get_Wyckoff_notation(site_sym_symbol,
											 position,
											 size,
											 rotations,
											 translations,
											 *multiplicity,
											 bravais_lattice,
											 hall_number,
											 symprec);
}
static int SPG_get_equivalent_atoms(int symmetry_num,int (*rot)[3][3],double (*tra)[3],double x,double y,double z,double symprec)
{
	SPG_tools spg;
	spg.num=0;
	
	int i,j;
	double symprec2=symprec*symprec;
	//get particle id from hallnumber
	
	int failed_flag;
	double vec[3];
	double dsp[3];
  
	for(i=0;i<symmetry_num;i++)
	{
		failed_flag=0;
		vec[0]=rot[i][0][0]*x+rot[i][0][1]*y+rot[i][0][2]*z+tra[i][0];
		vec[1]=rot[i][1][0]*x+rot[i][1][1]*y+rot[i][1][2]*z+tra[i][1];
		vec[2]=rot[i][2][0]*x+rot[i][2][1]*y+rot[i][2][2]*z+tra[i][2];
		for(j=0;j<spg.num;j++)
		{
			dsp[0]=vec[0]-spg.pos[j][0];
			dsp[1]=vec[1]-spg.pos[j][1];
			dsp[2]=vec[2]-spg.pos[j][2];
			while(dsp[0]<-0.5)dsp[0]+=1;
			while(dsp[1]<-0.5)dsp[1]+=1;
			while(dsp[2]<-0.5)dsp[2]+=1;
			while(dsp[0]>0.5)dsp[0]-=1;
			while(dsp[1]>0.5)dsp[1]-=1;
			while(dsp[2]>0.5)dsp[2]-=1;
			if((dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2])<symprec2)
			{
				failed_flag=1;
				break;
			}
		}
		if(!failed_flag)
			SPG_add_atom(&spg,1,vec[0],vec[1],vec[2]);
	}
	int num=spg.num;
	SPG_free(&spg);
	return num;
}
static int SPG_get_Wyckoff_notation(char site_sym_symbol[7],
									const double position[3],
									const int sym_size,
									int (*sym_rot)[3][3],
									const double (*sym_trans)[3],
									const int ref_multiplicity,
									SPGCONST double bravais_lattice[3][3],
									const int hall_number,
									const double symprec)
{
  int i, j, k, l, num_sitesym, multiplicity, wyckoff_letter;
  int indices_wyc[2];
  int rot[3][3];
  double trans[3], orbit[3];
  VecDBL *pos_rot;

  wyckoff_letter = -1;
  pos_rot = NULL;

  if ((pos_rot = mat_alloc_VecDBL(sym_size)) == NULL) {
    return -1;
  }

  for (i = 0; i < sym_size; i++) {
    mat_multiply_matrix_vector_id3(pos_rot->vec[i], sym_rot[i], position);
    for (j = 0; j < 3; j++) {
      pos_rot->vec[i][j] += sym_trans[i][j];
    }
  }

  ssmdb_get_wyckoff_indices(indices_wyc, hall_number);
  for (i = 0; i < indices_wyc[1]; i++) {
    /* (rot, trans) gives the first element of each Wyckoff position */
    /* of the 'Coordinates' in ITA */
    /* Example: (x,1/4,1/2)      */
    /* rot         trans         */
    /* [1, 0, 0]   [0, 1/4, 1/2] */
    /* [0, 0, 0]                 */
    /* [0, 0, 0]                 */
    multiplicity = ssmdb_get_coordinate(rot, trans, i + indices_wyc[0]);

    /* Effectively this iteration runs over all 'Coordinates' of each */
    /* Wyckoff position, i.e., works as looking for the first element. */
    for (j = 0; j < pos_rot->size; j++) {
      num_sitesym = 0;
      for (k = 0; k < pos_rot->size; k++) {
        if (cel_is_overlap(pos_rot->vec[j],
                           pos_rot->vec[k],
                           bravais_lattice,
                           symprec)) {
          mat_multiply_matrix_vector_id3(orbit, rot, pos_rot->vec[k]);
          for (l = 0; l < 3; l++) {
            orbit[l] += trans[l];
          }
          if (cel_is_overlap(pos_rot->vec[k],
                             orbit,
                             bravais_lattice,
                             symprec)) {
            num_sitesym++;
          }
        }
      }

      /* Consistency check */
      /* 1) num_sym == num_sitesym * m */
      /* 2) num_equiv_atoms in conventional cell == m */
      if ((num_sitesym * multiplicity == sym_size) &&
          (multiplicity == ref_multiplicity)) {
        /* Database is made reversed order, e.g., gfedcba. */
        /* wyckoff is set 0 1 2 3 4... for a b c d e..., respectively. */
        wyckoff_letter = indices_wyc[1] - i - 1;
        ssmdb_get_site_symmetry_symbol(site_sym_symbol, indices_wyc[0] + i);
        goto end;
      }
    }
  }

 end:
  mat_free_VecDBL(pos_rot);
  pos_rot = NULL;
  return wyckoff_letter;
}