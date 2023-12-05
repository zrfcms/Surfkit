#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"
#include "./../QSPG/spglib/cell.h"
#include "./../QSPG/spglib/mathfunc.h"
#include "./../QSPG/spglib/symmetry.h"
#include "./../QSPG/spglib/site_symmetry.h"
#include "./../QSPG/spglib/sitesym_database.h"
#define MAX_SITE_3D 4096
#define MAX_SITE_2D 256
static int SK_get_Wyckoff_notation(char site_sym_symbol[7],
									const double position[3],
									const int sym_size,
									int (*sym_rot)[3][3],
									const double (*sym_trans)[3],
									const int ref_multiplicity,
									SPGCONST double bravais_lattice[3][3],
									const int hall_number,
									const double symprec);
void SK_malloc_sitelist(int n,SK_site_list* list);
//loop space by 1/24 bravais lattice to search possible solute site
void SK_get_sitelist_3D(SK_site_list *list,
						const int sym_size,
						int (*rot)[3][3],
						const double (*trans)[3],
						SPGCONST double bravais_lattice[3][3],
						const int hall_number,
						double symprec)
{
	int i,j,k,n;
	char temp_symbol[7];
	double posible_pos[16]={
		  0 ,1.0/12.0,1.0/8.0,
	1.0/6.0 ,1.0/4.0 ,1.0/3.0,
	3.0/8.0 ,5.0/12.0,1.0/2.0,
	7.0/12.0,5.0/8.0 ,2.0/3.0,
	3.0/4.0 ,5.0/6.0 ,7.0/8.0,
	11.0/12.0};
	SK_malloc_sitelist(MAX_SITE_3D,list);
	n=0;
	for(i=0;i<16;i++)
	for(j=0;j<16;j++)
	for(k=0;k<16;k++)
	{
		list->pos[n][0]=posible_pos[i];
		list->pos[n][1]=posible_pos[j];
		list->pos[n][2]=posible_pos[k];
		list->equivalent_id[n]=-1;
		n++;
	}
	int equivalent_id=0;
	for(n=0;n<list->num;n++)
	{
		if(list->equivalent_id[n]<0)
		{
			double symprec2=symprec*symprec;
			//get particle id from hallnumber
			int failed_flag;
			double vec[3];
			double dsp[3];
			int equivalent_num=1;//skip self
			list->equivalent_id[n]=equivalent_id;
			for(i=0;i<sym_size;i++)
			{
				failed_flag=0;
				vec[0]=rot[i][0][0]*list->pos[n][0]+rot[i][0][1]*list->pos[n][1]+rot[i][0][2]*list->pos[n][2]+trans[i][0];
				vec[1]=rot[i][1][0]*list->pos[n][0]+rot[i][1][1]*list->pos[n][1]+rot[i][1][2]*list->pos[n][2]+trans[i][1];
				vec[2]=rot[i][2][0]*list->pos[n][0]+rot[i][2][1]*list->pos[n][1]+rot[i][2][2]*list->pos[n][2]+trans[i][2];
				while(vec[0]<-symprec)vec[0]+=1;
				while(vec[1]<-symprec)vec[1]+=1;
				while(vec[2]<-symprec)vec[2]+=1;
				while(vec[0]>1-symprec)vec[0]-=1;
				while(vec[1]>1-symprec)vec[1]-=1;
				while(vec[2]>1-symprec)vec[2]-=1;
				for(j=n;j<list->num;j++)
				{
					dsp[0]=vec[0]-list->pos[j][0];
					dsp[1]=vec[1]-list->pos[j][1];
					dsp[2]=vec[2]-list->pos[j][2];
					while(dsp[0]<-0.5)dsp[0]+=1;
					while(dsp[1]<-0.5)dsp[1]+=1;
					while(dsp[2]<-0.5)dsp[2]+=1;
					while(dsp[0]>0.5)dsp[0]-=1;
					while(dsp[1]>0.5)dsp[1]-=1;
					while(dsp[2]>0.5)dsp[2]-=1;
					if((dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2])<symprec2)
					{
						if(list->equivalent_id[j]!=equivalent_id)
						{
							equivalent_num++;
							list->equivalent_id[j]=equivalent_id;
						}
						failed_flag=1;
					}
				}
				if(!failed_flag)
				{
					equivalent_num++;
					SK_malloc_sitelist(list->num+1,list);
					list->equivalent_id[list->num-1]=equivalent_id;
					list->pos[list->num-1][0]=vec[0];
					list->pos[list->num-1][1]=vec[1];
					list->pos[list->num-1][2]=vec[2];
				}
			}
			list->wyckoff[n]=SK_get_Wyckoff_notation(temp_symbol,
									list->pos[n],
									sym_size,
									rot,
									trans,
									equivalent_num,
									bravais_lattice,
									hall_number,
									symprec);
			list->equivalent_num[n]=equivalent_num;
			for(j=n+1;j<list->num;j++)
			{
				if(list->equivalent_id[n]==list->equivalent_id[j])//copy wyckoff data
				{
					list->wyckoff[j]=list->wyckoff[n];
					list->equivalent_num[j]=equivalent_num;
				}
			}
			equivalent_id++;
		}
	}
}

//loop space by 1/24 bravais lattice to search possible solute site
void SK_get_sitelist_2D(SK_site_list *list,
						int v1,int v2,
						const int sym_size,
						int (*rot)[3][3],
						const double (*trans)[3],
						SPGCONST double bravais_lattice[3][3],
						const int hall_number,
						double symprec)
{
	int i,j,k,n;
	char temp_symbol[7];
	double posible_pos[16]={
		  0 ,1.0/12.0,1.0/8.0,
	1.0/6.0 ,1.0/4.0 ,1.0/3.0,
	3.0/8.0 ,5.0/12.0,1.0/2.0,
	7.0/12.0,5.0/8.0 ,2.0/3.0,
	3.0/4.0 ,5.0/6.0 ,7.0/8.0,
	11.0/12.0};
	SK_malloc_sitelist(MAX_SITE_2D,list);
	n=0;
	
	int v3=0;
	for(i=0;i<3;i++)if(i!=v1&&i!=v2)v3=i;
	
	for(i=0;i<16;i++)
	for(j=0;j<16;j++)
	{
		list->pos[n][v1]=posible_pos[i];
		list->pos[n][v2]=posible_pos[j];
		list->pos[n][v3]=0;
		list->equivalent_id[n]=-1;
		n++;
	}
	int equivalent_id=0;
	for(n=0;n<list->num;n++)
	{
		if(list->equivalent_id[n]<0)
		{
			double symprec2=symprec*symprec;
			//get particle id from hallnumber
			int failed_flag;
			double vec[3];
			double dsp[3];
			int equivalent_num=1;//skip self
			list->equivalent_id[n]=equivalent_id;
			for(i=0;i<sym_size;i++)
			{
				failed_flag=0;
				vec[0]=rot[i][0][0]*list->pos[n][0]+rot[i][0][1]*list->pos[n][1]+rot[i][0][2]*list->pos[n][2]+trans[i][0];
				vec[1]=rot[i][1][0]*list->pos[n][0]+rot[i][1][1]*list->pos[n][1]+rot[i][1][2]*list->pos[n][2]+trans[i][1];
				vec[2]=rot[i][2][0]*list->pos[n][0]+rot[i][2][1]*list->pos[n][1]+rot[i][2][2]*list->pos[n][2]+trans[i][2];
				while(vec[0]<-symprec)vec[0]+=1;
				while(vec[1]<-symprec)vec[1]+=1;
				while(vec[2]<-symprec)vec[2]+=1;
				while(vec[0]>1-symprec)vec[0]-=1;
				while(vec[1]>1-symprec)vec[1]-=1;
				while(vec[2]>1-symprec)vec[2]-=1;
				for(j=n;j<list->num;j++)
				{
					dsp[0]=vec[0]-list->pos[j][0];
					dsp[1]=vec[1]-list->pos[j][1];
					dsp[2]=vec[2]-list->pos[j][2];
					while(dsp[0]<-0.5)dsp[0]+=1;
					while(dsp[1]<-0.5)dsp[1]+=1;
					while(dsp[2]<-0.5)dsp[2]+=1;
					while(dsp[0]>0.5)dsp[0]-=1;
					while(dsp[1]>0.5)dsp[1]-=1;
					while(dsp[2]>0.5)dsp[2]-=1;
					if((dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2])<symprec2)
					{
						if(list->equivalent_id[j]!=equivalent_id)
						{
							equivalent_num++;
							list->equivalent_id[j]=equivalent_id;
						}
						failed_flag=1;
					}
				}
				if(!failed_flag)
				{
					equivalent_num++;
					SK_malloc_sitelist(list->num+1,list);
					list->equivalent_id[list->num-1]=equivalent_id;
					list->pos[list->num-1][0]=vec[0];
					list->pos[list->num-1][1]=vec[1];
					list->pos[list->num-1][2]=vec[2];
				}
			}
			list->wyckoff[n]=SK_get_Wyckoff_notation(temp_symbol,
									list->pos[n],
									sym_size,
									rot,
									trans,
									equivalent_num,
									bravais_lattice,
									hall_number,
									symprec);
			list->equivalent_num[n]=equivalent_num;
			for(j=n+1;j<list->num;j++)
			{
				if(list->equivalent_id[n]==list->equivalent_id[j])//copy wyckoff data
				{
					list->wyckoff[j]=list->wyckoff[n];
					list->equivalent_num[j]=equivalent_num;
				}
			}
			equivalent_id++;
		}
	}
}

void SK_malloc_sitelist(int n,SK_site_list* list)
{
	if(n<=0)
		return;
	if(list->num==0)
	{
		list->pos=(double(*)[3])malloc(3*n*sizeof(double));
		list->wyckoff=(int*)malloc(n*sizeof(int));
		list->equivalent_num=(int*)malloc(n*sizeof(int));
		list->equivalent_id=(int*)malloc(n*sizeof(int));
	}
	else
	{
		list->pos=(double(*)[3])realloc(list->pos,3*n*sizeof(double));
		list->wyckoff=(int*)realloc(list->wyckoff,n*sizeof(int));
		list->equivalent_num=(int*)realloc(list->equivalent_num,n*sizeof(int));
		list->equivalent_id=(int*)realloc(list->equivalent_id,n*sizeof(int));
	}
	list->num=n;
}

static int SK_get_Wyckoff_notation(char site_sym_symbol[7],
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