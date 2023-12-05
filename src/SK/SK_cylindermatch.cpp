#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"
extern double SK_limitd;
extern double SK_limitbox;
extern double SK_limitdt;
extern double SK_mismatch;
extern double SK_mingamma;
extern double SK_maxgamma;
extern void   SK_rotate_Cmatrix(double C_mat[3][3],double vv[3], double angle);
void SK_cylinder_scan(double C_mat1[3][3],double C_mat2[3][3],SK_scaned_data** rotate_list,
					  int selected,double thickness1,double thickness2)
{
	int i;
	double T_mat1[3][3],T_mat2[3][3];
	SK_copy_matrix(T_mat1,C_mat1);
	SK_copy_matrix(T_mat2,C_mat2);
	double R_axis[3]={0,0,1};
	SK_rotate_Cmatrix(T_mat1,R_axis,-(*rotate_list)[selected].theta1);
	SK_rotate_Cmatrix(T_mat2,R_axis,-(*rotate_list)[selected].theta2);
	
	int c1=(int)(thickness1/T_mat1[2][2])+1;
	int c2=(int)(thickness2/T_mat2[2][2])+1;
	(*rotate_list)[selected].trans1[2][2]=c1;
	(*rotate_list)[selected].trans2[2][2]=c2;
	
	for(i=0;i<3;i++)
	{
		T_mat1[2][i]*=c1;
		T_mat2[2][i]*=c2;
	}
	
	int	   i0,j0,k0,l0;
	if((T_mat1[1][0]*T_mat1[0][1]-T_mat1[0][0]*T_mat1[1][1])!=0)
	{
		j0=(int)((-T_mat1[2][0]*T_mat1[0][1]+T_mat1[2][1]*T_mat1[0][0])/(T_mat1[1][0]*T_mat1[0][1]-T_mat1[0][0]*T_mat1[1][1]));
		i0=(int)(( T_mat1[2][0]*T_mat1[1][1]-T_mat1[2][1]*T_mat1[1][0])/(T_mat1[1][0]*T_mat1[0][1]-T_mat1[0][0]*T_mat1[1][1]));
	}
	else
	{
		i0=j0=0;
	}
	
	if((T_mat2[1][0]*T_mat2[0][1]-T_mat2[0][0]*T_mat2[1][1])!=0)
	{
		l0=(int)((-T_mat2[2][0]*T_mat2[0][1]+T_mat2[2][1]*T_mat2[0][0])/(T_mat2[1][0]*T_mat2[0][1]-T_mat2[0][0]*T_mat2[1][1]));
		k0=(int)(( T_mat2[2][0]*T_mat2[1][1]-T_mat2[2][1]*T_mat2[1][0])/(T_mat2[1][0]*T_mat2[0][1]-T_mat2[0][0]*T_mat2[1][1]));
	}
	else
	{
		k0=l0=0;
	}
	(*rotate_list)[selected].trans1[0][2]=i0;
	(*rotate_list)[selected].trans1[1][2]=j0;
	(*rotate_list)[selected].trans2[0][2]=k0;
	(*rotate_list)[selected].trans2[1][2]=l0;
}

/*
void SK_cylinder_scan(double C_mat1[3][3],double C_mat2[3][3],SK_scaned_data** rotate_list,
					  int selected,double thickness1,double thickness2)
{
	int i,j,k,l;
	
	double T_mat1[3][3],T_mat2[3][3];
	SK_copy_matrix(T_mat1,C_mat1);
	SK_copy_matrix(T_mat2,C_mat2);
	double R_axis[3]={0,0,1};
	SK_rotate_Cmatrix(T_mat1,R_axis,-(*rotate_list)[selected].theta1);
	SK_rotate_Cmatrix(T_mat2,R_axis,-(*rotate_list)[selected].theta2);
	
	int c1=(int)(thickness1/T_mat1[2][2])+1;
	int c2=(int)(thickness2/T_mat2[2][2])+1;
	(*rotate_list)[selected].trans1[2][2]=c1;
	(*rotate_list)[selected].trans2[2][2]=c2;
	
	for(i=0;i<3;i++)
	{
		T_mat1[2][i]*=c1;
		T_mat2[2][i]*=c2;
	}
	
	int limit1=(int)(0.5*SK_limitbox/sqrt(C_mat1[0][0]*C_mat1[0][0]+C_mat1[0][1]*C_mat1[0][1]))+1;
	int limit2=(int)(0.5*SK_limitbox/sqrt(C_mat1[1][0]*C_mat1[1][0]+C_mat1[1][1]*C_mat1[1][1]))+1;
	double mismatch_a;
	double mismatch_b;
	double mismatch_m=SK_mismatch;
	double inv_c1=1/T_mat1[2][2];
	double inv_c2=1/T_mat2[2][2];
	int    mat_out[2][2]={0,0,0,0};
	int	   i0,j0,k0,l0;
	if((T_mat1[1][0]*T_mat1[0][1]-T_mat1[0][0]*T_mat1[1][1])!=0)
	{
		j0=(int)((-T_mat1[2][0]*T_mat1[0][1]+T_mat1[2][1]*T_mat1[0][0])/(T_mat1[1][0]*T_mat1[0][1]-T_mat1[0][0]*T_mat1[1][1]));
		i0=(int)(( T_mat1[2][0]*T_mat1[1][1]-T_mat1[2][1]*T_mat1[1][0])/(T_mat1[1][0]*T_mat1[0][1]-T_mat1[0][0]*T_mat1[1][1]));
	}
	else
	{
		i0=j0=0;
	}
	
	if((T_mat2[1][0]*T_mat2[0][1]-T_mat2[0][0]*T_mat2[1][1])!=0)
	{
		l0=(int)((-T_mat2[2][0]*T_mat2[0][1]+T_mat2[2][1]*T_mat2[0][0])/(T_mat2[1][0]*T_mat2[0][1]-T_mat2[0][0]*T_mat2[1][1]));
		k0=(int)(( T_mat2[2][0]*T_mat2[1][1]-T_mat2[2][1]*T_mat2[1][0])/(T_mat2[1][0]*T_mat2[0][1]-T_mat2[0][0]*T_mat2[1][1]));
	}
	else
	{
		k0=l0=0;
	}
	for(i=i0-limit1+1;i<i0+limit1;i++) 
	for(j=j0-limit1+1;j<j0+limit1;j++)
	for(k=k0-limit2+1;k<k0+limit2;k++) 
	for(l=l0-limit2+1;l<l0+limit2;l++)
	{
		mismatch_a=fabs((i*T_mat1[0][0]+j*T_mat1[1][0]+T_mat1[2][0])*inv_c1
					   -(k*T_mat2[0][0]+l*T_mat2[1][0]+T_mat2[2][0])*inv_c2);
		mismatch_b=fabs((i*T_mat1[0][1]+j*T_mat1[1][1]+T_mat1[2][1])*inv_c1
					   -(k*T_mat2[0][1]+l*T_mat2[1][1]+T_mat2[2][1])*inv_c2);
		if(mismatch_a+mismatch_b<mismatch_m)
		{
			mat_out[0][0]=i;
			mat_out[0][1]=j;
			mat_out[1][0]=k;
			mat_out[1][1]=l;
			mismatch_m=mismatch_a+mismatch_b;
			(*rotate_list)[selected].mismatch[3]=mismatch_a;
			(*rotate_list)[selected].mismatch[4]=mismatch_b;				
			(*rotate_list)[selected].cylinder_flag=1;
		}
	}
	(*rotate_list)[selected].trans1[0][2]=mat_out[0][0];
	(*rotate_list)[selected].trans1[1][2]=mat_out[0][1];
	(*rotate_list)[selected].trans2[0][2]=mat_out[1][0];
	(*rotate_list)[selected].trans2[1][2]=mat_out[1][1];
}
*/