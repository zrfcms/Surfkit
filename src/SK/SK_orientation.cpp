#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"
extern int SK_SPG_overloop_atom(SPG_tools* spg,int n,double symprec);
void SK_SPG_wrap_box(SPG_tools* spg,double symprec);
//find a orientation
//it doesn't matter what it is
//we will standardize it with minimization
void SK_auto_Frot(int c[3],int F_rot[3][3])
{
	int x1[3]={0,0,0};
	int x2[3]={0,0,0};
	if(c[2]!=0)//safe to use any vector use to skip calculation
	{
		x1[0]=1;
		x2[1]=1;
	}
	else if(c[0]!=0)//a axis use bc
	{
		x1[1]=1;
		x2[2]=1;
	}
	else
	{
		x1[0]=1;
		x2[2]=1;
	}
	//F_rot[0]=c X x1
	F_rot[0][0]=c[1]*x1[2]-c[2]*x1[1];
	F_rot[1][0]=c[2]*x1[0]-c[0]*x1[2];
	F_rot[2][0]=c[0]*x1[1]-c[1]*x1[0];
	//F_rot[1]=c X x2
	F_rot[0][1]=c[1]*x2[2]-c[2]*x2[1];
	F_rot[1][1]=c[2]*x2[0]-c[0]*x2[2];
	F_rot[2][1]=c[0]*x2[1]-c[1]*x2[0];
	//F_rot[2]=c_str
	F_rot[0][2]=c[0];
	F_rot[1][2]=c[1];
	F_rot[2][2]=c[2];
}
extern void SK_SPG_affine_delete(SPG_tools* spg);
void SK_redefine(SPG_tools* spg,int F_rot[3][3],double symprec)
{
	int i,j,k,m,n;
	double mat[3][3];
	double(*temppos)[3];
	double inverse[3][3];
	if((F_rot[0][0]*F_rot[1][1]*F_rot[2][2]+F_rot[0][1]*F_rot[1][2]*F_rot[2][0]+F_rot[0][2]*F_rot[1][0]*F_rot[2][1]-
	    F_rot[0][0]*F_rot[1][2]*F_rot[2][1]-F_rot[0][2]*F_rot[1][1]*F_rot[2][0]-F_rot[0][1]*F_rot[1][0]*F_rot[2][2])==0)
		return;
	if(spg->num==0)return;
	
	//duplicate box
	int max_sum=(fabs(F_rot[0][0])+fabs(F_rot[0][1])+fabs(F_rot[0][2]))
			   *(fabs(F_rot[1][0])+fabs(F_rot[1][1])+fabs(F_rot[1][2]))
			   *(fabs(F_rot[2][0])+fabs(F_rot[2][1])+fabs(F_rot[2][2]))*spg->num;
	
	spg->pos=(double(*)[3])realloc(spg->pos,3*max_sum*sizeof(double));
	temppos=(double(*)[3])malloc(3*max_sum*sizeof(double));
	spg->type=(int*)realloc(spg->type,max_sum*sizeof(int));
	n=0;
	for(i=0;i<fabs(F_rot[0][0])+fabs(F_rot[0][1])+fabs(F_rot[0][2]);i++)
	for(j=0;j<fabs(F_rot[1][0])+fabs(F_rot[1][1])+fabs(F_rot[1][2]);j++)
	for(k=0;k<fabs(F_rot[2][0])+fabs(F_rot[2][1])+fabs(F_rot[2][2]);k++)
	for(k=0;k<fabs(F_rot[2][0])+fabs(F_rot[2][1])+fabs(F_rot[2][2]);k++)
	for(m=0;m<spg->num;m++)
	{
		temppos[n][0]=spg->pos[m][0]+(double)i;
		temppos[n][1]=spg->pos[m][1]+(double)j;
		temppos[n][2]=spg->pos[m][2]+(double)k;
		spg->type[n]=spg->type[m];
		n++;
	}
	spg->num=max_sum;
	
	//remap current position to new fractional coordinates 
	double det=(double)(F_rot[0][0] * (F_rot[1][1] * F_rot[2][2] - F_rot[1][2] * F_rot[2][1])
					  + F_rot[0][1] * (F_rot[1][2] * F_rot[2][0] - F_rot[1][0] * F_rot[2][2])
					  + F_rot[0][2] * (F_rot[1][0] * F_rot[2][1] - F_rot[1][1] * F_rot[2][0]));
	inverse[0][0] = (double)(F_rot[1][1] * F_rot[2][2] - F_rot[1][2] * F_rot[2][1]) / det;
	inverse[1][0] = (double)(F_rot[1][2] * F_rot[2][0] - F_rot[1][0] * F_rot[2][2]) / det;
	inverse[2][0] = (double)(F_rot[1][0] * F_rot[2][1] - F_rot[1][1] * F_rot[2][0]) / det;
	inverse[0][1] = (double)(F_rot[2][1] * F_rot[0][2] - F_rot[2][2] * F_rot[0][1]) / det;
	inverse[1][1] = (double)(F_rot[2][2] * F_rot[0][0] - F_rot[2][0] * F_rot[0][2]) / det;
	inverse[2][1] = (double)(F_rot[2][0] * F_rot[0][1] - F_rot[2][1] * F_rot[0][0]) / det;
	inverse[0][2] = (double)(F_rot[0][1] * F_rot[1][2] - F_rot[0][2] * F_rot[1][1]) / det;
	inverse[1][2] = (double)(F_rot[0][2] * F_rot[1][0] - F_rot[0][0] * F_rot[1][2]) / det;
	inverse[2][2] = (double)(F_rot[0][0] * F_rot[1][1] - F_rot[0][1] * F_rot[1][0]) / det;
  
	for(i=0;i<spg->num;i++)
		for(j=0;j<3;j++)
			spg->pos[i][j]=inverse[j][0]*temppos[i][0]+inverse[j][1]*temppos[i][1]+inverse[j][2]*temppos[i][2];	
	
	free(temppos);
	
	//reset box
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			mat[i][j]=spg->mat[i][j];	
			
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			spg->mat[i][j]=mat[i][0]*F_rot[0][j]+mat[i][1]*F_rot[1][j]+mat[i][2]*F_rot[2][j];		
			
	for(i=spg->num-1;i>0;i--)
	{
		if(SK_SPG_overloop_atom(spg,i,symprec))spg->type[i]=-1;
	}
	
	SK_SPG_wrap_box(spg,symprec);
	SK_SPG_affine_delete(spg);
}
