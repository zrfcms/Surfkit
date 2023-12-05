#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"
//original methods get minimize cell of current scene
int SK_SPG_overloop_vec(SPG_tools* spg,int type,double v_x,double v_y,double v_z,double symprec);
int SK_SPG_get_mintype(SPG_tools* spg);
int SK_SPG_overloop_atom(SPG_tools* spg,int n,double symprec);
void SK_SPG_affine_delete(SPG_tools* spg);
void SK_SPG_wrap_box(SPG_tools* spg,double symprec);
void SK_SPG_minimize_face(SPG_tools* spg,int v1,int v2,double symprec)//v1 v2 = (012) means abc
{
	//if(v1==v2) 
	//	printf("need two basis to determine a plane!!"); return
	int min_type=SK_SPG_get_mintype(spg);
	int start=0;
	int i,j,k;
	double mat1[3]={0,0,0};
	mat1[v1]=1;
	double mat2[3]={0,0,0};
	mat2[v2]=1;
	double l1=1;
	double l2=1;
	double trans[3][3];
	double tempp[3];//store position
	double temp_mat[3][3]={{1,0,0},{0,1,0},{0,0,1}};//store matrix
	double temp_mat2[3][3];
	double dot_product;
	double lenth;
	for(i=0;i<spg->num;i++)
		if(spg->type[i]==min_type)
		{
			start=i;
			break;
		}
	double vec[3];
	double vec2[3];
	int failed_flag;
	for(i=0;i<spg->num;i++)
	{
		if(spg->type[i]==min_type)
		if(i!=start)
		{
			//find two vector 
			vec[0]=spg->pos[i][0]-spg->pos[start][0];
			vec[1]=spg->pos[i][1]-spg->pos[start][1];
			vec[2]=spg->pos[i][2]-spg->pos[start][2];
			//get minimum vector
			if(vec[0]>0.5)vec[0]-=1;
			if(vec[1]>0.5)vec[1]-=1;
			if(vec[2]>0.5)vec[2]-=1;
			if(vec[0]<-0.5)vec[0]+=1;
			if(vec[1]<-0.5)vec[1]+=1;
			if(vec[2]<-0.5)vec[2]+=1;
			//make sure vector in plane
			if(fabs(vec[0])+fabs(vec[1])+fabs(vec[2])-fabs(vec[v1])-fabs(vec[v2])<symprec)//only search inside plane
			{	
				failed_flag=0;
				for(k=0;k<spg->num;k++)
				{
					vec2[0]=spg->pos[k][0]+vec[0];
					vec2[1]=spg->pos[k][1]+vec[1];
					vec2[2]=spg->pos[k][2]+vec[2];
					//check translational symmetry
					if(!SK_SPG_overloop_vec(spg,spg->type[k],vec2[0],vec2[1],vec2[2],symprec))
					{
						failed_flag=1;
						break;
					}
				}
				
				if(!failed_flag)//save the least two vector
				{
					lenth=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
					if(l1>lenth)
					{
						dot_product=mat2[0]*vec[0]+mat2[1]*vec[1]+mat2[2]*vec[2];
						//parallel is forbidden
						if(fabs(dot_product*dot_product/l2/lenth-1)>symprec)
						{
							for(j=0;j<3;j++)
								mat1[j]=vec[j];
							l1=lenth;
						}
					}
					if(l2>lenth)
					{
						dot_product=mat1[0]*vec[0]+mat1[1]*vec[1]+mat1[2]*vec[2];
						if(fabs(dot_product*dot_product/l1/lenth-1)>symprec)
						{
							for(j=0;j<3;j++)
								mat2[j]=vec[j];
							l2=lenth;
						}
					}
				}
			}
		}
	}
	for(i=0;i<3;i++)
	{
		temp_mat[i][v1]=mat1[i];
		temp_mat[i][v2]=mat2[i];
	}
	double det = temp_mat[0][0] * (temp_mat[1][1] * temp_mat[2][2] - temp_mat[1][2] * temp_mat[2][1])
			   + temp_mat[0][1] * (temp_mat[1][2] * temp_mat[2][0] - temp_mat[1][0] * temp_mat[2][2])
			   + temp_mat[0][2] * (temp_mat[1][0] * temp_mat[2][1] - temp_mat[1][1] * temp_mat[2][0]);
	trans[0][0] = (temp_mat[1][1] * temp_mat[2][2] - temp_mat[1][2] * temp_mat[2][1]) / det;
	trans[1][0] = (temp_mat[1][2] * temp_mat[2][0] - temp_mat[1][0] * temp_mat[2][2]) / det;
	trans[2][0] = (temp_mat[1][0] * temp_mat[2][1] - temp_mat[1][1] * temp_mat[2][0]) / det;
	trans[0][1] = (temp_mat[2][1] * temp_mat[0][2] - temp_mat[2][2] * temp_mat[0][1]) / det;
	trans[1][1] = (temp_mat[2][2] * temp_mat[0][0] - temp_mat[2][0] * temp_mat[0][2]) / det;
	trans[2][1] = (temp_mat[2][0] * temp_mat[0][1] - temp_mat[2][1] * temp_mat[0][0]) / det;
	trans[0][2] = (temp_mat[0][1] * temp_mat[1][2] - temp_mat[0][2] * temp_mat[1][1]) / det;
	trans[1][2] = (temp_mat[0][2] * temp_mat[1][0] - temp_mat[0][0] * temp_mat[1][2]) / det;
	trans[2][2] = (temp_mat[0][0] * temp_mat[1][1] - temp_mat[0][1] * temp_mat[1][0]) / det;
	//apply transformation
	for(i=0;i<spg->num;i++)
	{
		for(j=0;j<3;j++)
			tempp[j]=spg->pos[i][j];
		for(j=0;j<3;j++)
			spg->pos[i][j]=tempp[0]*trans[j][0]+tempp[1]*trans[j][1]+tempp[2]*trans[j][2];
	}
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			temp_mat2[i][j]=spg->mat[i][j];
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			spg->mat[i][j]=temp_mat2[i][0]*temp_mat[0][j]+temp_mat2[i][1]*temp_mat[1][j]+temp_mat2[i][2]*temp_mat[2][j];
		}
	//remove duplicate atom
	for(i=spg->num-1;i>0;i--)
	{
		if(SK_SPG_overloop_atom(spg,i,symprec))spg->type[i]=-1;
	}
	SK_SPG_affine_delete(spg);
	SK_SPG_wrap_box(spg,symprec);
	
}

void SK_SPG_minimize_cylinder(SPG_tools* spg,int v1,int v2,double symprec)//v1 v2 = (012) means abc
{
	//if(v1==v2) 
	//	printf("need two basis to determine a plane!!"); return
	int min_type=SK_SPG_get_mintype(spg);
	int start=0;
	int i,j,k;
	//only transform mat1
	int	   v3=2;
	for(i=0;i<3;i++)
		if(i!=v1&&i!=v2)v3=i;
	double mat1[3]={0,0,0};
	mat1[v3]=1;
	double trans[3][3];
	double tempp[3];//store position
	double temp_mat[3][3]={{1,0,0},{0,1,0},{0,0,1}};//store matrix
	double temp_mat2[3][3];
	
	for(i=0;i<spg->num;i++)
		if(spg->type[i]==min_type)
		{
			start=i;
			break;
		}
	
	double vec[3];
	double vec2[3];
	int failed_flag;
	for(i=0;i<spg->num;i++)
	{
		if(spg->type[i]==min_type)
		if(i!=start)
		{
			//find two vector 
			vec[0]=spg->pos[i][0]-spg->pos[start][0];
			vec[1]=spg->pos[i][1]-spg->pos[start][1];
			vec[2]=spg->pos[i][2]-spg->pos[start][2];
			//get minimum vector
			if(vec[0]>0.5)vec[0]-=1;
			if(vec[1]>0.5)vec[1]-=1;
			if(vec[2]>0.5)vec[2]-=1;
			if(vec[0]<-0.5)vec[0]+=1;
			if(vec[1]<-0.5)vec[1]+=1;
			if(vec[2]<-0.5)vec[2]+=1;
			//make sure vector positive
			if(vec[v3]<0)
				for(k=0;k<3;k++)vec[k]*=-1;
			//make sure vector in plane
			if(fabs(vec[0])+fabs(vec[1])+fabs(vec[2])-fabs(vec[v1])-fabs(vec[v2])>symprec)//only search out of plane
			{	
				failed_flag=0;
				for(k=0;k<spg->num;k++)
				{
					vec2[0]=spg->pos[k][0]+vec[0];
					vec2[1]=spg->pos[k][1]+vec[1];
					vec2[2]=spg->pos[k][2]+vec[2];
					//check translational symmetry
					if(!SK_SPG_overloop_vec(spg,spg->type[k],vec2[0],vec2[1],vec2[2],symprec))
					{
						failed_flag=1;
						break;
					}
				}
				if(!failed_flag)
				{
					if((mat1[0]*mat1[0]+mat1[1]*mat1[1]+mat1[2]*mat1[2])>fabs(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]))
					{
						for(j=0;j<3;j++)
							mat1[j]=vec[j];
					}
				}
			}
		}
	}
	
	for(i=0;i<3;i++)
		temp_mat[i][v3]=mat1[i];
	double det = temp_mat[0][0] * (temp_mat[1][1] * temp_mat[2][2] - temp_mat[1][2] * temp_mat[2][1])
			   + temp_mat[0][1] * (temp_mat[1][2] * temp_mat[2][0] - temp_mat[1][0] * temp_mat[2][2])
			   + temp_mat[0][2] * (temp_mat[1][0] * temp_mat[2][1] - temp_mat[1][1] * temp_mat[2][0]);
	trans[0][0] = (temp_mat[1][1] * temp_mat[2][2] - temp_mat[1][2] * temp_mat[2][1]) / det;
	trans[1][0] = (temp_mat[1][2] * temp_mat[2][0] - temp_mat[1][0] * temp_mat[2][2]) / det;
	trans[2][0] = (temp_mat[1][0] * temp_mat[2][1] - temp_mat[1][1] * temp_mat[2][0]) / det;
	trans[0][1] = (temp_mat[2][1] * temp_mat[0][2] - temp_mat[2][2] * temp_mat[0][1]) / det;
	trans[1][1] = (temp_mat[2][2] * temp_mat[0][0] - temp_mat[2][0] * temp_mat[0][2]) / det;
	trans[2][1] = (temp_mat[2][0] * temp_mat[0][1] - temp_mat[2][1] * temp_mat[0][0]) / det;
	trans[0][2] = (temp_mat[0][1] * temp_mat[1][2] - temp_mat[0][2] * temp_mat[1][1]) / det;
	trans[1][2] = (temp_mat[0][2] * temp_mat[1][0] - temp_mat[0][0] * temp_mat[1][2]) / det;
	trans[2][2] = (temp_mat[0][0] * temp_mat[1][1] - temp_mat[0][1] * temp_mat[1][0]) / det;
	//apply transformation
	for(i=0;i<spg->num;i++)
	{
		for(j=0;j<3;j++)
			tempp[j]=spg->pos[i][j];
		for(j=0;j<3;j++)
			spg->pos[i][j]=tempp[0]*trans[j][0]+tempp[1]*trans[j][1]+tempp[2]*trans[j][2];
	}
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			temp_mat2[i][j]=spg->mat[i][j];
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			spg->mat[i][j]=temp_mat2[i][0]*temp_mat[0][j]+temp_mat2[i][1]*temp_mat[1][j]+temp_mat2[i][2]*temp_mat[2][j];
		}
	//remove duplicate atom
	for(i=spg->num-1;i>0;i--)
	{
		if(SK_SPG_overloop_atom(spg,i,symprec))spg->type[i]=-1;
	}
	SK_SPG_affine_delete(spg);
	SK_SPG_wrap_box(spg,symprec);
	
}

int SK_SPG_pick_3Dlattice(SPG_tools* spg,int v1,int v2,double symprec)//v1 v2 = (012) means abc
{
	//if(v1==v2) 
	//	printf("need two basis to determine a plane!!"); return
	int i,j,k,ii;
	//only transform mat1
	int	   v3;
	for(i=0;i<3;i++)
		if(i!=v1&&i!=v2)v3=i;
	double mat1[3]={0,0,0};
	mat1[v3]=1;
	double trans[3][3];
	double tempp[3];//store position
	double temp_mat[3][3]={{1,0,0},{0,1,0},{0,0,1}};//store matrix
	double temp_mat2[3][3];
	double max_v3=0;
	double min_v3=1;
	int    start=-1;
	for(i=0;i<spg->num;i++)
	{
		if(spg->pos[i][v3]>max_v3)
			max_v3=spg->pos[i][v3];
		if(spg->pos[i][v3]<min_v3)
			min_v3=spg->pos[i][v3];
	}
	double vec[3];
	double vec2[3];
	int failed_flag;
	for(i=0;i<spg->num;i++)
	for(ii=i+1;ii<spg->num;ii++)
	{
		if(spg->type[i]==spg->type[ii])
		{
			//find two vector 
			vec[0]=spg->pos[i][0]-spg->pos[ii][0];
			vec[1]=spg->pos[i][1]-spg->pos[ii][1];
			vec[2]=spg->pos[i][2]-spg->pos[ii][2];
			//make sure vector in plane
			if(fabs(vec[0])+fabs(vec[1])+fabs(vec[2])-fabs(vec[v1])-fabs(vec[v2])>symprec)//only search out of plane
			{	
				failed_flag=0;
				for(k=0;k<spg->num;k++)
				{
					vec2[0]=spg->pos[k][0]+vec[0];
					vec2[1]=spg->pos[k][1]+vec[1];
					vec2[2]=spg->pos[k][2]+vec[2];
					//check translational symmetry
					//skip outside surface
					if(vec2[v3]<=max_v3&&vec2[v3]>=min_v3)
					if(!SK_SPG_overloop_vec(spg,spg->type[k],vec2[0],vec2[1],vec2[2],symprec))
					{
						failed_flag=1;
						break;
					}
				}
				if(!failed_flag)
				{
					if(fabs(mat1[v3])>fabs(vec[v3]))
					{
						for(j=0;j<3;j++)
							mat1[j]=vec[j];
						//start=ii;
					}
				}
			}
		}
	}
	for(i=0;i<3;i++)
		temp_mat[i][v3]=mat1[i];
	double det = temp_mat[0][0] * (temp_mat[1][1] * temp_mat[2][2] - temp_mat[1][2] * temp_mat[2][1])
			   + temp_mat[0][1] * (temp_mat[1][2] * temp_mat[2][0] - temp_mat[1][0] * temp_mat[2][2])
			   + temp_mat[0][2] * (temp_mat[1][0] * temp_mat[2][1] - temp_mat[1][1] * temp_mat[2][0]);
	trans[0][0] = (temp_mat[1][1] * temp_mat[2][2] - temp_mat[1][2] * temp_mat[2][1]) / det;
	trans[1][0] = (temp_mat[1][2] * temp_mat[2][0] - temp_mat[1][0] * temp_mat[2][2]) / det;
	trans[2][0] = (temp_mat[1][0] * temp_mat[2][1] - temp_mat[1][1] * temp_mat[2][0]) / det;
	trans[0][1] = (temp_mat[2][1] * temp_mat[0][2] - temp_mat[2][2] * temp_mat[0][1]) / det;
	trans[1][1] = (temp_mat[2][2] * temp_mat[0][0] - temp_mat[2][0] * temp_mat[0][2]) / det;
	trans[2][1] = (temp_mat[2][0] * temp_mat[0][1] - temp_mat[2][1] * temp_mat[0][0]) / det;
	trans[0][2] = (temp_mat[0][1] * temp_mat[1][2] - temp_mat[0][2] * temp_mat[1][1]) / det;
	trans[1][2] = (temp_mat[0][2] * temp_mat[1][0] - temp_mat[0][0] * temp_mat[1][2]) / det;
	trans[2][2] = (temp_mat[0][0] * temp_mat[1][1] - temp_mat[0][1] * temp_mat[1][0]) / det;
	//apply transformation
	if(start>0)
	{
		vec[0]=spg->pos[start][0];
		vec[1]=spg->pos[start][1];
		vec[2]=spg->pos[start][2];
	}
	else
	{
		vec[0]=vec[1]=vec[2]=0;
	}
	for(i=0;i<spg->num;i++)
	{
		for(j=0;j<3;j++)
			tempp[j]=spg->pos[i][j]-vec[j];
		for(j=0;j<3;j++)
			spg->pos[i][j]=tempp[0]*trans[j][0]+tempp[1]*trans[j][1]+tempp[2]*trans[j][2];
	}
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			temp_mat2[i][j]=spg->mat[i][j];
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			spg->mat[i][j]=temp_mat2[i][0]*temp_mat[0][j]+temp_mat2[i][1]*temp_mat[1][j]+temp_mat2[i][2]*temp_mat[2][j];
		}
	//remove duplicate atom
	for(i=spg->num-1;i>0;i--)
	{
		if(SK_SPG_overloop_atom(spg,i,symprec))spg->type[i]=-1;
	}
	for(i=0;i<spg->num;i++)
		if(spg->type[i]>0)
		{
			start=i;
			break;
		}
	SK_SPG_affine_delete(spg);
	SK_SPG_wrap_box(spg,symprec);
	if(start>=0)
		return start;
	return -1;//not found
}

int SK_SPG_overloop_vec(SPG_tools* spg,int type,double v_x,double v_y,double v_z,double symprec)
{
	int i,j;
	double symprec2=symprec*symprec;
	double vec[3];	
	for(i=0;i<spg->num;i++)
	{
		//if(spg->type[i]!=type)
		//	continue;
			
		vec[0]=spg->pos[i][0]-v_x;
		vec[1]=spg->pos[i][1]-v_y;
		vec[2]=spg->pos[i][2]-v_z;
		for(j=0;j<3;j++)
		{
			while(vec[j]<-0.5)vec[j]+=1;
			while(vec[j]>0.5)vec[j]-=1;
		}
		if((vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])<symprec2)
		{
		   return 1;
		}
	}
	return 0;
}

int SK_SPG_get_mintype(SPG_tools* spg)
{
	int i;
	int ct=spg->num;
	int rt=1;
	int type_num=0;
	//for safety check TypeNum first
	for(i=0;i<spg->num;i++)
	{
		if(type_num<spg->type[i])
			type_num=spg->type[i];
	}
	if(type_num<=1)return 1;
	int *type_list=(int*)malloc(sizeof(int)*type_num);
	for(i=0;i<type_num;i++)
		type_list[i]=0;
	
	for(i=0;i<spg->num;i++)
	{
		if(spg->type[i]>0)
			type_list[spg->type[i]-1]++;
	}
	for(i=0;i<type_num;i++)
		if(type_list[i]!=0)
			if(type_list[i]<ct)
			{
				ct=type_list[i];
				rt=i+1;
			}
	free(type_list);
    return rt;
}

void SK_SPG_wrap_box(SPG_tools* spg,double symprec)
{
	int i,j;
	for(i=0;i<spg->num;i++)
	for(j=0;j<3;j++)
	{
		if (spg->pos[i][j] < -symprec)
			spg->pos[i][j] = spg->pos[i][j] + 1.0 - (int)(spg->pos[i][j]+symprec);
		else
			spg->pos[i][j] = spg->pos[i][j] - (int) (spg->pos[i][j]+symprec);
	}
	
}

void SK_SPG_affine_delete(SPG_tools* spg)
{
	int i,j;
	int initial_id=spg->num;
	double tempd;
	int tempi;
	for(i=0;i<spg->num;i++)
	{
		if(spg->type[i]<0)
			for(j=spg->num-1;j>=0;j--)
			{
				if(i==j)
				{	
					spg->num--;
					break;
				}
				if(spg->type[j]<0)
				{
					spg->num--;
					continue;
				}
				else
				{
					tempd=spg->pos[i][0];
					spg->pos[i][0]=spg->pos[j][0];
					spg->pos[j][0]=tempd;
					tempd=spg->pos[i][1];
					spg->pos[i][1]=spg->pos[j][1];
					spg->pos[j][1]=tempd;
					tempd=spg->pos[i][2];
					spg->pos[i][2]=spg->pos[j][2];
					spg->pos[j][2]=tempd;
					tempi=spg->type[i];
					spg->type[i]=spg->type[j];
					spg->type[j]=tempi;
					spg->num--;
					break;
				}
			}
	}
	if(initial_id!=spg->num)
	{
		spg->pos=(double(*)[3])realloc(spg->pos,3*spg->num*sizeof(double));
		spg->type=(int*)realloc(spg->type,spg->num*sizeof(int));
	}
}

int SK_SPG_overloop_atom(SPG_tools* spg,int n,double symprec)
{
	int i,j;
	double symprec2=symprec*symprec;
	double vec[3];	
	for(i=0;i<n;i++)
	{
		if(spg->type[i]!=spg->type[n])
			continue;
		vec[0]=spg->pos[i][0]-spg->pos[n][0];
		vec[1]=spg->pos[i][1]-spg->pos[n][1];
		vec[2]=spg->pos[i][2]-spg->pos[n][2];
		for(j=0;j<3;j++)
		{
			while(vec[j]<-0.5)vec[j]+=1;
			while(vec[j]>0.5)vec[j]-=1;
		}
		if((vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])<symprec2)
		{
		   return 1;
		}
	}
	return 0;
}