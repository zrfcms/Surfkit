#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"
static double point_product(double vec1[3],double vec2[3])
{
	return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
}

static void switch_vec(double vec1[3],double vec2[3])
{
	double temp;
	for(int i=0;i<3;i++)
	{
		temp=vec1[i];
		vec1[i]=vec2[i];
		vec2[i]=temp;
	}
}

void Gaussian_lattice_reduction(double vec1[3],double vec2[3])
{
	double vec3[3];
	if(point_product(vec1,vec2)<0)
		for(int i=0;i<3;i++)
			vec2[i]*=-1;
	
	if(point_product(vec1,vec1)<point_product(vec2,vec2))
		switch_vec(vec1,vec2);
		
	do
	{
		double rate=floor(point_product(vec1,vec2)/point_product(vec2,vec2));
		for(int i=0;i<3;i++)
		{
			vec3[i]=vec1[i]-floor(rate)*vec2[i];
			vec1[i]=vec2[i];
			vec2[i]=vec3[i];
		}
	}
	while(point_product(vec1,vec1)<point_product(vec2,vec2));
	
	if(point_product(vec1,vec2)<0.5*point_product(vec2,vec2))
		switch_vec(vec1,vec2);
	else 
	{
		for(int i=0;i<3;i++)
			vec3[i]=vec1[i]-vec2[i];
		if(point_product(vec3,vec3)>point_product(vec2,vec2))
			for(int i=0;i<3;i++)
			{
				vec1[i]=vec2[i];
				vec2[i]=vec3[i];
			}
		else 
			for(int i=0;i<3;i++)
			{
				vec1[i]=vec3[i];
				vec2[i]=-vec2[i];
			}
	}
}

void Right_hand_reduction(double mat[3][3])
{
	double A=(mat[0][0]*mat[1][1]*mat[2][2]+mat[0][1]*mat[1][2]*mat[2][0]+mat[0][2]*mat[1][0]*mat[2][1]
	  -mat[0][0]*mat[1][2]*mat[2][1]-mat[0][1]*mat[1][0]*mat[2][2]-mat[0][2]*mat[1][1]*mat[2][0]);
	if(A<0)
	for(int i=0;i<3;i++)
		mat[1][i]*=-1;
}