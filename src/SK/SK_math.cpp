#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"
void SK_cross_product(double out[3],double in1[3],double in2[3])
{
	out[0]=in1[1]*in2[2]-in1[2]*in2[1];
	out[1]=in1[2]*in2[0]-in1[0]*in2[2];
	out[2]=in1[0]*in2[1]-in1[1]*in2[0];
}

double SK_min_height(double mat[3][3])
{
	double xproduct[3];
	double xin1[3]={mat[0][0],mat[1][0],mat[2][0]};	
	double xin2[3]={mat[0][1],mat[1][1],mat[2][1]};
	double xin3[3]={mat[0][2],mat[1][2],mat[2][2]};
	
	SK_cross_product(xproduct,xin1,xin2);
	double normal_c=sqrt(xproduct[0]*xproduct[0]
						+xproduct[1]*xproduct[1]
						+xproduct[2]*xproduct[2]);
	double volumn=xproduct[0]*mat[2][0]
				 +xproduct[1]*mat[2][1]
				 +xproduct[2]*mat[2][2];
			 
	double min_height=volumn/normal_c;
	
	SK_cross_product(xproduct,xin1,xin3);
	double normal_b=sqrt(xproduct[0]*xproduct[0]
						+xproduct[1]*xproduct[1]
						+xproduct[2]*xproduct[2]);
	if(volumn/normal_b<min_height)
		min_height=volumn/normal_b;
		
	SK_cross_product(xproduct,xin2,xin3);	
	double normal_a=sqrt(xproduct[0]*xproduct[0]
						+xproduct[1]*xproduct[1]
						+xproduct[2]*xproduct[2]);
	if(volumn/normal_a<min_height)
		min_height=volumn/normal_a;	
		
	return min_height;
}

void SK_rotate_matrix(double F_mat[3][3],double vv[3], double angle)
{
	int i;
    angle*=M_PI/180.0;
    double normalize;
    normalize=sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
    double u = vv[0]/normalize;
    double v = vv[1]/normalize;
    double w = vv[2]/normalize;
	for(i=0;i<3;i++)
	{
		double x = F_mat[0][i];
		double y = F_mat[1][i];
		double z = F_mat[2][i];
		F_mat[0][i]=x*(u*u+(v*v+w*w)*cos(angle))+y*(u*v*(1-cos(angle))+w*sin(angle))+z*(u*w*(1-cos(angle))-v*sin(angle));
		F_mat[1][i]=x*(u*v*(1-cos(angle))-w*sin(angle))+y*(v*v+(u*u+w*w)*cos(angle))+z*(v*w*(1-cos(angle))+u*sin(angle));
		F_mat[2][i]=x*(u*w*(1-cos(angle))+v*sin(angle))+y*(v*w*(1-cos(angle))-u*sin(angle))+z*(w*w+(u*u+v*v)*cos(angle));
	}
}

void SK_rotate_Cmatrix(double C_mat[3][3],double vv[3], double angle)
{
	int i;
    angle*=M_PI/180.0;
    double normalize;
    normalize=sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
    double u = vv[0]/normalize;
    double v = vv[1]/normalize;
    double w = vv[2]/normalize;
	for(i=0;i<3;i++)
	{
		double x = C_mat[i][0];
		double y = C_mat[i][1];
		double z = C_mat[i][2];
		C_mat[i][0]=x*(u*u+(v*v+w*w)*cos(angle))+y*(u*v*(1-cos(angle))+w*sin(angle))+z*(u*w*(1-cos(angle))-v*sin(angle));
		C_mat[i][1]=x*(u*v*(1-cos(angle))-w*sin(angle))+y*(v*v+(u*u+w*w)*cos(angle))+z*(v*w*(1-cos(angle))+u*sin(angle));
		C_mat[i][2]=x*(u*w*(1-cos(angle))+v*sin(angle))+y*(v*w*(1-cos(angle))-u*sin(angle))+z*(w*w+(u*u+v*v)*cos(angle));
	}
}

void SK_copy_matrix(double output[3][3],double input[3][3])
{
	int i,j;
	for(i=0;i<3;i++)
	for(j=0;j<3;j++)
		output[i][j]=input[i][j];
}

void SK_inverse_matrix(double input[3][3],double output[3][3])
{
	double temp[3][3];
	double A=input[0][0]*input[1][1]*input[2][2]+input[0][1]*input[1][2]*input[2][0]+input[0][2]*input[1][0]*input[2][1]
		    -input[0][0]*input[1][2]*input[2][1]-input[0][1]*input[1][0]*input[2][2]-input[0][2]*input[1][1]*input[2][0];
	if(A==0)
	{
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				temp[i][j]=input[i][j];
	}
	else 
	{
		A=1.0/A;
		temp[0][0]=A*( input[1][1]*input[2][2]-input[1][2]*input[2][1]);
		temp[1][0]=A*(-input[1][0]*input[2][2]+input[1][2]*input[2][0]);
		temp[2][0]=A*( input[1][0]*input[2][1]-input[1][1]*input[2][0]);
		temp[0][1]=A*(-input[0][1]*input[2][2]+input[0][2]*input[2][1]);
		temp[1][1]=A*( input[0][0]*input[2][2]-input[0][2]*input[2][0]);
		temp[2][1]=A*(-input[0][0]*input[2][1]+input[0][1]*input[2][0]);
		temp[0][2]=A*( input[0][1]*input[1][2]-input[0][2]*input[1][1]);
		temp[1][2]=A*(-input[0][0]*input[1][2]+input[0][2]*input[1][0]);
		temp[2][2]=A*( input[0][0]*input[1][1]-input[0][1]*input[1][0]);
	}	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			output[i][j]=temp[i][j];
}

void SK_transposition_matrix(double input[3][3],double output[3][3])
{
	double temp[3][3];
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			temp[j][i]=input[i][j];
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			output[i][j]=temp[i][j];		
}

void SK_matrix_multiplication(double input1[3][3],double input2[3][3],double output[3][3])
{
	double temp[3][3];
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		temp[i][j]=0;
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	for(int k=0;k<3;k++)
		temp[i][j]+=input1[i][k]*input2[k][j];
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		output[i][j]=temp[i][j];
}

void SK_matrix_left_multiplication(double input1[3][3],double input2[3],double output[3])
{
	double temp[3];
	for(int i=0;i<3;i++)
		temp[i]=0;
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		temp[i]+=input1[i][j]*input2[j];
	for(int i=0;i<3;i++)
		output[i]=temp[i];
}


int SK_getrp(int a, int b)
{
	if(a==0)
	{
		if(b>0)
			return b;
		return 1;
	}
	if(b==0)
	{
		if(a>0)
			return a;
		return 1;
	}
	if(a==1||b==1)
		return 1;
	while(1)
    {
		int t = a%b;
		if(t == 0) 
        {
            break;
        }
		else
        {
			a = b;
			b = t;
		}
	}
	if(b>1)	return b;
	else return 1;
}

void SK_matrix_print(double input[3][3])
{
	printf("Current Matrix is:\n");
	printf("%lf\t%lf\t%lf\n",input[0][0],input[0][1],input[0][2]);
	printf("%lf\t%lf\t%lf\n",input[1][0],input[1][1],input[1][2]);
	printf("%lf\t%lf\t%lf\n",input[2][0],input[2][1],input[2][2]);
}