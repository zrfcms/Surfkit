#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"
#define MAXDANGLE 0.5
#define MAXROTATE 361
bool SK_isrp(int a, int b);
int SK_install_face(double theta,double dis,int a,int b,SK_crystal_face** input,int n);
double SK_arc_tan(double x,double y);

double SK_limitd=20.0;
double SK_limitdmin=0.0;
double SK_limitbox=100.0;
double SK_limitdt=20;
double SK_mismatch=0.005;
double SK_mingamma=20;
double SK_maxgamma=160;
double SK_maxdangle=1;
void SK_mismatch_set(double max_d,double min_d,double max_box,double limitdt,double min_gamma,
					 double max_gamma,double max_mismatch,double max_dangle)
{
    SK_limitd=max_d;
	SK_limitdmin=min_d;
    SK_limitbox=max_box;
	SK_limitdt=limitdt;
    SK_mismatch=max_mismatch;
	SK_mingamma=min_gamma;
	SK_maxgamma=max_gamma;
	SK_maxdangle=max_dangle;
}

int SK_mismatch_scan(double C_mat1[3][3],double C_mat2[3][3],SK_scaned_data** rotate_list)
{
    SK_crystal_face* face1;
    int fn1=0;
    SK_crystal_face* face2;
    int fn2=0;
    //FILE *output;
    int N_list=0;
    int i,j,k,m1=1,m2=1,n1=1,n2=1,o,p;
    double temp1,temp2,mismatch1,mismatch2;
	double T_mat1[3][3],T_mat2[3][3];
	SK_copy_matrix(T_mat1,C_mat1);
	SK_copy_matrix(T_mat2,C_mat2);
    //66666666666666666666666666666666666666
    //printf("step II: build face list\n");
    //66666666666666666666666666666666666666
	o=(int)(SK_limitd/sqrt(T_mat1[0][0]*T_mat1[0][0]+T_mat1[0][1]*T_mat1[0][1]))+1;
	p=(int)(SK_limitd/sqrt(T_mat1[1][0]*T_mat1[1][0]+T_mat1[1][1]*T_mat1[1][1]))+1;
    for(i=-o+1;i<o;i++) 
	for(j=-p+1;j<p;j++) 
    {
		if(SK_isrp(abs(i),abs(j)))
	    {
	        temp1=T_mat1[0][0]*i+T_mat1[1][0]*j;
	        temp2=T_mat1[0][1]*i+T_mat1[1][1]*j;
			if(sqrt(temp1*temp1+temp2*temp2)<SK_limitd)
			if(sqrt(temp1*temp1+temp2*temp2)>SK_limitdmin)
				fn1=SK_install_face(SK_arc_tan(temp1,temp2),sqrt(temp1*temp1+temp2*temp2),i,j,&face1,fn1);
	    }
    }
    o=(int)(SK_limitd/sqrt(C_mat2[0][0]*C_mat2[0][0]+C_mat2[0][1]*C_mat2[0][1]))+1;
	p=(int)(SK_limitd/sqrt(C_mat2[1][0]*C_mat2[1][0]+C_mat2[1][1]*C_mat2[1][1]))+1;
    for(i=-o+1;i<o;i++) 
	for(j=-p+1;j<p;j++) 
    {
		if(SK_isrp(abs(i),abs(j)))
	    {
	        temp1=C_mat2[0][0]*i+C_mat2[1][0]*j;
	        temp2=C_mat2[0][1]*i+C_mat2[1][1]*j;
			if(sqrt(temp1*temp1+temp2*temp2)<SK_limitd)
			if(sqrt(temp1*temp1+temp2*temp2)>SK_limitdmin)
				fn2=SK_install_face(SK_arc_tan(temp1,temp2),sqrt(temp1*temp1+temp2*temp2),i,j,&face2,fn2);
	    }
    }
    //66666666666666666666666666666666666666
    //printf("step IV conbine two groups' data\n");
    //66666666666666666666666666666666666666
    //output=fopen("mismatch.log","w");
    //fprintf(output,"theta1\ttheta2\t");
    //fprintf(output,"dy1\tdy2\t");
    //fprintf(output,"dz1\tdz2\t");
    //fprintf(output,"ny1\tny2\tmismatchy\t");
    //fprintf(output,"nz1\tnz2\tmismatchz\n");
    for(i=0;i<fn1;i++)
	{
		printf(">>Pairing %d/Total %d\n",i,fn1);
	for(o=0;o<fn1;o++)
	if(o!=i)
	{
		face1[i].dtheta=face1[o].theta-face1[i].theta;
		face1[i].id_guest=o;
		if(face1[i].dtheta>SK_mingamma&&face1[i].dtheta<SK_maxgamma)
		{
			for(j=0;j<fn2;j++)			
			for(p=0;p<fn2;p++)
			if(j!=p)
			{
				face2[j].dtheta=face2[p].theta-face2[j].theta;
				face2[j].id_guest=p;
				if(fabs(face2[j].dtheta-face1[i].dtheta)<SK_maxdangle)
				{
					mismatch1=mismatch2=SK_mismatch;
					for(k=1;((double)k*face1[i].dis)<SK_limitbox;k++)
					{
						temp1=(((double)k*face1[i].dis)-face2[j].dis*(double)(int)(((double)k*face1[i].dis)/face2[j].dis))/((double)k*face1[i].dis);
						temp2=(((double)k*face1[i].dis)-face2[j].dis*(double)(1+(int)(((double)k*face1[i].dis)/face2[j].dis)))/((double)k*face1[i].dis);
						if(temp1>fabs(temp2))
							temp1=temp2;
						if(fabs(temp1)<mismatch1)
						{
							mismatch1=fabs(temp1);
							m1=k;
							if(temp1<0)
							m2=1.0+(int)(((double)k*face1[i].dis)/face2[j].dis);
							else m2=(int)(((double)k*face1[i].dis)/face2[j].dis);
					   
						}
						if(mismatch1<SK_mismatch)
						break;
					}
					for(k=1;((double)k*face1[face1[i].id_guest].dis)<SK_limitbox;k++)
					{
						temp1=(
						((double)k*face1[face1[i].id_guest].dis)-
						face2[face2[j].id_guest].dis*
						(double)(int)(((double)k*face1[face1[i].id_guest].dis)/face2[face2[j].id_guest].dis)
						)/((double)k*face1[face1[i].id_guest].dis);

						temp2=(((double)k*face1[face1[i].id_guest].dis)-
						face2[face2[j].id_guest].dis*
						(double)(1+(int)(((double)k*face1[face1[i].id_guest].dis)/face2[face2[j].id_guest].dis))
						)/((double)k*face1[face1[i].id_guest].dis);
						if(temp1>fabs(temp2))
							temp1=temp2;
						if(fabs(temp1)<mismatch2)
						{
							mismatch2=fabs(temp1);
							n1=k;
							if(temp1<0)
								n2=1.0+(int)(((double)k*face1[face1[i].id_guest].dis)/face2[face2[j].id_guest].dis);
							else n2=(int)(((double)k*face1[face1[i].id_guest].dis)/face2[face2[j].id_guest].dis);     
						}
						if(mismatch2<SK_mismatch)
							break;
					}
					if(mismatch1<SK_mismatch)
					if(mismatch2<SK_mismatch)
					if(fabs(face1[i].dtheta-face2[j].dtheta)<SK_maxdangle)
					//if(face1[i].theta<SK_limitdt)
					if(fabs(face2[j].theta-face1[i].theta)<SK_limitdt)
					{
						int saved_flag=0;
						for(k=0;k<N_list;k++)
							if(fabs((*rotate_list)[k].theta2-(*rotate_list)[k].theta1-face2[j].theta+face1[i].theta)<MAXDANGLE)
							{
								//save the smallest area 
								if(((*rotate_list)[k].sigma1>
									fabs(face1[i].vec[0]*face1[face1[i].id_guest].vec[1]-
										 face1[i].vec[1]*face1[face1[i].id_guest].vec[0])*m1*n1)||
								   ((*rotate_list)[k].sigma1==
									fabs(face1[i].vec[0]*face1[face1[i].id_guest].vec[1]-
										 face1[i].vec[1]*face1[face1[i].id_guest].vec[0])*m1*n1&&
									fabs((*rotate_list)[k].theta1)>fabs(face1[i].theta)))
									{
										(*rotate_list)[k].theta1=face1[i].theta;
										(*rotate_list)[k].theta2=face2[j].theta;
										(*rotate_list)[k].sigma1=fabs(face1[i].vec[0]*face1[face1[i].id_guest].vec[1]-
																	  face1[i].vec[1]*face1[face1[i].id_guest].vec[0])*m1*n1;
										(*rotate_list)[k].sigma2=fabs(face2[j].vec[0]*face2[face2[j].id_guest].vec[1]-
																	  face2[j].vec[1]*face2[face2[j].id_guest].vec[0])*m2*n2;
										
										(*rotate_list)[k].trans1[0][0]=face1[i].vec[0]*m1;
										(*rotate_list)[k].trans1[1][0]=face1[i].vec[1]*m1;
										(*rotate_list)[k].trans1[0][1]=face1[face1[i].id_guest].vec[0]*n1;
										(*rotate_list)[k].trans1[1][1]=face1[face1[i].id_guest].vec[1]*n1;
										(*rotate_list)[k].trans1[0][2]=(*rotate_list)[k].trans1[1][2]=
										(*rotate_list)[k].trans1[2][0]=(*rotate_list)[k].trans1[2][1]=0;
										(*rotate_list)[k].trans1[2][2]=1;
										
										(*rotate_list)[k].trans2[0][0]=face2[j].vec[0]*m2;
										(*rotate_list)[k].trans2[1][0]=face2[j].vec[1]*m2;
										(*rotate_list)[k].trans2[0][1]=face2[face2[j].id_guest].vec[0]*n2;
										(*rotate_list)[k].trans2[1][1]=face2[face2[j].id_guest].vec[1]*n2;
										(*rotate_list)[k].trans2[0][2]=(*rotate_list)[k].trans2[1][2]=
										(*rotate_list)[k].trans2[2][0]=(*rotate_list)[k].trans2[2][1]=0;
										(*rotate_list)[k].trans2[2][2]=1;
											
										(*rotate_list)[k].mismatch[0]=mismatch1;
										(*rotate_list)[k].mismatch[1]=mismatch2;
										(*rotate_list)[k].mismatch[2]=fabs(face1[i].dtheta-face2[j].dtheta);
										(*rotate_list)[k].mismatch[3]=0;
										(*rotate_list)[k].mismatch[4]=0;
										
										(*rotate_list)[k].cylinder_flag=0;
										saved_flag=1;
									}
								else
									saved_flag=1;
							}
						if(saved_flag==0)
						{
							if(N_list==0)
								(*rotate_list)=(SK_scaned_data*)malloc(sizeof(SK_scaned_data));
							else 
								(*rotate_list)=(SK_scaned_data*)realloc((*rotate_list),(N_list+1)*sizeof(SK_scaned_data));
							(*rotate_list)[N_list].theta1=face1[i].theta;
							(*rotate_list)[N_list].theta2=face2[j].theta;
							(*rotate_list)[N_list].sigma1=fabs(face1[i].vec[0]*face1[face1[i].id_guest].vec[1]-
														       face1[i].vec[1]*face1[face1[i].id_guest].vec[0])*m1*n1;
							(*rotate_list)[N_list].sigma2=fabs(face2[j].vec[0]*face2[face2[j].id_guest].vec[1]-
															   face2[j].vec[1]*face2[face2[j].id_guest].vec[0])*m2*n2;
							
							(*rotate_list)[N_list].trans1[0][0]=face1[i].vec[0]*m1;
							(*rotate_list)[N_list].trans1[1][0]=face1[i].vec[1]*m1;
							(*rotate_list)[N_list].trans1[0][1]=face1[face1[i].id_guest].vec[0]*n1;
							(*rotate_list)[N_list].trans1[1][1]=face1[face1[i].id_guest].vec[1]*n1;
							(*rotate_list)[N_list].trans1[0][2]=(*rotate_list)[N_list].trans1[1][2]=
							(*rotate_list)[N_list].trans1[2][0]=(*rotate_list)[N_list].trans1[2][1]=0;
							(*rotate_list)[N_list].trans1[2][2]=1;
							
							(*rotate_list)[N_list].trans2[0][0]=face2[j].vec[0]*m2;
							(*rotate_list)[N_list].trans2[1][0]=face2[j].vec[1]*m2;
							(*rotate_list)[N_list].trans2[0][0]=face2[face2[j].id_guest].vec[0]*n2;
							(*rotate_list)[N_list].trans2[1][1]=face2[face2[j].id_guest].vec[1]*n2;
							(*rotate_list)[N_list].trans2[0][2]=(*rotate_list)[N_list].trans2[1][2]=
							(*rotate_list)[N_list].trans2[2][0]=(*rotate_list)[N_list].trans2[2][1]=0;
							(*rotate_list)[N_list].trans2[2][2]=1;
							
							(*rotate_list)[N_list].mismatch[0]=mismatch1;
							(*rotate_list)[N_list].mismatch[1]=mismatch2;
							(*rotate_list)[N_list].mismatch[2]=fabs(face1[i].dtheta-face2[j].dtheta);
							(*rotate_list)[N_list].mismatch[3]=0;
							(*rotate_list)[N_list].mismatch[4]=0;
							
							(*rotate_list)[N_list].cylinder_flag=0;
							N_list++;
						}
					}
				}
			}
		}
	}}
    return N_list;
}

bool SK_isrp(int a, int b)
{
	if(a==0)
	{
		if(b==1)
			return true;
		return false;
	}
	if(b==0)
	{
		if(a==1)
			return true;
		return false;
	}
	if(a==1||b==1)
		return true;
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
	if(b>1)	return false;
	else return true;
}

//666666666666666666666666666666666666666666666666666
//end of main
//666666666666666666666666666666666666666666666666666


int SK_install_face(double theta,double dis,int a,int b,SK_crystal_face** input,int n)
{
    if(n==0)
    {
		*input=(SK_crystal_face*)malloc(sizeof(SK_crystal_face));
		(*input)[n].dis=dis;
    	(*input)[n].theta=theta;
		(*input)[n].vec[0]=a;
		(*input)[n].vec[1]=b;
		return n+1;
    } 
    (*input)=(SK_crystal_face*)realloc(input[0],(n+1)*sizeof(SK_crystal_face));
    (*input)[n].dis=dis;
    (*input)[n].theta=theta;
	(*input)[n].vec[0]=a;
	(*input)[n].vec[1]=b;
    return n+1;
}


//     | y+ axis
//     |
//-----+------z+ axis
//     |
//     |
double SK_arc_tan(double x,double y)
{
    if(x==0)
    {
	if(y>0)return 90;
	else return -90;
    }
    if(y==0)
    {
	if(x>0)return 0;
	else return 180;
    }
    if(x>0)
    if(y>0)
    return atan(fabs(y/x))/M_PI*180.0;
    if(y>0)
    if(x<0)
    return 180.0-atan(fabs(y/x))/M_PI*180.0;
    if(y<0)
    if(x<0)
    return -180.0+atan(fabs(y/x))/M_PI*180.0;
    //if(y<0)
    //if(x>0)
    return -atan(fabs(y/x))/M_PI*180.0;
}