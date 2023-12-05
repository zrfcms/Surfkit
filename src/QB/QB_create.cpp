/*===============================================================
 * This is a part of SPaMD code
 * This program is modified at BUAA by  Z. R. Liu, March. 11 2019
 * Copyright[c] 2017-2019, zrfbuaa group
*================================================================*/
#include"QB.h"
#define BOUNDARY_DELTA 10e-6
#define MAX_DRATE 0.001

double QB_create_fcc(QB_tools *QB,double a,double M11,double M12,double M13,double M21,double M22,double M23,double M31,double M32,double M33)
{
    int i;
    double temp;
    QB_vector fcc_check[3];
    fcc_check[0].x=-1;fcc_check[0].y=1;fcc_check[0].z=1;
    fcc_check[1].x=1;fcc_check[1].y=-1;fcc_check[1].z=1;
    fcc_check[2].x=1;fcc_check[2].y=1;fcc_check[2].z=-1;
    double dangerous_rate=0;
    for(i=0;i<3;i++)
    {
	temp=fcc_check[i].x*M11+fcc_check[i].y*M12+fcc_check[i].z*M13;
	dangerous_rate+=fabs(temp-(int)temp);
	temp=fcc_check[i].x*M21+fcc_check[i].y*M22+fcc_check[i].z*M23;
	dangerous_rate+=fabs(temp-(int)temp);
	temp=fcc_check[i].x*M31+fcc_check[i].y*M32+fcc_check[i].z*M33;
	dangerous_rate+=fabs(temp-(int)temp);
    }
    dangerous_rate+=fabs(M11*M21+M12*M22+M13*M23);
    dangerous_rate+=fabs(M11*M31+M12*M32+M13*M33);
    dangerous_rate+=fabs(M31*M21+M32*M22+M33*M23);

    if(MAX_DRATE<dangerous_rate)
	QB_printf("WARNING: non-orthogonal or non-periodic system\n");

    QB_vector fcc[4];
    fcc[0].x=0;fcc[0].y=0;fcc[0].z=0;
    fcc[1].x=0.5*a;fcc[1].y=0.5*a;fcc[1].z=0;
    fcc[2].x=0.5*a;fcc[2].y=0;fcc[2].z=0.5*a;
    fcc[3].x=0;fcc[3].y=0.5*a;fcc[3].z=0.5*a;
    
    QB_vector M[3];
    temp=sqrt(M11*M11+M12*M12+M13*M13);
    M[0].x=M11/temp;
    M[1].x=M12/temp;
    M[2].x=M13/temp;
    
    temp=sqrt(M21*M21+M22*M22+M23*M23);
    M[0].y=M21/temp;
    M[1].y=M22/temp;
    M[2].y=M23/temp;

    temp=sqrt(M31*M31+M32*M32+M33*M33);
    M[0].z=M31/temp;
    M[1].z=M32/temp;
    M[2].z=M33/temp;

    QB->startx=QB->starty=QB->startz=QB->zerox=QB->zeroy=QB->zeroz=0;
    QB->mat[0][0]=QB->boundx=a*sqrt(M11*M11+M12*M12+M13*M13);
    QB->mat[1][1]=QB->boundy=a*sqrt(M21*M21+M22*M22+M23*M23);
    QB->mat[2][2]=QB->boundz=a*sqrt(M31*M31+M32*M32+M33*M33);

    temp=sqrt(QB->boundx*QB->boundx+QB->boundy*QB->boundy+QB->boundz*QB->boundz);
    int max_model=(int)(temp/a)+1;
    double x,y,z;
    double tx,ty,tz;
    for(x=-max_model;x<max_model;x++)
    for(y=-max_model;y<max_model;y++)
    for(z=-max_model;z<max_model;z++)
    for(i=0;i<4;i++)
    {
	tx=(x*a+fcc[i].x)*M[0].x+(y*a+fcc[i].y)*M[1].x+(z*a+fcc[i].z)*M[2].x;
	ty=(x*a+fcc[i].x)*M[0].y+(y*a+fcc[i].y)*M[1].y+(z*a+fcc[i].z)*M[2].y;
	tz=(x*a+fcc[i].x)*M[0].z+(y*a+fcc[i].y)*M[1].z+(z*a+fcc[i].z)*M[2].z;
	if(tx>0-BOUNDARY_DELTA&&tx<QB->boundx-BOUNDARY_DELTA)
	if(ty>0-BOUNDARY_DELTA&&ty<QB->boundy-BOUNDARY_DELTA)
	if(tz>0-BOUNDARY_DELTA&&tz<QB->boundz-BOUNDARY_DELTA)
	QB_create_atom(QB,1,tx,ty,tz);
    }
    return dangerous_rate;
}

double QB_create_bcc(QB_tools *QB,double a,double M11,double M12,double M13,double M21,double M22,double M23,double M31,double M32,double M33)
{
    int i;
    double temp;
    QB_vector bcc_check[3];
    bcc_check[0].x=0;bcc_check[0].y=1;bcc_check[0].z=1;
    bcc_check[1].x=1;bcc_check[1].y=0;bcc_check[1].z=1;
    bcc_check[2].x=1;bcc_check[2].y=1;bcc_check[2].z=0;
    double dangerous_rate=0;
    for(i=0;i<3;i++)
    {
	temp=bcc_check[i].x*M11+bcc_check[i].y*M12+bcc_check[i].z*M13;
	dangerous_rate+=fabs(temp-(int)temp);
	temp=bcc_check[i].x*M21+bcc_check[i].y*M22+bcc_check[i].z*M23;
	dangerous_rate+=fabs(temp-(int)temp);
	temp=bcc_check[i].x*M31+bcc_check[i].y*M32+bcc_check[i].z*M33;
	dangerous_rate+=fabs(temp-(int)temp);
    }
    dangerous_rate+=fabs(M11*M21+M12*M22+M13*M23);
    dangerous_rate+=fabs(M11*M31+M12*M32+M13*M33);
    dangerous_rate+=fabs(M31*M21+M32*M22+M33*M23);

    if(MAX_DRATE<dangerous_rate)
	QB_printf("WARNING: non-orthogonal or non-periodic system\n");

    QB_vector bcc[2];
    bcc[0].x=0;bcc[0].y=0;bcc[0].z=0;
    bcc[1].x=0.5*a;bcc[1].y=0.5*a;bcc[1].z=0.5*a;
    
    QB_vector M[3];
    temp=sqrt(M11*M11+M12*M12+M13*M13);
    M[0].x=M11/temp;
    M[1].x=M12/temp;
    M[2].x=M13/temp;
    
    temp=sqrt(M21*M21+M22*M22+M23*M23);
    M[0].y=M21/temp;
    M[1].y=M22/temp;
    M[2].y=M23/temp;

    temp=sqrt(M31*M31+M32*M32+M33*M33);
    M[0].z=M31/temp;
    M[1].z=M32/temp;
    M[2].z=M33/temp;

	QB->startx=QB->starty=QB->startz=QB->zerox=QB->zeroy=QB->zeroz=0;
    QB->mat[0][0]=QB->boundx=a*sqrt(M11*M11+M12*M12+M13*M13);
    QB->mat[1][1]=QB->boundy=a*sqrt(M21*M21+M22*M22+M23*M23);
    QB->mat[2][2]=QB->boundz=a*sqrt(M31*M31+M32*M32+M33*M33);

    temp=sqrt(QB->boundx*QB->boundx+QB->boundy*QB->boundy+QB->boundz*QB->boundz);
    int max_model=(int)(temp/a)+1;
    double x,y,z;
    double tx,ty,tz;
    for(x=-max_model;x<max_model;x++)
    for(y=-max_model;y<max_model;y++)
    for(z=-max_model;z<max_model;z++)
    for(i=0;i<2;i++)
    {
	tx=(x*a+bcc[i].x)*M[0].x+(y*a+bcc[i].y)*M[1].x+(z*a+bcc[i].z)*M[2].x;
	ty=(x*a+bcc[i].x)*M[0].y+(y*a+bcc[i].y)*M[1].y+(z*a+bcc[i].z)*M[2].y;
	tz=(x*a+bcc[i].x)*M[0].z+(y*a+bcc[i].y)*M[1].z+(z*a+bcc[i].z)*M[2].z;
	if(tx>0-BOUNDARY_DELTA&&tx<QB->boundx-BOUNDARY_DELTA)
	if(ty>0-BOUNDARY_DELTA&&ty<QB->boundy-BOUNDARY_DELTA)
	if(tz>0-BOUNDARY_DELTA&&tz<QB->boundz-BOUNDARY_DELTA)
	QB_create_atom(QB,1,tx,ty,tz);
    }
    return dangerous_rate;
}

double QB_create_hcp(QB_tools *QB,double a,double c,double M11,double M12,double M13,double M21,double M22,double M23,double M31,double M32,double M33)
{
    int i;
    double temp;
    M11=-0.5*M12+M11;
    M21=-0.5*M22+M21;
    M31=-0.5*M32+M31;
    double dangerous_rate=0;
    dangerous_rate+=fabs((M11*M21+0.75*M12*M22)*a*a+M13*M23*c*c);
    dangerous_rate+=fabs((M11*M31+0.75*M12*M32)*a*a+M13*M33*c*c);
    dangerous_rate+=fabs((M31*M21+0.75*M32*M22)*a*a+M33*M23*c*c);
    //dangerous_rate+=fabs(M11-(int)M11)+fabs(M12-(int)M12)+fabs(M13-(int)M13);
    //dangerous_rate+=fabs(M21-(int)M21)+fabs(M22-(int)M22)+fabs(M23-(int)M23);
    //dangerous_rate+=fabs(M31-(int)M31)+fabs(M32-(int)M32)+fabs(M33-(int)M33);
    if(MAX_DRATE<dangerous_rate)
	QB_printf("WARNING: non-orthogonal or non-periodic system\n");

    QB_vector hcp[4];
    hcp[0].x=0;hcp[0].y=0;hcp[0].z=0;
    hcp[1].x=0.5*a;hcp[1].y=0.5*sqrt(3.00)*a;hcp[1].z=0;
    hcp[2].x=0;hcp[2].y=2/sqrt(3.00)*a;hcp[2].z=0.5*c;
    hcp[3].x=0.5*a;hcp[3].y=0.5/sqrt(3.00)*a;hcp[3].z=0.5*c;
    
    QB_vector M[3];
    temp=sqrt((M11*M11+0.75*M12*M12)*a*a+M13*M13*c*c);
    M[0].x=a*M11/temp;
    M[1].x=0.5*sqrt(3.00)*a*M12/temp;
    M[2].x=c*M13/temp;
    
    temp=sqrt((M21*M21+0.75*M22*M22)*a*a+M23*M23*c*c);
    M[0].y=a*M21/temp;
    M[1].y=0.5*sqrt(3.00)*a*M22/temp;
    M[2].y=c*M23/temp;

    temp=sqrt((M31*M31+0.75*M32*M32)*a*a+M33*M33*c*c);
    M[0].z=a*M31/temp;
    M[1].z=0.5*sqrt(3.00)*a*M32/temp;
    M[2].z=c*M33/temp;


	QB->startx=QB->starty=QB->startz=QB->zerox=QB->zeroy=QB->zeroz=0;
    QB->mat[0][0]=QB->boundx=sqrt(a*a*(M11*M11+0.75*M12*M12)+c*c*M13*M13);
    QB->mat[1][1]=QB->boundy=sqrt(a*a*(M21*M21+0.75*M22*M22)+c*c*M23*M23);
    QB->mat[2][2]=QB->boundz=sqrt(a*a*(M31*M31+0.75*M32*M32)+c*c*M33*M33);

    temp=sqrt(QB->boundx*QB->boundx+QB->boundy*QB->boundy+QB->boundz*QB->boundz);
    int max_model_x=(int)(temp/a)+1;
    int max_model_y=(int)(temp/(0.5*sqrt(3.00)*a))+1;
    int max_model_z=(int)(temp/c)+1;
    double x,y,z;
    double tx,ty,tz;
    for(x=-max_model_x;x<max_model_x;x++)
    for(y=-max_model_y;y<max_model_y;y++)
    for(z=-max_model_z;z<max_model_z;z++)
    for(i=0;i<4;i++)
    {
	tx=(x*a+hcp[i].x)*M[0].x+(y*sqrt(3.00)*a+hcp[i].y)*M[1].x+(z*c+hcp[i].z)*M[2].x;
	ty=(x*a+hcp[i].x)*M[0].y+(y*sqrt(3.00)*a+hcp[i].y)*M[1].y+(z*c+hcp[i].z)*M[2].y;
	tz=(x*a+hcp[i].x)*M[0].z+(y*sqrt(3.00)*a+hcp[i].y)*M[1].z+(z*c+hcp[i].z)*M[2].z;
	if(tx>0-BOUNDARY_DELTA&&tx<QB->boundx-BOUNDARY_DELTA)
	if(ty>0-BOUNDARY_DELTA&&ty<QB->boundy-BOUNDARY_DELTA)
	if(tz>0-BOUNDARY_DELTA&&tz<QB->boundz-BOUNDARY_DELTA)
	QB_create_atom(QB,1,tx,ty,tz);
    }
    return dangerous_rate;
}

double QB_create_any(QB_tools *QB,double a,double M11,double M12,double M13,double M21,double M22,double M23,double M31,double M32,double M33,int n,double single_list[][4])
{
    int i;
    double temp;
    QB_vector check[3];
    check[0].x=0;check[0].y=1;check[0].z=1;
    check[1].x=1;check[1].y=0;check[1].z=1;
    check[2].x=1;check[2].y=1;check[2].z=0;
    double dangerous_rate=0;
    for(i=0;i<3;i++)
    {
		temp=check[i].x*M11+check[i].y*M12+check[i].z*M13;
		dangerous_rate+=fabs(temp-(int)temp);
		temp=check[i].x*M21+check[i].y*M22+check[i].z*M23;
		dangerous_rate+=fabs(temp-(int)temp);
		temp=check[i].x*M31+check[i].y*M32+check[i].z*M33;
		dangerous_rate+=fabs(temp-(int)temp);
    }
    dangerous_rate+=fabs(M11*M21+M12*M22+M13*M23);
    dangerous_rate+=fabs(M11*M31+M12*M32+M13*M33);
    dangerous_rate+=fabs(M31*M21+M32*M22+M33*M23);

    if(MAX_DRATE<dangerous_rate)
	QB_printf("WARNING: non-orthogonal or non-periodic system\n");
    
    QB_vector M[3];
    temp=sqrt(M11*M11+M12*M12+M13*M13);
    M[0].x=M11/temp;
    M[1].x=M12/temp;
    M[2].x=M13/temp;
    
    temp=sqrt(M21*M21+M22*M22+M23*M23);
    M[0].y=M21/temp;
    M[1].y=M22/temp;
    M[2].y=M23/temp;

    temp=sqrt(M31*M31+M32*M32+M33*M33);
    M[0].z=M31/temp;
    M[1].z=M32/temp;
    M[2].z=M33/temp;

    QB->startx=QB->starty=QB->startz=QB->zerox=QB->zeroy=QB->zeroz=0;
    QB->mat[0][0]=QB->boundx=a*sqrt(M11*M11+M12*M12+M13*M13);
    QB->mat[1][1]=QB->boundy=a*sqrt(M21*M21+M22*M22+M23*M23);
    QB->mat[2][2]=QB->boundz=a*sqrt(M31*M31+M32*M32+M33*M33);
	QB->TypeNumber=1;
    temp=sqrt(QB->boundx*QB->boundx+QB->boundy*QB->boundy+QB->boundz*QB->boundz);
    int max_model=(int)(temp/a)+1;
    double x,y,z;
    double tx,ty,tz;
    for(x=-max_model;x<max_model;x++)
    for(y=-max_model;y<max_model;y++)
    for(z=-max_model;z<max_model;z++)
    for(i=0;i<n;i++)
    {
		tx=(x+single_list[i][1])*a*M[0].x+(y+single_list[i][2])*a*M[1].x+(z+single_list[i][3])*a*M[2].x;
		ty=(x+single_list[i][1])*a*M[0].y+(y+single_list[i][2])*a*M[1].y+(z+single_list[i][3])*a*M[2].y;
		tz=(x+single_list[i][1])*a*M[0].z+(y+single_list[i][2])*a*M[1].z+(z+single_list[i][3])*a*M[2].z;
		if(tx>0-BOUNDARY_DELTA&&tx<QB->boundx-BOUNDARY_DELTA)
		if(ty>0-BOUNDARY_DELTA&&ty<QB->boundy-BOUNDARY_DELTA)
		if(tz>0-BOUNDARY_DELTA&&tz<QB->boundz-BOUNDARY_DELTA)
		QB_create_atom(QB,(int)single_list[i][0],tx,ty,tz);
		if(QB->TypeNumber<(int)single_list[i][0])QB->TypeNumber=(int)single_list[i][0];
    }
    return dangerous_rate;
}

