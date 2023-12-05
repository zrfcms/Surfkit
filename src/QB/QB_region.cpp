/*===============================================================
 * This is a part of QB code
 * This program is modified at BUAA by  Z. R. Liu, March. 11 2019
 * Copyright[c] 2017-2019, zrfbuaa group
 * Region system located here
*================================================================*/

#include "QB.h"

QB_tools *INV_SYS;
int QB_invest_system_flag=0;
void QB_region_set_invest_system(QB_tools*QB)
{
	QB_invest_system_flag=1;
	INV_SYS=QB;
}
int DEAL_ALL(int id,atomdata pt){return 1;}
int DEAL_NONE(int id,atomdata pt){return 0;}

double VOL_BOX(){return 1;}//(Pmax_x-Pmin_x)*(Pmax_y-Pmin_y)*(Pmax_z-Pmin_z);}
double VOL_BLANK(){return 1;}

int QB_region_sum = 0;
int QB_region_get_sum(){return QB_region_sum;}
//tip: show which region you are investigating
int QB_region_invest;
int QB_region_get_invest(){return QB_region_invest;}

REGION_CTRL* QB_region_control;

int QB_region_add(int deal(int id,atomdata pt),double vol)
{
    if(QB_region_sum==0)
        QB_region_control=(REGION_CTRL*)malloc(sizeof(REGION_CTRL));
    else
        QB_region_control=(REGION_CTRL*)realloc(QB_region_control,(QB_region_sum+1)*sizeof(REGION_CTRL));
    QB_region_control[QB_region_sum].deal=deal;
    QB_region_control[QB_region_sum].vol=vol;
    QB_region_sum++;
    return QB_region_sum-1;
}

//change region's volume
void QB_region_vol(int region_id,double vol)
{
    if(region_id<0||region_id>=QB_region_sum)return;
    QB_region_control[region_id].vol=vol;
}

void QB_region_clean()
{
    if(QB_region_sum!=0)
        free(QB_region_control);
    QB_region_sum=0;
}

int region_blank_deal(atomdata pt){return -1;}
double blank_vol(int i){return 1;}

int QB_region_logic_sum = 0;

int QB_region_extract(int id,atomdata pt,int i)
{
    if(i>=QB_region_sum)return 0;
    QB_region_invest=i;
    return QB_region_control[i].deal(id,pt);
}

double QB_region_volume(int i)
{
    QB_region_invest=i;
    return QB_region_control[i].vol;
}

REGION_MATH* QB_region_math;
int QB_region_logic(int id,atomdata pt)
{
    int i=QB_region_invest;
    int rt=0;
    if(QB_region_math[i].mode==QB_REGION_AND)
    {
        if(QB_region_extract(id,pt,QB_region_math[i].id[0])&&QB_region_extract(id,pt,QB_region_math[i].id[1]))
            rt=1;
        else rt=0;
    }
    else if(QB_region_math[i].mode==QB_REGION_OR)
    {
        if(QB_region_extract(id,pt,QB_region_math[i].id[0])||QB_region_extract(id,pt,QB_region_math[i].id[1]))
            rt=1;
        else rt=0;
    }
    else if(QB_region_math[i].mode==QB_REGION_COMPLE)
    {
        if(QB_region_extract(id,pt,QB_region_math[i].id[0])&&QB_region_extract(id,pt,QB_region_math[i].id[1]))
            rt=0;
        else if(QB_region_extract(id,pt,QB_region_math[i].id[0])||QB_region_extract(id,pt,QB_region_math[i].id[1]))
            rt=1;
        else rt=0;
    }
    else if(QB_region_math[i].mode==QB_REGION_NEITHER)
    {
        if(QB_region_extract(id,pt,QB_region_math[i].id[0])||QB_region_extract(id,pt,QB_region_math[i].id[1]))
            rt=0;
        else rt=1;
    }
    else if(QB_region_math[i].mode==QB_REGION_NOT)
    {
        if(QB_region_extract(id,pt,QB_region_math[i].id[0]))
            rt=0;
        else rt=1;
    }
    QB_region_invest=i;
    return rt;
}

void QB_region_logic_memory()
{
    if(QB_region_logic_sum==0)
    {
        QB_region_math=(REGION_MATH*)malloc(QB_region_sum*sizeof(REGION_MATH));
        QB_region_logic_sum=QB_region_sum;
        return;
    }
    if(QB_region_logic_sum<QB_region_sum)
    {
        QB_region_math=(REGION_MATH*)realloc(QB_region_math,QB_region_sum*sizeof(REGION_MATH));
        QB_region_logic_sum=QB_region_sum;
        return;
    }
    if(QB_region_logic_sum==QB_region_sum)
    {
        return;
    }
    if(QB_region_logic_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_math);
        QB_region_math=(REGION_MATH*)malloc(QB_region_sum*sizeof(REGION_MATH));
        QB_region_logic_sum=QB_region_sum;
        return;
    }
    return;
}

int QB_region_and(int a,int b,double vol)
{
	//no need for this setting
	/*
    if(a>=QB_region_sum||b>=QB_region_sum)
    {
        QB_printf("Create region failed: dealing with empty region\n");
        return -1;
    }
	*/
    int i=QB_region_add(QB_region_logic,vol);
    QB_region_logic_memory();
    QB_region_math[QB_region_sum-1].id[0]=a;
    QB_region_math[QB_region_sum-1].id[1]=b;
    QB_region_math[QB_region_sum-1].mode=QB_REGION_AND;
    return i;
}

int QB_region_or(int a,int b,double vol)
{
	//no need for this setting
	/*
    if(a>=QB_region_sum||b>=QB_region_sum)
    {
        QB_printf("Create region failed: dealing with empty region\n");
        return -1;
    }
	*/
    int i=QB_region_add(QB_region_logic,vol);
    QB_region_logic_memory();
    QB_region_math[QB_region_sum-1].id[0]=a;
    QB_region_math[QB_region_sum-1].id[1]=b;
    QB_region_math[QB_region_sum-1].mode=QB_REGION_OR;
    return i;
}

int QB_region_comple(int a,int b,double vol)
{
	//no need for this setting
	/*
    if(a>=QB_region_sum||b>=QB_region_sum)
    {
        QB_printf("Create region failed: dealing with empty region\n");
        return -1;
    }
	*/
    int i=QB_region_add(QB_region_logic,vol);
    QB_region_logic_memory();
    QB_region_math[QB_region_sum-1].id[0]=a;
    QB_region_math[QB_region_sum-1].id[1]=b;
    QB_region_math[QB_region_sum-1].mode=QB_REGION_COMPLE;
    return i;
}

int QB_region_neither(int a,int b,double vol)
{
	//no need for this setting
	/*
    if(a>=QB_region_sum||b>=QB_region_sum)
    {
        QB_printf("Create region failed: dealing with empty region\n");
        return -1;
    }
	*/
    int i=QB_region_add(QB_region_logic,vol);
    QB_region_logic_memory();
    QB_region_math[QB_region_sum-1].id[0]=a;
    QB_region_math[QB_region_sum-1].id[1]=b;
    QB_region_math[QB_region_sum-1].mode=QB_REGION_NEITHER;
    return i;
}

int QB_region_not(int a,double vol)
{
	//no need for this setting
	/*
    if(a>=QB_region_sum)
    {
        QB_printf("Create region failed: dealing with empty region\n");
        return -1;
    }
	*/
    int i=QB_region_add(QB_region_logic,vol);
    QB_region_logic_memory();
    QB_region_math[QB_region_sum-1].id[0]=a;
    QB_region_math[QB_region_sum-1].id[1]=a;
    QB_region_math[QB_region_sum-1].mode=QB_REGION_NOT;
    return i;
}


int INBOX_deal(int id,atomdata pt)
{
	if(INV_SYS->box)
	{
		QB_vector 	face[3];
		double 	endface[9];
		double x,y,z;
		int point=0,lp=0,n;
        face[0].x=INV_SYS->mat[1][1]*INV_SYS->mat[2][2]-INV_SYS->mat[1][2]*INV_SYS->mat[2][1];
        face[0].y=INV_SYS->mat[1][2]*INV_SYS->mat[2][0]-INV_SYS->mat[1][0]*INV_SYS->mat[2][2];
        face[0].z=INV_SYS->mat[1][0]*INV_SYS->mat[2][1]-INV_SYS->mat[1][1]*INV_SYS->mat[2][0];
		if(face[0].x*INV_SYS->mat[0][0]+face[0].y*INV_SYS->mat[0][1]+face[0].z*INV_SYS->mat[0][2]<0){face[0].x*=-1;face[0].y*=-1;face[0].z*=-1;}
        face[1].x=INV_SYS->mat[2][1]*INV_SYS->mat[0][2]-INV_SYS->mat[2][2]*INV_SYS->mat[0][1];
        face[1].y=INV_SYS->mat[2][2]*INV_SYS->mat[0][0]-INV_SYS->mat[2][0]*INV_SYS->mat[0][2];
        face[1].z=INV_SYS->mat[2][0]*INV_SYS->mat[0][1]-INV_SYS->mat[2][1]*INV_SYS->mat[0][0];
		if(face[1].x*INV_SYS->mat[1][0]+face[1].y*INV_SYS->mat[1][1]+face[1].z*INV_SYS->mat[1][2]<0){face[1].x*=-1;face[1].y*=-1;face[1].z*=-1;}
		face[2].x=INV_SYS->mat[0][1]*INV_SYS->mat[1][2]-INV_SYS->mat[0][2]*INV_SYS->mat[1][1];
        face[2].y=INV_SYS->mat[0][2]*INV_SYS->mat[1][0]-INV_SYS->mat[0][0]*INV_SYS->mat[1][2];
        face[2].z=INV_SYS->mat[0][0]*INV_SYS->mat[1][1]-INV_SYS->mat[0][1]*INV_SYS->mat[1][0];
		if(face[2].x*INV_SYS->mat[2][0]+face[2].y*INV_SYS->mat[2][1]+face[2].z*INV_SYS->mat[2][2]<0){face[2].x*=-1;face[2].y*=-1;face[2].z*=-1;}
    //init boundary face
    //end face 0~5 means D in Ax+By+Cz=D 6->8 means cutoff*sqrt(A*A+B*B+C*C)
        for(n=0;n<3;n++)
        {
            endface[n]=face[n].x*INV_SYS->zerox+face[n].y*INV_SYS->zeroy+face[n].z*INV_SYS->zeroz;
            endface[n+3]=face[n].x*(INV_SYS->mat[0][0]+INV_SYS->mat[1][0]+INV_SYS->mat[2][0]+INV_SYS->zerox)+
						 face[n].y*(INV_SYS->mat[0][1]+INV_SYS->mat[1][1]+INV_SYS->mat[2][1]+INV_SYS->zeroy)+
						 face[n].z*(INV_SYS->mat[0][2]+INV_SYS->mat[1][2]+INV_SYS->mat[2][2]+INV_SYS->zeroz);
        }
		x=INV_SYS->atom[id].x;
		y=INV_SYS->atom[id].y;
		z=INV_SYS->atom[id].z;
		if(x*face[0].x+y*face[0].y+z*face[0].z<endface[0])return 0;
		if(x*face[0].x+y*face[0].y+z*face[0].z>=endface[3])return 0;
		if(x*face[1].x+y*face[1].y+z*face[1].z<endface[1])return 0;
		if(x*face[1].x+y*face[1].y+z*face[1].z>=endface[4])return 0;
		if(x*face[2].x+y*face[2].y+z*face[2].z<endface[2])return 0;
		if(x*face[2].x+y*face[2].y+z*face[2].z>=endface[5])return 0; 
	}
	else 
	{
		if(pt.x>(INV_SYS->boundx+INV_SYS->startx))
			return 0;
		if(pt.y>(INV_SYS->boundy+INV_SYS->starty))
			return 0;
		if(pt.z>(INV_SYS->boundz+INV_SYS->startz))
			return 0;
		if(pt.x<INV_SYS->startx)
			return 0;
		if(pt.y<INV_SYS->starty)
			return 0;
		if(pt.z<INV_SYS->startz)
			return 0;
	}
    return 1;
}
int QB_region_inbox()
{
    int i=QB_region_add(INBOX_deal,1);
    return i;
}

//===================================================================================================================
//some packed up QB_region system
//===================================================================================================================

REGION_CUBOID* QB_region_cuboiddata;
int QB_region_cuboid_sum=0;
int CUBOID_deal(int id,atomdata pt)
{
    int i=QB_region_invest;
	if(pt.x>=QB_region_cuboiddata[i].minx)
	if(pt.x<QB_region_cuboiddata[i].maxx)
	if(pt.y>=QB_region_cuboiddata[i].miny)
	if(pt.y<QB_region_cuboiddata[i].maxy)
	if(pt.z>=QB_region_cuboiddata[i].minz)
	if(pt.z<QB_region_cuboiddata[i].maxz)
	return 1;
	return 0;
}

int QB_region_cuboid(double minx,double miny,double minz,double maxx,double maxy,double maxz)
{
    double lx=maxx-minx;
	//if(lx>Pmax_x-Pmin_x)lx=Pmax_x-Pmin_x;

	double ly=maxy-miny;
	//if(ly>Pmax_y-Pmin_y)ly=Pmax_y-Pmin_y;

	double lz=maxz-minz;
	//if(lz>Pmax_z-Pmin_z)lz=Pmax_z-Pmin_z;

    int i=QB_region_add(CUBOID_deal,lx*ly*lz);
	if(QB_region_cuboid_sum==0)
    {
        QB_region_cuboiddata=(REGION_CUBOID*)malloc(QB_region_sum*sizeof(REGION_CUBOID));
    }
    else if(QB_region_cuboid_sum<QB_region_sum)
    {
        QB_region_cuboiddata=(REGION_CUBOID*)realloc(QB_region_cuboiddata,QB_region_sum*sizeof(REGION_CUBOID));
    }
    else if(QB_region_cuboid_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_cuboiddata);
        QB_region_cuboiddata=(REGION_CUBOID*)malloc(QB_region_sum*sizeof(REGION_CUBOID));
    }
    QB_region_cuboid_sum=QB_region_sum;
    QB_region_cuboiddata[QB_region_cuboid_sum-1].minx=minx;
    QB_region_cuboiddata[QB_region_cuboid_sum-1].miny=miny;
    QB_region_cuboiddata[QB_region_cuboid_sum-1].minz=minz;
    QB_region_cuboiddata[QB_region_cuboid_sum-1].maxx=maxx;
    QB_region_cuboiddata[QB_region_cuboid_sum-1].maxy=maxy;
    QB_region_cuboiddata[QB_region_cuboid_sum-1].maxz=maxz;
    return i;
}
///////////////////////////////////////////////////////////////////////////////////////////

REGION_SPHERE* QB_region_spheredata;
int QB_region_sphere_sum=0;
int SPHERE_deal(int id,atomdata pt)
{
    int i=QB_region_invest;
	double x,y,z;
/*	
    if(fabs(pt->r.x+Pmax_x-Pmin_x-QB_region_spheredata[i].X)<fabs(pt->r.x-QB_region_spheredata[i].X))
	   	x=fabs(pt->r.x+Pmax_x-Pmin_x-QB_region_spheredata[i].X);
	else if(fabs(pt->r.x-Pmax_x+Pmin_x-QB_region_spheredata[i].X)<fabs(pt->r.x-QB_region_spheredata[i].X))
	   	x=fabs(pt->r.x-Pmax_x+Pmin_x-QB_region_spheredata[i].X);
	else x=fabs(pt->r.x-QB_region_spheredata[i].X);

	if(fabs(pt->r.y+Pmax_y-Pmin_y-QB_region_spheredata[i].Y)<fabs(pt->r.y-QB_region_spheredata[i].Y))
	   	y=fabs(pt->r.y+Pmax_y-Pmin_y-QB_region_spheredata[i].Y);
	else if(fabs(pt->r.y-Pmax_y+Pmin_y-QB_region_spheredata[i].Y)<fabs(pt->r.y-QB_region_spheredata[i].Y))
	   	y=fabs(pt->r.y-Pmax_y+Pmin_y-QB_region_spheredata[i].Y);
	else y=fabs(pt->r.y-QB_region_spheredata[i].Y);

	if(fabs(pt->r.z+Pmax_z-Pmin_z-QB_region_spheredata[i].Z)<fabs(pt->r.z-QB_region_spheredata[i].Z))
	   	z=fabs(pt->r.z+Pmax_z-Pmin_z-QB_region_spheredata[i].Z);
	else if(fabs(pt->r.z-Pmax_z+Pmin_z-QB_region_spheredata[i].Z)<fabs(pt->r.z-QB_region_spheredata[i].Z))
	   	z=fabs(pt->r.z-Pmax_z+Pmin_z-QB_region_spheredata[i].Z);
	else z=fabs(pt->r.z-QB_region_spheredata[i].Z);
*/
    x=pt.x-QB_region_spheredata[i].X;
    y=pt.y-QB_region_spheredata[i].Y;
    z=pt.z-QB_region_spheredata[i].Z;
	if(x*x+y*y+z*z<QB_region_spheredata[i].R*QB_region_spheredata[i].R)	
		return 1;
	return 0;
}

int QB_region_sphere(double cenx,double ceny,double cenz,double r)
{
    int i=QB_region_add(SPHERE_deal,4.000/3.000*M_PI*r*r*r);
	if(QB_region_sphere_sum==0)
    {
        QB_region_spheredata=(REGION_SPHERE*)malloc(QB_region_sum*sizeof(REGION_SPHERE));
    }
    else if(QB_region_sphere_sum<QB_region_sum)
    {
        QB_region_spheredata=(REGION_SPHERE*)realloc(QB_region_spheredata,QB_region_sum*sizeof(REGION_SPHERE));
    }
    else if(QB_region_sphere_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_spheredata);
        QB_region_spheredata=(REGION_SPHERE*)malloc(QB_region_sum*sizeof(REGION_SPHERE));
    }
    QB_region_sphere_sum=QB_region_sum;
    QB_region_spheredata[QB_region_sphere_sum-1].R=r;
    QB_region_spheredata[QB_region_sphere_sum-1].X=cenx;
    QB_region_spheredata[QB_region_sphere_sum-1].Y=ceny;
    QB_region_spheredata[QB_region_sphere_sum-1].Z=cenz;
    return i;
}
///////////////////////////////////////////////////////////////////////////////////////

REGION_CYLINDER* QB_region_cylinderdata;
int QB_region_cylinder_sum=0;
int CYLINDER_deal(int id,atomdata pt)
{
    int i=QB_region_invest;
	double x,y,z;
    x=pt.x-QB_region_cylinderdata[i].X;
    y=pt.y-QB_region_cylinderdata[i].Y;
    z=pt.z-QB_region_cylinderdata[i].Z;

	if(QB_region_cylinderdata[i].axis==0)
		if(y*y+z*z<QB_region_cylinderdata[i].R*QB_region_cylinderdata[i].R)
		if(fabs(x)<QB_region_cylinderdata[i].H/2.00)	
			return 1;
	if(QB_region_cylinderdata[i].axis==1)
		if(x*x+z*z<QB_region_cylinderdata[i].R*QB_region_cylinderdata[i].R)
		if(fabs(y)<QB_region_cylinderdata[i].H/2.00)	
			return 1;
	if(QB_region_cylinderdata[i].axis==2)
		if(x*x+y*y<QB_region_cylinderdata[i].R*QB_region_cylinderdata[i].R)
		if(fabs(z)<QB_region_cylinderdata[i].H/2.00)	
			return 1;
	
	
    return 0;
}

int QB_region_cylinder(double cenx,double ceny,double cenz,double r,double h,int axis)
{
    int i=QB_region_add(CYLINDER_deal, M_PI*r*r*h);
	if(QB_region_cylinder_sum==0)
    {
        QB_region_cylinderdata=(REGION_CYLINDER*)malloc(QB_region_sum*sizeof(REGION_CYLINDER));
    }
    else if(QB_region_cylinder_sum<QB_region_sum)
    {
        QB_region_cylinderdata=(REGION_CYLINDER*)realloc(QB_region_cylinderdata,QB_region_sum*sizeof(REGION_CYLINDER));
    }
    else if(QB_region_cylinder_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_cylinderdata);
        QB_region_cylinderdata=(REGION_CYLINDER*)malloc(QB_region_sum*sizeof(REGION_CYLINDER));
    }
    QB_region_cylinder_sum=QB_region_sum;
    QB_region_cylinderdata[QB_region_cylinder_sum-1].R=r;
    QB_region_cylinderdata[QB_region_cylinder_sum-1].H=h;
    QB_region_cylinderdata[QB_region_cylinder_sum-1].X=cenx;
    QB_region_cylinderdata[QB_region_cylinder_sum-1].Y=ceny;
    QB_region_cylinderdata[QB_region_cylinder_sum-1].Z=cenz;
    QB_region_cylinderdata[QB_region_cylinder_sum-1].axis=axis;
    return i;
}
///////////////////////////////////////////////////////////////////////////////////////

REGION_PLANE* QB_region_planedata;
int QB_region_plane_sum=0;
int PLANE_deal(int id,atomdata pt)
{
    int i=QB_region_invest;
	if(QB_region_planedata[i].T<=0)
	{
		if(pt.x*QB_region_planedata[i].X+pt.y*QB_region_planedata[i].Y+pt.z*QB_region_planedata[i].Z>QB_region_planedata[i].B)
			return 1;
	}
	else 
	{
		if(fabs(pt.x*QB_region_planedata[i].X+pt.y*QB_region_planedata[i].Y+pt.z*QB_region_planedata[i].Z-QB_region_planedata[i].B)
			/(QB_region_planedata[i].X*QB_region_planedata[i].X+QB_region_planedata[i].Y*QB_region_planedata[i].Y+QB_region_planedata[i].Z*QB_region_planedata[i].Z)<QB_region_planedata[i].T)
			return 1;
	}
    return 0;
}

int QB_region_plane(double x,double y,double z,double px,double py,double pz,double thickness)
{
    int i=QB_region_add(PLANE_deal,VOL_BOX());//tip:PLANE is only used to select region 
	if(QB_region_plane_sum==0)
    {
        QB_region_planedata=(REGION_PLANE*)malloc(QB_region_sum*sizeof(REGION_PLANE));
    }
    else if(QB_region_plane_sum<QB_region_sum)
    {
        QB_region_planedata=(REGION_PLANE*)realloc(QB_region_planedata,QB_region_sum*sizeof(REGION_PLANE));
    }
    else if(QB_region_plane_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_planedata);
        QB_region_planedata=(REGION_PLANE*)malloc(QB_region_sum*sizeof(REGION_PLANE));
    }
    QB_region_plane_sum=QB_region_sum;
    QB_region_planedata[QB_region_plane_sum-1].X=x;
    QB_region_planedata[QB_region_plane_sum-1].Y=y;
    QB_region_planedata[QB_region_plane_sum-1].Z=z;
    QB_region_planedata[QB_region_plane_sum-1].B=x*px+y*py+z*pz;
	QB_region_planedata[QB_region_plane_sum-1].T=thickness;
    return i;
}
///////////////////////////////////////////////////////////////////////////////////////
//tip X+ X- ... Z+ Z- => 0-5

REGION_CONE* QB_region_conedata;
int QB_region_cone_sum=0;
int CONE_deal(int id,atomdata pt)
{
    int i=QB_region_invest;
    double x,y,z;
    x=pt.x-QB_region_conedata[i].X;
    y=pt.y-QB_region_conedata[i].Y;
    z=pt.z-QB_region_conedata[i].Z;

    if(QB_region_conedata[i].DIR==0)
    {
	    if(x>0&&x<QB_region_conedata[i].H)
        if(z*z+y*y<pow(QB_region_conedata[i].R*(1-x/QB_region_conedata[i].H),2))
            return 1;
    }
    else if(QB_region_conedata[i].DIR==1)
    {
	    if(x<0&&x>-QB_region_conedata[i].H)
        if(z*z+y*y<pow(QB_region_conedata[i].R*(1+x/QB_region_conedata[i].H),2))
            return 1;
    }
    else if(QB_region_conedata[i].DIR==2)
    {
	    if(y>0&&y<QB_region_conedata[i].H)
        if(z*z+x*x<pow(QB_region_conedata[i].R*(1-y/QB_region_conedata[i].H),2))
            return 1;
    }
    else if(QB_region_conedata[i].DIR==3)
    {
	    if(y<0&&y>-QB_region_conedata[i].H)
        if(z*z+x*x<pow(QB_region_conedata[i].R*(1+y/QB_region_conedata[i].H),2))
            return 1;
    }
    else if(QB_region_conedata[i].DIR==4)
    {
	    if(z>0&&z<QB_region_conedata[i].H)
        if(x*x+y*y<pow(QB_region_conedata[i].R*(1-z/QB_region_conedata[i].H),2))
            return 1;
    }
    else if(QB_region_conedata[i].DIR==5)
    {
	    if(z<0&&z>-QB_region_conedata[i].H)
        if(x*x+y*y<pow(QB_region_conedata[i].R*(1+z/QB_region_conedata[i].H),2))
            return 1;
    }
    return 0;
}

int QB_region_cone(double x,double y,double z,double r,double h,int dir)
{
    int i=QB_region_add(CONE_deal,M_PI*r*r*h/3.000);//tip:PLANE is only used to select region 

	if(QB_region_cone_sum==0)
    {
        QB_region_conedata=(REGION_CONE*)malloc(QB_region_sum*sizeof(REGION_CONE));
    }
    else if(QB_region_cone_sum<QB_region_sum)
    {
        QB_region_conedata=(REGION_CONE*)realloc(QB_region_conedata,QB_region_sum*sizeof(REGION_CONE));
    }
    else if(QB_region_cone_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_conedata);
        QB_region_conedata=(REGION_CONE*)malloc(QB_region_sum*sizeof(REGION_CONE));
    }
    QB_region_cone_sum=QB_region_sum;
    QB_region_conedata[QB_region_cone_sum-1].X=x;
    QB_region_conedata[QB_region_cone_sum-1].Y=y;
    QB_region_conedata[QB_region_cone_sum-1].Z=z;
    QB_region_conedata[QB_region_cone_sum-1].R=r;
    QB_region_conedata[QB_region_cone_sum-1].H=h;
    QB_region_conedata[QB_region_cone_sum-1].DIR=dir;
    return i;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//special buff functions, used on regions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REGION_ROTATE* QB_region_rotatedata;
int QB_region_rotate_sum=0;
int ROTATE_deal(int id,atomdata pt)
{
    double storage_x=pt.x;
    double storage_y=pt.y;
    double storage_z=pt.z;
    int i=QB_region_invest;
	if(QB_region_rotatedata[i].L_X==0&&QB_region_rotatedata[i].L_Y==0&&QB_region_rotatedata[i].L_Z==0)
		return 0;
    pt.x-=QB_region_rotatedata[i].C_X;
    pt.y-=QB_region_rotatedata[i].C_Y;
    pt.z-=QB_region_rotatedata[i].C_Z;
	double normalize=1.0/sqrt(QB_region_rotatedata[i].L_X*QB_region_rotatedata[i].L_X+
						  QB_region_rotatedata[i].L_Y*QB_region_rotatedata[i].L_Y+
						  QB_region_rotatedata[i].L_Z*QB_region_rotatedata[i].L_Z);
    double u = QB_region_rotatedata[i].L_X*normalize;
    double v = QB_region_rotatedata[i].L_Y*normalize;
    double w = QB_region_rotatedata[i].L_Z*normalize;
    double x = pt.x;
    double y = pt.y;
    double z = pt.z;
    double cos_angle=QB_region_rotatedata[i].cos_angle;
    double sin_angle=QB_region_rotatedata[i].sin_angle;

    pt.x=x*(u*u+(v*v+w*w)*cos_angle)+y*(u*v*(1-cos_angle)+w*sin_angle)+z*(u*w*(1-cos_angle)-v*sin_angle);
    pt.y=x*(u*v*(1-cos_angle)-w*sin_angle)+y*(v*v+(u*u+w*w)*cos_angle)+z*(v*w*(1-cos_angle)+u*sin_angle);
    pt.z=x*(u*w*(1-cos_angle)+v*sin_angle)+y*(v*w*(1-cos_angle)-u*sin_angle)+z*(w*w+(u*u+v*v)*cos_angle);
    pt.x+=QB_region_rotatedata[i].C_X;
    pt.y+=QB_region_rotatedata[i].C_Y;
    pt.z+=QB_region_rotatedata[i].C_Z;
    if(QB_region_extract(id,pt,QB_region_rotatedata[i].id))
    {
        pt.x=storage_x;
        pt.y=storage_y;
        pt.z=storage_z;
	    return 1;
    }
    pt.x=storage_x;
    pt.y=storage_y;
    pt.z=storage_z;
    return 0;
}
int QB_region_rotate(int region_id,double cx,double cy,double cz,double lx,double ly,double lz,double angle)
{
    if(region_id>=QB_region_sum)
    {
        QB_printf("Create region failed: dealing with empty region\n");
        return -1;
    }
    int i=QB_region_add(ROTATE_deal,QB_region_volume(region_id));
	if(QB_region_rotate_sum==0)
    {
        QB_region_rotatedata=(REGION_ROTATE*)malloc(QB_region_sum*sizeof(REGION_ROTATE));
    }
    else if(QB_region_rotate_sum<QB_region_sum)
    {
        QB_region_rotatedata=(REGION_ROTATE*)realloc(QB_region_rotatedata,QB_region_sum*sizeof(REGION_ROTATE));
    }
    else if(QB_region_rotate_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_rotatedata);
        QB_region_rotatedata=(REGION_ROTATE*)malloc(QB_region_sum*sizeof(REGION_ROTATE));
    }
    QB_region_rotate_sum=QB_region_sum;
    QB_region_rotatedata[QB_region_rotate_sum-1].C_X=cx;
    QB_region_rotatedata[QB_region_rotate_sum-1].C_Y=cy;
    QB_region_rotatedata[QB_region_rotate_sum-1].C_Z=cz;
    QB_region_rotatedata[QB_region_rotate_sum-1].L_X=lx;  
    QB_region_rotatedata[QB_region_rotate_sum-1].L_Y=ly;
    QB_region_rotatedata[QB_region_rotate_sum-1].L_Z=lz;
    angle=angle*M_PI/180.0;
    QB_region_rotatedata[QB_region_rotate_sum-1].cos_angle=cos(angle);
    QB_region_rotatedata[QB_region_rotate_sum-1].sin_angle=sin(angle);
    QB_region_rotatedata[QB_region_rotate_sum-1].id=region_id;
    return i;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//special buff functions, used on regions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REGION_MOVE* QB_region_movedata;
int QB_region_move_sum=0;
int MOVE_deal(int id,atomdata pt)
{
	int i=QB_region_invest;
	if(QB_region_movedata[i].unit==0)
	{
		double storage_x=pt.x;
		double storage_y=pt.y;
		double storage_z=pt.z;
		pt.x-=QB_region_movedata[i].X;
		pt.y-=QB_region_movedata[i].Y;
		pt.z-=QB_region_movedata[i].Z;
		if(QB_region_extract(id,pt,QB_region_movedata[i].id))
		{
			pt.x=storage_x;
			pt.y=storage_y;
			pt.z=storage_z;
			return 1;
		}
		pt.x=storage_x;
		pt.y=storage_y;
		pt.z=storage_z;
		return 0;
	}
	else
	{
		if(!QB_invest_system_flag)
			return 0;
		double storage_x=pt.x;
		double storage_y=pt.y;
		double storage_z=pt.z;
		if(INV_SYS->box==0)
		{
			pt.x=storage_x-QB_region_movedata[i].X*(INV_SYS->boundx);
			pt.y=storage_y-QB_region_movedata[i].Y*(INV_SYS->boundy);
			pt.z=storage_z-QB_region_movedata[i].Z*(INV_SYS->boundz);
		}
		else
		{
			pt.x=storage_x-QB_region_movedata[i].X*(INV_SYS->mat[0][0])-QB_region_movedata[i].Y*(INV_SYS->mat[1][0])-QB_region_movedata[i].Z*(INV_SYS->mat[2][0]);
			pt.y=storage_y-QB_region_movedata[i].X*(INV_SYS->mat[0][1])-QB_region_movedata[i].Y*(INV_SYS->mat[1][1])-QB_region_movedata[i].Z*(INV_SYS->mat[2][1]);
			pt.z=storage_z-QB_region_movedata[i].X*(INV_SYS->mat[0][2])-QB_region_movedata[i].Y*(INV_SYS->mat[1][2])-QB_region_movedata[i].Z*(INV_SYS->mat[2][2]);
		}
		if(QB_region_extract(id,pt,QB_region_movedata[i].id))
		{
			pt.x=storage_x;
			pt.y=storage_y;
			pt.z=storage_z;
			return 1;
		}
		pt.x=storage_x;
		pt.y=storage_y;
		pt.z=storage_z;
		return 0;
	}
	return 0;
}
int QB_region_move(int region_id,int unit,double x,double y,double z)
{
    if(region_id>=QB_region_sum)
    {
        QB_printf("Create region failed: dealing with empty region\n");
        return -1;
    }
    int i=QB_region_add(MOVE_deal,QB_region_volume(region_id));
	if(QB_region_move_sum==0)
    {
        QB_region_movedata=(REGION_MOVE*)malloc(QB_region_sum*sizeof(REGION_MOVE));
    }
    else if(QB_region_move_sum<QB_region_sum)
    {
        QB_region_movedata=(REGION_MOVE*)realloc(QB_region_movedata,QB_region_sum*sizeof(REGION_MOVE));
    }
    else if(QB_region_move_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_movedata);
        QB_region_movedata=(REGION_MOVE*)malloc(QB_region_sum*sizeof(REGION_MOVE));
    }
    QB_region_move_sum=QB_region_sum;
    QB_region_movedata[QB_region_move_sum-1].X=x;
    QB_region_movedata[QB_region_move_sum-1].Y=y;
    QB_region_movedata[QB_region_move_sum-1].Z=z;
    QB_region_movedata[QB_region_move_sum-1].id=region_id;
	QB_region_movedata[QB_region_move_sum-1].unit=unit;
    return i;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//special buff functions, used on regions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REGION_DEFORM* QB_region_deformdata;
int QB_region_matrix_reverse(double in[3][3],double out[3][3])
{
	double A=in[0][0]*in[1][1]*in[2][2]+in[0][1]*in[1][2]*in[2][0]+in[0][2]*in[1][0]*in[2][1]
			-in[0][0]*in[1][2]*in[2][1]-in[0][1]*in[1][0]*in[2][2]-in[0][2]*in[1][1]*in[2][0];
	if(fabs(A)<0.0001)return 0;
	else
	{
		A=1/A;
		out[0][0]=(in[1][1]*in[2][2]-in[1][2]*in[2][1])*A;
		out[0][1]=-(in[1][0]*in[2][2]-in[1][2]*in[2][0])*A;
		out[0][2]=(in[1][0]*in[2][1]-in[1][1]*in[2][0])*A;
		out[1][0]=-(in[0][1]*in[2][2]-in[0][2]*in[2][1])*A;
		out[1][1]=(in[0][0]*in[2][2]-in[0][2]*in[2][0])*A;
		out[1][2]=-(in[0][0]*in[2][1]-in[0][1]*in[2][0])*A;
		out[2][0]=(in[0][1]*in[1][2]-in[0][2]*in[1][1])*A;		
		out[2][1]=-(in[0][0]*in[1][2]-in[0][2]*in[1][0])*A;
		out[2][2]=(in[0][0]*in[1][1]-in[0][1]*in[1][0])*A;
	}
	return 1;
}
int QB_region_deform_sum=0;
int DEFORM_deal(int id,atomdata pt)
{
	int i=QB_region_invest;
	double inv_mat[3][3];
	if(QB_region_matrix_reverse(QB_region_deformdata[i].mat,inv_mat)==0)
		return 0;
    double storage_x=pt.x;
    double storage_y=pt.y;
    double storage_z=pt.z;
	double tempx=pt.x-QB_region_deformdata[i].X;
	double tempy=pt.y-QB_region_deformdata[i].Y;
	double tempz=pt.z-QB_region_deformdata[i].Z;
    pt.x=tempx*inv_mat[0][0]+tempy*inv_mat[1][0]+tempz*inv_mat[2][0];
	pt.y=tempx*inv_mat[0][1]+tempy*inv_mat[1][1]+tempz*inv_mat[2][1];
	pt.z=tempx*inv_mat[0][2]+tempy*inv_mat[1][2]+tempz*inv_mat[2][2];
    pt.x+=QB_region_deformdata[i].X;
    pt.y+=QB_region_deformdata[i].Y;
    pt.z+=QB_region_deformdata[i].Z;
    if(QB_region_extract(id,pt,QB_region_deformdata[i].id))
    {
        pt.x=storage_x;
        pt.y=storage_y;
        pt.z=storage_z;
	    return 1;
    }
    pt.x=storage_x;
    pt.y=storage_y;
    pt.z=storage_z;
    return 0;
}
int QB_region_deform(int region_id,double x,double y,double z,double xx,double xy,double xz,double yx,double yy,double yz,double zx,double zy,double zz)
{
	/*
    if(region_id>=QB_region_sum)
    {
        QB_printf("Create region failed: dealing with empty region\n");
        return -1;
    }
	*/
    int i=QB_region_add(DEFORM_deal,1);
	if(QB_region_deform_sum==0)
    {
        QB_region_deformdata=(REGION_DEFORM*)malloc(QB_region_sum*sizeof(REGION_DEFORM));
    }
    else if(QB_region_deform_sum<QB_region_sum)
    {
        QB_region_deformdata=(REGION_DEFORM*)realloc(QB_region_deformdata,QB_region_sum*sizeof(REGION_DEFORM));
    }
    else if(QB_region_deform_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_deformdata);
        QB_region_deformdata=(REGION_DEFORM*)malloc(QB_region_sum*sizeof(REGION_DEFORM));
    }
    QB_region_deform_sum=QB_region_sum;
    QB_region_deformdata[QB_region_deform_sum-1].X=x;
    QB_region_deformdata[QB_region_deform_sum-1].Y=y;
    QB_region_deformdata[QB_region_deform_sum-1].Z=z;
	QB_region_deformdata[QB_region_deform_sum-1].mat[0][0]=xx;
	QB_region_deformdata[QB_region_deform_sum-1].mat[0][1]=xy;
	QB_region_deformdata[QB_region_deform_sum-1].mat[0][2]=xz;
	QB_region_deformdata[QB_region_deform_sum-1].mat[1][0]=yx;
	QB_region_deformdata[QB_region_deform_sum-1].mat[1][1]=yy;
	QB_region_deformdata[QB_region_deform_sum-1].mat[1][2]=yz;
	QB_region_deformdata[QB_region_deform_sum-1].mat[2][0]=zx;
	QB_region_deformdata[QB_region_deform_sum-1].mat[2][1]=zy;
	QB_region_deformdata[QB_region_deform_sum-1].mat[2][2]=zz;
    QB_region_deformdata[QB_region_deform_sum-1].id=region_id;
    return i;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//special buff functions, used on regions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REGION_PERIODIC* QB_region_periodicdata;

int QB_region_periodic_sum=0;
int PERIODIC_deal(int id,atomdata pt)
{
	if(!QB_invest_system_flag)
		return 0;
    double storage_x=pt.x;
    double storage_y=pt.y;
    double storage_z=pt.z;
    int i=QB_region_invest;
    int j,k,l;
    double u[3]={-1,0,1};
    for(j=0;j<3;j++)
    for(k=0;k<3;k++)
    for(l=0;l<3;l++)
    {
		if(INV_SYS->box==0)
		{
			pt.x=storage_x+u[j]*(INV_SYS->boundx);
			pt.y=storage_y+u[k]*(INV_SYS->boundy);
			pt.z=storage_z+u[l]*(INV_SYS->boundz);
		}
		else
		{
			pt.x=storage_x+u[j]*(INV_SYS->mat[0][0])+u[k]*(INV_SYS->mat[1][0])+u[l]*(INV_SYS->mat[2][0]);
			pt.y=storage_y+u[j]*(INV_SYS->mat[0][1])+u[k]*(INV_SYS->mat[1][1])+u[l]*(INV_SYS->mat[2][1]);
			pt.z=storage_z+u[j]*(INV_SYS->mat[0][2])+u[k]*(INV_SYS->mat[1][2])+u[l]*(INV_SYS->mat[2][2]);
		}
        if(QB_region_extract(id,pt,QB_region_periodicdata[i].id))
        {
            pt.x=storage_x;
            pt.y=storage_y;
            pt.z=storage_z;
	        return 1;
        }
    }
    pt.x=storage_x;
    pt.y=storage_y;
    pt.z=storage_z;
    return 0;
}
int QB_region_periodic(int region_id)
{
    /*
	if(region_id>=QB_region_sum)
    {
        QB_printf("Create region failed: dealing with empty region\n");
        return -1;
    }
	*/
    int i=QB_region_add(PERIODIC_deal,QB_region_volume(region_id));
	if(QB_region_periodic_sum==0)
    {
        QB_region_periodicdata=(REGION_PERIODIC*)malloc(QB_region_sum*sizeof(REGION_PERIODIC));
    }
    else if(QB_region_periodic_sum<QB_region_sum)
    {
        QB_region_periodicdata=(REGION_PERIODIC*)realloc(QB_region_periodicdata,QB_region_sum*sizeof(REGION_PERIODIC));
    }
    else if(QB_region_periodic_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_periodicdata);
        QB_region_periodicdata=(REGION_PERIODIC*)malloc(QB_region_sum*sizeof(REGION_PERIODIC));
    }
    QB_region_periodic_sum=QB_region_sum;
    QB_region_periodicdata[QB_region_periodic_sum-1].id=region_id;
    return i;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//special buff functions, used on regions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REGION_EXPRESSION* QB_region_expressiondata;
int QB_region_expression_sum=0;
int EXPRESSION_deal(int id,atomdata pt)
{
	if(!QB_invest_system_flag)
		return 0;
	int i=QB_region_invest;
	double expression=QB_get_data(INV_SYS,id,QB_region_expressiondata[i].slot);
	if(expression>=QB_region_expressiondata[i].min)
		if(expression<=QB_region_expressiondata[i].max)
			return 1;
    return 0;
}
int QB_region_expression(const char slot[1024],double min,double max)
{
    int i=QB_region_add(EXPRESSION_deal,1);
	if(QB_region_expression_sum==0)
    {
        QB_region_expressiondata=(REGION_EXPRESSION*)malloc(QB_region_sum*sizeof(REGION_EXPRESSION));
    }
    else if(QB_region_expression_sum<QB_region_sum)
    {
        QB_region_expressiondata=(REGION_EXPRESSION*)realloc(QB_region_expressiondata,QB_region_sum*sizeof(REGION_EXPRESSION));
    }
    else if(QB_region_expression_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_expressiondata);
        QB_region_expressiondata=(REGION_EXPRESSION*)malloc(QB_region_sum*sizeof(REGION_EXPRESSION));
    }
    QB_region_expression_sum=QB_region_sum;
    strcpy(QB_region_expressiondata[QB_region_expression_sum-1].slot,slot);
    QB_region_expressiondata[QB_region_expression_sum-1].min=min;
	QB_region_expressiondata[QB_region_expression_sum-1].max=max;
    return i;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//special buff functions, used on regions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

REGION_POLYHEDRON* QB_region_polyhedrondata;
int QB_region_polyhedrondata_sum=0;
int POLYHEDRON_deal(int id,atomdata pt)
{
	int i=QB_region_invest;
	if(QB_region_polyhedrondata[i].dir_x==0&&QB_region_polyhedrondata[i].dir_y==0&&QB_region_polyhedrondata[i].dir_z==0)
		return 0;
	static double cut[144][3];
	static int cut_n=0;
	int j,k,l,m;
	double u[2]={-1,1};
	double temp;
	double dir_x;
	double dir_y;
	double dir_z;
	if(QB_region_polyhedrondata[i].mode==0)
	{
		cut_n=0;
		for(m=0;m<QB_region_polyhedrondata[i].num+1;m++)
		{
			temp=sqrt(QB_squaredlenths(QB_region_polyhedrondata[i].dir_x[m],QB_region_polyhedrondata[i].dir_y[m],QB_region_polyhedrondata[i].dir_z[m]));
			dir_x=QB_region_polyhedrondata[i].dir_x[m]/temp;
			dir_y=QB_region_polyhedrondata[i].dir_y[m]/temp;
			dir_z=QB_region_polyhedrondata[i].dir_z[m]/temp;
			u[0]=-1.0/QB_region_polyhedrondata[i].distance[m];
			u[1]= 1.0/QB_region_polyhedrondata[i].distance[m];
			for(j=0;j<2;j++)
			for(k=0;k<2;k++)
			for(l=0;l<2;l++)
			{
				cut[cut_n][0]=u[j]*dir_x;
				cut[cut_n][1]=u[k]*dir_y;
				cut[cut_n][2]=u[l]*dir_z;
				cut_n++;
			}
		}
	}
	else if(QB_region_polyhedrondata[i].mode==1)
	{
		cut_n=0;
		for(m=0;m<QB_region_polyhedrondata[i].num+1;m++)
		{
			temp=sqrt(QB_squaredlenths(QB_region_polyhedrondata[i].dir_x[m],QB_region_polyhedrondata[i].dir_y[m],QB_region_polyhedrondata[i].dir_z[m]));
			dir_x=QB_region_polyhedrondata[i].dir_x[m]/temp;
			dir_y=QB_region_polyhedrondata[i].dir_y[m]/temp;
			dir_z=QB_region_polyhedrondata[i].dir_z[m]/temp;
			u[0]=-1.0/QB_region_polyhedrondata[i].distance[m];
			u[1]= 1.0/QB_region_polyhedrondata[i].distance[m];
			for(j=0;j<2;j++)
			for(k=0;k<2;k++)
			for(l=0;l<2;l++)
			{
				cut[cut_n][0]=u[j]*dir_x;
				cut[cut_n][1]=u[k]*dir_y;
				cut[cut_n][2]=u[l]*dir_z;
				cut_n++;
				cut[cut_n][0]=u[j]*dir_x;
				cut[cut_n][1]=u[k]*dir_z;
				cut[cut_n][2]=u[l]*dir_y;
				cut_n++;
				cut[cut_n][0]=u[j]*dir_y;
				cut[cut_n][1]=u[k]*dir_x;
				cut[cut_n][2]=u[l]*dir_z;
				cut_n++;
				cut[cut_n][0]=u[j]*dir_y;
				cut[cut_n][1]=u[k]*dir_z;
				cut[cut_n][2]=u[l]*dir_x;
				cut_n++;
				cut[cut_n][0]=u[j]*dir_z;
				cut[cut_n][1]=u[k]*dir_x;
				cut[cut_n][2]=u[l]*dir_y;
				cut_n++;
				cut[cut_n][0]=u[j]*dir_z;
				cut[cut_n][1]=u[k]*dir_y;
				cut[cut_n][2]=u[l]*dir_x;
				cut_n++;
			}
		}
	}
	else if(QB_region_polyhedrondata[i].mode==2)
	{
		cut_n=4;
		cut[0][0]=0;				cut[0][1]=0;			cut[0][2]=1;
		cut[1][0]=-2/(3*sqrt(2)); 	cut[1][1]=-2/sqrt(6);	cut[1][2]=-1.0/3.0;
		cut[2][0]=-2/(3*sqrt(2)); 	cut[2][1]=2/sqrt(6);	cut[2][2]=-1.0/3.0;
		cut[3][0]=4/(3*sqrt(2)); 	cut[3][1]=0;			cut[3][2]=-1.0/3.0;
	}
	else if(QB_region_polyhedrondata[i].mode==3)
	{		
		cut_n=6;
		cut[0][0]=0;	cut[0][1]=0;	cut[0][2]=1;
		cut[1][0]=0;	cut[1][1]=0;	cut[1][2]=-1;
		cut[2][0]=0;	cut[2][1]=1;	cut[2][2]=0;
		cut[3][0]=0;	cut[3][1]=-1;	cut[3][2]=0;
		cut[4][0]=1;	cut[4][1]=0;	cut[4][2]=0;
		cut[5][0]=-1;	cut[5][1]=0;	cut[5][2]=0;
	}
	else if(QB_region_polyhedrondata[i].mode==4)
	{
		cut_n=8;
		cut[0][0]=sqrt(3)/3;	cut[0][1]=sqrt(3)/3;	cut[0][2]=sqrt(3)/3;
		cut[1][0]=sqrt(3)/3;	cut[1][1]=sqrt(3)/3;	cut[1][2]=-sqrt(3)/3;
		cut[2][0]=sqrt(3)/3;	cut[2][1]=-sqrt(3)/3;	cut[2][2]=sqrt(3)/3;
		cut[3][0]=sqrt(3)/3;	cut[3][1]=-sqrt(3)/3;	cut[3][2]=-sqrt(3)/3;
		cut[4][0]=-sqrt(3)/3;	cut[4][1]=sqrt(3)/3;	cut[4][2]=sqrt(3)/3;
		cut[5][0]=-sqrt(3)/3;	cut[5][1]=sqrt(3)/3;	cut[5][2]=-sqrt(3)/3;
		cut[6][0]=-sqrt(3)/3;	cut[6][1]=-sqrt(3)/3;	cut[6][2]=sqrt(3)/3;
		cut[7][0]=-sqrt(3)/3;	cut[7][1]=-sqrt(3)/3;	cut[7][2]=-sqrt(3)/3;
	}
	else if(QB_region_polyhedrondata[i].mode==5)
	{
		cut_n=12;
		double mm=sqrt(50-10*sqrt(5))/10;
		double nn=sqrt(50+10*sqrt(5))/10;
		cut[0][0]=mm;	cut[0][1]=0;	cut[0][2]=nn;
		cut[1][0]=mm;	cut[1][1]=0;	cut[1][2]=-nn;
		cut[2][0]=-mm;	cut[2][1]=0;	cut[2][2]=nn;
		cut[3][0]=-mm;	cut[3][1]=0;	cut[3][2]=-nn;
		cut[4][0]=0;	cut[4][1]=nn;	cut[4][2]=mm;
		cut[5][0]=0;	cut[5][1]=nn;	cut[5][2]=-mm;
		cut[6][0]=0;	cut[6][1]=-nn;	cut[6][2]=mm;
		cut[7][0]=0;	cut[7][1]=-nn;	cut[7][2]=-mm;
		cut[8][0]=nn;	cut[8][1]=mm;	cut[8][2]=0;
		cut[9][0]=nn;	cut[9][1]=-mm;	cut[9][2]=0;
		cut[10][0]=-nn;	cut[10][1]=mm;	cut[10][2]=0;
		cut[11][0]=-nn;	cut[11][1]=-mm;	cut[11][2]=0;
	}
	else if(QB_region_polyhedrondata[i].mode==6)
	{
		cut_n=20;
		double mm=1/(0.5*(sqrt(5)+1))/sqrt(3);
		double nn=1*0.5*(sqrt(5)+1)/sqrt(3);
		cut[0][0]=mm;	cut[0][1]=0;	cut[0][2]=nn;
		cut[1][0]=mm;	cut[1][1]=0;	cut[1][2]=-nn;
		cut[2][0]=-mm;	cut[2][1]=0;	cut[2][2]=nn;
		cut[3][0]=-mm;	cut[3][1]=0;	cut[3][2]=-nn;
		cut[4][0]=0;	cut[4][1]=nn;	cut[4][2]=mm;
		cut[5][0]=0;	cut[5][1]=nn;	cut[5][2]=-mm;
		cut[6][0]=0;	cut[6][1]=-nn;	cut[6][2]=mm;
		cut[7][0]=0;	cut[7][1]=-nn;	cut[7][2]=-mm;
		cut[8][0]=nn;	cut[8][1]=mm;	cut[8][2]=0;
		cut[9][0]=nn;	cut[9][1]=-mm;	cut[9][2]=0;
		cut[10][0]=-nn;	cut[10][1]=mm;	cut[10][2]=0;
		cut[11][0]=-nn;	cut[11][1]=-mm;	cut[11][2]=0;
		cut[12][0]=sqrt(3)/3;	cut[12][0]=sqrt(3)/3;	cut[12][0]=sqrt(3)/3;
		cut[13][0]=sqrt(3)/3;	cut[13][0]=sqrt(3)/3;	cut[13][0]=-sqrt(3)/3;
		cut[14][0]=sqrt(3)/3;	cut[14][0]=-sqrt(3)/3;	cut[14][0]=sqrt(3)/3;
		cut[15][0]=sqrt(3)/3;	cut[15][0]=-sqrt(3)/3;	cut[15][0]=-sqrt(3)/3;
		cut[16][0]=-sqrt(3)/3;	cut[16][0]=sqrt(3)/3;	cut[16][0]=sqrt(3)/3;
		cut[17][0]=-sqrt(3)/3;	cut[17][0]=sqrt(3)/3;	cut[17][0]=-sqrt(3)/3;
		cut[18][0]=-sqrt(3)/3;	cut[18][0]=-sqrt(3)/3;	cut[18][0]=sqrt(3)/3;
		cut[19][0]=-sqrt(3)/3;	cut[19][0]=-sqrt(3)/3;	cut[19][0]=-sqrt(3)/3;
	}
	if(QB_region_polyhedrondata[i].mode<2)
	{
		for(j=0;j<cut_n;j++)
		{
			if((pt.x-QB_region_polyhedrondata[i].x)*cut[j][0]+
			   (pt.y-QB_region_polyhedrondata[i].y)*cut[j][1]+
			   (pt.z-QB_region_polyhedrondata[i].z)*cut[j][2]>1)
					return 0;
		}
	}
	else
	{
		for(j=0;j<cut_n;j++)
		{
			if((pt.x-QB_region_polyhedrondata[i].x)*cut[j][0]+
			   (pt.y-QB_region_polyhedrondata[i].y)*cut[j][1]+
			   (pt.z-QB_region_polyhedrondata[i].z)*cut[j][2]>QB_region_polyhedrondata[i].distance[0])
					return 0;
		}
	}
    return 1;
}
int QB_region_polyhedron(int mode,int num,
			double x,double y,double z,
			double dx1,double dy1,double dz1,
			double dx2,double dy2,double dz2,
			double dx3,double dy3,double dz3,
			double distance1,double distance2,double distance3)
{
    int i=QB_region_add(POLYHEDRON_deal,1);
	if(QB_region_polyhedrondata_sum==0)
    {
        QB_region_polyhedrondata=(REGION_POLYHEDRON*)malloc(QB_region_sum*sizeof(REGION_POLYHEDRON));
    }
    else if(QB_region_polyhedrondata_sum<QB_region_sum)
    {
        QB_region_polyhedrondata=(REGION_POLYHEDRON*)realloc(QB_region_polyhedrondata,QB_region_sum*sizeof(REGION_POLYHEDRON));
    }
    else if(QB_region_polyhedrondata_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_polyhedrondata);
        QB_region_polyhedrondata=(REGION_POLYHEDRON*)malloc(QB_region_sum*sizeof(REGION_POLYHEDRON));
    }
    QB_region_polyhedrondata_sum=QB_region_sum;
    QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].mode=mode;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].num=num;
    QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].x=x;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].y=y;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].z=z;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].dir_x[0]=dx1;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].dir_y[0]=dy1;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].dir_z[0]=dz1;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].dir_x[1]=dx2;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].dir_y[1]=dy2;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].dir_z[1]=dz2;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].dir_x[2]=dx3;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].dir_y[2]=dy3;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].dir_z[2]=dz3;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].distance[0]=distance1;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].distance[1]=distance2;
	QB_region_polyhedrondata[QB_region_polyhedrondata_sum-1].distance[2]=distance3;
    return i;
}

REGION_PRISM* QB_region_prismdata;
int QB_region_prismdata_sum=0;

int PRISM_deal(int id,atomdata pt)
{
	int i=QB_region_invest;
	int j;
	double theta=M_PI*2.0/QB_region_prismdata[i].n;
    double x,y,z;
	x=pt.x-QB_region_prismdata[i].x;
	y=pt.y-QB_region_prismdata[i].y;
	z=pt.z-QB_region_prismdata[i].z;
	
	if(QB_region_prismdata[i].mode==0)
	{
		if(x<-0.5*QB_region_prismdata[i].h||x>0.5*QB_region_prismdata[i].h)
			return 0;
	}
	if(QB_region_prismdata[i].mode==1)
	{
		if(y<-0.5*QB_region_prismdata[i].h||y>0.5*QB_region_prismdata[i].h)
			return 0;
	}
	if(QB_region_prismdata[i].mode==2)
	{
		if(z<-0.5*QB_region_prismdata[i].h||z>0.5*QB_region_prismdata[i].h)
			return 0;
	}
	if(QB_region_prismdata[i].mode==0)
	for(j=0;j<QB_region_prismdata[i].n;j++)
	{
		if(y*cos(j*theta)+z*sin(j*theta)>QB_region_prismdata[i].distance)
			return 0;
	}
	if(QB_region_prismdata[i].mode==1)
	for(j=0;j<QB_region_prismdata[i].n;j++)
	{
		if(x*cos(j*theta)+z*sin(j*theta)>QB_region_prismdata[i].distance)
			return 0;
	}
	if(QB_region_prismdata[i].mode==2)
	for(j=0;j<QB_region_prismdata[i].n;j++)
	{
		if(x*cos(j*theta)+y*sin(j*theta)>QB_region_prismdata[i].distance)
			return 0;
	}
    return 1;
}

int QB_region_prism(double x,double y,double z,double distance,double h,int n,int axis)
{
    int i=QB_region_add(PRISM_deal,1);
	if(QB_region_prismdata_sum==0)
    {
        QB_region_prismdata=(REGION_PRISM*)malloc(QB_region_sum*sizeof(REGION_PRISM));
    }
    else if(QB_region_prismdata_sum<QB_region_sum)
    {
        QB_region_prismdata=(REGION_PRISM*)realloc(QB_region_prismdata,QB_region_sum*sizeof(REGION_PRISM));
    }
    else if(QB_region_prismdata_sum>QB_region_sum)//it seems that sum of region reset once, so we reset too
    {
        free(QB_region_prismdata);
        QB_region_prismdata=(REGION_PRISM*)malloc(QB_region_sum*sizeof(REGION_PRISM));
    }
    QB_region_prismdata_sum=QB_region_sum;
    QB_region_prismdata[QB_region_prismdata_sum-1].mode=axis;
    QB_region_prismdata[QB_region_prismdata_sum-1].x=x;
	QB_region_prismdata[QB_region_prismdata_sum-1].y=y;
	QB_region_prismdata[QB_region_prismdata_sum-1].z=z;
	QB_region_prismdata[QB_region_prismdata_sum-1].distance=distance;
	QB_region_prismdata[QB_region_prismdata_sum-1].h=h;
	QB_region_prismdata[QB_region_prismdata_sum-1].n=n;
    return i;
}

void QB_region_free()
{
	if(QB_region_sum!=0)
		free(QB_region_control);
	QB_region_sum=0;
	if(QB_region_logic_sum!=0)
		free(QB_region_math);
	QB_region_logic_sum=0;
	if(QB_region_cuboid_sum!=0)
		free(QB_region_cuboiddata);
	QB_region_cuboid_sum=0;
	if(QB_region_sphere_sum!=0)
		free(QB_region_spheredata);
	QB_region_sphere_sum=0;
	if(QB_region_cylinder_sum!=0)
		free(QB_region_cylinderdata);
	QB_region_cylinder_sum=0;
	if(QB_region_plane_sum!=0)
		free(QB_region_planedata);
	QB_region_plane_sum=0;
	if(QB_region_cone_sum!=0)
		free(QB_region_conedata);
	QB_region_cone_sum=0;
	if(QB_region_rotate_sum!=0)
		free(QB_region_rotatedata);
	QB_region_rotate_sum=0;
	if(QB_region_move_sum!=0)
		free(QB_region_movedata);
	QB_region_move_sum=0;
	if(QB_region_deform_sum!=0)
		free(QB_region_deformdata);
	QB_region_deform_sum=0;
	if(QB_region_periodic_sum!=0)
		free(QB_region_periodicdata);
	QB_region_periodic_sum=0;
	if(QB_region_expression_sum!=0)
		free(QB_region_expressiondata);
	QB_region_expression_sum=0;
	if(QB_region_polyhedrondata_sum!=0)
		free(QB_region_polyhedrondata);
	QB_region_polyhedrondata_sum=0;
	if(QB_region_prismdata_sum!=0)
		free(QB_region_prismdata);
	QB_region_prismdata_sum=0;
}
