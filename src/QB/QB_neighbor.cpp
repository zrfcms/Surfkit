/*===============================================================
 * This is a part of SPaMD code
 * This program is modified at BUAA by  Z. R. Liu, March. 11 2019
 * Copyright[c] 2017-2019, zrfbuaa group
*================================================================*/
#include"QB.h"
// a limiter working on maxium coordination number
// in case of too many neighbors slow down program from running
void QB_list_build(QB_tools *QB,int i,double cutoff)
{
    if(QB->list_flag!=QB->MCN)
	{
		QB_list_clear(QB);
		QB->neb.id=(int*)malloc(QB->MCN*sizeof(int));
		QB->neb.dis=(double*)malloc(QB->MCN*sizeof(double));
	}
    int     j,k,x0,y0,z0,p;//used to save id and int atom position
    double  dis; //used to save atom distance
    int     template_int;
    double  template_double;//save temporary data
    int     x,y,z;//coordinate data in network
    QB->neb.num=0;
    x=(QB->atom[i].x-QB->startx+cutoff)/QB->tbx;
    x=QB_Min(QB_Max(x,0),QB->xlot-1);
    y=(QB->atom[i].y-QB->starty+cutoff)/QB->tby;
    y=QB_Min(QB_Max(y,0),QB->ylot-1);
    z=(QB->atom[i].z-QB->startz+cutoff)/QB->tbz;
    z=QB_Min(QB_Max(z,0),QB->zlot-1);
	
//set max distance for each neighbor slot

    for(x0=QB_Max(0,x-1);x0<QB_Min(QB->xlot,x+2);x0++)
        for(y0=QB_Max(0,y-1);y0<QB_Min(QB->ylot,y+2);y0++)
            for(z0=QB_Max(0,z-1);z0<QB_Min(QB->zlot,z+2);z0++)
                for(p=0;p<QB->network_num[x0][y0][z0];p++)
                {
                    j=QB->network_ID[x0][y0][z0][p];
                    if(j!=i)
                    {
                        dis=QB_squaredlenths(QB->atom[i].x-QB->atom[j].x,
                        		  QB->atom[i].y-QB->atom[j].y,
                        		  QB->atom[i].z-QB->atom[j].z);
                        ///====================================
                        /// sort neighbor from small to big
                        ///====================================
                        if(QB->neb.num<QB->MCN)
                        {
                            if(dis<cutoff*cutoff)
                            {
                                QB->neb.num++;
                                QB->neb.id[QB->neb.num-1]=j;
                                QB->neb.dis[QB->neb.num-1]=dis;

                                for(k=QB->neb.num-1;k>0;k--)
                                    if(QB->neb.dis[k]<QB->neb.dis[k-1])
                                    {
                                        template_double=QB->neb.dis[k];
                                        QB->neb.dis[k]=QB->neb.dis[k-1];
                                        QB->neb.dis[k-1]=template_double;
                                        template_int=QB->neb.id[k];
                                        QB->neb.id[k]=QB->neb.id[k-1];
                                        QB->neb.id[k-1]=template_int;
                                    }
                            }
                        }
                        else if(dis<QB->neb.dis[QB->neb.num-1])
                        {
                            QB->neb.id[QB->neb.num-1]=j;
                            QB->neb.dis[QB->neb.num-1]=dis;
                            for(k=QB->neb.num-1;k>0;k--)
                                if(QB->neb.dis[k]<QB->neb.dis[k-1])
                                {
                                    template_double=QB->neb.dis[k];
                                    QB->neb.dis[k]=QB->neb.dis[k-1];
                                    QB->neb.dis[k-1]=template_double;
                                    template_int=QB->neb.id[k];
                                    QB->neb.id[k]=QB->neb.id[k-1];
                                    QB->neb.id[k-1]=template_int;
                                }
                        }
                    }
                }
        for(k=0;k<QB->neb.num;k++)
              	QB->neb.dis[k]=sqrt(QB->neb.dis[k]);
	QB->list_flag=1;
}

void QB_nearest_atom(QB_tools *QB,double r_x,double r_y,double r_z,double cutoff)
{
    if(QB->list_flag!=QB->MCN)
	{
		QB_list_clear(QB);
		QB->neb.id=(int*)malloc(QB->MCN*sizeof(int));
		QB->neb.dis=(double*)malloc(QB->MCN*sizeof(double));
	}
    int     j,k,x0,y0,z0,p;//used to save id and int atom position
    double  dis; //used to save atom distance
    int     template_int;
    double  template_double;//save temporary data
    int     x,y,z;//coordinate data in network
    QB->neb.num=0;
    x=(r_x-QB->startx+cutoff)/QB->tbx;
    x=QB_Min(QB_Max(x,0),QB->xlot-1);
    y=(r_y-QB->starty+cutoff)/QB->tby;
    y=QB_Min(QB_Max(y,0),QB->ylot-1);
    z=(r_z-QB->startz+cutoff)/QB->tbz;
    z=QB_Min(QB_Max(z,0),QB->zlot-1);
	
//set max distance for each neighbor slot

    for(x0=QB_Max(0,x-1);x0<QB_Min(QB->xlot,x+2);x0++)
        for(y0=QB_Max(0,y-1);y0<QB_Min(QB->ylot,y+2);y0++)
            for(z0=QB_Max(0,z-1);z0<QB_Min(QB->zlot,z+2);z0++)
                for(p=0;p<QB->network_num[x0][y0][z0];p++)
                {
                    j=QB->network_ID[x0][y0][z0][p];
                    dis=QB_squaredlenths(r_x-QB->atom[j].x,
										r_y-QB->atom[j].y,
										r_z-QB->atom[j].z);
                    ///====================================
                    /// sort neighbor from small to big
                    ///====================================
                    if(QB->neb.num<QB->MCN)
                    {
                        if(dis<cutoff*cutoff)
                        {
                            QB->neb.num++;
                            QB->neb.id[QB->neb.num-1]=j;
                            QB->neb.dis[QB->neb.num-1]=dis;

                            for(k=QB->neb.num-1;k>0;k--)
                                if(QB->neb.dis[k]<QB->neb.dis[k-1])
                                {
                                    template_double=QB->neb.dis[k];
                                    QB->neb.dis[k]=QB->neb.dis[k-1];
                                    QB->neb.dis[k-1]=template_double;
                                    template_int=QB->neb.id[k];
                                    QB->neb.id[k]=QB->neb.id[k-1];
                                    QB->neb.id[k-1]=template_int;
                                }
                        }
                    }
                    else if(dis<QB->neb.dis[QB->neb.num-1])
                    {
                        QB->neb.id[QB->neb.num-1]=j;
                        QB->neb.dis[QB->neb.num-1]=dis;
                        for(k=QB->neb.num-1;k>0;k--)
                            if(QB->neb.dis[k]<QB->neb.dis[k-1])
                            {
                                template_double=QB->neb.dis[k];
                                QB->neb.dis[k]=QB->neb.dis[k-1];
                                QB->neb.dis[k-1]=template_double;
                                template_int=QB->neb.id[k];
                                QB->neb.id[k]=QB->neb.id[k-1];
                                QB->neb.id[k-1]=template_int;
                            }
                    }
                }
    for(k=0;k<QB->neb.num;k++)
      	QB->neb.dis[k]=sqrt(QB->neb.dis[k]);
    QB->list_flag=1;
}


void QB_pbclist_build(QB_tools *QB,int i,double cutoff)
{
    if(QB->list_flag!=QB->MCN)
	{
		QB_list_clear(QB);
		QB->neb.id=(int*)malloc(QB->MCN*sizeof(int));
		QB->neb.dis=(double*)malloc(QB->MCN*sizeof(double));
	}
    int     j,k,x0,y0,z0,p,q,l,m,n,judge;//used to save id and int atom position
    double  dis; //used to save atom distance
    int     template_int;
    double  template_double;//save temporary data
    int     x,y,z;//coordinate data in network
	double  x1,y1,z1;
	double backup_x=QB->atom[i].x;
	double backup_y=QB->atom[i].y;
	double backup_z=QB->atom[i].z;
    QB->neb.num=0;
    double u[3]={0,1,-1};
	QB_vector 	face[3];	
	double 	endface[9];
//set max distance for each neighbor slot
	if(QB->box==1)
    {
        face[0].x=QB->mat[1][1]*QB->mat[2][2]-QB->mat[1][2]*QB->mat[2][1];
        face[0].y=QB->mat[1][2]*QB->mat[2][0]-QB->mat[1][0]*QB->mat[2][2];
        face[0].z=QB->mat[1][0]*QB->mat[2][1]-QB->mat[1][1]*QB->mat[2][0];
        face[1].x=QB->mat[2][1]*QB->mat[0][2]-QB->mat[2][2]*QB->mat[0][1];
        face[1].y=QB->mat[2][2]*QB->mat[0][0]-QB->mat[2][0]*QB->mat[0][2];
        face[1].z=QB->mat[2][0]*QB->mat[0][1]-QB->mat[2][1]*QB->mat[0][0];
		face[2].x=QB->mat[0][1]*QB->mat[1][2]-QB->mat[0][2]*QB->mat[1][1];
        face[2].y=QB->mat[0][2]*QB->mat[1][0]-QB->mat[0][0]*QB->mat[1][2];
        face[2].z=QB->mat[0][0]*QB->mat[1][1]-QB->mat[0][1]*QB->mat[1][0];
    //init boundary face
    //end face 0~5 means D in Ax+By+Cz=D 6->8 means cutoff*sqrt(A*A+B*B+C*C)
        for(n=0;n<3;n++)
        {
            endface[n]=face[n].x*QB->zerox+face[n].y*QB->zeroy+face[n].z*QB->zeroz;
            endface[n+3]=face[n].x*(QB->mat[0][0]+QB->mat[1][0]+QB->mat[2][0]+QB->zerox)+
						 face[n].y*(QB->mat[0][1]+QB->mat[1][1]+QB->mat[2][1]+QB->zeroy)+
						 face[n].z*(QB->mat[0][2]+QB->mat[1][2]+QB->mat[2][2]+QB->zeroz);
            endface[n+6]=sqrt(QB_squaredlenths(face[n].x,face[n].y,face[n].z))*cutoff;
        }
		x1=QB->atom[i].x*face[0].x+QB->atom[i].y*face[0].y+QB->atom[i].z*face[0].z;
        y1=QB->atom[i].x*face[1].x+QB->atom[i].y*face[1].y+QB->atom[i].z*face[1].z;
        z1=QB->atom[i].x*face[2].x+QB->atom[i].y*face[2].y+QB->atom[i].z*face[2].z;
	}
	else
	{
		x1=QB->atom[i].x;
		y1=QB->atom[i].y;
		z1=QB->atom[i].z;
		
	}
	
    for(k=0;k<=2*QB->px;k++)
        for(p=0;p<=2*QB->py;p++)
            for(q=0;q<=2*QB->pz;q++)
			{
				QB->atom[i].x=backup_x;
				QB->atom[i].y=backup_y;
				QB->atom[i].z=backup_z;
				if(QB->box==1)
				{
                    judge=0;
                    if(k==2&&fabs(x1-endface[3])<endface[6])
                        judge+=1;
                    if(p==2&&fabs(y1-endface[4])<endface[7])
                        judge+=1;
                    if(q==2&&fabs(z1-endface[5])<endface[8])
                        judge+=1;
                    if(k==1&&fabs(x1-endface[0])<endface[6])
                        judge+=1;
                    if(p==1&&fabs(y1-endface[1])<endface[7])
                        judge+=1;
                    if(q==1&&fabs(z1-endface[2])<endface[8])
                        judge+=1;
                    if(k==0)
                        judge+=1;
                    if(p==0)
                        judge+=1;
                    if(q==0)
                        judge+=1;
                    if(judge==3)
                    {
                        QB->atom[i].x=QB->atom[i].x+u[k]*QB->mat[0][0]+u[p]*QB->mat[1][0]+u[q]*QB->mat[2][0];
                        QB->atom[i].y=QB->atom[i].y+u[k]*QB->mat[0][1]+u[p]*QB->mat[1][1]+u[q]*QB->mat[2][1];
                        QB->atom[i].z=QB->atom[i].z+u[k]*QB->mat[0][2]+u[p]*QB->mat[1][2]+u[q]*QB->mat[2][2];
                    }
					else continue;
                }
				else
				{
					judge=0;
					if(k==2&&x1>QB->startx+QB->boundx-cutoff)
						judge+=1;
					if(p==2&&y1>QB->starty+QB->boundy-cutoff)
						judge+=1;
					if(q==2&&z1>QB->startz+QB->boundz-cutoff)
						judge+=1;
					if(k==1&&x1<QB->startx+cutoff)
						judge+=1;
					if(p==1&&y1<QB->starty+cutoff)
						judge+=1;
					if(q==1&&z1<QB->startz+cutoff)
						judge+=1;
					if(k==0)
						judge+=1;
					if(p==0)
						judge+=1;
					if(q==0)
						judge+=1;
				    if(judge==3)
					{
						QB->atom[i].x+=u[k]*QB->boundx;
						QB->atom[i].y+=u[p]*QB->boundy;
						QB->atom[i].z+=u[q]*QB->boundz;
					}
					else continue;
				}
				x=(QB->atom[i].x-QB->startx+cutoff)/QB->tbx;
				x=QB_Min(QB_Max(x,0),QB->xlot-1);
				y=(QB->atom[i].y-QB->starty+cutoff)/QB->tby;
				y=QB_Min(QB_Max(y,0),QB->ylot-1);
				z=(QB->atom[i].z-QB->startz+cutoff)/QB->tbz;
				z=QB_Min(QB_Max(z,0),QB->zlot-1);
				for(x0=QB_Max(0,x-1);x0<QB_Min(QB->xlot,x+2);x0++)
					for(y0=QB_Max(0,y-1);y0<QB_Min(QB->ylot,y+2);y0++)
						for(z0=QB_Max(0,z-1);z0<QB_Min(QB->zlot,z+2);z0++)
							for(l=0;l<QB->network_num[x0][y0][z0];l++)
							{
								j=QB->network_ID[x0][y0][z0][l];
								if(j!=i)
								{
									dis=QB_squaredlenths(QB->atom[i].x-QB->atom[j].x,
														 QB->atom[i].y-QB->atom[j].y,
														 QB->atom[i].z-QB->atom[j].z);
									///====================================
									/// sort neighbor from small to big
									///====================================
									if(QB->neb.num<QB->MCN)
									{
										if(dis<cutoff*cutoff)
										{
											QB->neb.num++;
											
											QB->neb.id[QB->neb.num-1]=j;
											QB->neb.dis[QB->neb.num-1]=dis;

											for(m=QB->neb.num-1;m>0;m--)
												if(QB->neb.dis[m]<QB->neb.dis[m-1])
												{
													template_double=QB->neb.dis[m];
													QB->neb.dis[m]=QB->neb.dis[m-1];
													QB->neb.dis[m-1]=template_double;
													template_int=QB->neb.id[m];
													QB->neb.id[m]=QB->neb.id[m-1];
													QB->neb.id[m-1]=template_int;
												}
										}
									}
									else if(dis<QB->neb.dis[QB->neb.num-1])
									{
										QB->neb.id[QB->neb.num-1]=j;
										QB->neb.dis[QB->neb.num-1]=dis;
										for(m=QB->neb.num-1;m>0;m--)
											if(QB->neb.dis[m]<QB->neb.dis[m-1])
											{
												template_double=QB->neb.dis[m];
												QB->neb.dis[m]=QB->neb.dis[m-1];
												QB->neb.dis[m-1]=template_double;
												template_int=QB->neb.id[m];
												QB->neb.id[m]=QB->neb.id[m-1];
												QB->neb.id[m-1]=template_int;
											}
									}
								}
							}
			}
    for(k=0;k<QB->neb.num;k++)
       	QB->neb.dis[k]=sqrt(QB->neb.dis[k]);
	QB->atom[i].x=backup_x;
	QB->atom[i].y=backup_y;
	QB->atom[i].z=backup_z;
	QB->list_flag=QB->MCN;
}


void QB_list_clear(QB_tools* QB)
{ 
    if(QB->list_flag!=0)
    {
    	free(QB->neb.id);
    	free(QB->neb.dis);
    }
    QB->list_flag=0;
}





