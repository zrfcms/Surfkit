/*===============================================================
 * This is a part of SPaMD code
 * This program is modified at BUAA by  Z. R. Liu, March. 11 2019
 * Copyright[c] 2017-2019, zrfbuaa group
*================================================================*/
#include"QB.h"
int QB_pbc_start=0;
int QB_pbc_end=0;
void QB_pbc(QB_tools *QB,double rmax)
{
    if(QB_pbc_start!=0)
	if(QB_pbc_end!=0)
	{
	    QB_printf("PBC init twice!\n");
	    return;
	}
    double  u[3]={0,1,-1};
    double  x,y,z;
    int     k,p,q;
    int     judge;
    int     n;
    int     now;
    int     counter=0;
	
	if(QB->box==1)
    {
		QB_vector face[3];
		double 	endface[9];
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
            endface[n+6]=sqrt(QB_squaredlenths(face[n].x,face[n].y,face[n].z))*rmax;
        }
    //first we need to to check how many atom is too near to boundary
    //per run will check one dirction
        for(now=0;now<QB->TotalNumber;now++)
        {
            x=QB->atom[now].x*face[0].x+QB->atom[now].y*face[0].y+QB->atom[now].z*face[0].z;
            y=QB->atom[now].x*face[1].x+QB->atom[now].y*face[1].y+QB->atom[now].z*face[1].z;
            z=QB->atom[now].x*face[2].x+QB->atom[now].y*face[2].y+QB->atom[now].z*face[2].z;
            for(k=0;k<=2*QB->px;k++)
                for(p=0;p<=2*QB->py;p++)
                    for(q=0;q<=2*QB->pz;q++)//px = predioc boundary condition in dirction x
                    {
                        judge=0;
                        if(k==2&&fabs(x-endface[3])<endface[6])
                            judge+=1;
                        if(p==2&&fabs(y-endface[4])<endface[7])
                            judge+=1;
                        if(q==2&&fabs(z-endface[5])<endface[8])
                            judge+=1;
                        if(k==1&&fabs(x-endface[0])<endface[6])
                            judge+=1;
                        if(p==1&&fabs(y-endface[1])<endface[7])
                            judge+=1;
                        if(q==1&&fabs(z-endface[2])<endface[8])
                            judge+=1;
                        if(k==0)
                            judge+=1;
                        if(p==0)
                            judge+=1;
                        if(q==0)
                            judge+=1;
                        if(q==0&&k==0&&p==0)
                            judge-=1;
                        if(judge==3)
                            counter++;
                    }
            }
        counter+=QB->TotalNumber;
        QB->atom=(atomdata*)realloc(QB->atom,counter*sizeof(atomdata));
        //as we doesn't know how many atom needs copied we checked before and we saved it in counter
        //now beigin copy atom data
        now=QB->TotalNumber;
        for(n=0;n<QB->TotalNumber;n++)
        {
            x=QB->atom[n].x*face[0].x+QB->atom[n].y*face[0].y+QB->atom[n].z*face[0].z;
            y=QB->atom[n].x*face[1].x+QB->atom[n].y*face[1].y+QB->atom[n].z*face[1].z;
            z=QB->atom[n].x*face[2].x+QB->atom[n].y*face[2].y+QB->atom[n].z*face[2].z;
            for(k=0;k<=2*QB->px;k++)
                for(p=0;p<=2*QB->py;p++)
                    for(q=0;q<=2*QB->pz;q++)
                    {
                        judge=0;
                        if(k==2&&fabs(x-endface[3])<endface[6])
                            judge+=1;
                        if(p==2&&fabs(y-endface[4])<endface[7])
                            judge+=1;
                        if(q==2&&fabs(z-endface[5])<endface[8])
                            judge+=1;
                        if(k==1&&fabs(x-endface[0])<endface[6])
                            judge+=1;
                        if(p==1&&fabs(y-endface[1])<endface[7])
                            judge+=1;
                        if(q==1&&fabs(z-endface[2])<endface[8])
                            judge+=1;
                        if(k==0)
                            judge+=1;
                        if(p==0)
                            judge+=1;
                        if(q==0)
                            judge+=1;
                        if(q==0&&k==0&&p==0)
                            judge-=1;
                        if(judge==3)
                        {
                            QB->atom[now].x=QB->atom[n].x+u[k]*QB->mat[0][0]+u[p]*QB->mat[1][0]+u[q]*QB->mat[2][0];
                            QB->atom[now].y=QB->atom[n].y+u[k]*QB->mat[0][1]+u[p]*QB->mat[1][1]+u[q]*QB->mat[2][1];
                            QB->atom[now].z=QB->atom[n].z+u[k]*QB->mat[0][2]+u[p]*QB->mat[1][2]+u[q]*QB->mat[2][2];
                            QB->atom[now].type=QB->atom[n].type;
							QB->atom[now].id=QB->atom[n].id;
                            now++;
                        }
                }
        }
    }
	else if(QB->box==0)//square box -> we can speed up copy computing by running a different methord
    {
		for(now=0;now<QB->TotalNumber;now++)//check out how many atoms we need to add to system
			for(k=0;k<=2*QB->px;k++)
				for(p=0;p<=2*QB->py;p++)
					for(q=0;q<=2*QB->pz;q++)
					{
						x=QB->atom[now].x;
						y=QB->atom[now].y;
						z=QB->atom[now].z;
						judge=0;
						if(k==2&&x>QB->startx+QB->boundx-rmax)
							judge+=1;
						if(p==2&&y>QB->starty+QB->boundy-rmax)
							judge+=1;
						if(q==2&&z>QB->startz+QB->boundz-rmax)
							judge+=1;
						if(k==1&&x<QB->startx+rmax)
							judge+=1;
						if(p==1&&y<QB->starty+rmax)
							judge+=1;
						if(q==1&&z<QB->startz+rmax)
							judge+=1;
						if(k==0)
							judge+=1;
						if(p==0)
							judge+=1;
						if(q==0)
							judge+=1;
						if(q==0&&k==0&&p==0)
							judge-=1;
						if(judge==3)
							counter++;
					}
			counter+=QB->TotalNumber;
			QB->atom=(atomdata *)realloc(QB->atom,counter*sizeof(atomdata));
			//build new atoms which are too near to boundaries
			for(n=0;n<QB->TotalNumber;n++)
				for(k=0;k<=2*QB->px;k++)//consider setting in LAMMPS normal boundary condition will not be copied
					for(p=0;p<=2*QB->py;p++)
						for(q=0;q<=2*QB->pz;q++)
						{
							x=QB->atom[n].x;
							y=QB->atom[n].y;
							z=QB->atom[n].z;
							judge=0;
							if(k==2&&x>QB->startx+QB->boundx-rmax)
							{
								x-=QB->boundx;
								judge+=1;
							}
							if(p==2&&y>QB->starty+QB->boundy-rmax)
							{
								y-=QB->boundy;
								judge+=1;
							}
							if(q==2&&z>QB->startz+QB->boundz-rmax)
							{
								z-=QB->boundz;
								judge+=1;
							}
							if(k==1&&x<QB->startx+rmax)
							{
								x+=QB->boundx;
								judge+=1;
							}
							if(p==1&&y<QB->starty+rmax)
							{
								y+=QB->boundy;
								judge+=1;
							}
							if(q==1&&z<QB->startz+rmax)
							{
								z+=QB->boundz;
								judge+=1;
							}
							if(k==0)
								judge+=1;
							if(p==0)
								judge+=1;
							if(q==0)
								judge+=1;
							if(q==0&&k==0&&p==0)
								judge-=1;
						   if(judge==3)
							{
								QB->atom[now].x=x;
								QB->atom[now].y=y;
								QB->atom[now].z=z;
								QB->atom[now].type=QB->atom[n].type;
								QB->atom[now].id=QB->atom[n].id;
								now++;
							}
						}
	}
    QB_pbc_start=QB->TotalNumber;
    QB_pbc_end=counter;
    QB->TotalNumber=counter;
}

void QB_pbc_clean(QB_tools *QB)
{
    int i;
    for(i=QB_pbc_start;i<QB_pbc_end;i++)
    {
        QB_delete_atom(QB,i);
    }
    QB_pbc_start=QB_pbc_end=0;
}

