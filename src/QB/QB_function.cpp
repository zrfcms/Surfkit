/*===============================================================
 * This is a part of SPaMD code
 * This program is modified at BUAA by  Z. R. Liu, March. 11 2019
 * Copyright[c] 2017-2019, zrfbuaa group
*================================================================*/
#include"QB.h"
///maths functions
int QB_read_file(QB_tools *QB,const char name[50])
{
    int rt;
    QB_free_atom(QB);
	rt=QB_read_lmp(QB,name);
    if(rt==1)return QB_S_LMP;
    else if(rt==-2)return QB_F_LMPC;

    QB_free_atom(QB);
    rt=QB_read_lmc(QB,name);
    if(rt==1)return QB_S_LMC;
    else if(rt==-2)return QB_F_LMCC;
    
    QB_free_exchannel(QB);
    QB_free_atom(QB);
	
	QB_free_atom(QB);
	rt=QB_read_vasp(QB,name);
    if(rt==1)return QB_S_VASP;
	else if(rt==-2)return QB_F_VASPC;
	
    return QB_F_UNKNOWN;
}
void QB_move_atom(QB_tools *QB,int i,double movement_x,double movement_y,double movement_z)
{
    QB->atom[i].x+=movement_x;
    QB->atom[i].y+=movement_y;
    QB->atom[i].z+=movement_z;
}

void QB_sort_atom(QB_tools *QB,int deal(QB_tools *QB2,int a,int b))
{
    int i,j;
    for(i=1;i<QB->TotalNumber;i++)
        for(j=0;j<QB->TotalNumber-i;j++)
        {
            if(deal(QB,j,j+1))
            	QB_swap(QB,j,j+1);
        }
}

void QB_combine_system(QB_tools *output,QB_tools *input)
{
	int i,j,n;
	int last_TotalNumber=output->TotalNumber;
	if(last_TotalNumber==0)
		output->atom=(atomdata*)malloc((output->TotalNumber+input->TotalNumber)*sizeof(atomdata));
	else
		output->atom=(atomdata*)realloc(output->atom,(output->TotalNumber+input->TotalNumber)*sizeof(atomdata));
	char name[1024];
	for(i=0;i<input->TotalNumber;i++)
    {
        output->atom[last_TotalNumber+i].id=last_TotalNumber+input->atom[i].id;
		QB_load_elename(input,input->atom[i].type-1,name);
        output->atom[last_TotalNumber+i].type=QB_get_elename(output,name);
        output->atom[last_TotalNumber+i].x=input->atom[i].x;
        output->atom[last_TotalNumber+i].y=input->atom[i].y;
        output->atom[last_TotalNumber+i].z=input->atom[i].z;    
    }
	for(i=5;i<QB_slot_sum(output);i++)
	{
		QB_slot_name(output,i,name);
		n=QB_slot_get(input,name);
		if(n>0)
		{
			for(j=0;j<input->TotalNumber;j++)
				QB_slot_save(output,i,last_TotalNumber+j,QB_slot_data(input,n,j));
		}
		else
		{
			for(j=0;j<input->TotalNumber;j++)
				QB_slot_save(output,i,last_TotalNumber+j,0);
		}
	}
    output->TotalNumber+=input->TotalNumber;
}

void QB_copy_system(QB_tools *output,QB_tools *input)
{
	char name[32];
	int i,j;
    QB_free_atom(output);
    QB_free_extra(output);
    QB_free_exchannel(output);
    output->TotalNumber=input->TotalNumber;
    output->TypeNumber=input->TypeNumber;
    output->startx=input->startx;
    output->starty=input->starty;
    output->startz=input->startz;
    output->boundx=input->boundx;
    output->boundy=input->boundy;
    output->boundz=input->boundz;
	output->zerox=input->zerox;
    output->zeroy=input->zeroy;
    output->zeroz=input->zeroz;
	output->px=input->px;
	output->py=input->py;
	output->pz=input->pz;
	output->MCN=input->MCN;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			output->mat[i][j]=input->mat[i][j];
	output->box=input->box;
	
	for(i=0;i<input->TypeNumber;i++)
	{
		 QB_load_elename(input,i,name);
		 QB_save_elename(output,i,name);
	}
	
	if(output->TotalNumber<=0)return;
    
	output->atom=(atomdata*)malloc(output->TotalNumber*sizeof(atomdata));
    for(i=0;i<output->TotalNumber;i++)
    {
        output->atom[i].id=input->atom[i].id;
        output->atom[i].type=input->atom[i].type;
        output->atom[i].x=input->atom[i].x;
        output->atom[i].y=input->atom[i].y;
        output->atom[i].z=input->atom[i].z;    
    }
	
    QB_free_exchannel(output);
    output->exchannel=input->exchannel;
    if(output->exchannel)
    {
        output->exch=(char**)malloc(output->exchannel*sizeof(char*));
        output->exdat=(double**)malloc(output->exchannel*sizeof(double*));
        for(i=0;i<output->exchannel;i++)
        {
            output->exch[i]=(char*)malloc(50*sizeof(char));
            strcpy(output->exch[i],input->exch[i]);
            output->exdat[i]=(double*)malloc(output->TotalNumber*sizeof(double));
            for(j=0;j<output->TotalNumber;j++)
                output->exdat[i][j]=input->exdat[i][j];
        }
    }
	
	QB_free_extra(output);
    output->num_exint=input->num_exint;
    if(output->num_exint)
    {
		output->name_exint=(char**)malloc(output->num_exint*sizeof(char*));
		output->par_exint=(int**)malloc(output->num_exint*sizeof(int*));
        for(i=0;i<output->num_exint;i++)
        {
	        output->name_exint[i]=(char*)malloc(50*sizeof(char));
            strcpy(output->name_exint[i],input->name_exint[i]);
	        output->par_exint[i]=(int*)malloc(output->TotalNumber*sizeof(int));
	        for(j=0;j<output->TotalNumber;j++)
                output->par_exint[i][j]=input->par_exint[i][j];
        }
    }
	
    output->num_exdouble=input->num_exdouble;
    if(output->num_exdouble)
    {
		output->name_exdouble=(char**)malloc(output->num_exdouble*sizeof(char*));
		output->par_exdouble=(double**)malloc(output->num_exdouble*sizeof(double*));
        for(i=0;i<output->num_exdouble;i++)
        {
	        output->name_exdouble[i]=(char*)malloc(50*sizeof(char));
            strcpy(output->name_exdouble[i],input->name_exdouble[i]);
	        output->par_exdouble[i]=(double*)malloc(output->TotalNumber*sizeof(double));
	        for(j=0;j<output->TotalNumber;j++)
                output->par_exdouble[i][j]=input->par_exdouble[i][j];
        }
    }
    
}

void QB_copy_box(QB_tools *output,QB_tools *input)
{
    output->startx=input->startx;
    output->starty=input->starty;
    output->startz=input->startz;
    output->boundx=input->boundx;
    output->boundy=input->boundy;
    output->boundz=input->boundz;
	output->zerox=input->zerox;
    output->zeroy=input->zeroy;
    output->zeroz=input->zeroz;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			output->mat[i][j]=input->mat[i][j];
	output->box=input->box; 
}


void QB_move_box(QB_tools *QB,double movement_x,double movement_y,double movement_z)
{
    int i;
    for(i=0;i<QB->TotalNumber;i++)
    {
        QB->atom[i].x+=movement_x;
        QB->atom[i].y+=movement_y;
        QB->atom[i].z+=movement_z;
    }
    QB->startx+=movement_x;
    QB->starty+=movement_y;
    QB->startz+=movement_z;
	QB->zerox+=movement_x;
    QB->zeroy+=movement_y;
    QB->zeroz+=movement_z;
}

void QB_duplicate(QB_tools *QB,int nx,int ny,int nz)
{
    int i,j,k,m,n;
	int max_id=QB->TotalNumber;
	for(m=0;m<QB->TotalNumber;m++)      
	if(QB->atom[m].id>max_id)
		max_id=QB->atom[m].id;
    QB->atom=(atomdata*)realloc(QB->atom,nx*ny*nz*QB->TotalNumber*sizeof(atomdata));
    //copy model according to input parameter
    n=0;
	if(QB->box)
	{
		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
			   for(k=0;k<nz;k++)
				   for(m=0;m<QB->TotalNumber;m++)      
				   {
						QB->atom[n].x=QB->atom[m].x+(double)i*QB->mat[0][0]+(double)j*QB->mat[1][0]+(double)k*QB->mat[2][0];
						QB->atom[n].y=QB->atom[m].y+(double)i*QB->mat[0][1]+(double)j*QB->mat[1][1]+(double)k*QB->mat[2][1];
						QB->atom[n].z=QB->atom[m].z+(double)i*QB->mat[0][2]+(double)j*QB->mat[1][2]+(double)k*QB->mat[2][2];
						QB->atom[n].type=QB->atom[m].type;
						if(QB->atom[m].id<0)QB->atom[n].id=-(-QB->atom[m].id+(k+nz*j+nz*ny*i)*max_id);
						else QB->atom[n].id=QB->atom[m].id+(k+nz*j+nz*ny*i)*max_id;
						n++;
				   }
	}
	else
	{
		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
			   for(k=0;k<nz;k++)
				   for(m=0;m<QB->TotalNumber;m++)      
				   {
						QB->atom[n].x=QB->atom[m].x+(double)i*QB->boundx;
						QB->atom[n].y=QB->atom[m].y+(double)j*QB->boundy;
						QB->atom[n].z=QB->atom[m].z+(double)k*QB->boundz;
						QB->atom[n].type=QB->atom[m].type;
						if(QB->atom[m].id<0)QB->atom[n].id=-(-QB->atom[m].id+(k+nz*j+nz*ny*i)*max_id);
						else QB->atom[n].id=QB->atom[m].id+(k+nz*j+nz*ny*i)*max_id;
						n++;
				   }
	}
    QB->boundx*=(double)nx;
    QB->boundy*=(double)ny;
    QB->boundz*=(double)nz;
	if(QB->box)
	{
		QB->startx*=(double)nx;
		QB->starty*=(double)ny;
		QB->startz*=(double)nz;
	}
	for(int i=0;i<3;i++)
	{
		QB->mat[0][i]*=nx;
		QB->mat[1][i]*=ny;
		QB->mat[2][i]*=nz;
	}
	int n0=QB->TotalNumber;
	int l;
    QB->TotalNumber*=(nx*ny*nz);
	QB_update_extra(QB);
	n=0;
	for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
			   for(k=0;k<nz;k++)
				   for(m=0;m<n0;m++)      
				   {
					    for(l=5;l<QB_slot_sum(QB);l++)
						{
							QB_slot_save(QB,l,n,QB_slot_data(QB,l,m));
						}
						n++;
				   }
	
}

void QB_rotate(QB_tools *QB,int i,QB_vector center,QB_vector vv, float angle)
{
    angle*=M_PI/180.0;
    double normalize;
    QB->atom[i].x-=center.x;
    QB->atom[i].y-=center.y;
    QB->atom[i].z-=center.z;
    normalize=sqrt(QB_squaredlenths(vv.x,vv.y,vv.z));
    double u = vv.x/normalize;
    double v = vv.y/normalize;
    double w = vv.z/normalize;
    double x = QB->atom[i].x;
    double y = QB->atom[i].y;
    double z = QB->atom[i].z;
    QB->atom[i].x=x*(u*u+(v*v+w*w)*cos(angle))+y*(u*v*(1-cos(angle))+w*sin(angle))+z*(u*w*(1-cos(angle))-v*sin(angle));
    QB->atom[i].y=x*(u*v*(1-cos(angle))-w*sin(angle))+y*(v*v+(u*u+w*w)*cos(angle))+z*(v*w*(1-cos(angle))+u*sin(angle));
    QB->atom[i].z=x*(u*w*(1-cos(angle))+v*sin(angle))+y*(v*w*(1-cos(angle))-u*sin(angle))+z*(w*w+(u*u+v*v)*cos(angle));
    QB->atom[i].x+=center.x;
    QB->atom[i].y+=center.y;
    QB->atom[i].z+=center.z;
}

void QB_mirror(QB_tools *QB,int i,QB_vector p0, double dis)
{
    double norm=sqrt(QB_squaredlenths(p0.x,p0.y,p0.z));
    QB_vector p;
    p.x=p0.x/norm;
    p.y=p0.y/norm;
    p.z=p0.z/norm;
    dis/=norm;
    double dot=QB->atom[i].x*p.x+QB->atom[i].y*p.y+QB->atom[i].z*p.z-dis;
    QB_create_atom(QB,QB->atom[i].type,QB->atom[i].x-2*p.x*dot,
                QB->atom[i].y-2*p.y*dot,
                QB->atom[i].z-2*p.z*dot);
}

void QB_swichxz(QB_tools *QB)
{
    double temp;
	int tempi;
    int i;
    for(i=0;i<QB->TotalNumber;i++)
    {
        temp=QB->atom[i].x;
        QB->atom[i].x=QB->atom[i].z;
        QB->atom[i].z=temp;
    }
    temp=QB->boundz;
    QB->boundz=QB->boundx;
    QB->boundx=temp;
    temp=QB->startz;
    QB->startz=QB->startx;
    QB->startx=temp;
	
	tempi=QB->px;
	QB->px=QB->pz;
	QB->pz=tempi;
	if(QB->box)
	{
		temp=QB->mat[0][0];
		QB->mat[0][0]=QB->mat[0][2];
		QB->mat[0][2]=temp;
		temp=QB->mat[1][0];
		QB->mat[1][0]=QB->mat[1][2];
		QB->mat[1][2]=temp;
		temp=QB->mat[2][0];
		QB->mat[2][0]=QB->mat[2][2];
		QB->mat[2][2]=temp;
		temp=QB->zerox;
		QB->zerox=QB->zeroz;
		QB->zeroz=temp;
	}
}

void QB_swichxy(QB_tools *QB)
{
    double temp;
	int tempi;
    int i;
    for(i=0;i<QB->TotalNumber;i++)
    {
        temp=QB->atom[i].x;
        QB->atom[i].x=QB->atom[i].y;
        QB->atom[i].y=temp;
    }
    temp=QB->boundy;
    QB->boundy=QB->boundx;
    QB->boundx=temp;
    temp=QB->starty;
    QB->starty=QB->startx;
    QB->startx=temp;
	
	tempi=QB->px;
	QB->px=QB->py;
	QB->py=tempi;
	if(QB->box)
	{
		temp=QB->mat[0][0];
		QB->mat[0][0]=QB->mat[0][1];
		QB->mat[0][1]=temp;
		temp=QB->mat[1][0];
		QB->mat[1][0]=QB->mat[1][1];
		QB->mat[1][1]=temp;
		temp=QB->mat[2][0];
		QB->mat[2][0]=QB->mat[2][1];
		QB->mat[2][1]=temp;
		temp=QB->zerox;
		QB->zerox=QB->zeroy;
		QB->zeroy=temp;
	}
	
}

void QB_swichyz(QB_tools *QB)
{
    double temp;
	int tempi;
    int i;
    for(i=0;i<QB->TotalNumber;i++)
    {
        temp=QB->atom[i].y;
        QB->atom[i].y=QB->atom[i].z;
        QB->atom[i].z=temp;
    }
    temp=QB->boundz;
    QB->boundz=QB->boundy;
    QB->boundy=temp;
    temp=QB->startz;
    QB->startz=QB->starty;
    QB->starty=temp;
	
	tempi=QB->py;
	QB->py=QB->pz;
	QB->pz=tempi;
	if(QB->box)
	{
		temp=QB->mat[0][1];
		QB->mat[0][1]=QB->mat[0][2];
		QB->mat[0][2]=temp;
		temp=QB->mat[1][1];
		QB->mat[1][1]=QB->mat[1][2];
		QB->mat[1][2]=temp;
		temp=QB->mat[2][0];
		QB->mat[2][1]=QB->mat[2][2];
		QB->mat[2][2]=temp;
		temp=QB->zeroy;
		QB->zeroy=QB->zeroz;
		QB->zeroz=temp;
	}
}

void QB_wrap(QB_tools *QB,int i)
{
	if(QB->box)
	{
		QB_vector 	face[3];
		double 	endface[9];
		double x,y,z;
		int point=0,lp=0,n;
        face[0].x=QB->mat[1][1]*QB->mat[2][2]-QB->mat[1][2]*QB->mat[2][1];
        face[0].y=QB->mat[1][2]*QB->mat[2][0]-QB->mat[1][0]*QB->mat[2][2];
        face[0].z=QB->mat[1][0]*QB->mat[2][1]-QB->mat[1][1]*QB->mat[2][0];
		if(face[0].x*QB->mat[0][0]+face[0].y*QB->mat[0][1]+face[0].z*QB->mat[0][2]<0){face[0].x*=-1;face[0].y*=-1;face[0].z*=-1;}
        face[1].x=QB->mat[2][1]*QB->mat[0][2]-QB->mat[2][2]*QB->mat[0][1];
        face[1].y=QB->mat[2][2]*QB->mat[0][0]-QB->mat[2][0]*QB->mat[0][2];
        face[1].z=QB->mat[2][0]*QB->mat[0][1]-QB->mat[2][1]*QB->mat[0][0];
		if(face[1].x*QB->mat[1][0]+face[1].y*QB->mat[1][1]+face[1].z*QB->mat[1][2]<0){face[1].x*=-1;face[1].y*=-1;face[1].z*=-1;}
		face[2].x=QB->mat[0][1]*QB->mat[1][2]-QB->mat[0][2]*QB->mat[1][1];
        face[2].y=QB->mat[0][2]*QB->mat[1][0]-QB->mat[0][0]*QB->mat[1][2];
        face[2].z=QB->mat[0][0]*QB->mat[1][1]-QB->mat[0][1]*QB->mat[1][0];
		if(face[2].x*QB->mat[2][0]+face[2].y*QB->mat[2][1]+face[2].z*QB->mat[2][2]<0){face[2].x*=-1;face[2].y*=-1;face[2].z*=-1;}
    //init boundary face
    //end face 0~5 means D in Ax+By+Cz=D 6->8 means cutoff*sqrt(A*A+B*B+C*C)
        for(n=0;n<3;n++)
        {
            endface[n]=face[n].x*QB->zerox+face[n].y*QB->zeroy+face[n].z*QB->zeroz;
            endface[n+3]=face[n].x*(QB->mat[n][0]+QB->zerox)+
						 face[n].y*(QB->mat[n][1]+QB->zeroy)+
						 face[n].z*(QB->mat[n][2]+QB->zeroz);
        }
    //first we need to to check how many atom is too near to boundary
    //per run will check one dirction
        x=QB->atom[i].x*face[0].x+QB->atom[i].y*face[0].y+QB->atom[i].z*face[0].z;
        y=QB->atom[i].x*face[1].x+QB->atom[i].y*face[1].y+QB->atom[i].z*face[1].z;
        z=QB->atom[i].x*face[2].x+QB->atom[i].y*face[2].y+QB->atom[i].z*face[2].z;
           while(point!=6&&lp<10)
		{
			point=0;
			x=QB->atom[i].x;
			y=QB->atom[i].y;
			z=QB->atom[i].z;
			if(x*face[0].x+y*face[0].y+z*face[0].z<endface[0]){QB->atom[i].x+=QB->mat[0][0];QB->atom[i].y+=QB->mat[0][1];QB->atom[i].z+=QB->mat[0][2];}else point++;
			if(x*face[0].x+y*face[0].y+z*face[0].z>=endface[3]){QB->atom[i].x-=QB->mat[0][0];QB->atom[i].y-=QB->mat[0][1];QB->atom[i].z-=QB->mat[0][2];}else point++;
			if(x*face[1].x+y*face[1].y+z*face[1].z<endface[1]){QB->atom[i].x+=QB->mat[1][0];QB->atom[i].y+=QB->mat[1][1];QB->atom[i].z+=QB->mat[1][2];}else point++;
			if(x*face[1].x+y*face[1].y+z*face[1].z>=endface[4]){QB->atom[i].x-=QB->mat[1][0];QB->atom[i].y-=QB->mat[1][1];QB->atom[i].z-=QB->mat[1][2];}else point++;
			if(x*face[2].x+y*face[2].y+z*face[2].z<endface[2]){QB->atom[i].x+=QB->mat[2][0];QB->atom[i].y+=QB->mat[2][1];QB->atom[i].z+=QB->mat[2][2];}else point++;
			if(x*face[2].x+y*face[2].y+z*face[2].z>=endface[5]){QB->atom[i].x-=QB->mat[2][0];QB->atom[i].y-=QB->mat[2][1];QB->atom[i].z-=QB->mat[2][2];}else point++;
			lp++;
		}         
	}
	else	
	{
		while(QB->atom[i].x>(QB->boundx+QB->startx))
			QB->atom[i].x-=QB->boundx;
		while(QB->atom[i].y>(QB->boundy+QB->starty))
			QB->atom[i].y-=QB->boundy;
		while(QB->atom[i].z>(QB->boundz+QB->startz))
			QB->atom[i].z-=QB->boundz;
		while(QB->atom[i].x<QB->startx)
			QB->atom[i].x+=QB->boundx;
		while(QB->atom[i].y<QB->starty)
			QB->atom[i].y+=QB->boundy;
		while(QB->atom[i].z<QB->startz)
			QB->atom[i].z+=QB->boundz;
	}
}

void QB_wrap_vec(QB_tools *QB,double *v_x,double *v_y,double *v_z)//wrap any point back to box
{
	if(QB->box)
	{
		QB_vector 	face[3];
		double 	endface[9];
		double x,y,z;
		int point=0,lp=0,n;
        face[0].x=QB->mat[1][1]*QB->mat[2][2]-QB->mat[1][2]*QB->mat[2][1];
        face[0].y=QB->mat[1][2]*QB->mat[2][0]-QB->mat[1][0]*QB->mat[2][2];
        face[0].z=QB->mat[1][0]*QB->mat[2][1]-QB->mat[1][1]*QB->mat[2][0];
		if(face[0].x*QB->mat[0][0]+face[0].y*QB->mat[0][1]+face[0].z*QB->mat[0][2]<0){face[0].x*=-1;face[0].y*=-1;face[0].z*=-1;}
        face[1].x=QB->mat[2][1]*QB->mat[0][2]-QB->mat[2][2]*QB->mat[0][1];
        face[1].y=QB->mat[2][2]*QB->mat[0][0]-QB->mat[2][0]*QB->mat[0][2];
        face[1].z=QB->mat[2][0]*QB->mat[0][1]-QB->mat[2][1]*QB->mat[0][0];
		if(face[1].x*QB->mat[1][0]+face[1].y*QB->mat[1][1]+face[1].z*QB->mat[1][2]<0){face[1].x*=-1;face[1].y*=-1;face[1].z*=-1;}
		face[2].x=QB->mat[0][1]*QB->mat[1][2]-QB->mat[0][2]*QB->mat[1][1];
        face[2].y=QB->mat[0][2]*QB->mat[1][0]-QB->mat[0][0]*QB->mat[1][2];
        face[2].z=QB->mat[0][0]*QB->mat[1][1]-QB->mat[0][1]*QB->mat[1][0];
		if(face[2].x*QB->mat[2][0]+face[2].y*QB->mat[2][1]+face[2].z*QB->mat[2][2]<0){face[2].x*=-1;face[2].y*=-1;face[2].z*=-1;}
    //init boundary face
    //end face 0~5 means D in Ax+By+Cz=D 6->8 means cutoff*sqrt(A*A+B*B+C*C)
        for(n=0;n<3;n++)
        {
            endface[n]=face[n].x*QB->zerox+face[n].y*QB->zeroy+face[n].z*QB->zeroz;
            endface[n+3]=face[n].x*(QB->mat[n][0]+QB->zerox)+
						 face[n].y*(QB->mat[n][1]+QB->zeroy)+
						 face[n].z*(QB->mat[n][2]+QB->zeroz);
        }
    //first we need to to check how many atom is too near to boundary
    //per run will check one dirction
        x=(*v_x)*face[0].x+(*v_y)*face[0].y+(*v_z)*face[0].z;
        y=(*v_x)*face[1].x+(*v_y)*face[1].y+(*v_z)*face[1].z;
        z=(*v_x)*face[2].x+(*v_y)*face[2].y+(*v_z)*face[2].z;
           while(point!=6&&lp<10)
		{
			point=0;
			x=(*v_x);
			y=(*v_y);
			z=(*v_z);
			if(x*face[0].x+y*face[0].y+z*face[0].z<endface[0]){(*v_x)+=QB->mat[0][0];(*v_y)+=QB->mat[0][1];(*v_z)+=QB->mat[0][2];}else point++;
			if(x*face[0].x+y*face[0].y+z*face[0].z>=endface[3]){(*v_x)-=QB->mat[0][0];(*v_y)-=QB->mat[0][1];(*v_z)-=QB->mat[0][2];}else point++;
			if(x*face[1].x+y*face[1].y+z*face[1].z<endface[1]){(*v_x)+=QB->mat[1][0];(*v_y)+=QB->mat[1][1];(*v_z)+=QB->mat[1][2];}else point++;
			if(x*face[1].x+y*face[1].y+z*face[1].z>=endface[4]){(*v_x)-=QB->mat[1][0];(*v_y)-=QB->mat[1][1];(*v_z)-=QB->mat[1][2];}else point++;
			if(x*face[2].x+y*face[2].y+z*face[2].z<endface[2]){(*v_x)+=QB->mat[2][0];(*v_y)+=QB->mat[2][1];(*v_z)+=QB->mat[2][2];}else point++;
			if(x*face[2].x+y*face[2].y+z*face[2].z>=endface[5]){(*v_x)-=QB->mat[2][0];(*v_y)-=QB->mat[2][1];(*v_z)-=QB->mat[2][2];}else point++;
			lp++;
		}         
	}
	else	
	{
		while((*v_x)>(QB->boundx+QB->startx))
			(*v_x)-=QB->boundx;
		while((*v_y)>(QB->boundy+QB->starty))
			(*v_y)-=QB->boundy;
		while((*v_z)>(QB->boundz+QB->startz))
			(*v_z)-=QB->boundz;
		while((*v_x)<QB->startx)
			(*v_x)+=QB->boundx;
		while((*v_y)<QB->starty)
			(*v_y)+=QB->boundy;
		while((*v_z)<QB->startz)
			(*v_z)+=QB->boundz;
	}
}

void QB_combine(QB_tools *QB1,QB_tools *QB2)
{
	int i,j;
	if((QB1->TotalNumber)!=(QB2->TotalNumber))
		return; 
	for(i=0;i<QB2->exchannel;i++)
		QB_add_exdouble(QB1,QB2->exch[i]);
	for(i=0;i<QB1->TotalNumber;i++)
		for(j=0;j<QB2->exchannel;j++)
			QB_save_data(QB1,i,QB2->exch[j],QB2->exdat[j][i]);
}

int QB_inbox(QB_tools *QB,int i)
{
	if(QB->box)
	{
		QB_vector 	face[3];
		double 	endface[9];
		double x,y,z;
		int n;
        face[0].x=QB->mat[1][1]*QB->mat[2][2]-QB->mat[1][2]*QB->mat[2][1];
        face[0].y=QB->mat[1][2]*QB->mat[2][0]-QB->mat[1][0]*QB->mat[2][2];
        face[0].z=QB->mat[1][0]*QB->mat[2][1]-QB->mat[1][1]*QB->mat[2][0];
		if(face[0].x*QB->mat[0][0]+face[0].y*QB->mat[0][1]+face[0].z*QB->mat[0][2]<0){face[0].x*=-1;face[0].y*=-1;face[0].z*=-1;}
        face[1].x=QB->mat[2][1]*QB->mat[0][2]-QB->mat[2][2]*QB->mat[0][1];
        face[1].y=QB->mat[2][2]*QB->mat[0][0]-QB->mat[2][0]*QB->mat[0][2];
        face[1].z=QB->mat[2][0]*QB->mat[0][1]-QB->mat[2][1]*QB->mat[0][0];
		if(face[1].x*QB->mat[1][0]+face[1].y*QB->mat[1][1]+face[1].z*QB->mat[1][2]<0){face[1].x*=-1;face[1].y*=-1;face[1].z*=-1;}
		face[2].x=QB->mat[0][1]*QB->mat[1][2]-QB->mat[0][2]*QB->mat[1][1];
        face[2].y=QB->mat[0][2]*QB->mat[1][0]-QB->mat[0][0]*QB->mat[1][2];
        face[2].z=QB->mat[0][0]*QB->mat[1][1]-QB->mat[0][1]*QB->mat[1][0];
		if(face[2].x*QB->mat[2][0]+face[2].y*QB->mat[2][1]+face[2].z*QB->mat[2][2]<0){face[2].x*=-1;face[2].y*=-1;face[2].z*=-1;}
    //init boundary face
    //end face 0~5 means D in Ax+By+Cz=D 6->8 means cutoff*sqrt(A*A+B*B+C*C)
        for(n=0;n<3;n++)
        {
            endface[n]=face[n].x*QB->zerox+face[n].y*QB->zeroy+face[n].z*QB->zeroz;
            endface[n+3]=face[n].x*(QB->mat[0][0]+QB->mat[1][0]+QB->mat[2][0]+QB->zerox)+
						 face[n].y*(QB->mat[0][1]+QB->mat[1][1]+QB->mat[2][1]+QB->zeroy)+
						 face[n].z*(QB->mat[0][2]+QB->mat[1][2]+QB->mat[2][2]+QB->zeroz);
        }
		x=QB->atom[i].x;
		y=QB->atom[i].y;
		z=QB->atom[i].z;
		if(x*face[0].x+y*face[0].y+z*face[0].z<endface[0])return 0;
		if(x*face[0].x+y*face[0].y+z*face[0].z>=endface[3])return 0;
		if(x*face[1].x+y*face[1].y+z*face[1].z<endface[1])return 0;
		if(x*face[1].x+y*face[1].y+z*face[1].z>=endface[4])return 0;
		if(x*face[2].x+y*face[2].y+z*face[2].z<endface[2])return 0;
		if(x*face[2].x+y*face[2].y+z*face[2].z>=endface[5])return 0; 
	}
	else 
	{
		if(QB->atom[i].x>(QB->boundx+QB->startx))
			return 0;
		if(QB->atom[i].y>(QB->boundy+QB->starty))
			return 0;
		if(QB->atom[i].z>(QB->boundz+QB->startz))
			return 0;
		if(QB->atom[i].x<QB->startx)
			return 0;
		if(QB->atom[i].y<QB->starty)
			return 0;
		if(QB->atom[i].z<QB->startz)
			return 0;
	}
    return 1;
}

void QB_bound_update(QB_tools *QB)
{
	QB->box=1;
	QB->boundx=fabs(QB->mat[0][0])+fabs(QB->mat[1][0])+fabs(QB->mat[2][0]);
    QB->boundy=fabs(QB->mat[0][1])+fabs(QB->mat[1][1])+fabs(QB->mat[2][1]);
	QB->boundz=fabs(QB->mat[0][2])+fabs(QB->mat[1][2])+fabs(QB->mat[2][2]);
	QB->startx=QB->zerox=QB->starty=0;
	if(QB->mat[0][0]<0)QB->startx+=QB->mat[0][0];
	if(QB->mat[1][0]<0)QB->startx+=QB->mat[1][0];
	if(QB->mat[2][0]<0)QB->startx+=QB->mat[2][0];
	if(QB->mat[0][1]<0)QB->starty+=QB->mat[0][1];
	if(QB->mat[1][1]<0)QB->starty+=QB->mat[1][1];
	if(QB->mat[2][1]<0)QB->starty+=QB->mat[2][1];
	if(QB->mat[0][2]<0)QB->startz+=QB->mat[0][2];
	if(QB->mat[1][2]<0)QB->startz+=QB->mat[1][2];
	if(QB->mat[2][2]<0)QB->startz+=QB->mat[2][2];
}


void QB_mat_update(QB_tools *QB)
{
	if(QB->box==0)
	{
		QB->mat[0][0]=QB->boundx;
		QB->mat[1][1]=QB->boundy;
		QB->mat[2][2]=QB->boundz;
		QB->mat[1][0]=QB->mat[2][0]=QB->mat[0][1]=QB->mat[2][1]=QB->mat[0][2]=QB->mat[1][2]=0;
	}
}
//==========================================================================

double QB_squaredlenths(double x,double y,double z)
{
    return x*x+y*y+z*z;
}

//==========================================================================

int QB_Max(int x,int y)
{
    if(x>y)
        return x;
    return y;
}

//==========================================================================

int QB_Min(int x,int y)
{
    if(x<y)
        return x;
    return y;
}


//==========================================================================
int QB_checkdat_err=0;
///get double data from strings
double QB_checkdat(char c[50])
{
    double ans;
    char *err;
    ans=strtod(c,&err);
    if(strlen(err)!=0)
    {
        printf(">>ERROR: error data type %s should be float or int type",err);
        QB_checkdat_err=1;
    }
    return ans;
}

//==========================================================================

///check whether a string is not a "number"
int QB_checkchar(char c[50])
{
    char *err;
    strtod(c,&err);
    if(strlen(err)!=0)
        return 1;
    return 0;
}

//==========================================================================

