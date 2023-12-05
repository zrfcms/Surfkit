#include"./../QB/QB.h"
#include"QSPG.h"

extern void SPG_add_symmetry_atom(SPG_tools *spg,int hall_number,int type,double x,double y,double z,double symprec);
void QSPG_getbox(int mode,int R_flag,double lattice[3][3],double ia,double ib,double ic,double ialpha,double ibeta,double igamma)
{ 
	double metric[3][3];
	metric[0][0]=ia*ia;
	metric[0][1]=ia*ib*cos(igamma/180*M_PI);
	metric[0][2]=ia*ic*cos(ibeta/180*M_PI);
	metric[1][0]=ia*ib*cos(igamma/180*M_PI);
	metric[1][1]=ib*ib;
	metric[1][2]=ib*ic*cos(ialpha/180*M_PI);
	metric[2][0]=ia*ic*cos(ibeta/180*M_PI);
	metric[2][1]=ib*ic*cos(ialpha/180*M_PI);
	metric[2][2]=ic*ic;

	if(mode==1)//set_tricli(lattice, metric);
	{ 
		double a, b, c, alpha, beta, gamma, cg, cb, ca, sg;

		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		alpha = acos(metric[1][2] / b / c);
		beta = acos(metric[0][2] / a / c);
		gamma = acos(metric[0][1] / a / b);
		cg = cos(gamma);
		cb = cos(beta);
		ca = cos(alpha);
		sg = sin(gamma);

		lattice[0][0] = a;
		lattice[0][1] = b * cg;
		lattice[0][2] = c * cb;
		lattice[1][1] = b * sg;
		lattice[1][2] = c * (ca - cb * cg) / sg;
		lattice[2][2] = c * sqrt(1 - ca * ca - cb * cb - cg * cg +
							   2 * ca * cb * cg) / sg;
	}

	if(mode==2)//set_monocli(lattice, metric);
	{ 
		/* Lattice is expected to be C centring */
		double a, b, c, beta;

		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = a;
		lattice[1][1] = b;
		beta = acos(metric[0][2] / a / c);
		lattice[0][2] = c * cos(beta);
		lattice[2][2] = c * sin(beta);
	}
	
	if(mode==3)//set_ortho(lattice, metric);
	{
		double a, b, c;
		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = a;
		lattice[1][1] = b;
		lattice[2][2] = c;
	}

	if(mode==4)//set_tetra(lattice, metric);
	{
		double a, b, c;
		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = (a + b) / 2;
		lattice[1][1] = (a + b) / 2;
		lattice[2][2] = c;
	}
	
	if(mode==5)//set_rhomb(lattice, metric);
	{
		if(R_flag)
		{
			double a, b, c, angle, ahex, chex;

			a = sqrt(metric[0][0]);
			b = sqrt(metric[1][1]);
			c = sqrt(metric[2][2]);
			angle = acos((metric[0][1] / a / b +
						  metric[0][2] / a / c +
						  metric[1][2] / b / c) / 3);

			/* Reference, http://cst-www.nrl.navy.mil/lattice/struk/rgr.html */
			ahex = 2 * (a+b+c)/3 * sin(angle / 2);
			chex = (a+b+c)/3 * sqrt(3 * (1 + 2 * cos(angle))) ;
			lattice[0][0] = ahex / 2;
			lattice[1][0] = -ahex / (2 * sqrt(3));
			lattice[2][0] = chex / 3;
			lattice[1][1] = ahex / sqrt(3);
			lattice[2][1] = chex / 3;
			lattice[0][2] = -ahex / 2;
			lattice[1][2] = -ahex / (2 * sqrt(3));
			lattice[2][2] = chex / 3;
		}
		else//set_trigo(lattice, metric);
		{
			double a, b, c;

			a = sqrt(metric[0][0]);
			b = sqrt(metric[1][1]);
			c = sqrt(metric[2][2]);
			lattice[0][0] = (a + b) / 2;
			lattice[0][1] = - (a + b) / 4;
			lattice[1][1] = (a + b) / 4 * sqrt(3);
			lattice[2][2] = c;
		}
	}
	if(mode==6)//set_trigo(lattice, metric);
	{
		double a, b, c;

		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = (a + b) / 2;
		lattice[0][1] = - (a + b) / 4;
		lattice[1][1] = (a + b) / 4 * sqrt(3);
		lattice[2][2] = c;
	}

	if(mode==7)//set_cubic(lattice, metric);
	{
		double a, b, c;

		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = (a + b + c) / 3;
		lattice[1][1] = (a + b + c) / 3;
		lattice[2][2] = (a + b + c) / 3;
	}
}

void QSPG_setbox(QB_tools* output,QB_tools* input,double lattice[3][3],int hall_number,double symprec)
{
	int i,j;
	SPG_tools spg;
	spg.num=0;
	QB2SPG(input,&spg);
	//reset lattice
	for(i=0;i<3;i++)
	for(j=0;j<3;j++)
		spg.mat[i][j]=lattice[i][j];

	SPG2QB(&spg,output);
	
	if(spg.num!=0)
		SPG_free(&spg);
}

void QSPG_reshapebox(QB_tools* output,QB_tools* input,int a1,int x1,int a2,int x2)
{
	int i,j;
	int a3,x3;
	SPG_tools spg;
	spg.num=spg.ele_n=0;
	double temp[3];
	QB2SPG(input,&spg);
	double a,b,c,ab,bb,ac,bc,cc;
	double vb[3];
	double vc[3];
	a=sqrt(spg.mat[0][0]*spg.mat[0][0]+spg.mat[1][0]*spg.mat[1][0]+spg.mat[2][0]*spg.mat[2][0]);
	b=sqrt(spg.mat[0][1]*spg.mat[0][1]+spg.mat[1][1]*spg.mat[1][1]+spg.mat[2][1]*spg.mat[2][1]);
	ab=(spg.mat[0][0]*spg.mat[0][1]+spg.mat[1][0]*spg.mat[1][1]+spg.mat[2][0]*spg.mat[2][1])/a;
	for(i=0;i<3;i++)
		vb[i]=spg.mat[i][1]-spg.mat[i][0]*ab/a;
	bb=sqrt(vb[0]*vb[0]+vb[1]*vb[1]+vb[2]*vb[2]);
	
	ac=(spg.mat[0][0]*spg.mat[0][2]+spg.mat[1][0]*spg.mat[1][2]+spg.mat[2][0]*spg.mat[2][2])/a;
	bc=(        vb[0]*spg.mat[0][2]        +vb[1]*spg.mat[1][2]        +vb[2]*spg.mat[2][2])/bb;
	
	for(i=0;i<3;i++)
		vc[i]=spg.mat[i][2]-spg.mat[i][0]*ac/a-vb[i]*bc/bb;

	if(spg.mat[0][2]*(spg.mat[1][0]*spg.mat[2][1]-spg.mat[2][0]*spg.mat[1][1])+
	   spg.mat[1][2]*(spg.mat[2][0]*spg.mat[0][1]-spg.mat[0][0]*spg.mat[2][1])+
	   spg.mat[2][2]*(spg.mat[0][0]*spg.mat[1][1]-spg.mat[1][0]*spg.mat[0][1])>0)
		cc=sqrt(vc[0]*vc[0]+vc[1]*vc[1]+vc[2]*vc[2]);
	else
		cc=-sqrt(vc[0]*vc[0]+vc[1]*vc[1]+vc[2]*vc[2]);
	for(i=0;i<3;i++)
	{
		if(i!=a1&&i!=a2)
			a3=i;
		if(i!=x1&&i!=x2)
			x3=i;
	}
			   
	spg.mat[x1][a1] = a;
	spg.mat[x1][a2] = ab;
	spg.mat[x1][a3] = ac;
	spg.mat[x2][a1] = 0;
	spg.mat[x2][a2] = bb;
	spg.mat[x2][a3] = bc;
	spg.mat[x3][a1] = 0;
	spg.mat[x3][a2] = 0;
	spg.mat[x3][a3] = cc;
	
	for(i=0;i<spg.num;i++)
	{
		temp[0]=spg.pos[i][0];
		temp[1]=spg.pos[i][1];
		temp[2]=spg.pos[i][2];
		spg.pos[i][0]=temp[a1];
		spg.pos[i][1]=temp[a2];
		spg.pos[i][2]=temp[a3];
	}

	SPG2QB(&spg,output);
	if(spg.num!=0)
		SPG_free(&spg);
}