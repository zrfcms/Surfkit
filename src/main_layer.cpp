/*===============================================================
*   Copyright[c] 2022-2023, Z. R. Liu and R. F. Zhang
*
*   This file is part of the
*   Surfkit - An atomic toolkit for surface modelling with molecular adsorption
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*================================================================*/
#include"./QB/QB.h"
#include"./QSPG/QSPG.h"
#include"./SK/SurfaceKit.h"

int main_layer(int argc,char*argv[])
{
	if(strcmp(argv[1],"--space")==0)
	{
		QB_tools input1;
		QB_init(&input1); 
		QSPG_reshapebox(&input1,&input1,0,0,1,1);
		QB_read_file(&input1,argv[2]);
		double h=QB_checkdat(argv[4]);
		double r=4.0;
		int max_step=1000;
		if(argc<4)
		{
			printf(">>ERROR: Not enough input for adjusting space between layers\n");
			printf(">>Please use '--space [file(Slab)] [file2(Output)] [Height] {[Max step] [Distance/auto]}' instead\n");
			exit(0);
		}
		if(argc>5)
			max_step=(int)QB_checkdat(argv[5]);
		if(argc<7)
		{
			r=0.5*SK_rigid_average_r(&input1,4.0,input1.startz,h)
			 +0.5*SK_rigid_average_r(&input1,4.0,h,input1.startz+input1.boundz);
		}
		else if(!strcmp(argv[6],"auto"))
		{
			double r1=0.5*SK_rigid_average_r(&input1,4.0,h,input1.startz+input1.boundz);
			double r2=0.5*SK_rigid_average_r(&input1,4.0,input1.startz,h);
			printf("Average bond length is %lf\n",r);
			r=r1+r2;
		}
		else 
			r=QB_checkdat(argv[6]);
		SK_rigid_alias(&input1,h,r,max_step,10e-4);
		QB_dump_vasp_Cartesian(&input1,argv[3]);
		return 1;
	}
	if(strcmp(argv[1],"--match")==0)
	{
		if(argc<4)
		{
			printf(">>ERROR: Not enough input for searching twist angle\n");
			printf(">>Please use '–-match [file1(Slab)] [file2(Slab)] {[limit] [parameter]}…' instead\n");
			exit(0);
		}
		printf(">>Matching layers please wait\n");
		QB_tools input1;
		QB_init(&input1); 
		QB_tools input2;
		QB_init(&input2);
		QB_read_file(&input1,argv[2]);
		QB_read_file(&input2,argv[3]);
	
		SK_scaned_data* rotate_list;
		QSPG_reshapebox(&input1,&input1,0,0,1,1);
		QSPG_reshapebox(&input2,&input2,0,0,1,1);
		
		double limitd=20.0;
		double limitdmin=0.0;
		double limitbox=40.0;
		double limitdt=20;
		double mismatch=0.005;
		double mingamma=20;
		double maxgamma=160;
		double maxdangle=1;
		for(int i=4;i<argc-1;i+=2)
		{
			if(strcmp(argv[i],"maln")==0)
				limitd=QB_checkdat(argv[i+1]);
			if(strcmp(argv[i],"miln")==0)
				limitdmin=QB_checkdat(argv[i+1]);
			if(strcmp(argv[i],"mbox")==0)
				limitbox=QB_checkdat(argv[i+1]);
			if(strcmp(argv[i],"mtwi")==0)
				limitdt=QB_checkdat(argv[i+1]);
			if(strcmp(argv[i],"mvec")==0)
				mismatch=QB_checkdat(argv[i+1]);
			if(strcmp(argv[i],"mang")==0)
				maxdangle=QB_checkdat(argv[i+1]);
			if(strcmp(argv[i],"magm")==0)
				maxgamma=QB_checkdat(argv[i+1]);
			if(strcmp(argv[i],"migm")==0)
				mingamma=QB_checkdat(argv[i+1]);
		}
		SK_mismatch_set(limitd,limitdmin,limitbox,limitdt,mingamma,maxgamma,mismatch,maxdangle);
		int N_list=SK_mismatch_scan(input1.mat,input2.mat,&rotate_list);
	
		printf(">>%d Pairs of Twist Interface Detected\n",N_list);
		printf(">>Result Saved to mismatch.csv\n");
		SK_layerdata_save("mismatch.csv",N_list,rotate_list);
		return 1;
	}
	if(strcmp(argv[1],"--twist")==0)
	{
		if(argc<6)
		{
			printf(">>ERROR: Not enough input for build twist bilayer mode\n");
			printf(">>Please use '--twist [file1(Slab)] [file2(Slab)] [file3(Output)] [twist angle]' instead\n");
			exit(0);
		}
		printf(">>Merging twisted layers please wait\n");
		int N_list=0;
		SK_scaned_data* rotate_list;
		SK_layerdata_load("mismatch.csv",&N_list,&rotate_list);
		int selected=0;
		double twist_angle=QB_checkdat(argv[5]);
		double selected_angle=fabs(twist_angle-rotate_list[0].theta1+rotate_list[0].theta2);
		double current_angle;
		for(int i=1;i<N_list;i++)
		{
			current_angle=fabs(twist_angle-rotate_list[i].theta1+rotate_list[i].theta2);
			if(selected_angle>current_angle)
			{
				selected_angle=current_angle;
				selected=i;
			}
		}
		QB_tools input1;
		QB_init(&input1); 
		QB_tools input2;
		QB_init(&input2);
		QB_read_file(&input1,argv[2]);
		QB_read_file(&input2,argv[3]);
	
		QSPG_reshapebox(&input1,&input1,0,0,1,1);
		QSPG_reshapebox(&input2,&input2,0,0,1,1);
		
		SPG_tools spg1;
		spg1.num=spg1.ele_n=0;
		SPG_tools spg2;
		spg2.num=spg2.ele_n=0;
		
		double R_axis[3]={0,0,1};
		
		QB2SPG(&input1,&spg1);
		SK_redefine(&spg1,rotate_list[selected].trans1,10e-4);
		SK_rotate_matrix(spg1.mat,R_axis,rotate_list[selected].theta1);
			
		QB2SPG(&input2,&spg2);
		SK_redefine(&spg2,rotate_list[selected].trans2,10e-4);
		SK_rotate_matrix(spg2.mat,R_axis,rotate_list[selected].theta2);
		
		double ave;
		double lenth[2];
		for(int j=0;j<2;j++)
		{
			lenth[0]=sqrt(spg1.mat[0][j]*spg1.mat[0][j]+spg1.mat[1][j]*spg1.mat[1][j]);
			lenth[1]=sqrt(spg2.mat[0][j]*spg2.mat[0][j]+spg2.mat[1][j]*spg2.mat[1][j]);
			for(int k=0;k<2;k++)
			{
				ave=(spg1.mat[k][j]/lenth[0]+spg2.mat[k][j]/lenth[1])*0.25;
				spg1.mat[k][j]=ave*(lenth[0]+lenth[1]);
				spg2.mat[k][j]=ave*(lenth[0]+lenth[1]);
			}
		}
		
		SPG2QB(&spg1,&input1);
		SPG2QB(&spg2,&input2);
		
		for(int j=0;j<input2.TotalNumber;j++)
		{
			char elename[1024];
			QB_load_elename(&input2,input2.atom[j].type-1,elename);
			int k=QB_get_elename(&input1,elename);
			QB_create_atom(&input1,k,input2.atom[j].x+input1.mat[2][0],
									   input2.atom[j].y+input1.mat[2][1],
									   input2.atom[j].z+input1.mat[2][2]);
		}
		input1.mat[2][0]+=input2.mat[2][0];
		input1.mat[2][1]+=input2.mat[2][1];
		input1.mat[2][2]+=input2.mat[2][2];
		
		QB_dump_vasp_Cartesian(&input1,argv[4]);
		printf(">>Merging twisted layers complete! File dumped to %s\n",argv[4]);
		return 1;
	}
	//test function
	//simple merge layers together
	if(strcmp(argv[1],"--merge")==0)
	{
		QB_tools input1;
		QB_init(&input1); 
		QB_tools input2;
		QB_init(&input2);
		QB_read_file(&input1,argv[2]);
		QB_read_file(&input2,argv[3]);
	
		QSPG_reshapebox(&input1,&input1,0,0,1,1);
		QSPG_reshapebox(&input2,&input2,0,0,1,1);
		
		SPG_tools spg1;
		spg1.num=spg1.ele_n=0;
		SPG_tools spg2;
		spg2.num=spg2.ele_n=0;
		
		QB2SPG(&input1,&spg1);
		QB2SPG(&input2,&spg2);

		double ave;
		double lenth[2];
		for(int j=0;j<2;j++)
		for(int k=0;k<2;k++)
		{
			lenth[0]=sqrt(spg1.mat[j][0]*spg1.mat[j][0]+spg1.mat[j][1]*spg1.mat[j][1]);
			lenth[1]=sqrt(spg2.mat[j][0]*spg2.mat[j][0]+spg2.mat[j][1]*spg2.mat[j][1]);
			ave=(spg1.mat[j][k]/lenth[0]+spg2.mat[j][k]/lenth[1])*0.25;
			spg1.mat[j][k]=ave*(lenth[0]+lenth[1]);
			spg2.mat[j][k]=ave*(lenth[0]+lenth[1]);
		}
		
		SPG2QB(&spg1,&input1);
		SPG2QB(&spg2,&input2);
		
		for(int j=0;j<input2.TotalNumber;j++)
		{
			char elename[1024];
			QB_load_elename(&input2,input2.atom[j].type-1,elename);
			int k=QB_get_elename(&input1,elename);
			QB_create_atom(&input1,k,input2.atom[j].x+input1.mat[2][0],
									   input2.atom[j].y+input1.mat[2][1],
									   input2.atom[j].z+input1.mat[2][2]);
		}
		
		input1.mat[2][0]+=input2.mat[2][0];
		input1.mat[2][1]+=input2.mat[2][1];
		input1.mat[2][2]+=input2.mat[2][2];
		
		QB_dump_vasp_Cartesian(&input1,argv[4]);
		return 1;
	}
	return 0;
}
