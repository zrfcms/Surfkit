#include"QB/QB.h"
#include"QSPG/QSPG.h"
#include"SK/SurfaceKit.h"
#include"./QSPG/spglib/mathfunc.h"
int main_slab(int argc,char* argv[])
{
	char alphabet[30]="abcdefghijklmnopqrstuvwxyz";
	//get layer list 
	if(strcmp(argv[1],"--shift")==0)
	{
		if(argc<5)
		{
			printf(">>EROOR: Not enough input for shifting cell mode\n");
			printf(">>Please use '--shift file(INPOS) file(OUTPOS) [shift height]' instead\n");
			exit(0);
		}
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		double shift_height=QB_checkdat(argv[4]);
		double vec[2];
		if(input.mat[2][2]<0)
		{
			vec[0]=-input.mat[2][0]/input.mat[2][2];
			vec[1]=-input.mat[2][1]/input.mat[2][2];
		}
		else if(input.mat[2][2]>0)
		{
			vec[0]=input.mat[2][0]/input.mat[2][2];
			vec[1]=input.mat[2][1]/input.mat[2][2];
		}
		else 
		{
			printf(">>EROOR: C axis in xy plane,Please align model first\n");
			printf(">>Please use '--align file1(INPOS) file2(OUTPOS) [abc] [xyz] [abc] [xyz]' first\n");
			exit(0);
		}
		for(int i=0;i<input.TotalNumber;i++)
		{
			input.atom[i].x+=vec[0]*shift_height;
			input.atom[i].y+=vec[1]*shift_height;
			input.atom[i].z+=shift_height;
			QB_wrap(&input,i);
		}
		QB_dump_vasp_Cartesian(&input,argv[3]);
		return 1;
	}
	//get layer list 
	if(strcmp(argv[1],"--scan")==0)
	{
		if(argc<5)
		{
			printf(">>EROOR: Not enough input for scan atomic layer and non-polar center mode\n");
			printf(">>Please use '--scan file(INPOS) [min_height] [max_height]' instead\n");
			exit(0);
		}
		printf(">>Scanning atomic layer please wait\n");
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		double min_height=QB_checkdat(argv[3]);
		double max_height=QB_checkdat(argv[4]);
		if(input.mat[2][2]<0)
		{
			for(int i=0;i<3;i++)
				input.mat[2][i]*=-1;
			for(int i=0;i<input.TotalNumber;i++)
				QB_wrap(&input,i);
		}
		//switch to Fractional coordinate system
		SPG_tools spg;
		spg.num=spg.ele_n=0;
		QB2SPG(&input,&spg);

		int list_n=0;
		SK_layer_list *list;
		
		SK_polar_center(&list_n,&list,min_height,max_height,&spg,10e-4); 
		SK_atomic_layers(&list_n,&list,min_height,max_height,&input,10e-4);
		SK_layer_sort(list_n,list);
		SK_layer_dump("Layer_Scan.csv",list_n,list);
		return 1;		
	}
	//get unit cell 
	if(strcmp(argv[1],"--conv")==0)
	{
		if(argc<4)
		{
			printf(">>EROOR: Not enough input for conventional cell mode\n");
			printf(">>Please use '--conv file1(INPOS) file2(OUTPOS)' instead\n");
			exit(0);
		}
		printf(">>Finding conventional cell please wait\n");
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		
		QSPG_refined(&input,&input,10e-4);
			
		QB_dump_vasp_Cartesian(&input,argv[3]);
		QB_free_atom(&input);
		printf(">>>Finding conventional cell complete! File dumped to %s\n",argv[3]);
		return 1;
	}
	
	//get unit cell in plane
	if(strcmp(argv[1],"--psurf")==0)
	{
		if(argc<5)
		{
			printf(">>EROOR: Not enough input for primitive slab mode\n");
			printf(">>Please use '--psurf file1(INPOS) file2(OUTPOS) [ab ac bc]' instead\n");
			exit(0);
		}
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		
		//switch to Fractional coordinate system
		SPG_tools spg;
		spg.num=spg.ele_n=0;
		QB2SPG(&input,&spg);
		
		if(strcmp(argv[4],"ab")==0)
			SK_SPG_minimize_face(&spg,0,1,10e-4);
		else if(strcmp(argv[4],"ac")==0)
			SK_SPG_minimize_face(&spg,0,2,10e-4);
		else if(strcmp(argv[4],"bc")==0)
			SK_SPG_minimize_face(&spg,1,2,10e-4);
		else 
		{
			printf("WARNING: plane unselected, use ab as default\n");
			SK_SPG_minimize_face(&spg,0,1,10e-4);
		}
		SPG2QB(&spg,&input);
		QB_dump_vasp_Cartesian(&input,argv[3]);
		SPG_free(&spg);
		QB_free_atom(&input);
		return 1;
	}
	
	//get unit cell in one direction
	if(strcmp(argv[1],"--pslab")==0)
	{
		if(argc<5)
		{
			printf(">>ERROR: Not enough input for primitive slab mode\n");
			printf(">>Please use '--pslab file1(INPOS) file2(OUTPOS) [a b c]' instead\n");
			exit(0);
		}
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		
		//switch to Fractional coordinate system
		SPG_tools spg;
		spg.num=spg.ele_n=0;
		QB2SPG(&input,&spg);
		
		if(strcmp(argv[4],"c")==0)
			SK_SPG_minimize_cylinder(&spg,0,1,10e-4);
		else if(strcmp(argv[4],"b")==0)
			SK_SPG_minimize_cylinder(&spg,0,2,10e-4);
		else if(strcmp(argv[4],"a")==0)
			SK_SPG_minimize_cylinder(&spg,1,2,10e-4);
		else 
		{
			printf("WARNING: plane unselected, use c as default\n");
			SK_SPG_minimize_cylinder(&spg,0,1,10e-4);
		}
		SPG2QB(&spg,&input);
		QB_dump_vasp_Cartesian(&input,argv[3]);
		SPG_free(&spg);
		QB_free_atom(&input);
		return 1;
	}
	
	//get unit cell in one direction
	if(strcmp(argv[1],"--guass")==0)
	{
		if(argc<4)
		{
			printf(">>ERROR: Not enough input for Gussian reduction slab mode\n");
			printf(">>Please use '--guass file1(INPOS) file2(OUTPOS)' instead\n");
			exit(0);
		}
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		
		Gaussian_lattice_reduction(input.mat[0],input.mat[1]);
		for(int i=0;i<input.TotalNumber;i++)
			QB_wrap(&input,i);
		
		QB_dump_vasp_Cartesian(&input,argv[3]);
		QB_free_atom(&input);
		return 1;
	}
	
	//get unit cell in one direction
	if(strcmp(argv[1],"--prim2D")==0)
	{
		if(argc<4)
		{
			printf(">>ERROR: Not enough input for 2D primitive slab mode\n");
			printf(">>Please use '--prim2D file1(INPOS) file2(OUTPOS)' instead\n");
			exit(0);
		}
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		
		
		//switch to Fractional coordinate system
		SPG_tools spg;
		spg.num=spg.ele_n=0;
		QB2SPG(&input,&spg);
		SK_SPG_minimize_face(&spg,0,1,10e-4);	
		SK_SPG_minimize_cylinder(&spg,0,1,10e-4);
		SPG2QB(&spg,&input);
		Gaussian_lattice_reduction(input.mat[0],input.mat[1]);
		Right_hand_reduction(input.mat);
		for(int i=0;i<input.TotalNumber;i++)
			QB_wrap(&input,i);
		
		QB_dump_vasp_Cartesian(&input,argv[3]);
		QB_free_atom(&input);
		return 1;
	}
	//redefine slab by projection
	if(strcmp(argv[1],"--proj")==0)
	{
		if(argc<8)
		{
			printf(">>ERROR: Not enough input for projection mode\n");
			printf(">>Please use '--proj -auto [file1(In)] [file2(Out)] ch ck cl' instead\n");
			exit(0);
		}
		
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[3]);
		
		//switch to Fractional coordinate system
		SPG_tools spg;
		spg.num=spg.ele_n=0;
		QB2SPG(&input,&spg);
		
		int rot[3][3]={1,0,0,0,1,0,0,0,1};
		if(strcmp(argv[2],"-auto")==0)
		{
			if(argc<8)
			{
				printf(">>ERROR: Not enough input for projection mode\n");
				printf(">>Please use '--proj -auto [file1(In)] [file2(Out)] ch ck cl' instead\n");
				exit(0);
			}
			int c[3]={(int)QB_checkdat(argv[5]),(int)QB_checkdat(argv[6]),(int)QB_checkdat(argv[7])};
			if(c[0]==0&&c[1]==0&&c[2]==0)
			{
				printf(">>EROOR: Vector length can not be zero\n");
				exit(0);
			}
			printf(">>Redefining crystal cell please wait\n");
			SK_auto_Frot(c,rot);
			printf(">>Automatically Generated Matrix:\n");
			for(int i=0;i<3;i++)
				printf("%d\t%d\t%d\n",rot[0][i],rot[1][i],rot[2][i]);
		}
		if(strcmp(argv[2],"-mat")==0)
		{
			if(argc<14)
			{
				printf(">>ERROR: Not enough input for projection mode\n");
				printf(">>Please use '--proj -mat [file1(In)] [file2(Out)] ah ak al bh bk bl ch ck cl' instead\n");
				exit(0);
			}
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++)
					rot[j][i]=(int)QB_checkdat(argv[5+3*i+j]);
			if((rot[0][0]*rot[1][1]*rot[2][2]+rot[0][1]*rot[1][2]*rot[2][0]+rot[0][2]*rot[1][0]*rot[2][1]-
				rot[0][0]*rot[1][2]*rot[2][1]-rot[0][2]*rot[1][1]*rot[2][0]-rot[0][1]*rot[1][0]*rot[2][2])==0)
			{
				printf(">>ERROR: Rank of matrix can not be zero\n");
				exit(0);
			}
			printf(">>Redefining crystal cell please wait\n");
		}
		SK_redefine(&spg,rot,10e-4);
		SPG2QB(&spg,&input);
		QB_dump_vasp_Cartesian(&input,argv[4]);
		printf(">>>Redefining crystal cell complete! File dumped to %s\n",argv[4]);
		SPG_free(&spg);
		QB_free_atom(&input);
		return 1;
	}
	
	//duplicate slab
	if(strcmp(argv[1],"--dup")==0)
	{
		if(argc<7)
		{
			printf(">>ERROR: Not enough input for projection mode\n");
			printf(">>Please use '-–dup [file1(In)] [file2(Out)] ah bk cl' instead\n");
			exit(0);
		}
		
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		
		int n[3]={(int)QB_checkdat(argv[4]),(int)QB_checkdat(argv[5]),(int)QB_checkdat(argv[6])};
		if(n[0]<=0||n[1]<=0||n[2]<=0)
		{
			printf(">>WARNING: Volume of new system can not be zero\n");
			printf(">>Duplicate magnification will be set to at least 1 automatically\n");
			if(n[0]<=0)n[0]=1;
			if(n[1]<=0)n[1]=1;
			if(n[2]<=0)n[2]=1;
		}
		printf(">>Redefining crystal cell please wait\n");
		QB_duplicate(&input,n[0],n[1],n[2]);
		
		QB_dump_vasp_Cartesian(&input,argv[3]);
		printf(">>>Redefining crystal cell complete! File dumped to %s\n",argv[3]);
		return 1;
	}
	
	//redefine slab thickness 
	if(strcmp(argv[1],"--cleav")==0)
	{
		if(argc<5)
		{
			printf(">>ERROR: Not enough input for thickness mode\n");
			printf(">>Please use '--cleav file1(INPOS) file2(OUTPOS) [min height] [max height]' instead\n");
			exit(0);
		}
		printf(">>Cleaving crystal slab please wait\n");
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		double minh=QB_checkdat(argv[4]);
		double maxh=QB_checkdat(argv[5]);
		double xproduct[3];
		double xin1[3]={input.mat[0][0],input.mat[0][1],input.mat[0][2]};	
		double xin2[3]={input.mat[1][0],input.mat[1][1],input.mat[1][2]};
		SK_cross_product(xproduct,xin1,xin2);
		double normal_c=sqrt(xproduct[0]*xproduct[0]
							+xproduct[1]*xproduct[1]
							+xproduct[2]*xproduct[2]);
		for(int i=0;i<3;i++)
			xproduct[i]/=normal_c;
		double heightc=xproduct[0]*input.mat[2][0]
					  +xproduct[1]*input.mat[2][1]
					  +xproduct[2]*input.mat[2][2];
		
		if(heightc<0)
		{
			for(int i=0;i<3;i++)
				xproduct[i]=-xproduct[i];
			heightc=-heightc;
		}
		
		int dc=floor(maxh/heightc)-floor(minh/heightc)+1;
		double dz=floor(minh/heightc);
		QB_duplicate(&input,1,1,dc);
		
		minh-=(double)dz*heightc;
		maxh-=(double)dz*heightc;
		
		for(int i=0;i<input.TotalNumber;i++)
		{
			double trace_c=(input.atom[i].x-input.zerox)*xproduct[0]
						  +(input.atom[i].y-input.zeroy)*xproduct[1]
						  +(input.atom[i].z-input.zeroz)*xproduct[2];
			if(trace_c<minh)
				QB_delete_atom(&input,i);
			if(trace_c>maxh)
				QB_delete_atom(&input,i);
		}
		
		for(int i=0;i<3;i++)
			input.mat[2][i]=(double)dc*heightc*xproduct[i];
		
		for(int i=0;i<input.TotalNumber;i++)
			QB_wrap(&input,i);
		
		QB_dump_vasp_Cartesian(&input,argv[3]);
		printf(">>Cleaving crystal slab complete! File dumped to %s\n",argv[3]);
		return 1;
	}
	
	//add void layer for surface
	if(strcmp(argv[1],"--vacuum")==0)
	{
		if(argc<5)
		{
			printf(">>ERROR: Not enough input for vacuum layer mode\n");
			printf(">>Please use '--vacuum file1(INPOS) file2(OUTPOS) [thickness]' instead\n");
			exit(0);
		}
		printf(">>Cleaving crystal slab please wait\n");
		
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		
		double Thickness=QB_checkdat(argv[4]);
		double xproduct[3];
		double xin1[3]={input.mat[0][0],input.mat[0][1],input.mat[0][2]};	
		double xin2[3]={input.mat[1][0],input.mat[1][1],input.mat[1][2]};
		SK_cross_product(xproduct,xin1,xin2);
		double normal_c=sqrt(xproduct[0]*xproduct[0]
							+xproduct[1]*xproduct[1]
							+xproduct[2]*xproduct[2]);
		for(int i=0;i<3;i++)
			xproduct[i]/=normal_c;
		double heightc=xproduct[0]*input.mat[2][0]
					  +xproduct[1]*input.mat[2][1]
					  +xproduct[2]*input.mat[2][2];
					  
		//switch to Fractional coordinate system
		SPG_tools spg;
		spg.num=spg.ele_n=0;
		QB2SPG(&input,&spg);
		
		//get min & max fractional coordinate of c
		double min=10,max=-10;
		for(int i=0;i<spg.num;i++)
		{
			if(min>spg.pos[i][2])min=spg.pos[i][2];
			if(max<spg.pos[i][2])max=spg.pos[i][2];
		}
		double disp[3]={(-min*heightc+0.5*Thickness)*xproduct[0],
						(-min*heightc+0.5*Thickness)*xproduct[1],
						(-min*heightc+0.5*Thickness)*xproduct[2]};
		for(int i=0;i<input.TotalNumber;i++)
		{
			input.atom[i].x+=disp[0];
			input.atom[i].y+=disp[1];
			input.atom[i].z+=disp[2];
		}
		
		for(int i=0;i<3;i++)
			input.mat[2][i]=(Thickness+(max-min)*heightc)*xproduct[i];
		
		for(int i=0;i<input.TotalNumber;i++)
			QB_wrap(&input,i);
		
		QB_dump_vasp_Cartesian(&input,argv[3]);
		printf(">>Cleaving crystal slab complete! File dumped to %s\n",argv[3]);
		return 1;
	}
	//function in test
	if(strcmp(argv[1],"--align")==0)
	{
		if(argc<8)
		{
			printf(">>ERROR: Not enough input for align layer mode\n");
			printf(">>Please use '--align file1(INPOS) file2(OUTPOS) [abc] [xyz] [abc] [xyz]' instead\n");
			exit(0);
		}
		
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		
		int a1[2]={0,1};
		int x1[2]={0,1};
		if(strcmp("a",argv[4])==0)a1[0]=0;
		if(strcmp("b",argv[4])==0)a1[0]=1;
		if(strcmp("c",argv[4])==0)a1[0]=2;
		if(strcmp("x",argv[5])==0)x1[0]=0;
		if(strcmp("y",argv[5])==0)x1[0]=1;
		if(strcmp("z",argv[5])==0)x1[0]=2;
		if(strcmp("a",argv[6])==0)a1[1]=0;
		if(strcmp("b",argv[6])==0)a1[1]=1;
		if(strcmp("c",argv[6])==0)a1[1]=2;
		if(strcmp("x",argv[7])==0)x1[1]=0;
		if(strcmp("y",argv[7])==0)x1[1]=1;
		if(strcmp("z",argv[7])==0)x1[1]=2;
		
		if(a1[0]==a1[1])
		{
			printf(">>WARNING: Same vector selected for align layer\n");
			a1[1]=(a1[0]+1)%3;
		}
		
		if(x1[0]==x1[1])
		{
			printf(">>WARNING: Same vector selected for align layer\n");
			x1[1]=(x1[0]+1)%3;
		}
		
		QSPG_reshapebox(&input,&input,a1[0],x1[0],a1[1],x1[1]);
		QB_dump_vasp_Cartesian(&input,argv[3]);
		return 1;
	}
	return 0;
}