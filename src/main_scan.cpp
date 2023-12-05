#include"QB/QB.h"
#include"QSPG/QSPG.h"
#include"SK/SurfaceKit.h"
#include"./QSPG/spglib/mathfunc.h"
void SKA_dump_status(QB_tools* QB)
{
	char name[1024];
	FILE* status=fopen("Status.csv","w");
	printf(">>Total %d Type of Sites found\n",QB->TypeNumber);
	fprintf(status,"Type,Name,Sum\n");
	for(int i=0;i<QB->TypeNumber;i++)
	{
		QB_load_elename(QB,i,name);
		int number=0;
		for(int j=0;j<QB->TotalNumber;j++)
		{
			if(QB->atom[j].type==i+1)
				number++;
		}
		printf(">>Type %d named %s\t%d Sites in total\n",i+1,name,number);
		fprintf(status,"%d,%s,%d\n",i+1,name,number);
	}
	fclose(status);
}
void SKA_list_sort(SK_site_list& SKA_list)
{
	int i,j;
	double td;
	int ti;
	for(i=0;i<SKA_list.num-1;i++)
	for(j=0;j<SKA_list.num-i-1;j++)
	if(SKA_list.wyckoff[j]>SKA_list.wyckoff[j+1])
	{
		td=SKA_list.pos[j][0];
		SKA_list.pos[j][0]=SKA_list.pos[j+1][0];
		SKA_list.pos[j+1][0]=td;
		
		td=SKA_list.pos[j][1];
		SKA_list.pos[j][1]=SKA_list.pos[j+1][1];
		SKA_list.pos[j+1][1]=td;
		
		td=SKA_list.pos[j][2];
		SKA_list.pos[j][2]=SKA_list.pos[j+1][2];
		SKA_list.pos[j+1][2]=td;
		
		ti=SKA_list.wyckoff[j];
		SKA_list.wyckoff[j]=SKA_list.wyckoff[j+1];
		SKA_list.wyckoff[j+1]=ti;
		
		ti=SKA_list.equivalent_id[j];
		SKA_list.equivalent_id[j]=SKA_list.equivalent_id[j+1];
		SKA_list.equivalent_id[j+1]=ti;
		
		ti=SKA_list.equivalent_num[j];
		SKA_list.equivalent_num[j]=SKA_list.equivalent_num[j+1];
		SKA_list.equivalent_num[j+1]=ti;
	}
}

int main_scan(int argc,char* argv[])
{
	char alphabet[30]="abcdefghijklmnopqrstuvwxyz";
	int i,j;
	if(strcmp(argv[1],"--site2D")==0)
	{
		if(argc<6)
		{
			printf(">>ERROR: Not enough input for adsorption scan mode\n");
			printf(">>Please use '--site2D [file1(Slab)] [file2(Site List)] [Distance(2D)] [Max multiplicity]' instead\n");
			exit(0);
		}
		printf(">>Scanning adsorption please wait\n");
		SK_site_list SKA_list;
		QB_tools SKA_site;
		QB_init(&SKA_site);
		
		SKA_list.num=0;
		SKA_site.TotalNumber=0;
		
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		QB_copy_box(&SKA_site,&input);
		
		SPG_tools spg;
		spg.num=spg.ele_n=0;
		QB2SPG(&input,&spg);
		
		SPG_tools temp_spg;
		temp_spg.num=temp_spg.ele_n=0;
		int origin_id=0;
		//get Symmetry from spglib
		SpglibDataset *Dataset=spg_get_dataset(spg.mat,spg.pos,spg.type,spg.num,10e-4);
		int site_rot[192][3][3];
		double site_trans[192][3];
		int site_sym=spg_get_symmetry_from_database(site_rot,site_trans,Dataset->hall_number);
		
		//get site list from surfkit
		SK_get_sitelist_2D(&SKA_list,0,1,site_sym,site_rot,site_trans,Dataset->std_lattice,Dataset->hall_number,10e-4);
		SKA_list_sort(SKA_list);
		
		//get absorbtion list from sitelist
		int j,k,l;
		for(j=0;j<3;j++)
		for(k=0;k<3;k++)
		{
			temp_spg.mat[j][k]=Dataset->std_lattice[j][k];
		}
		int flag_num=-2;
		int meq=(int)QB_checkdat(argv[5]);
		
		for(j=0;j<SKA_list.num;j++)
		if(meq>=SKA_list.equivalent_num[j])
		{
			extern void SPG_add_atom(SPG_tools* spg,int type,double x,double y,double z);
			SPG_add_atom(&temp_spg,SKA_list.equivalent_id[j],
			SKA_list.pos[j][0],
			SKA_list.pos[j][1],
			SKA_list.pos[j][2]);
			if(SKA_list.equivalent_id[j]>flag_num)
				flag_num=SKA_list.equivalent_id[j];
		}
		//get mapping list
		flag_num++;
		int n=1;
		int *flag_list;
		if(flag_num>0)
		{
			flag_list=(int*)malloc(flag_num*sizeof(int));
			for(j=0;j<flag_num;j++)flag_list[j]=0;
			for(j=0;j<SKA_list.num;j++)
			if(meq>=SKA_list.equivalent_num[j])
			{
				if(flag_list[SKA_list.equivalent_id[j]]==0)
				{
					flag_list[SKA_list.equivalent_id[j]]=n;
					char label[256];
					sprintf(label,"%c%d",alphabet[SKA_list.wyckoff[j]],SKA_list.equivalent_num[j]);
					QB_save_elename(&SKA_site,n-1,label);
					n++;
				}
			}
		}
		//get rotation matrix
		double tmp_mat[3][3];
		mat_inverse_matrix_d3(tmp_mat, Dataset->std_rotation_matrix, 10e-4);
		mat_multiply_matrix_d3(temp_spg.mat,tmp_mat,Dataset->std_lattice);
		
		//set anker position
		int std_id;
		for(j=0;j<Dataset->n_std_atoms;j++)
		if((Dataset->mapping_to_primitive[origin_id])==(Dataset->std_mapping_to_primitive[j]))
		{
			std_id=j;
			break;
		}
		double anker_pos[3]={Dataset->std_positions[std_id][0],Dataset->std_positions[std_id][1],Dataset->std_positions[std_id][2]};
		
		//Mapping site list to backup model
		double displace_s=QB_checkdat(argv[4]);
		double max_c=0;
		for(j=0;j<spg.num;j++)
		{
			if(spg.pos[j][2]>max_c)max_c=spg.pos[j][2];
		}
		
		double axb[3];
		SK_cross_product(axb,input.mat[0],input.mat[1]);
		displace_s=displace_s*sqrt(axb[0]*axb[0]+axb[1]*axb[1]+axb[2]*axb[2])
				  /fabs(axb[0]*spg.mat[0][2]+axb[1]*spg.mat[1][2]+axb[2]*spg.mat[2][2]);
		max_c+=displace_s;
		
		double min_height=SK_min_height(Dataset->std_lattice);
		double max_height=sqrt(input.mat[0][0]*input.mat[0][0]
						  +input.mat[0][1]*input.mat[0][1]+input.mat[0][2]*input.mat[0][2])+
						  sqrt(input.mat[1][0]*input.mat[1][0]
						  +input.mat[1][1]*input.mat[1][1]+input.mat[1][2]*input.mat[1][2])+
						  sqrt(input.mat[2][0]*input.mat[2][0]
						  +input.mat[2][1]*input.mat[2][1]+input.mat[2][2]*input.mat[2][2]);
		
		for(j=0;j<temp_spg.num;j++)
		if(flag_list[temp_spg.type[j]])
		for(k=-0.5*max_height/min_height-1;k<0.5*max_height/min_height+1;k++)
		for(l=-0.5*max_height/min_height-1;l<0.5*max_height/min_height+1;l++)
		{
			double c_pos[3];
			double f_pos[3];
			c_pos[0]=input.atom[origin_id].x+(temp_spg.pos[j][0]+k-anker_pos[0])*temp_spg.mat[0][0]
												   +(temp_spg.pos[j][1]+l-anker_pos[1])*temp_spg.mat[0][1]
												   +max_c*spg.mat[0][2]-anker_pos[2]*temp_spg.mat[0][2];
											
			c_pos[1]=input.atom[origin_id].y+(temp_spg.pos[j][0]+k-anker_pos[0])*temp_spg.mat[1][0]
												   +(temp_spg.pos[j][1]+l-anker_pos[1])*temp_spg.mat[1][1]
												   +max_c*spg.mat[1][2]-anker_pos[2]*temp_spg.mat[1][2];
			
			c_pos[2]=input.atom[origin_id].z+(temp_spg.pos[j][0]+k-anker_pos[0])*temp_spg.mat[2][0]
												   +(temp_spg.pos[j][1]+l-anker_pos[1])*temp_spg.mat[2][1]
												   +max_c*spg.mat[2][2]-anker_pos[2]*temp_spg.mat[2][2];
						
			extern void SPG_add_atom(SPG_tools* spg,int type,double x,double y,double z);
			QB_c2f(&input,c_pos,f_pos);
				
			if(SK_SPG_overloop_vec(&spg,temp_spg.type[j],f_pos[0],f_pos[1],f_pos[2],10e-4)==0)
			{
				SPG_add_atom(&spg,SKA_list.equivalent_id[j],f_pos[0],f_pos[1],f_pos[2]);
				QB_create_atom(&SKA_site,flag_list[temp_spg.type[j]],c_pos[0],c_pos[1],c_pos[2]);
				QB_wrap(&SKA_site,SKA_site.TotalNumber-1);
			}
		}
		free(flag_list);
		QB_dump_vasp_Direct(&SKA_site,argv[3]);
		SKA_dump_status(&SKA_site);
		return 1;
	}
	
	if(strcmp(argv[1],"--site3D")==0)
	{
		if(argc<6)
		{
			printf(">>ERROR: Not enough input for adsorption scan mode\n");
			printf(">>Please use '–site3D [file1(Slab)] [file2(Site List)] [Distance(3D)] [Max multiplicity]' instead\n");
			exit(0);
		}
		printf(">>Scanning adsorption please wait\n");
		SK_site_list SKA_list;
		QB_tools SKA_site;
		QB_init(&SKA_site);
		
		SKA_list.num=0;
		SKA_site.TotalNumber=0;
		
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[2]);
		QB_copy_box(&SKA_site,&input);
		
		SPG_tools spg;
		spg.num=spg.ele_n=0;
		QB2SPG(&input,&spg);
		
		SPG_tools temp_spg;
		temp_spg.num=temp_spg.ele_n=0;
		
		//backtrack 3D lattice from 2D crystal
		int origin_id=SK_SPG_pick_3Dlattice(&spg,0,1,10e-4);
		
		//get symmetry from spglib
		SpglibDataset *Dataset=spg_get_dataset(spg.mat,spg.pos,spg.type,spg.num,10e-4);
		int site_rot[192][3][3];
		double site_trans[192][3];
		int site_sym=spg_get_symmetry_from_database(site_rot,site_trans,Dataset->hall_number);
		
		//get site list from surfkit
		SK_get_sitelist_3D(&SKA_list,site_sym,site_rot,site_trans,Dataset->std_lattice,Dataset->hall_number,10e-4);
		SKA_list_sort(SKA_list);
		int j,k,l,m;
		for(j=0;j<3;j++)
		for(k=0;k<3;k++)
		{
			temp_spg.mat[j][k]=Dataset->std_lattice[j][k];
		}
		
		//get absorbtion list from sitelist
		int flag_num=-2;
		int meq=(int)QB_checkdat(argv[5]);
		for(j=0;j<SKA_list.num;j++)
		if(meq>SKA_list.equivalent_num[j])
		{
			extern void SPG_add_atom(SPG_tools* spg,int type,double x,double y,double z);
			SPG_add_atom(&temp_spg,SKA_list.equivalent_id[j],
			SKA_list.pos[j][0],
			SKA_list.pos[j][1],
			SKA_list.pos[j][2]);
			if(SKA_list.equivalent_id[j]>flag_num)
				flag_num=SKA_list.equivalent_id[j];
		}
		int n=1;
		flag_num++;
		int *flag_list;
		//update site list to screen
		if(flag_num>0)
		{
			flag_list=(int*)malloc(flag_num*sizeof(int));
			for(j=0;j<flag_num;j++)flag_list[j]=0;
			for(j=0;j<SKA_list.num;j++)
			if(meq>SKA_list.equivalent_num[j])
			{
				if(flag_list[SKA_list.equivalent_id[j]]==0)
				{
					flag_list[SKA_list.equivalent_id[j]]=n;
					char label[256];
					sprintf(label,"%c%d",alphabet[SKA_list.wyckoff[j]],SKA_list.equivalent_num[j]);
					QB_save_elename(&SKA_site,n-1,label);
					n++;
				}
			}
		}
		//get rotation matrix
		double tmp_mat[3][3];
		mat_inverse_matrix_d3(tmp_mat, Dataset->std_rotation_matrix, 10e-4);
		mat_multiply_matrix_d3(temp_spg.mat,tmp_mat,Dataset->std_lattice);

		double min_height=SK_min_height(Dataset->std_lattice);
		double max_height=sqrt(input.mat[0][0]*input.mat[0][0]+
						  input.mat[0][1]*input.mat[0][1]+input.mat[0][2]*input.mat[0][2])+
						  sqrt(input.mat[1][0]*input.mat[1][0]+
						  input.mat[1][1]*input.mat[1][1]+input.mat[1][2]*input.mat[1][2])+
						  sqrt(input.mat[2][0]*input.mat[2][0]+
						  input.mat[2][1]*input.mat[2][1]+input.mat[2][2]*input.mat[2][2]);
		
		//set anker position		
		int std_id;
		for(j=0;j<Dataset->n_std_atoms;j++)
		if((Dataset->mapping_to_primitive[origin_id])==(Dataset->std_mapping_to_primitive[j]))
		{
			std_id=j;
			break;
		}
		double anker_pos[3]={Dataset->std_positions[std_id][0],Dataset->std_positions[std_id][1],Dataset->std_positions[std_id][2]};

		//get back 2D structure
		QB2SPG(&input,&spg);
		
		//Mapping site list to backup model
		double displace_s=QB_checkdat(argv[4]);
		double axb[3];
		SK_cross_product(axb,input.mat[0],input.mat[1]);
		displace_s=displace_s*sqrt(axb[0]*axb[0]+axb[1]*axb[1]+axb[2]*axb[2])
				  /fabs(axb[0]*spg.mat[0][2]+axb[1]*spg.mat[1][2]+axb[2]*spg.mat[2][2]);
		
		double max_c=0;
		for(j=0;j<spg.num;j++)
			if(spg.pos[j][2]>max_c)
				max_c=spg.pos[j][2];
				
		for(j=0;j<temp_spg.num;j++)
		if(flag_list[temp_spg.type[j]])
		for(k=-0.5*max_height/min_height-1;k<0.5*max_height/min_height+1;k++)
		for(l=-0.5*max_height/min_height-1;l<0.5*max_height/min_height+1;l++)
		for(m=-0.5*max_height/min_height-1;m<0.5*max_height/min_height+1;m++)
		{
			double c_pos[3];
			double f_pos[3];
			c_pos[0]=input.atom[origin_id].x+(temp_spg.pos[j][0]+k-anker_pos[0])*temp_spg.mat[0][0]
											+(temp_spg.pos[j][1]+l-anker_pos[1])*temp_spg.mat[0][1]
											+(temp_spg.pos[j][2]+m-anker_pos[2])*temp_spg.mat[0][2];
											
			c_pos[1]=input.atom[origin_id].y+(temp_spg.pos[j][0]+k-anker_pos[0])*temp_spg.mat[1][0]
											+(temp_spg.pos[j][1]+l-anker_pos[1])*temp_spg.mat[1][1]
											+(temp_spg.pos[j][2]+m-anker_pos[2])*temp_spg.mat[1][2];
			
			c_pos[2]=input.atom[origin_id].z+(temp_spg.pos[j][0]+k-anker_pos[0])*temp_spg.mat[2][0]
											+(temp_spg.pos[j][1]+l-anker_pos[1])*temp_spg.mat[2][1]
											+(temp_spg.pos[j][2]+m-anker_pos[2])*temp_spg.mat[2][2];
											
			extern void SPG_add_atom(SPG_tools* spg,int type,double x,double y,double z);
			QB_c2f(&input,c_pos,f_pos);
			
			if(f_pos[2]>max_c&&f_pos[2]<(max_c+displace_s))
			if(SK_SPG_overloop_vec(&spg,temp_spg.type[j],f_pos[0],f_pos[1],f_pos[2],10e-4)==0)
			{
				SPG_add_atom(&spg,temp_spg.type[j],f_pos[0],f_pos[1],f_pos[2]);
				QB_create_atom(&SKA_site,flag_list[temp_spg.type[j]],c_pos[0],c_pos[1],c_pos[2]);
				QB_wrap(&SKA_site,SKA_site.TotalNumber-1);
			}
		}
		free(flag_list);
		QB_dump_vasp_Direct(&SKA_site,argv[3]);
		SKA_dump_status(&SKA_site);
		return 1;
	}
	
	if(strcmp(argv[1],"--adsorb")==0)
	{
		//read input
		QB_tools input;
		QB_init(&input); 
		QB_read_vasp(&input,argv[3]);
			
		QB_tools SKA_site;
		QB_init(&SKA_site); 
		QB_read_vasp(&SKA_site,argv[4]);
			
		int found_flag=0;
		int type=(int)QB_checkdat(argv[6]);
		int target_id=(int)QB_checkdat(argv[7]);
		double target_pos[3]={0,0,0};
		int current_id=1;
		for(int i=0;i<SKA_site.TotalNumber;i++)
		{
			if(SKA_site.atom[i].type==type)
			{
				if(current_id==target_id)
				{
					target_pos[0]=SKA_site.atom[i].x;
					target_pos[1]=SKA_site.atom[i].y;
					target_pos[2]=SKA_site.atom[i].z;
					found_flag=1;
					break;
				}
				else current_id++;
			}
		}
		if(!found_flag)
		{
			printf(">>ERROR: Requested site not found\n");
			exit(0);
		}
		if(strcmp(argv[2],"-ele")==0)
		{
			if(argc<8)
			{
				printf(">>EROOR: Not enough input for add adsorption mode\n");
				printf(">>Please use '--adsorb -ele [file1(Slab)] [file2(Site List)] [file3(Output)] [Type] [ID] [Element]' instead\n");
				exit(0);
			}
			printf(">>Adding adsorption please wait\n");
			int target_type=QB_get_elename(&input,argv[8]);
			QB_create_atom(&input,target_type,target_pos[0],target_pos[1],target_pos[2]);
		}
		
		if(strcmp(argv[2],"-file")==0)
		{
			if(argc<10)
			{
				printf(">>EROOR: Not enough input for adsorption scan mode\n");
				printf(">>Please use '--adsorb -file [file1(Slab)] [file2(Site List)] [file3(Output)] [Type] [ID] [file4(Molecular)] [Element] [Element ID]' instead\n");
				exit(0);
			}
			printf(">>Adding adsorption please wait\n");
			QB_tools adsorb;
			QB_init(&adsorb); 
			QB_read_vasp(&adsorb,argv[8]);
			int found_flag=0;
			int ad_type=QB_get_elename(&adsorb,argv[9]);
			int ad_target_id=(int)QB_checkdat(argv[10]);
			double ad_target_pos[3]={0,0,0};
			int ad_current_id=1;
			char ele_name[256];
			for(int i=0;i<adsorb.TotalNumber;i++)
			{
				if(adsorb.atom[i].type==ad_type)
				{
					if(ad_current_id==ad_target_id)
					{
						ad_target_pos[0]=adsorb.atom[i].x;
						ad_target_pos[1]=adsorb.atom[i].y;
						ad_target_pos[2]=adsorb.atom[i].z;
						found_flag=1;
						break;
					}
					else ad_current_id++;
				}
			}
			//get rotate matrix
			double rotate_mat[3][3]={{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
			for(int i=0;(i<argc-11)&&(i<9);i++)
				rotate_mat[i/3][i%3]=atof(argv[i+11]);
			//get unit matrix
			for(int i=0;i<3;i++)
			{
				double temp=sqrt(rotate_mat[i][0]*rotate_mat[i][0]+rotate_mat[i][1]*rotate_mat[i][1]+rotate_mat[i][2]*rotate_mat[i][2]);
				rotate_mat[i][0]/=temp;
				rotate_mat[i][1]/=temp;
				rotate_mat[i][2]/=temp;
			}
				
			for(int i=0;i<adsorb.TotalNumber;i++)
			{
				QB_load_elename(&adsorb,adsorb.atom[i].type-1,ele_name);
				int target_type=QB_get_elename(&input,ele_name);
				double atom_pos[3];
				for(int j=0;j<3;j++)
				{
					atom_pos[j]=rotate_mat[j][0]*(adsorb.atom[i].x-ad_target_pos[0])+
								rotate_mat[j][1]*(adsorb.atom[i].y-ad_target_pos[1])+
								rotate_mat[j][2]*(adsorb.atom[i].z-ad_target_pos[2]);
				}
				QB_create_atom(&input,target_type,atom_pos[0]+target_pos[0],
												  atom_pos[1]+target_pos[1],
												  atom_pos[2]+target_pos[2]);
			}
		}
		QB_dump_vasp_Direct(&input,argv[5]);
		printf(">>Adding adsorption complete! File dumped to %s\n",argv[5]);
		return 1;
	}
	return 0;
}