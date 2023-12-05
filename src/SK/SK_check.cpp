#include"./../QB/QB.h"
#include"./../QSPG/QSPG.h"
#include"SurfaceKit.h"

static void SK_layer_add(int *list_n,SK_layer_list **list,int type,double height,char name[1024],char integration[1024])
{
	if((*list_n)==0)
		(*list)=(SK_layer_list*)malloc(sizeof(SK_layer_list));
	else
		(*list)=(SK_layer_list*)realloc(*list,(*list_n+1)*sizeof(SK_layer_list));
	(*list)[*list_n].type=type;
	(*list)[*list_n].height=height;
	strcpy((*list)[*list_n].name,name);
	strcpy((*list)[*list_n].integration,integration);
	(*list_n)++;
}

void SK_layer_sort(int list_n,SK_layer_list *list)
{
	double height;
	int type;
	char name[1024];
	char integration[1024];
	for(int i=0;i<list_n;i++)
	for(int j=0;j<list_n-1;j++)
	{
		if(list[j].height<list[j+1].height)
		{
			type=list[j].type;
			list[j].type=list[j+1].type;
			list[j+1].type=type;
			
			height=list[j].height;
			list[j].height=list[j+1].height;
			list[j+1].height=height;
			
			strcpy(name,list[j].name);
			strcpy(list[j].name,list[j+1].name);
			strcpy(list[j+1].name,name);
			
			strcpy(integration,list[j].integration);
			strcpy(list[j].integration,list[j+1].integration);
			strcpy(list[j+1].integration,integration);
		}
	}
}
//dump layer data to screen
void SK_layer_dump(char output_name[1024],int list_n,SK_layer_list *list)
{
	FILE* output=fopen(output_name,"w");
	printf("%s\t%s\t%s\t%s\n","Layer_Type","Height","Contains","Integration");
	fprintf(output,"%s,%s,%s,%s\n","Layer_Type","Height","Contains","Integration");
	char type_name[2][32]={"Atomic_Layer","Nonpolar_Center"};
	for(int i=0;i<list_n;i++)
	{
		
		if(list[i].type==SK_ATOMIC_LAYER)
		{
			printf("%s\t%lf\t%s\t%s\n",type_name[list[i].type],list[i].height,list[i].name,list[i].integration);	
			fprintf(output,"%s,%lf,%s,%s\n",type_name[list[i].type],list[i].height,list[i].name,list[i].integration);	
		}
		else
		{	
			printf("%s\t%lf\n",type_name[list[i].type],list[i].height);	
			fprintf(output,"%s,%lf\n",type_name[list[i].type],list[i].height);	
		}
	}
	fclose(output);
}

void SK_polar_center(int *list_n,SK_layer_list **list,double min_height,double max_height,SPG_tools* spg,double symprec)
{
	double xproduct[3];
	double xin1[3]={spg->mat[0][0],spg->mat[1][0],spg->mat[2][0]};	
	double xin2[3]={spg->mat[0][1],spg->mat[1][1],spg->mat[2][1]};
	double xin3[3]={spg->mat[0][2],spg->mat[1][2],spg->mat[2][2]};
	
	SK_cross_product(xproduct,xin1,xin2);
	double normal_c=sqrt(xproduct[0]*xproduct[0]
						+xproduct[1]*xproduct[1]
						+xproduct[2]*xproduct[2]);
	double height=(xproduct[0]*xin3[0]
				  +xproduct[1]*xin3[1]
				  +xproduct[2]*xin3[2])/normal_c;
				 
	SpglibDataset *dataset=spg_get_dataset(spg->mat,spg->pos,spg->type,spg->num,symprec);
	for(int i=0;i<dataset->n_operations;i++)
	{
		if(fabs(dataset->rotations[i][0][2])+fabs(dataset->rotations[i][1][2])+fabs(dataset->rotations[i][2][2]+1)<symprec)
		{
			double temp_height=0.5*dataset->translations[i][2];
			while(temp_height<0)temp_height++;
			while(temp_height>=1)temp_height--;
			temp_height=height*temp_height;
			for(int j=2*(int)floor(min_height/height)-2;j<2*(int)floor(max_height/height)+3;j++)
			{
				double out_height=0.5*height*j+temp_height;
				if(out_height>=min_height)
				if(out_height<=max_height)
				{
					//evade same layer
					int flag=0;
					for(int k=0;k<(*list_n);k++)
					{
						if((*list)[k].type==SK_POLAR_CENTER)
						if(fabs((*list)[k].height-out_height)<symprec)
						{
							flag=1;
							break;
						}
					}
					if(!flag)
						SK_layer_add(list_n,list,SK_POLAR_CENTER,out_height,"None","None");
				}
			}	
		}
	}
}

static void SK_atomic_layer(char name[1024],double min_height,double max_height,QB_tools* QB)
{
	strcpy(name,"");
	char ele_name[1024]="";
	for(int j=0;j<QB->TypeNumber;j++)
	{
		int type_num=0;
		QB_load_elename(QB,j,ele_name);
		for(int i=0;i<QB->TotalNumber;i++)
		{
			if(QB->atom[i].z<=max_height)
			if(QB->atom[i].z>min_height)
			if(QB->atom[i].type==j+1)
				type_num++;
		}
		if(type_num)
			sprintf(name,"%s%s%d",name,ele_name,type_num);
	}
}

void SK_atomic_layers(int *list_n,SK_layer_list **list,double min_height,double max_height,QB_tools *QB,double symprec)
{
	QSPG_reshapebox(QB,QB,0,0,1,1);
	int low_layer=(int)floor(min_height/QB->mat[2][2]);
	int high_layer=(int)floor(max_height/QB->mat[2][2])+1;
	double delta=(double)low_layer*QB->mat[2][2];
	double test_min_height=min_height-delta;
	double test_max_height=max_height-delta;
	
	QB_duplicate(QB,1,1,high_layer-low_layer);
	QB_add_exint(QB,"Mark");	
	for(int i=0;i<QB->TotalNumber;i++)
		QB_save_data(QB,i,"Mark",0);
	
	for(int i=0;i<QB->TotalNumber;i++)
	{
		if(QB->atom[i].z<=test_max_height)
		if(QB->atom[i].z>test_min_height)
		{
			int data=(int)QB_get_data(QB,i,"Mark");
			if(!data)
			{
				for(int j=0;j<QB->TotalNumber;j++)
				{
					if(fabs(QB->atom[j].z-QB->atom[i].z)<symprec)
						QB_save_data(QB,j,"Mark",1);
					
				}
				char name[1024];
				SK_atomic_layer(name,QB->atom[i].z-symprec,QB->atom[i].z+symprec,QB);
				char integration[1024];
				SK_atomic_layer(integration,test_min_height,QB->atom[i].z+symprec,QB);
				SK_layer_add(list_n,list,SK_ATOMIC_LAYER,QB->atom[i].z+delta,
					name,integration);
			}
		}
	}
		
}
