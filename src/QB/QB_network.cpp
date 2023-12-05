/*===============================================================
 * This is a part of SPaMD code
 * This program is modified at BUAA by  Z. R. Liu, March. 11 2019
 * Copyright[c] 2017-2019, zrfbuaa group
*================================================================*/
#include"QB.h"
void QB_network_init(QB_tools *QB,double cutoff)
{
    if(cutoff<3)
    	cutoff=3;
    if(QB->network_flag==1)
	QB_network_free(QB);
    QB->xlot=(int)(QB->boundx/cutoff)+2;
    QB->ylot=(int)(QB->boundy/cutoff)+2;
    QB->zlot=(int)(QB->boundz/cutoff)+2;
    int i,j,k;
    QB->tbx=(QB->boundx+2.000*cutoff)/(double)(QB->xlot);
    QB->tby=(QB->boundy+2.000*cutoff)/(double)(QB->ylot);
    QB->tbz=(QB->boundz+2.000*cutoff)/(double)(QB->zlot);
    int x,y,z;
    QB->network_ID=(int****)malloc(QB->xlot*sizeof(int***));
    QB->network_num=(int***)malloc(QB->xlot*sizeof(int**));
    for(i=0;i<QB->xlot;i++)
    {
        QB->network_ID[i]=(int***)malloc(QB->ylot*sizeof(int**));
        QB->network_num[i]=(int**)malloc(QB->ylot*sizeof(int*));
    }
    for(i=0;i<QB->xlot;i++)
        for(j=0;j<QB->ylot;j++)
        {
            QB->network_ID[i][j]=(int**)malloc(QB->zlot*sizeof(int*));
            QB->network_num[i][j]=(int*)malloc(QB->zlot*sizeof(int));
        }
    //know how many atoms there are in each box
    for(i=0;i<QB->xlot;i++)
        for(j=0;j<QB->ylot;j++)
            for(k=0;k<QB->zlot;k++)
            {
                QB->network_num[i][j][k]=0;
            }

    for(i=0;i<QB->TotalNumber;i++)
    {
        x=(QB->atom[i].x-QB->startx+cutoff)/(QB->tbx);
        x=QB_Min(QB_Max(x,0),QB->xlot-1);
        y=(QB->atom[i].y-QB->starty+cutoff)/(QB->tby);
        y=QB_Min(QB_Max(y,0),QB->ylot-1);
        z=(QB->atom[i].z-QB->startz+cutoff)/(QB->tbz);
        z=QB_Min(QB_Max(z,0),QB->zlot-1);
        QB->network_num[x][y][z]++;
    }
    //give space to every pointer
    for(i=0;i<(QB->xlot);i++)
        for(j=0;j<(QB->ylot);j++)
            for(k=0;k<(QB->zlot);k++)
            {
                QB->network_ID[i][j][k]=(int*)malloc(QB_Max(QB->network_num[i][j][k],1)*sizeof(int));
            }
    //reset number of atoms in each box so network_num can show how many atoms is set in that box
    for(i=0;i<(QB->xlot);i++)
        for(j=0;j<(QB->ylot);j++)
            for(k=0;k<(QB->zlot);k++)
                QB->network_num[i][j][k]=0;
                
    for(i=0;i<QB->TotalNumber;i++)
    {
        x=(QB->atom[i].x-QB->startx+cutoff)/(QB->tbx);
        x=QB_Min(QB_Max(x,0),QB->xlot-1);
        y=(QB->atom[i].y-QB->starty+cutoff)/(QB->tby);
        y=QB_Min(QB_Max(y,0),QB->ylot-1);
        z=(QB->atom[i].z-QB->startz+cutoff)/(QB->tbz);
        z=QB_Min(QB_Max(z,0),QB->zlot-1);
        QB->network_ID[x][y][z][QB->network_num[x][y][z]]=i;
        QB->network_num[x][y][z]++;
    }
    QB->network_flag=1;
}
void QB_network_update(QB_tools *QB,double cutoff,int i)
{
    if(QB->network_flag==0)return;
    int x,y,z;
    int x0,y0,z0;
    int temp;
    int j;
    x=(QB->atom[i].x-QB->startx+cutoff)/(QB->tbx);
    x=QB_Min(QB_Max(x,0),QB->xlot-1);
    y=(QB->atom[i].y-QB->starty+cutoff)/(QB->tby);
    y=QB_Min(QB_Max(y,0),QB->ylot-1);
    z=(QB->atom[i].z-QB->startz+cutoff)/(QB->tbz);
    z=QB_Min(QB_Max(z,0),QB->zlot-1);
    

    for(j=0;j<QB->network_num[x][y][z];j++)
    {
        if(i==QB->network_ID[x][y][z][j])
        {
            return;//no need to add
        }
    }//delete it from old list 
	
   
    QB->network_ID[x][y][z]=(int*)realloc(QB->network_ID[x][y][z],(QB->network_num[x][y][z]+1)*sizeof(int));
    QB->network_ID[x][y][z][QB->network_num[x][y][z]]=i;
    QB->network_num[x][y][z]++;
}
//free pointer
void QB_network_free(QB_tools *QB)
{
    int i,j,k;
    for(i=0;i<(QB->xlot);i++)
        for(j=0;j<(QB->ylot);j++)
            for(k=0;k<(QB->zlot);k++)
                free(QB->network_ID[i][j][k]);
    for(i=0;i<QB->xlot;i++)
        for(j=0;j<QB->ylot;j++)
        {
            free(QB->network_ID[i][j]);
            free(QB->network_num[i][j]);
        }
    for(i=0;i<QB->xlot;i++)
    {
        free(QB->network_ID[i]);
        free(QB->network_num[i]);
    }
    QB->network_flag=0;
}