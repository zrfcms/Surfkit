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
#include"main.h"
#include"stdio.h"
#include"string.h"
void help()
{
	printf(">>To find conventional cell:\n");
	printf("--conv   [file1(In)] [file2(Out)]\n");
	printf(">>To find 2D primitive cell:\n");
	printf("--prim2D [file1(In)] [file2(Out)]\n");
	printf(">>To redefinition of crystal cell:\n");
	printf("--proj -mat  [file1(In)] [file2(Out)] ah ak al bh bk bl ch ck cl\n");
	printf("--proj -auto [file1(In)] [file2(Out)] ch ck cl\n");
	printf("-–dup [file1(In)] [file2(Out)] ah bk cl\n");
	printf(">>To scan atomic layer and non-polar center mode\n");
	printf("--scan file(INPOS) [min_height] [max_height]\n");
	printf(">>To cleavage of crystal slab:\n");
	printf("--shift  [file1(In)] [file2(Out)] [shift height]\n");
	printf("--cleav  [file1(In)] [file2(Out)] [min_height] [max_height]\n");
	printf("--vacuum [file1(In)] [file2(Out)] [vacuum thickness]\n");
	printf(">>To search twist angles with small strain:\n");
	printf("–-match [file1(Slab)] [file2(Slab)] {[limit] [parameter]}...\n");
	printf(">>To merge slabs with twist angles\n");
	printf("--twist [file1(Slab)] [file2(Slab)] [file3(Output)] [twist angle]\n");
	printf(">>To adjust space between layers\n");
	printf("--space [file(Slab)] [file2(Output)] [Height] {[Max step] [Distance/auto]}\n");		
	printf(">>To search for adsorption site\n");
	printf("--site2D [file1(Slab)] [file2(Site List)] [Distance(2D)] [Max multiplicity]\n");
	printf("-–site3D [file1(Slab)] [file2(Site List)] [Distance(3D)] [Max multiplicity]\n");
	printf(">>To absorb molecule to site\n");
	printf("--adsorb -ele  [file1(Slab)] [file2(Site List)] [file3(Output)] [Type] [ID] [Element]\n");
	printf("--adsorb -file [file1(Slab)] [file2(Site List)] [file3(Output)] [Type] [ID] [file4(Molecular)] [Element] [Element ID]\n");
	printf("--adsorb -file [file1(Slab)] [file2(Site List)] [file3(Output)] [Type] [ID] [file4(Molecular)] [Element] [Element ID] [M11] [M12] [M13] [M21] [M22] [M23] [M31] [M32] [M33]\n");	
}

int main(int argc,char* argv[])
{
	printf("**************************************************\n");
    printf("***   This program is designed at BUAA by      ***\n");
    printf("***       Z. R. Liu and R. F. Zhang            ***\n");
    printf("***   Provided 'as is' without any warranty    ***\n");
    printf("***   Copyright[c] 2022-2023, zrfbuaa group    ***\n");
    printf("**************************************************\n");
	
	if(argc<2)
		help();
	else if(!strcmp(argv[1],"-h"))
		help();
	else 
	{
		int status=0;
		status+=main_slab(argc,argv);
		status+=main_layer(argc,argv);
		status+=main_scan(argc,argv);
		if(!status)
		{
			printf(">>ERROR: Syntax not found\n");
			help();
		}
	}

		
	return 0;
}