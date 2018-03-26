# include "xtool.h"

Molecule::Molecule()
{
	Clear();
}

Molecule::~Molecule()
{
	if(atom) delete [] atom; atom=NULL;
	if(bond) delete [] bond; bond=NULL;

	ring.clear();
	vol_dot.clear();
	sur_dot.clear();
}

Molecule::Molecule(int max_atom_num, int max_bond_num)
{
	Clear();

	atom=new Atom[max_atom_num];
	if(atom==NULL) Memory_Allocation_Error();

	bond=new Bond[max_bond_num];
	if(bond==NULL) Memory_Allocation_Error();
}

Molecule::Molecule(const Molecule &original)
{
	xtool_format=original.xtool_format;
	id=original.id;
	valid=original.valid;
	strcpy(name,original.name);
	weight=original.weight;
	strcpy(formula,original.formula);
	logp=original.logp;
	surface=original.surface;
	nsur=original.nsur;
	psur=original.psur;
	volume=original.volume;
	num_hb_atom=original.num_hb_atom;
	num_rotor=original.num_rotor;

	num_subst=original.num_subst;
	num_feature=original.num_feature;
	num_set=original.num_set;
	strcpy(mol_type,original.mol_type);
	strcpy(charge_type,original.charge_type);

	int i,n;

	num_atom=original.num_atom;

	if(original.atom!=NULL)
	{
	 atom=new Atom[num_atom]; if(atom==NULL) Memory_Allocation_Error();
	 for(i=0;i<num_atom;i++) atom[i]=original.atom[i];
	}

	num_bond=original.num_bond;

	if(original.bond!=NULL)
	{
	 bond=new Bond[num_bond]; if(bond==NULL) Memory_Allocation_Error();
	 for(i=0;i<num_bond;i++) bond[i]=original.bond[i];
	}

	num_ring=original.num_ring;
	ring.clear(); n=original.ring.size(); 
	for(i=0;i<n;i++) ring.push_back(original.ring[i]);

	vol_dot.clear(); n=original.vol_dot.size();
	for(i=0;i<n;i++) vol_dot.push_back(original.vol_dot[i]);

	sur_dot.clear(); n=original.sur_dot.size();
	for(i=0;i<n;i++) sur_dot.push_back(original.sur_dot[i]);
}

Molecule& Molecule::operator = (const Molecule &original)
{
	if(this==&original) return *this;

	xtool_format=original.xtool_format;
	id=original.id;
	valid=original.valid;
	strcpy(name,original.name);
	weight=original.weight;
	strcpy(formula,original.formula);
	logp=original.logp;
	surface=original.surface;
	nsur=original.nsur;
	psur=original.psur;
	volume=original.volume;
	num_hb_atom=original.num_hb_atom;
	num_rotor=original.num_rotor;

	num_subst=original.num_subst;
	num_feature=original.num_feature;
	num_set=original.num_set;
	strcpy(mol_type,original.mol_type);
	strcpy(charge_type,original.charge_type);

	int i,n;

	// you need to delete atom and bond first before re-create them
	// this is slightly different from the copy constructor

	num_atom=original.num_atom; 
	if(atom) delete [] atom; atom=NULL;

	if(original.atom!=NULL)
	{
	 atom=new Atom[num_atom]; if(atom==NULL) Memory_Allocation_Error();
	 for(i=0;i<num_atom;i++) atom[i]=original.atom[i];
	}

	num_bond=original.num_bond; 
	if(bond) delete [] bond; bond=NULL; 

	if(original.bond!=NULL)
	{
	 bond=new Bond[num_bond]; if(bond==NULL) Memory_Allocation_Error();
	 for(i=0;i<num_bond;i++) bond[i]=original.bond[i];
	}

        num_ring=original.num_ring;
        ring.clear(); n=original.ring.size();
        for(i=0;i<n;i++) ring.push_back(original.ring[i]);

       	vol_dot.clear(); n=original.vol_dot.size();
        for(i=0;i<n;i++) vol_dot.push_back(original.vol_dot[i]);

        sur_dot.clear(); n=original.sur_dot.size();
        for(i=0;i<n;i++) sur_dot.push_back(original.sur_dot[i]);

	return *this;
}

void Molecule::Clear()
{
	xtool_format=false;
	id=0;
	valid=0;
        strcpy(name,"");
        strcpy(formula,"");
        weight=0.000;
	logp=0.000;
	surface=0.000;
	nsur=psur=0.000;
	volume=0.000;
	num_hb_atom=0;
	num_rotor=0;

        num_atom=0; atom=NULL;
        num_bond=0; bond=NULL;

	num_subst=1; num_feature=0; num_set=0;
	strcpy(mol_type,"SMALL"); strcpy(charge_type,"USER_CHARGES");

	num_ring=0; ring.clear();

	vol_dot.clear(); sur_dot.clear();

	return;
}

void Molecule::Show_Contents() const
{
	printf("Molecule: %s\t", name);
	printf("%s\t", formula);
	printf("MW = %6.1f\n", weight);

	printf("HB atom = %d\t", num_hb_atom);
	printf("Rotor = %d\n", num_rotor);

	printf("LogP = %6.2f\t", logp);
	printf("Surface = %6.1f\t", surface);
	printf("Volume = %6.1f\n", volume);

	printf("Number of substructs = %d\n", num_subst);
	printf("Number of features = %d\n", num_feature);
	printf("Number of sets = %d\n", num_set);
	printf("Molecule type = %s\n", mol_type);
	printf("Charge type = %s\n", charge_type);

	Show_Atoms();
	Show_Bonds();
	Show_Rings();

	return;
}

void Molecule::Show_Atoms() const
{
	int i;

	printf("\nTotal number of atoms in this molecule = %d\n",num_atom);

        for(i=0;i<num_atom;i++)
                {
                 printf("%-3d ", atom[i].id);
		 printf("%-4s ", atom[i].name);
                 printf("%8.3f ", atom[i].coor[0]);
                 printf("%8.3f ", atom[i].coor[1]);
                 printf("%8.3f ", atom[i].coor[2]);
		 printf("=%1d ", atom[i].valid);
                 printf("%-5s ", atom[i].type);
		 printf("%4.2f ", atom[i].R);
		 printf("%5.2f ", atom[i].q);
		 printf("%6.1f ", atom[i].weight);
		 printf("%-10s ", atom[i].xtype);
                 printf("%-2s ",  atom[i].hb);
		 printf("%5.2f ", atom[i].logp);
		 printf("(%1d) ", atom[i].ring);
		 printf("<%1d> ", atom[i].part);
		 printf("%7.2f ", atom[i].root[0]);
		 printf("%7.2f ", atom[i].root[1]);
		 printf("%7.2f ", atom[i].root[2]);
		 printf("%1d ",  atom[i].num_neib);
                 printf("%-3d ", atom[i].neib[0]);
                 printf("%-3d ", atom[i].neib[1]);
                 printf("%-3d ", atom[i].neib[2]);
                 printf("%-3d ", atom[i].neib[3]);
                 printf("\n");
                }

	return;
}

void Molecule::Show_Bonds() const
{
	int i;

	printf("\nTotal number of bonds in this molecule = %d\n",num_bond);

        for(i=0;i<num_bond;i++)
                {
                 printf("%-3d ", i+1); 
		 printf("%-3d ", bond[i].id);
		 printf("=%1d ", bond[i].valid);
		 printf("%-3d ", bond[i].atom_1);
		 printf("%-3d ", bond[i].atom_2);
		 printf("%-2s ", bond[i].type); 
                 printf("(%1d) ", bond[i].ring);
		 printf("<%1d> ", bond[i].part);
                 printf("%1d ", bond[i].num_neib);
                 printf("%-3d ", bond[i].neib[0]);
                 printf("%-3d ", bond[i].neib[1]);
                 printf("%-3d ", bond[i].neib[2]);
                 printf("%-3d ", bond[i].neib[3]);
                 printf("%-3d ", bond[i].neib[4]);
                 printf("%-3d ", bond[i].neib[5]);
		 printf("\n");
		}

	return;
}

void Molecule::Show_Rings() const
{
	int i,j;

	printf("\nTotal number of rings in this molecule = %d\n", num_ring);

	for(i=0;i<num_ring;i++)
                {
                 printf("%-2d ", i+1);
                 printf("=%1d ", ring[i].valid);
                 printf("%1d | ", ring[i].type);
		 printf("%2d ", ring[i].num_member);

		 for(j=0;j<ring[i].num_member;j++)
			{
			 printf("%5d ",ring[i].atom_id[j]);
			}

		 printf("%8.3f%8.3f%8.3f", 
			ring[i].centroid[0],
			ring[i].centroid[1],
			ring[i].centroid[2]);

		 printf("\n");
                }

	return;
}

// *************************************************************************
// it will return TRUE if read successfully
// it will return FALSE otherwise
// notice: this function will open a file and then read it
// so it is not suitable for reading a multi-Mol2 file
// latest update: 12/04/2003
// *************************************************************************
int Molecule::Read_From_Mol2(char *filename, int position)
{
	FILE *fp;
	int i,tmp,count;
	char buf[256],head[256];

	if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

	count=0; this->xtool_format=false;

	for(;;)
		{
		 if(fgets(buf,256,fp)==NULL) return FALSE;

		 strcpy(head,"");  sscanf(buf,"%s",head);

		 if(!strcmp(head,"###"))
			{
			 if(strstr(buf,"X-TOOL")||
			    strstr(buf,"X-Tool")||
			    strstr(buf,"XTOOL")||
			    strstr(buf,"XTool"))
				{
				 this->xtool_format=true;
				}
			}

		 if(strcmp(head,"@<TRIPOS>MOLECULE")) continue;
		 else count++;

		 if(count==position) break;	// the right one to be read
		 else continue;
		}

	Molecule_Section:

	// get the name of the molecule

       	fgets(buf,256,fp); tmp=strlen(buf);

	for(i=tmp-1;i>=0;i--)
		{
		 if(buf[i]==' '||buf[i]=='\n') buf[i]='\0';
		 else break;
		}

/*
	tmp=strlen(buf);

	for(i=tmp-1;i>=0;i--)
		{
		 if(buf[i]!=' ') continue;
		 else buf[i]='_';
		}
*/

	strcpy(name,buf);

	// now get the number of atoms and the number of bonds

	fgets(buf,256,fp); sscanf(buf,"%d%d",&num_atom,&num_bond);
	fgets(buf,256,fp); sscanf(buf,"%s",mol_type);
	fgets(buf,256,fp); sscanf(buf,"%s",charge_type);

	if(atom) delete [] atom; atom=NULL; 
	atom=new Atom[num_atom];
	if(atom==NULL) Memory_Allocation_Error();

	if(bond) delete [] bond; bond=NULL; 
	bond=new Bond[num_bond];
	if(bond==NULL) Memory_Allocation_Error();

	for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) return FALSE;
                 else {strcpy(head,""); sscanf(buf,"%s",head);}

                 if(!strcmp(head,"@<TRIPOS>MOLECULE")) return FALSE;
                 else if(!strcmp(head,"@<TRIPOS>ATOM")) goto Atom_Section;
                 else if(!strcmp(head,"@<TRIPOS>BOND")) return FALSE;
                 else continue;
                }

	Atom_Section:

	for(i=0;i<num_atom;i++)
		{
		 if(fgets(buf,256,fp)==NULL) return FALSE;
                 else {strcpy(head,""); sscanf(buf,"%s",head);}

                 if(!strcmp(head,"@<TRIPOS>MOLECULE")) return FALSE;
                 else if(!strcmp(head,"@<TRIPOS>ATOM")) return FALSE;
                 else if(!strcmp(head,"@<TRIPOS>BOND")) return FALSE;
                 else 
			{
		 	 sscanf(buf,"%d%s%f%f%f%s%*d%*s%f", 
				&atom[i].id, 
				 atom[i].name,
				&atom[i].coor[0], 
				&atom[i].coor[1], 
				&atom[i].coor[2], 
				 atom[i].type,
				&atom[i].q);
		 	 atom[i].valid=1; atom[i].origin=1;
			}
		}

	for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) return FALSE;
                 else {strcpy(head,""); sscanf(buf,"%s",head);}

                 if(!strcmp(head,"@<TRIPOS>MOLECULE")) return FALSE;
                 else if(!strcmp(head,"@<TRIPOS>ATOM")) return FALSE;
                 else if(!strcmp(head,"@<TRIPOS>BOND")) goto Bond_Section;
                 else continue;
                }

	Bond_Section:

	for(i=0;i<num_bond;i++)
		{
		 if(fgets(buf,256,fp)==NULL) return FALSE;
                 else {strcpy(head,""); sscanf(buf,"%s",head);}

		 if(!strcmp(head,"@<TRIPOS>MOLECULE")) return FALSE;
                 else if(!strcmp(head,"@<TRIPOS>ATOM")) return FALSE;
                 else if(!strcmp(head,"@<TRIPOS>BOND")) return FALSE;
                 else
			{
		 	 sscanf(buf,"%hd%hd%hd%s",
			 	&bond[i].id, 
				&bond[i].atom_1,
				&bond[i].atom_2, 
				 bond[i].type);
		 	 bond[i].valid=1;
			}
		}

	fclose(fp);

	this->valid=1; return TRUE;
}

// *************************************************************************
// notice: this function will read a file from the current position
// so it is suitable for reading multi-Mol2
// it will return TRUE if read successfully
// it will return FALSE otherwise
// *************************************************************************
int Molecule::Read_From_Mol2(FILE *fp, int position)
{
	int i,tmp,count;
	char buf[256],head[256];

	count=0; this->xtool_format=false;

	for(;;)
		{
		 if(fgets(buf,256,fp)==NULL) return FALSE;

		 strcpy(head,"");  sscanf(buf,"%s",head);

		 if(!strcmp(head,"###"))
			{
			 if(strstr(buf,"X-TOOL")||
			    strstr(buf,"X-Tool")||
			    strstr(buf,"XTOOL")||
			    strstr(buf,"XTool"))
				{
				 this->xtool_format=true;
				}
			}

		 if(strcmp(head,"@<TRIPOS>MOLECULE")) continue;
		 else count++;

		 if(count==position) break;	// the right one to be read
		 else continue;
		}

	Molecule_Section:

	// get the name of the molecule

       	fgets(buf,256,fp); tmp=strlen(buf);

	for(i=tmp-1;i>=0;i--)
		{
		 if(buf[i]==' '||buf[i]=='\n') buf[i]='\0';
		 else break;
		}

/*
	tmp=strlen(buf);

	for(i=tmp-1;i>=0;i--)
		{
		 if(buf[i]!=' ') continue;
		 else buf[i]='_';
		}
*/

	strcpy(name,buf);

	// now get the number of atoms and the number of bonds

	fgets(buf,256,fp); sscanf(buf,"%d%d",&num_atom,&num_bond);
	fgets(buf,256,fp); sscanf(buf,"%s",mol_type);
	fgets(buf,256,fp); sscanf(buf,"%s",charge_type);

	if(atom) delete [] atom; atom=NULL; 
	atom=new Atom[num_atom];
	if(atom==NULL) Memory_Allocation_Error();

	if(bond) delete [] bond; bond=NULL; 
	bond=new Bond[num_bond];
	if(bond==NULL) Memory_Allocation_Error();

	for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) return FALSE;
                 else {strcpy(head,""); sscanf(buf,"%s",head);}

                 if(!strcmp(head,"@<TRIPOS>MOLECULE")) goto Error_Section;
                 else if(!strcmp(head,"@<TRIPOS>ATOM")) goto Atom_Section;
                 else if(!strcmp(head,"@<TRIPOS>BOND")) goto Error_Section;
                 else continue;
                }

	Atom_Section:

	for(i=0;i<num_atom;i++)
		{
		 if(fgets(buf,256,fp)==NULL) return FALSE;
                 else {strcpy(head,""); sscanf(buf,"%s",head);}

                 if(!strcmp(head,"@<TRIPOS>MOLECULE")) goto Error_Section;
                 else if(!strcmp(head,"@<TRIPOS>ATOM")) goto Error_Section;
                 else if(!strcmp(head,"@<TRIPOS>BOND")) goto Error_Section;
                 else 
			{
		 	 sscanf(buf,"%d%s%f%f%f%s%*d%*s%f", 
				&atom[i].id, 
				 atom[i].name,
				&atom[i].coor[0], 
				&atom[i].coor[1], 
				&atom[i].coor[2], 
				 atom[i].type,
				&atom[i].q);
		 	 atom[i].valid=1; atom[i].origin=1;
			}
		}

	for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) return FALSE;
                 else {strcpy(head,""); sscanf(buf,"%s",head);}

                 if(!strcmp(head,"@<TRIPOS>MOLECULE")) goto Error_Section;
                 else if(!strcmp(head,"@<TRIPOS>ATOM")) goto Error_Section;
                 else if(!strcmp(head,"@<TRIPOS>BOND")) goto Bond_Section;
                 else continue;
                }

	Bond_Section:

	for(i=0;i<num_bond;i++)
		{
		 if(fgets(buf,256,fp)==NULL) return FALSE;
                 else {strcpy(head,""); sscanf(buf,"%s",head);}

		 if(!strcmp(head,"@<TRIPOS>MOLECULE")) goto Error_Section;
                 else if(!strcmp(head,"@<TRIPOS>ATOM")) goto Error_Section;
                 else if(!strcmp(head,"@<TRIPOS>BOND")) goto Error_Section;
                 else
			{
		 	 sscanf(buf,"%hd%hd%hd%s",
			 	&bond[i].id, 
				&bond[i].atom_1,
				&bond[i].atom_2, 
				 bond[i].type);
		 	 bond[i].valid=1;
			}
		}

	Right_Section:

	this->valid=1; return TRUE;

	Error_Section:

	tmp=strlen(buf); fseek(fp,-tmp,SEEK_CUR);

	this->valid=0; return FALSE;
}

int Molecule::Write_Out_Mol2(char *filename, char *mode) const
// suitable for single molecule output
{
	FILE *fp;
	int i;

	if((fp=fopen(filename,mode))==NULL) Open_File_Error(filename);

	fprintf(fp,"### \n"); 
	fprintf(fp,"### Created by X-TOOL on %s\n", Get_Time()); 
	fprintf(fp,"### \n\n");

	fprintf(fp,"@<TRIPOS>MOLECULE\n");
	fprintf(fp,"%s\n",name);
	fprintf(fp,"%5d %5d %5d %5d %5d\n",
		num_atom,num_bond,num_subst,num_feature,num_set);
	fprintf(fp,"%s\n",mol_type);
	fprintf(fp,"%s\n",charge_type);
	fprintf(fp,"\n\n");

	fprintf(fp,"@<TRIPOS>ATOM\n");

	for(i=0;i<num_atom;i++)
		{
		 if(!strcmp(atom[i].type,"H.spc")) strcpy(atom[i].type,"H");

		 fprintf(fp,"%7d %-8s %9.4f %9.4f %9.4f %-5s ",
			 atom[i].id, 
			 atom[i].name, 
			 atom[i].coor[0],
			 atom[i].coor[1], 
			 atom[i].coor[2], 
			 atom[i].type);
		 fprintf(fp,"%5d %3s%-6d %7.3f\n",
			 atom[i].part,
			 "SUB",
			 atom[i].part,
			 atom[i].q); 
		}
	
	fprintf(fp,"@<TRIPOS>BOND\n");

	for(i=0;i<num_bond;i++)
		{
		 fprintf(fp,"%6d %4d %4d %-3s\n", 
			 bond[i].id, 
			 bond[i].atom_1, 
			 bond[i].atom_2, 
			 bond[i].type);
		}

	int *root;

	root=new int[num_subst];
	if(root==NULL) Memory_Allocation_Error();

	for(i=0;i<num_subst;i++) root[i]=0;

	for(i=0;i<num_atom;i++)
		{
		 if(root[atom[i].part-1]!=0) continue;
		 else root[atom[i].part-1]=atom[i].id;
		}

	fprintf(fp,"@<TRIPOS>SUBSTRUCTURE\n");

	for(i=0;i<num_subst;i++)
                {
                 fprintf(fp,"%6d %3s%-3d %6d %-s\n",
			 i+1,"SUB",i+1,root[i],"GROUP");
                }

	fprintf(fp,"\n");

	if(root) delete [] root;

	fclose(fp); return TRUE;
}

int Molecule::Write_Out_Mol2(FILE *fp) const
// suitable for multiple molecule output
{
	int i;

	fprintf(fp,"### \n"); 
	fprintf(fp,"### Created by X-TOOL on %s\n", Get_Time()); 
	fprintf(fp,"### \n\n");

	fprintf(fp,"@<TRIPOS>MOLECULE\n");
	fprintf(fp,"%s\n",name);
	fprintf(fp,"%5d %5d %5d %5d %5d\n",
		num_atom,num_bond,num_subst,num_feature,num_set);
	fprintf(fp,"%s\n",mol_type);
	fprintf(fp,"%s\n",charge_type);
	fprintf(fp,"\n\n");

	fprintf(fp,"@<TRIPOS>ATOM\n");

	for(i=0;i<num_atom;i++)
		{
		 if(!strcmp(atom[i].type,"H.spc")) strcpy(atom[i].type,"H");

		 fprintf(fp,"%7d %-8s %9.4f %9.4f %9.4f %-5s ",
			 atom[i].id, 
			 atom[i].name, 
			 atom[i].coor[0],
			 atom[i].coor[1], 
			 atom[i].coor[2], 
			 atom[i].type);
		 fprintf(fp,"%5d %3s%-6d %7.3f\n",
			 atom[i].part,
			 "SUB",
			 atom[i].part,
			 atom[i].q); 
		}
	
	fprintf(fp,"@<TRIPOS>BOND\n");

	for(i=0;i<num_bond;i++)
		{
		 fprintf(fp,"%6d %4d %4d %-3s\n", 
			 bond[i].id, 
			 bond[i].atom_1, 
			 bond[i].atom_2, 
			 bond[i].type);
		}

	int *root;

	root=new int[num_subst];
	if(root==NULL) Memory_Allocation_Error();

	for(i=0;i<num_subst;i++) root[i]=0;

	for(i=0;i<num_atom;i++)
		{
		 if(root[atom[i].part-1]!=0) continue;
		 else root[atom[i].part-1]=atom[i].id;
		}

	fprintf(fp,"@<TRIPOS>SUBSTRUCTURE\n");

	for(i=0;i<num_subst;i++)
                {
                 fprintf(fp,"%6d %3s%-3d %6d %-s\n",
			 i+1,"SUB",i+1,root[i],"GROUP");
                }

	fprintf(fp,"\n");

	if(root) delete [] root;

	fflush(fp); return TRUE;
}

int Molecule::Value_Atom()
{
	extern ForceField *ff;
	int i;
	bool mark,success;
	char xtype[20];

	success=true;

	this->Detect_Connections();

	// check if this molecule is valid or not

	// if(this->xtool_format==false) success=this->Check_Atom_Type();
	success=this->Check_Atom_Type();

	// now detect the ring systems in this molecule

	this->Detect_Rings(); 
	this->Detect_Aromatic_Rings();

	// now assign XTOOL atom type and corresponding parameters

	for(i=0;i<num_atom;i++)
		{
	 	 mark=this->Get_XTOOL_Type(atom[i].id,xtype);

	 	 if(mark==true) 
			{
			 strcpy(atom[i].xtype, xtype); 
			 atom[i].valid=1;
			}
	 	 else 
			{
			 printf("Warning: unknown atom type: %d %s\n", 
				atom[i].id, atom[i].type);
			 strcpy(atom[i].xtype, xtype);
			 atom[i].valid=0; success=false; 
			}
		}

	for(i=0;i<num_atom;i++) 
		{
	 	 if(atom[i].valid<=0) continue;

	 	 mark=ff->Assign_Atom_Parameters(atom[i]);

	 	 if(mark==false) {atom[i].valid=0; success=false;} 
	 	 else continue;
		}

	// get other necessary information

	this->Get_Formula(this->formula);
	this->weight=this->Get_Weight();

	// now return value

	if(success==true) return TRUE;
	else return FALSE;
}

void Molecule::Assign_Apparent_Charge()
{
	int i;

	strcpy(this->charge_type,"FORMAL_CHARGES");

	for(i=0;i<num_atom;i++)
		{
		 if(!strcmp(atom[i].xtype,"O.co2")) atom[i].q=-0.500;
		 else if(!strcmp(atom[i].xtype,"N.4")) atom[i].q=1.000;
		 else atom[i].q=0.000;
		}
	
	for(i=0;i<num_atom;i++)
		{
		 if(strcmp(atom[i].xtype,"N.pl3.h")) continue;
		 else if(atom[i].num_nonh!=1) continue;
		 else if(strcmp(atom[atom[i].neib[0]-1].type,"C.cat")) continue;
		 else atom[i].q=0.500;
		}

	return;
}

void Molecule::Assign_MMFF94_Charge()  // to be done 
{
	return;
}

// *****************************************************************************
// Detect and automatically correct some common errors in atom type 
// and bond type assignment for a Mol2 file 
// Latest update 12/01/2003 
// *****************************************************************************
bool Molecule::Check_Atom_Type()
{
	int i,j,id,count,temp;
	bool mark;
	Group group;

	// first of all, check carbon atoms
	// this is needed by the following hydrogen atom check

	mark=true;

	for(i=0;i<num_atom;i++)
		{
	 	 if(!strcmp(atom[i].type,"C.2"))
			{
		 	 group=Find_A_Group(atom[i].id);

		 	 if((group.num_db==2)&&
                            (group.num_nonh==2)&&
                    	    (group.num_neib==2))  // =C=
				{
			 	 strcpy(atom[i].type,"C.1");
				}
			 else if(group.num_db<1)
				{
			         mark=false; break;
				}
			}
		 else if(!strcmp(atom[i].type,"C.1"))
			{
			 group=Find_A_Group(atom[i].id);

			 if(group.num_db<1&&group.num_tb<1)
				{
				 mark=false; break;
				}
			}
		 else continue;
		}

	 if(mark==false)
		{
		 puts("Warning: some carbon atoms have wrong types.");
		 return false;
		}

	// first check if hydrogen atoms have been added to this molecule
	// this is done by checking all of the carbon atoms

	count=0; temp=0;

	for(i=0;i<num_atom;i++)
		{
		 if(!strcmp(atom[i].type,"H")) 
			{
			 if(atom[atom[i].neib[0]-1].type[0]!='C') continue;
			 else count++;
			}

		 if(!strcmp(atom[i].type,"C.3"))
			{
			 temp+=(4-atom[i].num_nonh);
			}
		 else if(!strcmp(atom[i].type,"C.2"))
			{
			 temp+=(3-atom[i].num_nonh);
			}
		 else if(!strcmp(atom[i].type,"C.ar"))
			{
			 temp+=(3-atom[i].num_nonh);
			}
		 else if(!strcmp(atom[i].type,"C.cat"))
			{
			 temp+=(3-atom[i].num_nonh);
			}
		 else if(!strcmp(atom[i].type,"C.1"))
			{
			 temp+=(2-atom[i].num_nonh);
			}
		 else continue;
		}

	if(count!=temp)
		{
	 	 puts("Warning: hydrogen atoms may not be correctly added"); 
		 // printf("theory = %d  reality = %d\n", temp, count); 
		 return false;
		}

	// now check atom by atom

	for(i=0;i<num_atom;i++)
        {
         group=Find_A_Group(atom[i].id); 

	 if(group.center.type[0]=='O')
                {
                 if(group.num_neib==1&&group.num_nonh==1)   // =O 
                        {
			 if(!strcmp(group.center.type,"O.co2"))
				{
				 strcpy(group.bond[0].type,"2");
				}
                         else if(!strcmp(group.center.type,"O.2"))
				{
				 strcpy(group.bond[0].type,"2");
				}
			 else
                                {
                                 strcpy(group.center.type,"O.2");
				 strcpy(group.bond[0].type,"2");
                                }
                        }
		 else if(group.num_neib>=2)
			{
			 strcpy(group.center.type,"O.3");
			 for(j=0;j<group.num_nonh;j++)
				{
				 strcpy(group.bond[j].type,"1");
				}
			}
		}

	 if(group.center.type[0]=='P') strcpy(group.center.type,"P.3");

	 if(group.center.type[0]=='S')
		{
		 if((group.db_type==3)&&(group.num_neib==3)) // -SO-
			{
			 strcpy(group.center.type,"S.o");
			}
		 else if((group.db_type==3)&&(group.num_neib==4)) // -SO2-
			{
			 strcpy(group.center.type,"S.o2");
			}
		 else if(group.num_nonh>=4)	// e.g. -SF5
			{
                         strcpy(group.center.type,"S.3");
			 for(j=0;j<group.num_nonh;j++)
				{
				 strcpy(group.bond[j].type,"1");
				}
			}
		 else if(group.num_nonh==1&&group.num_neib==1)   // =S
			{
                         strcpy(group.center.type,"S.2");
			 strcpy(group.bond[0].type,"2");
			}
		}

	 if(!strcmp(group.center.type,"N.2"))
		{
                 if((group.num_db==2)&&(group.num_nonh==2)&&
		    (group.num_neib==2))  // =N=
                        {
                         strcpy(group.center.type,"N.1");
                        }
		}
	 else if(group.center.type[0]=='N')
		{
		 if(group.db_type==3&&group.num_nonh>=2)  // -NO, -NO2 
			{
			 strcpy(group.center.type,"N.2");
			}
		 else if(!strcmp(group.center.type,"N.3")||
 			 !strcmp(group.center.type,"N.4"))
                        {
			 if(group.amide>0) // -NH-SO-, -NH-PO-, or -NH-CO-
				{
				 strcpy(group.center.type,"N.am");
				 strcpy(group.bond[group.amide-1].type,"am");
				}
                         else if(group.num_pi==group.num_nonh)
                                {
                                 strcpy(group.center.type,"N.pl3");
                                }
                        }
		 else if(!strcmp(group.center.type,"N.pl3"))
			{
			 if(group.amide>0) // -NH-SO-, -NH-PO-, or -NH-CO-
                                {
                                 strcpy(group.center.type,"N.am");
				 strcpy(group.bond[group.amide-1].type,"am");
                                }
			}
		}

	 // now correct the atom type and bond type for this atom

	 strcpy(atom[i].type, group.center.type);

	 for(j=0;j<group.num_nonh;j++)
		{
		 id=group.bond[j].id;
		 strcpy(bond[id-1].type,group.bond[j].type);
		}
	}

	return true;
}

int Molecule::Get_Atom_Hybridizing_Type(char *type) const
{
	int mark; 	// sp=1, sp2=2; sp3=3; none=4;

        if(!strcmp(type,"C.3")) mark=3;
        else if(!strcmp(type,"C.2"))   mark=2;
        else if(!strcmp(type,"C.1"))   mark=1; 
        else if(!strcmp(type,"C.cat")) mark=2;
        else if(!strcmp(type,"C.ar"))  mark=2;
        else if(!strcmp(type,"H"))     mark=4;
        else if(!strcmp(type,"H.spc")) mark=4;
        else if(!strcmp(type,"N.4"))   mark=3;
        else if(!strcmp(type,"N.3"))   mark=3;
        else if(!strcmp(type,"N.2"))   mark=2;
        else if(!strcmp(type,"N.1"))   mark=1; 
        else if(!strcmp(type,"N.ar"))  mark=2;
        else if(!strcmp(type,"N.pl3")) mark=2;
        else if(!strcmp(type,"N.am"))  mark=2;
        else if(!strcmp(type,"O.3"))   mark=3;
        else if(!strcmp(type,"O.2"))   mark=2;
        else if(!strcmp(type,"O.co2")) mark=2;
        else if(!strcmp(type,"P.3"))   mark=3;
        else if(!strcmp(type,"S.3"))   mark=3;
        else if(!strcmp(type,"S.2"))   mark=2;
        else if(!strcmp(type,"S.o"))   mark=3;
        else if(!strcmp(type,"S.o2"))  mark=3;
        else if(!strcmp(type,"F"))     mark=4;
        else if(!strcmp(type,"Cl"))    mark=4;
        else if(!strcmp(type,"Br"))    mark=4;
        else if(!strcmp(type,"I"))     mark=4;
	else if(!strcmp(type,"Si"))    mark=3;
        else mark=4;
	
	return mark;
}

// *************************************************************************
// function 1: detect all the fragments in the molecule
// function 2: set up the connection tables for atoms and bonds
// latest update: 11/25/2003
// *************************************************************************
void Molecule::Detect_Connections()
{
	int i,j,k,mark,tmp1,tmp2,old_part,new_part,id1,id2;

	// task 1: detect all the fragments and determine num_subst 

	// first, detect how many fragments are in the molecule
	// initialize the molecule as the assembly of isolated atoms

	for(i=0;i<num_atom;i++) atom[i].part=atom[i].id;

	//  then merge the atoms into fragments 

	for(i=0;i<num_bond;i++)
		{
		 if(bond[i].valid<=0) continue;

		 tmp1=bond[i].atom_1; tmp2=bond[i].atom_2;
		 if(atom[tmp1-1].part>=atom[tmp2-1].part) 
			{
			 old_part=atom[tmp1-1].part; 
			 new_part=atom[tmp2-1].part;
			}
		 else  
			{
			 old_part=atom[tmp2-1].part; 
			 new_part=atom[tmp1-1].part;
			}

		 for(j=0;j<num_atom;j++)
			{
			 if(atom[j].part!=old_part) continue;
			 else atom[j].part=new_part;
			} 
		}

	// then count the number of fragments 
	// also re-arrange the id of fragments to make them continuous

	int *part_list;

	part_list=new int[num_atom];
	if(part_list==NULL) Memory_Allocation_Error();

	for(i=0;i<num_atom;i++) part_list[i]=0;

	num_subst=0;

	for(i=0;i<num_atom;i++)
		{
		 mark=FALSE;

		 for(j=0;j<num_subst;j++)
			{
			 if(atom[i].part!=part_list[j]) continue;
			 else {mark=TRUE; break;}
			}

		 if(mark==TRUE) {atom[i].part=j+1;}	// an existing fragment 
		 else					// a new fragment 
			{
			 part_list[num_subst]=atom[i].part;
			 num_subst++;
			 atom[i].part=num_subst;	// 1,2,3 ...
			}
		}

	if(part_list) delete [] part_list;

	// also define all the bonds

	for(i=0;i<num_bond;i++) 
		{
		 if(bond[i].valid<=0) bond[i].part=0;
		 else bond[i].part=atom[bond[i].atom_1-1].part;
		}

	// task 2: set up the connection tables for atoms and bonds

	// assign atomic weight to atoms, which is necessary to rank 
	// the neighboring atoms for a central atom

        for(i=0;i<num_atom;i++)
                {
		 if(!strcmp(atom[i].type,"F")) atom[i].weight=19.00;
                 else if(!strcmp(atom[i].type,"Cl")) atom[i].weight=35.45;
                 else if(!strcmp(atom[i].type,"Br")) atom[i].weight=79.90;
		 else if(!strcmp(atom[i].type,"I")) atom[i].weight=126.90;
		 else if(!strcmp(atom[i].type,"Si")) atom[i].weight=28.09;
		 else if(atom[i].type[0]=='C') atom[i].weight=12.01;
		 else if(atom[i].type[0]=='H') atom[i].weight=1.00;
		 else if(atom[i].type[0]=='N') atom[i].weight=14.01;
		 else if(atom[i].type[0]=='O') atom[i].weight=16.00;
		 else if(atom[i].type[0]=='P') atom[i].weight=30.97;
		 else if(atom[i].type[0]=='S') atom[i].weight=32.06;
		 else atom[i].weight=0.00;
                }

	// now detect the neighboring atoms for each atom

	for(i=0;i<num_atom;i++) atom[i].num_neib=atom[i].num_nonh=0;

	for(i=0;i<num_bond;i++)
		{
		 if(bond[i].valid<=0) continue;

		 id1=bond[i].atom_1; id2=bond[i].atom_2;

		 tmp1=atom[id1-1].num_neib; 

		 if(tmp1<MAX_ATOM_NEIB) 
			{
		 	 atom[id1-1].neib[tmp1]=id2; 
			 atom[id1-1].bond[tmp1]=i+1; 
		 	 atom[id1-1].num_neib++;
			}

		 tmp2=atom[id2-1].num_neib; 

		 if(tmp2<MAX_ATOM_NEIB)
			{
		 	 atom[id2-1].neib[tmp2]=id1; 
			 atom[id2-1].bond[tmp2]=i+1;
		 	 atom[id2-1].num_neib++;
			}
		}

	// now arrange the neighboring atoms in a decreasing order
	// according to atomic weights; bonds are also re-arranged

	for(i=0;i<num_atom;i++)
		{
		 if(atom[i].valid==0) continue;
		 else if(atom[i].num_neib<=1) continue;

		 for(j=0;j<atom[i].num_neib-1;j++)
		 for(k=j+1;k<atom[i].num_neib;k++)
			{
			 tmp1=atom[i].neib[j]; tmp2=atom[i].neib[k];
			 if(atom[tmp1-1].weight>=atom[tmp2-1].weight) continue;
			 else
				{
				 mark=atom[i].neib[j];
				 atom[i].neib[j]=atom[i].neib[k];
				 atom[i].neib[k]=mark;

				 mark=atom[i].bond[j];
                                 atom[i].bond[j]=atom[i].bond[k];
                                 atom[i].bond[k]=mark;
				}
			}
		}

	for(i=0;i<num_atom;i++)
		{
		 for(j=0;j<atom[i].num_neib;j++)
			{
			 if(atom[atom[i].neib[j]-1].type[0]=='H') continue;
		 	 else atom[i].num_nonh++;
			}
		}

	// now detect the neighboring bonds

	for(i=0;i<num_bond;i++) bond[i].num_neib=0;

	for(i=0;i<num_bond-1;i++)
	for(j=i+1;j<num_bond;j++)
		{
		 if(bond[i].valid<=0||bond[j].valid<=0) continue;

		 if(Two_Bonds_Connection_Check(bond[i],bond[j]))
			{
			 tmp1=bond[i].num_neib;
			 if(tmp1<MAX_BOND_NEIB)
				{
			 	 bond[i].neib[tmp1]=bond[j].id; 
				 bond[i].num_neib++;
				}

			 tmp2=bond[j].num_neib;
			 if(tmp2<MAX_BOND_NEIB)
				{
			 	 bond[j].neib[tmp2]=bond[i].id; 
				 bond[j].num_neib++;
				}
			}

		 else continue;
		}

	return;
}

void Molecule::Translate(float move[3])
{
	int i;

	for(i=0;i<num_atom;i++)
		{
		 atom[i].coor[0]+=move[0];
		 atom[i].coor[1]+=move[1];
		 atom[i].coor[2]+=move[2];
		}

	return;
}

void Molecule::Rotate(float angle,float axis[3],float origin[3])
{
	int i,j,k;
	float old_vector[3],new_vector[3];
	float mat[3][3];	
	float a,b;

	Normal_Vector(axis);

	a=sin((double)(angle*(PI)/180.0));
        b=cos((double)(angle*(PI)/180.0));

	mat[0][0]=axis[0]*axis[0]+(1-axis[0]*axis[0])*b;
        mat[1][0]=axis[0]*axis[1]*(1-b)-axis[2]*a;
        mat[2][0]=axis[0]*axis[2]*(1-b)+axis[1]*a;
        mat[0][1]=axis[0]*axis[1]*(1-b)+axis[2]*a;
        mat[1][1]=axis[1]*axis[1]+(1-axis[1]*axis[1])*b;
        mat[2][1]=axis[1]*axis[2]*(1-b)-axis[0]*a;
        mat[0][2]=axis[0]*axis[2]*(1-b)-axis[1]*a;
        mat[1][2]=axis[1]*axis[2]*(1-b)+axis[0]*a;
        mat[2][2]=axis[2]*axis[2]+(1-axis[2]*axis[2])*b;

	for(i=0;i<num_atom;i++)
		{
		 for(j=0;j<3;j++) 
			{
			 old_vector[j]=atom[i].coor[j]-origin[j];
			 new_vector[j]=0.000;
			}

		 for(j=0;j<3;j++)
			{
			 for(k=0;k<3;k++) 
				new_vector[j]+=(old_vector[k]*mat[k][j]);
			}

		 for(j=0;j<3;j++)
			{
			 atom[i].coor[j]=new_vector[j]+origin[j];
			}
		}

	return;
}

int Molecule::Get_XTOOL_Type(int atom_id, char *type) const
{
	Group group;
	Atom atom;

	atom=this->atom[atom_id-1]; group=Find_A_Group(atom.id);

	strcpy(type,"Un");

	if(!strcmp(atom.type,"H")||!strcmp(atom.type,"H.spc"))
		{
	 	 if(group.neib[0].type[0]=='O') strcpy(type,"H.hb");
	 	 else if(group.neib[0].type[0]=='N') strcpy(type,"H.hb");
	 	 else strcpy(type,"H"); 
		}

	if(!strcmp(atom.type,"C.3"))
		{
		 if(group.num_hetero==0) strcpy(type,"C.3");
		 else if(group.num_hetero>0) strcpy(type,"C.3.x");
		 else strcpy(type,"C.3.un");
		}

	if(!strcmp(atom.type,"C.2")&&(atom.ring!=2))
		{
		 if(group.num_hetero==0) strcpy(type,"C.2");
		 else if(group.num_hetero>0) strcpy(type,"C.2.x");
		 else strcpy(type,"C.2.un");
		}

	if(!strcmp(atom.type,"C.ar")||(!strcmp(atom.type,"C.2")&&atom.ring==2))
		{
		 if(group.num_hetero==0) strcpy(type,"C.ar");
		 else if(group.num_hetero>0) strcpy(type,"C.ar.x");
		 else strcpy(type,"C.ar.un");
		}

	if(!strcmp(atom.type,"C.1"))
		{
		 if(group.num_hetero==0) strcpy(type,"C.1");
		 else if(group.num_hetero>0) strcpy(type,"C.1.x");
		 else strcpy(type,"C.1.un");
		}

	if(!strcmp(atom.type,"C.cat")) strcpy(type,"C.cat");

	if(!strcmp(atom.type,"N.4")||!strcmp(atom.type,"N.3"))
		{
		 if(group.num_nonh<=2) strcpy(type,"N.4");
		 else if(group.num_nonh==3) strcpy(type,"N.3");
		 else strcpy(type,"N.3.un");
		}

	if((!strcmp(atom.type,"N.am")&&atom.ring!=2)||
	   (!strcmp(atom.type,"N.pl3")&&atom.ring!=2))
		{
		 if(group.num_nonh==1) strcpy(type,"N.pl3.h");
		 else if(group.num_nonh==2) strcpy(type,"N.pl3.h");
		 else if(group.num_nonh==3) strcpy(type,"N.pl3");
		 else strcpy(type,"N.pl3.un");
		}

	if(!strcmp(atom.type,"N.2")&&atom.ring!=2)
		{
		 if(group.num_nonh==1) strcpy(type,"N.2.h");
		 else if(group.num_nonh==2) strcpy(type,"N.2");
		 else strcpy(type,"N.2.un");
		}

	if(!strcmp(atom.type,"N.ar")||
	  (!strcmp(atom.type,"N.2")&&atom.ring==2)||
	  (!strcmp(atom.type,"N.pl3")&&atom.ring==2)||
	  (!strcmp(atom.type,"N.am")&&atom.ring==2))
		{
		 if(group.num_h==1) strcpy(type,"N.ar.h");
		 else if(group.num_h==0) strcpy(type,"N.ar");
		 else strcpy(type,"N.ar.un");
		}

	if(!strcmp(atom.type,"N.1")) 
		{
		 if(group.num_nonh==1) strcpy(type,"N.1");
		 else strcpy(type,"N.1.un");
		}

	if(!strcmp(atom.type,"O.3"))
		{
		 if(group.num_nonh==1) strcpy(type,"O.3.h");
		 else if(group.num_nonh==2) strcpy(type,"O.3");
		 else strcpy(type,"O.3.un");
		}

	if(!strcmp(atom.type,"O.2")) strcpy(type,"O.2");

	if(!strcmp(atom.type,"O.co2")) strcpy(type,"O.co2");

	if(!strcmp(atom.type,"S.3"))
		{
		 if(group.num_nonh==1) strcpy(type,"S.3.h");
		 else if(group.num_nonh==2) strcpy(type,"S.3");
		 else strcpy(type,"S.3.un");
		}

	if(!strcmp(atom.type,"S.2")) strcpy(type,"S.2");

	if(!strcmp(atom.type,"S.o")) strcpy(type,"S.o");

	if(!strcmp(atom.type,"S.o2")) strcpy(type,"S.o");

	if(!strcmp(atom.type,"P.3")) strcpy(type,"P.3");

	if(!strcmp(atom.type,"F")) strcpy(type,"F");

	if(!strcmp(atom.type,"Cl")) strcpy(type,"Cl");

	if(!strcmp(atom.type,"Br")) strcpy(type,"Br");

	if(!strcmp(atom.type,"I")) strcpy(type,"I");

	if(!strcmp(atom.type,"Si")) strcpy(type,"Si");

	if(!strcmp(type,"Un")) return FALSE; 
	else return TRUE;
}

Group Molecule::Find_A_Group(int atom_id) const
{
        int i,j,id,num;
	bool mark;
        Group group;

        // define the center

        group.center=atom[atom_id-1];

        // find the center's neighbours

        num=group.center.num_neib;

        for(i=0;i<num;i++)
                {
                 group.neib[i]=this->atom[group.center.neib[i]-1];
		 group.bond[i]=this->bond[group.center.bond[i]-1];
                }

        // count the necessary parameters

	group.num_neib=group.center.num_neib;
        group.num_h=0; group.num_nonh=0;

        for(i=0;i<num;i++)
                {
                 if(group.neib[i].type[0]=='H') group.num_h++;
                 else group.num_nonh++;
                }

        group.num_hetero=0;

        for(i=0;i<group.num_nonh;i++)
                {
                 if(!strcmp(group.neib[i].type,"F")) group.num_hetero++;
		 else if(!strcmp(group.neib[i].type,"Cl")) group.num_hetero++;
		 else if(!strcmp(group.neib[i].type,"Br")) group.num_hetero++; 
		 else if(!strcmp(group.neib[i].type,"I")) group.num_hetero++;
		 else if(!strcmp(group.neib[i].type,"Si")) continue; 
                 else if(group.neib[i].type[0]=='N') group.num_hetero++;
                 else if(group.neib[i].type[0]=='O') group.num_hetero++;
                 else if(group.neib[i].type[0]=='P') group.num_hetero++;
                 else if(group.neib[i].type[0]=='S') group.num_hetero++;
                 else continue;
                }

        group.db_type=0; group.num_db=0; group.num_tb=0;

        for(i=0;i<group.num_nonh;i++)
                {
                 if(!strcmp(group.bond[i].type,"2"))
                        {
			 group.num_db++;
                         if(group.neib[i].type[0]=='C') group.db_type=1;
                         else if(group.neib[i].type[0]=='N') group.db_type=2;
                         else if(group.neib[i].type[0]=='O') group.db_type=3;
                         else if(group.neib[i].type[0]=='S') group.db_type=4;
                         else if(group.neib[i].type[0]=='P') group.db_type=5;
                         else continue;
                        }
		 else if(!strcmp(group.bond[i].type,"1")&&
			 !strcmp(group.neib[i].type,"O.co2"))
			{
			 group.db_type=3; group.num_db++;
			}
		 else if(!strcmp(group.bond[i].type,"ar")&&
                         !strcmp(group.neib[i].type,"O.co2"))
                        {
                         group.db_type=3; group.num_db++;
                        }
		 else if(!strcmp(group.bond[i].type,"3"))
                        {
                         group.num_tb++;
                        }
                 else continue;
                }

        group.num_pi=0;

        for(i=0;i<group.num_nonh;i++)
                {
                 if(!strcmp(group.bond[i].type,"2")) continue;
                 else if(!strcmp(group.bond[i].type,"3")) continue;
                 else if(!strcmp(group.neib[i].type,"C.ar")) group.num_pi++;
                 else if(!strcmp(group.neib[i].type,"C.2")) group.num_pi++;
                 else if(!strcmp(group.neib[i].type,"C.1"))group.num_pi++;
                 else if(!strcmp(group.neib[i].type,"C.cat")) group.num_pi++;
                 else if(!strcmp(group.neib[i].type,"N.2")) group.num_pi++;
                 else if(!strcmp(group.neib[i].type,"N.1")) group.num_pi++;
                 else if(!strcmp(group.neib[i].type,"N.ar")) group.num_pi++;
                 else continue;
                }

	// check if the central atom is adjacent to any -SO-, -PO-, or -CO- 

	mark=false;

        for(i=0;i<group.num_nonh;i++)
                {
                 if(!strcmp(group.neib[i].type,"P.3")||
		    !strcmp(group.neib[i].type,"S.o")||
		    !strcmp(group.neib[i].type,"S.o2")||
		    !strcmp(group.neib[i].type,"C.2"))
			{
			 num=group.neib[i].num_nonh; 

			 for(j=0;j<num;j++)
				{
				 id=this->atom[group.neib[i].id-1].neib[j];
				 if(id==group.center.id) continue;
				 else if(!strcmp(this->atom[id-1].type,"O.2"))
					{
					 mark=true; break;
					}
				 else if(!strcmp(this->atom[id-1].type,"O.co2"))
                                        {
                                         mark=true; break;
                                        }
				 else continue;
				}

			 if(mark==true) break;
			 else continue;
			}
                 else continue;
                }

	if(mark==false) group.amide=0;
	else group.amide=j+1;   // assign the value of the amide bond

        group.valid=1; return group;
}

int Molecule::Two_Bonds_Connection_Check(Bond bond1, Bond bond2) const
// this function returns the ID of the joint atom
{
        int id;

        id=bond1.atom_1;

        if(id==bond2.atom_1) return id;
        else if(id==bond2.atom_2) return id;

        id=bond1.atom_2;

        if(id==bond2.atom_1) return id;
        else if(id==bond2.atom_2) return id;

        return 0;       // two bonds are not connected
}

int Molecule::Connection_In_Ring_Check(int id1, int id2) const
{
	int i,tmp;

	for(i=0;i<atom[id1-1].num_neib;i++)
		{
		 if(atom[id1-1].neib[i]!=id2) continue;
		 else 
			{
			 tmp=atom[id1-1].bond[i]-1;
			 if(bond[tmp].ring==0) continue;
			 else return TRUE;
			}
		}

	return FALSE;
}

int Molecule::Inner_Collision_Check() const
{
	int i,j;
	float d,d0;

	for(i=0;i<num_atom-1;i++)
	for(j=i+1;j<num_atom;j++)
		{
		 d=Distance(atom[i].coor,atom[j].coor);	
		 d0=atom[i].r+atom[j].r;

		 if(d>(d0-0.50)) continue;	// no bump;
		 else if(Connection_1_2_Check(atom[i].id, atom[j].id)==TRUE) 
			continue; 
		 else if(Connection_1_3_Check(atom[i].id, atom[j].id)==TRUE) 
			continue;
		 else if(Connection_1_4_Check(atom[i].id, atom[j].id)==TRUE) 
			continue;
		 else return TRUE; 
		}

	return FALSE;
}

int Molecule::Connection_1_2_Check(int id1, int id2) const
{
	int i;

	if(id1==id2) return FALSE;

	for(i=0;i<atom[id2-1].num_neib;i++)
		{
		 if(id1==atom[id2-1].neib[i]) return TRUE; 
		 else continue;
		}
	
	return FALSE;
}

int Molecule::Connection_1_3_Check(int id1, int id2) const
{
	int i,j;

	if(id1==id2) return FALSE;
	if(Connection_1_2_Check(id1,id2)==TRUE) return FALSE;

	id1--; id2--;

	for(i=0;i<atom[id1].num_neib;i++)
	for(j=0;j<atom[id2].num_neib;j++)
		{
		 if(atom[id1].neib[i]==atom[id2].neib[j]) return TRUE;
		 else continue;
		}

	return FALSE;
}

int Molecule::Connection_1_4_Check(int id1, int id2) const
{
	int i,j;

	if(id1==id2) return FALSE;
	if(Connection_1_2_Check(id1,id2)==TRUE) return FALSE;
	if(Connection_1_3_Check(id1,id2)==TRUE) return FALSE;

	id1--; id2--;

	for(i=0;i<atom[id1].num_neib;i++)
	for(j=0;j<atom[id2].num_neib;j++)
		{
		 if(Connection_1_2_Check(atom[id1].neib[i],
		    atom[id2].neib[j])==TRUE) return TRUE;
		 else continue;
		}

	return FALSE;
}

int Molecule::Connection_1_5_Check(int id1, int id2) const
{
        int i,j;

        if(id1==id2) return FALSE;
        if(Connection_1_2_Check(id1,id2)==TRUE) return FALSE;
        if(Connection_1_3_Check(id1,id2)==TRUE) return FALSE;
        if(Connection_1_4_Check(id1,id2)==TRUE) return FALSE;

        id1--; id2--;

        for(i=0;i<atom[id1].num_neib;i++)
        for(j=0;j<atom[id2].num_neib;j++)
                {
                 if(Connection_1_3_Check(atom[id1].neib[i],
                    atom[id2].neib[j])==TRUE) return TRUE;
                 else continue;
                }

        return FALSE;
}

int Molecule::Connection_1_6_Check(int id1, int id2) const
{
        int i,j;

        if(id1==id2) return FALSE;
        if(Connection_1_2_Check(id1,id2)==TRUE) return FALSE;
        if(Connection_1_3_Check(id1,id2)==TRUE) return FALSE;
        if(Connection_1_4_Check(id1,id2)==TRUE) return FALSE;
        if(Connection_1_5_Check(id1,id2)==TRUE) return FALSE;

        id1--; id2--;

        for(i=0;i<atom[id1].num_neib;i++)
        for(j=0;j<atom[id2].num_neib;j++)
                {
                 if(Connection_1_4_Check(atom[id1].neib[i],
                    atom[id2].neib[j])==TRUE) return TRUE;
                 else continue;
                }

        return FALSE;
}

void Molecule::Detect_Rings()
// find the ring systems in the given molecule
// define all the atoms in a ring as atom.ring=1; otherwise atom.ring=0;
{
	int i,j,mark,count;
	Group group;

	// first of all, check whether there is any ring in the molecule

	if((num_bond-num_atom+num_subst)==0)
		{
		 for(i=0;i<num_atom;i++) atom[i].ring=0;
		 for(i=0;i<num_bond;i++) bond[i].ring=0;
		 return;
		}

	// detect terminal atoms

	for(i=0;i<num_atom;i++) 	
		{
		 if(atom[i].num_neib<=1) atom[i].ring=0;
		 else atom[i].ring=-1;	// still ambiguous here 
		}

	// collapse the structure to eliminate branches

	do		
		{
		 mark=0;
		 for(i=0;i<num_atom;i++)
			{
			 if(atom[i].ring!=-1) continue;

			 count=0;	// count possible ring neighbors	
		 	 for(j=0;j<atom[i].num_neib;j++)
			{
			 if(atom[atom[i].neib[j]-1].ring==0) continue;
			 else count++;
			}

			 if(count<=1) {atom[i].ring=0; mark++;}
			}
		}
	while(mark);

	// detect branching bonds

	mark=0;

	for(i=0;i<num_bond;i++) 
		{
		 if(bond[i].valid<=0) continue; 
		 else if(atom[bond[i].atom_1-1].ring==0) bond[i].ring=0;
		 else if(atom[bond[i].atom_2-1].ring==0) bond[i].ring=0;
		 else {bond[i].ring=-1; mark++;}	// ambiguous
		}

	if(mark==0) return; 	// finished, no ring detected

	// now comes the ring perception algorithm

	int *atom_path=NULL; atom_path=new int[num_bond];
	if(atom_path==NULL) Memory_Allocation_Error();

	int *bond_path=NULL; bond_path=new int[num_bond];
	if(bond_path==NULL) Memory_Allocation_Error();

	for(i=0;i<num_bond;i++)	
		{
		 if(bond[i].ring!=-1) continue;

		 for(j=0;j<num_bond;j++) atom_path[j]=bond_path[j]=0;

		 count=Look_For_A_Ring(bond[i].id,atom_path,bond_path);

		 if(count==FALSE) {bond[i].ring=0; continue;}

		 for(j=0;j<count;j++) 
			{
			 bond[bond_path[j]-1].ring=1;
			 atom[atom_path[j]-1].ring=1;
			}
		}

	// define all the left atoms which are not in ring

	for(i=0;i<num_atom;i++)
		{
		 if(atom[i].ring!=-1) continue;
		 else atom[i].ring=0;
		}

	if(atom_path) delete [] atom_path; 
	if(bond_path) delete [] bond_path;

	return;
}

void Molecule::Reset_Choices(int bond_id, int choice[]) const
// find all the neighboring bonds for the given bond 
// it is exclusively used by the ring perception algorithm
{
	int i,tmp;

	for(i=0;i<MAX_BOND_NEIB;i++) choice[i]=0;

	for(i=0;i<bond[bond_id-1].num_neib;i++)
		{
		 tmp=bond[bond_id-1].neib[i];
		 if(bond[tmp-1].ring==0) choice[i]=0;
		 else if(tmp>num_bond) choice[i]=0;
		 else choice[i]=tmp; 
		}

	return;
}

void Molecule::Clean_Choices(int tried_choice, int choice[]) const
{
        int i;

        for(i=0;i<MAX_BOND_NEIB;i++)
                {
                 if(choice[i]==0) continue;
                 else if(choice[i]==tried_choice) choice[i]=0;
                 else continue;
                }

        return;
}

int Molecule::Get_Choices(int choice[], int wrong_end, 
                          int current_path_length, int current_path[]) const
// notice that wrong_end helps to define the direction of the searching path 
{
	int i,j,result;
	bool mark;

	result=0;

	for(i=0;i<MAX_BOND_NEIB;i++)
	{
	 if(choice[i]==0) continue;
	 else if(choice[i]>num_bond) continue;
	 else if(bond[choice[i]-1].ring==0) continue;
	 else if(bond[choice[i]-1].atom_1==wrong_end) continue;
         else if(bond[choice[i]-1].atom_2==wrong_end) continue;

	// Is this bond already included in the current path?

	 mark=false;

	 for(j=0;j<current_path_length;j++)
		{
		 if(current_path[j]!=choice[i]) continue;
		 else {mark=true; break;}
		}

	 if(mark==true) continue;
	 else {result=choice[i]; break;}
	}

	return result;
}

int Molecule::Look_For_A_Ring(int bond_id, int atom_path[], int bond_path[],
                              int required_size) const
// this is a bond-based algorithm to detect rings in a molecule.
// it returns the size of the ring (if zero, the given bond is not in a ring)
// also the corresponding atom path and bond path
{
	int i,j,ring_size=0;
	int p,next;

	// initialize the variables first

	int **choice=NULL; choice=new int*[num_bond]; 
	if(choice==NULL) Memory_Allocation_Error();
	for(i=0;i<num_bond;i++) choice[i]=new int[MAX_BOND_NEIB];

	for(i=0;i<num_bond;i++)
	for(j=0;j<MAX_BOND_NEIB;j++) choice[i][j]=0;

	// note that atom_path[] and bond_path[] is not cleaned here
	// they are supposed to be cleaned at the upper level 

	for(i=0;i<num_bond;i++) 
		{
		 if(bond[i].valid<=0) continue;
		 Reset_Choices(bond[i].id,choice[i]);
		 // printf("Choices for bond %d: ", bond[i].id);
		 // for(j=0;j<MAX_BOND_NEIB;j++) printf("%d ", choice[i][j]);
		 // printf("\n");
		}

	// put the given bond at the beginning of the searching chain 

	p=0; 
	bond_path[0]=bond_id; atom_path[0]=bond[bond_id-1].atom_1;

	for(;;)
	{
	 // search the next possible step 

	 for(;;) 
		{
		 next=Get_Choices(choice[bond_path[p]-1],atom_path[p],p+1,bond_path);

		 if(next==FALSE)  	// next step is not available 
			{
			 if(p==0) 	// the given bond is not in a ring.
				{
				 ring_size=0;
				 goto Look_For_A_Ring_Exit;
				}

			 // trace back to the previous step 

			 Reset_Choices(bond_path[p],choice[bond_path[p]-1]);
			 bond_path[p]=0;
			 atom_path[p]=0;
			 p--;
			 continue;
			}
		 else break;
		} 

	 // define the next step 

	 p++;
	 bond_path[p]=next;
	 atom_path[p]=Two_Bonds_Connection_Check(bond[bond_path[p-1]-1],
				                 bond[bond_path[p]-1]);

	 // record the searching history, prevent repeating  

	 Clean_Choices(bond_path[p],choice[bond_path[p-1]-1]);

	 // now check if a ring has formed

	 if(Two_Bonds_Connection_Check(bond[bond_path[0]-1],
            bond[bond_path[p]-1])==atom_path[0])
		{
		 if(required_size<=0) // a ring is found
		{
                 // printf("Ring size %d: ", p+1);
                 // for(i=0;i<=p;i++) printf("%d ", atom_path[i]);
                 // printf("\n");
		 ring_size=p+1; goto Look_For_A_Ring_Exit;
                }
		 else if((p+1)==required_size)  // a ring with required size
		{
		 // printf("Ring size %d: ", p+1);
		 // for(i=0;i<=p;i++) printf("%d ", atom_path[i]);
		 // printf("\n"); 
		 ring_size=p+1; goto Look_For_A_Ring_Exit;	 
		}
		 else   // pre-mature ring, trace back
		{
		 Reset_Choices(bond_path[p], choice[bond_path[p]-1]);
		 bond_path[p]=0;
		 atom_path[p]=0;
		 p--;
		 continue;
		}
		}
	 else   // no ring is formed at this step
		{
		 if(required_size<=0) continue;
		 else if((p+1)==required_size) // no ring on this path, trace back
		{
		 Reset_Choices(bond_path[p], choice[bond_path[p]-1]);
		 bond_path[p]=0;
		 atom_path[p]=0;
		 p--;
		 continue;
		}
		 else continue;
		}
	}

	Look_For_A_Ring_Exit: 

	for(i=0;i<num_bond;i++) 
		{
		 if(choice[i]) delete [] choice[i];
		}

	return ring_size;
}

// *************************************************************************
// find all the atoms in aromatic rings (6-member and 5-member)
// define them as atom.ring=2 
// *************************************************************************
void Molecule::Detect_Aromatic_Rings()
{
	int i,j;
	Ring tmp_ring;

	int *atom_path=NULL; atom_path=new int[6];
	if(atom_path==NULL) Memory_Allocation_Error();

	int *bond_path=NULL; bond_path=new int[6];
	if(bond_path==NULL) Memory_Allocation_Error();

	this->ring.clear();

	// check all the 6-membered rings first 

	for(i=0;i<num_bond;i++)
	{
 	 if(bond[i].ring!=1) continue;

	 for(j=0;j<6;j++) atom_path[j]=bond_path[j]=0;

	 if(Look_For_A_Ring(bond[i].id,atom_path,bond_path,6)==0) continue;
	 if(Aromatic_Ring_Check_6(atom_path,bond_path)==0) continue;

	 tmp_ring.Clear();

	 for(j=0;j<6;j++)
		{
		 atom[atom_path[j]-1].ring=2;
		 bond[bond_path[j]-1].ring=2;

		 tmp_ring.atom_id.push_back(atom_path[j]);
		 tmp_ring.bond_id.push_back(bond_path[j]);
		 tmp_ring.centroid[0]+=atom[atom_path[j]-1].coor[0];
		 tmp_ring.centroid[1]+=atom[atom_path[j]-1].coor[1];
		 tmp_ring.centroid[2]+=atom[atom_path[j]-1].coor[2];
		}

	 tmp_ring.num_member=tmp_ring.atom_id.size();

	 if(tmp_ring.num_member>0)
		{
	 	 tmp_ring.centroid[0]/=tmp_ring.num_member;
	 	 tmp_ring.centroid[1]/=tmp_ring.num_member;
	 	 tmp_ring.centroid[2]/=tmp_ring.num_member;

	 	 tmp_ring.valid=1; tmp_ring.type=2; 
		 this->ring.push_back(tmp_ring); 
		}
	}

	// then check all the 5-membered rings

	for(i=0;i<num_bond;i++)
	{
 	 if(bond[i].ring!=1) continue;

	 for(j=0;j<5;j++) atom_path[j]=bond_path[j]=0;

	 if(Look_For_A_Ring(bond[i].id,atom_path,bond_path,5)==0) continue;
	 if(Aromatic_Ring_Check_5(atom_path,bond_path)==0) continue;

	 tmp_ring.Clear();

	 for(j=0;j<5;j++)
		{
		 atom[atom_path[j]-1].ring=2;
		 bond[bond_path[j]-1].ring=2;

		 tmp_ring.atom_id.push_back(atom_path[j]);
		 tmp_ring.bond_id.push_back(bond_path[j]);
		 tmp_ring.centroid[0]+=atom[atom_path[j]-1].coor[0];
		 tmp_ring.centroid[1]+=atom[atom_path[j]-1].coor[1];
		 tmp_ring.centroid[2]+=atom[atom_path[j]-1].coor[2];
		}

	 tmp_ring.num_member=tmp_ring.atom_id.size();

	 if(tmp_ring.num_member>0)
		{
         	 tmp_ring.centroid[0]/=tmp_ring.num_member;
         	 tmp_ring.centroid[1]/=tmp_ring.num_member;
         	 tmp_ring.centroid[2]/=tmp_ring.num_member;

	 	 tmp_ring.valid=1; tmp_ring.type=2; 
		 this->ring.push_back(tmp_ring); 
		}
	}

	this->num_ring=this->ring.size();

	if(atom_path) delete [] atom_path; 
	if(bond_path) delete [] bond_path; 

	return;
}

int Molecule::Aromatic_Ring_Check_6(int atom_path[],int bond_path[]) const
{
	int i,count;
	bool mark,mark1,mark2;

	// check the number of pi electrons 

	mark=true; count=0;

	for(i=0;i<6;i++)
		{
		 if(!strcmp(atom[atom_path[i]-1].type,"C.ar")) 
			{
			 count++; 
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"N.ar")) 
			{
			 count++; 
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"C.2")) 
			{
			 count++; 
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"N.2")) 
			{
			 count++; 
			}
		 else 
			{
			 mark=false; break;
			}
		}

	if(mark==false||count!=6) return FALSE;

	// check if there are two continuous single bonds on the path
	// this algorithm is superior to check -=-=-=
	
	mark=false;

	for(i=0;i<5;i++)
		{
		 if(!strcmp(bond[bond_path[i]-1].type,"1")||
                    !strcmp(bond[bond_path[i]-1].type,"am")) mark1=true;
		 else mark1=false;

		 if(!strcmp(bond[bond_path[i+1]-1].type,"1")||
                    !strcmp(bond[bond_path[i+1]-1].type,"am")) mark2=true;
		 else mark2=false;
	
		 if(mark1&&mark2) {mark=true; break;}
		 else continue;
		}

	if(mark==true) return FALSE;
	else return TRUE;
}

int Molecule::Aromatic_Ring_Check_5(int atom_path[],int bond_path[]) const
{
	int i,count;
	bool mark,mark1,mark2;
	int pi_path[5];

	// check the number of pi electrons 

	mark=true; count=0;

	for(i=0;i<5;i++)
		{
		 if(!strcmp(atom[atom_path[i]-1].type,"C.ar")) 
			{
			 count++; pi_path[i]=1;
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"N.ar")) 
			{
			 count++; pi_path[i]=1;
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"C.2")) 
			{
			 count++; pi_path[i]=1;
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"N.2")) 
			{
			 count++; pi_path[i]=1;
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"N.pl3")) 
			{
			 count+=2; pi_path[i]=2;
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"N.am")) 
			{
			 count+=2; pi_path[i]=2;
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"O.3")) 
			{
			 count+=2; pi_path[i]=2;
			}
		 else if(!strcmp(atom[atom_path[i]-1].type,"S.3")) 
			{
			 count+=2; pi_path[i]=2;
			}
		 else 
			{
			 mark=false; break;
			}
		}

	if(mark==false||count!=6) return FALSE;

	// check if there are two continuous single bonds on the path
	// but it is okay if these two single bonds explain the special atom
	
	mark=false;

	for(i=0;i<4;i++)
		{
		 if(!strcmp(bond[bond_path[i]-1].type,"1")||
                    !strcmp(bond[bond_path[i]-1].type,"am")) mark1=true;
		 else mark1=false;

		 if(!strcmp(bond[bond_path[i+1]-1].type,"1")||
                    !strcmp(bond[bond_path[i+1]-1].type,"am")) mark2=true;
		 else mark2=false;
	
		 if(mark1!=true||mark2!=true) continue;
		 else if(pi_path[i+1]==2) continue;
		 else {mark=true; break;}
		}

	if(mark==true) return FALSE;
	else return TRUE;

	return mark;
}

float Molecule::Get_Weight() const
{
        int i;
	float sum=0.000;

        for(i=0;i<num_atom;i++) 
		{
		 if(atom[i].valid==0) continue;
		 else sum+=(atom[i].weight);
		}

	return sum;
}

void Molecule::Get_Formula(char *string) const
{
	int i;
        int num[12]={0,0,0,0,0,0,0,0,0,0,0,0}; // C,H,N,O,P,S,F,Cl,Br,I,Si,Un
	char element[12][3];
	char tmp[10],tmp_formula[256];

	strcpy(element[0],"C");
	strcpy(element[1],"H");
	strcpy(element[2],"N");
	strcpy(element[3],"O");
	strcpy(element[4],"P");
	strcpy(element[5],"S");
	strcpy(element[6],"F");
	strcpy(element[7],"Cl");
	strcpy(element[8],"Br");
	strcpy(element[9],"I");
	strcpy(element[10],"Si");
	strcpy(element[11],"Un");

	for(i=0;i<num_atom;i++)
		{
		 if(atom[i].valid==0) num[11]++;
		 else if(!strcmp(atom[i].type,"C.3")) num[0]++;
		 else if(!strcmp(atom[i].type,"C.2")) num[0]++;
		 else if(!strcmp(atom[i].type,"C.1")) num[0]++;
		 else if(!strcmp(atom[i].type,"C.cat")) num[0]++;
                 else if(!strcmp(atom[i].type,"C.ar")) num[0]++;
                 else if(!strcmp(atom[i].type,"H")) num[1]++;
                 else if(!strcmp(atom[i].type,"H.spc")) num[1]++;
		 else if(!strcmp(atom[i].type,"N.4")) num[2]++;
                 else if(!strcmp(atom[i].type,"N.3")) num[2]++;
                 else if(!strcmp(atom[i].type,"N.2")) num[2]++;
		 else if(!strcmp(atom[i].type,"N.1")) num[2]++;
                 else if(!strcmp(atom[i].type,"N.ar")) num[2]++;
                 else if(!strcmp(atom[i].type,"N.pl3")) num[2]++;
		 else if(!strcmp(atom[i].type,"N.am")) num[2]++;
                 else if(!strcmp(atom[i].type,"O.3")) num[3]++;
		 else if(!strcmp(atom[i].type,"O.2")) num[3]++;
                 else if(!strcmp(atom[i].type,"O.co2")) num[3]++;
                 else if(!strcmp(atom[i].type,"P.3")) num[4]++;
		 else if(!strcmp(atom[i].type,"S.3")) num[5]++;
                 else if(!strcmp(atom[i].type,"S.2")) num[5]++;
                 else if(!strcmp(atom[i].type,"S.o")) num[5]++;
		 else if(!strcmp(atom[i].type,"S.o2")) num[5]++;
                 else if(!strcmp(atom[i].type,"F")) num[6]++;
                 else if(!strcmp(atom[i].type,"Cl")) num[7]++;
                 else if(!strcmp(atom[i].type,"Br")) num[8]++;
		 else if(!strcmp(atom[i].type,"I")) num[9]++;
		 else if(!strcmp(atom[i].type,"Si")) num[10]++;
                 else num[11]++;
		}

	strcpy(tmp_formula,"");

	for(i=0;i<12;i++)
		{
		 if(num[i]==0) continue;
		 else
			{
			 strcat(tmp_formula,element[i]);
			 if(num[i]==1) continue;
			 else 
				{
				 sprintf(tmp,"%d", num[i]);
				 strcat(tmp_formula,tmp);
				}
			}
		}

	strcpy(string,tmp_formula);

        return;
}

int Molecule::Get_Num_HB_Atom() const
{
	int i,num=0;

        for(i=0;i<num_atom;i++)
                {
                 if(atom[i].valid==0) continue;
                 else if(!strcmp(atom[i].hb,"D")) num++;
		 else if(!strcmp(atom[i].hb,"DA")) num++;
		 else if(!strcmp(atom[i].hb,"A")) num++;
                 else continue;
                }

	return num;
}

int Molecule::Get_Num_Heavy_Atom() const
{
	int i,num=0;

	for(i=0;i<num_atom;i++)
		{
		 if(atom[i].valid==0) continue;
		 else if(!strcmp(atom[i].type,"H")) continue;
		 else num++;
		}

	return num;
}

float Molecule::Count_Rotor()
// if a single bond is normal, bond.valid=1;
// if a single bond is a rotor, bond.valid=2;
{
	int i,mark,tmp;
	int id1,id2;
	float sum;

	// clean the variables

	for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

	// eliminate all the bonds in_ring and all the non-single bonds

	for(i=0;i<num_bond;i++)
		{
		 if(bond[i].ring>0) continue;
		 else if(strcmp(bond[i].type,"1")) continue;
		 else bond[i].valid=2;
		}

	// eliminate all the R-H, R-X, R-OH, R-NH2, R-CH3 bonds 

	for(i=0;i<num_bond;i++)
		{
		 if(bond[i].valid!=2) continue;

		 id1=bond[i].atom_1; id2=bond[i].atom_2;

		 if(atom[id1-1].num_nonh==1||
		    atom[id2-1].num_nonh==1) bond[i].valid=1;
		 else continue;
		}

	// sp2-sp2 rotors

	for(i=0;i<num_bond;i++)
		{
		 if(bond[i].valid!=2) continue;

		 id1=bond[i].atom_1; id2=bond[i].atom_2;

		 mark=0; 

		 tmp=Get_Atom_Hybridizing_Type(atom[id1-1].type);
		 if(tmp==1||tmp==2) mark++;

		 tmp=Get_Atom_Hybridizing_Type(atom[id2-1].type);
                 if(tmp==1||tmp==2) mark++;

		 if(mark==2) {bond[i].valid=1; continue;}
		}

	// eliminate terminal rotors, e.g. -PO3, -CF3, -CMe3, -NMe3

	for(i=0;i<num_bond;i++)
                {
                 if(bond[i].valid!=2) continue;

		 id1=bond[i].atom_1; id2=bond[i].atom_2;

		 if(Judge_Terminal_Atom(atom[id1-1],id2)==TRUE)
			{
			 bond[i].valid=1; 
			}
		 else if(Judge_Terminal_Atom(atom[id2-1],id1)==TRUE)
			{
			 bond[i].valid=1;
			}
		 else continue; 
		}

	// eliminate abnormal rotors

	for(i=0;i<num_bond;i++)
                {
                 if(bond[i].valid!=2) continue;

                 id1=bond[i].atom_1; id2=bond[i].atom_2;

		 if(atom[id1-1].valid<=0) bond[i].valid=1;
		 else if(atom[id2-1].valid<=0) bond[i].valid=1;
		 else continue;
		}

	// now count the frozen rotors, all the rotors have been labeled as 2

	sum=0.000;

	for(i=0;i<num_bond;i++)
                {
                 if(bond[i].valid!=2) continue;
		 else sum+=1.00;
                }

	return sum;
}

int Molecule::Judge_Terminal_Atom(Atom &atm, int partner) const
{
	int i,j;

	if(Get_Atom_Hybridizing_Type(atm.type)!=3) return FALSE;

	if(atm.num_nonh!=4) return FALSE;

	for(i=0;i<atm.num_nonh;i++)
		{
		 if(atm.neib[i]==partner) continue;
		 else if(this->atom[atm.neib[i]-1].num_nonh!=1) return FALSE;
		 else continue;
		}

	for(i=0;i<atm.num_nonh-1;i++)
                {
                 if(atm.neib[i]==partner) continue;

		 for(j=i+1;j<atm.num_nonh;j++)
			{
			 if(atm.neib[j]==partner) continue;
			 else if(!strcmp(this->atom[atm.neib[j]-1].xtype,
				 this->atom[atm.neib[i]-1].xtype)) continue;
			 else return FALSE;
			}
                }

	return TRUE;
}

void Molecule::Generate_Volume_Dots(float spacing, float probe_r)
// if probe_r=0, then it represents the VDW volume;
// if probe_r=1.40, then it represents the solvent-inaccessible volume;
// note that it is a unit-atom model
{
	int i,j,num,atom_id;
	float max_x,min_x,max_y,min_y,max_z,min_z;
	float coor[3],d,margin,dmin,unit;
	bool mark;
	Dot tmp_dot;

	vol_dot.clear();	// clear all the existing point, if any

	max_x=max_y=max_z=-9999.000;
	min_x=min_y=min_z= 9999.000;

	for(i=0;i<this->num_atom;i++)
	{
 	 if(this->atom[i].valid<=0) continue;
	 else if(this->atom[i].type[0]=='H') continue;

 	 if(this->atom[i].coor[0]>max_x) max_x=ceil(this->atom[i].coor[0]);
 	 if(this->atom[i].coor[0]<min_x) min_x=floor(this->atom[i].coor[0]);
 	 if(this->atom[i].coor[1]>max_y) max_y=ceil(this->atom[i].coor[1]);
       	 if(this->atom[i].coor[1]<min_y) min_y=floor(this->atom[i].coor[1]);
 	 if(this->atom[i].coor[2]>max_z) max_z=ceil(this->atom[i].coor[2]);
       	 if(this->atom[i].coor[2]<min_z) min_z=floor(this->atom[i].coor[2]);
	}

	// add a margin to the box, 3.00 is the largest possible atomic radius

	margin=probe_r+3.00;

	max_x+=margin; min_x-=margin; 
	max_y+=margin; min_y-=margin; 
	max_z+=margin; min_z-=margin; 

	// now check each dot

	unit=spacing*spacing*spacing;  // volume of each dot

	for(coor[0]=min_x;coor[0]<max_x;coor[0]+=spacing)
	for(coor[1]=min_y;coor[1]<max_y;coor[1]+=spacing)
	for(coor[2]=min_z;coor[2]<max_z;coor[2]+=spacing)
		{
		 mark=false; dmin=9999.000; atom_id=0;

		 for(i=0;i<this->num_atom;i++)
			{	
			 if(this->atom[i].valid<=0) continue;
			 else if((this->atom[i].type[0]=='H')&&
                                 (!strcmp(this->atom[i].hb,"N"))) continue;
			 // non-polar hydrogen are excluded

			 d=Distance(coor,this->atom[i].coor);

			 if(d>(this->atom[i].R+probe_r)) continue;
			 else if(d>=dmin) continue; 
			 else {dmin=d; atom_id=i+1; mark=true;} 
			}

		 if(mark==false) continue;

		 tmp_dot.valid=atom_id;
		 strcpy(tmp_dot.type,this->atom[atom_id-1].type);
		 tmp_dot.coor[0]=coor[0];
		 tmp_dot.coor[1]=coor[1];
		 tmp_dot.coor[2]=coor[2];

		 // correct contributions of the surfacial dots

		 if((this->atom[atom_id-1].R+probe_r-dmin)<(spacing/2.0))
			{		
			 tmp_dot.unit=unit/2.0;
			}
		 else
			{	
		 	 tmp_dot.unit=unit;
			}

		 vol_dot.push_back(tmp_dot);
		}

	// now finally sort all the points according to atom_id

	num=vol_dot.size();

	for(i=0;i<num-1;i++)
	for(j=i+1;j<num;j++)
		{
		 if(vol_dot[i].valid<=vol_dot[j].valid) continue;
		 else 
			{
			 tmp_dot=vol_dot[i];
			 vol_dot[i]=vol_dot[j];
			 vol_dot[j]=tmp_dot;
			}
		}

	return;
}

float Molecule::Calculate_Volume() const 
{
	if(vol_dot.size()==0) 
		{
		 puts("Volume dots have not been generated yet.");
		 puts("Program stops at Calculate_Volume().");
		 exit(1);
		}

	int i,num;
	float volume;

	num=vol_dot.size(); volume=0.000;

	for(i=0;i<num;i++)
		{
		 if(vol_dot[i].valid==0) continue;
		 else volume+=vol_dot[i].unit;
		}

	return volume;
}

void Molecule::Show_Volume_Dots(char *filename, char *show) const
// output volume dots to a PDB file
{
        FILE *fp;
	int i,num;
	float *property=NULL;

        if((fp=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	fprintf(fp,"COMPND    Molecular Volume Dots\n");
	fprintf(fp,"REMARK Creation time: %s\n", Get_Time());

	num=vol_dot.size();

	property=new float[num];
	if(property==NULL) Memory_Allocation_Error();

	if(!strcasecmp(show,"unit"))
		{
		 for(i=0;i<num;i++) property[i]=vol_dot[i].unit;
		}
	else
		{
		 for(i=0;i<num;i++) property[i]=vol_dot[i].score;
		}

	for(i=0;i<num;i++)
		{
		 if(vol_dot[i].valid==0) continue;
                 else if(!strcmp(vol_dot[i].type,"F")||
                         !strcmp(vol_dot[i].type,"Cl")||
                         !strcmp(vol_dot[i].type,"Br")||
                         !strcmp(vol_dot[i].type,"I"))
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "F");
                         fprintf(fp,"VOL   %-5d   ", vol_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 vol_dot[i].coor[0], 
				 vol_dot[i].coor[1], 
				 vol_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(!strcmp(vol_dot[i].type,"Si"))
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "SI");
                         fprintf(fp,"VOL   %-5d   ", vol_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 vol_dot[i].coor[0], 
				 vol_dot[i].coor[1], 
				 vol_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(vol_dot[i].type[0]=='C')
			{
			 fprintf(fp,"HETATM%6d %-3s ", i+1, "C");
                         fprintf(fp,"VOL   %-5d   ", vol_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 vol_dot[i].coor[0], 
				 vol_dot[i].coor[1], 
				 vol_dot[i].coor[2], 1.00, property[i]);
			}
		 else if(vol_dot[i].type[0]=='N')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "N");
                         fprintf(fp,"VOL   %-5d   ", vol_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 vol_dot[i].coor[0], 
				 vol_dot[i].coor[1], 
				 vol_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(vol_dot[i].type[0]=='O')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "O");
                         fprintf(fp,"VOL   %-5d   ", vol_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 vol_dot[i].coor[0], 
				 vol_dot[i].coor[1], 
				 vol_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(vol_dot[i].type[0]=='S')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "S");
                         fprintf(fp,"VOL   %-5d   ", vol_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 vol_dot[i].coor[0], 
				 vol_dot[i].coor[1], 
				 vol_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(vol_dot[i].type[0]=='P')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "P");
                         fprintf(fp,"VOL   %-5d   ", vol_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 vol_dot[i].coor[0], 
				 vol_dot[i].coor[1], 
				 vol_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(vol_dot[i].type[0]=='H')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "H");
                         fprintf(fp,"VOL   %-5d   ", vol_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 vol_dot[i].coor[0], 
				 vol_dot[i].coor[1], 
				 vol_dot[i].coor[2], 1.00, property[i]);
                        }
                 else
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "ZN");
                         fprintf(fp,"VOL   %-5d   ", vol_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 vol_dot[i].coor[0], 
				 vol_dot[i].coor[1], 
				 vol_dot[i].coor[2], 1.00, property[i]);
                        }
		}

        fclose(fp); if(property) delete [] property;

	return;
}

// ***************************************************************************
// generate dot surface of the given molecule (unit-atom model)
// if probe_r=0.00, it is vdw surface;
// if probe_r=WATER_R, it is solvent accessible surface.
// Latest update: 11/17/2003 
// ***************************************************************************
void Molecule::Generate_Surface_Dots(float probe_r)
{
	extern ForceField *ff;
	int i,j,k;
	bool mark;
	float d,dd,dmin;
	DotSet tmp_set;
	Dot tmp_dot;

	sur_dot.clear();

	int *ligand_check_list;

	ligand_check_list=new int[this->num_atom];
	if(ligand_check_list==NULL) Memory_Allocation_Error();	

	for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid<=0) continue;
		 else if(!strcmp(this->atom[i].xtype,"H")) continue;

		 for(j=0;j<this->num_atom;j++)
			{
			 if(i==j||this->atom[j].valid<=0) 
				{
				 ligand_check_list[j]=0; continue;
				}
			 else if(!strcmp(this->atom[j].xtype,"H")) 
				{
				 ligand_check_list[j]=0; continue;
				}

			 d=Distance(this->atom[i].coor,this->atom[j].coor);

			 if(d>(this->atom[i].R+this->atom[j].R+2*probe_r))
				{
				 ligand_check_list[j]=0;
				}
			 else
				{
				 ligand_check_list[j]=1;
				}
			}

		 tmp_set=ff->Get_Surface_Dot(this->atom[i],probe_r);

		 for(j=0;j<tmp_set.num_dot;j++) 
			{
			 // check whether this dot is on the surface or not 

			 mark=true; dmin=1000.0;

			 for(k=0;k<this->num_atom;k++)
			{
			 if(ligand_check_list[k]==0) continue;

			 d=Distance(tmp_set.dot[j].coor,this->atom[k].coor);
			 dd=d/(this->atom[k].R+probe_r)-1.00;

			 if(dd<0.000) {mark=false; break;} 
			 else if(dd<dmin) {dmin=dd; continue;}
			 else continue; 
			}

			 if(mark==false) tmp_set.dot[j].valid=0;
			 else if(dmin>=0.10) tmp_set.dot[j].valid=1; // regular 
			 else tmp_set.dot[j].valid=2;	// dots at edge
                        }

		 // now record the surface dots of the current atom

		 for(j=0;j<tmp_set.num_dot;j++)
			{
			 if(tmp_set.dot[j].valid==0) continue;

			 tmp_dot=tmp_set.dot[j]; tmp_dot.valid=i+1;

			 if(tmp_set.dot[j].valid==1)
				{
				 tmp_dot.unit=tmp_set.unit;
				}
			 else  // correct the overlapping of edging dots
				{
				 tmp_dot.unit=tmp_set.unit*0.500;
				}
			 
			 sur_dot.push_back(tmp_dot);
			}
		}

	if(ligand_check_list) delete [] ligand_check_list;

	return;
}

float Molecule::Calculate_Surface() const
{
	if(sur_dot.size()==0)
		{
		 puts("Surface dots have not been generated yet.");
		 puts("Program stops at Calculate_Surface().");
		 exit(1);
		}

	int i,num,type,id;
	float surface,nsur,psur;

	surface=psur=nsur=0.000;

	num=sur_dot.size();

	for(i=0;i<num;i++)
		{
		 if(sur_dot[i].valid==0) continue;

		 // first, determine whether it is non-polar or polar
		 // type=1, non-polar
		 // type=2, polar
		 // type=3, others

		 id=sur_dot[i].valid;

		 if(!strcmp(this->atom[id-1].hb,"N")) type=1;
		 else if(!strcmp(this->atom[id-1].hb,"DH")) type=2;
		 else if(!strcmp(this->atom[id-1].hb,"H")) type=1;
		 else if(!strcmp(this->atom[id-1].hb,"P")) type=1;
		 else if(!strcmp(this->atom[id-1].hb,"D")) type=2;
		 else if(!strcmp(this->atom[id-1].hb,"A")) type=2;
		 else if(!strcmp(this->atom[id-1].hb,"DA")) type=2;
		 else type=3;

		 switch(type)
			{
			 case 1: {nsur+=sur_dot[i].unit; break;}
			 case 2: {psur+=sur_dot[i].unit; break;}
			 default: {nsur+=sur_dot[i].unit; break;}
			}

		 surface+=sur_dot[i].unit;
		}

	return surface;
}

// *************************************************************************
// output surface dots to a PDB file
// the "show" flag could be either "unit" or "score"
// Latest update: 11/17/2003
// *************************************************************************
void Molecule::Show_Surface_Dots(char *filename, char *show) const
{
        FILE *fp;
	int i,num;
	float *property=NULL;

        if((fp=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	fprintf(fp,"COMPND    Molecular Surface Dots\n");
	fprintf(fp,"REMARK Creation time: %s\n", Get_Time());

	num=sur_dot.size();

	property=new float[num];
	if(property==NULL) Memory_Allocation_Error();

	if(!strcasecmp(show,"unit"))
		{
		 for(i=0;i<num;i++) property[i]=sur_dot[i].unit;
		}
	else
		{
		 for(i=0;i<num;i++) property[i]=sur_dot[i].score;
		}

	for(i=0;i<num;i++)
		{
		 if(sur_dot[i].valid==0) continue;
                 else if(!strcmp(sur_dot[i].type,"F")||
                         !strcmp(sur_dot[i].type,"Cl")||
                         !strcmp(sur_dot[i].type,"Br")||
                         !strcmp(sur_dot[i].type,"I"))
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "F");
                         fprintf(fp,"SUR   %-5d   ", sur_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 sur_dot[i].coor[0], 
				 sur_dot[i].coor[1], 
				 sur_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(!strcmp(sur_dot[i].type,"Si"))
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "SI");
                         fprintf(fp,"SUR   %-5d   ", sur_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 sur_dot[i].coor[0], 
				 sur_dot[i].coor[1], 
				 sur_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(sur_dot[i].type[0]=='C')
			{
			 fprintf(fp,"HETATM%6d %-3s ", i+1, "C");
                         fprintf(fp,"SUR   %-5d   ", sur_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 sur_dot[i].coor[0], 
				 sur_dot[i].coor[1], 
				 sur_dot[i].coor[2], 1.00, property[i]);
			}
		 else if(sur_dot[i].type[0]=='N')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "N");
                         fprintf(fp,"SUR   %-5d   ", sur_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 sur_dot[i].coor[0], 
				 sur_dot[i].coor[1], 
				 sur_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(sur_dot[i].type[0]=='O')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "O");
                         fprintf(fp,"SUR   %-5d   ", sur_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 sur_dot[i].coor[0], 
				 sur_dot[i].coor[1], 
				 sur_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(sur_dot[i].type[0]=='S')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "S");
                         fprintf(fp,"SUR   %-5d   ", sur_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 sur_dot[i].coor[0], 
				 sur_dot[i].coor[1], 
				 sur_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(sur_dot[i].type[0]=='P')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "P");
                         fprintf(fp,"SUR   %-5d   ", sur_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 sur_dot[i].coor[0], 
				 sur_dot[i].coor[1], 
				 sur_dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(sur_dot[i].type[0]=='H')
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "H");
                         fprintf(fp,"SUR   %-5d   ", sur_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 sur_dot[i].coor[0], 
				 sur_dot[i].coor[1], 
				 sur_dot[i].coor[2], 1.00, property[i]);
                        }
                 else 
                        {
                         fprintf(fp,"HETATM%6d %-3s ", i+1, "ZN");
                         fprintf(fp,"SUR   %-5d   ", sur_dot[i].valid);
                         fprintf(fp,"%7.3f %7.3f %7.3f %5.2f %5.2f\n",
                                 sur_dot[i].coor[0], 
				 sur_dot[i].coor[1], 
				 sur_dot[i].coor[2], 1.00, property[i]);
                        }
		}

        fclose(fp); if(property) delete [] property;

	return;
}

