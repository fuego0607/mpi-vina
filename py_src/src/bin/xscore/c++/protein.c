# include "xtool.h"

Protein::Protein()
{
	Clear();
}

Protein::~Protein()
{
	atom.clear(); chain.clear(); ring.clear();
	vol_dot.clear(); sur_dot.clear();
}

Protein::Protein(const Protein &original)
{
 this->xtool_format=original.xtool_format;
 strcpy(this->name,original.name);
 this->surface=original.surface;
 this->bnsur=original.bnsur;
 this->bpsur=original.bpsur;

 int i,n;

 this->num_atom=original.num_atom;
 this->atom.clear(); n=original.num_atom;
 for(i=0;i<n;i++) this->atom.push_back(original.atom[i]);

 this->num_chain=original.num_chain;
 this->chain.clear(); n=original.num_chain;
 for(i=0;i<n;i++) this->chain.push_back(original.chain[i]);

 this->num_ring=original.num_ring;
 this->ring.clear(); n=original.num_ring;
 for(i=0;i<n;i++) this->ring.push_back(original.ring[i]);

 this->vol_dot.clear(); n=original.vol_dot.size();
 for(i=0;i<n;i++) this->vol_dot.push_back(original.vol_dot[i]);

 this->sur_dot.clear(); n=original.sur_dot.size();
 for(i=0;i<n;i++) this->sur_dot.push_back(original.sur_dot[i]);
}

Protein& Protein::operator = (const Protein &original)
{
 if(this==&original) return *this;

 this->xtool_format=original.xtool_format;
 strcpy(this->name,original.name);
 this->surface=original.surface;
 this->bnsur=original.bnsur;
 this->bpsur=original.bpsur;

 int i,n;

 this->num_atom=original.num_atom;
 this->atom.clear(); n=original.num_atom;
 for(i=0;i<n;i++) this->atom.push_back(original.atom[i]);

 this->num_chain=original.num_chain;
 this->chain.clear(); n=original.num_chain;
 for(i=0;i<n;i++) this->chain.push_back(original.chain[i]);

 this->num_ring=original.num_ring;
 this->ring.clear(); n=original.num_ring;
 for(i=0;i<n;i++) this->ring.push_back(original.ring[i]);

 this->vol_dot.clear(); n=original.vol_dot.size();
 for(i=0;i<n;i++) this->vol_dot.push_back(original.vol_dot[i]);

 this->sur_dot.clear(); n=original.sur_dot.size();
 for(i=0;i<n;i++) this->sur_dot.push_back(original.sur_dot[i]);

 return *this;
}

void Protein::Clear()
{
	strcpy(name,"");
	surface=bnsur=bpsur=0.000;
	num_atom=0; atom.clear();
	num_chain=0; chain.clear();
	num_ring=0; ring.clear();
	vol_dot.clear(); sur_dot.clear();
}

void Protein::Show_Contents() const
{
	printf("\n%s\n", this->name);

	Show_Chains();
	Show_Atoms();
	Show_Rings();

	return;
}

void Protein::Show_Chains() const
{
	int i,j,k,temp,count;

	printf("\nTotal number of chains in this protein = %d\n", num_chain);

	for(i=0;i<this->num_chain;i++)
	{
	 count=(int)ceil(this->chain[i].length/13.0);

	 for(j=1;j<=count;j++)
		{
		 printf("SEQRES %3d %c %4d  ", 
			j, this->chain[i].label, this->chain[i].length);

		 for(k=1;k<=13;k++)
			{
			 temp=(j-1)*13+k;
			 if(temp>this->chain[i].length) break;
			 else printf("%3s ", this->chain[i].residue[temp-1].name);
			}
		
		 printf("\n");
		}
	}

	return;
}

void Protein::Show_Atoms() const
{
	int i;

	printf("\nTotal number of atoms in this protein = %d\n",num_atom);

        for(i=0;i<num_atom;i++)
                {
                 printf("%-5d ", atom[i].id);
		 printf("%-4s ", atom[i].name);
		 printf("%-3s ", atom[i].residue);
		 printf("%4s ", atom[i].res_id);
		 printf("%c ", atom[i].chain);
		 printf("%7.2f ", atom[i].coor[0]);
		 printf("%7.2f ", atom[i].coor[1]);
		 printf("%7.2f ", atom[i].coor[2]);
		 printf("%-5s ", atom[i].type);
		 printf("%4.2f ", atom[i].R);
		 printf("%5.2f ", atom[i].q);
		 printf("=%1d ", atom[i].valid);
		 printf("%-7s ", atom[i].xtype);
		 printf("%-2s ", atom[i].hb);
		 printf("%5.2f ", atom[i].logp);
		 printf("(%1d) ", atom[i].ring);
		 printf("%7.2f ", atom[i].root[0]);
                 printf("%7.2f ", atom[i].root[1]);
                 printf("%7.2f ", atom[i].root[2]);
		 printf("%-5d ", atom[i].neib[0]);
		 printf("\n");                
		}

	return;
}

void Protein::Show_Rings() const
{
	int i,j;

	printf("\nTotal number of rings in this protein = %d\n", num_ring);

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

// ****************************************************************************
// Read information from a given PDB file
// 'flag' determines if the CONECT information is read or not
// latest update: 10/28/2003
// ****************************************************************************
void Protein::Read_From_PDB(char *filename)
{
	FILE *fp;
	int i,j,n,count,len;
	char line[256],head[80],temp[256];
	bool mark,mark1;

	if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

	Clear();	// clear this protein

	// first, check the PDB file 

	count=0; this->xtool_format=false;	 

	mark=false;	// SEQRES information indicator
	mark1=false;	// COMPND information indicator 

	for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else
                        {
                         sscanf(line,"%s", head);
                         if(!strcmp(head,"ATOM")) count++;
			 else if(!strcmp(head,"HETATM")) count++;
			 else if(!strcmp(head,"SEQRES")) mark=true;
			 else if(!strcmp(head,"COMPND")&&mark1==false)
				{
				 len=strlen(line);

				 if(len<70) n=len;
				 else n=70;

				 j=0;
				 for(i=10;i<n;i++) 
					{
					 if(line[i]=='\n') continue;
					 else {temp[j]=line[i];j++;}
					}
				 temp[j]='\0';

				 strcpy(this->name,temp); mark1=true;

				 // note that sometimes COMPND spans more than
				 // one line, but we only take the first line.
				 // why the hell the name should be that long?
				}
			 else if(!strcmp(head,"REMARK"))
				{
				 if(strstr(line,"X-Tool")||
                                    strstr(line,"X-TOOL")||
				    strstr(line,"XTool")||
				    strstr(line,"XTOOL"))
					{
					 this->xtool_format=true;
					}
				}
                         else continue;
                        }
                }

        rewind(fp);

	if(count<5) PDB_Format_Error(filename);

	// read the ATOM/HETATM section 

	this->Read_ATOM(fp);

	// read the input sequence information

	if(mark==true) this->Read_SEQRES(fp);
	else this->Analyze_Sequence();

	// read the CONECT section if it is generated by X-Tool itself

	if(this->xtool_format==true) this->Read_CONECT(fp);

	fclose(fp); return;
}

// ****************************************************************************
// Read the sequence information from a given PDB file
// if the sequence is broken (often seen in SYBYL outputs), it can also fix it. 
// note that the length of each chain is re-counted according to the real
// number of residues
// latest update: 03/04/2004
// ****************************************************************************
void Protein::Read_SEQRES(FILE *fp)
{
	char line[256],head[256],temp[256],label,residue[10];
	int i,id,len;
	char *p;
	Residue tmp_residue;
	Chain tmp_chain;

	// clear the record first

	this->chain.clear(); this->num_chain=0;

	while(fgets(line,256,fp))
	{
	 strcpy(temp,line); temp[6]='\0'; 
	 strcpy(head,""); sscanf(temp,"%s",head);

	 if(strcmp(head,"SEQRES")) continue;

	 label=line[11]; line[11]=' ';

	 // now check if this line starts a new chain

	 id=-1;

	 for(i=0;i<this->num_chain;i++)
		{
		 if(label!=this->chain[i].label) continue;
		 else {id=i; break;}
		}

	 if(id>=0)  // an existing chain 
		{
		 for(i=0;i<=17;i++) line[i]=' '; line[70]='\0';

		 for(;;)
			{
			 strcpy(residue,"XXX"); sscanf(line,"%s", residue);
			 if(!strcmp(residue,"XXX")) break; 
			 else
				{
				 tmp_residue.Clear();
				 tmp_residue.valid=1;
				 strcpy(tmp_residue.name,residue);
				 tmp_residue.chain=label;
				 strcpy(tmp_residue.id,"0"); // not known yet

				 this->chain[id].residue.push_back(tmp_residue);
				 this->chain[id].length++;

				 // erase this residue from the line

				 len=strlen(residue);
			 	 p=strstr(line,residue);
				 for(i=1;i<=len;i++) {*p=' '; p++;}
				}
			}
		}
	 else  // a new chain
		{
		 tmp_chain.Clear();
		 tmp_chain.label=label;
		 tmp_chain.valid=1;

		 this->chain.push_back(tmp_chain);
		 this->num_chain=this->chain.size();
		 id=this->chain.size()-1;

		 for(i=0;i<=17;i++) line[i]=' '; line[70]='\0';

                 for(;;)
                        {
                         strcpy(residue,"XXX"); sscanf(line,"%s", residue);
                         if(!strcmp(residue,"XXX")) break;
                         else
                                {
				 tmp_residue.Clear();
                                 tmp_residue.valid=1;
                                 strcpy(tmp_residue.name,residue);
                                 tmp_residue.chain=label;
                                 strcpy(tmp_residue.id,"0");  // not known yet

                                 this->chain[id].residue.push_back(tmp_residue);
				 this->chain[id].length++;

                                 // erase this residue from the line

                                 len=strlen(residue);
                                 p=strstr(line,residue);
                                 for(i=1;i<=len;i++) {*p=' '; p++;}
                                }
                        }
		}
	}

	rewind(fp); return;
}

// ****************************************************************************
// Generate sequence information based on atom information 
// latest update: 03/04/2004
// ****************************************************************************
void Protein::Analyze_Sequence()
{
	int i,j,id,total;
	bool mark;
	Residue tmp_residue;
	Chain tmp_chain;

	this->chain.clear(); this->num_chain=0;

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].part>1) continue;   // hetero atom

	 // first, check if this atom belongs to a known chain

	 mark=false; id=-1; total=this->chain.size();

	 for(j=0;j<total;j++)
		{
		 if(this->atom[i].chain!=this->chain[j].label) continue;
		 else {id=j; mark=true; break;}
		}

	 if(mark==true)
		{
		 // check if this atom belongs to a known residue

		 mark=false; total=this->chain[id].residue.size();

		 for(j=0;j<total;j++)
			{
	 	 	 if(strcmp(atom[i].res_id, chain[id].residue[j].id)) continue;
	 	 	 else if(strcmp(atom[i].residue,chain[id].residue[j].name)) continue;
			 else {mark=true; break;}
			}

		 if(mark==true) continue;
		 else 
			{
			 // add a new residue to this chain

			 tmp_residue.Clear();
			 tmp_residue.valid=1;
			 strcpy(tmp_residue.name,atom[i].residue);
			 strcpy(tmp_residue.id,atom[i].res_id);
			 tmp_residue.chain=chain[id].label;

			 chain[id].residue.push_back(tmp_residue);
			 chain[id].length++;
			}
		}
	 else
		{
		 id=this->chain.size();

		 // add a new chain
		 tmp_chain.Clear();
		 tmp_chain.valid=1;
		 tmp_chain.label=this->atom[i].chain;
		 tmp_chain.residue.clear();

		 this->chain.push_back(tmp_chain);

		 // add the first residue to this chain
		 tmp_residue.Clear();
		 tmp_residue.valid=1;
		 strcpy(tmp_residue.name,atom[i].residue);
		 strcpy(tmp_residue.id,atom[i].res_id);
		 tmp_residue.chain=chain[id].label;

		 chain[id].residue.push_back(tmp_residue);
		 chain[id].length++;
		}
	}

	this->num_chain=this->chain.size();
	for(i=0;i<num_chain;i++) chain[i].length=chain[i].residue.size();

	return;
}

void Protein::Read_ATOM(FILE *fp)
{
	int i,len;
        char line[256],head[80],temp[256];
        char name[80],id[80],residue[80],chain,res_id[80];
	float x,y,z,occupancy,bfactor;
        Atom tmp_atom;

	while(fgets(line,256,fp))
	{
	 // get the heading and erase that section 

	 strcpy(temp,line); temp[6]='\0'; 
	 strcpy(head,""); sscanf(temp,"%s",head);
	 for(i=0;i<6;i++) line[i]=' ';

	 if(!strcmp(head,"END")) break;
	 else if(!strcmp(head,"ATOM")||!strcmp(head,"HETATM"))
		{
		 // get the atom id and erase that section

		 strcpy(temp,line); temp[11]='\0';
		 sscanf(temp,"%s",id);
		 for(i=0;i<11;i++) line[i]=' ';
		 
		 // get the atom name and erase that section 

		 strcpy(temp,line); temp[17]='\0';
		 sscanf(temp,"%s",name);
		 for(i=0;i<17;i++) line[i]=' ';

		 if(strstr(name,"LP")) continue;  // skip lone pairs 

		 // get the residue name and erase that section

		 strcpy(temp,line); temp[20]='\0';
		 sscanf(temp,"%s",residue);
		 for(i=0;i<20;i++) line[i]=' ';

		 // get the chain label and erase that section 

                 chain=line[21]; 
		 for(i=0;i<22;i++) line[i]=' ';

		 // get the residue id and erase that section

		 sscanf(line,"%s", res_id);
		 for(i=0;i<27;i++) line[i]=' ';

		 len=strlen(res_id);
		 if(isdigit(res_id[len-1])) strcat(res_id," ");

		 // get coordinate x and erase that section

		 strcpy(temp,line); temp[38]='\0';
		 sscanf(temp,"%f",&x);
		 for(i=0;i<38;i++) line[i]=' ';

		 // get coordinate y and erase that section

		 strcpy(temp,line); temp[46]='\0';
		 sscanf(temp,"%f",&y);
		 for(i=0;i<46;i++) line[i]=' ';

		 // get coordinate z and erase that section

		 strcpy(temp,line); temp[54]='\0';
		 sscanf(temp,"%f",&z);
		 for(i=0;i<54;i++) line[i]=' ';

		 // get occupancy probability and erase that section
		 // note that not all PDB files have this field

		 len=strlen(line);

		 if(len>=60)  
			{
			 strcpy(temp,line); temp[60]='\0';
			 sscanf(temp,"%f",&occupancy);
			 for(i=0;i<60;i++) line[i]=' ';
			}

		 // get B-factor and erase that section
		 // note that not all PDB files have this field

		 len=strlen(line);

		 if(len>=66)
			{
			 strcpy(temp,line); temp[66]='\0';
			 sscanf(temp,"%f",&bfactor);
			 for(i=0;i<66;i++) line[i]=' ';
			}
			
		 // now summarize

		 sscanf(id,"%d", &tmp_atom.id);
		 sscanf(name,"%s", tmp_atom.name);
		 sscanf(residue,"%s", tmp_atom.residue);
		 tmp_atom.chain=chain;
		 strcpy(tmp_atom.res_id,res_id);
		 tmp_atom.coor[0]=x;
		 tmp_atom.coor[1]=y;
		 tmp_atom.coor[2]=z;
		 tmp_atom.occupancy=occupancy;
		 tmp_atom.bfactor=bfactor;

		 if(!strcmp(head,"ATOM")) tmp_atom.part=1;  // regular atoms
		 else tmp_atom.part=2;  // hetero-atoms

		 tmp_atom.valid=1; tmp_atom.origin=2; // protein atom

		 atom.push_back(tmp_atom);
		}
	 else continue;
	}

	num_atom=atom.size();

	rewind(fp); return;
}

// ****************************************************************************
// This function corrects the common mistakes in PDB files regarding atom names
// and residue names
// Latest update: 12/09/2003
// ****************************************************************************
void Protein::Check_Atom_Type()
{
	int i,j,n,len;
	Atom atm;
	char tmp_name[10];

	for(i=0;i<this->num_atom;i++)
	{
	 atm=this->atom[i];

	// first, correct residue names -------------------------------------

	if(!strcmp(atm.residue,"WAT")) strcpy(atm.residue,"HOH");
	if(!strcmp(atm.residue,"WTR")) strcpy(atm.residue,"HOH");
	if(!strcmp(atm.residue,"CYX")) strcpy(atm.residue,"CYS");
	if(!strcmp(atm.residue,"HID")) strcpy(atm.residue,"HIS");
	if(!strcmp(atm.residue,"HIE")) strcpy(atm.residue,"HIS");
	if(!strcmp(atm.residue,"HIP")) strcpy(atm.residue,"HIS");

	// the following three mistakes are often seen in SYBYL-generated 
	// PDB files

	if(!strcmp(atm.residue,"SO")&&
	   (!strcmp(atm.name,"S")|| !strcmp(atm.name,"O1")||
            !strcmp(atm.name,"O2")|| !strcmp(atm.name,"O3")||
            !strcmp(atm.name,"O4")))
		{
		 strcpy(atm.residue,"SO4");
		}

	if(!strcmp(atm.residue,"PO")&&
	   (!strcmp(atm.name,"P")|| !strcmp(atm.name,"O1")||
            !strcmp(atm.name,"O2")|| !strcmp(atm.name,"O3")||
            !strcmp(atm.name,"O4")))
		{
		 strcpy(atm.residue,"PO4");
		}

	if(!strcmp(atm.residue,"CA")&&!strcmp(atm.name,"O"))
		{
		 strcpy(atm.residue,"HOH");
		}

	// now correct atom names -------------------------------------------

 	if(isdigit(atm.name[0]))  
		{
		 len=strlen(atm.name); n=0;
		 for(j=1;j<len;j++)
			{
			 tmp_name[n]=atm.name[j]; n++;
			}
		 tmp_name[n]=atm.name[0]; n++; tmp_name[n]='\0';
		 strcpy(atm.name,tmp_name);
		}	

	if(!strcmp(atm.name,"OCT")) strcpy(atm.name,"OXT");  // all residue
	if(!strcmp(atm.name,"HN")) strcpy(atm.name,"H");     // all residue

	if(!strcmp(atm.residue,"ACE"))
                {
                 if(!strcmp(atm.name,"CH3")) strcpy(atm.name,"CA");
                }
	else if(!strcmp(atm.residue,"ASN"))
		{
		 if(!strcmp(atm.name,"AD1")) strcpy(atm.name,"OD1");
		 else if(!strcmp(atm.name,"AD2")) strcpy(atm.name,"ND2");
		}
	else if(!strcmp(atm.residue,"GLN"))
                {
                 if(!strcmp(atm.name,"AE1")) strcpy(atm.name,"OE1");
                 else if(!strcmp(atm.name,"AE2")) strcpy(atm.name,"NE2");
                }
	else if(!strcmp(atm.residue,"LEU"))
		{
		 if(!strcmp(atm.name,"CD")) strcpy(atm.name,"CD1");
		 else if(!strcmp(atm.name,"CE")) strcpy(atm.name,"CD2");
		}

	 this->atom[i]=atm;
	}

	return;
}

// ****************************************************************************
// this function reads the CONECT information from the input PDB file
// latest update: 10/28/2003
// ****************************************************************************
void Protein::Read_CONECT(FILE *fp)
{
	char line[256],head[256],temp[256];
	int i,j,id,atom_id,len,count;
	bool mark;

	while(fgets(line,256,fp))
	{
	 strcpy(temp,line); temp[6]='\0'; 
	 strcpy(head,""); sscanf(temp,"%s",head);

	 if(strcmp(head,"CONECT")) continue;

	 for(i=0;i<6;i++) line[i]=' '; // remove the header

	 // read the atom id first

	 strcpy(temp,line); temp[11]='\0'; 
	 sscanf(temp,"%d",&atom_id);
	 for(i=0;i<11;i++) line[i]=' ';

	 mark=false; 

	 for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid==0) continue;
		 else if(this->atom[i].part==1) continue;
		 else if(this->atom[i].id!=atom_id) continue;
		 else {id=i; mark=true; break;}
		}

	 if(mark==false) continue;	// no such atom in hetero atom list

	 for(i=1;i<=MAX_ATOM_NEIB;i++) 	
		{
		 len=11+i*5; strcpy(temp,line); temp[len]='\0';
		 if(!isdigit(temp[len-1])) break;  // no more numbers  
		 else
			{
			 sscanf(temp,"%d",&atom_id);
			 for(j=0;j<len;j++) line[j]=' ';

			 count=this->atom[id].num_neib;
			 this->atom[id].neib[count]=atom_id;
			 this->atom[id].num_neib++;
			}
		}

/*
	 printf("CONECT read: %d ", this->atom[id].id);
	 for(i=0;i<this->atom[id].num_neib;i++)
		printf("%d ", this->atom[id].neib[i]);
	 printf("\n");
*/
	}

	rewind(fp); return;
}

// ****************************************************************************
// detect the connection table of this protein molecule, including ATOMs,
// metal ions, SO4s, and PO4s
// latest update: 11/04/2003
// ****************************************************************************
void Protein::Detect_Connections()
{
	extern ForceField *ff;
	int i,j,k,count,id;
	bool mark;
	float cutoff,d;

	// build the connections for regular atoms first

	cutoff=2.00;  // covalent bond distance cutoff

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid==0) continue;
	 else if(this->atom[i].part!=1) continue;
	 else if(this->atom[i].num_neib>0) continue; // already done 

	 count=0;

	 // find internal connections inside the same residue first

	 for(j=0;j<this->num_atom;j++)
		{
		 if(i==j) continue;
		 else if(this->atom[j].valid==0) continue;
		 else if(this->atom[j].part!=1) continue;
		 else if(atom[i].chain!=atom[j].chain) continue;
		 else if(strcmp(atom[i].residue,atom[j].residue)) continue;
		 else if(strcmp(atom[i].res_id,atom[j].res_id)) continue;

		 // there is no bond between hydrogen atoms

		 if(!strcmp(atom[i].type,"H")&&
		    !strcmp(atom[j].type,"H")) continue;

		 // now check whether atom i and atom j are covalently bound

		 if(ff->Patom_Connection_Test(atom[i],atom[j])==true)
			{
			 if(count>=MAX_ATOM_NEIB) continue;  // already full
			 else
				{
			 	 this->atom[i].neib[count]=j+1; count++;
				}
			}
		 else continue;
		}

	 // now detect peptide amide bonds

	 if(!strcmp(atom[i].name,"N")&&!strcmp(atom[i].type,"N.am"))
		{
		 for(j=0;j<this->num_atom;j++)
			{
			 if(this->atom[j].valid==0) continue;
			 else if(this->atom[j].part!=1) continue;
			 else if(atom[i].chain!=atom[j].chain) continue;
			 else if(!strcmp(atom[i].res_id,atom[j].res_id)) continue;
			 else if(strcmp(atom[j].type,"C.2")) continue;
			 else if(strcmp(atom[j].name,"C")) continue;

			 d=Distance(atom[i].coor,atom[j].coor);

			 if(d>cutoff) continue;
			 else if(count>=MAX_ATOM_NEIB) continue; 
			 else
				{
			 	 this->atom[i].neib[count]=j+1; 
				 count++;
				}
			}
		}
	 else if(!strcmp(atom[i].name,"C")&&!strcmp(atom[i].type,"C.2"))
		{
		 for(j=0;j<this->num_atom;j++)
                        {
                         if(this->atom[j].valid==0) continue;
                         else if(this->atom[j].part!=1) continue;
                         else if(atom[i].chain!=atom[j].chain) continue;
                         else if(!strcmp(atom[i].res_id,atom[j].res_id)) continue;
                         else if(strcmp(atom[j].type,"N.am")) continue;
                         else if(strcmp(atom[j].name,"N")) continue;

                         d=Distance(atom[i].coor,atom[j].coor);

                         if(d>cutoff) continue;
			 else if(count>=MAX_ATOM_NEIB) continue; 
                         else
                                {
                                 this->atom[i].neib[count]=j+1; 
                                 count++;
                                }
                        }
		}

	 this->atom[i].num_neib=count;
	}

	// now find the patoms bound to the metal ions

	cutoff=3.00;  // M-bond distance cutoff

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid==0) continue;
	 else if(this->atom[i].part==1) continue;
	 else if(strcmp(this->atom[i].xtype,"M+")) continue;
	 else if(this->atom[i].num_neib>0) continue;  // already done

	 count=0;

	 for(j=0;j<this->num_atom;j++)
		{
		 if(this->atom[j].valid==0) continue;
		 else if(this->atom[j].part!=1) continue;
		 else if(strcmp(this->atom[j].hb,"A")&&
                         strcmp(this->atom[j].hb,"DA")) continue;

		 d=Distance(this->atom[i].coor,this->atom[j].coor);

		 if(d>cutoff) continue;
		 else if(count>=MAX_ATOM_NEIB) continue;
		 else 
			{
			 this->atom[i].neib[count]=j+1;  // note this
			 count++;
			 // here we do not add connections to atom[j]
			}
		}

	 this->atom[i].num_neib=count;
	}

	// now build the connections for SO4s and PO4s

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid==0) continue;
	 else if(this->atom[i].part==1) continue;
	 else if(strcmp(this->atom[i].residue,"SO4")&&
		 strcmp(this->atom[i].residue,"PO4")) continue;
	 else if(strcmp(this->atom[i].name,"S")&&
		 strcmp(this->atom[i].name,"P")) continue;

	 if(this->atom[i].num_neib>0) continue;  // done 

	 // now check the satellite atoms for this P or S atom

	 for(j=0;j<this->num_atom;j++)
		{
		 if(this->atom[j].valid==0) continue;
		 else if(this->atom[j].part==1) continue;
		 else if(strcmp(this->atom[j].residue,"SO4")&&
                         strcmp(this->atom[j].residue,"PO4")) continue;
		 else if(strcmp(this->atom[j].res_id,
				this->atom[i].res_id)) continue;
		 else if(strstr(this->atom[j].name,"S")) continue;
		 else if(strstr(this->atom[j].name,"P")) continue;

		 count=this->atom[i].num_neib; mark=false;

		 for(k=0;k<count;k++)
			{
			 if(this->atom[i].neib[k]==(j+1))
				{
				 mark=true; break;
				}
			 else continue;
			}

		 if(mark==false)
			{
			 this->atom[i].neib[count]=j+1;  // note this
			 this->atom[i].num_neib++;
			}

		 count=this->atom[j].num_neib; mark=false;

		 for(k=0;k<count;k++)
			{
			 if(this->atom[j].neib[k]==(i+1))
				{
				 mark=true; break;
				}
			 else continue;
			}

		 if(mark==false)
			{
			 this->atom[j].neib[count]=i+1;  // note this
			 this->atom[j].num_neib++;
			}
		}
	}

/*
	// now list what we get

	for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid<=0) continue;

		 printf("ATOM %d %s %s %s ", 
			this->atom[i].id, 
			this->atom[i].name,
			this->atom[i].residue,
			this->atom[i].res_id);

		 for(j=0;j<this->atom[i].num_neib;j++)
			{
			 id=this->atom[i].neib[j];
			 printf("%d ", id);
			}

		 printf("\n");
		}
*/

	// now determine number of heavy atoms for each valid atom

	for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid<=0) continue;

		 count=0;

		 for(j=0;j<this->atom[i].num_neib;j++)
			{
			 id=this->atom[i].neib[j];
			 if(!strcmp(this->atom[id-1].type,"H")) continue;
			 else count++;
			}

		 this->atom[i].num_nonh=count;
		}

	return;
}

// ****************************************************************************
// Write out a standard PDB file, which includes HEADER, SEQRES, ATOM, HETATM,
// and CONECT sections. 
// latest update: 11/04/2003
// ****************************************************************************
void Protein::Write_Out_PDB(char *filename)
{
	FILE *fp;
	int i,j,k,len,temp,count;
	bool mark;

	if((fp=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	// convert the name into upper case

	len=strlen(this->name);
	for(i=0;i<len;i++) this->name[i]=toupper(this->name[i]);

	// now start to output

	fprintf(fp,"HEADER    %s\n", this->name);
	fprintf(fp,"COMPND    %s\n", this->name);
	fprintf(fp,"REMARK    GENERATED BY X-TOOL on %s\n", Get_Time());

	// write out the SEQRES section

	for(i=0;i<this->num_chain;i++)
	{
	 count=(int)ceil(this->chain[i].length/13.0);

	 for(j=1;j<=count;j++)
		{
		 fprintf(fp,"SEQRES %3d %c %4d  ", 
			j, this->chain[i].label, this->chain[i].length);

		 for(k=1;k<=13;k++)
			{
			 temp=(j-1)*13+k;
			 if(temp>this->chain[i].length) break;
			 else fprintf(fp,"%3s ", this->chain[i].residue[temp-1].name);
			}
		
		 fprintf(fp,"\n");
		}
	}

	// write out the ATOM section 

	for(i=0;i<this->num_chain;i++)
	{
	 for(j=0;j<this->num_atom;j++)
		{
		 if(this->atom[j].valid<=0) continue;
		 else if(this->atom[j].part!=1) continue;
		 else if(this->atom[j].chain!=this->chain[i].label) continue;

		 // the atom section

		 fprintf(fp,"ATOM  %5d ",this->atom[j].id);

		 if(strlen(this->atom[j].name)>=4) 
			{ 
			 fprintf(fp,"%-4s ",this->atom[j].name);
			}
		 else
			{
			 fprintf(fp," %-3s ",this->atom[j].name);
			}

		 // the residue section

		 fprintf(fp,"%3s %c%5s   ",
			 this->atom[j].residue,
			 this->atom[j].chain,
			 this->atom[j].res_id);

		 // the coordinate section

		 fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
			 this->atom[j].coor[0],
			 this->atom[j].coor[1],
			 this->atom[j].coor[2],
			 this->atom[j].occupancy,
			 this->atom[j].bfactor);
		}

	 fprintf(fp,"TER\n");
	}

	// now write out the HETATM section, including cofactor in any 

	for(j=0;j<this->num_atom;j++)
		{
		 if(this->atom[j].valid<=0) continue;
		 else if(this->atom[j].part<=1) continue;

		 // the atom section

		 fprintf(fp,"HETATM%5d ",this->atom[j].id);

		 if(strlen(this->atom[j].name)>4) this->atom[j].name[4]='\0';
		 
		 if(!strcmp(this->atom[j].hb,"M")||
                    strlen(this->atom[j].name)>=4)
			{
			 fprintf(fp,"%-4s ",this->atom[j].name);
			}
		 else
			{
			 fprintf(fp," %-3s ",this->atom[j].name);
			}

		 // the residue section

		 fprintf(fp,"%3s %c%5s   ", 
			 this->atom[j].residue, 
			 ' ',
			 this->atom[j].res_id);

		 // the coordinate section

		 fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
			 this->atom[j].coor[0],
			 this->atom[j].coor[1],
			 this->atom[j].coor[2],
			 this->atom[j].occupancy, 
			 this->atom[j].bfactor);
		}

	// write out the CONECT section, note that XTOOL will only consider
        // the cofactor, metal ions, SO4 and PO4 for this purpose. 

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid==0) continue;
	 else if(this->atom[i].part==1) continue;
	 else if(this->atom[i].num_neib==0) continue;

	 // check eligibility

	 if(!strcmp(this->atom[i].residue,"COF")) mark=true;
	 else if(!strcmp(this->atom[i].residue,"SO4")) mark=true;
	 else if(!strcmp(this->atom[i].residue,"PO4")) mark=true;
	 else if(!strcmp(this->atom[i].xtype,"M+")) mark=true;
	 else mark=false;

	 if(mark==false) continue;
	 else
		{
	 	 fprintf(fp,"CONECT%5d",this->atom[i].id);
	 	 for(j=0;j<this->atom[i].num_neib;j++) 
			{
		 	 fprintf(fp,"%5d",this->atom[i].neib[j]);
			}
	 	 fprintf(fp,"\n");
		}
	}

	// the end of the PDB file

	fprintf(fp,"END\n");

	fclose(fp); return;
}

// ***************************************************************************
// prepare a PDB file that meets X-Tool's format, which only requires 
// polar hydrogen atoms. Metal ions, PO4, SO4, CL, and water molecules
// are included in the HETATM section
// latest update: 08/02/2003
// ***************************************************************************
void Protein::Write_Out_PDB_XTool(char *filename)
{
	int i;

	// filter out unwanted components

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(strstr(this->atom[i].name,"LP")) this->atom[i].valid=0;
	 else if((this->atom[i].part==1)&&
                 (!strcmp(this->atom[i].xtype,"H"))) this->atom[i].valid=0;
	 else if((this->atom[i].part>1)&&
                 (!strcmp(this->atom[i].xtype,"H"))) this->atom[i].valid=0;
	 else if((this->atom[i].part>1)&&
                 (!strcmp(this->atom[i].xtype,"H.hb"))) this->atom[i].valid=0;
	 else continue;
	}

	this->Rearrange_IDs();

	// output the PDB file

	this->Write_Out_PDB(filename);

	return;
}

// ***************************************************************************
// prepare a PDB file that meets DrugScore's requirements, which does not
// include any hudrogen atom and water molecule
// latest update: 08/02/2003
// ***************************************************************************
void Protein::Write_Out_PDB_DrugScore(char *filename)
{
	FILE *fin,*fout;
	int i,len,count;
	char line[81];

	// filter out unwanted components

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].type,"H")) this->atom[i].valid=0;
	 else if(!strcmp(this->atom[i].residue,"HOH")) this->atom[i].valid=0;
	 else if(!strcmp(this->atom[i].residue,"WAT")) this->atom[i].valid=0;
	 else if(!strcmp(this->atom[i].residue,"PO4")) this->atom[i].valid=0;
	 else if(!strcmp(this->atom[i].residue,"SO4")) this->atom[i].valid=0;
	 else if(!strcmp(this->atom[i].residue,"CL")) this->atom[i].valid=0;
	 else if(strstr(this->atom[i].name,"LP")) this->atom[i].valid=0;
	 else continue;
	}

	this->Rearrange_IDs();

	// save the result in a temparary file

	Write_Out_PDB("temp_drugscore.pdb");

	// now write out the desired PDB file

	if((fin=fopen("temp_drugscore.pdb","r"))==NULL) Open_File_Error("temp_drugscore.pdb");
	if((fout=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	count=1;

	while(fgets(line,80,fin)) 
	{
	 len=strlen(line);

	 for(i=80;i>=72;i--) line[i]=' ';
	 if(len<=72) {for(i=len-1;i<=71;i++) line[i]=' ';}
	 line[72]='\0';

	 fprintf(fout,"%72s%4s%4d\n",line,"XXXX",count);
	 count++;
	}

	fclose(fin); fclose(fout); system("rm temp_drugscore.pdb");

	return;
}

// ***************************************************************************
// the flag controls if the HB-roots are computed or not
// latest update: 03/10/2003
// ***************************************************************************
void Protein::Value_Atom(int flag)
{
	extern ForceField *ff;
	int i;

	// note that within this subroutine the protein is considered to be
	// all-atom model. this is because this subroutine is shared by
	// all of the functions related to protein

	// correct atom names and residue names first.

	this->Check_Atom_Type();

	// assign the parameters for all atoms

	for(i=0;i<num_atom;i++) ff->Assign_Patom_Parameters(atom[i]);

	// build connection tables, must assign parameters first

	this->Detect_Connections();

	// calculate H-bond root

	if(flag!=0) Calculate_HB_Root();

	return;
}

// ***********************************************************************
// for each HB atom, calculate the position of its 'root', and also locate 
// its neighboring atoms.
// latest update: 11/06/2003
// ***********************************************************************
void Protein::Calculate_HB_Root()
{
	int i,j,id,count;
	float tmpx,tmpy,tmpz;

	// (1) we do calculate root for polar hydrogens 
	// (2) we do NOT calculate root for non-HB atoms
	// (3) we do NOT calculate root for metal ions
	// (4) we do NOT calculate root for waters   
	// (5) we do calculate root for the oxygen atoms on SO4 and PO4

	for(i=0;i<num_atom;i++)
	{
	 if(atom[i].valid<=0) continue;
	 else if(!strcmp(atom[i].xtype,"H")) continue;
	 else if(!strcmp(atom[i].xtype,"O.w")) continue;
	 else if(!strcmp(atom[i].hb,"N")) continue;
	 else if(!strcmp(atom[i].hb,"H")) continue;
	 else if(!strcmp(atom[i].hb,"P")) continue;
	 else if(!strcmp(atom[i].hb,"M")) continue;

	 if(atom[i].num_neib==0) continue; 

	 tmpx=tmpy=tmpz=0.000; count=0;

	 for(j=0;j<atom[i].num_neib;j++)
		{
		 id=atom[i].neib[j]-1;
		 if(!strcmp(atom[id].type,"H")) continue;
		 else
			{
		 	 tmpx+=(atom[id].coor[0]);
		 	 tmpy+=(atom[id].coor[1]);
		 	 tmpz+=(atom[id].coor[2]);
			 count++;
			}
		}

	 if(count==0) continue;

	 atom[i].root[0]=tmpx/count;
	 atom[i].root[1]=tmpy/count;
	 atom[i].root[2]=tmpz/count;
	}

	return;
}


// *****************************************************************************
// find all neigboring residues around a given ligand within a distance cutoff
// including metal ions and water molecules
// lastest update 01/08/2003
// *****************************************************************************
void Protein::Define_Pocket(const Ligand *ligand, float cutoff)
{
	extern ForceField *ff;
	int i,j,mark,count;
	float d;
	int num_pocket_res; 
        Residue *pocket_res=NULL;

	num_pocket_res=0;

	pocket_res=new Residue[300];
        if(pocket_res==NULL) Memory_Allocation_Error();

	// find pocket atoms on the protein and label their valid as '2',
	// including the waters close to the ligand 

	for(i=0;i<num_atom;i++)
	{
	 if(atom[i].valid<=0) continue;

	 // filter out all of the non-polar hydrogen atoms
	 // since X-Score uses unit-atom model

	 if(!strcmp(atom[i].xtype,"H")&&strcmp(atom[i].residue,"COF")) 
		{
		 atom[i].valid=0; continue;
		}

	 // check if this atom is close to the ligand

	 mark=FALSE;

	 for(j=0;j<ligand->num_atom;j++)
		{
		 if(ligand->atom[j].valid<=0) continue;
		 else if(!strcmp(ligand->atom[j].type,"H")) continue;
		 else
			{
			 d=Distance(atom[i].coor,ligand->atom[j].coor);
			 if(d>cutoff) continue;
                         else {mark=TRUE; break;}
			}
		}

	 if(mark==FALSE) continue;	// not a pocket atom

	 atom[i].valid=2;

	 // if this atom is a water or a metal ion, do not check its residue

	 if(!strcmp(atom[i].xtype,"O.w")) continue;
	 if(!strcmp(atom[i].xtype,"M+")) continue;

	 // for regular atoms, check if it has been included in a 
	 // known pocket residue; if not, add that newly found pocket residue	

	 mark=FALSE;

	 for(j=0;j<num_pocket_res;j++)
       		{
       		 if(!strcmp(atom[i].residue,pocket_res[j].name)&& 
		    !strcmp(atom[i].res_id,pocket_res[j].id)&&
		    (atom[i].chain==pocket_res[j].chain)) {mark=TRUE;break;}
		 else continue;
		}

	 if(mark==FALSE)  // new pocket residue found
		{
		 strcpy(pocket_res[num_pocket_res].name,atom[i].residue);
		 strcpy(pocket_res[num_pocket_res].id,atom[i].res_id);
		 pocket_res[num_pocket_res].chain=atom[i].chain;
		 pocket_res[num_pocket_res].valid=1;
		 num_pocket_res++;
		}
	}

	// now check if the binding pocket is well defined

	count=0;

	for(i=0;i<num_pocket_res;i++)
	{
	 if(!strcmp(pocket_res[i].name,"HET")) continue;
	 else if(!strcmp(pocket_res[i].name,"COF")) continue;
	 else if(!strcmp(pocket_res[i].name,"WAT")) continue;
	 else if(!strcmp(pocket_res[i].name,"HOH")) continue;
	 else count++;
	}

	if(count<3) // there should be at least 3 residues as pocket!
	{
	 puts("Error: cannot find binding pocket residues on the protein.");
	 puts("Probably the ligand has not been docked with the protein.");
	 exit(1);
	}

	// define all the left atoms in pocket residues as pocket atoms 

	for(i=0;i<num_atom;i++)
	{
	 if(atom[i].valid<=0) continue;
	 else if(atom[i].valid==2) continue;
	 else if(!strcmp(atom[i].type,"O.w")) continue;
	 else if(!strcmp(atom[i].hb,"M")) continue;

 	 for(j=0;j<num_pocket_res;j++)
		{
		 if(strcmp(atom[i].residue,pocket_res[j].name)) continue;
		 else if(strcmp(atom[i].res_id,pocket_res[j].id)) continue;
	 	 else if(atom[i].chain!=pocket_res[j].chain) continue;
	 	 else {atom[i].valid=2; break;}
		}
	}

	// now detect all of the aromatic rings within pocket residues
	// this information is need for later scoring purposes

	Ring tmp_ring;

	this->ring.clear(); this->num_ring=0;

	for(i=0;i<num_pocket_res;i++)
	{
	 if(strcmp(pocket_res[i].name,"PHE")&&
            strcmp(pocket_res[i].name,"TYR")&&
            strcmp(pocket_res[i].name,"HIS")&&
            strcmp(pocket_res[i].name,"TRP")) continue;

	 tmp_ring.Clear();

	 for(j=0;j<num_atom;j++)
		{
		 if(atom[j].valid!=2) continue;
		 else if(atom[j].part!=1) continue;
		 else if(atom[j].ring!=2) continue;
		 else if(atom[j].chain!=pocket_res[i].chain) continue;
		 else if(strcmp(atom[j].residue,pocket_res[i].name)) continue;
		 else if(strcmp(atom[j].res_id,pocket_res[i].id)) continue;

		 tmp_ring.atom_id.push_back(j+1);   
		 tmp_ring.centroid[0]+=atom[j].coor[0];
		 tmp_ring.centroid[1]+=atom[j].coor[1];
		 tmp_ring.centroid[2]+=atom[j].coor[2];
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

	if(pocket_res) if(pocket_res) delete [] pocket_res; return;
}

// *****************************************************************************
// find all neigboring residues around a given origin and a radius 
// including metal ions and water molecules
// lastest update 09/25/2003
// *****************************************************************************
void Protein::Define_Pocket(float origin[], float radius)
{
	// need to be finished 
	return;
}

// ***************************************************************************
// determine which water to keep (valid=2) and which to disregard (valid=1).
// latest update 03/04/03
// ***************************************************************************
void Protein::Define_Water(const Ligand *ligand, bool w_flag)
{
        int i,j,mark;
        float d;

	if(w_flag==false)	// water molecules are not considered at all
	{
	 for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid<2) continue;
		 else if(strcmp(this->atom[i].type,"O.w")) continue;
		 else this->atom[i].valid=1;
		}
	 return;
	}

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<2) continue;
	 else if(strcmp(this->atom[i].type,"O.w")) continue;

	 // first, this water should be close enough to the ligand
	 // only keep the first solvation shell

	 mark=false;

	 for(j=0;j<ligand->num_atom;j++)
		{
		 if(ligand->atom[j].valid<=0) continue;
		 else if(!strcmp(ligand->atom[j].type,"H")) continue;

		 d=Distance(this->atom[i].coor,ligand->atom[j].coor);

		 if(d<(ligand->atom[j].R+2*WATER_R))
			{
			 mark=true; break;
			}
		 else continue;
		}

	 if(mark==false) {this->atom[i].valid=1; continue;}

	 // second, buried status check 

	 mark=Check_Buried_Ratio(this->atom[i].coor);

	 if(mark==true) this->atom[i].valid=2;
	 else this->atom[i].valid=1;
	}

	return;
}

void Protein::Write_Pocket_PDB(char *filename) const
{
	int i,count;
        FILE *fp;

        if((fp=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	fprintf(fp,"HEADER    BINDING POCKET\n");
	fprintf(fp,"COMPND    BINDING POCKET\n");
	fprintf(fp,"REMARK    CREATED BY XTOOL  %s\n", Get_Time());

	count=1;

	for(i=0;i<num_atom;i++)
		{
		 if(atom[i].valid!=2) continue;
		 // note that the format of metal is slightly different
		 // from other hetero-atoms
		 else if(!strcmp(atom[i].hb,"M"))
			{
		 	 fprintf(fp,"HETATM");
			 fprintf(fp,"%5d ", count);
                 	 fprintf(fp,"%-4s ", atom[i].name);
                 	 fprintf(fp,"%3s ", atom[i].residue);
                 	 fprintf(fp,"%c", ' ');
                 	 fprintf(fp,"%5s", atom[i].res_id);
			 fprintf(fp,"   ");
                 	 fprintf(fp,"%8.3f", atom[i].coor[0]);
                 	 fprintf(fp,"%8.3f", atom[i].coor[1]);
                 	 fprintf(fp,"%8.3f", atom[i].coor[2]);
                 	 fprintf(fp,"%6.2f", 1.00);
			 fprintf(fp,"%6.2f", atom[i].q);
			 fprintf(fp,"\n");
			 count++;
			}
		 else if(atom[i].part>1)       // hetero atoms
                        {
			 fprintf(fp,"HETATM");
			 fprintf(fp,"%5d  ", count);
                         fprintf(fp,"%-4s", atom[i].name);
                         fprintf(fp,"%3s ", atom[i].residue);
                         fprintf(fp,"%c", ' ');
			 fprintf(fp,"%5s", atom[i].res_id);
                         fprintf(fp,"   ");
                         fprintf(fp,"%8.3f", atom[i].coor[0]);
                         fprintf(fp,"%8.3f", atom[i].coor[1]);
                         fprintf(fp,"%8.3f", atom[i].coor[2]);
                         fprintf(fp,"%6.2f", 1.00);
                         fprintf(fp,"%6.2f", atom[i].q);
                         fprintf(fp,"\n");
                         count++;
                        }
		 else
			{
			 fprintf(fp,"ATOM  ");
			 fprintf(fp,"%5d ", count);

			 if(strlen(atom[i].name)>=4)
				{
				  fprintf(fp,"%-4s ", atom[i].name);
				}
			 else
				{
				 fprintf(fp," %-3s ", atom[i].name);
				}
	
			 fprintf(fp,"%3s ", atom[i].residue);
                         fprintf(fp,"%c ", ' ');
                         fprintf(fp,"%4s", atom[i].res_id);
                         fprintf(fp,"   ");
                         fprintf(fp,"%8.3f", atom[i].coor[0]);
                         fprintf(fp,"%8.3f", atom[i].coor[1]);
                         fprintf(fp,"%8.3f", atom[i].coor[2]);
                         fprintf(fp,"%6.2f", 1.00);
                         fprintf(fp,"%6.2f", atom[i].q);
                         fprintf(fp,"\n");
                         count++;
			}
		}

	fprintf(fp,"END\n");

	fclose(fp);

	return;
}

// ***************************************************************************
// output the parameters of the already defined binding pocket
// latest update: 12/04/2003
// ***************************************************************************
void Protein::Write_Pocket_XTOOL(char *filename) const
{
	FILE *fp;
	int i,j,count;
	int *id=NULL;

	id=new int[this->num_atom];
	if(id==NULL) Memory_Allocation_Error();

	if((fp=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	fprintf(fp,"#\n");
	fprintf(fp,"# This file stores the binding pocket atoms.\n");
	fprintf(fp,"# Creation time: %s\n", Get_Time());
	fprintf(fp,"#\n");

	fprintf(fp,"# column 1: atom ID\n");
	fprintf(fp,"# column 2: atom name\n");
	fprintf(fp,"# column 3: residue name\n");
	fprintf(fp,"# column 4: residue ID\n");
	fprintf(fp,"# column 5: X coordinate\n");
	fprintf(fp,"# column 6: Y coordinate\n");
	fprintf(fp,"# column 7: Z coordinate\n");
	fprintf(fp,"# column 8: SYBYL type\n");
	fprintf(fp,"# column 9: vdw radius\n");
	fprintf(fp,"# column 10: vdw weight factor\n");
	fprintf(fp,"# column 11: atomic charge\n");
	fprintf(fp,"# column 12: X-TOOL type\n");
	fprintf(fp,"# column 13: H-bonding character\n");
	fprintf(fp,"# column 14: atomic logP\n");
	fprintf(fp,"# column 15: ring classfication\n");
	fprintf(fp,"# column 16: H-bond root X coordinate\n"); 
	fprintf(fp,"# column 17: H-bond root Y coordinate\n"); 
	fprintf(fp,"# column 18: H-bond root Z coordinate\n");
	fprintf(fp,"# column 19: primary neighbor ID\n");
	fprintf(fp,"#\n");

	// re-assign new ids for output atoms

	count=0;

        for(i=0;i<num_atom;i++)
                {
                 if(atom[i].valid<2) {id[i]=0;}
                 else {count++; id[i]=count;}
                }

	// now write out the atom section

	fprintf(fp,"<START>\n");
	fprintf(fp,"<ATOM> %d\n", count);

	for(i=0;i<num_atom;i++)
		{
		 if(atom[i].valid<2) continue;

		 fprintf(fp,"%-4d ", id[i]);
		 fprintf(fp,"%-4s ", atom[i].name);
		 fprintf(fp,"%-4s ", atom[i].residue);
		 fprintf(fp,"%-4s ", atom[i].res_id);
		 fprintf(fp,"%8.3f ", atom[i].coor[0]);
		 fprintf(fp,"%8.3f ", atom[i].coor[1]);
		 fprintf(fp,"%8.3f ", atom[i].coor[2]);
		 fprintf(fp,"%-5s ", atom[i].type);
		 fprintf(fp,"%4.2f ", atom[i].R);
		 fprintf(fp,"%4.2f ", atom[i].eps);
		 fprintf(fp,"%5.2f ", atom[i].q);
		 fprintf(fp,"%-7s ", atom[i].xtype);
		 fprintf(fp,"%-2s ", atom[i].hb);
		 fprintf(fp,"%5.2f ", atom[i].logp);
		 fprintf(fp,"%1d ", atom[i].ring);
		 fprintf(fp,"%8.3f ", atom[i].root[0]);
                 fprintf(fp,"%8.3f ", atom[i].root[1]);
                 fprintf(fp,"%8.3f ", atom[i].root[2]);

		 // note that neighboring atom is only used for H.hb 

		 if(strcmp(atom[i].hb,"DH")||atom[i].neib[0]==0)
			{
		  	 fprintf(fp,"%-d ", 0);
			}
		 else 
			{
			 fprintf(fp,"%-d", id[atom[i].neib[0]-1]);
			}

		 fprintf(fp,"\n");
		}

	// now write out the ring section

	fprintf(fp,"<RING> %d\n", this->num_ring);

	for(i=0;i<this->num_ring;i++)
		{
		 fprintf(fp, "%-2d ", i+1);
		 fprintf(fp, "%1d ", this->ring[i].type);
		 fprintf(fp, "%-2d ", this->ring[i].num_member); 

		 for(j=0;j<this->ring[i].num_member;j++)
			{
			 fprintf(fp, "%4d ", id[ring[i].atom_id[j]-1]);
			}

		 fprintf(fp,"\n");
		}

	// now finish it

	fprintf(fp,"<END>\n");

	fclose(fp); if(id) delete [] id; return;
}

// ***************************************************************************
// read the parameters of pre-defined binding pocket
// latest update: 12/04/2003
// ***************************************************************************
void Protein::Read_Pocket_XTOOL(char *filename)
{
	FILE *fp;
	int i,j,k,n,id;
        char line[256],head[256],format[256];
	Atom tmp_atom;
	Ring tmp_ring;

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        do
        {
         if(fgets(line,256,fp)==NULL) Read_File_Error(filename);
         else sscanf(line,"%s",head);
        } while(strcasecmp(head,"<START>"));

	// now read the <ATOM> setction 

	do
        {
         if(fgets(line,256,fp)==NULL) Read_File_Error(filename);
         else sscanf(line,"%s",head);
        } while(strcasecmp(head,"<ATOM>"));

	sscanf(line,"%*s%d",&this->num_atom);

	this->atom.clear();

	for(i=0;i<this->num_atom;i++)
		{
		 fgets(line,256,fp); tmp_atom.Clear();

		 sscanf(line,"%d%s%s%s%f%f%f%s%f%f%f%s%s%f%hd%f%f%f%d", 
			&tmp_atom.id, 
			 tmp_atom.name,
			 tmp_atom.residue,
			 tmp_atom.res_id,
			&tmp_atom.coor[0], 
			&tmp_atom.coor[1], 
			&tmp_atom.coor[2], 
			 tmp_atom.type,
			&tmp_atom.R,
			&tmp_atom.eps,
			&tmp_atom.q, 
			 tmp_atom.xtype, 
			 tmp_atom.hb,
			&tmp_atom.logp, 
			&tmp_atom.ring,
			&tmp_atom.root[0], 
			&tmp_atom.root[1], 
			&tmp_atom.root[2],
			&tmp_atom.neib[0]);

		 tmp_atom.valid=2; tmp_atom.origin=2;

		 this->atom.push_back(tmp_atom);
		}

	if(this->num_atom!=this->atom.size())
	{
	 printf("Warning: number of atoms does not match in %s\n", filename);
	}

	// now read the <RING> section

	do
        {
         if(fgets(line,256,fp)==NULL) Read_File_Error(filename);
         else sscanf(line,"%s",head);
        } while(strcasecmp(head,"<RING>"));

        sscanf(line,"%*s%d",&this->num_ring);

	this->ring.clear();

	for(i=0;i<this->num_ring;i++)
		{
		 fgets(line,256,fp); tmp_ring.Clear();
		 sscanf(line,"%*d%d%d",&tmp_ring.type,&tmp_ring.num_member);

		 for(j=0;j<tmp_ring.num_member;j++)
			{
			 n=3+j; strcpy(format,"");
			 for(k=1;k<=n;k++) strcat(format,"%*d");
			 strcat(format,"%d");
			 
			 sscanf(line,format,&id);

			 tmp_ring.atom_id.push_back(id);
			 tmp_ring.centroid[0]+=this->atom[id-1].coor[0];
			 tmp_ring.centroid[1]+=this->atom[id-1].coor[1];
			 tmp_ring.centroid[2]+=this->atom[id-1].coor[2];
			}

		 tmp_ring.centroid[0]/=tmp_ring.num_member;
		 tmp_ring.centroid[1]/=tmp_ring.num_member;
		 tmp_ring.centroid[2]/=tmp_ring.num_member;
		 tmp_ring.valid=1;
		 this->ring.push_back(tmp_ring);
		}

	if(this->num_ring!=this->ring.size())
        {
         printf("Warning: number of rings does not match in %s\n", filename);
        }

	fclose(fp); return;
}

void Protein::Generate_Volume_Dots(float spacing, float probe_r)
// if probe_r=0, then it represents the VDW volume;
// if probe_r=1.40, then it represents the solvent-inaccessible volume;
// note that it is a unit-atom model
{
	int i,j,num,atom_id;
	float max_x,min_x,max_y,min_y,max_z,min_z;
	float coor[3],d,margin,dmin,unit;
	bool mark;
	Dot tmp_dot;

	vol_dot.clear(); // clear all the existing volume points, if any

	max_x=max_y=max_z=-9999.000;
	min_x=min_y=min_z= 9999.000;

	for(i=0;i<this->num_atom;i++)
	{
 	 if(this->atom[i].valid<2) continue;  // pocket atoms only
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
			 if(this->atom[i].valid<2) continue;
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

		 // correct the contributions of the edging dots

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

float Protein::Calculate_Volume() const 
{
	int i,num;
	float volume=0.000;

	num=vol_dot.size(); 

	for(i=0;i<num;i++)
		{
		 if(vol_dot[i].valid==0) continue;
		 else volume+=vol_dot[i].unit;
		}

	return volume;
}

void Protein::Show_Volume_Dots(char *filename, char *show) const
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

void Protein::Generate_Surface_Dots(const Ligand *ligand, float probe_r)
// generate the solvent-accessible surface of the given molecule
// if probe_r=0.00, it represents vdw surface;
// if probe_r=WATER_R, it represents solvent accessible surface.
// note that it is a unit-atom model
{
	extern ForceField *ff;
	int i,j,k;
	bool mark;
	float d,dd,dmin;
	DotSet tmp_set;
	Dot tmp_dot;

	sur_dot.clear();	// clear all the existing surface dots, if any

	int *protein_check_list;

	protein_check_list=new int[this->num_atom];
	if(protein_check_list==NULL) Memory_Allocation_Error();	

	for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid<2) continue;
		 else if(!strcmp(this->atom[i].xtype,"H")) continue;

		 // first, check if this atom is close enough to the ligand

		 mark=false;

		 for(j=0;j<ligand->num_atom;j++)
			{
			 if(ligand->atom[j].valid==0) continue;
			 else if(!strcmp(ligand->atom[j].xtype,"H")) continue;

			 d=Distance(ligand->atom[j].coor,this->atom[i].coor);

			 if(d<(ligand->atom[j].R+this->atom[i].R+2*probe_r)) 
				{
				 mark=true; break;
				}
			 else continue;
			}

		 if(mark==false) continue;

		 // then calculate the surface dots for this atom

		 for(j=0;j<this->num_atom;j++)
			{
			 if(i==j||this->atom[j].valid<2)
				{
				 protein_check_list[j]=0; continue;
				}
			 else if((this->atom[j].type[0]=='H')&&
				 (!strcmp(this->atom[j].hb,"N"))) 
				{
				 protein_check_list[j]=0; continue;
				}

			 d=Distance(this->atom[i].coor,this->atom[j].coor);

			 if(d>(this->atom[i].R+this->atom[j].R+2*probe_r))
				{
				 protein_check_list[j]=0;
				}
			 else
				{
				 protein_check_list[j]=1;
				}
			}

		 tmp_set=ff->Get_Surface_Dot(this->atom[i],probe_r);

		 for(j=0;j<tmp_set.num_dot;j++) 
			{
			 // check whether this dot is on the surface or not 

			 mark=true; dmin=9999.0;

			 for(k=0;k<this->num_atom;k++)
			{
			 if(protein_check_list[k]==0) continue;

			 d=Distance(tmp_set.dot[j].coor,this->atom[k].coor);
			 dd=d-(this->atom[k].R+probe_r);

			 if(dd<0.000) {mark=false; break;} 
			 else if(dd<dmin) {dmin=dd; continue;}
			 else continue; 
			}

			 dd=sqrt(tmp_set.unit);

			 if(mark==false) tmp_set.dot[j].valid=0;
			 else if(dmin>=dd) tmp_set.dot[j].valid=1; // regular 
			 else tmp_set.dot[j].valid=2;	       // dots at edge
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
				 tmp_dot.unit=tmp_set.unit*0.667;
				}
			 
			 sur_dot.push_back(tmp_dot);
			}
		}

	if(protein_check_list) delete [] protein_check_list;

	return;
}

float Protein::Calculate_Surface() const
{
	int i,num;
	float surface=0.000;

	num=sur_dot.size();

	for(i=0;i<num;i++)
		{
		 if(sur_dot[i].valid==0) continue;
		 else surface+=sur_dot[i].unit;
		}

	return surface;
}

void Protein::Show_Surface_Dots(char *filename, char *show) const
// output surface dots to a PDB file
{
        FILE *fp;
	int i,num;
	float *property=NULL;

        if((fp=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	fprintf(fp,"HEADER    Molecular Surface Dots\n");
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

void Protein::Merge_Cofactor(const Ligand &cofactor)
{
	int i,j,last_id;
	Atom tmp_atom;

	last_id=this->num_atom;

	for(i=0;i<cofactor.num_atom;i++)
	{
	 tmp_atom=cofactor.atom[i];

	 // re-assign atom IDs

	 tmp_atom.id+=last_id;

	 for(j=0;j<tmp_atom.num_neib;j++) 
		{
		 if(tmp_atom.neib[j]<=0) tmp_atom.neib[j]=0;
		 else tmp_atom.neib[j]+=last_id;
		}

	 strcpy(tmp_atom.residue,"COF");  // special name for cofactor
	 strcpy(tmp_atom.res_id,"0 ");
	 tmp_atom.chain=' ';
	 tmp_atom.part=2;  	          // hetero-atoms
	 // note that tmp_atom.origin is still 1

	 this->atom.push_back(tmp_atom);
	}

	this->num_atom=this->atom.size();

	return;
}

// ***************************************************************************
// given a water, check its buried SAS inside protein's SAS and ligand's SAS 
// latest update: 03/07/2003
// ***************************************************************************
int Protein::Check_Buried_Ratio(float coor[]) const
{
	extern ForceField *ff;
	extern Ligand *ligand;
	int j,k,count;
	float d,buried_ratio;
	bool mark;
	DotSet tmp_set;

	// notice that only protein atoms are considered in this check
	// current criterion is 90% 

	buried_ratio=0.90;

	tmp_set=ff->Get_Surface_Dot(2*WATER_R,coor[0],coor[1],coor[2]);

	count=0;

	for(j=0;j<tmp_set.num_dot;j++)
		{
                 mark=false;	// check the protein

                 for(k=0;k<this->num_atom;k++)
                        {
                         if(this->atom[k].valid<2) continue;
			 else if(!strcmp(this->atom[k].type,"H")) continue;
                         else if(!strcmp(this->atom[k].type,"O.w")) continue;

                         d=Distance(tmp_set.dot[j].coor,this->atom[k].coor);

                         if(d>(WATER_R+this->atom[k].R)) continue;
                         else {mark=true; break;}
                        }

		 if(mark==true) {count++; tmp_set.dot[j].valid=0; continue;}

		 mark=false;	// check the ligand

		 for(k=0;k<ligand->num_atom;k++)
                        {
                         if(ligand->atom[k].valid<=0) continue;
                         else if(!strcmp(ligand->atom[k].type,"H")) continue;

                         d=Distance(tmp_set.dot[j].coor,ligand->atom[k].coor);

                         if(d>(WATER_R+ligand->atom[k].R)) continue;
                         else {mark=true; break;}
                        }

                 if(mark==true) {count++; tmp_set.dot[j].valid=0;}
		 else continue;
		}

	 // now determine this water is buried or not 

	 if(((float)count/tmp_set.num_dot)<buried_ratio) return FALSE; 
	 else return TRUE;
}

// ****************************************************************************
// this function processs an original PDB file (which typically comes from 
// manual process from Sybyl) and write out a standard PDB file meeting XTOOL
// criteria (all-atom, metal ions, waters, PO4 and SO4s). 
// this function is used by the PDBbind project
// latest update: 10/28/2003
// ****************************************************************************
void Protein::Fix_PDB_File(char *input, char *output)
{
	// read the original PDB file

        this->Read_From_PDB(input);   // do not read the CONECT section
        this->Value_Atom(0);	// do not calculate HB root

	// check structure and format 

	this->Refine_Structure();
	this->Rearrange_IDs();

	// output the new one 

	this->Write_Out_PDB(output);

	return;
}

// ****************************************************************************
// this function reads in an all-atom PDB file and write out an united-atom 
// PDB file meeting the XTOOL standard 
// this function is used by the PDBbind project
// latest update: 10/28/2003
// ****************************************************************************
void Protein::Eliminate_H_in_PDB(char *input, char *output)
{
	int i;

	this->Read_From_PDB(input);

	this->Value_Atom(0);

	for(i=0;i<this->num_atom;i++)
	{
	 if(strcmp(this->atom[i].xtype,"H")) continue;
	 else this->atom[i].valid=0;
	}

	this->Rearrange_IDs();

	this->Write_Out_PDB(output);

	return;
}

// ****************************************************************************
// this function checks the PDB, eliminate unwanted hydrogens, waters, metal
// ions, PO4, and SO4. 
// latest update: 07/30/2003
// ****************************************************************************
void Protein::Refine_Structure()
{
	extern Ligand *ligand;
	int i,j,count;
	bool mark;
	float d,cutoff;

	// first, eliminate all the hydrogen atoms on the hetero molecules
	// but the hydrogen atoms on the cofactor should remain

	count=0;

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid==0) continue;
	 else if(this->atom[i].part==1)
		{
		 if(!strcmp(this->atom[i].xtype,"H.hb")) count++;
		 else continue;
		}
	 else if(!strcmp(this->atom[i].residue,"COF")) continue;
	 else if(!strcmp(this->atom[i].xtype,"H")) this->atom[i].valid=0;
	 else if(!strcmp(this->atom[i].xtype,"H.hb")) this->atom[i].valid=0;
	 else continue;
	}

	if(count<5)
	{
	 printf("Warning: hydrogen atoms may not have been added on %s\n", 
		this->name);
	}

	// eliminate "suspended" waters, metal ions, SO4 and PO4

	cutoff=12.0; count=0;

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid==0) continue;
	 else if(this->atom[i].part==1) continue;

	 // check if this atom is close to any protein atom

	 mark=false; 

	 for(j=0;j<this->num_atom;j++)
		{
		 if(this->atom[j].valid==0) continue;
		 else if(this->atom[j].part>1) continue;
		 else if(!strcmp(this->atom[j].type,"H")) continue;

		 d=Distance(this->atom[i].coor, this->atom[j].coor);

		 if(d>cutoff) continue;
		 else {mark=true; break;}
		}

	 // check if this atom is close to any ligand atom, if ligand is given

	 if((mark==false)&&(ligand!=NULL)) 
	{
	 for(j=0;j<ligand->num_atom;j++)
		{
		 if(ligand->atom[j].valid==0) continue;
		 else if(!strcmp(ligand->atom[j].type,"H")) continue;

		 d=Distance(this->atom[i].coor, ligand->atom[j].coor);

		 if(d>cutoff) continue;
		 else {mark=true; break;}
		}
	}

	 if(mark==true) continue;
	 else {this->atom[i].valid=0; count++;}
	}

	// if any part of a PO4 or SO4 is within the cutoff, then keep it

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid==0) continue;
	 else if(this->atom[i].part==1) continue;
	 else if(strcmp(this->atom[i].residue,"PO4")&&
		 strcmp(this->atom[i].residue,"SO4")) continue;

	 mark=false;

	 for(j=0;j<this->num_atom;j++)
		{
		 if(this->atom[j].valid==0) continue;
		 else if(this->atom[j].part==1) continue;
		 else if(strcmp(this->atom[j].res_id,
                                this->atom[i].res_id)) continue;
		 else if(strcmp(this->atom[j].residue,
				this->atom[i].residue)) continue;

		 if(this->atom[j].valid==0) continue;
		 else {mark=true; break;}
		}

	 if(mark==true)
	{
	 for(j=0;j<this->num_atom;j++)
		{
		 if(this->atom[j].valid==0) continue;
		 else if(this->atom[j].part==1) continue;
		 else if(strcmp(this->atom[j].res_id,
				this->atom[i].res_id)) continue;
		 else if(strcmp(this->atom[j].residue,
				this->atom[i].residue)) continue;
		
		 if(this->atom[j].valid==0) 
			{
			 this->atom[j].valid=1;
			 count--;
			}
		}
	}
	}

	if(count>0) printf("%d hetero-atoms are eliminated\n", count);

	return;
}

// **************************************************************************
// this routine reassign continuous IDs for regular atoms and hetero-atoms
// to meet the requirments of X-Tool
// latest update: 07/29/03
// **************************************************************************
void Protein::Rearrange_IDs()
{
	int i,j,k,id1,id2,count;
	char res_id[10];
	Atom temp;
	vector <int> old_id;
	vector <int> new_id;

	// now re-assign residue IDs for hetero-atoms, started from 1
	// give metal ions, PO4, SO4 the priority 
	// note that res_id 0 is reserved for cofactor molecule

	count=0; strcpy(res_id,"0 ");

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid==0) continue;
	 else if(this->atom[i].part==1) continue;
	 else if(!strcmp(this->atom[i].xtype,"M+")||
		 !strcmp(this->atom[i].xtype,"F-")||
		 !strcmp(this->atom[i].xtype,"Cl-")||
		 !strcmp(this->atom[i].xtype,"Br-")||
		 !strcmp(this->atom[i].xtype,"I-"))
		{
		 strcpy(res_id,this->atom[i].res_id); count++;
		 sprintf(this->atom[i].res_id,"%d ", count);
		}
	 else if(!strcmp(this->atom[i].residue,"SO4")||
		 !strcmp(this->atom[i].residue,"PO4"))
		{
		 if(strcmp(this->atom[i].res_id,res_id))
			{
			 strcpy(res_id,this->atom[i].res_id); count++;
			 sprintf(this->atom[i].res_id,"%d ", count);
			}
		 else
			{
			 sprintf(this->atom[i].res_id,"%d ", count);
			}
		}
	 else continue;
	}

	// now assign residue ID to water molecules

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid==0) continue;
	 else if(this->atom[i].part==1) continue;
	 else if(strcmp(this->atom[i].xtype,"O.w")) continue;
	 else {count++; sprintf(this->atom[i].res_id,"%d ", count);}
	}

	// now arrange hetero-atoms according to their residue IDs

	for(i=0;i<this->num_atom-1;i++)
	{
	 if(this->atom[i].part==1) continue;

	 for(j=i+1;j<this->num_atom;j++)
		{
		 if(this->atom[j].part==1) continue;

		 sscanf(this->atom[i].res_id,"%d",&id1);
		 sscanf(this->atom[j].res_id,"%d",&id2);

		 if(id1<id2) continue;
		 else if(id1==id2)
			{
			 if(this->atom[i].id<=this->atom[j].id) continue;
			 else {SWAP(this->atom[i],this->atom[j]);}
			}
		 else {SWAP(this->atom[i],this->atom[j]);}
		}
	}

	// re-assign continous IDs for all atoms, started from 1

	count=0; old_id.clear(); new_id.clear();

	// first, the regular atoms

	for(i=0;i<this->num_chain;i++)
       		{
       	 	 for(j=0;j<this->num_atom;j++)
                {
                 if(this->atom[j].valid<=0) continue;
                 else if(this->atom[j].part!=1) continue;
                 else if(this->atom[j].chain!=this->chain[i].label) continue;
                 else 
			{
			 count++;
			 id1=j+1; old_id.push_back(id1);  
			 id2=count; new_id.push_back(id2); 
			 this->atom[j].id=id2;
			} 
		}
		}

	// second, the hetero atoms

	for(i=0;i<this->num_atom;i++)
                {
                 if(this->atom[i].valid<=0) continue;
                 else if(this->atom[i].part<=1) continue;
		 else
			{
			 count++;
			 id1=i+1; old_id.push_back(id1);  
			 id2=count; new_id.push_back(id2); 
			 this->atom[i].id=id2;
			} 
		}

	// also, update the neib IDs for each atom

	count=old_id.size();

	for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid<=0) continue;

		 for(j=0;j<this->atom[i].num_neib;j++)
			{
			 if(this->atom[i].neib[j]<=0) continue;
			 for(k=0;k<count;k++)
				{
				 if(this->atom[i].neib[j]!=old_id[k]) continue;
				 else {this->atom[i].neib[j]=new_id[k];break;}
				}
			}
		}

	return;
}

void Protein::Classify_Residues(char *filename)
{
	extern Ligand *ligand; 
	FILE *fp;
	int i,j,k,l,num;
	float cutoff,ratio,area_total,area_exposed,d;
	bool mark;
	Residue tmp_residue;

	this->Analyze_Sequence();

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].xtype,"H")) this->atom[i].valid=0;
	 else if(!strcmp(this->atom[i].xtype,"O.w")) this->atom[i].valid=0;
	 else continue;
	}

	if((fp=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	for(i=0;i<this->num_chain;i++)
	{
	 for(j=0;j<chain[i].length;j++)
	{
	 mark=tmp_residue.Get_Contents_From_Protein(chain[i].residue[j].name, 
						    chain[i].residue[j].id, 
						    chain[i].residue[j].chain);

	 if(!mark) continue; // no necessary information for this residue 

	 // calculate surface area first

	 tmp_residue.Generate_Surface_Dots(WATER_R); 

	 num=tmp_residue.sur_dot.size(); 
	 area_exposed=area_total=0.000;

	 for(k=0;k<num;k++)
		{
		 area_total+=tmp_residue.sur_dot[k].unit;

		 if(tmp_residue.sur_dot[k].valid<=0) continue;
		 else area_exposed+=tmp_residue.sur_dot[k].unit;
		}

	 ratio=area_exposed/area_total;

	 // determine if it is in close contact with the ligand 

	 mark=false;

	 for(k=0;k<tmp_residue.num_atom;k++)
		{
		 if(tmp_residue.atom[k].valid<=0) continue;

		 for(l=0;l<ligand->num_atom;l++)
			{
			 if(ligand->atom[l].valid<=0) continue;
			 else if(!strcmp(ligand->atom[l].type,"H")) continue;

			 d=Distance(tmp_residue.atom[k].coor,
				    ligand->atom[l].coor);

			 if(d>4.00) continue;
			 else {mark=true; break;}
			}

		 if(mark==true) break;
		}

	 // now do the classification and output the list

	 cutoff=0.33;

	 fprintf(fp,"%3s  ", tmp_residue.name); 

	 if(tmp_residue.chain==' ') fprintf(fp,"%c  ",'_');
	 else fprintf(fp,"%c  ",tmp_residue.chain);

	 fprintf(fp,"%4s  ", tmp_residue.id);

	 if(ratio>cutoff) fprintf(fp,"SURFACE  ");
	 else 
		{
		 fprintf(fp,"CORE     ");
/*
		 for(k=0;k<this->num_atom;k++)
		{
		 if(strcmp(this->atom[k].residue,tmp_residue.name)) continue;
		 else if(strcmp(this->atom[k].res_id,tmp_residue.id)) continue;
		 else if(this->atom[k].chain!=tmp_residue.chain) continue;
		 this->atom[k].valid=0;
		}
*/
		}

	 if(mark==true) fprintf(fp,"POCKET  ");
	 else fprintf(fp,"______  ");

	 fprintf(fp,"%6.2f  %6.2f  %4.2f\n", area_total, area_exposed, ratio);
	}
	}

	// Write_Out_PDB("temp_surface.pdb");

	fclose(fp); return;
}

