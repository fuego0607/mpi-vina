# include "xtool.h"

Ligand::Ligand()
{
	Clear();
}

Ligand::~Ligand()
{
	if(abs_inf) delete [] abs_inf;
}

Ligand::Ligand(int max_atom_num, int max_bond_num)
{
	Clear();

	atom=new Atom[max_atom_num];
	if(atom==NULL) Memory_Allocation_Error();

	bond=new Bond[max_bond_num];
	if(bond==NULL) Memory_Allocation_Error();
}

Ligand::Ligand(const Ligand &original)
{
	// this part is from Molecule

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

	// this part is for Ligand

        bind_score=original.bind_score;
        chem_score=original.chem_score;
	vdw=original.vdw;
	sb=original.sb;
	hb=original.hb;
	hp=original.hp;
	hm=original.hm;
	hs=original.hs;
	ar=original.ar;
	rt=original.rt;
	pmf=original.pmf;
	uhb=original.uhb;
	bnsur=original.bnsur;
	bpsur=original.bpsur;
	pkd1=original.pkd1; 
	pkd2=original.pkd2;
	pkd3=original.pkd3;

	if(original.abs_inf!=NULL)
	{
	 abs_inf=new ABS[num_atom]; if(abs_inf==NULL) Memory_Allocation_Error();
	 for(i=0;i<num_atom;i++) abs_inf[i]=original.abs_inf[i];
	}
}

Ligand& Ligand::operator = (const Ligand &original)
{
	if(this==&original) return *this;

	// this part is from Molecule

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

	// this part is for Ligand

        bind_score=original.bind_score;
        chem_score=original.chem_score;
	vdw=original.vdw;
	sb=original.sb;
	hb=original.hb;
	hp=original.hp;
	hm=original.hm;
	hs=original.hs;
	ar=original.ar;
	rt=original.rt;
	pmf=original.pmf;
	uhb=original.uhb;
	bnsur=original.bnsur;
	bpsur=original.bpsur;
	pkd1=original.pkd1;
	pkd2=original.pkd2;
	pkd3=original.pkd3;

	if(abs_inf) delete [] abs_inf; abs_inf=NULL;

	if(original.abs_inf!=NULL)
	{
	 abs_inf=new ABS[num_atom]; if(abs_inf==NULL) Memory_Allocation_Error();
	 for(i=0;i<num_atom;i++) abs_inf[i]=original.abs_inf[i];
	}

	return *this;
}

void Ligand::Clear()
{
	Molecule::Clear();

        bind_score=chem_score=0.000;
	vdw=0.000;
	sb=0.000;
	hb=0.000;
	hp=0.000; 
	hm=0.000;
	hs=0.000;
	ar=0.000;
	rt=0.000;
	pmf=0.000;
	uhb=0.000;
	bnsur=bpsur=0.000;
	pkd1=pkd2=pkd3=0.000;

	abs_inf=NULL;

	return;
}

void Ligand::Show_Contents() const
{
	Molecule::Show_Contents();

	return;
}

int Ligand::Value_Atom()
{
	int mark;

	mark=Molecule::Value_Atom();

	// now compute the H-bond root for each HB atom

	Calculate_HB_Root();

	// assign atomic logp values, which is needed by binding score
	// atomic logp is not read from ATOM_DEF_XTOOL!

	this->logp=Calculate_LogP(); 

	// some other molecular properties

	this->num_hb_atom=this->Get_Num_HB_Atom();
	this->num_rotor=this->Count_Rotor();

	return mark;
}

void Ligand::Calculate_HB_Root()
{
	int i,j,tmp,num_nonh;
	float tmpx,tmpy,tmpz;

	for(i=0;i<num_atom;i++)
		{
		 if(atom[i].valid==0) continue;
		 else if(!strcmp(atom[i].hb,"N")) continue;
		 else if(!strcmp(atom[i].hb,"H")) continue;
		 else if(!strcmp(atom[i].hb,"P")) continue;

		 tmpx=tmpy=tmpz=0.000; num_nonh=0;
		
		 for(j=0;j<atom[i].num_neib;j++)
			{
			 tmp=atom[i].neib[j]-1;
			 if(atom[tmp].type[0]=='H') continue;
			 else
				{
				 tmpx+=atom[tmp].coor[0];
				 tmpy+=atom[tmp].coor[1];
				 tmpz+=atom[tmp].coor[2];
				 num_nonh++;
				}
			}

		 if(num_nonh==0) strcpy(atom[i].hb,"P");
		 else
			{
		 	 tmpx/=num_nonh;
			 tmpy/=num_nonh;
			 tmpz/=num_nonh;
			 atom[i].root[0]=tmpx;
			 atom[i].root[1]=tmpy;
			 atom[i].root[2]=tmpz;
			}
		}

	return;
}

int Ligand::Chemical_Viability_Check() const
{
 extern Input *input;

 if(!strcmp(((SCORE_Input *)input)->apply_chemical_rules,"YES"))
{
 if(this->weight>((SCORE_Input *)input)->max_weight) return FALSE;
 else if(this->weight<((SCORE_Input *)input)->min_weight) return FALSE; 
 else if(this->logp>((SCORE_Input *)input)->max_logp) return FALSE; 
 else if(this->logp<((SCORE_Input *)input)->min_logp) return FALSE; 
 else if(this->num_hb_atom>((SCORE_Input *)input)->max_hb_atom) return FALSE; 
 else if(this->num_hb_atom<((SCORE_Input *)input)->min_hb_atom) return FALSE; 
}

 return TRUE;
}

float Ligand::Calculate_Buried_Surface() const
{
	int i;
	float sum,total,buried;

	sum=0.000;

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].xtype,"H")) continue; // non-polar H 

	 this->Atom_Buried_Surface(i+1,total,buried);

	 sum+=buried;
	}

	return sum;
}

float Ligand::Atom_Buried_Surface(int id, float &total, float &buried) const
{
 extern Protein *protein;
 int j,k,num;
 bool mark;
 float d,ratio;
 int *atom_check_list;

 atom_check_list=new int[protein->num_atom];
 if(atom_check_list==NULL) Memory_Allocation_Error();

 total=buried=0.000;

 for(j=0;j<protein->num_atom;j++)
       	{
         if((protein->atom[j].valid!=2)|| 
            (!strcmp(protein->atom[j].xtype,"H"))||
            (!strcmp(protein->atom[j].xtype,"O.w")))
		{
		 atom_check_list[j]=0; continue;
		}

         d=Distance(this->atom[id-1].coor,protein->atom[j].coor);

         if(d>(this->atom[id-1].R+protein->atom[j].R+2*WATER_R))
                {
                 atom_check_list[j]=0;
                }
         else
                {
                 atom_check_list[j]=1;
                }
        }

 num=this->sur_dot.size();

 for(j=0;j<num;j++) 
        {
	 if(this->sur_dot[j].valid!=id) continue;

	 // check if this dot is buried 

	 mark=false;

	 for(k=0;k<protein->num_atom;k++)
                 {
                  if(atom_check_list[k]==0) continue;

                  d=Distance(this->sur_dot[j].coor, protein->atom[k].coor);

                  if(d>(protein->atom[k].R+WATER_R)) continue;
                  else {mark=true; break;}
                 }

	 total+=this->sur_dot[j].unit;

 	 if(mark==false) continue;
	 else buried+=this->sur_dot[j].unit;
       	}

 if(atom_check_list) delete [] atom_check_list;

 if(total<=0.00) ratio=0.00;
 else ratio=buried/total; 

 return ratio;
}

// ***************************************************************************
// Check buried status of each ligand atom: if the buried ratio exceeds
// a cutoff value, then assign valid=2 for this atom.
// This functions also returns the total buried surface area.
// Latest update: 11/18/2003
// ***************************************************************************
float Ligand::Check_Buried_Status(float cutoff)
{
	int i;
	float sum,ratio,total,buried;

	sum=0.000;

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].xtype,"H")) continue; // non-polar H 

	 ratio=Atom_Buried_Surface(i+1,total,buried);
	 sum+=buried;

	 if(ratio>cutoff) this->atom[i].valid=2;  // buried atom
	 else continue;
	}

	return sum;
}

