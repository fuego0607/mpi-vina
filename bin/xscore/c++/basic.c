# include "xtool.h"

Atom::Atom() {Clear();}

Atom::~Atom() {}

void Atom::Clear()
{
	id=0; valid=0; part=1; origin=1; 

	strcpy(name,"Un"); strcpy(type,"Un"); 
	strcpy(xtype,"Un"); strcpy(type2,"Un"); 
	strcpy(residue,"Un"); strcpy(res_id,"0"); chain=' ';

	coor[0]=coor[1]=coor[2]=0.000;
	root[0]=root[1]=root[2]=0.000;

	r=0.000; eps=0.000; q=0.000; weight=0.000; R=0.000; 
	strcpy(hb,"N"); logp=0.000; solv=0.000; score=0.000;
	ring=0;

	occupancy=1.000; bfactor=0.000;

	num_neib=num_nonh=0; 

	int i;

	for(i=0;i<MAX_ATOM_NEIB;i++) neib[i]=0; 
	for(i=0;i<MAX_ATOM_NEIB;i++) bond[i]=0;

	temp=0;
}

void Atom::Show_Contents() const
{
	printf("Atom: ");
	printf("%-4d ",id);
	printf("%-5s ",name);
	printf("%8.3f ",coor[0]);
	printf("%8.3f ",coor[1]);
	printf("%8.3f ",coor[2]);
	printf("%-5s ",type);
	printf("\n");

	return;
}

int Atom::Get_Donor_Type() const
{
 if(strcmp(this->hb,"D")&&strcmp(this->hb,"DA")&&strcmp(this->hb,"M")) return 0;

 int type; 

 if(this->origin==2) // protein atoms
	{
 	 if(!strcmp(this->hb,"M"))
		{
	 	 type=1;
		}
 	 else if(!strcmp(this->type,"O.3")||
		!strcmp(this->type,"O.co2")||
         	!strcmp(this->type,"O.w"))
       		{
         	 type=2;
        	}
 	 else if(strstr(this->xtype,"N.4")||
		strstr(this->xtype,"N.3"))
        	{
         	 type=2;
        	}
 	 else if(strstr(this->xtype,"N.pl3")||
		strstr(this->xtype,"N.2")||
         	strstr(this->xtype,"N.ar"))
        	{
	 	 if(strstr(this->name,"NH")&&strstr(this->residue,"ARG")) type=2;
	 	 else if(strstr(this->name,"ND")&&strstr(this->residue,"ASN")) type=2;
	 	 else if(strstr(this->name,"NE")&&strstr(this->residue,"GLN")) type=2;
         	else type=1;
        	}
 	 else if(!strcmp(this->type,"S.3"))
        	{
         	 type=2;
        	}
 	 else type=2;
	}
 else  // ligand atoms
	{
 	 if(!strcmp(this->type,"O.3")||!strcmp(this->type,"O.co2"))
     		{
         	 type=2;
        	}
 	 else if(strstr(this->xtype,"N.4")||strstr(this->xtype,"N.3"))
        	{
         	 if(this->num_nonh<=2) type=2;
         	 else type=1;
        	}
 	 else if(strstr(this->xtype,"N.pl3")||strstr(this->xtype,"N.2")||
         	 strstr(this->xtype,"N.ar"))
        	{
         	 if(this->num_nonh<=1) type=2;
         	 else type=1;
        	}
 	 else if(!strcmp(this->type,"S.3")) type=2;
 	 else type=2;
	}

 return type;
}

int Atom::Get_Acceptor_Type() const
{
 if(strcmp(this->hb,"A")&&strcmp(this->hb,"DA")) return 0;

 int type;

 if(this->origin==2)  // protein atoms
	{
 	 if(!strcmp(this->type,"O.3")||
	    !strcmp(this->type,"O.2")||
    	    !strcmp(this->type,"O.co2")||
	    !strcmp(this->type,"O.w"))
       		{
         	 type=2;
        	}
 	 else if(strstr(this->xtype,"N.pl3")||
		 strstr(this->xtype,"N.2")||
         	 strstr(this->xtype,"N.ar"))
        	{
         	 type=1;
        	}
 	 else type=2;
	}
 else
	{
 	 if(!strcmp(this->type,"O.3")||
	    !strcmp(this->type,"O.2")||
    	    !strcmp(this->type,"O.co2"))
       		{
        	 type=2;
       		}
 	 else if(strstr(this->xtype,"N.pl3")||
		 strstr(this->xtype,"N.2")||
         	 strstr(this->xtype,"N.ar"))
       		{
		 if(this->num_nonh<=1) type=2;
        	 else type=1;
       		}
 	 else type=2;
	}

 return type;
}

Water::Water() {Clear();}

Water::~Water() {}

void Water::Clear()
{
	id=0; valid=0; 

	strcpy(name,"O"); 
	strcpy(type,"O.w"); 
	strcpy(xtype,"O.w"); 
	strcpy(residue,"HOH");
	strcpy(hb,"DA");

	coor[0]=coor[1]=coor[2]=0.000;

	r=1.770; q=0.000; eps=0.116; 
	logp=-0.500; depth=0.000; score=0.000;
}

void Water::Show_Contents() const
{
	printf("Water: ");
	printf("%-4d ",id);
	printf("%-5s ",name);
	printf("%8.3f ",coor[0]);
	printf("%8.3f ",coor[1]);
	printf("%8.3f ",coor[2]);
	printf("%-5s ",type);
	printf("\n");

	return;
}

Bond::Bond() {Clear();}

Bond::~Bond() {}

void Bond::Clear()
{
	int i;

	id=0; valid=0; part=1; ring=0;

	atom_1=atom_2=0;

	strcpy(type,"1"); length=0.000;

	num_neib=0;
	for(i=0;i<MAX_BOND_NEIB;i++) neib[i]=0;
}

void Bond::Show_Contents() const
{
	printf("Bond %s between: %d - %d\n", type, atom_1, atom_2);

	return;
}

Torsion::Torsion() {Clear();}

Torsion::~Torsion() {}

void Torsion::Clear()
{
	strcpy(type,"");
	angle=0; e=0.000; V=0.000; n=0; S=0;
}

void Torsion::Show_Contents() const
{
	printf("Torsion %s between:\n", type);
        printf("\tAtom 1\t%d\t%s\t%s\n", atom_1.id, atom_1.name, atom_1.type);
        printf("\tAtom 2\t%d\t%s\t%s\n", atom_2.id, atom_2.name, atom_2.type);
	printf("\tAtom 3\t%d\t%s\t%s\n", atom_3.id, atom_3.name, atom_3.type);
        printf("\tAtom 4\t%d\t%s\t%s\n", atom_4.id, atom_4.name, atom_4.type);

        return;
}

Group::Group() {Clear();}

Group::~Group() {}

void Group::Clear()
{
	valid=0;
	num_neib=0; num_nonh=0;	num_h=0; 
	num_hetero=0; num_pi=0;	num_car=0; num_nar=0;
	num_db=0; db_type=0; num_tb=0; 
	amide=0;

	int i; center.Clear();

	for(i=0;i<MAX_ATOM_NEIB;i++) neib[i].Clear(); 
	for(i=0;i<MAX_ATOM_NEIB;i++) bond[i].Clear(); 
}

void Group::Show_Contents() const
{
	int i;

	printf("Group center:\n");

	center.Show_Contents();

	for(i=0;i<num_neib;i++)
		{
		 printf("Neighbor %d\t", i+1);
		 printf("Atom %d\t%s\t%s\tvia Bond %s\n", 
			neib[i].id, neib[i].name, neib[i].type, bond[i].type);
		}

	return;
}

HBond::HBond() {Clear();}

HBond::~HBond() {}

void HBond::Clear()
{
	valid=latom=patom=0;
	type=d_type=a_type=0; 
	score=0.000; sb=false;
	d=a0=a1=a2=0.000;
	D.Clear(); H.Clear(); A.Clear();
}

void HBond::Show_Contents() const
{
	printf("HBond type: %d -------------------------------------\n", type);

	printf("Hydrogen: %d %s %8.3f %8.3f %8.3f\n", 
		H.id, H.xtype, H.coor[0], H.coor[1], H.coor[2]);

	printf("Donor: %d %s %s %s type=%d ", 
		D.id, D.type, D.residue, D.res_id, d_type);

	printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
		D.coor[0], D.coor[1], D.coor[2],
		D.root[0], D.root[1], D.root[2]);

	printf("Acceptor: %d %s %s %s type=%d ", 
		A.id, A.type, A.residue, A.res_id, a_type);

	printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
		A.coor[0], A.coor[1], A.coor[2],
		A.root[0], A.root[1], A.root[2]);

	printf("D-A distance = %5.2f\n", d);
	printf("D-H-A angle  = %5.1f\n", a0);
	printf("DR-D-A angle = %5.1f\n", a1);
	printf("D-A-AR angle = %5.1f\n", a2);
	printf("Score = %5.3f\n", score);
}

// *************************************************************************
// Evaluate a given h-bond, which relies on explicit hydrogen atoms
// latest update: 11/18/2003
// *************************************************************************
int HBond::Value_HBond()
{
	float d,d0,d1,d2,d3,d4;
	float a0,a1,a2,angle1,angle2,angle3,angle4;
	float tmp1,tmp2,tmp3,tmp4;
	int mark0,mark1,mark2;

	// first, calculate the necessary parameters

	d=Distance(this->D.coor,this->A.coor);

	// angle D-H-A 

	if(!strcmp(this->D.hb,"M"))
		{
		 a0=0.1; mark0=false;
		}
	else if(!strcmp(this->D.type,"O.3"))
		{
	 	 // align H to the optimal position and then recalculate <D-H-A 

	 	 float a,b,c,A,C;

	 	 b=0.98; // D-H length 
	 	 c=d;	 // D-A length 
	 	 A=fabs(Angle(this->D.root,this->D.coor,this->A.coor)-109.0); 
	 	 a=sqrt(b*b+c*c-2*b*c*cos(A/180.0));	// H-A length
		 C=acos((a*a+b*b-c*c)/(2*a*b))/PI*180.0;
	 	 a0=180.0-C;
		 mark0=true;
		}
	else 
		{
		 a0=180.0-Angle(this->D.coor,this->H.coor,this->A.coor);
		 mark0=true;
		}

	// angle DR-D-A

	if(!strcmp(this->D.hb,"M")||!strcmp(this->D.type,"O.w")) 
		{
		 a1=0.1; mark1=false; 
		}
	else 
		{
		 a1=180.0-Angle(this->D.root,this->D.coor,this->A.coor);
		 mark1=true;
		}

	// D-A-AR

	if(!strcmp(this->A.hb,"M")||!strcmp(this->A.type,"O.w")) 
		{
		 a2=0.1; mark2=false; 
		}
	else 
		{
		 a2=180.0-Angle(this->D.coor,this->A.coor,this->A.root);
		 mark2=true;
		}

	// now determine the geometry type of this h-bond

	this->Determine_DA_Type();

	// now compute the strength of this h-bond

	d0=this->D.R+this->A.R; d1=0.00; d2=1.00; d3=d0-0.40; d4=d0+0.10;

	if(d<d1) tmp1=0.000;
	else if(d<d2) tmp1=(d-d1)/(d2-d1);
	else if(d<d3) tmp1=1.000;
	else if(d<d4) tmp1=(d4-d)/(d4-d3);
	else tmp1=0.000;

	if(mark0==true)
		{
		 angle1=0.0; angle2=0.001; angle3=60.0; angle4=90.0;
		 if(a0<angle1) tmp2=0.000;
        	 else if(a0<angle2) tmp2=(a0-angle1)/(angle2-angle1);
        	 else if(a0<angle3) tmp2=1.000;
        	 else if(a0<angle4) tmp2=(angle4-a0)/(angle4-angle3);
        	 else tmp2=0.000;
		}
	else tmp2=1.000;

	if(mark1==true)
		{
		 if(d_type==1) 
			{angle1=0.0; angle2=0.001; angle3=30.0; angle4=60.0;}
		 else 
			{angle1=20.0; angle2=40.0; angle3=70.0; angle4=90.0;}

        	 if(a1<angle1) tmp3=0.000;
        	 else if(a1<angle2) tmp3=(a1-angle1)/(angle2-angle1);
        	 else if(a1<angle3) tmp3=1.000;
        	 else if(a1<angle4) tmp3=(angle4-a1)/(angle4-angle3);
        	 else tmp3=0.000;
		}
	else tmp3=1.000;

	if(mark2==true)
		{
		 if(a_type==1) 
                	{angle1=0.0; angle2=0.001; angle3=30.0; angle4=60.0;}
        	 else {angle1=0.0; angle2=5.0; angle3=60.0; angle4=90.0;}

        	 if(a2<angle1) tmp4=0.000;
        	 else if(a2<angle2) tmp4=(a2-angle1)/(angle2-angle1);
        	 else if(a2<angle3) tmp4=1.000;
        	 else if(a2<angle4) tmp4=(angle4-a2)/(angle4-angle3);
        	 else tmp4=0.000;
		}
	else tmp4=1.000;

	this->d=d; this->a0=a0; this->a1=a1; this->a2=a2;

	this->score=tmp1*tmp2*tmp3*tmp4;
	// this->score=tmp1*tmp3*tmp4;

	if(this->score<0.001) return FALSE;
	else return TRUE;
}

void HBond::Determine_DA_Type()
{
 this->d_type=this->D.Get_Donor_Type();
 this->a_type=this->A.Get_Acceptor_Type();

 return;
}

// *************************************************************************
// Evaluate a given h-bond, which does not rely on explicit hydrogen atoms
// latest update: 11/21/2003
// *************************************************************************
int HBond::Value_HBond_2()
{
	float d,d0,d1,d2,d3,d4;
	float a1,a2,angle1,angle2,angle3,angle4;
	float tmp1,tmp3,tmp4;
	int mark1,mark2;

	// first, calculate the necessary parameters

	d=Distance(this->D.coor,this->A.coor);

	// angle DR-D-A

	if(!strcmp(this->D.hb,"M")||!strcmp(this->D.type,"O.w")) 
		{
		 mark1=false; a1=0.0; 
		}
	else 
		{
		 a1=180.0-Angle(this->D.root,this->D.coor,this->A.coor);
		 mark1=true;
		}

	// D-A-AR

	if(!strcmp(this->A.hb,"M")||!strcmp(this->A.type,"O.w")) 
		{
		 mark2=false; a2=0.0; 
		}
	else 
		{
		 a2=180.0-Angle(this->D.coor,this->A.coor,this->A.root);
		 mark2=true;
		}

	// now determine the geometry type of this h-bond

	this->Determine_DA_Type();

	// now compute the strength of this h-bond

	d0=this->D.R+this->A.R; 
	d1=0.00; d2=1.00; d3=d0-0.40; d4=d0+0.20;

	if(d<d1) tmp1=0.000;
	else if(d<d2) tmp1=(d-d1)/(d2-d1);
	else if(d<d3) tmp1=1.000;
	else if(d<d4) tmp1=(d4-d)/(d4-d3);
	else tmp1=0.000;

	if(mark1==true)
		{
		 if(d_type==1) 
			{angle1=0.0; angle2=0.001; angle3=25.0; angle4=50.0;}
		 else 
			{angle1=25.0; angle2=50.0; angle3=75.0; angle4=100.0;}

        	 if(a1<angle1) tmp3=0.000;
        	 else if(a1<angle2) tmp3=(a1-angle1)/(angle2-angle1);
        	 else if(a1<angle3) tmp3=1.000;
        	 else if(a1<angle4) tmp3=(angle4-a1)/(angle4-angle3);
        	 else tmp3=0.000;
		}
	else tmp3=1.000;

	if(mark2==true)
		{
		 if(a_type==1) 
                	{angle1=0.0; angle2=0.001; angle3=30.0; angle4=55.0;}
        	 else {angle1=0.0; angle2=5.0; angle3=70.0; angle4=95.0;}

        	 if(a2<angle1) tmp4=0.000;
        	 else if(a2<angle2) tmp4=(a2-angle1)/(angle2-angle1);
        	 else if(a2<angle3) tmp4=1.000;
        	 else if(a2<angle4) tmp4=(angle4-a2)/(angle4-angle3);
        	 else tmp4=0.000;
		}
	else tmp4=1.000;

	this->d=d; this->a1=a1; this->a2=a2;

	if(tmp3>=tmp4) this->score=tmp1*tmp4;
	else this->score=tmp1*tmp3;

	// this->score=tmp1*tmp3*tmp4;

	if(this->score<0.001) return FALSE;
	else return TRUE;
}

int HBond::Value_SBond()
{
	float d,d0,d1,d2,d3,d4;
	float a1,a2,angle1,angle2,angle3,angle4;
	float tmp1,tmp3,tmp4;
	int mark1,mark2;

	// first, calculate the necessary parameters

	d=Distance(this->D.coor,this->A.coor);

	// angle DR-D-A

	if(!strcmp(this->D.hb,"M")||!strcmp(this->D.type,"O.w")) 
		{
		 mark1=false; 
		}
	else 
		{
		 a1=180.0-Angle(this->D.root,this->D.coor,this->A.coor);
		 mark1=true;
		}

	// D-A-AR

	if(!strcmp(this->A.hb,"M")||!strcmp(this->A.type,"O.w")) 
		{
		 mark2=false; 
		}
	else 
		{
		 a2=180.0-Angle(this->D.coor,this->A.coor,this->A.root);
		 mark2=true;
		}

	// now determine the geometry type of this h-bond

	this->Determine_DA_Type();

	// now compute the strength of this h-bond

	d0=this->D.R+this->A.R; 
	d1=0.00; d2=1.00; d3=d0-0.50; d4=d0+0.50;

	if(d<d1) tmp1=0.000;
	else if(d<d2) tmp1=(d-d1)/(d2-d1);
	else if(d<d3) tmp1=1.000;
	else if(d<d4) tmp1=(d4-d)/(d4-d3);
	else tmp1=0.000;

	if(mark1==true)
		{
		 if(d_type==1) 
			{angle1=0.0; angle2=0.001; angle3=30.0; angle4=60.0;}
		 else 
			{angle1=25.0; angle2=45.0; angle3=80.0; angle4=90.0;}

        	 if(a1<angle1) tmp3=0.000;
        	 else if(a1<angle2) tmp3=(a1-angle1)/(angle2-angle1);
        	 else if(a1<angle3) tmp3=1.000;
        	 else if(a1<angle4) tmp3=(angle4-a1)/(angle4-angle3);
        	 else tmp3=0.000;
		}
	else tmp3=1.000;

	if(mark2==true)
		{
		 if(a_type==1) 
                	{angle1=0.0; angle2=0.001; angle3=30.0; angle4=55.0;}
        	 else {angle1=0.0; angle2=5.0; angle3=70.0; angle4=80.0;}

        	 if(a2<angle1) tmp4=0.000;
        	 else if(a2<angle2) tmp4=(a2-angle1)/(angle2-angle1);
        	 else if(a2<angle3) tmp4=1.000;
        	 else if(a2<angle4) tmp4=(angle4-a2)/(angle4-angle3);
        	 else tmp4=0.000;
		}
	else tmp4=1.000;

	this->d=d; this->a1=a1; this->a2=a2;

	if(tmp3<=tmp4) this->score=tmp1*tmp3;
	else this->score=tmp1*tmp4;

	if((this->D.q*this->A.q)>0.000) this->score*=(-1.00);  // bad SB	
	else this->score*=(1.00);  // good SB

	if(fabs(this->score)<0.001) return FALSE;
	else return TRUE;
}

Ring::Ring() {Clear();}

Ring::~Ring() {atom_id.clear(); bond_id.clear();}

Ring::Ring(const Ring &original)
{
 this->valid=original.valid;
 this->type=original.type;
 this->num_member=original.num_member;
 this->centroid[0]=original.centroid[0];
 this->centroid[1]=original.centroid[1];
 this->centroid[2]=original.centroid[2];

 int i,n;

 this->atom_id.clear(); n=original.atom_id.size();
 for(i=0;i<n;i++) this->atom_id.push_back(original.atom_id[i]);

 this->bond_id.clear(); n=original.bond_id.size();
 for(i=0;i<n;i++) this->bond_id.push_back(original.bond_id[i]);
}

Ring& Ring::operator = (const Ring &original)
{
 if(this==&original) return *this;

 this->valid=original.valid;
 this->type=original.type;
 this->num_member=original.num_member;
 this->centroid[0]=original.centroid[0];
 this->centroid[1]=original.centroid[1];
 this->centroid[2]=original.centroid[2];

 int i,n;

 this->atom_id.clear(); n=original.atom_id.size();
 for(i=0;i<n;i++) this->atom_id.push_back(original.atom_id[i]);

 this->bond_id.clear(); n=original.bond_id.size();
 for(i=0;i<n;i++) this->bond_id.push_back(original.bond_id[i]);

 return *this;
}

void Ring::Clear()
{
	valid=0; type=0;
	num_member=0; atom_id.clear(); bond_id.clear();

	centroid[0]=centroid[1]=centroid[2]=0.000;
}

void Ring::Show_Contents() const
{
	int i;

	printf("Ring: type = %d ------------------------\n",type);

	for(i=0;i<num_member;i++)
		{
		 printf("Member %d: %d\n", i+1, atom_id[i]); 
		}

	printf("Centroid: %8.3f %8.3f %8.3f\n", 
		centroid[0], centroid[1], centroid[2]);

	printf("----------------------------------------\n");
}

Dot::Dot() {Clear();}

Dot::~Dot() {}

void Dot::Clear()
{
	valid=0; strcpy(type,"Un");
	coor[0]=coor[1]=coor[2]=0.000;
	unit=score=0.000;
}

DotSet::DotSet() {Clear();}

DotSet::~DotSet() {dot.clear();}

DotSet::DotSet(const DotSet &original)
{
 int i,n;

 this->num_dot=original.num_dot;
 this->dot.clear(); n=original.num_dot;
 for(i=0;i<n;i++) this->dot.push_back(original.dot[i]);

 this->r=original.r;
 strcpy(this->type,original.type);
 this->unit=original.unit;
 this->total=original.total;
}

DotSet& DotSet::operator = (const DotSet &original)
{
 if(this==&original) return *this;

 int i,n;

 this->num_dot=original.num_dot;
 this->dot.clear(); n=original.num_dot;
 for(i=0;i<n;i++) this->dot.push_back(original.dot[i]);

 this->r=original.r;
 strcpy(this->type,original.type);
 this->unit=original.unit;
 this->total=original.total;

 return *this;
}

void DotSet::Clear()
{
	r=0.000; unit=0.000; strcpy(type,"Un");
	num_dot=0; dot.clear();
}

void DotSet::Show_Contents() const
{
        int i;

        printf("DOTSET: %4.2f  ", r);
        printf("%3d  ", num_dot);
	printf("%5.3f  ", unit);
	printf("%6.2f  ", total);
        printf("%s\n", type);

        for(i=0;i<num_dot;i++)
        {
         printf("%6.3f %6.3f %6.3f %1d\n", 
		dot[i].coor[0], dot[i].coor[1], dot[i].coor[2], dot[i].valid);
        }

        return;
}

// *****************************************************************************
// output a set of dots to a PDB file
// *****************************************************************************
void DotSet::Show_Dots(char *filename, char *header, char *show) const
{
        FILE *fp;
	int i,num;
	float *property=NULL;

        if((fp=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	fprintf(fp,"HEADER    %s\n", header);
	fprintf(fp,"REMARK Creation time: %s\n", Get_Time());

	num=dot.size();

	property=new float[num];
	if(property==NULL) Memory_Allocation_Error();

	if(!strcasecmp(show,"unit"))
		{
		 for(i=0;i<num;i++) property[i]=dot[i].unit;
		}
	else
		{
		 for(i=0;i<num;i++) property[i]=dot[i].score;
		}

	for(i=0;i<num;i++)
		{
		 if(dot[i].valid==0) continue;
                 else if(!strcmp(dot[i].type,"F"))
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "F");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
                 else if(!strcmp(dot[i].type,"Cl"))
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "CL");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
                 else if(!strcmp(dot[i].type,"Br"))
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "BR");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
                 else if(!strcmp(dot[i].type,"I"))
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "I");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(!strcmp(dot[i].type,"Si"))
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "SI");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(dot[i].type[0]=='C')
			{
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "C");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(dot[i].type[0]=='N')
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "N");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(dot[i].type[0]=='O')
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "O");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(dot[i].type[0]=='S')
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "S");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(dot[i].type[0]=='P')
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "P");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
		 else if(dot[i].type[0]=='H')
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "H");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
                 else
                        {
                         fprintf(fp,"HETATM%5d  %-3s ", i+1, "ZN");
                         fprintf(fp,"DOT %5d    ", dot[i].valid);
                         fprintf(fp,"%8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                 dot[i].coor[0], dot[i].coor[1], 
				 dot[i].coor[2], 1.00, property[i]);
                        }
		}

        fclose(fp); if(property) delete [] property;

	return;
}

Residue::Residue() {Clear();}

Residue::~Residue() {atom.clear(); vol_dot.clear(); sur_dot.clear();}

Residue::Residue(const Residue &original)
{
 this->valid=original.valid;
 strcpy(this->name,original.name);
 this->chain=original.chain;
 strcpy(this->id,original.id);

 int i,n;

 this->num_atom=original.num_atom;
 this->atom.clear(); n=original.num_atom;
 for(i=0;i<n;i++) this->atom.push_back(original.atom[i]);

 this->vol_dot.clear(); n=original.vol_dot.size();
 for(i=0;i<n;i++) this->vol_dot.push_back(original.vol_dot[i]);

 this->sur_dot.clear(); n=original.sur_dot.size();
 for(i=0;i<n;i++) this->sur_dot.push_back(original.sur_dot[i]);
}

Residue& Residue::operator = (const Residue &original)
{
 if(this==&original) return *this;

 this->valid=original.valid;
 strcpy(this->name,original.name);
 this->chain=original.chain;
 strcpy(this->id,original.id);

 int i,n;

 this->num_atom=original.num_atom;
 this->atom.clear(); n=original.num_atom;
 for(i=0;i<n;i++) this->atom.push_back(original.atom[i]);

 this->vol_dot.clear(); n=original.vol_dot.size();
 for(i=0;i<n;i++) this->vol_dot.push_back(original.vol_dot[i]);

 this->sur_dot.clear(); n=original.sur_dot.size();
 for(i=0;i<n;i++) this->sur_dot.push_back(original.sur_dot[i]);

 return *this;
}

void Residue::Clear()
{
	valid=0; strcpy(name,""); chain='\0'; strcpy(id,"0");
	num_atom=0; atom.clear();
	vol_dot.clear(); sur_dot.clear();
	return;
}

void Residue::Show_Contents() const
{
	printf("Residue %s %c %s  valid=%d\n", name,chain,id,valid);
	printf("contains %d atoms:\n", num_atom);

	int i; for(i=0;i<num_atom;i++) atom[i].Show_Contents();

	return;
}

int Residue::Get_Contents_From_Protein(const char *name, const char *id, const char chain)
{
	extern Protein *protein;
	int i;

	Clear();

	for(i=0;i<protein->num_atom;i++)
	{
	 if(protein->atom[i].valid<=0) continue;
	 else if(protein->atom[i].chain!=chain) continue;
	 else if(strcmp(protein->atom[i].residue,name)) continue;
	 else if(strcmp(protein->atom[i].res_id,id)) continue;
	 else this->atom.push_back(protein->atom[i]);
	}

	this->num_atom=this->atom.size();

	if(this->num_atom==0) return FALSE;
	else
	{
	 strcpy(this->name,this->atom[0].residue);
	 strcpy(this->id,this->atom[0].res_id);
	 this->chain=this->atom[0].chain;
	 this->valid=1;

	 return TRUE;
	}
}

// notice, surface dots are generated in context of the protein
void Residue::Generate_Surface_Dots(float probe_r)
{
	extern ForceField *ff;
	extern Protein *protein;
	int i,j,k;
	bool mark;
	float d,dd,dmin;
	DotSet tmp_set;
	Dot tmp_dot;

	sur_dot.clear();

	// first, generate the surface for this residue

	int *protein_check_list;

	protein_check_list=new int[protein->num_atom];
	if(protein_check_list==NULL) Memory_Allocation_Error();	

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].xtype,"H")) continue;

	 // generate the dots

	 tmp_set=ff->Get_Surface_Dot(this->atom[i],probe_r);

	 // determine if the dots are on the surface of this residue

	 for(j=0;j<protein->num_atom;j++)
		{
		 protein_check_list[j]=0;

		 if(protein->atom[j].valid<=0) continue;
		 else if(!strcmp(protein->atom[j].xtype,"H")) continue;
		 else if(strcmp(protein->atom[j].residue,this->name)) continue;
		 else if(strcmp(protein->atom[j].res_id,this->id)) continue;
		 else if(protein->atom[j].chain!=this->chain) continue;
		 else if(protein->atom[j].id==this->atom[i].id) continue;
		 else protein_check_list[j]=1;
		}

	 for(j=0;j<tmp_set.num_dot;j++) 
		{
		 mark=true; dmin=9999.0;

		 for(k=0;k<protein->num_atom;k++)
			{
			 if(protein_check_list[k]==0) continue;

			 d=Distance(tmp_set.dot[j].coor,protein->atom[k].coor);
			 dd=d-(protein->atom[k].R+probe_r);

			 if(dd<0.000) {mark=false; break;} 
			 else if(dd<dmin) {dmin=dd; continue;}
			 else continue; 
			}

		 dd=sqrt(tmp_set.unit);

		 if(mark==false) tmp_set.dot[j].valid=0;
		 else if(dmin>=dd) tmp_set.dot[j].valid=1; // regular dots 
		 else tmp_set.dot[j].valid=2;	           // dots at edge
                }

	 // now determine if the surface dots are exposed 

	 for(j=0;j<protein->num_atom;j++)
		{
		 protein_check_list[j]=0;

		 if(protein->atom[j].valid<=0) continue;
		 else if(!strcmp(protein->atom[j].xtype,"H")) continue;
		 else if(protein->atom[j].id==this->atom[i].id) continue;

		 d=Distance(this->atom[i].coor,protein->atom[j].coor); 

		 if(d>(this->atom[i].R+protein->atom[j].R+2*probe_r)) continue;
		 else protein_check_list[j]=1;
		}

	 for(j=0;j<tmp_set.num_dot;j++)
		{
		 if(tmp_set.dot[j].valid==0) continue;

		 tmp_dot=tmp_set.dot[j];

		 if(tmp_dot.valid==1) tmp_dot.unit=tmp_set.unit;
		 else tmp_dot.unit=tmp_set.unit*0.500;

		 mark=true;

		 for(k=0;k<protein->num_atom;k++)
			{
			 if(protein_check_list[k]==0) continue;

			 d=Distance(tmp_dot.coor,protein->atom[k].coor);
			 dd=d-(protein->atom[k].R+probe_r);

			 if(dd<0.000) {mark=false; break;} 
			 else continue; 
			}

		 if(mark==false) tmp_dot.valid=0;
		 else tmp_dot.valid=1; 

		 this->sur_dot.push_back(tmp_dot);
		}
	}

	if(protein_check_list) delete [] protein_check_list;

	return;
}

Chain::Chain() {Clear();}

Chain::~Chain() {residue.clear();}

Chain::Chain(const Chain &original)
{
 this->valid=original.valid;
 this->label=original.label;
 this->length=original.length;

 int i,n;

 this->residue.clear(); n=original.length;
 for(i=0;i<n;i++) this->residue.push_back(original.residue[i]);
}

Chain& Chain::operator = (const Chain &original)
{
 if(this==&original) return *this;

 this->valid=original.valid;
 this->label=original.label;
 this->length=original.length;

 int i,n;

 this->residue.clear(); n=original.length;
 for(i=0;i<n;i++) this->residue.push_back(original.residue[i]);

 return *this;
}

void Chain::Clear()
{
	valid=0; label='\0';
	length=0; residue.clear(); 
	return;
}

void Chain::Show_Contents() const
{
	printf("Chain %c: valid=%d\n", label, valid);
	printf("contains %d residues:\n", length);

	int i; for(i=0;i<length;i++) residue[i].Show_Contents();

	return;
}

