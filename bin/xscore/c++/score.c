# include "xtool.h"

int Xtool_Score()
{
	extern Input *input;
        extern ForceField *ff;
        extern Protein *protein;
        extern Ligand *ligand;
	Ligand reference,cofactor;

	printf("Now reading parameters from '%s' ...\n",input->parameter_dir);
        ff=new ForceField(input->parameter_dir);
        if(ff==NULL) Memory_Allocation_Error();

	if(strcmp(input->reference_file,"none"))  // reference molecule given
	{
	 printf("Now reading the reference ligand from '%s' ... \n", input->reference_file);
	 reference.Read_From_Mol2(input->reference_file);
	 reference.Value_Atom();
	}
	else  // use the ligand molecule as the reference
	{
	 printf("Now reading the reference ligand from '%s' ... \n", input->ligand_file);
	 reference.Read_From_Mol2(input->ligand_file);
	 reference.Value_Atom();
	}
	
        printf("Now reading the protein from '%s' ... \n",input->receptor_file);
	protein=new Protein; 
	if(protein==NULL) Memory_Allocation_Error();
        protein->Read_From_PDB(input->receptor_file);
        protein->Value_Atom();

	if(strcmp(input->cofactor_file,"none"))
	{
	 printf("Now reading the cofactor from '%s' ... \n",input->cofactor_file);
	 cofactor.Read_From_Mol2(input->cofactor_file);
	 cofactor.Value_Atom();
	 protein->Merge_Cofactor(cofactor);
	}

	printf("Now defining the binding pocket ...\n");
	protein->Define_Pocket(&reference,10.0);

	printf("Now scoring the ligands in '%s'...\n", input->ligand_file);
	Population population(input->ligand_file);
	population.Score_Members();
	population.Rank_Members();
	population.Output_Results();

	return TRUE;
}

int Xtool_Score_Shortcut()
{
        extern Input *input;
        extern ForceField *ff;
        extern Protein *protein;
        extern Ligand *ligand;
	int i,num_molecule;
	FILE *fin,*fout;
	Ligand cofactor;

	num_molecule=Check_Mol2_File(input->ligand_file);

	if(num_molecule==0)
        {
	 puts("Error: no valid ligand molecule is given.");
	 exit(1);
        }

	ligand=new Ligand; if(ligand==NULL) Memory_Allocation_Error();
	protein=new Protein; if(protein==NULL) Memory_Allocation_Error();

	printf("Now reading parameters from '%s' ...\n",input->parameter_dir);
        ff=new ForceField(input->parameter_dir);
        if(ff==NULL) Memory_Allocation_Error();

	printf("Now reading the ligand from '%s' ...\n", input->ligand_file);
	ligand->Read_From_Mol2(input->ligand_file);
	ligand->Value_Atom();
	// ligand->Write_Out_Mol2("temp.mol2");

	printf("Now reading the protein from '%s' ... \n",input->receptor_file);
	protein->Read_From_PDB(input->receptor_file);
	protein->Value_Atom();	

	if(strcmp(input->cofactor_file,"none"))
        {
         printf("Now reading the cofactor from '%s' ... \n",input->cofactor_file);
         cofactor.Read_From_Mol2(input->cofactor_file);
         cofactor.Value_Atom();
         protein->Merge_Cofactor(cofactor);
        }

        printf("Now defining the binding pocket ...\n");
        protein->Define_Pocket(ligand,10.0);

	if(num_molecule==1)	// single mol2 file
	{
	 ligand->Calculate_Binding_Score();

	 printf("\n***********************************************\n");

	 printf("HPScore -log(Kd) = %-5.2f\n", ligand->pkd1);
	 printf("HMScore -log(Kd) = %-5.2f\n", ligand->pkd2);
	 printf("HSScore -log(Kd) = %-5.2f\n", ligand->pkd3);
	 printf("Predicted average -log(Kd) = %-5.2f\n", ligand->bind_score);

	 printf("Predicted binding energy = %-6.2f kcal/mol\n", 
		ligand->bind_score*(-1.364)); 

	 printf("***********************************************\n\n");

	 printf("Now writing scoring information to '%s' ...\n", input->log_file);
	 ligand->Write_Out_Log(input->log_file);
	}

	else	// multi-mol2 file
	{	
	 printf("ID: HPSCORE HMSCORE HSSCORE AVE_SCORE BIND_ENERGY NAME\n");

	 if((fin=fopen(input->ligand_file,"r"))==NULL) Open_File_Error(input->ligand_file);

	 if((fout=fopen(input->log_file,"w"))==NULL) Open_File_Error(input->log_file);

	 for(i=1;i<=num_molecule;i++)
		{
		 ligand->Clear();
		 ligand->Read_From_Mol2(fin);
        	 ligand->Value_Atom();
		 ligand->Calculate_Binding_Score();

		 printf("Molecule %6d: ", i);
		 printf("%5.2f  ", ligand->pkd1);
		 printf("%5.2f  ", ligand->pkd2);
		 printf("%5.2f  ", ligand->pkd3);
		 printf("%5.2f  ", ligand->bind_score);
		 printf("%6.2f  ", ligand->bind_score*(-1.364));
		 printf("%s\n", ligand->name);

		 ligand->Write_Out_Log(fout);
		}

	 fclose(fin); fclose(fout);
	}

	return TRUE;
}

float Ligand::Calculate_Binding_Score()
{
	extern Input *input;
	int i;
	Ligand bak;

	if(input->num_method==0) return 0.000;
	else bak=(*this);  // backup the atomic charge information

	// these scoring functions may need formal charges

	this->Assign_Apparent_Charge();

	// prepare for atomic binding score

	if(this->abs_inf) delete [] this->abs_inf; 
	this->abs_inf=new ABS[this->num_atom];
	if(abs_inf==NULL) Memory_Allocation_Error();

	for(i=0;i<this->num_atom;i++)
		{
		 abs_inf[i].pkd1=0.000;
		 abs_inf[i].pkd2=0.000;
		 abs_inf[i].pkd3=0.000;
		 abs_inf[i].vdw=0.000;
		 abs_inf[i].hb=0.000;
		 abs_inf[i].hm=0.000;
		 abs_inf[i].hp=0.000;
		 abs_inf[i].hs=0.000;
		 abs_inf[i].rt=0.000;
		 abs_inf[i].score=0.000;
		}

	// clean other variables

	this->pkd1=this->pkd2=this->pkd3=0.000;
	this->vdw=this->hb=this->hm=this->hp=this->hs=this->rt=0.000;

	// first, generate dot surface

	this->Generate_Surface_Dots(WATER_R);

	// now calculate the van der Waals term

	this->vdw=this->Calculate_VDW();
	for(i=0;i<this->num_atom;i++) 
		{
		 abs_inf[i].vdw=this->atom[i].score;
		}

	// now calculate the H-bond term

	this->hb=this->Calculate_HB();
	for(i=0;i<this->num_atom;i++) 
		{
		 abs_inf[i].hb=this->atom[i].score;
		}

	// now calculate the rotor term

	this->rt=this->Calculate_RT();
	for(i=0;i<this->num_atom;i++) 
		{
		 abs_inf[i].rt=this->atom[i].score;
		}

	// now calculate the hydrophobic term

	if(!strcmp(input->apply_hpscore,"YES")) 
	{
	 this->hp=this->Calculate_HP();

	 for(i=0;i<this->num_atom;i++) 
		{
		 abs_inf[i].hp=this->atom[i].score;
		}

	 this->pkd1=input->hpscore_c0+
                    input->hpscore_cvdw*this->vdw+ 
                    input->hpscore_chb*this->hb+ 
                    input->hpscore_chp*this->hp+
                    input->hpscore_crt*this->rt;
	}

	if(!strcmp(input->apply_hmscore,"YES")) 
	{
	 this->hm=this->Calculate_HM();

	 for(i=0;i<this->num_atom;i++) 
		{
		 abs_inf[i].hm=this->atom[i].score;
		}

	 this->pkd2=input->hmscore_c0+
                    input->hmscore_cvdw*this->vdw+ 
                    input->hmscore_chb*this->hb+ 
                    input->hmscore_chm*this->hm+
              	    input->hmscore_crt*this->rt;

	}

	if(!strcmp(input->apply_hsscore,"YES")) 
	{
	 this->hs=this->Calculate_HS();

	 for(i=0;i<this->num_atom;i++) 
		{
		 abs_inf[i].hs=this->atom[i].score;
		}

	 this->pkd3=input->hsscore_c0+
                    input->hsscore_cvdw*this->vdw+ 
                    input->hsscore_chb*this->hb+ 
                    input->hsscore_chs*this->hs+
                    input->hsscore_crt*this->rt;
	}

	this->bind_score=(pkd1+pkd2+pkd3)/(input->num_method);

	// now compute atomic binding scores

	int num_nonh=0; num_nonh=this->Get_Num_Heavy_Atom(); 

	if(!strcmp(input->apply_hpscore,"YES")) 
		{
		 for(i=0;i<this->num_atom;i++)
			{
			 if(this->atom[i].valid<=0) continue;
			 else if(!strcmp(this->atom[i].type,"H")) continue;

	 	 	 abs_inf[i].pkd1=input->hpscore_c0/num_nonh;
			 abs_inf[i].pkd1+=input->hpscore_cvdw*abs_inf[i].vdw;
			 abs_inf[i].pkd1+=input->hpscore_chb*abs_inf[i].hb;
			 abs_inf[i].pkd1+=input->hpscore_chp*abs_inf[i].hp;
			 abs_inf[i].pkd1+=input->hpscore_crt*abs_inf[i].rt;
			}
		}

	if(!strcmp(input->apply_hmscore,"YES")) 
		{
		 for(i=0;i<this->num_atom;i++)
			{
			 if(this->atom[i].valid<=0) continue;
			 else if(!strcmp(this->atom[i].type,"H")) continue;

	 	 	 abs_inf[i].pkd2=input->hmscore_c0/num_nonh;
			 abs_inf[i].pkd2+=input->hmscore_cvdw*abs_inf[i].vdw;
			 abs_inf[i].pkd2+=input->hmscore_chb*abs_inf[i].hb;
			 abs_inf[i].pkd2+=input->hmscore_chm*abs_inf[i].hm;
			 abs_inf[i].pkd2+=input->hmscore_crt*abs_inf[i].rt;
			}
		}

	if(!strcmp(input->apply_hsscore,"YES")) 
		{
		 for(i=0;i<this->num_atom;i++)
			{
			 if(this->atom[i].valid<=0) continue;
			 else if(!strcmp(this->atom[i].type,"H")) continue;

	 	 	 abs_inf[i].pkd3=input->hsscore_c0/num_nonh;
			 abs_inf[i].pkd3+=input->hsscore_cvdw*abs_inf[i].vdw;
			 abs_inf[i].pkd3+=input->hsscore_chb*abs_inf[i].hb;
			 abs_inf[i].pkd3+=input->hsscore_chs*abs_inf[i].hs;
			 abs_inf[i].pkd3+=input->hsscore_crt*abs_inf[i].rt;
			}
		}

	float sum=0.000;

	for(i=0;i<this->num_atom;i++)
		{
		 abs_inf[i].score=0.000;
		 abs_inf[i].score+=abs_inf[i].pkd1;
		 abs_inf[i].score+=abs_inf[i].pkd2;
		 abs_inf[i].score+=abs_inf[i].pkd3;
		 abs_inf[i].score/=input->num_method;
		 this->atom[i].score=abs_inf[i].score;
		 sum+=this->atom[i].score;
		}

/*
	// check if anything wrong in atomic binding score calculation

	if(fabs(this->bind_score-sum)>0.01)
		{
	 	 puts("Warning: sum of ABS != total binding score.");
		}
*/

	if(!strcmp(input->show_abs,"YES"))  // use abs as charges
		{
	 	 strcpy(this->charge_type,"USER_CHARGES");

	 	 for(i=0;i<this->num_atom;i++) 
			{
			 this->atom[i].q=abs_inf[i].score;
			}
		}
	else		// restore original charges
		{
	 	 strcpy(this->charge_type, bak.charge_type);

	 	 for(i=0;i<this->num_atom;i++)
			{
		 	 this->atom[i].q=bak.atom[i].q;
			}
		}

	return this->bind_score;
}

// **************************************************************************
// Computation of the VDW term in X-Score v1.2
// Latest update: 11/21/2003
// **************************************************************************
float Ligand::Calculate_VDW()
{
	extern Protein *protein;
	int i,j;
	float cutoff,d0,d,tmp,tmp1,tmp2,sum,asum;

	// clear the variables

	cutoff=DIST_CUTOFF; sum=0.000;
	for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

	// now calculate the P-L vdw interaction 

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].type,"H")) continue;

	 asum=0.000;

	 for(j=0;j<protein->num_atom;j++)
		{
		 if(protein->atom[j].valid!=2) continue;
		 else if(!strcmp(protein->atom[j].type,"H")) continue;
		 else if(!strcmp(protein->atom[j].type,"O.w")) continue;

		 d0=this->atom[i].R+protein->atom[j].R;
		 d=Distance(this->atom[i].coor,protein->atom[j].coor);

		 if(d>cutoff) continue;

		 tmp1=d0/d; 
		 tmp1=tmp1*tmp1*tmp1*tmp1; tmp2=tmp1*tmp1;
		 tmp=tmp2-2.0*tmp1;

		 asum+=tmp;
		}

	 // revert the sign

	 asum*=(-1.00);

	 // if this atom has unfavorable contribution then neglect it

	 if(asum<0.00) continue;
	 else {sum+=asum; this->atom[i].score=asum;}
	}

	return sum;
}

// **************************************************************************
// Computation of the HB term in X-Score v1.2 (including M-bonds)
// Latest update 11/21/2003
// **************************************************************************
float Ligand::Calculate_HB(char *pdb_entry) 
{
	int i,j;
	float sum;
	HBond candidate[1000],temp;
	int num_candidate;

	// clear the variables

	sum=0.000; num_candidate=0;
	for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

	// first, get the HB pairs between protein and ligand

	num_candidate=this->Get_HBond_Pair_PL(candidate);
	this->Sum_HBonds(num_candidate,candidate);

	// then sum their contributions 

	for(i=0;i<num_candidate;i++)
	{
	 if(fabs(candidate[i].score<0.01)) continue;
	 else
		{
	 	 sum+=candidate[i].score;
		 this->atom[candidate[i].latom-1].score+=candidate[i].score;
		}
	}

	// output H-bonds for analysis

	if(strcmp(pdb_entry,"unknown"))
	{
	 // arrange them according to ligand atom id 

	 for(i=0;i<num_candidate-1;i++)
	 for(j=i+1;j<num_candidate;j++)
		{
		 if(candidate[i].latom<=candidate[j].latom) continue;
		 else {SWAP(candidate[i],candidate[j]);}
		}

	 for(i=0;i<num_candidate;i++)
		{
		 if(fabs(candidate[i].score)<=0.000) continue;
		 else candidate[i].Show_Contents();
		}
	}

	return sum;
}

void Ligand::Sum_HBonds(int num_candidate, HBond candidate[]) const
{
 extern Protein *protein;
 int i,j,k,count,limit,latom,patom;
 float angle,v1[3],v2[3];
 HBond temp;
 Group tmp_group;

 // Step 1: rank candidates according their strength in decreasing order
 // note "fabs" is applied because SB strength could be negative

 for(i=0;i<num_candidate-1;i++)
 for(j=i+1;j<num_candidate;j++)
        {
         if(fabs(candidate[i].score)>=fabs(candidate[j].score)) continue;
         else {SWAP(candidate[i],candidate[j]);}
        }

 // Step 2: check the angular limit: the angle between any two H-bonds 
 // on the same atom must be larger than 45 degrees
 // note this this filter could be applied to both ligand and protein 

 for(i=0;i<num_candidate-1;i++)
 for(j=i+1;j<num_candidate;j++)
        {
         if(candidate[i].latom!=candidate[j].latom) continue;

	 for(k=0;k<3;k++) 
	{
	 latom=candidate[i].latom; patom=candidate[i].patom;
	 v1[k]=protein->atom[patom-1].coor[k]-this->atom[latom-1].coor[k];
	 latom=candidate[j].latom; patom=candidate[j].patom;
	 v2[k]=protein->atom[patom-1].coor[k]-this->atom[latom-1].coor[k];
	}

         angle=fabs(Angle_of_Two_Vectors(v1,v2));

         if(angle<45.0) candidate[j].score=0.000;
         else continue;
        }

/*
 for(i=0;i<num_candidate-1;i++)
 for(j=i+1;j<num_candidate;j++)
        {
         if(candidate[i].patom!=candidate[j].patom) continue;

         for(k=0;k<3;k++)
        {
	 latom=candidate[i].latom; patom=candidate[i].patom;
	 v1[k]=this->atom[latom-1].coor[k]-protein->atom[patom-1].coor[k];
	 latom=candidate[j].latom; patom=candidate[j].patom;
	 v2[k]=this->atom[latom-1].coor[k]-protein->atom[patom-1].coor[k];
        }

         angle=fabs(Angle_of_Two_Vectors(v1,v2));

         if(angle<45.0) candidate[j].score=0.000;
         else continue;
        }
*/

 // Step 3: an donor atom shall not form more H-bonds than the H atoms it has.
 // note that this filter is applied only to the ligand side

 for(i=0;i<num_candidate-1;i++)
	{
	 if(candidate[i].type!=1) continue;

	 count=1; latom=candidate[i].latom;

	 if(!strcmp(this->atom[latom-1].hb,"DA")) 
		{
		 limit=1; 
		}
	 else
		{
		 tmp_group=this->Find_A_Group(latom); 
		 limit=tmp_group.num_h;
		}

	 for(j=i+1;j<num_candidate;j++)
		{
		 if(candidate[j].type!=1) continue;
		 else if(candidate[i].latom!=candidate[j].latom) continue;

		 count++;

		 if(count<=limit) continue;
		 else candidate[j].score=0.000;
		}
	}

 // Step 4: an acceptor atom shall not form more H-bonds than its LPs 
 // note that this filter is applied only to the ligand side

 for(i=0;i<num_candidate-1;i++)
	{
	 if(candidate[i].type!=2&&candidate[i].type!=3) continue;

	 count=1; latom=candidate[i].latom;

	 if(this->atom[latom-1].type[0]=='O') limit=2;
	 else if(this->atom[latom-1].type[0]=='N') limit=1;
	 else if(this->atom[latom-1].type[0]=='S') limit=2;

	 for(j=i+1;j<num_candidate;j++)
		{
		 if(candidate[j].type!=2&&candidate[j].type!=3) continue;
		 else if(candidate[i].latom!=candidate[j].latom) continue;

		 count++;

		 if(count<=limit) continue;
		 else candidate[j].score=0.000;
		}
	}

 return;
}

// ************************************************************************
// Get all of the potential H-bond pairs between protein and ligand
// sb_flag determines if SB should be treated differently from HB
// Latest update: 11/21/2003
// ************************************************************************
int Ligand::Get_HBond_Pair_PL(HBond candidate[], bool sb_flag) const
{
 extern Protein *protein;
 int i,j,num,type;
 bool sb;
 float d,cutoff;
 HBond tmp_candidate;

 num=0; cutoff=5.00;

 // note that the following H-bond algorithm is not based on any
 // explicit hydrogen atom at all, and this is the right thing to do.

 for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].type,"H")) continue;
	 else if(!strcmp(this->atom[i].hb,"N")) continue;
	 else if(!strcmp(this->atom[i].hb,"H")) continue;
	 else if(!strcmp(this->atom[i].hb,"P")) continue;
	 else if(!strcmp(this->atom[i].hb,"DH")) continue;

	 for(j=0;j<protein->num_atom;j++)
		{
		 if(protein->atom[j].valid!=2) continue;
		 else if(!strcmp(protein->atom[j].type,"H")) continue;
		 else if(!strcmp(protein->atom[j].type,"O.w")) continue;
		 else if(!strcmp(protein->atom[j].hb,"H")) continue;
		 else if(!strcmp(protein->atom[j].hb,"P")) continue;
		 else if(!strcmp(protein->atom[j].hb,"N")) continue;

		 // determine the type of this H-bond first
		 // type=0, no H-bond
		 // type=1, latom is the donor, patom is the acceptor;
		 // type=2, latom is the acceptor, patom is the donor;
		 // type=3, latom bound with metal ion; 

		 if(!strcmp(this->atom[i].hb,"D"))
		 	{
		 	 if(!strcmp(protein->atom[j].hb,"D")) type=0;
		 	 else if(!strcmp(protein->atom[j].hb,"A")) type=1;
		 	 else if(!strcmp(protein->atom[j].hb,"DA")) type=1;
		 	 else if(!strcmp(protein->atom[j].hb,"M")) type=0;
		 	 else type=0;
			}
		 else if(!strcmp(this->atom[i].hb,"A"))
			{
		 	 if(!strcmp(protein->atom[j].hb,"D")) type=2;
		 	 else if(!strcmp(protein->atom[j].hb,"A")) type=0;
                 	 else if(!strcmp(protein->atom[j].hb,"DA")) type=2;
		 	 else if(!strcmp(protein->atom[j].hb,"M")) type=3;
                 	 else type=0;
			}
		 else if(!strcmp(this->atom[i].hb,"DA"))
			{
		 	 if(!strcmp(protein->atom[j].hb,"A")) type=1;
		 	 else if(!strcmp(protein->atom[j].hb,"D")) type=2;
		 	 else if(!strcmp(protein->atom[j].hb,"DA")) type=1;
		 	 else if(!strcmp(protein->atom[j].hb,"M")) type=3;
		 	 else type=0;
			}
		 else type=0;

		 if(type==0) continue;  // no H-bond

		 // a crude distance check

		 d=Distance(this->atom[i].coor,protein->atom[j].coor);
		 if(d>cutoff) continue;

		 // this section is used to differentiate HB and SB

		 if(!strcmp(this->atom[i].type,"O.co2")&&
		    !strcmp(protein->atom[j].type,"O.co2")) 
			{
			 sb=false;
			}
	 	 else if((fabs(this->atom[i].q)>0.01)&&
		         (fabs(protein->atom[j].q)>0.01)) 
			{
			 sb=true; 
			}
		 else 
			{
			 sb=false;
			}

		 // now handle this h-bond

		 tmp_candidate.Clear();
		 tmp_candidate.type=type;
		 tmp_candidate.sb=sb;
		 tmp_candidate.latom=i+1;
		 tmp_candidate.patom=j+1;

		 if(type==1)
			{
			 tmp_candidate.D=this->atom[i];
			 tmp_candidate.A=protein->atom[j];
			}
		 else if(type==2)
			{
			 tmp_candidate.A=this->atom[i];
			 tmp_candidate.D=protein->atom[j];
			}
		 else if(type==3)
			{
			 tmp_candidate.A=this->atom[i];
			 tmp_candidate.D=protein->atom[j];
			}

		 // if need to treat charged H-bonds differently

		 if(sb_flag==true&&sb==true) tmp_candidate.Value_SBond();  
		 else tmp_candidate.Value_HBond_2();

		 if(fabs(tmp_candidate.score)>0.000)
			{
			 candidate[num]=tmp_candidate; num++;
			}
		 else continue;
		}
	}

 return num;
}

// ************************************************************************
// Get all of the potential H-bond pairs inside the protein 
// Latest update: 11/18/2003
// ************************************************************************
int Ligand::Get_HBond_Pair_PP(HBond candidate[]) const
{
 extern Protein *protein;
 int i,j,num,type;
 HBond tmp_candidate;

 num=0;

 // note that the following H-bond algorithm is not based on any
 // hydrogen atom at all, and this is the right thing to do.

 for(i=0;i<protein->num_atom-1;i++)
	{
	 if(protein->atom[i].valid<2) continue;
	 else if(!strcmp(protein->atom[i].type,"H")) continue;
	 else if(!strcmp(protein->atom[i].type,"O.w")) continue;
	 else if(!strcmp(protein->atom[i].hb,"H")) continue;
	 else if(!strcmp(protein->atom[i].hb,"P")) continue;
	 else if(!strcmp(protein->atom[i].hb,"N")) continue;
	 else if(!strcmp(protein->atom[i].hb,"DH")) continue;

	 for(j=i+1;j<protein->num_atom;j++)
		{
		 if(protein->atom[j].valid<2) continue;
		 else if(!strcmp(protein->atom[j].type,"H")) continue;
		 else if(!strcmp(protein->atom[j].type,"O.w")) continue;
		 else if(!strcmp(protein->atom[j].hb,"H")) continue;
		 else if(!strcmp(protein->atom[j].hb,"P")) continue;
		 else if(!strcmp(protein->atom[j].hb,"N")) continue;
		 else if(!strcmp(protein->atom[j].hb,"DH")) continue;

		 // atoms in the same residue are not supposed to form H-bond

		 if(!strcmp(protein->atom[i].res_id,protein->atom[j].res_id)&&
		    !strcmp(protein->atom[i].residue,protein->atom[j].residue))
			 continue;

		 // type=0, no H-bond
		 // type=1, atom i is the donor, atom j is the acceptor;
		 // type=2, atom i is the acceptor, atom j is the donor;

		 if(!strcmp(protein->atom[i].hb,"M"))
			{
			 if(!strcmp(protein->atom[j].hb,"D")) type=0;
			 else if(!strcmp(protein->atom[j].hb,"A")) type=1;
			 else if(!strcmp(protein->atom[j].hb,"DA")) type=1;
			 else if(!strcmp(protein->atom[j].hb,"M")) type=0;
			 else type=0;
			}
		 else if(!strcmp(protein->atom[i].hb,"D"))
		 	{
		 	 if(!strcmp(protein->atom[j].hb,"D")) type=0;
		 	 else if(!strcmp(protein->atom[j].hb,"A")) type=1;
		 	 else if(!strcmp(protein->atom[j].hb,"DA")) type=1;
		 	 else if(!strcmp(protein->atom[j].hb,"M")) type=0;
		 	 else type=0;
			}
		 else if(!strcmp(protein->atom[i].hb,"A"))
			{
		 	 if(!strcmp(protein->atom[j].hb,"D")) type=2;
		 	 else if(!strcmp(protein->atom[j].hb,"A")) type=0;
                 	 else if(!strcmp(protein->atom[j].hb,"DA")) type=2;
                 	 else if(!strcmp(protein->atom[j].hb,"M")) type=2;
                 	 else type=0;
			}
		 else if(!strcmp(protein->atom[i].hb,"DA"))
			{
		 	 if(!strcmp(protein->atom[j].hb,"A")) type=1;
		 	 else if(!strcmp(protein->atom[j].hb,"D")) type=2;
		 	 else if(!strcmp(protein->atom[j].hb,"DA")) type=1;
		 	 else if(!strcmp(protein->atom[j].hb,"M")) type=2;
		 	 else type=0;
			}
		 else type=0;

		 if(type==0) continue;

		 // now handle this h-bond

		 tmp_candidate.Clear();
		 tmp_candidate.type=5;
		 tmp_candidate.sb=false;

		 if(type==1)
			{
			 tmp_candidate.D=protein->atom[i];
			 tmp_candidate.A=protein->atom[j];
			}
		 else if(type==2)
			{
			 tmp_candidate.A=protein->atom[i];
			 tmp_candidate.D=protein->atom[j];
			}

		 tmp_candidate.Value_HBond_2();

		 if(tmp_candidate.score>0.000)
			{
			 candidate[num]=tmp_candidate; num++;
			}
		 else continue;
		}
	}

 return num;
}

// *************************************************************************
// Computation of the HP term in X-Score v1.2
// Latest update: 11/21/2003
// *************************************************************************
float Ligand::Calculate_HP()
{
	extern Protein *protein;
	int i,j;
	float d,d1,d2,cutoff;
	float tmp,asum,sum;

	sum=0.000;
	for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

	cutoff=DIST_CUTOFF;
	
	for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid<=0) continue;
		 else if(strcmp(this->atom[i].hb,"H")) continue;

		 asum=0.000;

		 for(j=0;j<protein->num_atom;j++)
			{
			 if(protein->atom[j].valid!=2) continue;
			 else if(!strcmp(protein->atom[j].type,"H")) continue;
			 else if(!strcmp(protein->atom[j].type,"O.w")) continue;
			 else if(strcmp(protein->atom[j].hb,"H")) continue;

			 d=Distance(this->atom[i].coor, protein->atom[j].coor);
			 if(d>=cutoff) continue;

			 d1=this->atom[i].R+protein->atom[j].R+0.500;
			 d2=this->atom[i].R+protein->atom[j].R+2.200;

                         if(d<d1) tmp=1.000;
                         else if(d<d2) tmp=(1/(d1-d2))*(d-d2);
                         else tmp=0.000;

                         asum+=tmp;
			}

		 this->atom[i].score=asum; sum+=asum;
		}

	return sum;
}

// ***************************************************************************
// Computation of the HM term in X-Score v1.2
// Latest update: 11/21/2002
// ***************************************************************************
float Ligand::Calculate_HM()
{
 extern Protein *protein;
 int i,j;
 float asum,sum,tmp,d,d1,d2,total,cutoff;

 for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

 sum=0.000; cutoff=DIST_CUTOFF;

 for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].type,"H")) continue;
	 else if(strcmp(this->atom[i].hb,"H")) continue;
	 else if(this->atom[i].logp<=0.00) continue;

	 total=0.000; asum=0.000;

	 for(j=0;j<protein->num_atom;j++)
		{
		 if(protein->atom[j].valid!=2) continue;
		 else if(!strcmp(protein->atom[j].type,"H")) continue;
		 else if(!strcmp(protein->atom[j].type,"O.w")) continue;

		 d=Distance(this->atom[i].coor, protein->atom[j].coor);
		 if(d>cutoff) continue;

		 d1=this->atom[i].R+protein->atom[j].R+0.50;
		 d2=this->atom[i].R+protein->atom[j].R+2.20;

		 if(d<d1) tmp=1.000; 
                 else if(d<d2) tmp=(1/(d1-d2))*(d-d2);
                 else tmp=0.000;
			
		 total+=(protein->atom[j].logp*tmp);
		}

	 if(this->atom[i].logp>=0.50) asum=this->atom[i].logp;
 	 else if(total>-0.50) asum=this->atom[i].logp;
	 else asum=0.000;

	 this->atom[i].score=asum; sum+=asum;
	}

 return sum;
}

// ***********************************************************
// Computation of the HS term in X-Score v1.2 
// Latest update: 11/21/2003
// ***********************************************************
float Ligand::Calculate_HS()
{
 int i;
 float sum,total,buried;

 // clear the variables first

 sum=0.000;
 for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

 // then get the buried surface atom by atom

 for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].type,"H")) continue;
	 else if(!strcmp(this->atom[i].hb,"DH")) continue;
	 else if(!strcmp(this->atom[i].hb,"D")) continue;
	 else if(!strcmp(this->atom[i].hb,"A")) continue;
	 else if(!strcmp(this->atom[i].hb,"DA")) continue;
	 else if(!strcmp(this->atom[i].hb,"P")) continue;
	 else if(!strcmp(this->atom[i].hb,"N")) continue;

	 this->Atom_Buried_Surface(i+1,total,buried);

	 sum+=buried; this->atom[i].score=buried;

	 // sum+=(buried*this->atom[i].logp);
	 // this option works *slightly* better
	}

 return sum;
}

// ****************************************************************************
// Computation of the RT term in X-Score v1.2
// Latest update: 11/21/2003
// ****************************************************************************
float Ligand::Calculate_RT()
{
	int i,j,mark,tmp;
	int id1,id2;
	float sum;

	// if a single bond is normal, bond.valid=1;
	// if a single bond is a rotor, bond.valid=2;

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

	// notice: the following part is different from Molecule::Count_Rotor()

	sum=0.000;

	for(i=0;i<num_atom;i++)
		{
		 if(atom[i].valid<=0) continue;

		 mark=0;

		 for(j=0;j<atom[i].num_neib;j++)
			{
			 tmp=atom[i].bond[j];
			 if(tmp==0) continue;
			 else if(bond[tmp-1].valid!=2) continue;
			 else mark++;
			}

		 if(mark==1) this->atom[i].score+=0.50;
		 else if(mark==2) this->atom[i].score+=1.00;
		 else if(mark>=3) this->atom[i].score+=0.50;

		 sum+=this->atom[i].score;
		}

	return sum;
}

// ******************************************************************
// calculate unpaired H-bonding pairs between the complex 
// latest update 11/21/2003
// ******************************************************************
float Ligand::Calculate_UHB() const 
{
	extern Protein *protein;
	int i,j;
	float sum,d,d1,d2,d0,tmp;

	// **********************************************
	// now count the unpaired HB atoms on the ligand
	// **********************************************

	sum=0.000;

	// first, count the pairs between the HB atoms on the ligand 
  	// and the H/P atoms on the protein

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].type,"H")) continue;
	 else if(!strcmp(this->atom[i].hb,"N")) continue;
	 else if(!strcmp(this->atom[i].hb,"H")) continue;
	 else if(!strcmp(this->atom[i].hb,"P")) continue;

	 for(j=0;j<protein->num_atom;j++)
		{
		 if(protein->atom[j].valid<2) continue;
		 else if(!strcmp(protein->atom[j].type,"H")) continue;
		 else if(!strcmp(protein->atom[j].type,"O.w")) continue;
		 else if(!strcmp(protein->atom[j].hb,"N")) continue;
		 else if(!strcmp(protein->atom[j].hb,"DH")) continue;
		 else if(!strcmp(protein->atom[j].hb,"D")) continue;
		 else if(!strcmp(protein->atom[j].hb,"A")) continue;
		 else if(!strcmp(protein->atom[j].hb,"DA")) continue;

		 d0=this->atom[i].R+protein->atom[j].R;
		 d1=d0; d2=d0+0.50;
		 d=Distance(this->atom[i].coor,protein->atom[j].coor);

		 if(d<d1) tmp=1.000;
		 else if(d<d2) tmp=(1/(d1-d2))*(d-d2);
		 else tmp=0.000;

		 sum+=tmp;
		}
	}

	// second, count the pairs between the H/P atoms on the ligand
	// and the HB atoms on the protein

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].type,"H")) continue;
	 else if(!strcmp(this->atom[i].hb,"DH")) continue;
	 else if(!strcmp(this->atom[i].hb,"D")) continue;
	 else if(!strcmp(this->atom[i].hb,"A")) continue;
	 else if(!strcmp(this->atom[i].hb,"DA")) continue;

	 for(j=0;j<protein->num_atom;j++)
		{
		 if(protein->atom[j].valid<2) continue;
		 else if(!strcmp(protein->atom[j].type,"H")) continue;
		 else if(!strcmp(protein->atom[j].type,"O.w")) continue;
		 else if(!strcmp(protein->atom[j].hb,"DH")) continue;
		 else if(!strcmp(protein->atom[j].hb,"N")) continue;
		 else if(!strcmp(protein->atom[j].hb,"H")) continue;
		 else if(!strcmp(protein->atom[j].hb,"P")) continue;

		 d0=this->atom[i].R+protein->atom[j].R;
		 d1=d0; d2=d0+0.50;
		 d=Distance(this->atom[i].coor,protein->atom[j].coor);

		 if(d<d1) tmp=1.000;
		 else if(d<d2) tmp=(1/(d1-d2))*(d-d2);
		 else tmp=0.000;

		 sum+=tmp;
		}
	}

	// third, count the pairs between the HB atoms on the ligand
	// and the HB atoms on the protein

	for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;
	 else if(!strcmp(this->atom[i].type,"H")) continue;
	 else if(!strcmp(this->atom[i].hb,"N")) continue;
	 else if(!strcmp(this->atom[i].hb,"H")) continue;
	 else if(!strcmp(this->atom[i].hb,"P")) continue;

	 for(j=0;j<protein->num_atom;j++)
		{
		 if(protein->atom[j].valid<2) continue;
		 else if(!strcmp(protein->atom[j].type,"H")) continue;
		 else if(!strcmp(protein->atom[j].type,"O.w")) continue;
		 else if(!strcmp(protein->atom[j].hb,"N")) continue;
		 else if(!strcmp(protein->atom[j].hb,"H")) continue;
		 else if(!strcmp(protein->atom[j].hb,"P")) continue;

		 d0=this->atom[i].R+protein->atom[j].R;
		 d1=d0; d2=d0+0.500;
		 d=Distance(this->atom[i].coor,protein->atom[j].coor);

		 if(d<d1) tmp=1.000;
		 else if(d<d2) tmp=(1/(d1-d2))*(d-d2);
		 else tmp=0.000;

		 sum+=tmp;
		}
	}

	return sum;
}

// ***************************************************************************
// compute the stacking of aromatic rings
// latest update: 12/04/2003
// ***************************************************************************
float Ligand::Calculate_AR() const
{
	extern Protein *protein;
	int i,j,id11,id12,id13,id21,id22,id23;
	float d,d1,d2,angle,angle1,angle2,tmp,tmp1,tmp2,sum;

	// aromatic rings on ligand are pre-defined in Molecule::Value_Atom()
        // aromatic rings on protein are pre-defined in Protein::Define_Pocket()

	// note that a ring system like indole is defined as one ring on protein
	// while defined as two rings on ligand
	// thus the following procedure is based primarily on protein rings

	sum=0.000;

	for(i=0;i<protein->num_ring;i++)
	{
	 if(protein->ring[i].valid==0) continue;

	 tmp=0.000;

	 for(j=0;j<this->num_ring;j++)
		{
		 if(this->ring[j].valid==0) continue;

		 d=Distance(protein->ring[i].centroid,this->ring[j].centroid);

		 id11=protein->ring[i].atom_id[0];
		 id12=protein->ring[i].atom_id[1];
		 id13=protein->ring[i].atom_id[2];

		 id21=this->ring[j].atom_id[0];
		 id22=this->ring[j].atom_id[1];
		 id23=this->ring[j].atom_id[2];

		 angle=Angle_of_Two_Planes(protein->atom[id11-1].coor,
                                           protein->atom[id12-1].coor,
                                           protein->atom[id13-1].coor,
                                           this->atom[id21-1].coor,
                                           this->atom[id22-1].coor,
                                           this->atom[id23-1].coor);

		 d1=4.00; d2=5.00;

		 if(d<d1) tmp1=1.000;
                 else if(d<d2) tmp1=(d-d2)/(d1-d2);
                 else tmp1=0.000;

		 angle1=30.0; angle2=60.0;

		 if(angle<angle1) tmp2=1.000;
                 else if(angle<angle2) tmp2=(angle-angle2)/(angle1-angle2);
                 else tmp2=0.000;

		 if((tmp1*tmp2)>tmp) tmp=tmp1*tmp2;
		 else continue;
		}

	 sum+=tmp;
	}

	return sum;
}

// *************************************************************************
// Output detailed score information of each atom
// This function is suitable for outputing a single molecule
// Latest update: 11/21/2003
// *************************************************************************
void Ligand::Write_Out_Log(char *filename) const
{
	extern Input *input;
	FILE *fp;
	int i;

	if((fp=fopen(filename,"w"))==NULL) Open_File_Error(filename);

	fprintf(fp,"### X-Score scoring information\n");
	fprintf(fp,"### Created by X-Score on %s\n\n", Get_Time());

	fprintf(fp,"<MOLECULE> %s\n", this->name);
	fprintf(fp,"<FORMULA> %s\n", this->formula);
	fprintf(fp,"<WEIGHT> %6.1f\n", this->weight);

	fprintf(fp,"<ATOM> %d\n", this->num_atom);
	fprintf(fp,"--------------------------------------------------------------\n");
	fprintf(fp,"ID   Type     VDW     HB     HP     HM      HS     RT   Score\n"); 
	fprintf(fp,"--------------------------------------------------------------\n");

	for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid<=0) continue;

		 fprintf(fp,"%-3d  ", this->atom[i].id);
		 fprintf(fp,"%-5s  ", this->atom[i].type);
		 fprintf(fp,"%6.1f  ", this->abs_inf[i].vdw);
		 fprintf(fp,"%5.1f  ", this->abs_inf[i].hb);
		 fprintf(fp,"%5.1f  ", this->abs_inf[i].hp);
		 fprintf(fp,"%5.2f  ", this->abs_inf[i].hm);
		 fprintf(fp,"%6.1f  ", this->abs_inf[i].hs);
		 fprintf(fp,"%5.1f  ", this->abs_inf[i].rt);
		 fprintf(fp,"%5.2f\n", this->abs_inf[i].score);
		}

	fprintf(fp,"--------------------------------------------------------------\n");
	fprintf(fp,"Total       ");
	fprintf(fp,"%6.1f  ", this->vdw);
	fprintf(fp,"%5.1f  ", this->hb);
	fprintf(fp,"%5.1f  ", this->hp);
	fprintf(fp,"%5.2f  ", this->hm);
	fprintf(fp,"%6.1f  ", this->hs);
	fprintf(fp,"%5.1f  ", this->rt);
	fprintf(fp,"%5.2f\n", this->bind_score);

	fprintf(fp,"--------------------------------------------------------------\n\n");

	fprintf(fp,"Scoring functions used in computation:\n\n");

	fprintf(fp,"HPSCORE = %5.3f + ", input->hpscore_c0);
	fprintf(fp,"(%5.3f)*VDW + ", input->hpscore_cvdw);
	fprintf(fp,"(%5.3f)*HB + ", input->hpscore_chb); 
	fprintf(fp,"(%5.3f)*HP + ", input->hpscore_chp);
	fprintf(fp,"(%6.3f)*RT = ", input->hpscore_crt);
	fprintf(fp,"%5.2f\n", this->pkd1);

	fprintf(fp,"HMSCORE = %5.3f + ", input->hmscore_c0);
        fprintf(fp,"(%5.3f)*VDW + ", input->hmscore_cvdw);
        fprintf(fp,"(%5.3f)*HB + ", input->hmscore_chb);
        fprintf(fp,"(%5.3f)*HM + ", input->hmscore_chm);
        fprintf(fp,"(%6.3f)*RT = ", input->hmscore_crt);
	fprintf(fp,"%5.2f\n", this->pkd2);

	fprintf(fp,"HSSCORE = %5.3f + ", input->hsscore_c0);
        fprintf(fp,"(%5.3f)*VDW + ", input->hsscore_cvdw);
        fprintf(fp,"(%5.3f)*HB + ", input->hsscore_chb);
        fprintf(fp,"(%5.3f)*HS + ", input->hsscore_chs);
        fprintf(fp,"(%6.3f)*RT = ", input->hsscore_crt);
	fprintf(fp,"%5.2f\n", this->pkd3);

	fprintf(fp,"\n<END>\n\n");

	fclose(fp); return;
}

// *************************************************************************
// Output detailed score information of each atom
// This function is suitable for outputing multiple molecules
// Latest update: 11/21/2003
// *************************************************************************
void Ligand::Write_Out_Log(FILE *fp) const
{
	extern Input *input;
	int i;

	fprintf(fp,"### X-Score scoring information\n");
	fprintf(fp,"### Created by X-Score on %s\n\n", Get_Time());

	fprintf(fp,"<MOLECULE> %s\n", this->name);
	fprintf(fp,"<FORMULA> %s\n", this->formula);
	fprintf(fp,"<WEIGHT> %6.1f\n", this->weight);

	fprintf(fp,"<ATOM> %d\n", this->num_atom);
	fprintf(fp,"--------------------------------------------------------------\n");
	fprintf(fp,"ID   Type     VDW     HB     HP     HM      HS     RT   Score\n"); 
	fprintf(fp,"--------------------------------------------------------------\n");

	for(i=0;i<this->num_atom;i++)
		{
		 if(this->atom[i].valid<=0) continue;

		 fprintf(fp,"%-3d  ", this->atom[i].id);
		 fprintf(fp,"%-5s  ", this->atom[i].type);
		 fprintf(fp,"%6.1f  ", this->abs_inf[i].vdw);
		 fprintf(fp,"%5.1f  ", this->abs_inf[i].hb);
		 fprintf(fp,"%5.1f  ", this->abs_inf[i].hp);
		 fprintf(fp,"%5.2f  ", this->abs_inf[i].hm);
		 fprintf(fp,"%6.1f  ", this->abs_inf[i].hs);
		 fprintf(fp,"%5.1f  ", this->abs_inf[i].rt);
		 fprintf(fp,"%5.2f\n", this->abs_inf[i].score);
		}

	fprintf(fp,"--------------------------------------------------------------\n");
	fprintf(fp,"Total       ");
	fprintf(fp,"%6.1f  ", this->vdw);
	fprintf(fp,"%5.1f  ", this->hb);
	fprintf(fp,"%5.1f  ", this->hp);
	fprintf(fp,"%5.2f  ", this->hm);
	fprintf(fp,"%6.1f  ", this->hs);
	fprintf(fp,"%5.1f  ", this->rt);
	fprintf(fp,"%5.2f\n", this->bind_score);

	fprintf(fp,"--------------------------------------------------------------\n\n");

	fprintf(fp,"Scoring functions used in computation:\n\n");

	fprintf(fp,"HPSCORE = %5.3f + ", input->hpscore_c0);
	fprintf(fp,"(%5.3f)*VDW + ", input->hpscore_cvdw);
	fprintf(fp,"(%5.3f)*HB + ", input->hpscore_chb); 
	fprintf(fp,"(%5.3f)*HP + ", input->hpscore_chp);
	fprintf(fp,"(%6.3f)*RT = ", input->hpscore_crt);
	fprintf(fp,"%5.2f\n", this->pkd1);

	fprintf(fp,"HMSCORE = %5.3f + ", input->hmscore_c0);
        fprintf(fp,"(%5.3f)*VDW + ", input->hmscore_cvdw);
        fprintf(fp,"(%5.3f)*HB + ", input->hmscore_chb);
        fprintf(fp,"(%5.3f)*HM + ", input->hmscore_chm);
        fprintf(fp,"(%6.3f)*RT = ", input->hmscore_crt);
	fprintf(fp,"%5.2f\n", this->pkd2);

	fprintf(fp,"HSSCORE = %5.3f + ", input->hsscore_c0);
        fprintf(fp,"(%5.3f)*VDW + ", input->hsscore_cvdw);
        fprintf(fp,"(%5.3f)*HB + ", input->hsscore_chb);
        fprintf(fp,"(%5.3f)*HS + ", input->hsscore_chs);
        fprintf(fp,"(%6.3f)*RT = ", input->hsscore_crt);
	fprintf(fp,"%5.2f\n", this->pkd3);

	fprintf(fp,"\n<END>\n\n");

	fflush(fp); return;
}

