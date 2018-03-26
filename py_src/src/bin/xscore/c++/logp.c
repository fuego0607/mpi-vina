# include "xtool.h"

const float LOGP_HYDROPHOBIC_CARBON = 0.211;
const float LOGP_INTERNAL_HBOND = 0.429;
const float LOGP_HALOGEN_PAIR = 0.137;
const float LOGP_NAR_PAIR = 0.485;
const float LOGP_O3_PAIR = -0.268;
const float LOGP_ACCEPTOR_PAIR = 0.580;
const float LOGP_AMINO_ACID = -2.166;
const float LOGP_SALICYLIC_ACID = 0.554;
const float LOGP_SULFONIC_ACID = -0.501;

struct LogP_Factor
	{
	 char symbol[30];
	 float num;
	 float coeff;
	} logp_factor[10];

// *****************************************************************************
// calculate the chemical properties of the input molecule(s)
// also dissect the LogP value into each atom
// Latest update: 04/11/2003
// *****************************************************************************
int Xtool_LogP()
{
        extern Input *input;
        extern ForceField *ff;
	extern LogP_Factor logp_factor[10];
	FILE *fin, *fout;
	int i,count,num;
	Molecule mol;

	printf("Now reading parameters from '%s' ...\n",input->parameter_dir);
        ff=new ForceField(input->parameter_dir);
        if(ff==NULL) Memory_Allocation_Error();

	printf("Now reading the ligand from '%s' ...\n", input->input_file);
	num=Check_Mol2_File(input->input_file);

	if((fin=fopen(input->input_file,"r"))==NULL) 
		{
		 Open_File_Error(input->input_file);
		}

	if(strcmp(input->output_file,"none"))
		{
		 if((fout=fopen(input->output_file,"w"))==NULL) 
			{
			 Open_File_Error(input->output_file);
			}
		}
	else
		{
		 fout=NULL;
		}

	count=0;

	for(i=0;i<num;i++)
	{
	 mol.Clear();
	 mol.Read_From_Mol2(fin);
	 mol.Value_Atom();

	 if(mol.num_subst>1||mol.num_subst==0)
		{
		 printf("Invalid molecule: %s ... skipped\n", mol.name);
		 continue;
		}
	 else
		{
		 printf("Molecule: ");
		}

	 if(!strcmp(((LOGP_Input *)input)->calculate_mw,"YES"))
		{
	 	 // mol.weight=mol.Get_Weight();
	 	 printf("MW= %6.1f ", mol.weight);
		}

	 if(!strcmp(((LOGP_Input *)input)->calculate_logp,"YES"))
		{
		 mol.logp=mol.Calculate_LogP();
	 	 printf("logP= %6.2f ", mol.logp);
		}

	 if(!strcmp(((LOGP_Input *)input)->count_hb_atom,"YES"))
		{
		 mol.num_hb_atom=mol.Get_Num_HB_Atom();
	 	 printf("HB#= %2d ", mol.num_hb_atom);
		}

	 if(!strcmp(((LOGP_Input *)input)->count_rotor,"YES"))
                {
		 mol.num_rotor=(int)mol.Count_Rotor();
	 	 printf("rotor= %2d ", mol.num_rotor);
                }

	 printf("formula= %-s ", mol.formula);
	 printf("name= %-s\n", mol.name);

	 if(fout!=NULL) mol.Write_Out_Log(fout);

	 count++;
	}	

	printf("%d molecules are processed; %d are valid.\n", num, count);

	fclose(fin); fclose(fout);

	return TRUE;
}

// *****************************************************************************
// The shortcut form of Xtool_LogP()
// Latest update: 11/17/2003
// *****************************************************************************
int Xtool_LogP_Shortcut()
{
 extern Input *input;
 extern ForceField *ff;
 extern LogP_Factor logp_factor[10];
 FILE *fin, *fout;
 int i,count,num;
 Molecule mol;

 printf("Now reading parameters from '%s' ...\n",input->parameter_dir);
 ff=new ForceField(input->parameter_dir);
 if(ff==NULL) Memory_Allocation_Error();

 printf("Now reading the ligand from '%s' ...\n", input->input_file);
 num=Check_Mol2_File(input->input_file);

 if((fin=fopen(input->input_file,"r"))==NULL) 
	{
	 Open_File_Error(input->input_file);
	}

 if(strcmp(input->output_file,"none"))
	{
	 if((fout=fopen(input->output_file,"w"))==NULL) 
		{
		 Open_File_Error(input->output_file);
		}
	}
 else
	{
	 fout=NULL;
	}

 count=0;

 puts("");
 puts("**********************************************************************");

 for(i=0;i<num;i++)
	{
 	 mol.Clear(); mol.Read_From_Mol2(fin);
	 mol.Value_Atom();

	 if(mol.num_subst>1||mol.num_subst==0)
		{
		 printf("Invalid molecule: %s ... skipped\n", mol.name);
		 continue;
		}

	 printf("Molecule: ");

	 // mol.weight=mol.Get_Weight(); // already calculated in Value_Atom()
	 printf("MW = %-6.1f ", mol.weight);

	 mol.logp=mol.Calculate_LogP();
	 printf("logP = %-6.2f ", mol.logp);

	 // mol.Get_Formula(mol.formula); // already got in Value_Atom()
	 printf("formula = %-s ", mol.formula);

	 printf("name = %-s\n", mol.name);

	 if(fout!=NULL) mol.Write_Out_Log(fout);

	 count++;
	}	

 puts("**********************************************************************");

 printf("%d molecules are processed; %d are valid.\n", num, count);

 fclose(fin); fclose(fout);

 return TRUE;
}

// *****************************************************************************
// calculate the LogP value for the molecule
// also dissect the LogP value into each atom
// Latest update: 11/13/2002
// *****************************************************************************
float Molecule::Calculate_LogP()
{
	extern ForceField *ff;
        int i,flag;
        float xlogp;
	char type[20];
	extern LogP_Factor logp_factor[10];

	strcpy(logp_factor[0].symbol,"Hydrophobic carbon");
        strcpy(logp_factor[1].symbol,"Internal H-bond");
        strcpy(logp_factor[2].symbol,"Halogen 1-3 pair");
        strcpy(logp_factor[3].symbol,"Aromatic nitrogen 1-4 pair");
        strcpy(logp_factor[4].symbol,"Ortho sp3 oxygen pair");
        strcpy(logp_factor[5].symbol,"Acceptor 1-5 pair");
        strcpy(logp_factor[6].symbol,"Paralleled donor pair");
        strcpy(logp_factor[7].symbol,"Alpha amino acid");
        strcpy(logp_factor[8].symbol,"Salicylic acid");
        strcpy(logp_factor[9].symbol,"P-amino sulfonic acid");

	logp_factor[0].coeff = 0.211;
	logp_factor[1].coeff = 0.429;
        logp_factor[2].coeff = 0.137;
        logp_factor[3].coeff = 0.485;
        logp_factor[4].coeff =-0.268;
        logp_factor[5].coeff = 0.580;
        logp_factor[6].coeff =-0.423;
        logp_factor[7].coeff =-2.166;
        logp_factor[8].coeff = 0.554;
        logp_factor[9].coeff =-0.501;

	xlogp=0.000;

        for(i=0;i<num_atom;i++)
                {
                 if(atom[i].valid==0) continue;

		 Get_XLOGP_Type(atom[i].id,type);
		 atom[i].logp=ff->Get_Atom_LogP(type);

                 xlogp+=atom[i].logp;
                }

	for(i=0;i<10;i++) logp_factor[i].num=0;

	// this flag determines if correction factors are also distributed
	// onto each atoms or not

	flag=FALSE;

        logp_factor[0].num=Count_Hydrophobic_Carbon(flag);
        logp_factor[1].num=Count_Internal_HBond(flag);
        logp_factor[2].num=Count_Halogen_1_3_Pair(flag);
        logp_factor[3].num=Count_Nar_1_4_Pair(flag);
        logp_factor[4].num=Count_O3_1_4_Pair(flag);
        logp_factor[5].num=Count_Acceptor_1_5_Pair(flag);
	logp_factor[6].num=0.000;
        logp_factor[7].num=Count_Amino_Acid(flag);
        logp_factor[8].num=Count_Salicylic_Acid(flag);
	logp_factor[9].num=Count_Sulfonic_Acid(flag);

	for(i=0;i<10;i++) xlogp+=(logp_factor[i].num*logp_factor[i].coeff);

        return xlogp;
}

int Molecule::Get_XLOGP_Type(int atom_id, char *type) const
{
	Group group;
	Atom atom;

	atom=this->atom[atom_id-1]; group=Find_X_Group(atom.id);

	strcpy(type,"Un"); 

	if(!strcmp(atom.type,"H")||!strcmp(atom.type,"H.spc"))
	{
	 if(group.neib[0].type[0]=='O') strcpy(type,"H.hb");
	 else if(group.neib[0].type[0]=='N') strcpy(type,"H.hb");
	 // else if(!strcmp(group.neib[0].type,"S.3")) strcpy(type,"H.hb");
	 else strcpy(type,"H"); 
	}

	 if(!strcmp(atom.type,"C.3"))
		{
		 if(group.num_nonh==1)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0) strcpy(type,"C.3.h3.pi=0");
				 else strcpy(type,"C.3.h3.pi=1");
				}
			 else
				{
                                 strcpy(type,"C.3.h3.x");
				}
			}
		 else if(group.num_nonh==2)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0)
					{
					 strcpy(type,"C.3.h2.pi=0");
					}
				 else if(group.num_pi==1) 
					{ 
					 strcpy(type,"C.3.h2.pi=1"); 
					}
				 else
					{
					 strcpy(type,"C.3.h2.pi=2");
					}
				}
			 else 
				{
				 if(group.num_pi==0)
                                        {
                                         strcpy(type,"C.3.h2.x.pi=0");
                                        }
                                 else if(group.num_pi==1)
                                        {
                                         strcpy(type,"C.3.h2.x.pi=1");
                                        }
                                 else
                                        {
                                         strcpy(type,"C.3.h2.x.pi=2");
                                        }
				}
			}
		 else if(group.num_nonh==3)
			{
			 if(group.num_hetero==0)
				{
                                 if(group.num_pi==0)
                                        {
                                         strcpy(type,"C.3.h.pi=0");
                                        }
                                 else if(group.num_pi==1)
                                        {
                                         strcpy(type,"C.3.h.pi=1");
                                        }
                                 else
                                        {
                                         strcpy(type,"C.3.h.pi>1");
                                        }
				}
			 else 
				{
                                 if(group.num_pi==0)
                                        {
                                         strcpy(type,"C.3.h.x.pi=0");
                                        }
                                 else if(group.num_pi==1)
                                        {
                                         strcpy(type,"C.3.h.x.pi=1");
                                        }
                                 else
                                        {
                                         strcpy(type,"C.3.h.x.pi>1");
                                        }
				}
			}
		 else if(group.num_nonh==4)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0)
                                       	{
                                         strcpy(type,"C.3.pi=0");
                                       	}
                                 else if(group.num_pi==1) 
					{ 
					 strcpy(type,"C.3.pi=1");
                                       	}
                                 else
                                       	{
                                         strcpy(type,"C.3.pi>1");
                                       	}
				}
			 else 
				{
				 if(group.num_pi==0) strcpy(type,"C.3.x.pi=0");
                                 else strcpy(type,"C.3.x.pi>0");
				}
			}
		 else	
			{
			 strcpy(type,"C.3.unknown");
			}
		}

	if(!strcmp(atom.type,"C.2"))
		{
		 if(group.num_nonh==1)
			{
			 strcpy(type,"C.2.h2");	
			}
		 else if(group.num_nonh==2)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0) strcpy(type,"C.2.h.pi=0");
				 else strcpy(type,"C.2.h.pi=1");
				}
			 else
				{
				 if(group.num_pi==0) 
					{
					 strcpy(type,"C.2.h.x.pi=0");
					}
                                 else strcpy(type,"C.2.h.x.pi=1");
				}
			}
		 else if(group.num_nonh==3)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0) strcpy(type,"C.2.pi=0"); 
				 else strcpy(type,"C.2.pi>0");
				}
			 else if(group.num_hetero==1)
				{
				 if(group.num_pi==0) strcpy(type,"C.2.x.pi=0"); 
                                 else strcpy(type,"C.2.x.pi>0");
				}
			 else
				{
				 if(group.num_pi==0) strcpy(type,"C.2.x2.pi=0");
                                 else strcpy(type,"C.2.x2.pi>0");
				}
			}
		 else
			{
			 strcpy(type,"C.2.unknown");
			}
		}

	 if(!strcmp(atom.type,"C.cat"))
		{
		 strcpy(type,"C.2.x2.pi>0");
		}

	 if(!strcmp(atom.type,"C.ar"))
		{
		 if(group.num_nonh==2)
			{
			 if(group.num_nar==0) strcpy(type,"C.ar.h");
			 else strcpy(type,"C.ar.h.(X)");
			}
		 else if(group.num_nonh==3)
			{
			 if(group.num_nar==0)
				{
				 if(group.num_hetero==0) strcpy(type,"C.ar");
				 else strcpy(type,"C.ar.x");
				}
			 else 
				{
				 if(group.num_hetero==0)
                                        {
                                         strcpy(type,"C.ar.(X)");
                                        }
                                 else 
                                        {
                                         strcpy(type,"C.ar.(X).x");
                                        }
				}
			}
		 else	
			{
			 strcpy(type,"C.ar.unknown");
			}
		}

	if(!strcmp(atom.type,"C.1"))
		{
		 if(group.db_type!=0) strcpy(type,"C.1.==");
		 else if(group.num_nonh==1) strcpy(type,"C.1.h");
		 else if(group.num_nonh==2) strcpy(type,"C.1");
		 else strcpy(type,"C.1.unknown");
		}

	if(!strcmp(atom.type,"N.4")||
	   !strcmp(atom.type,"N.3")||
	   !strcmp(atom.type,"N.pl3"))
		{
		 if(group.num_nonh==1)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0) 
					{
					 strcpy(type,"N.3.h2.pi=0");
					}
				 else 
					{
					 strcpy(type,"N.3.h2.pi=1");
					}
				}
			 else 
				{
				 strcpy(type,"N.3.h2.x");
				}
			}
		 else if(group.num_nonh==2)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0)
					{
					 strcpy(type,"N.3.h.pi=0");
					}
				 else if(group.num_pi==1)
					{
					 strcpy(type,"N.3.h.pi>0");
					}
				 else if(atom.ring==0)
					{
					 strcpy(type,"N.3.h.pi>0");
					}
				 else strcpy(type,"N.3.h.ring");
				}
                         else 
                                {
				 if(group.num_pi<=1)
                                        {
                                         strcpy(type,"N.3.h.x");
                                        }
                                 else if(atom.ring==0)
                                        {
                                         strcpy(type,"N.3.h.x");
                                        }
                                 else
                                        {
                                         strcpy(type,"N.3.h.x.ring");
                                        }
                                }
			}
		 else if(group.num_nonh==3)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0) 
					{
					 strcpy(type,"N.3.pi=0");
					}
				 else if(group.num_pi==1) 
					{
					 strcpy(type,"N.3.pi>0");
					}
				 else if(atom.ring==0)
					{
					 strcpy(type,"N.3.pi>0");
					}
				 else
					{
					 strcpy(type,"N.3.ring");
					}
				}
                         else 
                                {
                                 if(group.num_pi<=1)
                                        {
                                         strcpy(type,"N.3.x");
                                        }
                                 else if(atom.ring==0)
                                        {
                                         strcpy(type,"N.3.x");
                                        }
                                 else
                                        {
                                         strcpy(type,"N.3.x.ring");
                                        }
                                }
			}
		 else
			{
			 strcpy(type,"N.3.unknown");
			}
		}

	if(!strcmp(atom.type,"N.am"))
		{
		 if(group.num_nonh==1) strcpy(type,"N.am.h2");
		 else if(group.num_nonh==2)
			{
			 if(group.num_hetero==0) strcpy(type,"N.am.h");
			 else strcpy(type,"N.am.h.x");
			}
		 else if(group.num_nonh==3)
			{
			 if(group.num_hetero==0) strcpy(type,"N.am");
			 else strcpy(type,"N.am.x");
			}
		 else strcpy(type,"N.am.unknown");
		}

	if(!strcmp(atom.type,"N.2"))
		{
		 // N=C, N=S, N=P
		 if(group.db_type==1||group.db_type==4||group.db_type==5)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0) 
					{ 
					 strcpy(type,"N.2.(=C).pi=0");
					}
				 else
					{
					 strcpy(type,"N.2.(=C).pi=1");
					}
				}
			 else
				{
				 if(group.num_pi==0)
                                      {
                                       strcpy(type,"N.2.(=C).x.pi=0");
                                      }
                                 else
                                      {
                                       strcpy(type,"N.2.(=C).x.pi=1");
                                      }
				}
			}
		 else if(group.db_type==2)
			{
			 if(group.num_hetero==0) strcpy(type,"N.2.(=N)");
			 else strcpy(type,"N.2.(=N).x");
			}
		 else if(group.db_type==3)
			{
			 if(group.num_nonh==2) strcpy(type,"N.2.o");
			 else if(group.num_nonh==3) strcpy(type,"N.2.o2");
			 else strcpy(type,"N.2.o");
			}
		 else
			{
			 strcpy(type,"N.2.unknown");
			}
		}

	if(!strcmp(atom.type,"N.ar")) strcpy(type,"N.ar");

	if(!strcmp(atom.type,"N.1")) strcpy(type,"N.1");

	if(!strcmp(atom.type,"O.3"))
		{
		 if(group.num_nonh==1)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0) strcpy(type,"O.3.h.pi=0");
				 else strcpy(type,"O.3.h.pi=1");
				}
			 else 
				{
                                 strcpy(type,"O.3.h.x");
				}
			}
		 else if(group.num_nonh==2)
			{
			 if(group.num_hetero==0)
				{
				 if(group.num_pi==0) strcpy(type,"O.3.pi=0");
				 else strcpy(type,"O.3.pi>0");
				}
			 else
				{
				 strcpy(type,"O.3.x");
				}
			}
		 else
			{
			 strcpy(type,"O.3.unknown");
			}
		}

	if(!strcmp(atom.type,"O.2")) strcpy(type,"O.2");

	if(!strcmp(atom.type,"O.co2")) strcpy(type,"O.co2");

	if(!strcmp(atom.type,"S.3"))
		{
		 if(group.num_nonh==1) strcpy(type,"S.3.h");
		 else if(group.num_nonh==2) strcpy(type,"S.3");
		 else strcpy(type,"S.3.unknown");
		}

	if(!strcmp(atom.type,"S.2")) strcpy(type,"S.2");

	if(!strcmp(atom.type,"S.o")) strcpy(type,"S.o");

	if(!strcmp(atom.type,"S.o2")) strcpy(type,"S.o2");

	if(!strcmp(atom.type,"P.3"))
		{
		 if(group.db_type==3) strcpy(type,"P.3.(=O)");
		 else if(group.db_type==4) strcpy(type,"P.3.(=S)");
		 else strcpy(type,"P.3.unknown");
		}

	if(!strcmp(atom.type,"F"))
		{
                 if(group.num_pi==0) strcpy(type,"F.pi=0");
                 else if(group.num_pi==1) strcpy(type,"F.pi=1");
		 else strcpy(type,"F.unknown");
		}

	if(!strcmp(atom.type,"Cl"))
		{
                 if(group.num_pi==0) strcpy(type,"Cl.pi=0");
                 else if(group.num_pi==1) strcpy(type,"Cl.pi=1");
		 else strcpy(type,"Cl.unknown");
		}

	if(!strcmp(atom.type,"Br"))
		{
                 if(group.num_pi==0) strcpy(type,"Br.pi=0");
                 else if(group.num_pi==1) strcpy(type,"Br.pi=1");
		 else strcpy(type,"Br.unknown");
		}

	if(!strcmp(atom.type,"I"))
		{
		 if(group.num_pi==0) strcpy(type,"I.pi=0");
		 else if(group.num_pi==1) strcpy(type,"I.pi=1");
		 else strcpy(type,"I.unknown");
		}

	if(!strcmp(atom.type,"Si")) strcpy(type,"Si");

	if(!strcmp(type,"Un")) return FALSE;
	return TRUE;
}

Group Molecule::Find_X_Group(int atom_id) const
{
        int i,num;
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

        group.num_neib=num;

        group.num_h=0; group.num_nonh=0;

        for(i=0;i<num;i++)
                {
                 if(group.neib[i].type[0]=='H') group.num_h++;
                 else group.num_nonh++;
                }

        group.num_hetero=0;

        for(i=0;i<group.num_nonh;i++)
                {
                 if(!strcmp(group.bond[i].type,"2")) continue;
                 else if(!strcmp(group.bond[i].type,"3")) continue;
                 else if(!strcmp(group.bond[i].type,"ar")) continue;
                 else if(group.neib[i].type[0]=='N') group.num_hetero++;
                 else if(group.neib[i].type[0]=='O') group.num_hetero++;
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

        group.num_nar=group.num_car=0;

        for(i=0;i<group.num_nonh;i++)
                {
                 if(strcmp(group.bond[i].type,"ar")) continue;
                 else if(!strcmp(group.neib[i].type,"N.ar")) group.num_nar++;
                 else if(!strcmp(group.neib[i].type,"C.ar")) group.num_car++;
                 else if(group.neib[i].type[0]=='C') group.num_car++;
                 else group.num_nar++;
                }

        group.db_type=0;

        for(i=0;i<group.num_nonh;i++)
                {
                 if(!strcmp(group.center.type,"O.co2")||
                    !strcmp(group.center.type,"O.2")||
                    !strcmp(group.center.type,"S.2"))
                        {
                         if(group.neib[i].type[0]=='C') group.db_type=1;
                         else if(group.neib[i].type[0]=='N') group.db_type=2;
                         else if(group.neib[i].type[0]=='O') group.db_type=3;
                         else if(group.neib[i].type[0]=='S') group.db_type=4;
                         else if(group.neib[i].type[0]=='P') group.db_type=5;
                         else continue;
                        }
                 else if(!strcmp(group.bond[i].type,"2")||
                         !strcmp(group.neib[i].type,"O.co2")||
                         !strcmp(group.neib[i].type,"O.2")||
                         !strcmp(group.neib[i].type,"S.2"))
                        {
                         if(group.neib[i].type[0]=='C') group.db_type=1;
                         else if(group.neib[i].type[0]=='N') group.db_type=2;
                         else if(group.neib[i].type[0]=='O') group.db_type=3;
                         else if(group.neib[i].type[0]=='S') group.db_type=4;
                         else if(group.neib[i].type[0]=='P') group.db_type=5;
                         else continue;
                        }
                 else continue;
                }

        group.valid=1; return group;
}

float Molecule::Count_Hydrophobic_Carbon(int flag)
{
        int i,num;

        num=0;

        for(i=0;i<num_atom;i++)
                {
                 if(strcmp(atom[i].type,"C.3")&&
                    strcmp(atom[i].type,"C.2")) continue;
                 else if(Hydrophobic_Neighbor_Check(atom[i].id)
                         ==FALSE) continue;

                 num++; 

		 if(flag) atom[i].logp+=(LOGP_HYDROPHOBIC_CARBON);
		 else continue;
                }

        if(num>=10) num/=2;

        return (float)num;
}

int Molecule::Hydrophobic_Neighbor_Check(int id) const
{
        int i,mark;

        id--;

        mark=TRUE;

        for(i=0;i<num_atom;i++)
                {
                 if(i==id) continue;
                 else if(!strcmp(atom[i].type,"F")||
                         !strcmp(atom[i].type,"Cl")||
                         !strcmp(atom[i].type,"Br")||
                         !strcmp(atom[i].type,"I")||
			 !strcmp(atom[i].type,"Si")||
                         atom[i].type[0]=='N'||
                         atom[i].type[0]=='O'||
                         atom[i].type[0]=='S'||
                         atom[i].type[0]=='P')
                        {
			 if(Connection_1_2_Check(atom[id].id,atom[i].id)==TRUE)
				{mark=FALSE; break;}
			 else if(Connection_1_3_Check(atom[id].id,atom[i].id)==TRUE)
				{mark=FALSE; break;}
			 else if(Connection_1_4_Check(atom[id].id,atom[i].id)==TRUE)
				{mark=FALSE; break;}
                         else continue;
                        }
                 else continue;
                }

        return mark;
}

float Molecule::Count_Internal_HBond(int flag)
{
	int i,j,mark1,mark2;
	float num;
	int *record;

	record=new int[num_atom];
	if(record==NULL) Memory_Allocation_Error();

	for(i=0;i<num_atom;i++) record[i]=0;

	num=0;

	for(i=0;i<num_atom;i++)
		{
		 if(strcmp(atom[i].hb,"D")&&strcmp(atom[i].hb,"DA")) continue;
		 else if(atom[i].ring!=0) continue; // not allowed in ring

		 if(Adjacent_Ring_Check(atom[i].id)==FALSE) mark1=FALSE;
		 else mark1=TRUE;

		 for(j=0;j<num_atom;j++)
			{
			 if(i==j) continue;
			 else if(strcmp(atom[j].hb,"A")&&strcmp(atom[j].hb,"DA")) continue;
			 else if(!strcmp(atom[j].type,"O.3")&&!strcmp(atom[j].hb,"A")) continue;
			 else if(!strcmp(atom[j].type,"N.2")) continue;
			 else if(!strcmp(atom[j].type,"N.ar")) continue;
			 else if(atom[j].ring!=0) continue; // not in ring 

			 if(Adjacent_Ring_Check(atom[j].id)==FALSE) mark2=FALSE;
			 else mark2=TRUE;

			 if(mark1==TRUE&&mark2==TRUE)
				{
				 if(Connection_1_4_Check(atom[i].id,
				    atom[j].id)==FALSE) continue;
				 else 
					{
					 record[i]++; record[j]++; 
					}
				}
			 else if(mark1==TRUE&&mark2==FALSE)
			 	{
				 if(Connection_1_5_Check(atom[i].id,
                                    atom[j].id)==FALSE) continue;
				 else
                                        {
                                         record[i]++; record[j]++;
                                        }
				}
			 else if(mark2==FALSE&&mark2==TRUE)
				{
				 if(Connection_1_5_Check(atom[i].id,
                                    atom[j].id)==FALSE) continue;
				 else
                                        {
                                         record[i]++; record[j]++;
                                        }
				}
			 else continue;
			}
		}

	for(i=0;i<num_atom;i++)
		{
		 if(record[i]==0) continue;

		 num+=0.500;

		 if(flag) atom[i].logp+=(LOGP_INTERNAL_HBOND/2.0);
		 else continue;
		}

	if(record) delete [] record;

	return num;
}

int Molecule::Adjacent_Ring_Check(int id) const
{
        int i,num,tmp,mark;

        id--;

        mark=FALSE; num=atom[id].num_neib;

        for(i=0;i<num;i++)
                {
                 tmp=atom[id].neib[i];
                 if(atom[tmp-1].ring!=0) {mark=tmp; break;}
                 else continue;
                }

        return mark;
}

float Molecule::Count_Halogen_1_3_Pair(int flag)
{
        int i,j;
        int num1,num2;

        num1=num2=0;

        for(i=0;i<num_atom-1;i++)
                {
                 if(strcmp(atom[i].type,"F")) continue;

                 for(j=i+1;j<num_atom;j++)
                        {
                         if(strcmp(atom[j].type,"F")) continue;
                         else if(Connection_1_3_Check(atom[i].id,atom[j].id)
                                 ==FALSE) continue;

                         num1++;

			 if(flag)
				{
                                 atom[i].logp+=((LOGP_HALOGEN_PAIR)/2.0);
                                 atom[j].logp+=((LOGP_HALOGEN_PAIR)/2.0);
				}
			 else continue;
                        }
                }

        for(i=0;i<num_atom-1;i++)
                {
                 if(strcmp(atom[i].type,"Cl")&&
                    strcmp(atom[i].type,"Br")&&
                    strcmp(atom[i].type,"I")) continue;

                 for(j=i+1;j<num_atom;j++)
                        {
                         if(strcmp(atom[j].type,"Cl")&&
                            strcmp(atom[j].type,"Br")&&
                            strcmp(atom[j].type,"I")) continue;
                         else if(Connection_1_3_Check(atom[i].id,atom[j].id)
				==FALSE) continue;

                         num2++;

			 if(flag)
				{
                                 atom[i].logp+=((LOGP_HALOGEN_PAIR)/2.0);
                                 atom[j].logp+=((LOGP_HALOGEN_PAIR)/2.0);
                                }
			 else continue;
                        }
                }

	return (float)(num1+num2);
}

float Molecule::Count_Nar_1_4_Pair(int flag)
{
        int i,j,num,tmp1,tmp2,tmp3,tmp4;

        num=0;

        for(i=0;i<num_atom-1;i++)
                {
                 if(strcmp(atom[i].type,"N.ar")) continue;

                 tmp1=atom[i].neib[0]; tmp2=atom[i].neib[1];

                 for(j=i+1;j<num_atom;j++)
                        {
                         if(strcmp(atom[j].type,"N.ar")) continue;
                         else if(Connection_1_4_Check(atom[i].id, atom[j].id)
                                 ==FALSE) continue;
                         else
                                {
                                 tmp3=atom[j].neib[0]; tmp4=atom[j].neib[1];

                                 if(Connection_1_2_Check(tmp1,tmp3)==TRUE)
                                        {
                                         if(Connection_1_2_Check(tmp2,tmp4)==TRUE)
                                                {
                                                 num++;
						 if(flag)
						{
                                                 atom[i].logp+=((LOGP_NAR_PAIR)/2.0);
                                                 atom[j].logp+=((LOGP_NAR_PAIR)/2.0);
						}
                                                }
                                         else continue;
                                        }
                                 else if(Connection_1_2_Check(tmp1,tmp4)==TRUE)
                                        {
                                         if(Connection_1_2_Check(tmp2,tmp3)
                                            ==TRUE)
                                                {
                                                 num++;
						 if(flag)
						{
                                                 atom[i].logp+=((LOGP_NAR_PAIR)/2.0);
                                                 atom[j].logp+=((LOGP_NAR_PAIR)/2.0);
						}
                                                }
                                         else continue;
                                        }
                                 else continue;
                                }
                        }
                }

        return (float)num;
}

float Molecule::Count_O3_1_4_Pair(int flag)
{
	int i,j;
	float num;
	int *record;

	record=new int[num_atom];
	if(record==NULL) Memory_Allocation_Error();

	for(i=0;i<num_atom;i++) record[i]=0;

	num=0;

	for(i=0;i<num_atom-1;i++)
		{
		 if(strcmp(atom[i].type,"O.3")) continue;
		 else if(!strcmp(atom[i].hb,"DA")) continue;
		 else if(atom[i].ring!=0) continue;
		 else if(Adjacent_Aromatic_Check(atom[i].id)==FALSE) continue;

		 for(j=i+1;j<num_atom;j++)
			{
			 if(strcmp(atom[j].type,"O.3")) continue;
			 else if(!strcmp(atom[j].hb,"DA")) continue;
			 else if(atom[j].ring!=0) continue;
                 	 else if(Adjacent_Aromatic_Check(atom[j].id)==FALSE) 
				continue;
			 else if(Connection_1_4_Check(atom[i].id,atom[j].id)
				 ==FALSE) continue;
			 else 
				{
				 record[i]++; record[j]++; 
				}
			}
		}

	for(i=0;i<num_atom;i++)
		{
		 if(record[i]==0) continue;

		 num+=0.500;

		 if(flag) atom[i].logp+=(LOGP_O3_PAIR/2.0);
		 else continue;
		}

	if(record) delete [] record;

	return (float)num;
}

float Molecule::Count_Acceptor_1_5_Pair(int flag)
{
        int i,j,num,tmp1,tmp2;

        num=0;

        for(i=0;i<num_atom-1;i++)
                {
                 if(strcmp(atom[i].hb,"A")) continue;
		 else if(!strcmp(atom[i].type,"O.3")) continue;
		 else if(!strcmp(atom[i].type,"N.2")) continue;
		 else if(!strcmp(atom[i].type,"N.ar")) continue;

                 tmp1=atom[i].neib[0]-1;
                 if(atom[tmp1].type[0]=='S') continue;
                 else if(atom[tmp1].type[0]=='P') continue;
                 else if(atom[tmp1].ring!=0) continue;

                 for(j=i+1;j<num_atom;j++)
                        {
                         if(strcmp(atom[j].hb,"A")) continue;
			 else if(!strcmp(atom[j].type,"O.3")) continue;
			 else if(!strcmp(atom[j].type,"N.2")) continue;
			 else if(!strcmp(atom[j].type,"N.ar")) continue;

                         tmp2=atom[j].neib[0]-1;
                         if(atom[tmp2].type[0]=='S') continue;
                         else if(atom[tmp2].type[0]=='P') continue;
                         else if(atom[tmp2].ring!=0) continue;

                         if(Connection_1_5_Check(atom[i].id, atom[j].id)
                                 ==FALSE) continue;

                         num++;

			 if(flag)
				{
                                 atom[i].logp+=((LOGP_ACCEPTOR_PAIR)/2.0);
                                 atom[j].logp+=((LOGP_ACCEPTOR_PAIR)/2.0);
				}
			 else continue;
                        }
                }

        return (float)num;
}

float Molecule::Count_Salicylic_Acid(int flag)
{
        int i,j,num,tmp,mark;

        num=0;

        for(i=0;i<num_atom;i++)
                {
                 if(strcmp(atom[i].type,"O.2")) continue;

                 tmp=atom[i].neib[0]-1;
                 if(atom[tmp].type[0]!='C') continue;
                 else if(atom[tmp].ring!=0) continue;
                 else if(Adjacent_Aromatic_Check(atom[tmp].id)==FALSE) continue;
                 mark=FALSE;

                 for(j=0;j<num_atom;j++)
                        {
                         if(i==j) continue;
                         else if(strcmp(atom[j].type,"O.3")) continue;
                         else if(atom[j].ring!=0) continue;
                         else if(Connection_1_3_Check(atom[i].id, atom[j].id)
                                 ==FALSE) continue;
                         else {mark=TRUE; break;}
                        }

                 if(mark==FALSE) continue;

                 mark=FALSE;

                 for(j=0;j<num_atom;j++)
                        {
                         if(i==j) continue;
                         else if(strcmp(atom[j].type,"O.3")) continue;
                         else if(strcmp(atom[j].hb,"DA")) continue;
                         else if(Adjacent_Aromatic_Check(atom[j].id)==FALSE)
                                continue;
                         else if(Connection_1_5_Check(atom[i].id, atom[j].id)
                                 ==FALSE) continue;

                         num++;

			 if(flag)
				{
                                 atom[i].logp+=((LOGP_SALICYLIC_ACID)/2.0);
                                 atom[j].logp+=((LOGP_SALICYLIC_ACID)/2.0);
				}

                         break;
                        }

                 if(num!=0) break;
                }

        return (float)num;
}

int Molecule::Adjacent_Aromatic_Check(int id) const
{
        int i,num,tmp,mark;

        id--;

        mark=FALSE; num=atom[id].num_neib;

        for(i=0;i<num;i++)
                {
                 tmp=atom[id].neib[i];
                 if(!strcmp(atom[tmp-1].type,"C.ar")) {mark=tmp; break;}
                 else continue;
                }

        return mark;
}

float Molecule::Count_Amino_Acid(int flag)
{
        int i,j,tmp,num,mark;

        num=0;

        for(i=0;i<num_atom;i++)
                {
                 if(strcmp(atom[i].type,"O.2")) continue;

                 tmp=atom[i].neib[0]-1;
                 if(atom[tmp].type[0]!='C') continue;
                 else if(atom[tmp].ring!=0) continue;

                 mark=FALSE;

                 for(j=0;j<num_atom;j++)
                        {
                         if(i==j) continue;
                         else if(strcmp(atom[j].type,"O.3")) continue;
                         else if(strcmp(atom[j].hb,"DA")) continue;
                         else if(atom[j].ring!=0) continue;
                         else if(Connection_1_3_Check(atom[i].id, atom[j].id)
                                 ==FALSE) continue;
                         else {mark=TRUE; break;}
                        }

                 if(mark==FALSE) continue;

                 for(j=0;j<num_atom;j++)
                        {
                         if(strcmp(atom[j].xtype,"N.3.h2.pi=0")) continue;
                         else if(Connection_1_4_Check(atom[i].id, atom[j].id)
                                 ==FALSE) continue;

                         num++;

			 if(flag)
				{
                                 atom[i].logp+=((LOGP_AMINO_ACID)/2.0);
                                 atom[j].logp+=((LOGP_AMINO_ACID)/2.0);
                                }

                         break;
                        }

                 if(num!=0) break;

                 for(j=0;j<num_atom;j++)
                        {
                         if(strcmp(atom[j].type,"N.ar")) continue;
                         else if(Connection_1_4_Check(atom[i].id, atom[j].id)
                                 ==FALSE) continue;

                         num++;

			 if(flag)
				{
                                 atom[i].logp+=((LOGP_AMINO_ACID)/2.0);
                                 atom[j].logp+=((LOGP_AMINO_ACID)/2.0);
                                }

                         break;
                        }

                 if(num!=0) break;
                }

        return (float)num;
}

float Molecule::Count_Sulfonic_Acid(int flag)
{
	int i,j,num,tmp1,tmp2;

	num=0;

	for(i=0;i<num_atom;i++)
		{
		 if(strcmp(atom[i].type,"S.o2")) continue;
		 else if(atom[i].ring!=0) continue;

		 tmp1=Adjacent_Aromatic_Check(atom[i].id);
		 if(tmp1==FALSE) continue;

		 for(j=0;j<num_atom;j++)
			{
			 if(strcmp(atom[j].xtype,"N.3.h2.pi=1")) continue;

			 tmp2=Adjacent_Aromatic_Check(atom[j].id);
			 if(tmp2==FALSE) continue;

			 if(Connection_1_6_Check(atom[i].id,atom[j].id)==FALSE) continue;
			 else if(Connection_1_4_Check(tmp1,tmp2)==FALSE) continue;

			 num++; 

			 if(flag)
				{
			 	 atom[i].logp+=(LOGP_SULFONIC_ACID/2.0);
			 	 atom[j].logp+=(LOGP_SULFONIC_ACID/2.0);
				}

			 break;
			}
		}

	return (float)num;
}

void Molecule::Write_Out_Log(FILE *fp) const
{
	extern LogP_Factor logp_factor[10];
        int i,mark;

	fprintf(fp,"#\n");
	fprintf(fp,"# XLOGP v2.0 logP calculation: %s\n", Get_Time());
	fprintf(fp,"#\n");
	fprintf(fp,"# Name of molecule : %s\n", this->name);
	fprintf(fp,"# Molecular formula: %s\n", this->formula); 
	fprintf(fp,"# Molecular weight : %-7.2f\n", this->weight);
	fprintf(fp,"#\n");

	fprintf(fp,"----------------------------------------------\n");
	fprintf(fp,"ID   atom type                  contribution\n");
	fprintf(fp,"----------------------------------------------\n");

	for(i=0;i<num_atom;i++)
		{
		 if(!strcmp(atom[i].type,"H")) continue;
		 
		 fprintf(fp,"%-3d  ", atom[i].id);
		 fprintf(fp,"%-10s  ", atom[i].xtype);
		 fprintf(fp,"%-13s  ", "");
		 fprintf(fp,"%6.3f  ", atom[i].logp);
		 fprintf(fp,"\n");
		}

	fprintf(fp,"----------------------------------------------\n");

	mark=FALSE;

	for(i=0;i<10;i++)
		{
		 if(logp_factor[i].num==0) continue;
		 
		 fprintf(fp,"%-30s\t%6.3f\n", 
			 logp_factor[i].symbol, 
			 logp_factor[i].num*logp_factor[i].coeff);
		 mark++;
		}

	if(mark!=FALSE)
		fprintf(fp,"----------------------------------------------\n");

	fprintf(fp,"Calculated LogP =%6.2f\n\n", this->logp);

	fflush(fp);

	return;
}

