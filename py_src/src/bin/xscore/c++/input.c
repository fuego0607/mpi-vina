# include "xtool.h"

Input::Input()
{
 strcpy(function,"none");

 strcpy(input_file,"none");
 strcpy(output_file,"none");
 strcpy(log_file,"none");

 strcpy(receptor_file,"none");
 strcpy(cofactor_file,"none");
 strcpy(reference_file,"none");
 strcpy(ligand_file,"none");

 num_method=3;
 strcpy(apply_hpscore,"YES");
 strcpy(apply_hmscore,"YES");
 strcpy(apply_hsscore,"YES");
 strcpy(apply_pmfscore,"NO");
 strcpy(show_abs,"NO");

 hpscore_cvdw=0.004;
 hpscore_chb=0.054;
 hpscore_chp=0.009;
 hpscore_crt=-0.061;
 hpscore_c0=3.441;

 hmscore_cvdw=0.004;
 hmscore_chb=0.101;
 hmscore_chm=0.387;
 hmscore_crt=-0.097;
 hmscore_c0=3.567;

 hsscore_cvdw=0.004;
 hsscore_chb=0.073;
 hsscore_chs=0.004;
 hsscore_crt=-0.090;
 hsscore_c0=3.328;

 if(getenv("XSCORE_PARAMETER")==NULL) 
	{
	 puts("Warning: XSCORE_PARAMETER is not set ... use default setting"); 
	 strcpy(parameter_dir,"../parameter/");
	}
 else strcpy(parameter_dir, getenv("XSCORE_PARAMETER"));

 Check_Directory(parameter_dir);
}

Input::~Input()
{
 // destructor
}

void Input::Read_Inputs(char *filename)
{
        FILE *fp;
        char line[256],head[256];

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else sscanf(line,"%s",head);

                 if(!strcmp(head,"RECEPTOR_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",receptor_file);
                        }
                 else if(!strcmp(head,"LIGAND_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",ligand_file);
                        }
                 else if(!strcmp(head,"COFACTOR_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",cofactor_file);
                        }
		 else continue;
		}

	fclose(fp);

        return;
}

void Input::Show_Contents() const
{
        printf("FUNCTION = %s\n", this->function);
	printf("RECEPTOR_PDB_FILE = %s\n", this->receptor_file);
	printf("LIGAND_MOL2_FILE = %s\n", this->ligand_file);
	printf("COFACTOR_MOL2_FILE = %s\n", this->cofactor_file);

	return;
}

void Input::Missing_Parameter_Error(const char *name) const
{
	printf("\n");
	printf("Error: You have probably forgot telling me %s.\n", name);
	return;
}

void Input::Invalid_Parameter_Error(const char *name) const
{
	printf("\n");
        printf("Error: %s has an invalid value.\n", name);
	return;
}

// ***********************************************************************
// add '/' to the directory name automatically
// if the directory is set for output, the 'flag' needs to be non-zero
// latest update: 11/21/2003
// ***********************************************************************
void Input::Check_Directory(char *dirname, int flag)
{
	int len;
        FILE *fp;
        char filename[256],command[256];

	len=strlen(dirname);

	if(dirname[len-1]!='/') strcat(dirname,"/");

	if(flag==FALSE) return;

	// check if the directory exists

        strcpy(filename,dirname); strcat(filename,".");

        if((fp=fopen(filename,"r"))==NULL)
        	{
         	 strcpy(command, "mkdir ");
         	 strcat(command, dirname);
         	 system(command);
         	 printf("Directory '%s' does not exist ... ", dirname);
         	 printf("it is created.\n");
	 	 return;
        	}
	else 
		{
		 fclose(fp); return;
		}
}

DOCK_Input::DOCK_Input()
{
        strcpy(function,"DOCK");
}

DOCK_Input::~DOCK_Input()
{
}

void DOCK_Input::Read_Inputs(char *filename)
{
        FILE *fp;
        char line[256],head[256],entry[5];
	int count;

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

	count=0;

        for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else sscanf(line,"%s",head);

                 if(!strcmp(head,"WORK_DIRECTORY"))
                        {
                         sscanf(line,"%*s%s",work_directory);
                        }
                 else if(!strcmp(head,"NUM_OF_ENTRY"))
                        {
                         sscanf(line,"%*s%d",&num_entry);
                        }
		 else if(!strcmp(head,"PDB"))
			{
			 sscanf(line,"%*s%s",entry);
			 strcpy(pdb_entry[count],entry);
			 count++;
			}
		 else continue;
		}

	fclose(fp);

	if(count!=num_entry)
	{
	 printf("Warning: something wrong with the number of PDB entries\n");
	 exit(1);
	}

	Check_Directory(work_directory);

        return;
}

void DOCK_Input::Show_Contents() const
{
	int i;

	printf("WORK_DIRECTORY = %s\n", work_directory);
	printf("NUM_OF_ENTRY = %d\n", num_entry);

	for(i=0;i<num_entry;i++)
		{
		 printf("PDB %s\n", pdb_entry[i]);
		}

	return;
}

SCORE_Input::SCORE_Input()
{
 strcpy(function,"SCORE");

 num_hits=0; strcpy(hits_dir,"none");

 strcpy(apply_chemical_rules,"NO");
 max_weight=600.0; min_weight=200.0;
 max_logp=6.00; min_logp=1.00;
 max_hb_atom=8; min_hb_atom=2;
}

SCORE_Input::~SCORE_Input()
{
}

void SCORE_Input::Read_Inputs(char *filename)
{
        FILE *fp;
        char line[256],head[256];
        int num_error;

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else sscanf(line,"%s",head);

                 if(!strcmp(head,"RECEPTOR_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",receptor_file);
                        }
                 else if(!strcmp(head,"COFACTOR_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",cofactor_file);
                        }
                 else if(!strcmp(head,"REFERENCE_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",reference_file);
                        }
                 else if(!strcmp(head,"LIGAND_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",ligand_file);
                        }
                 else if(!strcmp(head,"OUTPUT_TABLE_FILE"))
                        {
                         sscanf(line,"%*s%s",output_file);
                        }
                 else if(!strcmp(head,"OUTPUT_LOG_FILE"))
                        {
                         sscanf(line,"%*s%s",log_file);
                        }
		 else if(!strcmp(head,"NUMBER_OF_HITS"))
                        {
                         sscanf(line,"%*s%d",&num_hits);
                        }
		 else if(!strcmp(head,"HITS_DIRECTORY"))
                        {
                         sscanf(line,"%*s%s",hits_dir);
                        }

		 if(!strcmp(head,"APPLY_HPSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_hpscore);
                        }
		 else if(!strcmp(head,"HPSCORE_CVDW"))
			{
			 sscanf(line,"%*s%f", &hpscore_cvdw);
			}
		 else if(!strcmp(head,"HPSCORE_CHB"))
			{
			 sscanf(line,"%*s%f", &hpscore_chb);
			}
		 else if(!strcmp(head,"HPSCORE_CHP"))
			{
			 sscanf(line,"%*s%f", &hpscore_chp);
			}
		 else if(!strcmp(head,"HPSCORE_CRT"))
			{
			 sscanf(line,"%*s%f", &hpscore_crt);
			}
		 else if(!strcmp(head,"HPSCORE_C0"))
			{
			 sscanf(line,"%*s%f", &hpscore_c0);
			}

                 if(!strcmp(head,"APPLY_HMSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_hmscore);
                        }
		 else if(!strcmp(head,"HMSCORE_CVDW"))
			{
			 sscanf(line,"%*s%f", &hmscore_cvdw);
			}
		 else if(!strcmp(head,"HMSCORE_CHB"))
			{
			 sscanf(line,"%*s%f", &hmscore_chb);
			}
		 else if(!strcmp(head,"HMSCORE_CHM"))
			{
			 sscanf(line,"%*s%f", &hmscore_chm);
			}
		 else if(!strcmp(head,"HMSCORE_CRT"))
			{
			 sscanf(line,"%*s%f", &hmscore_crt);
			}
		 else if(!strcmp(head,"HMSCORE_C0"))
			{
			 sscanf(line,"%*s%f", &hmscore_c0);
			}

		 if(!strcmp(head,"APPLY_HSSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_hsscore);
                        }
		 else if(!strcmp(head,"HSSCORE_CVDW"))
			{
			 sscanf(line,"%*s%f", &hsscore_cvdw);
			}
		 else if(!strcmp(head,"HSSCORE_CHB"))
			{
			 sscanf(line,"%*s%f", &hsscore_chb);
			}
		 else if(!strcmp(head,"HSSCORE_CHS"))
			{
			 sscanf(line,"%*s%f", &hsscore_chs);
			}
		 else if(!strcmp(head,"HSSCORE_CRT"))
			{
			 sscanf(line,"%*s%f", &hsscore_crt);
			}
		 else if(!strcmp(head,"HSSCORE_C0"))
			{
			 sscanf(line,"%*s%f", &hsscore_c0);
			}

		 if(!strcmp(head,"SHOW_ATOM_BIND_SCORE"))
			{
			 sscanf(line,"%*s%s",show_abs);
			}

		 if(!strcmp(head,"APPLY_CHEMICAL_RULES"))
                        {
                         sscanf(line,"%*s%s",apply_chemical_rules);
                        }
		 else if(!strcmp(head,"MAXIMAL_MOLECULAR_WEIGHT"))
                        {
                         sscanf(line,"%*s%f",&max_weight);
                        }
                 else if(!strcmp(head,"MINIMAL_MOLECULAR_WEIGHT"))
                        {
                         sscanf(line,"%*s%f",&min_weight);
                        }
                 else if(!strcmp(head,"MAXIMAL_LOGP"))
                        {
                         sscanf(line,"%*s%f",&max_logp);
                        }
                 else if(!strcmp(head,"MINIMAL_LOGP"))
                        {
                         sscanf(line,"%*s%f",&min_logp);
                        }
                 else if(!strcmp(head,"MAXIMAL_HB_ATOM"))
                        {
                         sscanf(line,"%*s%d",&max_hb_atom);
                        }
                 else if(!strcmp(head,"MINIMAL_HB_ATOM"))
                        {
                         sscanf(line,"%*s%d",&min_hb_atom);
                        }
		 else continue;
		}

	fclose(fp);

	num_error=0;

        if(!strcmp(receptor_file,"none"))
                {
                 Missing_Parameter_Error("RECEPTOR_PDB_FILE");
                 num_error++;
                }
        if(!strcmp(ligand_file,"none"))
                {
                 Missing_Parameter_Error("LIGAND_MOL2_FILE");
                 num_error++;
                }
        if(!strcmp(output_file,"none"))
                {
                 Missing_Parameter_Error("OUTPUT_TABLE_FILE");
                 num_error++;
                }
	if(!strcmp(parameter_dir,"none"))
                {
                 Missing_Parameter_Error("PARAMETER_DIRECTORY");
                 num_error++;
                }

	if(!strncasecmp(show_abs,"Y",1))
		{
                 strcpy(show_abs,"YES");
		}
        else if(!strncasecmp(show_abs,"N",1))
		{
                 strcpy(show_abs,"NO");
		}
        else
                {
                 Invalid_Parameter_Error("SHOW_ATOM_BIND_SCORE");
                 num_error++;
                }

	if((num_hits>0)&&(!strcmp(hits_dir,"none")))
		{
		 Missing_Parameter_Error("HITS_DIRECTORY");
		 num_error++;
		}

	if(num_hits>0) Check_Directory(hits_dir,TRUE);
	else Check_Directory(hits_dir,FALSE);

	num_method=0;

	if(!strncasecmp(apply_hpscore,"Y",1)) 
		{
		 strcpy(apply_hpscore,"YES");
		 num_method++;
		}
        else if(!strncasecmp(apply_hpscore,"N",1)) 
		{
		 strcpy(apply_hpscore,"NO");
		}
        else
                {
                 Invalid_Parameter_Error("APPLY_HPSCORE");
                 num_error++;
                }

	if(!strncasecmp(apply_hmscore,"Y",1)) 
                {
                 strcpy(apply_hmscore,"YES");
                 num_method++;
                }
        else if(!strncasecmp(apply_hmscore,"N",1)) 
                {
                 strcpy(apply_hmscore,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("APPLY_HMSCORE");
                 num_error++;
                }

	if(!strncasecmp(apply_hsscore,"Y",1))
                {
                 strcpy(apply_hsscore,"YES");
                 num_method++;
                }
        else if(!strncasecmp(apply_hsscore,"N",1))
                {
                 strcpy(apply_hsscore,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("APPLY_HSSCORE");
                 num_error++;
                }

	if(num_method==0)
	{
	 puts("Error: no scoring function has been chosen.");
	 num_error++;
	}

	if(!strncasecmp(apply_chemical_rules,"Y",1))
		{
                 strcpy(apply_chemical_rules,"YES");
		}
        else if(!strncasecmp(apply_chemical_rules,"N",1))
		{
                 strcpy(apply_chemical_rules,"NO");
		}
        else
                {
                 Invalid_Parameter_Error("APPLY_CHEMICAL_RULES");
                 num_error++;
                }

	if(max_weight<min_weight)
                {
                 puts("Error: MAXIMAL_MOLECULAR_WEIGHT < MINIMAL_MOLECULAR_WEIGHT");
                 num_error++;
                }

        if(max_logp<min_logp)
                {
                 puts("Error: MAXIMAL_LOGP < MINIMAL_LOGP");
                 num_error++;
                }

        if(max_hb_atom<min_hb_atom)
                {
                 puts("Error: MAXIMAL_HB_ATOM < MINIMAL_HB_ATOM");
                 num_error++;
                }

	if(num_error!=0)
               	{
                 printf("\n");
                 printf("%d errors have been detected in the parameter file.\n",num_error);
                 printf("Please correct them and try again.\n");
                 exit(1);
                }

        return;
}

void SCORE_Input::Show_Contents() const
{
        printf("RECEPTOR_PDB_FILE = %s\n",receptor_file);
	printf("COFACTOR_MOL2_FILE = %s\n",cofactor_file);
	printf("REFERENCE_MOL2_FILE = %s\n",reference_file);
        printf("LIGAND_MOL2_FILE = %s\n",ligand_file);
        printf("OUTPUT_TABLE_FILE = %s\n",output_file);
	printf("OUTPUT_LOG_FILE = %s\n",log_file);

        printf("PARAMETER_DIRECTORY = %s\n",parameter_dir);

	printf("NUMBER_OF_HITS = %d\n",num_hits);
	printf("HITS_DIRECTORY = %s\n",hits_dir);
	printf("SHOW_ATOM_BIND_SCORE = %s\n",show_abs);

        printf("APPLY_HPSCORE = %s\n",apply_hpscore);
        printf("APPLY_HMSCORE = %s\n",apply_hmscore);
	printf("APPLY_HSSCORE = %s\n",apply_hsscore);
	printf("Number of methods = %d\n", num_method);
	
        printf("APPLY_CHEMICAL_RULES = %s\n",apply_chemical_rules);
        printf("MAXIMAL_MOLECULAR_WEIGHT = %6.1f\n",max_weight);
        printf("MINIMAL_MOLECULAR_WEIGHT = %6.1f\n",min_weight);
        printf("MAXIMAL_LOGP = %-6.2f\n",max_logp);
        printf("MINIMAL_LOGP = %-6.2f\n",min_logp);
        printf("MAXIMAL_HB_ATOM = %d\n",max_hb_atom);
        printf("MINIMAL_HB_ATOM = %d\n",min_hb_atom);

	return;
}

LOGP_Input::LOGP_Input()
{
	strcpy(function,"LOGP");

	strcpy(calculate_logp,"YES");
	strcpy(calculate_mw,"YES");
	strcpy(count_hb_atom,"YES");
	strcpy(count_rotor,"YES");
}

LOGP_Input::~LOGP_Input()
{
}

void LOGP_Input::Read_Inputs(char *filename)
{
        FILE *fp;
        char line[256],head[256];
        int num_error;

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else sscanf(line,"%s",head);

                 if(!strcmp(head,"INPUT_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s", input_file);
                        }
                 else if(!strcmp(head,"OUTPUT_LOG_FILE"))
                        {
                         sscanf(line,"%*s%s", output_file);
                        }
		 else if(!strcmp(head,"CALCULATE_LOGP"))
                        {
                         sscanf(line,"%*s%s", calculate_logp);
                        }
		 else if(!strcmp(head,"CALCULATE_MW"))
                        {
                         sscanf(line,"%*s%s", calculate_mw);
                        }
		 else if(!strcmp(head,"COUNT_HB_ATOM"))
                        {
                         sscanf(line,"%*s%s", count_hb_atom);
                        }
		 else if(!strcmp(head,"COUNT_ROTOR"))
                        {
                         sscanf(line,"%*s%s", count_rotor);
                        }
		 else continue;
		}

	fclose(fp);

	num_error=0;

        if(!strcmp(input_file,"none"))
                {
                 Missing_Parameter_Error("INPUT_MOL2_FILE");
                 num_error++;
                }

	if(!strcmp(parameter_dir,"none"))
                {
                 Missing_Parameter_Error("PARAMETER_DIRECTORY");
                 num_error++;
                }

	if(!strncasecmp(calculate_logp,"N",1))
		 strcpy(calculate_logp,"NO");
	else strcpy(calculate_logp,"YES");

	if(!strncasecmp(calculate_mw,"N",1))
                 strcpy(calculate_mw,"NO");
        else strcpy(calculate_mw,"YES");

	if(!strncasecmp(count_hb_atom,"N",1)) 
		strcpy(count_hb_atom,"NO");
        else strcpy(count_hb_atom,"YES"); 

	if(!strncasecmp(count_rotor,"N",1)) 
		strcpy(count_rotor,"NO");
        else strcpy(count_rotor,"YES");

	if(num_error!=0)
        {
         printf("\n");
       	 printf("%d errors have been detected in the input file.\n",num_error);
         printf("Please correct them and try again.\n");
         exit(1);
        }

        return;
}

void LOGP_Input::Show_Contents() const
{
        printf("INPUT_MOL2_FILE = %s\n",input_file);
        printf("OUTPUT_LOG_FILE = %s\n",output_file);
        printf("PARAMETER_DIRECTORY = %s\n",parameter_dir);

	return;
}

POCKET_Input::POCKET_Input()
{
	strcpy(function,"POCKET");

	strcpy(pocket_txt_file,"none");
	strcpy(grid_txt_file,"none");
	strcpy(pocket_pdb_file,"none");
	strcpy(box_pdb_file,"none");
	strcpy(shape_pdb_file,"none");
	strcpy(grid_pdb_file,"none");
	strcpy(pharmacophore_pdb_file,"none");

	probe_radius=4.00;
	shell_depth=8.00;
	min_feature_distance=3.50;
}

POCKET_Input::~POCKET_Input()
{
}

void POCKET_Input::Read_Inputs(char *filename)
{
        FILE *fp;
        char line[256],head[256];
        int num_error;

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else sscanf(line,"%s",head);

                 if(!strcmp(head,"RECEPTOR_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",receptor_file);
                        }
                 else if(!strcmp(head,"LIGAND_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",ligand_file);
                        }
                 else if(!strcmp(head,"SHELL_DEPTH"))
                        {
                         sscanf(line,"%*s%f",&shell_depth);
                        }
                 else if(!strcmp(head,"PROBE_RADIUS"))
                        {
                         sscanf(line,"%*s%f",&probe_radius);
                        }
		 else if(!strcmp(head,"POCKET_TXT_FILE"))
                        {
                         sscanf(line,"%*s%s",pocket_txt_file);
                        }
		 else if(!strcmp(head,"GRID_TXT_FILE"))
                        {
                         sscanf(line,"%*s%s",grid_txt_file);
                        }
		 else if(!strcmp(head,"POCKET_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",pocket_pdb_file);
                        }
		 else if(!strcmp(head,"BOX_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",box_pdb_file);
                        }
		 else if(!strcmp(head,"SHAPE_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",shape_pdb_file);
                        }
		 else if(!strcmp(head,"GRID_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",grid_pdb_file);
                        }
		 else if(!strcmp(head,"PHARMACOPHORE_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",pharmacophore_pdb_file);
                        }
		 else if(!strcmp(head,"MINIMAL_FEATURE_DISTANCE"))
                        {
                         sscanf(line,"%*s%f",&min_feature_distance);
                        }
		 else continue;
		}

	fclose(fp);

	num_error=0;

        if(!strcmp(receptor_file,"none"))
                {
                 Missing_Parameter_Error("RECEPTOR_PDB_FILE");
                 num_error++;
                }
        if(!strcmp(ligand_file,"none"))
                {
                 Missing_Parameter_Error("LIGAND_MOL2_FILE");
                 num_error++;
                }
	if(!strcmp(parameter_dir,"none"))
                {
                 Missing_Parameter_Error("PARAMETER_DIRECTORY");
                 num_error++;
                }

	if(probe_radius<1.50||probe_radius>10.0)
		{
		 Invalid_Parameter_Error("PROBE_RADIUS");
		 num_error++;
		}

	if(shell_depth<0.00||shell_depth>100.00)
		{
		 Invalid_Parameter_Error("SHELL_DEPTH");
		 num_error++;
		}

	if(min_feature_distance<0.00)
		{
		 Invalid_Parameter_Error("MINIMAL_FEATURE_DISTANCE");
                 num_error++;
		}

	if(num_error!=0)
        {
         printf("\n");
         printf("%d errors detected in the parameter file.\n",num_error);
         printf("Please correct them and try again.\n");
         exit(1);
        }

        return;
}

void POCKET_Input::Show_Contents() const
{
        printf("RECEPTOR_PDB_FILE = %s\n",receptor_file);
        printf("LIGAND_MOL2_FILE = %s\n",ligand_file);
        printf("PARAMETER_DIRECTORY = %s\n",parameter_dir);
	printf("PROBE_RADIUS = %f\n",probe_radius);
	printf("SHELL_DEPTH = %f\n",shell_depth);
	printf("POCKET_TXT_FILE = %s\n", pocket_txt_file);
        printf("GRID_TXT_FILE = %s\n", grid_txt_file);
	printf("POCKET_PDB_FILE = %s\n", pocket_pdb_file);
        printf("BOX_PDB_FILE = %s\n", box_pdb_file);
	printf("SHAPE_PDB_FILE = %s\n", shape_pdb_file);
	printf("GRID_PDB_FILE = %s\n", grid_pdb_file);
	printf("PHARMACOPHORE_PDB_FILE = %s\n", pharmacophore_pdb_file);
	printf("MINIMAL_FEATURE_DISTANCE =%f\n", min_feature_distance);
	
	return;
}

AUTODOCK_Input::AUTODOCK_Input()
{
 strcpy(function,"AUTODOCK");

 strcpy(autodock_mdb_location,"none");
 strcpy(autodock_job_name,"none");
 strcpy(autodock_pdbq_file,"none");
 strcpy(autodock_dlg_file,"none");
 autodock_num_conformation=0;
 autodock_rmsd_tolerance=2.00;
 strcpy(autodock_rank_conformation,"YES");
 strcpy(calculate_rmsd_distribution,"NO");
}

AUTODOCK_Input::~AUTODOCK_Input()
{
}

void AUTODOCK_Input::Read_Inputs(char *filename)
{
        FILE *fp;
        char line[256],head[256];
        int num_error;

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else sscanf(line,"%s",head);

                 if(!strcmp(head,"RECEPTOR_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",receptor_file);
                        }
                 else if(!strcmp(head,"LIGAND_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",ligand_file);
                        }
		 else if(!strcmp(head,"AUTODOCK_PDBQ_FILE"))
                        {
                         sscanf(line,"%*s%s",autodock_pdbq_file);
                        }
		 else if(!strcmp(head,"AUTODOCK_DLG_FILE"))
                        {
                         sscanf(line,"%*s%s",autodock_dlg_file);
                        }
		 else if(!strcmp(head,"AUTODOCK_JOB_NAME"))
                        {
                         sscanf(line,"%*s%s",autodock_job_name);
                        }
		 else if(!strcmp(head,"NUMBER_OF_CONFORMATION"))
                        {
                         sscanf(line,"%*s%d",&autodock_num_conformation);
                        }
		 else if(!strcmp(head,"CLUSTER_RMSD_TOLERANCE"))
                        {
                         sscanf(line,"%*s%f",&autodock_rmsd_tolerance);
                        }
		 else if(!strcmp(head,"APPLY_HPSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_hpscore);
                        }
                 else if(!strcmp(head,"APPLY_HMSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_hmscore);
                        }
                 else if(!strcmp(head,"APPLY_HSSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_hsscore);
                        }
                 else if(!strcmp(head,"APPLY_PMFSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_pmfscore);
                        }
		 else if(!strcmp(head,"SHOW_ATOM_BIND_SCORE"))
                        {
                         sscanf(line,"%*s%s",show_abs);
                        }
                 else if(!strcmp(head,"OUTPUT_CLUSTER_TABLE"))
                        {
                         sscanf(line,"%*s%s",output_file);
                        }
		 else if(!strcmp(head,"RANK_CONFORMATION"))
                        {
                         sscanf(line,"%*s%s",autodock_rank_conformation);
                        }
		 else if(!strcmp(head,"AUTODOCK_MDB_LOCATION"))
                        {
                         sscanf(line,"%*s%s",autodock_mdb_location);
                        }
		 else if(!strcmp(head,"CALCULATE_RMSD_DISTRIBUTION"))
                        {
                         sscanf(line,"%*s%s",calculate_rmsd_distribution);
                        }
		 else continue;
		}

	fclose(fp);

	num_error=0;

        if(!strcmp(receptor_file,"none"))
                {
                 Missing_Parameter_Error("RECEPTOR_PDB_FILE");
                 num_error++;
                }

        if(!strcmp(ligand_file,"none"))
                {
                 Missing_Parameter_Error("LIGAND_MOL2_FILE");
                 num_error++;
                }

	if(!strcmp(autodock_pdbq_file,"none"))
                {
                 Missing_Parameter_Error("AUTODOCK_PDBQ_FILE");
                 num_error++;
                }

	if(!strcmp(autodock_dlg_file,"none"))
                {
                 Missing_Parameter_Error("AUTODOCK_DLG_FILE");
                 num_error++;
                }

	if(!strcmp(autodock_job_name,"none"))
                {
                 Missing_Parameter_Error("AUTODOCK_JOB_NAME");
                 num_error++;
                }

	if(autodock_num_conformation==0)
		{
		 Invalid_Parameter_Error("NUMBER_OF_CONFORMATION");
		 num_error++;
		}

	if(autodock_rmsd_tolerance<=0.00)
		{
		 Invalid_Parameter_Error("CLUSTER_RMSD_TOLERANCE");
		 num_error++;
		}

	if(!strcmp(parameter_dir,"none"))
                {
                 Missing_Parameter_Error("PARAMETER_DIRECTORY");
                 num_error++;
                }

	num_method=0;

	if(!strncasecmp(apply_hpscore,"Y",1)) 
		{
		 strcpy(apply_hpscore,"YES");
		 num_method++;
		}
        else if(!strncasecmp(apply_hpscore,"N",1)) 
		{
		 strcpy(apply_hpscore,"NO");
		}
        else
                {
                 Invalid_Parameter_Error("APPLY_HPSCORE");
                 num_error++;
                }

	if(!strncasecmp(apply_hmscore,"Y",1)) 
                {
                 strcpy(apply_hmscore,"YES");
                 num_method++;
                }
        else if(!strncasecmp(apply_hmscore,"N",1)) 
                {
                 strcpy(apply_hmscore,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("APPLY_HMSCORE");
                 num_error++;
                }

	if(!strncasecmp(apply_hsscore,"Y",1))
                {
                 strcpy(apply_hsscore,"YES");
                 num_method++;
                }
        else if(!strncasecmp(apply_hsscore,"N",1))
                {
                 strcpy(apply_hsscore,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("APPLY_HSSCORE");
                 num_error++;
                }

	if(!strncasecmp(apply_pmfscore,"Y",1))
                {
                 strcpy(apply_pmfscore,"YES");
                }
        else if(!strncasecmp(apply_pmfscore,"N",1))
                {
                 strcpy(apply_pmfscore,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("APPLY_PMFSCORE");
                 num_error++;
                }

	if(!strncasecmp(show_abs,"Y",1))
                {
                 strcpy(show_abs,"YES");
                }
        else if(!strncasecmp(show_abs,"N",1))
                {
                 strcpy(show_abs,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("SHOW_ATOM_BIND_SCORE");
                 num_error++;
                }

	if(!strncasecmp(autodock_rank_conformation,"Y",1))
                {
                 strcpy(autodock_rank_conformation,"YES");
                }
        else if(!strncasecmp(autodock_rank_conformation,"N",1))
                {
                 strcpy(autodock_rank_conformation,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("RANK_CONFORMATION");
                 num_error++;
                }

	if(!strncasecmp(calculate_rmsd_distribution,"Y",1))
                {
                 strcpy(calculate_rmsd_distribution,"YES");
                }
        else if(!strncasecmp(calculate_rmsd_distribution,"N",1))
                {
                 strcpy(calculate_rmsd_distribution,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("CALCULATE_RMSD_DISTRIBUTION");
                 num_error++;
                }

        if(!strcmp(output_file,"none"))
                {
                 Missing_Parameter_Error("OUTPUT_CLUSTER_TABLE");
                 num_error++;
                }

	if(!strcmp(autodock_mdb_location,"none"))
                {
                 Missing_Parameter_Error("AUTODOCK_MDB_LOCATION");
                 num_error++;
                }

	if(autodock_num_conformation>0) Check_Directory(autodock_mdb_location,TRUE);
	else Check_Directory(autodock_mdb_location,FALSE);

	if(num_error!=0)
                {
                 printf("\n");
                 printf("%d errors have been detected in the parameter file.\n",num_error);
                 printf("Please correct them and try again.\n");
                 exit(1);
                }

        return;
}

void AUTODOCK_Input::Show_Contents() const
{
        printf("RECEPTOR_PDB_FILE = %s\n", receptor_file);
        printf("LIGAND_MOL2_FILE = %s\n", ligand_file);
	printf("AUTODOCK_PDBQ_FILE = %s\n", autodock_pdbq_file);
	printf("AUTODOCK_DLG_FILE = %s\n", autodock_dlg_file);
	printf("AUTODOCK_JOB_NAME = %s\n", autodock_job_name);
	printf("NUMBER_OF_CONFORMATION = %d\n", autodock_num_conformation);
	printf("CLUSTER_RMSD_TOLERANCE = %5.2f\n", autodock_rmsd_tolerance);
        printf("PARAMETER_DIRECTORY = %s\n", parameter_dir);
	printf("APPLY_HPSCORE = %s\n", apply_hpscore);
	printf("APPLY_HMSCORE = %s\n", apply_hmscore);
	printf("APPLY_HSSCORE = %s\n", apply_hsscore);
	printf("SHOW_ATOM_BIND_SCORE = %s\n", show_abs);
	printf("APPLY_PMFSCORE = %s\n", apply_pmfscore);
        printf("OUTPUT_CLUSTER_TABLE = %s\n", output_file);
	printf("RANK_CONFORMATION = %s\n", autodock_rank_conformation);
	printf("AUTODOCK_MDB_LOCATION = %s\n", autodock_mdb_location);
	printf("CALCULATE_RMSD_DISTRIBUTION = %s\n", calculate_rmsd_distribution);

	return;
}

FLEXX_Input::FLEXX_Input()
{
 strcpy(function,"FLEXX");

 strcpy(flexx_mdb_location,"none");
 strcpy(flexx_job_name,"none");
 flexx_num_conformation=0;
 strcpy(flexx_energy_table,"none");
 flexx_rmsd_tolerance=2.00;
 strcpy(flexx_rank_conformation,"YES");
}

FLEXX_Input::~FLEXX_Input()
{
}

void FLEXX_Input::Read_Inputs(char *filename)
{
        FILE *fp;
        char line[256],head[256];
        int num_error;

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else sscanf(line,"%s",head);

                 if(!strcmp(head,"RECEPTOR_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",receptor_file);
                        }
		 else if(!strcmp(head,"LIGAND_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",ligand_file);
                        }
                 else if(!strcmp(head,"FLEXX_MDB_LOCATION"))
                        {
                         sscanf(line,"%*s%s",flexx_mdb_location);
                        }
                 else if(!strcmp(head,"FLEXX_JOB_NAME"))
                        {
                         sscanf(line,"%*s%s",flexx_job_name);
                        }
		 else if(!strcmp(head,"NUMBER_OF_CONFORMATION"))
                        {
                         sscanf(line,"%*s%d",&flexx_num_conformation);
                        }
		 else if(!strcmp(head,"FLEXX_ENERGY_TABLE"))
                        {
                         sscanf(line,"%*s%s",flexx_energy_table);
                        }
		 else if(!strcmp(head,"CLUSTER_RMSD_TOLERANCE"))
                        {
                         sscanf(line,"%*s%f",&flexx_rmsd_tolerance);
                        }
		 else if(!strcmp(head,"APPLY_HPSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_hpscore);
                        }
                 else if(!strcmp(head,"APPLY_HMSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_hmscore);
                        }
                 else if(!strcmp(head,"APPLY_HSSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_hsscore);
                        }
		 else if(!strcmp(head,"APPLY_PMFSCORE"))
                        {
                         sscanf(line,"%*s%s",apply_pmfscore);
                        }
                 else if(!strcmp(head,"SHOW_ATOM_BIND_SCORE"))
                        {
                         sscanf(line,"%*s%s",show_abs);
                        }
                 else if(!strcmp(head,"OUTPUT_CLUSTER_TABLE"))
                        {
                         sscanf(line,"%*s%s",output_file);
                        }
		 else if(!strcmp(head,"RANK_CONFORMATION"))
                        {
                         sscanf(line,"%*s%s",flexx_rank_conformation);
                        }
		 else continue;
		}

	fclose(fp);

	num_error=0;

        if(!strcmp(receptor_file,"none"))
                {
                 Missing_Parameter_Error("RECEPTOR_PDB_FILE");
                 num_error++;
                }

	if(!strcmp(ligand_file,"none"))
                {
                 Missing_Parameter_Error("LIGAND_MOL2_FILE");
                 num_error++;
                }

        if(!strcmp(flexx_mdb_location,"none"))
                {
                 Missing_Parameter_Error("FLEXX_MDB_LOCATION");
                 num_error++;
                }

	Check_Directory(flexx_mdb_location);

        if(!strcmp(flexx_job_name,"none"))
                {
                 Missing_Parameter_Error("FLEXX_JOB_NAME");
                 num_error++;
                }

	if(flexx_num_conformation==0)
		{
		 Invalid_Parameter_Error("NUMBER_OF_CONFORMATION");
                 num_error++;
                }

	if(!strcmp(flexx_energy_table,"none"))
                {
                 Missing_Parameter_Error("FLEXX_ENERGY_TABLE");
                 num_error++;
                }

	if(flexx_rmsd_tolerance<=0.00)
		{
		 Invalid_Parameter_Error("CLUSTER_RMSD_TOLERANCE");
		 num_error++;
		}

	if(!strcmp(parameter_dir,"none"))
                {
                 Missing_Parameter_Error("PARAMETER_DIRECTORY");
                 num_error++;
                }

	num_method=0;

	if(!strncasecmp(apply_hpscore,"Y",1)) 
		{
		 strcpy(apply_hpscore,"YES");
		 num_method++;
		}
        else if(!strncasecmp(apply_hpscore,"N",1)) 
		{
		 strcpy(apply_hpscore,"NO");
		}
        else
                {
                 Invalid_Parameter_Error("APPLY_HPSCORE");
                 num_error++;
                }

	if(!strncasecmp(apply_hmscore,"Y",1)) 
                {
                 strcpy(apply_hmscore,"YES");
                 num_method++;
                }
        else if(!strncasecmp(apply_hmscore,"N",1)) 
                {
                 strcpy(apply_hmscore,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("APPLY_HMSCORE");
                 num_error++;
                }

	if(!strncasecmp(apply_hsscore,"Y",1))
                {
                 strcpy(apply_hsscore,"YES");
                 num_method++;
                }
        else if(!strncasecmp(apply_hsscore,"N",1))
                {
                 strcpy(apply_hsscore,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("APPLY_HSSCORE");
                 num_error++;
                }

	if(!strncasecmp(apply_pmfscore,"Y",1))
                {
                 strcpy(apply_pmfscore,"YES");
                }
        else if(!strncasecmp(apply_pmfscore,"N",1))
                {
                 strcpy(apply_pmfscore,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("APPLY_PMFSCORE");
                 num_error++;
                }

	if(!strncasecmp(show_abs,"Y",1))
		{
                 strcpy(show_abs,"YES");
		}
        else if(!strncasecmp(show_abs,"N",1))
		{
                 strcpy(show_abs,"NO");
		}
        else
                {
                 Invalid_Parameter_Error("SHOW_ATOM_BIND_SCORE");
                 num_error++;
                }

	if(!strncasecmp(flexx_rank_conformation,"Y",1))
                {
                 strcpy(flexx_rank_conformation,"YES");
                }
        else if(!strncasecmp(flexx_rank_conformation,"N",1))
                {
                 strcpy(flexx_rank_conformation,"NO");
                }
        else
                {
                 Invalid_Parameter_Error("RANK_CONFORMATION");
                 num_error++;
                }

        if(!strcmp(output_file,"none"))
                {
                 Missing_Parameter_Error("OUTPUT_CLUSTER_TABLE");
                 num_error++;
                }

	if(num_error!=0)
                {
                 printf("\n");
                 printf("%d errors have been detected in the parameter file.\n",num_error);
                 printf("Please correct them and try again.\n");
                 exit(1);
                }

	// this->Show_Contents();

        return;
}

void FLEXX_Input::Show_Contents() const
{
        printf("RECEPTOR_PDB_FILE = %s\n", receptor_file);
	printf("LIGAND_MOL2_FILE = %s\n", ligand_file);
        printf("FLEXX_MDB_LOCATION = %s\n", flexx_mdb_location);
        printf("FLEXX_JOB_NAME = %s\n", flexx_job_name);
	printf("NUMBER_OF_CONFORMATION = %d\n", flexx_num_conformation);
	printf("FLEXX_ENERGY_TABLE = %s\n", flexx_energy_table);
	printf("CLUSTER_RMSD_TOLERANCE = %4.2f\n", flexx_rmsd_tolerance);
        printf("PARAMETER_DIRECTORY = %s\n", parameter_dir);
        printf("APPLY_HPSCORE = %s\n", apply_hpscore);
        printf("APPLY_HMSCORE = %s\n", apply_hmscore);
	printf("APPLY_HSSCORE = %s\n", apply_hsscore);
	printf("Number of methods = %d\n", num_method);
	printf("APPLY_PMFSCORE = %s\n", apply_pmfscore);
	printf("SHOW_ATOM_BIND_SCORE = %s\n", show_abs);
        printf("OUTPUT_CLUSTER_TABLE = %s\n", output_file);
	printf("RANK_CONFORMATION = %s\n", flexx_rank_conformation);

	return;
}

SURFACE_Input::SURFACE_Input()
{
	strcpy(function,"SURFACE");
	probe_radius=1.50;
	dot_density=4.0;
}

SURFACE_Input::~SURFACE_Input()
{
}

void SURFACE_Input::Read_Inputs(char *filename)
{
        FILE *fp;
        char line[256],head[256];
        int num_error;

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else sscanf(line,"%s",head);

                 if(!strcmp(head,"INPUT_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",input_file);
                        }
		 else if(!strcmp(head,"SURFACE_PDB_FILE"))
                        {
                         sscanf(line,"%*s%s",output_file);
                        }
		 else if(!strcmp(head,"PROBE_RADIUS"))
                        {
                         sscanf(line,"%*s%f",&probe_radius);
                        }
		 else if(!strcmp(head,"DOT_DENSITY"))
                        {
                         sscanf(line,"%*s%f",&dot_density);
                        }
		 else continue;
		}

	fclose(fp);

	num_error=0;

        if(!strcmp(input_file,"none"))
                {
                 Missing_Parameter_Error("LIGAND_MOL2_FILE");
                 num_error++;
                }

        if(!strcmp(output_file,"none"))
                {
                 Missing_Parameter_Error("SURFACE_PDB_FILE");
                 num_error++;
                }

	if(num_error!=0)
       {
        printf("\n");
        printf("%d errors in the parameter file.\n",num_error);
        printf("Please correct them and try again.\n");
        exit(1);
       }

        return;
}

void SURFACE_Input::Show_Contents() const
{
        printf("INPUT_MOL2_FILE = %s\n", input_file);
        printf("SURFACE_PDB_FILE = %s\n", output_file);
        printf("PARAMETER_DIRECTORY = %s\n", parameter_dir);
	printf("PROBE_RADIUS = %f\n", probe_radius);
	printf("DOT_DENSITY = %f\n", dot_density);

	return;
}

DATABASE_Input::DATABASE_Input()
{
        strcpy(function,"DATABASE");

	strcpy(output_00,"database_00.mol2");
	strcpy(output_01,"database_01.mol2");
	strcpy(output_02,"database_02.mol2");
	strcpy(output_03,"database_03.mol2");
	strcpy(output_04,"database_04.mol2");
	strcpy(output_05,"database_05.mol2");
	strcpy(output_06,"database_06.mol2");
	strcpy(output_07,"database_07.mol2");
	strcpy(output_08,"database_08.mol2");
	strcpy(output_09,"database_09.mol2");
	strcpy(output_10,"database_10.mol2");
	strcpy(output_error,"database_error.mol2");
}

DATABASE_Input::~DATABASE_Input()
{
}

void DATABASE_Input::Read_Inputs(char *filename)
{
        FILE *fp;
        char line[256],head[256];
        int num_error;

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        for(;;)
                {
                 if(fgets(line,256,fp)==NULL) break;
                 else if(line[0]=='#') continue;
                 else if(Blank_Line_Check(line)==TRUE) continue;
                 else sscanf(line,"%s",head);

                 if(!strcmp(head,"INPUT_MOL2_FILE"))
                        {
                         sscanf(line,"%*s%s",input_file);
                        }
                 else if(!strcmp(head,"OUTPUT_0_100_FILE"))
                        {
                         sscanf(line,"%*s%s",output_00);
                        }
		 else if(!strcmp(head,"OUTPUT_100_200_FILE"))
                        {
                         sscanf(line,"%*s%s",output_01);
                        }
		 else if(!strcmp(head,"OUTPUT_200_300_FILE"))
                        {
                         sscanf(line,"%*s%s",output_02);
                        }
		 else if(!strcmp(head,"OUTPUT_300_400_FILE"))
                        {
                         sscanf(line,"%*s%s",output_03);
                        }
		 else if(!strcmp(head,"OUTPUT_400_500_FILE"))
                        {
                         sscanf(line,"%*s%s",output_04);
                        }
		 else if(!strcmp(head,"OUTPUT_500_600_FILE"))
                        {
                         sscanf(line,"%*s%s",output_05);
                        }
		 else if(!strcmp(head,"OUTPUT_600_700_FILE"))
                        {
                         sscanf(line,"%*s%s",output_06);
                        }
		 else if(!strcmp(head,"OUTPUT_700_800_FILE"))
                        {
                         sscanf(line,"%*s%s",output_07);
                        }
		 else if(!strcmp(head,"OUTPUT_800_900_FILE"))
                        {
                         sscanf(line,"%*s%s",output_08);
                        }
		 else if(!strcmp(head,"OUTPUT_900_1000_FILE"))
                        {
                         sscanf(line,"%*s%s",output_09);
                        }
		 else if(!strcmp(head,"OUTPUT_1000+_FILE"))
                        {
                         sscanf(line,"%*s%s",output_10);
                        }
		 else if(!strcmp(head,"OUTPUT_ERROR_FILE"))
                        {
                         sscanf(line,"%*s%s",output_error);
                        }
		 else continue;
		}

	fclose(fp);

	num_error=0;

	if(!strcmp(input_file,"none"))
                {
                 Missing_Parameter_Error("INPUT_MOL2_FILE");
                 num_error++;
                }

	if(!strcmp(parameter_dir,"none"))
                {
                 Missing_Parameter_Error("PARAMETER_DIRECTORY");
                 num_error++;
                }

	if(num_error!=0)
       {
        printf("\n");
        printf("%d errors have been detected in the parameter file.\n",num_error);
        printf("Please correct them and try again.\n");
        exit(1);
       }

        return;
}

void DATABASE_Input::Show_Contents() const
{
	printf("INPUT_MOL2_FILE = %s\n", input_file);
        printf("PARAMETER_DIRECTORY = %s\n", parameter_dir);
	printf("OUTPUT_0_100_FILE = %s\n", output_00);
	printf("OUTPUT_100_200_FILE = %s\n", output_01);
	printf("OUTPUT_200_300_FILE = %s\n", output_02);
	printf("OUTPUT_300_400_FILE = %s\n", output_03);
	printf("OUTPUT_400_500_FILE = %s\n", output_04);
	printf("OUTPUT_500_600_FILE = %s\n", output_05);
	printf("OUTPUT_600_700_FILE = %s\n", output_06);
	printf("OUTPUT_700_800_FILE = %s\n", output_07);
	printf("OUTPUT_800_900_FILE = %s\n", output_08);
	printf("OUTPUT_900_1000_FILE = %s\n", output_09);
	printf("OUTPUT_1000+_FILE = %s\n", output_10);
	printf("OUTPUT_ERROR_FILE = %s\n", output_error);

        return;
}

