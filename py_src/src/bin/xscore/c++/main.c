/******************************************************************************

#     #         ####### ####### ####### #
 #   #             #    #     # #     # #
  # #              #    #     # #     # #
   #     #####     #    #     # #     # #
  # #              #    #     # #     # #
 #   #             #    #     # #     # #
#     #            #    ####### ####### #######  v.2004.7

Structure-Based Drug Design Toolkits, developed by Dr. Renxiao Wang
******************************************************************************/

# include "xtool.h"

Input *input=NULL;		// input parameters
ForceField *ff=NULL;            // the force field
Protein *protein=NULL;		// the receptor
Ligand *ligand=NULL;		// the ligand

char timestring[256];           // global string for storing local time

void Check_Function(char *filename, char *function);

int main(int argc,char *argv[])
{
	extern Input *input;
	extern ForceField *ff;
	extern Protein *protein;
	extern Ligand *ligand;

	printf("\nX-Score starts to run ... %s\n\n", Get_Time());

        if(argc==2) 	       		// major functions:  
		{
		 Xtool_Functions(argc,argv);
		}
        else if(argc>2) 	  	// misc utilities: 
		{
		 Xtool_Utilities(argc,argv);
		}
	else
                {
                 puts("Error: wrong usage"); 
		 puts("Synopsis: xscore input_file");
		 puts("       or xscore -function ... ... ...");
		 exit(1);
                }

	delete input; delete ff; delete protein; delete ligand; 

	printf("\nDone ... %s\n\n", Get_Time());

	return 0;
}

// ***************************************************************************
// perform major functions implemented in X-TOOL
// latest update: 11/06/2003
// ***************************************************************************
int Xtool_Functions(int argc, char *argv[])
{
 extern Input *input;
 char function[80];

 Check_Function(argv[1], function);

 if(!strcasecmp(function,"SCORE"))
	{
	 input=new SCORE_Input; if(input==NULL) Memory_Allocation_Error();
	}
 else if(!strcasecmp(function,"LOGP"))
        {
         input=new LOGP_Input; if(input==NULL) Memory_Allocation_Error();
        }
 else 
	{
	 puts("Error: no specific function assigned in the input file.");
	 puts("The program needs to know what to do.");
	 return FALSE;
	}

 // now read the information in the input file

 printf("Now reading the input file '%s' ...\n", argv[1]);
 input->Read_Inputs(argv[1]);
 // input->Show_Contents();

 if(!strcasecmp(function,"SCORE")) {Xtool_Score(); return TRUE;}
 else if(!strcasecmp(function,"LOGP")) {Xtool_LogP(); return TRUE;}
 else return FALSE;
}

void Check_Function(char *filename, char *function)
{
	FILE *fp;
	char buf[256],head[256];

	strcpy(function,"none");

	if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
                 else if(Blank_Line_Check(buf)==TRUE) continue;

                 strcpy(head,""); sscanf(buf,"%s",head);

                 if(!strcasecmp(head,"FUNCTION"))
                        {
                         sscanf(buf,"%*s%s", function);
                        }
		 else continue;
		}

	fclose(fp);

	return;
}

// ***************************************************************************
// perform misc functions
// latest update: 11/06/2003
// ***************************************************************************
int Xtool_Utilities(int argc, char *argv[])
{
 extern Input *input;
 int i,len;
 char function[80];

 strcpy(function,argv[1]); len=strlen(function); 
 for(i=0;i<len;i++) function[i]=tolower(function[i]);

 if(strstr(function,"score"))
	{
	 if(argc>=4)
		{
	 	 input=new Input; 
		 if(input==NULL) Memory_Allocation_Error();

         	 strcpy(input->receptor_file,argv[2]);

         	 if(argc==4)
			{
		 	 strcpy(input->cofactor_file,"none");
		 	 strcpy(input->ligand_file,argv[3]);
			}
	 	 else
			{
		 	 strcpy(input->cofactor_file,argv[3]);
		 	 strcpy(input->ligand_file,argv[4]);
			}

		 strcpy(input->log_file,"xscore.log");
         	 strcpy(input->apply_hpscore,"YES");
         	 strcpy(input->apply_hmscore,"YES");
         	 strcpy(input->apply_hsscore,"YES");
         	 input->num_method=3;
	 	 strcpy(input->apply_pmfscore,"NO");

	 	 Xtool_Score_Shortcut(); return TRUE;
		}
	 else
		{
		 puts("Error: wrong usage!");
		 puts("Synopsis: xscore -score protein_file [cofactor_file] ligand_file");
		 return FALSE; 
		}
	}

 if(strstr(function,"logp"))
	{
	 if(argc>=3)
		{
	 	 input=new Input; 
		 if(input==NULL) Memory_Allocation_Error();

         	 strcpy(input->input_file,argv[2]);
		 strcpy(input->output_file,"xlogp.log");

	 	 Xtool_LogP_Shortcut(); return TRUE;
		}
	 else
		{
		 puts("Error: wrong usage!");
		 puts("Synopsis: xscore -logp Mol2_file");
		 return FALSE; 
		}
	}

 if(strstr(function,"fixpdb"))
	{
	 if(argc!=4)
		{
		 puts("Error: wrong usage!");
		 puts("Synopsis: xscore -fixpdb input_file output_file");
		 return FALSE; 
		}
	 else
		{
		 input=new Input; if(input==NULL) Memory_Allocation_Error();
		 strcpy(input->input_file,argv[2]);
		 strcpy(input->output_file,argv[3]);
		 
		 Xtool_Fix_PDB(); return TRUE;
		}
	}

 if(strstr(function,"fixmol2"))
	{
	 if(argc!=4)
		{
		 puts("Error: wrong usage!");
		 puts("Synopsis: xscore -fixmol2 input_file output_file");
		 return FALSE; 
		}
	 else
		{
		 input=new Input; if(input==NULL) Memory_Allocation_Error();
		 strcpy(input->input_file,argv[2]);
		 strcpy(input->output_file,argv[3]);

		 Xtool_Fix_Mol2(); return TRUE;
		}
	}

 puts("Error: No valid function is assigned.");
 puts("Synopsis: xscore -function ... ...");

 return FALSE;
}


