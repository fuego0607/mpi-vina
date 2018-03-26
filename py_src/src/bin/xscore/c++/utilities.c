# include "xtool.h"

// ***************************************************************************
// read in a PDB file and write out a XTOOL/PDB file (polar-H only)
// latest update: 11/06/2003
// ***************************************************************************
int Xtool_Fix_PDB()
{
 extern Input *input;
 extern ForceField *ff;
 extern Protein *protein;
 int i;

 ff=new ForceField(input->parameter_dir);
 if(ff==NULL) Memory_Allocation_Error("Xtool_Fix_PDB()");

 protein=new Protein;
 if(protein==NULL) Memory_Allocation_Error("Xtool_Fix_PDB()");

 printf("Now reading the input PDB file from '%s' ...\n", input->input_file);
 protein->Read_From_PDB(input->input_file); 

 printf("Now analyzing the structure ...\n");

 protein->Value_Atom(0);                      // do not calculate HB root

 // disable non-polar hydrogen atoms

 for(i=0;i<protein->num_atom;i++)
        {
         if(strcmp(protein->atom[i].xtype,"H")) continue;
         else protein->atom[i].valid=0;
        }

 protein->Refine_Structure();
 protein->Rearrange_IDs();

 printf("Now writing the fixed PDB file to '%s' ...\n", input->output_file);
 protein->Write_Out_PDB(input->output_file);

 return TRUE;
}

// ***************************************************************************
// read in a Mol2 file and write out a XTOOL/Mol2 file
// latest update: 11/06/2003
// ***************************************************************************
int Xtool_Fix_Mol2()
{
 extern Input *input;
 extern ForceField *ff;
 FILE *fin,*fout;
 int i,total,valid;
 Molecule mol;

 ff=new ForceField(input->parameter_dir);
 if(ff==NULL) Memory_Allocation_Error("Xtool_Fix_Mol2()");

 printf("Now reading the input Mol2 file from '%s' ...\n", input->input_file);

 total=Check_Mol2_File(input->input_file);

 if(total==0)
	{
	 puts("Error: no valid ligand molecule is given.");
	 exit(1);
	}

 if((fin=fopen(input->input_file,"r"))==NULL) 
	Open_File_Error(input->input_file);

 if((fout=fopen(input->output_file,"w"))==NULL)
	Open_File_Error(input->output_file);

 valid=0;

 for(i=0;i<total;i++)
	{
	 if(mol.Read_From_Mol2(fin)==FALSE) continue;

	 printf("Now processing %s ... \n", mol.name);

	 if(mol.Value_Atom()==FALSE)
		{
		 printf("This molecule is skipped.\n");
		}
	 else 
		{
		 mol.Write_Out_Mol2(fout); valid++;
		}
	}

 printf("%d molecules have been processed; %d are valid.\n", total, valid);

 printf("A new Mol2 file '%s' has been created.\n", input->output_file);

 fclose(fin); fclose(fout); return TRUE;
}

