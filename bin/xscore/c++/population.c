# include "xtool.h"

Population::Population(char *filename)
{
	extern Input *input;
	int i;

        num_member=Check_Mol2_File(filename);
	if(num_member==0) Mol2_Format_Error(filename);

	record=new Record[num_member];
	if(record==NULL) Memory_Allocation_Error();

	if(((SCORE_Input *)input)->num_hits>num_member) num_output=num_member;
	else num_output=((SCORE_Input *)input)->num_hits;

	member=NULL; member=new Ligand[num_output];
	if(member==NULL) Memory_Allocation_Error();

	for(i=0;i<num_output;i++) member[i].valid=0;
}

Population::~Population()
{
	if(member) delete [] member;
	if(record) delete [] record;
}

void Population::Score_Members()
{
	extern Input *input;
	FILE *fp;
	int i,j,id;
	float worst;
	Ligand mol;

	if((fp=fopen(input->ligand_file,"r"))==NULL) Open_File_Error(input->ligand_file);

        for(i=0;i<num_member;i++)
                {
		 mol.Clear();
		 if(mol.Read_From_Mol2(fp)==FALSE) continue;

		 printf("Now processing %s ...\n", mol.name);

                 // if(mol.Value_Atom()==FALSE) continue;
		 mol.Value_Atom();

                 record[i].id=i+1;
                 strcpy(record[i].name,mol.name);
                 strcpy(record[i].formula,mol.formula);
                 record[i].weight=mol.weight;
                 record[i].logp=mol.logp;
                 record[i].num_hb_atom=mol.num_hb_atom;

                 if(mol.Chemical_Viability_Check()==FALSE)
                        {
			 record[i].hpscore=0.000;
			 record[i].hmscore=0.000;
			 record[i].hsscore=0.000;
			 record[i].score=0.000;
                         record[i].valid=0;
                         continue;
                        }

               	 mol.Calculate_Binding_Score();

	 	 record[i].hpscore=mol.pkd1;
	 	 record[i].hmscore=mol.pkd2;
	 	 record[i].hsscore=mol.pkd3;
	 	 record[i].score=mol.bind_score;
	 	 record[i].valid=1;

		 if(num_output<=0) continue;

		 // now check if the new molecule is qualified for final output 
		 // this is a dynamic stack

		 worst=mol.bind_score; id=0;

		 for(j=0;j<num_output;j++)
			{
			 if(member[j].bind_score>=worst) continue;
			 else {worst=member[j].bind_score; id=j+1;}
			}

		 if(id==0) continue;
		 else {mol.valid=1; member[id-1]=mol;}
                }

	fclose(fp);

	printf("Totally %d molecules are processed.\n", num_member);

	return;
}

void Population::Rank_Members()
{
	int i,j;
	Record temp;
	Ligand tmp_mol;

	// rank the records first

	for(i=0;i<num_member-1;i++)
	for(j=i+1;j<num_member;j++)
		{
		 if(record[i].score>=record[j].score) continue;
		 else {SWAP(record[i],record[j]);}
		}

	for(i=0;i<num_member;i++) record[i].rank=i+1;

	// rank the top molecules then
	// note that there are only num_output real molecules 

	for(i=0;i<num_output-1;i++)
        for(j=i+1;j<num_output;j++)
                {
                 if(member[i].bind_score>=member[j].bind_score) continue;
                 else
                        {
                         tmp_mol=member[i];
                         member[i]=member[j];
                         member[j]=tmp_mol;
                        }
                }

	return;
}

void Population::Output_Results()
{
	extern Input *input;
	FILE *fp;
	int i;
	char filename[256];

	// output the table

	if((fp=fopen(input->output_file,"w"))==NULL) Open_File_Error(input->output_file);

	printf("Now writing the results to '%s' ...\n", input->output_file);

	fprintf(fp,"%s,", "Rank");
	fprintf(fp,"%s,", "Formula");
	fprintf(fp,"%s,", "MW");
	fprintf(fp,"%s,", "LogP");
	fprintf(fp,"%s,", "HPScore");
	fprintf(fp,"%s,", "HMScore");
	fprintf(fp,"%s,", "HSScore");
	fprintf(fp,"%s,","Average");
	fprintf(fp,"%s\n", "Name");

	for(i=0;i<num_member;i++)
		{
		 fprintf(fp,"%-5d ", record[i].rank);
                 fprintf(fp,"%-16s ", record[i].formula);
		 fprintf(fp,"%6.1f  ", record[i].weight);
		 fprintf(fp,"%6.2f  ", record[i].logp);
		 fprintf(fp,"%5.2f  ", record[i].hpscore);
		 fprintf(fp,"%5.2f  ", record[i].hmscore);
		 fprintf(fp,"%5.2f  ", record[i].hsscore);
		 fprintf(fp,"%5.2f  ", record[i].score);
		 fprintf(fp,"%s\n", record[i].name);
		}

	fclose(fp);

	// output top candidates ----------------------------------------

	if(((SCORE_Input *)input)->num_hits<=0) return;

	printf("Now extracting top hits to '%s' ...\n", 
		((SCORE_Input *)input)->hits_dir);

	for(i=0;i<num_output;i++)
	{
	 if(member[i].valid==0) continue;
	 sprintf(filename,"%sNo%d.mol2",((SCORE_Input *)input)->hits_dir,i+1);
	 member[i].Write_Out_Mol2(filename);
	}

	if((fp=fopen(input->log_file,"w"))==NULL) Open_File_Error(input->log_file);

	printf("Now writing scoring information of top hits to '%s' ...\n", input->log_file);

	for(i=0;i<num_output;i++)
		{
		 if(member[i].valid==0) continue;
		 member[i].Write_Out_Log(fp);
		}

	fclose(fp); return;
}

