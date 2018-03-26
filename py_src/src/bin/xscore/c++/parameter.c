# include "xtool.h"

ForceField::ForceField(char *dirname)
{
	char filename[256];

	num_restype=0; residue=NULL;
	num_atomtype=0;	atom=NULL;
	num_xatomtype=0; xatom=NULL;
	num_bondtype=0; bond=NULL;
	num_torstype=0; torsion=NULL;
	num_pmftype=0; pmf=NULL;
	num_sdot_type=0; sdot=NULL;
	num_vdot_type=0; vdot=NULL;
	num_hbtype=0; hb_pmf=NULL;
	num_wpmftype=0; wpmf=NULL;

	strcpy(filename,dirname); strcat(filename,"RESIDUE_DEF_XTOOL");
	Read_RESIDUE_DEF(filename);

	strcpy(filename,dirname); strcat(filename,"ATOM_DEF_XTOOL");
	Read_ATOM_DEF(filename);

	strcpy(filename,dirname); strcat(filename,"ATOM_DEF_XLOGP");
        Read_XATOM_DEF(filename);

	// strcpy(filename,dirname); strcat(filename,"BOND_DEF_SYBYL");
	// Read_BOND_DEF(filename);

	// strcpy(filename,dirname); strcat(filename,"TORSION_DEF_SYBYL");
	// Read_TORSION_DEF(filename);

	// strcpy(filename,dirname); strcat(filename,"PMF_DEF_XTOOL");
	// Read_PMF_DEF(filename);

	strcpy(filename,dirname); strcat(filename,"SURFACE_DEF_XTOOL");
        Read_SURFACE_DEF(filename);

	strcpy(filename,dirname); strcat(filename,"VOLUME_DEF_XTOOL");
	Read_VOLUME_DEF(filename);
}

ForceField::~ForceField()
{
	if(residue) delete [] residue;
	if(atom) delete [] atom;
	if(xatom) delete [] xatom;
	if(bond) delete [] bond;
	if(torsion) delete [] torsion;
	if(pmf) delete [] pmf;
	if(sdot) delete [] sdot;
	if(vdot) delete [] vdot;
	if(hb_pmf) delete [] hb_pmf;
	if(wpmf) delete [] wpmf;
}

void ForceField::Read_RESIDUE_DEF(char *filename)
{
	FILE *fp;
	int i,num,count;
	char buf[256],head[256];

	if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        // first, determine num_restype

        count=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
                 else if(Blank_Line_Check(buf)==true) continue;
		 else
			{
			 sscanf(buf,"%s", head);
			 if(strcmp(head,"RESI")) continue;
			 else count++;
			}
                }

        rewind(fp);

        num_restype=count;

        residue=new Residue_Def[num_restype];
        if(residue==NULL) Memory_Allocation_Error();

	// second, read the residue templates

	count=0;

	for(;;)
	{
	 if(fgets(buf,256,fp)==NULL) break;
	 else if(Blank_Line_Check(buf)==true) continue;
	 else if(buf[0]=='#') continue;
		 
	 sscanf(buf,"%s",head);

	 if(!strcmp(head,"END")) break;
	 else if(!strcmp(head,"RESI"))
		{
		 num=0;
		 sscanf(buf,"%*s%s%d%d",
			residue[count].name,
		        &residue[count].num_atom,
			&residue[count].num_bond);

		 // read the ATOM information first

		 num=residue[count].num_atom;

		 for(i=0;i<num;i++)
			{
		 	 fgets(buf,256,fp);
		 	 sscanf(buf,"%*s%s%s%s%f%f%f%s%f%f%d%s",
			 	residue[count].atom[i].name,
				residue[count].atom[i].type,
				residue[count].atom[i].xtype,
		       	       &residue[count].atom[i].r,
		               &residue[count].atom[i].eps,
                               &residue[count].atom[i].q,
				residue[count].atom[i].hb,
		       	       &residue[count].atom[i].logp,
		       	       &residue[count].atom[i].solv,
		       	       &residue[count].atom[i].ring,
				residue[count].atom[i].pmftype);
			}

		 // then read the BOND information

		 num=residue[count].num_bond;

		 if(num>0)
		{
		 for(i=0;i<num;i++)
                        {
                         fgets(buf,256,fp);
                         sscanf(buf,"%*s%s%s%s",
                                residue[count].bond[i].atom1,
                                residue[count].bond[i].atom2,
                                residue[count].bond[i].type);
                        }
		}

		 count++;
		}
	}

	fclose(fp);

        if(count!=num_restype) Read_File_Error(filename);

	return;
}

void ForceField::Read_ATOM_DEF(char *filename)
{
	FILE *fp;
	int count;
	char buf[256];

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

	// first, determine num_atomtype

        count=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
		 else if(Blank_Line_Check(buf)==true) continue;
		 else count++;
		}

	rewind(fp);

	num_atomtype=count;

	atom=new Atom_Def[num_atomtype];
	if(atom==NULL) Memory_Allocation_Error();

	// second, read parameters for each atom type

        count=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
		 else if(Blank_Line_Check(buf)==true) continue;
		 else
			{
			 sscanf(buf,"%*d%s%f%f%f%f%s", 
				 atom[count].type, 
				&atom[count].weight, 
				&atom[count].r,
				&atom[count].eps,
				&atom[count].q,
				 atom[count].hb);
			 count++;
			}
		}

	fclose(fp);

	if(num_atomtype!=count) Read_File_Error(filename);

	return;
}

void ForceField::Read_XATOM_DEF(char *filename)
{
	FILE *fp;
	int count;
	char buf[256];

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

	// first, determine num_xatomtype

        count=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
		 else if(Blank_Line_Check(buf)==true) continue;
		 else count++;
		}

	rewind(fp);

	num_xatomtype=count;

	xatom=new Xatom_Def[num_xatomtype];
	if(xatom==NULL) Memory_Allocation_Error();

	// second, read parameters for each atom type

        count=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
		 else if(Blank_Line_Check(buf)==true) continue;
		 else
			{
			 sscanf(buf,"%*d%s%*s%f", 
				xatom[count].type, 
			       &xatom[count].logp);
			 count++;
			}
		}

	fclose(fp);

	if(num_xatomtype!=count) Read_File_Error(filename);

	return;
}

void ForceField::Read_BOND_DEF(char *filename)
{
	FILE *fp;
	int tmp;
	char buf[256];

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

	// first, determine num_bondtype

        tmp=0;

	for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
		 else if(Blank_Line_Check(buf)==true) continue;
		 else tmp++;
		}

	rewind(fp);

	num_bondtype=tmp;

	bond=new Bond_Def[num_bondtype];
	if(bond==NULL) Memory_Allocation_Error();

	// second, read bond parameters

	tmp=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
		 else if(Blank_Line_Check(buf)==true) continue;
		 else 
			{
			 sscanf(buf,"%s%s%s%f", 
				 bond[tmp].atom_1, 
				 bond[tmp].atom_2,
				 bond[tmp].type, 
				&bond[tmp].length);
			 tmp++;
			}
		}

	fclose(fp);

	if(num_bondtype!=tmp) Read_File_Error(filename);

	return;
}

void ForceField::Read_TORSION_DEF(char *filename)
{
	FILE *fp;
	int tmp,mark;
	char buf[256];

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

	// first, determine num_torstype
 
	tmp=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
		 else if(Blank_Line_Check(buf)==true) continue;
                 else tmp++;
                }

	rewind(fp);

	num_torstype=tmp;

	torsion=new Tors_Def[num_torstype];
	if(torsion==NULL) Memory_Allocation_Error();

	// second, read torsion parameters

        tmp=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
		 else if(Blank_Line_Check(buf)==true) continue;
                 else
                        {
                         sscanf(buf,"%s%s%s%s%s%f%d",
                                 torsion[tmp].atom_1, 
				 torsion[tmp].atom_2,
				 torsion[tmp].atom_3, 
				 torsion[tmp].atom_4,
                                 torsion[tmp].type, 
				&torsion[tmp].V, 
				&mark);

			 if(mark>0) {torsion[tmp].n=mark; torsion[tmp].S=1;}
			 else {torsion[tmp].n=-mark; torsion[tmp].S=-1;}
                         tmp++;
                        }
                }

        fclose(fp);

        if(num_torstype!=tmp) Read_File_Error(filename);

	return;
}

void ForceField::Read_PMF_DEF(char *filename)
{
	FILE *fp;
	int i,tmp;
	char buf[256],head[256];

	if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        // first, determine num_restype

        tmp=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
                 else if(Blank_Line_Check(buf)==true) continue;
		 else
			{
			 sscanf(buf,"%s", head);
			 if(strcmp(head,"PAIR")) continue;
			 else tmp++;
			}
                }

        rewind(fp);

        num_pmftype=tmp;

        pmf=new PMF_Def[num_pmftype];
        if(pmf==NULL) Memory_Allocation_Error();

	// second, read the residue templates

	tmp=0;

	for(;;)
		{
		 if(fgets(buf,256,fp)==NULL) break;
		 else if(Blank_Line_Check(buf)==true) continue;
		 else if(buf[0]=='#') continue;
		 
		 sscanf(buf,"%s",head);

		 if(!strcmp(head,"END")) break;
		 else if(!strcmp(head,"PAIR"))
			{
			 sscanf(buf,"%*s%s%s",pmf[tmp].ptype,pmf[tmp].ltype);
			 for(i=0;i<60;i++)
				{
				 fgets(buf,256,fp);
				 sscanf(buf,"%f%f", 
					&pmf[tmp].d[i], &pmf[tmp].p[i]);
				}
			 tmp++;
			}
		}

	fclose(fp);

        if(tmp!=num_pmftype) Read_File_Error(filename);

	return;
}

void ForceField::Read_HBOND_DEF(char *filename)
{
	FILE *fp;
	int i,count;
	char buf[256],head[256];

	if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        // first, determine num_hbtype

        count=0;

	while(fgets(buf,256,fp)!=NULL)
                {
                 if(buf[0]=='#') continue;
                 else if(Blank_Line_Check(buf)==true) continue;
		 else
			{
			 sscanf(buf,"%s", head);
			 if(strcmp(head,"TYPE")) continue;
			 else 
				{
				 count++;
				}
			}
                }

        rewind(fp);

        num_hbtype=count;

        hb_pmf=new HB_Def[num_hbtype];
        if(hb_pmf==NULL) Memory_Allocation_Error();

	// second, read the residue templates

	count=0;

	while(fgets(buf,256,fp)!=NULL)
		{
		 if(Blank_Line_Check(buf)==true) continue;
		 else if(buf[0]=='#') continue;
		 
		 sscanf(buf,"%s",head);

		 if(!strcmp(head,"END")) break;
		 else if(!strcmp(head,"TYPE"))
		{
		 sscanf(buf,"%*s%s", hb_pmf[count].type);

		 fgets(buf,256,fp);
		 sscanf(buf,"%*s%f%*f%*f", &hb_pmf[count].low_cutoff_d);
		 sscanf(buf,"%*s%*f%f%*f", &hb_pmf[count].high_cutoff_d);
		 sscanf(buf,"%*s%*f%*f%f", &hb_pmf[count].step_d);
		 hb_pmf[count].num_bin_d=(int)((hb_pmf[count].high_cutoff_d-
                                          hb_pmf[count].low_cutoff_d)/
					  hb_pmf[count].step_d);

		 fgets(buf,256,fp);
		 sscanf(buf,"%*s%f%*f%*f", &hb_pmf[count].low_cutoff_a1);
		 sscanf(buf,"%*s%*f%f%*f", &hb_pmf[count].high_cutoff_a1);
		 sscanf(buf,"%*s%*f%*f%f", &hb_pmf[count].step_a1);
		 hb_pmf[count].num_bin_a1=(int)((hb_pmf[count].high_cutoff_a1-
                                          hb_pmf[count].low_cutoff_a1)/
					  hb_pmf[count].step_a1);

		 fgets(buf,256,fp);
		 sscanf(buf,"%*s%f%*f%*f", &hb_pmf[count].low_cutoff_a2);
		 sscanf(buf,"%*s%*f%f%*f", &hb_pmf[count].high_cutoff_a2);
		 sscanf(buf,"%*s%*f%*f%f", &hb_pmf[count].step_a2);
		 hb_pmf[count].num_bin_a2=(int)((hb_pmf[count].high_cutoff_a2-
                                          hb_pmf[count].low_cutoff_a1)/
					  hb_pmf[count].step_a2);

		 hb_pmf[count].num_bin_total=hb_pmf[count].num_bin_d*
                                             hb_pmf[count].num_bin_a1*
                                             hb_pmf[count].num_bin_a2;

		 hb_pmf[count].pmf=new float[hb_pmf[count].num_bin_total];
		 if(hb_pmf[count].pmf==NULL) Memory_Allocation_Error();

		 for(i=0;i<hb_pmf[count].num_bin_total;i++)
			{
			 fgets(buf,256,fp);
			 sscanf(buf,"%*f%*f%*f%*s%*f%*f%*s%*f%*f%*s%f",
				&hb_pmf[count].pmf[i]);
			}

		 count++;
		}
		}

	fclose(fp);

        if(count!=num_hbtype) Read_File_Error(filename);

	return;
}

void ForceField::Read_WPMF_DEF(char *filename)
{
	FILE *fp;
	int i,count;
	char buf[256],head[256];

	if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        // first, determine num_hbtype

        count=0;

	while(fgets(buf,256,fp)!=NULL)
                {
                 if(buf[0]=='#') continue;
                 else if(Blank_Line_Check(buf)==true) continue;
		 else
			{
			 sscanf(buf,"%s", head);
			 if(strcmp(head,"PAIR")) continue;
			 else {count++;}
			}
                }

        rewind(fp);

        num_wpmftype=count;

	wpmf=new WPMF_Def[num_wpmftype];
	if(wpmf==NULL) Memory_Allocation_Error();

	// read the bin definition

	fgets(buf,256,fp);

	for(i=0;i<num_wpmftype;i++)
		{
		 sscanf(buf,"%*s%f%f%f", 
			&wpmf[i].low_cutoff,
			&wpmf[i].high_cutoff,
			&wpmf[i].step);
		 wpmf[i].num_bin=(int)floor((wpmf[i].high_cutoff-
                                        wpmf[i].low_cutoff)/wpmf[i].step);

		 wpmf[i].pmf=new float[wpmf[i].num_bin];
		 if(wpmf[i].pmf==NULL) Memory_Allocation_Error();
		}

	// read the PMFs 

	count=0;

	while(fgets(buf,256,fp)!=NULL)
		{
		 if(Blank_Line_Check(buf)==true) continue;
		 else if(buf[0]=='#') continue;
		 
		 sscanf(buf,"%s",head);

		 if(!strcmp(head,"END")) break;
		 else if(!strcmp(head,"PAIR"))
		{
		 sscanf(buf,"%*s%*s%s", wpmf[count].type);

		 for(i=0;i<wpmf[count].num_bin;i++)
			{
			 fgets(buf,256,fp);
			 sscanf(buf,"%*s%*s%*s%*d%*f%*s%*d%*f%*s%*f%f",
				&wpmf[count].pmf[i]);
			}

		 count++;
		}
		}

	fclose(fp);

        if(count!=num_wpmftype) Read_File_Error(filename);

	return;
}

void ForceField::Read_SURFACE_DEF(char *filename)
{
	FILE *fp;
	int i,count;
	char buf[256],head[256];
	Dot tmp_dot;

	if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        // first, determine num_sdot_type;

        count=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
                 else if(Blank_Line_Check(buf)==true) continue;
		 else
			{
			 sscanf(buf,"%s", head);
			 if(strcmp(head,"DOTSET")) continue;
			 else count++;
			}
                }

        rewind(fp);

        num_sdot_type=count;

        sdot=new DotSet[num_sdot_type];
        if(sdot==NULL) Memory_Allocation_Error();

	// now read pre-calculated volume dots 

	count=0;

	for(;;)
		{
		 if(fgets(buf,256,fp)==NULL) break;
		 else if(Blank_Line_Check(buf)==true) continue;
		 else if(buf[0]=='#') continue;
		 
		 sscanf(buf,"%s",head);

		 if(!strcmp(head,"END")) break;
		 else if(!strcmp(head,"DOTSET"))
			{
			 sscanf(buf,"%*s%f%d%f%f", 
				&sdot[count].r,
				&sdot[count].num_dot,
				&sdot[count].unit,
				&sdot[count].total);

			 strcpy(sdot[count].type,"Un");

			 sdot[count].dot.clear();

			 for(i=0;i<sdot[count].num_dot;i++)
				{
				 fgets(buf,256,fp);
				 sscanf(buf,"%f%f%f", 
					&tmp_dot.coor[0],
					&tmp_dot.coor[1],
					&tmp_dot.coor[2]);

				 tmp_dot.valid=1;
				 tmp_dot.unit=sdot[count].unit;
				 strcpy(tmp_dot.type,"Un");

				 sdot[count].dot.push_back(tmp_dot);
				}

			 if(sdot[count].dot.size()!=sdot[count].num_dot)
				Read_File_Error(filename);

			 count++;
			}
		}

	fclose(fp);

        if(count!=num_sdot_type) Read_File_Error(filename);

	return;
}

void ForceField::Read_VOLUME_DEF(char *filename)
{
	FILE *fp;
	int i,count;
	char buf[256],head[256];
	Dot tmp_dot;

	if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        // first, determine num_vdot_type;

        count=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break;
                 else if(buf[0]=='#') continue;
                 else if(Blank_Line_Check(buf)==true) continue;
		 else
			{
			 sscanf(buf,"%s", head);
			 if(strcmp(head,"DOTSET")) continue;
			 else count++;
			}
                }

        rewind(fp);

        num_vdot_type=count;

        vdot=new DotSet[num_vdot_type];
        if(vdot==NULL) Memory_Allocation_Error();

	// now read pre-calculated volume dots 

	count=0;

	for(;;)
		{
		 if(fgets(buf,256,fp)==NULL) break;
		 else if(Blank_Line_Check(buf)==true) continue;
		 else if(buf[0]=='#') continue;
		 
		 sscanf(buf,"%s",head);

		 if(!strcmp(head,"END")) break;
		 else if(!strcmp(head,"DOTSET"))
			{
			 sscanf(buf,"%*s%f%d%f%f", 
				&vdot[count].r,
				&vdot[count].num_dot,
				&vdot[count].unit,
				&vdot[count].total);

			 strcpy(vdot[count].type,"Un");

			 vdot[count].dot.clear();

			 for(i=0;i<vdot[count].num_dot;i++)
				{
				 fgets(buf,256,fp);
				 sscanf(buf,"%f%f%f", 
					&tmp_dot.coor[0],
					&tmp_dot.coor[1],
					&tmp_dot.coor[2]);

				 tmp_dot.valid=1;
				 tmp_dot.unit=vdot[count].unit;
				 strcpy(tmp_dot.type,"Un");

				 vdot[count].dot.push_back(tmp_dot);
				}

			 if(vdot[count].dot.size()!=vdot[count].num_dot)
				Read_File_Error(filename);

			 count++;
			}
		}

	fclose(fp);

        if(count!=num_vdot_type) Read_File_Error(filename);

	return;
}

void ForceField::Show_Contents() const
{
	int i,j;

	if(num_restype>0)
	{
	 printf("\nTotal number of residue templates is: %d\n",num_restype);

         for(i=0;i<num_restype;i++)
                {
                 printf("RESI %3s %d %d\n",
                        residue[i].name,
			residue[i].num_atom,
			residue[i].num_bond);
                 for(j=0;j<residue[i].num_atom;j++)
                        {
                         printf("ATOM %-4s ", residue[i].atom[j].name);
			 printf("%-5s ", residue[i].atom[j].type);
			 printf("%-8s ", residue[i].atom[j].xtype);
			 printf("%5.3f ", residue[i].atom[j].r);
			 printf("%6.3f ", residue[i].atom[j].q);
			 printf("%-2s ", residue[i].atom[j].hb);
			 printf("%6.3f ", residue[i].atom[j].logp);
			 printf("%1d ", residue[i].atom[j].ring);
			 printf("%2s\n", residue[i].atom[j].pmftype);
                        }
		 for(j=0;j<residue[i].num_bond;j++)
			{
			 printf("BOND %-4s ", residue[i].bond[j].atom1);
			 printf("%-5s ", residue[i].bond[j].atom2);
			 printf("%-2s\n", residue[i].bond[j].type);
			}
                }
	}

	if(num_atomtype>0)
	{
        printf("\nTotal number of atom types is %d\n", num_atomtype);

        for(i=0;i<num_atomtype;i++)
                {
                 printf("%-10s ", atom[i].type);  
		 printf("%6.2f ", atom[i].weight);
		 printf("%5.3f ", atom[i].r); 
		 printf("%-2s\n", atom[i].hb);
                }

	printf("\nTotal number of X types is %d\n", num_xatomtype);

        for(i=0;i<num_xatomtype;i++)
                {
                 printf("%-20s ", xatom[i].type);
                 printf("%6.3f\n",xatom[i].logp);
                }
	}

	if(num_bondtype>0)
	{
        printf("\nTotal number of bond types is %d\n", num_bondtype);

        for(i=0;i<num_bondtype;i++)
                {
                 printf("%10s %10s %5s %5.3f\n",
                        bond[i].atom_1,
                        bond[i].atom_2,
                        bond[i].type,
                        bond[i].length);
                }
	}

	if(num_torstype>0)
	{
        printf("\nTotal number of torsion angle types is %d\n", num_bondtype);

        for(i=0;i<num_torstype;i++)
                {
                 printf("%10s %10s %10s %10s %5s %5.3f %2d %2d\n",
                        torsion[i].atom_1,
                        torsion[i].atom_2,
                        torsion[i].atom_3,
                        torsion[i].atom_4,
                        torsion[i].type,
                        torsion[i].V,
                        torsion[i].n,
                        torsion[i].S);
                }
	}

	if(num_pmftype>0)
	{
	printf("\nTotal number of PMF pairs is %d\n", num_pmftype);

        for(i=0;i<num_pmftype;i++)
                {
                 printf("PAIR %s %s\n", pmf[i].ptype, pmf[i].ltype);
                 for(j=0;j<60;j++)
                        {
                         printf("%8.3f %8.3f\n", pmf[i].d[j], pmf[i].p[j]);
                        }
                }
	}

	if(num_sdot_type>0)
	{
	 printf("\nTotal number of surface dot sets is %d\n", num_sdot_type);
	 for(i=0;i<num_sdot_type;i++) sdot[i].Show_Contents();
	}

	if(num_vdot_type>0)
	{
	 printf("\nTotal number of volume dot sets is %d\n", num_vdot_type);
	 for(i=0;i<num_vdot_type;i++) vdot[i].Show_Contents();
	}

	if(num_hbtype>0)
	{
	 int num,d,a1,a2;

	 printf("\nTotal number of HB_PMF profiles is %d\n", num_hbtype);

         for(i=0;i<num_hbtype;i++)
                {
                 printf("TYPE %s\n", hb_pmf[i].type);

                 for(j=0;j<hb_pmf[i].num_bin_total;j++)
                {
		 num=j;
		 d=num/(hb_pmf[i].num_bin_a1*hb_pmf[i].num_bin_a2);
		 num-=(d*hb_pmf[i].num_bin_a1*hb_pmf[i].num_bin_a2);
		 a1=num/hb_pmf[i].num_bin_a2;
		 num-=(a1*hb_pmf[i].num_bin_a2);
		 a2=num;

		 printf("%5.2f  ",hb_pmf[i].low_cutoff_d+d*hb_pmf[i].step_d);
		 printf("%5.1f  ",hb_pmf[i].low_cutoff_a1+a1*hb_pmf[i].step_a1);
		 printf("%5.1f  ",hb_pmf[i].low_cutoff_a2+a2*hb_pmf[i].step_a2);
		 printf("%6.3f\n", hb_pmf[i].pmf[j]);
                }
                }
	}

	if(num_wpmftype>0)
	{
	 printf("\nTotal number of WPMF profiles is %d\n", num_wpmftype);

	 for(i=0;i<num_wpmftype;i++)
		{
		 printf("PAIR  O.w  %-s\n", wpmf[i].type);
		 for(j=0;j<wpmf[i].num_bin;j++)
			{
			 printf("%4.2f -- %4.2f: ",
				wpmf[i].step*j+wpmf[i].low_cutoff,	
				wpmf[i].step*(j+1)+wpmf[i].low_cutoff);
			 printf("%7.3f\n", wpmf[i].pmf[j]);
			}
		}
	}

	return;
}

// ****************************************************************************
// find parameters for protein atoms from the residue template file
// latest update: 11/04/2003
// ****************************************************************************
bool ForceField::Assign_Patom_Parameters(Atom &atm) const
{
        int i,j;

	// First, check the standard residue templates to map the atom 

	for(i=0;i<num_restype;i++)
	{
         if(strcmp(atm.residue,residue[i].name)) continue; 

         for(j=0;j<residue[i].num_atom;j++) // residue template found
                {
		 if(strcmp(residue[i].atom[j].name,atm.name)) continue;

                 strcpy(atm.type,residue[i].atom[j].type);
		 strcpy(atm.xtype,residue[i].atom[j].xtype);
		 strcpy(atm.type2,residue[i].atom[j].pmftype);
                 atm.r=residue[i].atom[j].r;
		 atm.eps=residue[i].atom[j].eps;
                 atm.q=residue[i].atom[j].q;
		 atm.R=residue[i].atom[j].r;
                 strcpy(atm.hb,residue[i].atom[j].hb);
                 atm.logp=residue[i].atom[j].logp;
		 atm.solv=residue[i].atom[j].solv;
	 	 atm.ring=residue[i].atom[j].ring;
                 atm.valid=1; return true;
                }
	}

	// check if they are terminal atoms

	for(i=0;i<num_restype;i++)
	{
         if(strcmp(residue[i].name,"TER")) continue; 

         for(j=0;j<residue[i].num_atom;j++)  // residue template found
                {
		 if(strcmp(residue[i].atom[j].name,atm.name)) continue;

                 strcpy(atm.type,residue[i].atom[j].type);
		 strcpy(atm.xtype,residue[i].atom[j].xtype);
		 strcpy(atm.type2,residue[i].atom[j].pmftype);
                 atm.r=residue[i].atom[j].r;
		 atm.eps=residue[i].atom[j].eps;
                 atm.q=residue[i].atom[j].q;
		 atm.R=residue[i].atom[j].r;
                 strcpy(atm.hb,residue[i].atom[j].hb);
                 atm.logp=residue[i].atom[j].logp;
		 atm.solv=residue[i].atom[j].solv;
	 	 atm.ring=residue[i].atom[j].ring;
                 atm.valid=1; return true;
                }
	}

	// check if they are valid hetero atoms 

	for(i=0;i<num_restype;i++)
        {
         if(strcmp(residue[i].name,"HET")) continue;

         for(j=0;j<residue[i].num_atom;j++)  // residue template found
                {
                 if(strcmp(residue[i].atom[j].name,atm.name)) continue;

                 strcpy(atm.type,residue[i].atom[j].type);
                 strcpy(atm.xtype,residue[i].atom[j].xtype);
                 strcpy(atm.type2,residue[i].atom[j].pmftype);
                 atm.r=residue[i].atom[j].r;
                 atm.eps=residue[i].atom[j].eps;
                 atm.q=residue[i].atom[j].q;
                 atm.R=residue[i].atom[j].r;
                 strcpy(atm.hb,residue[i].atom[j].hb);
                 atm.logp=residue[i].atom[j].logp;
                 atm.solv=residue[i].atom[j].solv;
                 atm.ring=residue[i].atom[j].ring;
                 atm.valid=1; return true;
                }
        }

	printf("Warning: no parameter for patom %s(%s,%s) ... ",
		atm.name, atm.residue, atm.res_id);
	printf("this atom is ignored\n");

	atm.valid=0; return false;
}

bool ForceField::Assign_Atom_Parameters(Atom &atm) const
{
	int i;

	// notice that partial charge is not assigned here
	// in order not to overwrite the original charge on this atom

	for(i=0;i<num_atomtype;i++)
                {
                 if(strcmp(atm.xtype,atom[i].type)) continue;
                 else
                        {
			 atm.weight=atom[i].weight;
			 atm.R=atom[i].r;
			 atm.r=atom[i].r;
			 atm.eps=atom[i].eps;
                         strcpy(atm.hb,atom[i].hb);
			 atm.valid=1; return true;
                        }
                }

	printf("Warning: no parameter for atom %d (%s,%s) ... ",
                atm.id, atm.type, atm.xtype);
        printf("this atom is ignored\n");

	atm.valid=0; return false;
}

bool ForceField::Patom_Connection_Test(Atom atm1, Atom atm2) const
{
	if(atm1.chain!=atm2.chain) return false;
	if(strcmp(atm1.residue,atm2.residue)) return false;
	if(strcmp(atm1.res_id,atm2.res_id)) return false;
	if(atm1.id==atm2.id) return false;

	int i,j;

	for(i=0;i<this->num_restype;i++)
	{
	 if(strcmp(this->residue[i].name,atm1.residue)) continue;
	 else if(strcmp(this->residue[i].name,atm2.residue)) continue;

	 if(this->residue[i].num_bond==0) break;

	 for(j=0;j<this->residue[i].num_bond;j++)
		{
		 if(!strcmp(this->residue[i].bond[j].atom1,atm1.name)&&
                    !strcmp(this->residue[i].bond[j].atom2,atm2.name))
			{
			 return true;
			}
		 else if(!strcmp(this->residue[i].bond[j].atom1,atm2.name)&&
                         !strcmp(this->residue[i].bond[j].atom2,atm1.name))
                        {
                         return true;
                        }
		 else continue;
		}

	 break;
	}

	return false;
}

float ForceField::Get_Atom_LogP(char *type) const
{
	int i;
	bool mark;
	float logp;

	mark=false;

	for(i=0;i<num_xatomtype;i++)
	{
	 if(strcmp(type,xatom[i].type)) continue;
	 else {logp=xatom[i].logp; mark=true; break;}
	}

	if(mark==true) return logp;
	else 	
	{
	 printf("Warning: no LogP parameter for atom type %s ... ", type);
	 printf("zero value assigned\n");
	 return 0.000;
	}
}

float ForceField::Get_Bond_Length(char *a1, char *a2, char *bond_type) const
{
        int i;
	bool mark;
	float length;
	char atom_1[20], atom_2[20];

	strcpy(atom_1,a1); strcpy(atom_2,a2);

        if(!strcmp(atom_1,"H.spc")) strcpy(atom_1,"H");
        if(!strcmp(atom_2,"H.spc")) strcpy(atom_2,"H");

        mark=false;

        for(i=0;i<num_bondtype;i++)
        {
         if(strcmp(bond_type,bond[i].type)) continue;
         else if(!strcmp(atom_1,bond[i].atom_1)&&!strcmp(atom_2,bond[i].atom_2))
                {length=bond[i].length; mark=true; break;}
         else if(!strcmp(atom_2,bond[i].atom_1)&&!strcmp(atom_1,bond[i].atom_2))
                {length=bond[i].length; mark=true; break;}
         else continue;
        }

	if(mark==true) return length;
	else
                {
                 printf("Warning: no bond parameter for ");
                 printf("%s-%s-%s ... ", atom_1, bond_type, atom_2);
		 printf("default value assigned\n");
		 return 1.54;  // default value for C.3-C.3
                }
}

bool ForceField::Assign_Torsion_Parameters(Torsion &tors) const
{
	int i,mark,tmp;

	if(!strcmp(tors.atom_1.type,"H.spc")) strcpy(tors.atom_1.type,"H");
	if(!strcmp(tors.atom_2.type,"H.spc")) strcpy(tors.atom_2.type,"H");
	if(!strcmp(tors.atom_3.type,"H.spc")) strcpy(tors.atom_3.type,"H");
	if(!strcmp(tors.atom_4.type,"H.spc")) strcpy(tors.atom_4.type,"H");
	
	if(!strcmp(tors.atom_1.type,"N.4")) strcpy(tors.atom_1.type,"N.3");
	if(!strcmp(tors.atom_2.type,"N.4")) strcpy(tors.atom_2.type,"N.3");
	if(!strcmp(tors.atom_3.type,"N.4")) strcpy(tors.atom_3.type,"N.3");
	if(!strcmp(tors.atom_4.type,"N.4")) strcpy(tors.atom_4.type,"N.3");

	if(!strcmp(tors.atom_1.type,"O.co2")) strcpy(tors.atom_1.type,"O.2");
	if(!strcmp(tors.atom_2.type,"O.co2")) strcpy(tors.atom_2.type,"O.2");
	if(!strcmp(tors.atom_3.type,"O.co2")) strcpy(tors.atom_3.type,"O.2");
	if(!strcmp(tors.atom_4.type,"O.co2")) strcpy(tors.atom_4.type,"O.2");

	mark=false;

	for(i=0;i<num_torstype;i++)
		{
		 if(strcmp(tors.type,torsion[i].type)) continue;
		 else if(!strcmp(tors.atom_2.type,torsion[i].atom_2)&&
			 !strcmp(tors.atom_3.type,torsion[i].atom_3))
			{
			 // find the most suitable torsion parameter
			 tmp=2;
			 if(!strcmp(tors.atom_1.type,torsion[i].atom_1)) tmp++;
			 if(!strcmp(tors.atom_4.type,torsion[i].atom_4)) tmp++;
			 if(tmp>mark)
				{
				 mark=tmp;
				 tors.V=torsion[i].V;
				 tors.n=torsion[i].n;
				 tors.S=torsion[i].S;
				 if(mark==4) break;
				}
			}
		 else if(!strcmp(tors.atom_3.type,torsion[i].atom_2)&&
                         !strcmp(tors.atom_2.type,torsion[i].atom_3))
                        {
			 // find the most suitable torsion parameter
                         tmp=2;
                         if(!strcmp(tors.atom_4.type,torsion[i].atom_1)) tmp++;
                         if(!strcmp(tors.atom_1.type,torsion[i].atom_4)) tmp++;
                         if(tmp>mark)
                                {
                                 mark=tmp;
                                 tors.V=torsion[i].V;
                                 tors.n=torsion[i].n;
                                 tors.S=torsion[i].S;
				 if(mark==4) break;
                                }

                        }
		 else continue;
		}

	if(mark==false)
		{
                 printf("Warning: no torsion parameter for ");
		 printf("%s-%s-%s-%s-%s ... ", 
			 tors.atom_1.type, tors.atom_2.type, tors.type,
			 tors.atom_3.type, tors.atom_4.type);
		 printf("default value assigned.\n");
               	 tors.V=0.5; tors.n=3; tors.S=1;  // C.3-C.3-C.3-C.3
		 return false;
		}
	else return true;
}

float ForceField::Calculate_Torsion_Energy(const Torsion &tors) const
{
	float e,tmp1,tmp2;

	// TORSION_E = 1/2*V*[1+S*cos(n*angle)] 

	tmp1=0.50*tors.V;
	tmp2=tors.angle/180.0*(PI)*tors.n;

	e=tmp1*(1+tors.S*cos((double)tmp2));
	
	return e;
}

float ForceField::Calculate_VDW_Energy(const Atom &a1, const Atom &a2, int mark_1_4) const
{
	float d0,d,eps,e;
	float tmp1,tmp2,tmp3,tmp4;

	// VDW_E = SQRT(e1*e2)*[(d0/d)^12-2.0*(d0/d)^6] 

	e=0.000;

	d=Distance(a1.coor,a2.coor);
	if(d>8.00) return 0.000;

	// I do not use DIST_CUTOFF to replace 8.00 here in order to
	// keep this function independent.

	d0=a1.r+a2.r;

	eps=sqrt((double)(a1.eps*a2.eps));

	tmp1=d0/d; tmp2=tmp1*tmp1*tmp1;
	tmp3=tmp2*tmp2; tmp4=tmp3*tmp3;

	e=eps*(tmp4-2.0*tmp3);

	if(mark_1_4) e*=(0.50);		// VDW_1_4_FACTOR

	if(e>1000.0) e=1000.0;	

	return e;
}

float ForceField::Get_PMF(char *ptype, char *ltype, float d) const
{
	int i,tmp1,tmp2,tmp3;
	float cutoff,sum;

	// set the dual-cutoff

	if((ptype[0]=='C'||ptype[0]=='c')&&
	   ((ltype[0]=='C'||ltype[0]=='c')&&strcmp(ltype,"Cl"))) cutoff=6.00;
	else cutoff=9.00;

	if(d>cutoff) return 0.000;

	// now calculate the PMF score

	sum=0.000;

	for(i=0;i<num_pmftype;i++)
		{
		 if(strcmp(ptype,pmf[i].ptype)) continue;
		 else if(strcmp(ltype,pmf[i].ltype)) continue;
		 else
			{
			 tmp2=(int)(d/0.200); tmp1=tmp2-1; tmp3=tmp2+1; 

			 if(tmp1<0) tmp1=0;
                         else if(tmp1>59) tmp1=59;

			 if(tmp2<0) tmp2=0;
                         else if(tmp2>59) tmp2=59;

			 if(tmp3<0) tmp3=0;
			 else if(tmp3>59) tmp3=59;

			 // no smoothing algorithm
			 // sum+=pmf[i].p[tmp2];

			 // smoothing algorithm

			 sum+=(0.2*pmf[i].p[tmp1]);
			 sum+=(0.6*pmf[i].p[tmp2]);
			 sum+=(0.2*pmf[i].p[tmp3]);

			 break;
			}
		}

	return sum;
}

// **************************************************************************
// latest update: 11/06/2003
// **************************************************************************
float ForceField::Get_WPMF(char *type, float d) const
{
 int i,bin;

 for(i=0;i<num_wpmftype;i++)
	{
	 if(strcmp(type,wpmf[i].type)) continue;

	 if(d<wpmf[i].low_cutoff) return 0.000;
	 else if(d>wpmf[i].high_cutoff) return 0.000;
	 else
		{
	 	 bin=(int)floor((d-wpmf[i].low_cutoff)/wpmf[i].step);
	 	 return wpmf[i].pmf[bin];
		}
	}

 return 0.000;  // type not found
}

float ForceField::Get_HB_PMF(int type, float d, float a1, float a2) const
{
	int i,id,tmp1,tmp2,tmp3,tmp;
	char hbtype[5];
	float pmf;

	if(type==1) strcpy(hbtype,"D1");
	else if(type==2) strcpy(hbtype,"D2");
	else return 0.000;

	id=0;

	for(i=0;i<num_hbtype;i++)
		{
		 if(strcmp(hb_pmf[i].type,hbtype)) continue;
		 else {id=i; break;} 
		}

	if(d<hb_pmf[id].low_cutoff_d) return 0.000;
	if(d>hb_pmf[id].high_cutoff_d) return 0.000;
	if(a1<hb_pmf[id].low_cutoff_a1) return 0.000;
	if(a1>hb_pmf[id].high_cutoff_a1) return 0.000;
	if(a2<hb_pmf[id].low_cutoff_a2) return 0.000;
	if(a2>hb_pmf[id].high_cutoff_a2) return 0.000;

	tmp1=(int)floor((d-hb_pmf[id].low_cutoff_d)/hb_pmf[id].step_d);
	tmp2=(int)floor((a1-hb_pmf[id].low_cutoff_a1)/hb_pmf[id].step_a1);
	tmp3=(int)floor((a2-hb_pmf[id].low_cutoff_a2)/hb_pmf[id].step_a2);

	tmp=tmp1*hb_pmf[id].num_bin_a1*hb_pmf[id].num_bin_a2;
	tmp+=tmp2*hb_pmf[id].num_bin_a2;
	tmp+=tmp3;

	pmf=hb_pmf[id].pmf[tmp];

	if(pmf<0.000) return 0.000;
	else return pmf;
}

DotSet ForceField::Get_Surface_Dot(Atom atom, float r) const
{
	int i;
	float R;
	DotSet tmp_set;

	R=atom.R+r;

	tmp_set=Get_Surface_Dot(R,atom.coor[0],atom.coor[1],atom.coor[2]);

	strcpy(tmp_set.type, atom.type); 

	for(i=0;i<tmp_set.num_dot;i++) 
		{
		 tmp_set.dot[i].valid=atom.id;
		 strcpy(tmp_set.dot[i].type, atom.type);
		}

	return tmp_set;
}

DotSet ForceField::Get_Surface_Dot(float R, float x, float y, float z) const
{
	int i,num;
	bool mark;
	float tmp,theta,phi,theta_step,phi_step;
	float r,d,total;
	DotSet tmp_set;
	Dot tmp_dot;

	// check the pre-calculated surface dot sets

        mark=false;

        for(i=0;i<num_sdot_type;i++)
                {
                 if(fabs(R-sdot[i].r)>0.025) continue;
                 else {tmp_set=sdot[i]; mark=true; break;}
                }

	if(mark==true) 
	{
	 for(i=0;i<tmp_set.num_dot;i++)
		{
		 tmp_set.dot[i].coor[0]+=x;
		 tmp_set.dot[i].coor[1]+=y;
		 tmp_set.dot[i].coor[2]+=z;
		}
	 return tmp_set;
	}

        // if it is not pre-calculated, calculate it now 

	total=(4.000*PI*R*R);

	// d=0.500;		// spacing between two dots

	d=sqrt(total/300);
	if(d<0.500) d=0.500;

	tmp=(int)(PI*R/d+0.500); theta_step=PI/tmp;

	num=0; tmp_set.dot.clear();

	for(theta=0.00;theta<PI;theta+=theta_step)
		{
		 r=R*sin(theta);
		 tmp=(int)(2*PI*r/d+0.500); phi_step=2*PI/tmp;

		 for(phi=0.00;phi<(2*PI);phi+=phi_step)
                        {
                         tmp_dot.coor[0]=r*cos(phi);
                         tmp_dot.coor[1]=r*sin(phi);
                         tmp_dot.coor[2]=R*cos(theta);
			 tmp_dot.valid=1;
			 strcpy(tmp_dot.type,"Un");
			 tmp_set.dot.push_back(tmp_dot);
                         num++;
                        }
		}

	strcpy(tmp_set.type, "Un"); tmp_set.r=R; 
        tmp_set.num_dot=num; tmp_set.unit=total/num;
	tmp_set.total=total;

	for(i=0;i<tmp_set.num_dot;i++) 
		{
		 tmp_set.dot[i].unit=tmp_set.unit;
		 tmp_set.dot[i].coor[0]+=x;
		 tmp_set.dot[i].coor[1]+=y;
		 tmp_set.dot[i].coor[2]+=z;
		}

	return tmp_set;
}

DotSet ForceField::Get_Volume_Dot(float R, float x, float y, float z) const
{
	int i,num;
	bool mark;
	float r,total,spacing;
	float x_min,x_max,y_min,y_max,z_min,z_max;
	float probe[3];
	DotSet tmp_set;
	Dot tmp_dot;

	// check the pre-calculated surface dot sets

        mark=false;

        for(i=0;i<num_vdot_type;i++)
                {
                 if(fabs(R-vdot[i].r)>0.025) continue;
                 else {tmp_set=vdot[i]; mark=true; break;}
                }

	if(mark==true) 
	{
	 for(i=0;i<tmp_set.num_dot;i++)
		{
		 tmp_set.dot[i].coor[0]+=x;
		 tmp_set.dot[i].coor[1]+=y;
		 tmp_set.dot[i].coor[2]+=z;
		}
	 return tmp_set;
	}

        // if it is not pre-calculated, calculate it now 

	total=4.00*PI*R*R*R/3.00;
	spacing=powf(8.0*R*R*R/200,1/3.0);
	if(spacing<0.25) spacing=0.25;

	x_max=ceil(2*R/spacing)*spacing/2.0; x_min=-x_max;
	y_max=ceil(2*R/spacing)*spacing/2.0; y_min=-y_max;
	z_max=ceil(2*R/spacing)*spacing/2.0; z_min=-z_max;

	num=0; tmp_set.dot.clear();
	r=R;
	// r=R-spacing/4.0;

	for(probe[0]=x_min; probe[0]<(x_max+0.01); probe[0]+=spacing)
	for(probe[1]=y_min; probe[1]<(y_max+0.01); probe[1]+=spacing)
	for(probe[2]=z_min; probe[2]<(z_max+0.01); probe[2]+=spacing)
		{
		 if(sqrt(probe[0]*probe[0]+
                         probe[1]*probe[1]+
                         probe[2]*probe[2])>r) continue;
                 else
                        {
                         tmp_dot.coor[0]=probe[0];
                         tmp_dot.coor[1]=probe[1];
                         tmp_dot.coor[2]=probe[2];
			 tmp_dot.valid=1;
			 strcpy(tmp_dot.type,"Un");
			 tmp_set.dot.push_back(tmp_dot);
                         num++;
                        }
		}

	strcpy(tmp_set.type, "Un"); tmp_set.r=R; 
        tmp_set.num_dot=num; tmp_set.unit=total/num;
	tmp_set.total=total;

	for(i=0;i<tmp_set.num_dot;i++) 
		{
		 tmp_set.dot[i].unit=tmp_set.unit;
		 tmp_set.dot[i].coor[0]+=x;
		 tmp_set.dot[i].coor[1]+=y;
		 tmp_set.dot[i].coor[2]+=z;
		}

	return tmp_set;
}

