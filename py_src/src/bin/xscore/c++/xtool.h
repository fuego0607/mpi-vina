# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <string>
# include <math.h>
# include <ctype.h>
# include <time.h>
# include <vector.h>
# define TRUE 1
# define FALSE 0
# define PI 3.1416
# define TINY 1.0e-6
# define LARGE 1.0e+6
# define DIST_CUTOFF 8.00
# define WATER_R 1.40
# define POCKET_DEPTH 4.00
# define LAYER_DEPTH 3.00
# define MAX_ATOM_NEIB 6 
# define MAX_BOND_NEIB 10 
# define SWAP(a,b) temp=(a); (a)=(b); (b)=temp
# define MAX(a,b) (a)>=(b) ? (a):(b) 
# define MIN(a,b) (a)<=(b) ? (a):(b)
# define SQR(a) (a)==0.0 ? 0.0:(a)*(a)
# define SHIFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d)
# define SIGN(a,b) ((b)>=0.0 ? fabs(a):-fabs(a))

// ************************* functions ********************************

int Xtool_Functions(int argc, char *argv[]);	// main.c
	int Xtool_Score();			// score.c
	int Xtool_LogP();			// logp.c
int Xtool_Utilities(int argc, char *argv[]);	// main.c
	int Xtool_Score_Shortcut();		// score.c
	int Xtool_LogP_Shortcut();		// logp.c
	int Xtool_Fix_PDB();			// utilities.c
	int Xtool_Fix_Mol2();			// utilities.c

// ************************* input.c ********************************

class Input 
{
 public:
	char function[256];

	char input_file[256];
	char output_file[256];
	char log_file[256];

	char receptor_file[256];
	char cofactor_file[256];
	char reference_file[256];
	char ligand_file[256];

	char parameter_dir[256];

	int num_method;
	char apply_hpscore[10];
	char apply_hmscore[10];
	char apply_hsscore[10];
	char apply_pmfscore[10];
	char show_abs[10];

	float hpscore_cvdw;
	float hpscore_chb;
	float hpscore_chp;
	float hpscore_crt;
	float hpscore_c0;

	float hmscore_cvdw;
	float hmscore_chb;
	float hmscore_chm;
	float hmscore_crt;
	float hmscore_c0;

	float hsscore_cvdw;
	float hsscore_chb;
	float hsscore_chs;
	float hsscore_crt;
	float hsscore_c0;

	Input(); ~Input();

	virtual void Read_Inputs(char *filename);
	virtual void Show_Contents() const;

	void Missing_Parameter_Error(const char *name) const;
	void Invalid_Parameter_Error(const char *name) const;
	void Check_Directory(char *dirname, int flag=FALSE);
};

class SCORE_Input: public Input
{
 public:
	int num_hits;
	char hits_dir[256];
	
	char apply_chemical_rules[10];
	float max_weight;
       	float min_weight;
       	float max_logp;
       	float min_logp;
       	int max_hb_atom;
       	int min_hb_atom;

	SCORE_Input(); ~SCORE_Input();

	void Read_Inputs(char *filename);
	void Show_Contents() const;
}; 

class LOGP_Input: public Input
{
 public:
	char calculate_logp[10];
	char calculate_mw[10];
	char count_hb_atom[10];
	char count_rotor[10];

        LOGP_Input(); ~LOGP_Input();

        void Read_Inputs(char *filename);
        void Show_Contents() const;
};

class POCKET_Input: public Input
{
 public:
	char pocket_txt_file[256];
	char grid_txt_file[256];	
	char pocket_pdb_file[256];
	char box_pdb_file[256];
	char shape_pdb_file[256];
	char grid_pdb_file[256];
	char pharmacophore_pdb_file[256];
	float probe_radius;
	float shell_depth;
	float min_feature_distance;

        POCKET_Input(); ~POCKET_Input();

        void Read_Inputs(char *filename);
        void Show_Contents() const;
};

class AUTODOCK_Input: public Input
{
 public:
	char autodock_pdbq_file[256];
	char autodock_dlg_file[256];
	char autodock_job_name[256];
	int autodock_num_conformation;
	float autodock_rmsd_tolerance;
	char autodock_mdb_location[256];
	char autodock_rank_conformation[256];
	char calculate_rmsd_distribution[256];

        AUTODOCK_Input(); ~AUTODOCK_Input();

        void Read_Inputs(char *filename);
        void Show_Contents() const;
};

class FLEXX_Input: public Input
{
 public:
	char flexx_mdb_location[256];
	char flexx_job_name[256];
	int flexx_num_conformation;
	char flexx_energy_table[256];
	float flexx_rmsd_tolerance;
	char flexx_rank_conformation[256];

        FLEXX_Input(); ~FLEXX_Input();

        void Read_Inputs(char *filename);
        void Show_Contents() const;
};

class DOCK_Input: public Input
{
 public:
	char work_directory[256];
	int num_entry;
	char pdb_entry[1000][5];

	DOCK_Input(); ~DOCK_Input();

	void Read_Inputs(char *filename);
	void Show_Contents() const;
};

class SURFACE_Input: public Input
{
 public:
	float probe_radius;
	float dot_density;

        SURFACE_Input(); ~SURFACE_Input();

        void Read_Inputs(char *filename);
        void Show_Contents() const;
};

class DATABASE_Input: public Input
{
 public:
	char output_00[256];
	char output_01[256];
	char output_02[256];
	char output_03[256];
	char output_04[256];
	char output_05[256];
	char output_06[256];
	char output_07[256];
	char output_08[256];
	char output_09[256];
	char output_10[256];
	char output_error[256];

        DATABASE_Input(); ~DATABASE_Input();

        void Read_Inputs(char *filename);
        void Show_Contents() const;
};

// ************************* basic.c ********************************

class Atom	// full version
{
 public:
        int id;  	        // atom id 
        short int valid; 	// valid indicator
        char name[10];          // atom name
        char type[10];          // basic atom type
        char xtype[10];         // xtool atom type
	char type2[10];         // another atom type when necessary
	char residue[10];       // residue name
	char res_id[10];       	// residue id: must be a string for PDB
	char chain;             // chain label
        float coor[3];          // coordinates
        float root[3];          // HB root's coordinates 
        float weight;       	// atomic weight
        float r;                // vdw radius
        float eps;              // vdw epsilon value
	float q;		// partial atomic charge
	float R;		// X radius
        float logp;             // atomic hydrophobic scale
	float solv;		// atomic solvation parameter
        char hb[3];             // HB property
	float occupancy;	// occupancy probability, for protein atoms
	float bfactor;		// B-factor, for protein atoms
        float score;            // atomic binding score
        short int ring;         // ring indicator: 1=normal; 2=aromatic
	short int origin;	// atom origin indicator: 1=ligand; 2=protein 
        short int part;         // component indicator: 
			        // for protein atoms: 1=ATOM; 2=HETATM
        short int num_neib;     // number of neighboring atoms
        short int num_nonh;     // number of non-H neighboring atoms
        int neib[MAX_ATOM_NEIB];  // ID of neighboring atoms
	int bond[MAX_ATOM_NEIB];  // ID of neighboring bonds
	short int temp;		// for misc uses

        Atom(); ~Atom();              

	void Clear();
        void Show_Contents() const;

	int Get_Donor_Type() const;
	int Get_Acceptor_Type() const;
};

class Water
{
 public:
        short int id;      	// atom id
        short int valid;        // valid indicator
        char name[10];          // atom name
        char type[10];          // atom type
        char xtype[10];         // another atom type
        char residue[10];       // residue name
        float coor[3];          // coordinates
        float r;                // vdw radius
        float eps;              // vdw epsilon value
        float q;                // partial atomic charge
        float logp;             // atomic hydrophobic scale
        char hb[3];             // HB property
	float depth;		// buried depth
	float score;		// score value

        Water(); ~Water();      

	void Clear();
        void Show_Contents() const;
};

class Bond	// full version
{
 public:
        short int id;
        short int valid;        // valid indicator
        short int atom_1;       // ID of the atom_1
        short int atom_2;       // ID of the atom_2
        char type[3];           // bond type
        short int part;         // ID of the component
        short int ring;         // ring indicator
        float length;           // bond length
        short int num_neib;     // number of neighboring bonds
        short int neib[MAX_BOND_NEIB];     // ID of neighboring bonds

        Bond(); ~Bond(); 

	void Clear();
        void Show_Contents() const;
};

class Torsion	// full version
{
 public:
        Atom atom_1;
        Atom atom_2;
        Atom atom_3;
        Atom atom_4;
        char type[3];           // type of the torsion between 2-3
        short int angle;        // torsion angle, in degree
        float V;                // potential barrier
        short int n;            // periodicity
        short int S;            // sign
        float e;                // torsion energy

        Torsion(); ~Torsion();    

	void Clear();
        void Show_Contents() const;
};

class Group	// full version
{
 public:
        short int valid;
        short int num_neib;    	// number of neighboring atoms
        short int num_nonh;    	// number of neighboring non-h atoms
        short int num_h;       	// number of neighboring hydrogen atoms
        short int num_hetero;  	// number of neighboring heteroatoms
        short int num_pi;      	// number of neighboring pi atoms
        short int num_car;     	// number of neighboring C.ar
        short int num_nar;     	// number of neighboring N.ar
	short int num_db;	// number of double bonds
	short int num_tb;	// number of triple bonds
        short int db_type;      // double bond type
	short int amide;	// indicator for amide group

        Atom center;            // center atom of the group
        Atom neib[MAX_ATOM_NEIB];   // neighboring atoms
        Bond bond[MAX_ATOM_NEIB];   // bonds

        Group(); ~Group(); 

	void Clear();
        void Show_Contents() const;
};

class HBond
{
 public:
	int valid;
	int type;	// 1: donor=latom, acceptor=patom
			// 2: donor=patom, acceptor=latom
			// 3: donor=metal, acceptor=latom/patom
			// 4: donor=latom, acceptor=latom 
			// 5: donor=patom, acceptor=patom
			// 0: invalid, no H-bond
	int d_type;	// donor type: 1=straight; 2=angled
	int a_type;	// acceptor type: 1=straight; 2=angled
	bool sb;	// flag for neutral HB or SB

	int latom;	// id of the ligand atom
	int patom;	// id of the protein atom

	Atom H,D,A;

	float d;	// D-A distance
	float a0;	// D-H-A angle
	float a1;	// DR-D-A angle
	float a2;	// D-A-AR angle

	float score;	// strength of this HBond

	HBond(); ~HBond();

	void Clear();
	void Show_Contents() const;

	int Value_HBond();
		void Determine_DA_Type();
	int Value_HBond_2();
	int Value_SBond();
};

class Ring
{
 public:
	int valid;	// 0 = invalid; 1 = valid
	int type;	// 1 = normal; 2 = aromatic

	int num_member;
	vector <int> atom_id;   // atom id in this ring
	vector <int> bond_id;   // bond id in this ring

	float centroid[3];

	Ring(); ~Ring();

	Ring(const Ring &original);
        Ring& operator = (const Ring &original);

	void Clear();
	void Show_Contents() const;
};

class Dot
{
 public: 
	int valid;		// status indicator
	char type[10];		// type
	float coor[3];		// coordinates
	float unit;		// contribution, in either A^2 or A^3
	float score;		// score on this dot

	Dot(); ~Dot();

	void Clear();
};

class DotSet 
{
 public:
	int num_dot;
	vector <Dot> dot;
	float r;
	char type[10];
	float unit;		// default contribution of each dot to total
	float total;		// total volume or surface

	DotSet(); ~DotSet();

	DotSet(const DotSet &original);
        DotSet& operator = (const DotSet &original);

	void Clear();
	void Show_Contents() const;
	void Show_Dots(char *filename, char *header, char *show="unit") const;
};

class Residue
{
 public: 
	int valid;
	char name[10];		// for PDB files it is a 3-letter string
	char chain;
	char id[10];

	int num_atom;
	vector <Atom> atom;

	vector <Dot> vol_dot;
	vector <Dot> sur_dot;

	Residue(); ~Residue();

	Residue(const Residue &original);
        Residue& operator = (const Residue &original);

	void Clear();
	void Show_Contents() const;

	int Get_Contents_From_Protein(const char *name, 
                                      const char *id, const char chain);
	void Generate_Surface_Dots(float probe_r=0.00);
};

class Chain
{
 public:
	int valid;
	char label;

	int length;
	vector <Residue> residue;

	Chain(); ~Chain();

	Chain(const Chain &original);
        Chain& operator = (const Chain &original);

	void Clear();
	void Show_Contents() const;
};

// ************************* forcefield.c ********************************

class ForceField	// full version
{
 private:
	int num_restype;
        struct Residue_Def
		{
                 char name[10];
                 int num_atom, num_bond;
                 struct
                        {
                         char name[10];
                         char type[10];
			 char xtype[10];
                         float r;
			 float eps;
                         float q;
                         char hb[3];
                         float logp;
			 float solv;
			 int ring;
			 char pmftype[3]; 
                        } atom[50];
		 struct
			{
			 char atom1[10];
			 char atom2[10];
			 char type[3];
			} bond[50];
                }; 
	Residue_Def *residue;

        int num_atomtype;
        struct Atom_Def
                {
                 char type[10];
                 float weight;
                 float r;
                 float eps;
		 float q;
                 char hb[3];
                };
        Atom_Def *atom;

	int num_xatomtype;
	struct Xatom_Def
		{
		 char type[20];
		 float logp;
		};
	Xatom_Def *xatom;

        int num_bondtype;
        struct Bond_Def
                {
                 char atom_1[10];
                 char atom_2[10];
                 char type[3];
                 float length;
                };
        Bond_Def *bond;

        int num_torstype;
        struct Tors_Def
                {
                 char atom_1[10];
                 char atom_2[10];
                 char atom_3[10];
                 char atom_4[10];
                 char type[3];
                 float V;               // twisting force constant
                 int n;                 // periodicity
                 int S;                 // sign of torsion angle type
                };
        Tors_Def *torsion;

	int num_pmftype;
        struct PMF_Def
                {
                 char ltype[3];
                 char ptype[3];
                 float d[60];		// distance
                 float p[60];		// PMF potential
                };
        PMF_Def *pmf;

	int num_sdot_type;		// pre-calculated surface dots
	DotSet *sdot;

	int num_vdot_type;		// pre-calculated volume dots
	DotSet *vdot;

	int num_hbtype;
	struct HB_Def
		{
		 char type[10];
		 float low_cutoff_d, high_cutoff_d, step_d;
		 float low_cutoff_a1, high_cutoff_a1, step_a1;
		 float low_cutoff_a2, high_cutoff_a2, step_a2;
		 int num_bin_d, num_bin_a1, num_bin_a2, num_bin_total;
		 float *pmf;
		};
	HB_Def *hb_pmf;

	int num_wpmftype;
	struct WPMF_Def
		{
		 char type[10];
		 float low_cutoff, high_cutoff, step;
		 int num_bin;
		 float *pmf;
		};
	WPMF_Def *wpmf;

 public:
        ForceField(char *dirname); ~ForceField();                 

	void Read_RESIDUE_DEF(char *filename);
        void Read_ATOM_DEF(char *filename);
        void Read_XATOM_DEF(char *filename);
        void Read_BOND_DEF(char *filename);
        void Read_TORSION_DEF(char *filename);
	void Read_PMF_DEF(char *filename);
	void Read_SURFACE_DEF(char *filename);
	void Read_VOLUME_DEF(char *filename);
	void Read_HBOND_DEF(char *filename);
	void Read_WPMF_DEF(char *filename);
        void Show_Contents() const;

	bool Assign_Patom_Parameters(Atom &atm) const;
        bool Assign_Atom_Parameters(Atom &atm) const;
        bool Assign_Torsion_Parameters(Torsion &tors) const;

	bool Patom_Connection_Test(Atom atm1, Atom atm2) const;

	float Get_Atom_LogP(char *type) const;
        float Get_Bond_Length(char *a1,char *a2,char *bond_type) const;
        float Calculate_Torsion_Energy(const Torsion &tors) const;
        float Calculate_VDW_Energy(const Atom &a1, const Atom &a2,
				     int mark_1_4=FALSE) const;
	float Get_PMF(char *ptype, char *ltype, float d) const;
	float Get_HB_PMF(int type, float d, float a1, float a2) const;
	float Get_WPMF(char *type, float d) const;

	DotSet Get_Surface_Dot(float R, 
                               float x=0.000, 
                               float y=0.000, 
                               float z=0.000) const;	
	DotSet Get_Surface_Dot(Atom atom, float r=0.000) const;
	DotSet Get_Volume_Dot(float R, 
                              float x=0.000, 
                              float y=0.000, 
                              float z=0.000) const;
};

// *************************** molecule.c ********************************** 

class Molecule
{
 public:
	bool xtool_format;
	int id;
        int valid;
        char name[256];
        char formula[256];
        float weight;
	int num_hb_atom;
	int num_rotor;
        float logp;
	float surface,nsur,psur;
	float volume;

        int num_atom;
	Atom *atom;

        int num_bond;
        Bond *bond;

        int num_subst;
	int num_feature;
	int num_set;
	char mol_type[256];
	char charge_type[256]; 

	int num_ring;
	vector <Ring> ring;

	vector <Dot> vol_dot;
	vector <Dot> sur_dot;

        Molecule(); ~Molecule();
        Molecule(int max_atom_num, int max_bond_num);

        Molecule(const Molecule &original);
        Molecule& operator = (const Molecule &original);

        void Clear();
	void Show_Contents() const;
		void Show_Atoms() const;
		void Show_Bonds() const;
		void Show_Rings() const;

        int Read_From_Mol2(char *filename,int position=1);
	int Read_From_Mol2(FILE *fp,int position=1);
        int Write_Out_Mol2(char *filename, char *mode="w") const;
	int Write_Out_Mol2(FILE *fp) const;

        int Value_Atom();
            	void Detect_Connections();
            	bool Check_Atom_Type();	
            	void Detect_Rings();
                	void Reset_Choices(int id, int choice[]) const;
		 	void Clean_Choices(int tried_choice, 
                                           int choice[]) const;
		 	int Get_Choices(int choice[], int wrong_end, 
                                        int current_path_length, 
                                        int current_path[]) const; 
		 	int Look_For_A_Ring(int bond_id, 
                                            int atom_path[],
                                            int bond_path[], 
                                            int required_size=0) const;
	    	void Detect_Aromatic_Rings();
			int Aromatic_Ring_Check_5(int atom_path[],
                                                  int bond_path[]) const;
			int Aromatic_Ring_Check_6(int atom_path[],
						  int bond_path[]) const;
        	void Get_Formula(char *string) const;
        	float Get_Weight() const;
		int Get_XTOOL_Type(int atom_id, char *type) const;

	Group Find_A_Group(int atom_id) const;

	int Get_Num_HB_Atom() const;
	int Get_Num_Heavy_Atom() const;
	float Count_Rotor();

	void Assign_Apparent_Charge();	
	void Assign_MMFF94_Charge();

	void Generate_Volume_Dots(float spacing, float probe_r=0.00);
	float Calculate_Volume() const;
	void Show_Volume_Dots(char *filename, char *show="unit") const; 

	void Generate_Surface_Dots(float probe_r=0.00);
	float Calculate_Surface() const;
	void Show_Surface_Dots(char *filename, char *show="unit") const;

        int Get_Atom_Hybridizing_Type(char *type) const;
        int Inner_Collision_Check() const;
	int Connection_1_2_Check(int id1, int id2) const;
        int Connection_1_3_Check(int id1, int id2) const;
        int Connection_1_4_Check(int id1, int id2) const;
        int Connection_1_5_Check(int id1, int id2) const;
	int Connection_1_6_Check(int id1, int id2) const;
        int Two_Bonds_Connection_Check(Bond b1, Bond b2) const;
	int Connection_In_Ring_Check(int id1, int id2) const;
	int Judge_Terminal_Atom(Atom &atm, int partner) const;

        void Translate(float move[3]);
        void Rotate(float angle,float axis[3],float origin[3]);

	// logp.c

	float Calculate_LogP();
		int Get_XLOGP_Type(int atom_id, char *type) const;
        	Group Find_X_Group(int atom_id) const;
        	float Count_Hydrophobic_Carbon(int flag);
        	float Count_Internal_HBond(int flag);
        	float Count_Halogen_1_3_Pair(int flag);
        	float Count_Nar_1_4_Pair(int flag);
        	float Count_O3_1_4_Pair(int flag);
        	float Count_Acceptor_1_5_Pair(int flag);
		// float Count_Remote_Donor_Pair(int flag);
        	float Count_Amino_Acid(int flag);
        	float Count_Salicylic_Acid(int flag);
		float Count_Sulfonic_Acid(int flag);
        int Hydrophobic_Neighbor_Check(int id) const;
        int Adjacent_Ring_Check(int id) const;
        int Adjacent_Aromatic_Check(int id) const;
	void Write_Out_Log(FILE *fp) const;
};

// ************************ ligand.c *****************************

class Ligand: public Molecule
{
 public:
        float bind_score,chem_score;
	float vdw,sb,hb,hp,hm,hs,ar,rt,pmf,uhb,bnsur,bpsur;
	float pkd1,pkd2,pkd3;

	struct ABS
                {
                 float pkd1,pkd2,pkd3;
                 float vdw,hb,hm,hp,hs,rt;
                 float score;
                };

	ABS *abs_inf;

        Ligand(); ~Ligand();
        Ligand(int max_atom_num, int max_bond_num);

        Ligand(const Ligand &original);
        Ligand& operator = (const Ligand &original);

        void Clear();				// override Molecule::()
	void Show_Contents() const;		// override Molecule::()

        int Value_Atom();			// override Molecule::()
	void Calculate_HB_Root();

	int Chemical_Viability_Check() const;
	float Check_Buried_Status(float cutoff=0.90);
	float Atom_Buried_Surface(int id, float &total, float &buried) const;
	float Calculate_Buried_Surface() const;

	// score.c

	float Calculate_Binding_Score();
		float Calculate_VDW();
		float Calculate_HB(char *pdb_entry="unknown");
			void Sum_HBonds(int num_candidate, 
				        HBond candidate[]) const;
			int Get_HBond_Pair_PL(HBond candidate[], 
				      	      bool sb_flag=false) const; 
			int Get_HBond_Pair_PP(HBond candidate[]) const;
		float Calculate_HP();
		float Calculate_HM();
		float Calculate_HS();
		float Calculate_RT();

		float Calculate_UHB() const;
		float Calculate_AR() const;

	void Write_Out_Log(char *filename) const;
	void Write_Out_Log(FILE *fp) const;
};

// **************************** protein.c *******************************

class Protein	// full version
{
 public:
	bool xtool_format;
	char name[256];
	float surface,bnsur,bpsur;

	int num_atom;
	vector <Atom> atom;

	int num_chain;
	vector <Chain> chain;

	int num_ring;
	vector <Ring> ring;

	vector <Dot> vol_dot;
	vector <Dot> sur_dot;

        Protein(); ~Protein();

        Protein(const Protein &original);
        Protein& operator = (const Protein &original);

	void Clear();
	void Show_Contents() const;
		void Show_Chains() const;
		void Show_Atoms() const;
		void Show_Rings() const;

        void Read_From_PDB(char *filename);
		void Read_ATOM(FILE *fp);
		void Read_SEQRES(FILE *fp);
		void Analyze_Sequence();
		void Read_CONECT(FILE *fp);

        void Value_Atom(int flag=1);  
		void Check_Atom_Type();
		void Detect_Connections();
        	void Calculate_HB_Root();

	void Fix_PDB_File(char *input, char *output);
		void Refine_Structure();

	void Rearrange_IDs();
	void Merge_Cofactor(const Ligand &cofactor);
	void Eliminate_H_in_PDB(char *input, char *output);

	void Write_Out_PDB(char *filename);
	void Write_Out_PDB_XTool(char *filename);
	void Write_Out_PDB_DrugScore(char *filename);

	void Define_Water(const Ligand *ligand, bool w_flag=false);
	void Define_Pocket(const Ligand *ligand, float cutoff=(DIST_CUTOFF));
	void Define_Pocket(float origin[], float radius);

	void Write_Pocket_PDB(char *filename) const;
	void Write_Pocket_XTOOL(char *filename) const;
	void Read_Pocket_XTOOL(char *filename);

	int Check_Buried_Ratio(float coor[]) const;

	void Generate_Volume_Dots(float spacing, float probe_r=0.00);
	float Calculate_Volume() const;
	void Show_Volume_Dots(char *filename, char *show="unit") const; 

	void Generate_Surface_Dots(const Ligand *ligand, float probe_r=0.00);
	float Calculate_Surface() const;
	void Show_Surface_Dots(char *filename, char *show="unit") const;

	void Classify_Residues(char *filename);
};

// ************************** population.c *********************************

class Population		// in multiple Mol2 format 
{
 private:
	Ligand *member;

        struct Record
		{
		 int id;
                 int valid;
                 char name[256];
                 char formula[256];
                 float weight;
                 float logp;
		 int num_hb_atom;
		 float hpscore;
		 float hmscore;
		 float hsscore;
                 float score;
		 int rank;
                };
	Record *record;

 public:
	int num_member;
	int num_output;

	Population(char *filename); ~Population();

	void Score_Members();
	void Rank_Members();
	void Output_Results();
};

// ************************ misc.c ********************************

float Distance(const float a[3], const float b[3]);
float Distance2(const float a[3], const float b[3]);
void Translate_Point(const float start[3], const float move[3], float end[3]);
void Rotate_Point(float theta, float axis[3], const float origin[3],
                  const float start[3], float end[3]);
void Cross_Multiply(const float v1[3], const float v2[3], float result[3]);
float Dot_Multiply(const float v1[3], const float v2[3]);
float Triple_Multiply(const float v1[3], const float v2[3], const float v3[3]);
float Torsion_Angle(const float p1[3], const float p2[3], 
		    const float p3[3], const float p4[3]);
void Unit_Vector(const float start[3], const float end[3], float result[3]);
void Normal_Vector(float v[3]);
float Angle(const float a[3], const float b[3], const float c[3]);
float Angle_of_Two_Vectors(const float v1[3], const float v2[3]);
float Angle_of_Two_Planes(const float p1[3], const float p2[3],
			  const float p3[3], const float q1[3],
			  const float q2[3], const float q3[3]);
float Distance_From_Point_To_Line(const float origin[3],
				  const float axis[3],
				  const float point[3]);
float Angle_of_Point_and_Plane(const float point[3],
			       const float centroid[3],
			       const float p1[3], 
			       const float p2[3],
			       const float p3[3]);	 
void Get_Sp3_Coordinates(const float v1[3], const float v2[3], 
			 float v3[3], float v4[3]);
void Get_Sp2_Coordinates(const float v1[3], const float v2[3], float v3[3]);

void Set_Random_Number();
char *Get_Time();
char *Chomp(char *string);
int Blank_Line_Check(const char *line);

int Pow(int x, int y);
void Memory_Allocation_Error(char *position="undefined");
void Open_File_Error(const char *filename);
void Read_File_Error(const char *filename);
void PDB_Format_Error(const char *filename);
void Mol2_Format_Error(const char *filename);
void Lig_Format_Error(const char *filename);
void Int_To_Char(int number, char string[], int digit=3);

int Check_Mol2_File(const char *filename);

float Calculate_RMSD(const Molecule &mol1, const Molecule &mol2);

void Sort_Array(int n, double array[], char *order="increasing");

void Simplex_Minimization(int ndim, int mpts, float **p, float y[], float ftol,
                          float (*funk)(int, float []), int &nfunk);
float Simplex_Generation(int ndim, int mpts, float **p, float y[], float psum[],
                         float (*funk)(int, float []), int ihi, float factor);

void Powell_Minimization(int n, float p[], float ftol,
			 int &iter, float &fret, float (*func)(int, float[]));
	void linmin(float p[], float xi[], int n, float &fret, 
	    	    float (*func)(int, float []));
	float f1dim(float x);
	void mnbrak(float &ax, float &bx, float &cx, 
		    float &fa, float &fb, float &fc,
            	    float (*func)(float));
	float brent(float ax, float bx, float cx, 
            	    float (*f)(float), float tol, float &xmin);

