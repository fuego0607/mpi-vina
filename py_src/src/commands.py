from subprocess import call
import os

def runVina(vina_exe, conf_file, receptor, ligand, scores_out):
	#"${p}/vina --config ${conf_file} --ligand ${ligand} --out ${pdbqt} > ${stdout} 2> ${stderr}
	command = "{0} --config {1} --receptor {2} --ligand {3} --out {4} > /dev/null 2>&1".format(vina_exe, conf_file, receptor, ligand, scores_out)
	call(command, shell=True)

def convertToMol2(vinaScores_pdbqt, scores_mol2):
	#"obabel -i pdbqt path_to_pdbqt -o mol2 -O path_to_new_mol2 -h"
	command = "obabel -i pdbqt {0} -o mol2 -O {1} -h > /dev/null 2>&1".format(vinaScores_pdbqt, scores_mol2)

	call(command, shell=True)

def runDSXScore(dsx_exe, receptor_pdb, vinaScore_mol2, potentials_dir, dsx_out):
	#"dsx_linux_64.lnx -P ${receptor_pdb} -L ${mol2} -D ${potentials} -F ${dsx_dir}/${out_name}.out"
	command = "{0} -P {1} -L {2} -D {3} -F {4} > /dev/null 2>&1".format(dsx_exe, receptor_pdb, vinaScore_mol2, potentials_dir, dsx_out)

	call(command, shell=True)

def runXScore(xscore_exe, receptor_pdb, pdb_fix, vinaScore_mol2, mol2_fix, output_file):
	#"xscore -fixpdb ${receptor_pdb} ${pdb_fix}"
	#convert receptor from .pdb to .pdb_fixed using xscore
	pdb_command = "{0} -fixpdb {1} {2} > /dev/null 2>&1".format(xscore_exe, receptor_pdb, pdb_fix)

	call(pdb_command, shell=True)
	#"xscore -fixmol2 ${mol2} ${mol2_fix}"
	#convert vina score from .mol2 to .mol2_fixed using xscore
	mol2_command = "{0} -fixmol2 {1} {2} > /dev/null 2>&1".format(xscore_exe, vinaScore_mol2, mol2_fix)

	call(mol2_command, shell=True)

	#update and write config file to disk

	command = "{0} -score {1} {2} > {3}".format(xscore_exe, pdb_fix, mol2_fix, output_file)

	call(command, shell = True)

def runNNScore(nns_exe, receptor_pdbqt, vinaScore_pdbqt, networks_dir, outfile):
	#specify output
	#specif neural networks directory
	#python NNScore.py -receptor {receptor} -ligand {ligand} -networks_dir > outfile
	command = "python {0} -receptor {1} -ligand {2} -networks_dir {3} > {4}".format(nns_exe, receptor_pdbqt, vinaScore_pdbqt, networks_dir, outfile)

	call(command, command=True)

def runRFScore(rfs_exe, receptor, vinaScores_pdbqt, outfile):
	#callRFScore
	#vina pdbqt == ligands
	#rf-score-vs --receptor protein.pdb ligands.sdf -o csv --field "name" --field "RFScoreVS_v2"
	command = "{0} --receptor {1} {2} -o csv --field \"name\" --field \"RFScoreVS_V2\" > {3}".format(rfs_exe, receptor, vinaScores_pdbqt, outfile)

	call(command, command=True)




