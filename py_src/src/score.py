from subprocess import call
import os

def runVina(vinaConf, ligand):
	call(["vina", "--config", vinaConf, "--ligand", ligand, "--out", conf[vinaOutput], ">", "${stdout} 2> ${stderr}"])

def runNNScore(receptor_pdbqt, vinaScore_pdbqt, conf):
	#specify output
	#specif neural networks directory
	call(["python", "NNScore.py", "-receptor", receptor_pdbqt, "-ligand", vinaScore_pdbqt, "-networks_dir", "./networks/top_3_networks/"])

def runDSXScore(receptor_pdb, vinaScore_mol2, conf):
	#specify output
	potentials = "/pdb_pot_0511/"
	call(["dsx_linux_64.lnx", "-P", receptor_pdb, "-L", vinaScore_mol2, "-D", potentials, "-F", conf["output"]])

def runRFScore():
	#callRFScore

def runXScore(receptor_pdb, vinaScore_mol2, conf, input_file):
	#convert receptor from .pdb to .pdb_fixed using xscore
	call(["xscore", "-fixpdb", receptor, "pdb_fix_directory"])
	#convert vina score from .mol2 to .mol2_fixed using xscore
	call(["xscore", "-fixmol2", score, "mol2_fix_directory"])
	call(["xscore", input_file])

def convertToMol2(vinaScores_pdbqt):
	#need to define mol2 directory
	call(["obabel", "-i", "pdbqt", score, "-o", "mol2", "-O", "mol2_directory", "-h"])
