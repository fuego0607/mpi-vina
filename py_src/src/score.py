from subprocess import call
import os

def runNNScore(receptors_pdbqt, vinaScores_pdbqt):
	for index, score in enumerate(vinaScores_pdbqt):
		call(["python", "NNScore.py", "-receptor", receptors_pdbqt[index], "-ligand", score, "-vina_executable", "/Vina/vina"])

def runDSXScore(receptors_pdb, vinaScores_mol2):
	#build dsx_output_directory for result output
	potentials = "/pdb_pot_0511/"

	for index, score in enumerate(vinaScores_mol2):
		call(["dsx_linux_64.lnx", "-P", receptors_pdb[index], "-L", score, "-D", potentials, "-F", "dsx_output_directory"])

def runRFScore():
	#callRFScore

def runXScore(receptors_pdb, vinaScores_mol2):
	#convert receptor from .pdb to .pdb_fixed using xscore
	for receptor in receptors_pdb:
		call(["xscore", "-fixpdb", receptor, "pdb_fix_directory"])
	
	#convert vina score from .mol2 to .mol2_fixed using xscore
	for score in vinaScores_mol2:
		call(["xscore", "-fixmol2", score, "mol2_fix_directory"])
	
	#converted_score = ** get list of all mol2_fixed scores in mol2_fix_directory **
	#for converted_score in converted_scores:
		#create score specific config file using base config file template and path names
		#https://github.com/elmotec/massedit seems like a good replacement for sed
		#run xscore with input file
		#call(["xscore", input_file])

def convertToMol2(vinaScores_pdbqt):
	#need to define mol2 directory
	for score in vinaScores_pdbqt:
		call(["obabel", "-i", "pdbqt", score, "-o", "mol2", "-O", "mol2_directory", "-h"])


if __name__ == '__main__':
	#get receptor pdbqt dir, need to change this to a pdbqt-specific receptor directory
	receptors_pdbqt = os.listdir("/receptors")
	#get receptor pdb dir, need to establish a receptor pdb file
	#get vina scores pdbqt dir
	vinaScores_pdbqt = os.listdir("/ligands")
	convertToMol2(vinaScores_pdbqt)
	#vinaScores_mol2 = os.listdir("/ligands_mol2")  *** need to establish ligands mol2 directory

	#runNNScore(receptors_pdbqt, vinaScores_pdbqt)
	#runDSXScore(receptors_pdb, vinaScores_mol2)
	#runRFScore
	#runXScore(receptors_pdb, vinaScores_mol2)

	#all scoring methods assume receptors and lignd/vinaScore lists are built such that
	#the elements at the same index correspond
	#THIS IS A WRONG ASSUMPTION
	#For loops inside scoring functions need to be adjusted
	#need to map each ligand/vinaScore to its receptor
