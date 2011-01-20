from os import system, path, remove
import shutil

def run_modelgenerator(modelgenerator_path, alignment_path, num_gamma_cat = 4):
		
	aln_name = alignment_path.split("/")[-1]
	output_filename = "%s_modelgenerator.out" %(alignment_path)
	
	#1. run modelgenerator, if we haven't already run it for this alignment
	if not(path.isfile(output_filename)):
		print "Running ModelGenerator on %s" %(alignment_path)
		system("java -jar %s %s %d" %(modelgenerator_path, alignment_path, num_gamma_cat))
	
		#2. move and rename output file, and delete the other files 
		shutil.move("modelgenerator0.out", output_filename)
		remove("%s_phyml.sh" %(aln_name))
		remove("%s_phymlBoot.sh" %(aln_name))
		remove("%s_puzzleBoot.sh" %(aln_name))
		remove("%s_treePuzzle.sh" %(aln_name))
		
	else:
		print "modelgenerator has already been run for the file:\n'%s'\nthe program won't re-run it, but if you'd like to rerun it, delete the modelgenerator outputfile:'%s'" %(aln_name, output_filename)

	return output_filename

def extract_results_modelgenerator(filepath):

	#extract modelgenerator results as lists of tuples
	find_this = "\tAIC1\t\t\tModel\t\tLn\t\t\tAIC2\t\t\tModel\t\tLn\t\t\tBIC\t\t\tModel\t\tLn\n"
	file = open(filepath, "r")
	lines = file.readlines()
	for line in lines:
		if line == find_this:
			results = lines[lines.index(find_this)+3:lines.index(find_this)+60]
			break

	AIC_results = {}
	AICc_results = {}
	BIC_results = {}
	lnL_results = {}

	for line in results:
		line = line.split("\t")
		while line.count(""): #sometimes there are double tabs in the line, which gives empty strings in results
			line.remove("")

		#stored as: [modelname] = score
		AIC_results[str(line[2])] = float(line[1])
		AICc_results[str(line[5])] = float(line[4])
		BIC_results[str(line[8])] = float(line[7])
		lnL_results[str(line[2])] = float(line[3])
		
	return AIC_results, AICc_results, BIC_results, lnL_results

if __name__ == "__main__":
	
	curdir = path.dirname(path.abspath(__file__))
	print curdir
	
	alignment_path = "%s/testfiles/test.fas" %(curdir)
	modelgenerator_path = "%s/programs/modelgenerator.jar" %(curdir)

	results_file = run_modelgenerator(modelgenerator_path, alignment_path)	

	print alignment_path, modelgenerator_path
	
	test = extract_results_modelgenerator("/Users/Rob/Dropbox/Current_work/06Partition_Finder/v0.2/partition_alignments/1_models.out")

	for metric in test:
		for entry in test[metric]:
			print entry