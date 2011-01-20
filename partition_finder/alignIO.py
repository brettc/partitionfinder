def import_fasta_as_dict(alignment_filename):
	######################        Import fasta file as dict        #################################
	alignment_file = open(alignment_filename, 'r')
	#and set up the alignment as a dictionary
	seqdict = {}
	#get the seqs from the alignment file
	start = 0 # a cheap trick to do this simply
	for line in alignment_file:
		if line.startswith(">"): #it's a new sequence
			#first, stash the old sequence into the dictionary
			if start == 1:
				seqdict[seq_name] = seq
			seq_name = line
			seq = '' #an empty sequence
		else: #we must have a sequence, maybe with a newline on the end of it
			seq_part = line.strip("\n") #remove newlines as we go
			seq = ''.join([seq, seq_part])
		start=1
	#add the last sequence, which you'll have missed because I'm too lazy to write this properly
	seqdict[seq_name] = seq	
	alignment_file.close()
	
	return seqdict
