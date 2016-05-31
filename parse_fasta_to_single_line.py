#
# Builds new text files for fasta files
# Creates text file with fasta sequence in single line without header
#

#Function to read in fasta and flatten it into a single line txt file
def flatten_fasta(id, organism):
	#Reads fasta file
	file_dir = "Data/"
	fasta_file = file_dir + id + ".fasta"
	f = open(fasta_file, "r")

	#Array to take in data
	ref = []

	#Append all the lines from fasta to array
	for line in f:
		ref.append(line.rstrip())

	#Close connection
	f.close()

	#Remove the header
	ref = ref[1:len(ref)]

	#Join array pieces together to single line
	ref = "".join(ref)

	#Output single line fasta
	output_file = "Output/" + organism + "_raw.txt"
	f = open(output_file, 'w')
	f.write(ref)
	f.close()
