#import sys
import os
import random
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

#working directory
os.chdir('/home/henrike/Documents/PD_AS/teaching/advanced_topics_binf/random_mut_exercise/')

#function for translating
def translate_DNA(DNA):
	protein = ''
	for i in range(0,len(DNA),3):
		protein += codon_table[DNA[i:i+3]]
	
	return(protein)

#mutate one position only
# ~ def mutate(DNA):
	# ~ #pick a random position:
	# ~ pos = random.randint(0,len(DNA)-1)
	# ~ #pick a random base to mutate into
	# ~ new_base = random.choice(change_base[DNA[pos]])
	# ~ #substitute old base
	# ~ new_DNA = DNA[:pos] + new_base + DNA[pos+1:]
	# ~ return(new_DNA)

#mutate chosen number of positions	
def mutate(DNA, nr_muts=1):
	new_bases = {}
	for i in range(nr_muts):
		#we don't want to chose the same position seveal times (though this of course happens in nature, but I actually want x number of mutated positions)
		while True:
			#pick a random position:
			pos = random.randint(0,len(DNA)-1)
			if not pos in new_bases:
				#pick a random base to mutate into
				new_bases[pos] = random.choice(change_base[DNA[pos]])
				break
			else:
				continue	
		
	#sort dict by positions and substitute in the changes
	sorted_positions = [int(k) for k, v in sorted(new_bases.items())]
	
	new_DNA = DNA[:sorted_positions[0]]
	for i in range(len(sorted_positions)-1):
		#mind the placement of the 1 inside the [] (meaning the next chosen position!) or outside (meaning the current chosen position + 1)
		new_DNA += new_bases[sorted_positions[i]] + DNA[sorted_positions[i]+1:sorted_positions[i+1]]
	#add last position and go until end of original seq
	new_DNA += new_bases[sorted_positions[-1]] + DNA[sorted_positions[-1]+1:]
	
	return(new_DNA)

def prot_diff(WT, new_prot):
	changes = []
	for pos,aa in enumerate(WT):
		if WT[pos] != new_prot[pos]:
			#changes.add(aa+str[pos]+new_prot[pos])	
			changes.append((aa,pos,new_prot[pos]))	
			
	return(changes)

def print_table(muts, filename):
	#you can probable find a smarter way to do this, I just hacked it
	with open(filename, 'w') as OUT:
		print(','.join(order_AA), file = OUT)
		for AA_WT in order_AA:
			print(AA_WT, end = '', file = OUT)
			for AA_mut in order_AA:
				if AA_mut == AA_WT:
					counts = 'NA'
				else:
					counts = muts[AA_WT+AA_mut] if str(AA_WT+AA_mut) in muts else 0
				print(',', counts, sep = '', end = '', file = OUT)
			print('\n', end = '', file = OUT)

def make_heatmap(table_name, output_name):
	my_df = pd.read_csv(table_name, sep = ',')
	#https://seaborn.pydata.org/generated/seaborn.heatmap.html
	ax = sns.heatmap(my_df, cmap='YlOrRd')
	#https://stackoverflow.com/questions/34232073/seaborn-heatmap-y-axis-reverse-order
	ax.invert_yaxis()
	plt.savefig(output_name)
	#clear the plot so we can make another one
	plt.clf()
	#plt.show()

#codon table: adapted from https://pythonforbiologists.com/dictionaries
codon_table = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}

#bases = ['A', 'T', 'G', 'C']
#one level up: force to choose a different base than the current when mutating
#two levels up: actually transitions (A<->G, C<->T) occur at a much higher freq than transversions so you could choose with some prior
change_base = {
	'A' : ['T', 'C', 'G'],
	'T' : ['A', 'C', 'G'],
	'G' : ['A', 'T', 'C'],
	'C' : ['A', 'T', 'G']		
}

order_AA = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']

#read in file:
WT_nt = ""
with open('native_DNA.fa', 'r') as IN:
	for line in IN:
		if not line.startswith('>'):
			#make sure DNa is in upper case since we use that in the dict and remove the trailing newline
			WT_nt += line.rstrip().upper()

#translate to protein:
WT_prot = translate_DNA(WT_nt)
print(WT_nt)
print(WT_prot)

#introduce 1 single mutation (on nucleotide level, that's where mutations happen!) a 1000 times
n_iter = 1000
muts = {}
for i in range(n_iter):
	#mutate
	new_DNA = mutate(WT_nt)
	#translate
	new_prot = translate_DNA(new_DNA)
	#save the changes to a dict
	for change in prot_diff(WT_prot,new_prot):
		if not str(change[0]+change[2]) in muts:
			muts[str(change[0]+change[2])] = 1
		else:
			muts[str(change[0]+change[2])] += 1

#write to file
print_table(muts, 'single_muts2.csv')

#heatmap
make_heatmap('single_muts2.csv','single_muts2.png')

#introduce 3 new mutations
for i in range(n_iter):
	#mutate
	new_DNA = mutate(WT_nt, 3)
	#translate
	new_prot = translate_DNA(new_DNA)
	#save the changes to a dict
	for change in prot_diff(WT_prot,new_prot):
		if not str(change[0]+change[2]) in muts:
			muts[str(change[0]+change[2])] = 1
		else:
			muts[str(change[0]+change[2])] += 1
			
#write to file
print_table(muts, 'three_muts.csv')
#heatmap
make_heatmap('three_muts.csv','three_muts.png')	
	
#new_DNA = mutate(WT_nt, 3)
#print(new_DNA)
#print(translate_DNA(new_DNA))

#introduce 5
