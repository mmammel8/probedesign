#>Rosalind_79
#CAGCACTTGGATTCTCGG
#>Rosalind_98
#CAGCGTGG
INNAME = "rosalind_smgb9.txt"

def read():
	lineno = 0
	with open(INNAME) as f:
		for line in f: #single line sequences FASTA
			if lineno == 1:
				dna1 = line.strip()
			elif lineno == 3:
				dna2 = line.strip()				
			lineno += 1
	f.close()
	return dna1, dna2

def swater(dna1, dna2):
	len1 = len(dna1) + 1
	len2 = len(dna2) + 1
	table = [ [0 for _ in range(len1)] for _ in range (len2)]
	optbl = [ ["" for _ in range(len1)] for _ in range (len2)]
	scores = [1,-1,-1,-1,0]

	print(len1, len2)
	for i in range(len1): #initialize upper row
		table[0][i] = 0
		optbl[0][i] = "k"
	for j in range(len2): #init left column
		table[j][0] = 0
		optbl[j][0] = "j"
	j = 1
	while j < len2:
		i = 1
		while i < len1:
			maxv = 0
			c1 = 0
			if dna1[i-1] == dna2[j-1]:
				maxv = table[j-1][i-1] + scores[0] #copy
			else:
				maxv = table[j-1][i-1] + scores[1] #replace
			op = "c"
			if i < len1 - 1:
				c1 = table[j-1][i] + scores[2] #insert
				if c1 > maxv:
					maxv = c1
					op = "i"
			else:
				c1 = table[j-1][i] #kill
				if c1 > maxv:
					maxv = c1
					op = "k"
			if j < len2 - 1:
				c1 = table[j][i-1] + scores[3] #delete
				if c1 > maxv:
					maxv = c1
					op = "d"
			else:
				c1 = table[j][i-1] #kill
				if c1 > maxv:
					maxv = c1
					op = "j"
			table[j][i] = maxv
			optbl[j][i] = op
			i += 1
		j += 1

	startx = len1 - 1
	starty = len2 - 1
	aln1 = ""
	aln2 = ""
	maxv = table[starty][startx] #last entry
	#print("score =", maxv)
	savemin = maxv
	i = startx
	j = starty #backtrack
	mismatch = 0
	matched_seq = ""
	while i > 0 and j > 0:
		#find which operation was used to get here
		op = optbl[j][i]
		if op == "c":
			#copy/replace
			i -= 1
			j -= 1
			aln1 = dna1[i] + aln1
			aln2 = dna2[j] + aln2
			matched_seq = dna1[i] + matched_seq
			if dna1[i] != dna2[j]:
				mismatch += 1
		elif op =="i":
			#insert
			j -= 1
			aln1 = '-' + aln1
			aln2 = dna2[j] + aln2
			matched_seq = '-' + matched_seq			
			mismatch += 1
		elif op == "k":
			j -= 1	
		elif op== "d":
			#delete
			i -= 1
			aln1 = dna1[i] + aln1
			aln2 = '-' + aln2
			matched_seq = dna1[i] + matched_seq
			mismatch += 1
		elif op == "j":
			i -= 1					
 
	print(aln1)
	print(aln2)
	out_file = open("table.txt", 'w')
	out = ""
	for j in range(len2):
		out = ""	
		for i in range(len1):
			#out += str(table[j][i]) + "\t"
			out += optbl[j][i] + "\t"
		out_file.write(out + "\n")
	out_file.close()
	return mismatch, matched_seq

dna1, dna2 = read()
mm, aln2 = swater(dna1, dna2)
print("mismatches", mm)
print(aln2)

