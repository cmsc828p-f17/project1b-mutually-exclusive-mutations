import math
import argparse
import random
import operator

genes = set()
patients = set()
patient2genes = {}
genes2patients = {}

def weight(M):
	gamma_of_g = 0
	for g in M:
		if g in genes2patients:
			gamma_of_g += len(genes2patients[g])
	gamma_of_m = 0
	for patient in patient2genes:
		to_add = False
		for gene in patient2genes[patient]:
			if gene in M:
				to_add = True
				break
		if to_add:
			gamma_of_m += 1

	return 2*gamma_of_m - gamma_of_g

def read_matrix(name):
	with open(name,'r') as f:
		for line in f:
			attrs = line.split()
			if len(attrs) >= 2:
				patient2genes[attrs[0]] = set(attrs[1:])
				for gene in attrs[1:]:
					genes.add(gene)
					if gene not in genes2patients:
						genes2patients[gene] = set()
					genes2patients[gene].add(attrs[0])
				patients.add(attrs[0])




def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-m','--matrix',help="Mutation Matrix for Patients")
	parser.add_argument('-k','--size',help="Size of the gene set sampled")
	parser.add_argument('-i','--iteration',help="Number of Iteration")
	parser.add_argument('-c','--exponent',help="Exponent parameter", default=0.5,required=False)
	parser.add_argument('-n','--num',help="Number of sets to output", default=100,required=False)
	parser.add_argument('-o','--prefix',help="output prefix ", default="OUT",required=False)
	
	
	args = parser.parse_args()
	iterations = int(args.iteration)
	read_matrix(args.matrix)
	c = float(args.exponent)
	curr_answer = set(random.sample(genes,int(args.size)))

	solutions_visited = {}

	for i in xrange(0,iterations):
		w = random.sample(genes,1)

		if w[0] in curr_answer:
			v = set(w)
		else:
			v = set(random.sample(curr_answer,1))	
		to_test = (curr_answer.difference(v)).union(w)
		# print to_test
		# print len(to_test)
		to_test_weight = weight(to_test)
		curr_answer_weight = weight(curr_answer)
		exp_prob = math.exp(c*(to_test_weight - curr_answer_weight))
		# print exp_prob
		P = min(1.0,exp_prob)
		samp_prob = random.random()
		if samp_prob < P:
			curr_answer = to_test

		fset = frozenset(curr_answer)
		if fset not in solutions_visited:
			solutions_visited[fset] = 0
		solutions_visited[fset] += 1


	sorted_solutions_visited = sorted(solutions_visited.items(),key=operator.itemgetter(1),reverse=True)
	
	solutions_weight = {}
	for key in solutions_visited:
		solutions_weight[frozenset(key)] = weight(key)

	sorted_solutions_weights = sorted(solutions_weight.items(),key=operator.itemgetter(1),reverse=True)

	to_output = int(args.num)
	count = 0
	ofile = open(args.prefix+'_byNumber.txt','w')
	for key in sorted_solutions_visited:
		for each in key[0]:
			ofile.write(each+'\t')
		ofile.write(str(key[1]))
		ofile.write('\n')

	ofile.close()

	count = 0
	ofile = open(args.prefix+'_byWeight.txt','w')
	for key in sorted_solutions_weights:
		for each in key[0]:
			ofile.write(each+'\t')
		ofile.write(str(key[1]))
		ofile.write('\n')

	ofile.close()

if __name__ == '__main__':
	main()