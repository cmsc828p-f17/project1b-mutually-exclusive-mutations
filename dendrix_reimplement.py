import math
import argparse
import random
import operator
import sys

genes = set()
patients = set()
patient2genes = {}
genes2patients = {}
c = 0.5

'''
This function takes 3 arguments:
1. Gene set (M)
2. patient to mutated genes mapping(patient2genes)
3. genes to patients in which that gene is mutated mapping (genes2patients)

Returns the weight based on the formula in the paper
'''
def weight(M,patient2genes,genes2patients):
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

'''
This function reads the input file containing mutated genes information.
Stores information as:
genes: all unique genes
patients: all unique patients
patients2genes: patient to mutated genes mapping
genes2patients: genes to patients in which that gene is mutated mapping 
'''
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


'''
Standalone MCMC function, that takes following arguments:
patients2genes: patient to mutated genes mapping
genes2patients: genes to patients in which that gene is mutated mapping 
iterations: Number of iterations for MCMC
k: gene set size
weight_cutoff: weight being tested in case of the permutation test
'''
def MCMC(patient2genes_sample,genes2patients_sample,iterations,k,weight_cutoff=0):
    '''
    Initial gene set
    '''    
    curr_answer = set(random.sample(genes,k))
    solution_map = {}
    solution_map_all = {}
    for i in xrange(0,iterations):
        '''
        Randomly choose a gene from all genes
        '''
        w = random.sample(genes,1)

        if w[0] in curr_answer:
            v = set(w)
        else:
            v = set(random.sample(curr_answer,1))

        '''
        to_test is the set by adding the gene to test w and removing a random gene v from the gene set
        '''
        to_test = (curr_answer.difference(v)).union(w)
        to_test_weight = weight(to_test,patient2genes_sample,genes2patients_sample)
        curr_answer_weight = weight(curr_answer,patient2genes_sample,genes2patients_sample)
        '''
        Compute the probability using the formula in the paper
        '''
        exp_prob = math.exp(c*(to_test_weight - curr_answer_weight))

        '''
        Choose the minimum of 1 and the computed probability
        '''
        P = min(1.0,exp_prob)
        samp_prob = random.random()
        '''
        If P is highm then change the current set to the set being tested
        '''
        if samp_prob < P:
            curr_answer = to_test
            curr_weight = weight(curr_answer,patient2genes,genes2patients)
            if curr_weight >= weight_cutoff:
                if frozenset(curr_answer) not in solution_map:
                    solution_map[frozenset(curr_answer)] = 0
                
                '''
                Count the number of times a set is visited
                '''
                solution_map[frozenset(curr_answer)] += 1
            fset = frozenset(curr_answer)
            if fset not in solution_map_all:
                solution_map_all[fset] = 0
            solution_map_all[fset] += 1
    

    return solution_map,solution_map_all

def permutation_test(rounds,weight_cutoff,k,rank,iterations):

    consider = 0

    for i in xrange(0,rounds):
        '''
        Randomly permute the patients for each mutated genes, keeping the number of patients same
        Also compute other data structure based on this 
        '''
        patient2genes_sample = {}
        genes2patients_sample = {}
        for gene in genes2patients:
            rand_subset = random.sample(patients,len(genes2patients[gene]))
            genes2patients_sample[gene] = set(rand_subset)
            for patient in rand_subset:
                if patient not in patient2genes_sample:
                    patient2genes_sample[patient] = set()
                patient2genes_sample[patient].add(gene)

        '''
        Call MCMC function
        '''
        solution_map,_ = MCMC(patient2genes,genes2patients,iterations,k,weight_cutoff)

        '''
        Counts the number of time the set has higher weight than the rank being tested
        '''
        if sum(solution_map.values()) >= rank:
            consider += 1

        

    pvalue = consider*1.0/rounds
    # print "PValue for the permutation test = " + str(pvalue)
    return pvalue





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--matrix',help="Mutation Matrix for Patients")
    parser.add_argument('-k','--size',help="Size of the gene set sampled")
    parser.add_argument('-i','--iteration',help="Number of Iteration")
    parser.add_argument('-c','--exponent',help="Exponent parameter", default=0.5,required=False)
    parser.add_argument('-n','--num',help="Number of sets to output", default=100,required=False)
    parser.add_argument('-o','--prefix',help="output prefix ", default="OUT",required=False)
    parser.add_argument('-p','--permut',help="set this to true to perform permutation test",default=False,required=False)
    parser.add_argument('-r','--rounds',help="Number of rounds to be performed for the permutation test",default=100,required=False)
    parser.add_argument('-w','--weight',help="Weight to be used for the permutation test",default = 10, required=False)
    parser.add_argument('-x','--rank',help="Rank to be used for the permutation test",default = 10, required=False)

    args = parser.parse_args()
    iterations = int(args.iteration)
    '''
    Function to load the data
    '''
    read_matrix(args.matrix)
    c = float(args.exponent)
    curr_answer = set(random.sample(genes,int(args.size)))
    k = int(args.size)

    solutions_visited = {}

    '''
    If permutation test needs to be performed, call the method.
    '''
    if args.permut:
       permutation_test(int(args.rounds),int(args.weight),int(args.size),int(args.rank),int(args.iteration))
       sys.exit(1)

    '''
    Call MCMC function
    '''
    _,solutions_visited = MCMC(patient2genes,genes2patients,iterations,k)


    sorted_solutions_visited = sorted(solutions_visited.items(),key=operator.itemgetter(1),reverse=True)

    solutions_weight = {}
    for key in solutions_visited:
        solutions_weight[frozenset(key)] = weight(key,patient2genes,genes2patients)

    sorted_solutions_weights = sorted(solutions_weight.items(),key=operator.itemgetter(1),reverse=True)

    '''
    Output the sets visited in the decreasing order of the number of times they were visited during MCMC
    '''
    to_output = int(args.num)
    count = 0
    ofile = open(args.prefix+'_byNumber.txt','w')
    for key in sorted_solutions_visited:
        for each in key[0]:
            ofile.write(each+'\t')
        ofile.write(str(key[1]))
        ofile.write('\n')

    ofile.close()

    '''
    Output the sets visited in the decreasing order of the number of weights
    '''
    count = 0
    ofile = open(args.prefix+'_byWeight.txt','w')
    for key in sorted_solutions_weights:
        for each in key[0]:
            ofile.write(each+'\t')
        ofile.write(str(key[1]))
        ofile.write('\n')

    ofile.close()

    '''
    Output permutation test results for gene sets
    Foramat: Gene set, weight, sampling frequency, p-value
    '''
    ofile = open(args.prefix+'_'+str(k)+'_stats.txt','w',1)
    rnk = 1
    for key in sorted_solutions_weights:
        if rnk <= 20:
            print 'set with rank '+ str(rnk) + ' testing!'
            pvalue = permutation_test(1000,key[1],k,rnk,iterations)
            line = ''
            for each in key[0]:
                line += each+'\t'
            line += str(key[1])+'\t'
            line += str(solutions_visited[key[0]])+'\t'+str(pvalue)+'\n'
            ofile.write(line)
            print 'set with rank '+ str(rnk) + ' done!'
            rnk += 1
        else:
            break
    ofile.close()

if __name__ == '__main__':
    main()
