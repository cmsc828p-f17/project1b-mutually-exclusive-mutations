# project1b_permutation
# doing the permutation test for a given weight/rank of a function
# rank is order on weight list,

# basic idea - with permutation, null is we expect to see the same amount of 
# runs where we find at least one set with that weight when sampling ervery x# of iterations
# so if we do not see a solution of that weight in the run, then we do not add that to the count

# run as python project1b_permutation.py data\examples\example-mutation-matrix.tsv 3 100000 100 47 1
# or, python project1b_permutation.py project1b_permutation.py nputdata setsize numiter numruns weight rank

# inputdata - the mutation input data text file
# setsize - parameter for the size of the mutually exclusive gene set to test, usually 2 or 3
# numiter - number of iterations to run Metropolis-Hasting algorithm
# numruns - number of runs you want to consecutively run this algorithm for, so for diffeent random initializations
# weight - weight of statistic we wish to test
# rank - how many times we expect to see this weight, as a null hypothesis
# rank comes from the postion on the list of highest weights of sets seen
# we assume, say, if the rank was one, if that was truly the optimal solution and the null hypothesis is wrong,
# then we shoud lnever see a weight that high are randomly permuting the samples with a given gene mutated

# note no step size, we sample every iteration now

# copyright stuff again
#Copyright 2010,2011,2012,2013 Brown University, Providence, RI.

                         #All Rights Reserved

#Permission to use, copy, modify, and distribute this software and its
#documentation for any purpose other than its incorporation into a
#commercial product is hereby granted without fee, provided that the
#above copyright notice appear in all copies and that both that
#copyright notice and this permission notice appear in supporting
#documentation, and that the name of Brown University not be used in
#advertising or publicity pertaining to distribution of the software
#without specific, written prior permission.

#BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
#INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
#PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
#ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
#WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
#ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
#OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#http://cs.brown.edu/people/braphael/software.html


# import packages
import sys
import os

import random
import math


# define weight functiion
def weight(solution):
    lambdaG = 0
    lambdaM = 0
    examples_seen = set()
    for gene in solution:
        # add to coverage overlap
        lambdaG += len(genes_setexamples[gene])
        for example in genes_setexamples[gene]:
            # add to coverage (only unique patients)
            if example not in examples_seen:
                examples_seen.add(example)
                lambdaM += 1
                
    return float(2*lambdaM - lambdaG)

if __name__ == "__main__":
    # read in input arguments
    inputfile = sys.argv[1]
    setsize = int(sys.argv[2])
    numiter = int(sys.argv[3])
    numruns = int(sys.argv[4])
    testweight = int(sys.argv[5])
    rank = int(sys.argv[6])
    
    # creating sets/dictionaries to use later
    genes_setexamples = dict()
    examples_setgenes = dict()
    numvisits = dict()
    exampleslist = set()
    geneslist = set()
    
    # upload the example files
    mutationdata = open(inputfile,'r')
    # parse examples
    for example in mutationdata:
        line = example.split()
        samplepatient = line[0]
        exampleslist.add(samplepatient)
        examples_setgenes[samplepatient] = set()
        
        for i in range(len(line)-1):
            gene = line[i+1]
        
            #add this gene to the list of mutated genes for that example
            examples_setgenes[samplepatient].add(gene)
            
            # add the patient to the set of patients with the given gene mutated
            if gene not in genes_setexamples:
                genes_setexamples[gene] = set()
                geneslist.add(gene)
                
            genes_setexamples[gene].add(samplepatient)
        # end of i
    # end of example        
    mutationdata.close()
# run the Metroplis-Hasting experiment
    random.seed()
# setting c =-0.5 like they empirically found in the paper
    c = 0.5
    
# opening file to write now
    
# this set-up above is the same as before, now we will do the permutation test
# not sure if we can cut out making the examples_setgenes, I do not think so since we need the # of patients
# but I keep it since this is how a formulate the gene_Setexamples set, so I have to keep it for now
    runs_greaterweight = 0
    for runs in range(numruns):
        examples_setgenes = dict()
        num_greaterweight = 0
        # randomly permute the patients with each gene mutated
        for gene in genes_setexamples:
            permset = set(random.sample(exampleslist,len(genes_setexamples[gene])))
            genes_setexamples[gene] = permset
            # update the example_setgenes list as well for new set of pateints
            for example in permset:
                # adding key to dictionary
                if example not in examples_setgenes:
                    examples_setgenes[example] = set()
                examples_setgenes[example].add(gene) # adding the gene to the lsit of genes mutated in each patient
                # since we need to update both sets
        
        # initialize a random set to start
        state = random.sample(geneslist, setsize)
        state = set(state)
        # now run same Metropolis-Hasting algorithm
        for ii in range(numiter):
            # NOTE: THIS IS DIFFERENT IN THE RANDOM PERMUTATION TEST FILE FROM THE PAPER
            # While I do not know why, I will do the same thing. I will talk to you about it on slack
#             geneswapin = random.sample(geneslist,1)
#             if geneswapin[0] in state:
#               geneswapout = geneswapin
#             else:
#                geneswapout = random.sample(state,1)
    
            # doing the switching of genes like they do in the sample permutation test, although is this 
            # a smart thing to do? I guess so, since it insinuates a new set to test each time
            # meaning we can test more potential gene sets?
             genes_notinstate = geneslist.difference(state)
             geneswapin = random.sample(genes_notinstate,1)
             geneswapout = random.sample(state,1)
                
            # swap that given gene 
             newstate = state.difference(geneswapout)
             newstate = newstate.union(geneswapin)
        
        # assess the measure of both old solution and new solution
             stateweight = weight(state)
             newstateweight = weight(newstate)
        
        # compute exponent
             weightprob = math.exp(c*(newstateweight-stateweight))
        # paper states they take a minium of this and 1, so we will do it as well
             prob = min(1.0,weightprob)
        
        # do random prob, switch if the random prob is less that the probabilty found above
             switch = random.random()
             if switch <= prob: # should not matter <= or < due to float number
                state = newstate
                stateweight = newstateweight
                
            # checking to see if the new state has a higher weight than the testweight
             if stateweight >= testweight:
                num_greaterweight += 1
            # end iterations loop
            
        if num_greaterweight >= rank:
             runs_greaterweight += 1
             
        print('Done Run %s' % str(runs))
        
        # done runs loop
    # finding the empirical p-value which is the fraction of runs where we saw a gene set with at least the
    # same weight as the testweight in mroe number if iterations than the rank we had
    # this tests the hypothesis of if that weight would be seen randomly
    # still iffy in the validity of this statistic, since it allows for double-counting 
    # of a solution with high weight if seen in a lot of iterations
    empirical_pvalue = float(runs_greaterweight/numruns)
    # writing to output file
    out_perm = open('permutationtest_weight'+str(testweight)+'_rank'+str(rank)+'.txt','w')
    out_perm.write('Weight\t'+str(testweight)+'\tRank\t'+str(rank)+'\tNumber of Runs\t'+str(numruns)+'\n')
    out_perm.write('P-value:\t'+str(empirical_pvalue)+'\n')
    out_perm.close()
                
