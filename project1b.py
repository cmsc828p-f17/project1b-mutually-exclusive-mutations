# Project 1b
# Graham Antoszewski

# run as project1b.py inputdata setsize numiter numexper stepsize
# use command line inputs to mirror the previous Dendrix file, but will happily fix this for simplicity if you want

# Example of how it is run on my personal computer:
# python project1b.py data\examples\example-mutation-matrix.tsv 3 100000 1 1000

# inputdata - the mutation input data text file
# setsize - parameter for the size of the mutually exclusive gene set to test, usually 2 or 3
# numiter - number of iterations to run Metropolis-Hasting algorithm
# numruns - number of runs you want to consecutively run this algorithm for, so for diffeent random initializations
# stepsize - number of iteration between counting the current state of the Markov chain


# putting this because of copyright
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
    

# define main
if __name__ == "__main__":
    # read in input arguments
    inputfile = sys.argv[1]
    setsize = int(sys.argv[2])
    numiter = int(sys.argv[3])
    numruns = int(sys.argv[4])
    stepsize = int(sys.argv[5])
    
    # creating sets/dictionaries to use later
    genes_setexamples = dict()
    examples_setgenes = dict()
    numvisits = dict()
    # examples list is meant for the permutation test for later to randomly sample, doing it now for copy and paste purposes
    exampleslist = set()
    # this is used for both to create random sample
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
# for loop over the number of runs you want to do
   # print(len(geneslist))
    for runs in range(numruns):
    # randomly select a sample to start with
        state = random.sample(geneslist, setsize)
        state = set(state)
    # for loop for number of iterations
        for ii in range(numiter):
    
        # randomly select a new gene to swap
            geneswapin = random.sample(geneslist,1)
            # cannot do this right away, what if the swapin is in the gene set
            # but the swap out is a different gene than that one?
            # doing the swap functions below would then reduce the set size to 2 or 1
           # geneswapout = random.sample(state,1)
            if geneswapin[0] in state:
               geneswapout = geneswapin
            else:
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
        # update solutions, amount of times a soltuon has been visited
            if (ii+1) % (stepsize)==0:
                frzstate = frozenset(state)
                if frzstate not in numvisits:
                    numvisits[frzstate] = 0
                numvisits[frzstate] += 1
        
    # end for loop for iterations
    
    #OUTPUT TO FILE - NOTE, SINCE WE ARE TRYING TO REPLICATE RESULTS IN VANDIN ET AL,
    # I HAVE THE OUTPUT IN THE SAME FORMAT AS THEIR OUTPUT FILES
    # this included only outputting the top 1000 for most-frequent and highest-weight files
        print('Done simulation')
    # find the most-visited sets of genes

        # for ease of sorting
        sorted_visits = list()
        sorted_weights = list()
        for geneset in numvisits:
            sorted_visits.append([numvisits[geneset], geneset])
            
        sorted_visits.sort()
    # write to file
        out_mostvisit = open('highestfreq_genesets_run'+str(runs)+'.txt','w')
        out_mostvisit.write('Top 1000 sets with highest frequency seen during simulation, number_visits set_of_genes state_weight \n')
        for i in range(len(sorted_visits)):
#             if i < 1000:
#                out_mostvisit.write([str(sorted_visits[(-i-1)][0])+'\t']) 
                
             gene_sort = list(sorted_visits[(-i-1)][1])
             # sorting the genes in alphabetical order like their output
             gene_sort.sort()
             gene_str = ''
             for gene in gene_sort:
                    gene_str = gene_str+str(gene)+'\t'

             # creating weight function value for each set, to be used for sorting by weight later                
             stateweight = weight(gene_sort)
             if i < 1000:
                 # frequency
                 out_mostvisit.write(str(sorted_visits[(-i-1)][0])+'\t')
                 # genes
                 out_mostvisit.write(gene_str)
                 # weight
                 out_mostvisit.write(str(stateweight)+'\n')
                 
             sorted_weights.append([stateweight, gene_str, str(sorted_visits[(-i-1)][0])])
             
        out_mostvisit.close()

    # find the sets visited with the highest weight seen
        sorted_weights.sort()
        out_mostweight = open('highestweight_genesets_run'+str(runs)+'.txt','w')
        out_mostweight.write('Top 1000 sets with highest weight seen during simulation, state_weights set_of_genes number_visits \n')
       
        # can we assume, for purposes of our project, we have at least 1000 unique sets seen?
        # not for initial testing at least
        for j in range(min(len(numvisits),1000)):
    # write to file
            out_mostweight.write(str(sorted_weights[(-j-1)][0])+'\t')
            out_mostweight.write(sorted_weights[(-j-1)][1])
            out_mostweight.write(sorted_weights[(-j-1)][2]+'\n')
        out_mostweight.close()
        print('Done Writing')
        


