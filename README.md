# Project 1b: Mutually Exclusive Mutations

For this project, you will implement and run the Markov chain Monte Carlo algorithm from [Vandin, et al. (Genome Research, 2012)](http://genome.cshlp.org/content/22/2/375.full).

### Data

The input data for your algorithm is a binary mutation matrix. The mutation matrices are stored in a tab-separated text file, where each line lists the genes with mutations in a single patient. The patient name will be stored in the first column, with subsequent columns containing gene names.

#### Example data

You can find a small example dataset for your project in [data/examples](https://github.com/cmsc828p-f17/project1b-mutually-exclusive-mutations/blob/master/data/examples). When run on this example dataset, your program should identify the implanted set of four "genes" (Gene-1, Gene-2, Gene-3, Gene-4) with significantly mutually exclusive mutations.

#### Real data

You will need to download real data for your project. You can find the lung cancer data originally used in the Vandin, et al. (Genome Research, 2012) study at [http://compbio.cs.brown.edu/projects/dendrix/](http://compbio.cs.brown.edu/projects/dendrix/).
