# TCC
Self adaptive Differential Evolution codes for my undergraduate thesis in Computer Science (TCC)

Each folder contains the codes for an algorithm, with the exception of the folder "Fitness", which contains the codes for the benchmark functions.
The codes are using the CEC-2013 benchmark set (28 functions). The CEC codes were modified to allow the simultaneous evaluation of multiple individuals. The codes are running with OpenMP, and the number of threads used are equal to the maximum number supported by the system minus one.
The codes were written and run in a Ubuntu Linux environment and in a Windows Subsystem for Linux (also Ubuntu).

**Compilation**

Use the command "make" in a terminal to compile de codes for an algorithm.
The output file will be called "main".

**Running**

Run with "./main <function number> <number of runs>".
For example, "./main 6 30" will execute 30 independent runs of function 6, and the outputs will give the means of the 30 runs.

**Parameters input**

Each algorithm folder has an "input.txt" file, which holds the parameters. The complete list of parameters is:

- np: population size;
- npi(x d): initial population size to be multiplied by the dimensionality (L and A-SHADE algorithms);
- npf: final population size (L and A-SHADE algorithms);
- d: problem dimensionality;
- aval: number of function evaluations divided by the dimensionality;
- li: lower bound for the value of each gene;
- ui: upper bound for the value of each gene;
- f: mutation factor (only needed for the canonical DE);
- cr: crossover rate (only needed for the canonical DE);
- h: memory size (for SHADE based algorithms);
- p: percentage of best individuals used in the pbest mutation;
- rArc: this multiplied by the current population size will be the maximum size of the external archive (L and A-SHADE algorithms);
- c: learning rate used by EB algorithms;
- mut\_sel: probability for a target vector to undergo non-uniform mutation (EB-A-SHADE-M only);
- mut\_p: probability for a given gene in the target vector to undergo non-uniform mutation (EB-A-SHADE-M only);
- mut\_b: value of 'b' used in the non-uniform mutation formula (EB-A-SHADE-M only). It is kept at 5 by default.

Only the parameters used by a given algorithm will be in their input file.

**Output**

The execution outputs will be saved in the folder "tests". All algorithms will create, for each combination of function and dimensionality, an output file, a fitness file, and a diversity file. EB algorithms will also output a mutProb file.
The contents of these files are:

- output: this file will contain the best and mean result for each run, as well as the best result overall, means for both the best and mean results, and the total execution time.
- fitness: this file has five columns. They are, in order:
    - generation number.
    - best fitness averaged over all runs.
    - mean fitness averaged over all runs.
    - mean fitness averaged over all runs + standard deviation.
    - mean fitness averaged over all runs - standard deviation.
- diversity: this file has two columns. They are, in order:
    - generation number.
    - population diversity averaged over all runs.
- mutProb: this file has one plus h columns. The first one is the generation number, and each other one is entry for the list of mutation probabilities used in the EB algorithms, averaged over all runs.