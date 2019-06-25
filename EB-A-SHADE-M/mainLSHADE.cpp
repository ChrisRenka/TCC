#include"LSHADE.hpp"

int main(int argc, char **argv){
    
    Population *pop = NULL;
    FILE *input = NULL;
    FILE *fitFile = NULL;
    FILE *divFile = NULL;
    FILE *output = NULL;
    FILE *mutProbFile = NULL;
    char nameFitFile[50];
    char nameDivFile[50];
    char nameOutput[50];
    char nameMutProbFile[50];
    int bestIndividual;
    int runs;
    int gens = 0;
    int func;
    int nThreads;
    int currentAval;
    int currentGen;
    double bestFitness;
    double meanFitGen;
    double meanFitRuns;
    double *runsBestFitness;
    double sdBestFitness;
    double start_time, end_time;
    bool alloc = false;
    
    if ((argc != 3)) {
        printf("To run: %s <func_num> <runs>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
   
    func = atoi(argv[1]);
    runs = atoi(argv[2]);
    
    vector<vector<double>> averageFitness(runs);
    vector<double> averageBestFitness;
    vector<double> sdFitness;
    vector<double> averageDiversity;
    vector<vector<double>> averageMutProb;
    
    input = fopen("input.txt", "r");
    
    nThreads = thread::hardware_concurrency() - 1;
    printf("nThreads: %d\n", nThreads);
    
    start_time = wtime();
    
    for(int r=0;r<runs;r++){
    
        pop = initializePop(input);
        
        if(!alloc){
            sprintf(nameFitFile, "tests/L-SHADE_f%d_D%d_fitness.txt", func, pop->d);
            sprintf(nameDivFile, "tests/L-SHADE_f%d_D%d_diversity.txt", func, pop->d);
            sprintf(nameOutput, "tests/L-SHADE_f%d_D%d_output.txt", func, pop->d);
            sprintf(nameMutProbFile, "tests/L-SHADE_f%d_D%d_mutProb.txt", func, pop->d);
    
            fitFile = fopen(nameFitFile, "w+");
            divFile = fopen(nameDivFile, "w+");
            output = fopen(nameOutput, "w+");
            mutProbFile = fopen(nameMutProbFile, "w+");
            
            fprintf(output, "L-SHADE, function number %d, %d dimensions\n", func, pop->d);
            fprintf(output, "average of %d runs\n", runs);
            fprintf(output, "\n");
            fprintf(output, "Best fitness for each run:\n");
            
            for(int i=0;i<runs;i++){
                averageFitness[i].reserve((double)pop->aval/(double)pop->np);
            }
            averageBestFitness.reserve((double)pop->aval/(double)pop->np);
            sdFitness.reserve((double)pop->aval/(double)pop->np);
            averageDiversity.reserve((double)pop->aval/(double)pop->np);
            
            runsBestFitness = (double*)malloc(runs*sizeof(double));
            
            averageMutProb.reserve(pop->h);
            for(int i=0;i<pop->h;i++){
                averageMutProb[i].reserve((double)pop->aval/(double)pop->np);
            }
            
            int np = pop->np;
            double div = (double)pop->npf/(double)pop->npi;
            currentAval = np;
            while(currentAval + np < pop->aval){
                gens += 1;
                currentAval += np;
                np = round(pop->npi*(pow(div, (currentAval/(double)pop->aval))));
                if(np < pop->npf){
                    np = pop->npf;
                }
            }
            
            alloc = true;
            setupCEC2013(pop->d, func, nThreads);
        }
        
        cec2013Fitness(pop, true, func, nThreads);
        meanFitGen = meanFitness(pop);
        
        bestIndividual = bestFitnessIndex(pop);
        averageFitness[r].push_back(meanFitGen);
        if(r==0){
            averageBestFitness.push_back(pop->ind[bestIndividual]->fit);
            averageDiversity.push_back(populationDiversity(pop));
            for(int i=0;i<pop->h;i++){
                averageMutProb[i].push_back(pop->mfcp[i]);
            }
        }
        else{
            averageBestFitness[0] += pop->ind[bestIndividual]->fit;
            averageDiversity[0] += populationDiversity(pop);
            for(int i=0;i<pop->h;i++){
                averageMutProb[i][0] += pop->mfcp[i];
            }
        }
        
        currentAval = pop->np;
        currentGen = 1;
        while(currentAval + pop->np < pop->aval){
            generateTrialPop(pop, gens, currentGen);
            
            cec2013Fitness(pop, false, func, nThreads);
            
            selectionMin(pop);
            meanFitGen = meanFitness(pop);
            
            bestIndividual = bestFitnessIndex(pop);
            averageFitness[r].push_back(meanFitGen);
            if(r==0){
                averageBestFitness.push_back(pop->ind[bestIndividual]->fit);
                averageDiversity.push_back(populationDiversity(pop));
                for(int i=0;i<pop->h;i++){
                    averageMutProb[i].push_back(pop->mfcp[i]);
                }
            }
            else{
                averageBestFitness[currentGen] += pop->ind[bestIndividual]->fit;
                averageDiversity[currentGen] += populationDiversity(pop);
                for(int i=0;i<pop->h;i++){
                    averageMutProb[i][currentGen] += pop->mfcp[i];
                }
            }
            
            currentGen += 1;
            currentAval += pop->np;
            
            updateMeans(pop);
            alternativePopulationReduction(pop, currentAval);
            adjustSettings(pop);
        }
        if(r == 0){
            bestFitness = pop->ind[bestIndividual]->fit;
        }
        else{
            if(pop->ind[bestIndividual]->fit < bestFitness){
                bestFitness = pop->ind[bestIndividual]->fit;
            }
        }
        
        runsBestFitness[r] = pop->ind[bestIndividual]->fit;
        printf("Best fitness run %d: %25.15e\n", r+1, pop->ind[bestIndividual]->fit);
        fprintf(output, "%25.15e\n", pop->ind[bestIndividual]->fit);
        
        destroyPopulation(pop);
    }
    end_time = wtime();
    
    fprintf(output, "\nMean fitness for each run:\n");
    for(int i=0;i<runs;i++){
        fprintf(output, "%25.15e\n", averageFitness[i][gens-1]);
    }
    fprintf(output, "\n\n");
    for(int i=0;i<gens;i++){
        meanFitRuns = 0.0;
        for(int j=0;j<runs;j++){
            meanFitRuns += averageFitness[j][i];
        }
        meanFitRuns = meanFitRuns/(double)runs;
        sdFitness.push_back(standardDeviationRuns(averageFitness, meanFitRuns, runs, i));
        
        fprintf(fitFile, "%10d %25.15e %25.15e %25.15e %25.15e\n", i, averageBestFitness[i]/runs, meanFitRuns, 
                                                                  (meanFitRuns+sdFitness[i]), (meanFitRuns-sdFitness[i]));
        fprintf(divFile, "%10d %25.15e\n", i, averageDiversity[i]/runs);
        fprintf(mutProbFile, "%10d ", i);
        for(int j=0;j<pop->h;j++){
            fprintf(mutProbFile, "%12.6f ", averageMutProb[j][i]/runs);
        }
        fprintf(mutProbFile, "\n");
    }
    
    sdBestFitness = standardDeviationArray(runsBestFitness, runs);
    
    printf("Best fitness overall: %25.15e\n", bestFitness);
    printf("Average best fitness: %25.15e +- %25.15e\n", averageBestFitness[gens-1]/runs, sdBestFitness);
    printf("Average fitness: %25.15e +- %25.15e\n", meanFitRuns, sdFitness[gens-1]);
    printf("Total time: %8.3f s\n", end_time - start_time);
    
    fprintf(output, "Best fitness overall: %25.15e\n", bestFitness);
    fprintf(output, "Average best fitness: %25.15e +- %25.15e\n", averageBestFitness[gens-1]/runs, sdBestFitness);
    fprintf(output, "Average fitness: %25.15e +- %25.15e\n", meanFitRuns, sdFitness[gens-1]);
    fprintf(output, "Total time: %8.3f s\n", end_time - start_time);
    
    fclose(input);
    fclose(fitFile);
    fclose(divFile);
    fclose(output);
    fclose(mutProbFile);
    
    free(runsBestFitness);
    
    return 0;
    
}
