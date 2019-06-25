#ifndef _DE_HPP
#define _DE_HPP

#include<stdio.h>
#include<omp.h>
#include<sys/time.h>
#include<random>
#include<vector>
#include<thread>
#include"../Fitness/test_func_parallel.hpp"

using namespace std;

typedef struct{
    vector<double> data;
    double fit;
} Individual;

typedef struct{
    vector<Individual*> ind;
    vector<Individual*> trial;
    int np;
    int d;
    int aval;
    double li;
    double ui;
    double cr;
    double f;
} Population;

//Main functions
Population *initializePop(FILE *input);
void mainPopulationFitness(Population *pop, double (*p) (Individual*, int));
void trialPopulationFitness(Population *pop, double (*p) (Individual*, int));
void generateTrialPop(Population *pop);
void mutation(Population *pop, int ind_parent, int *ind_mut, void (*p) (Population*, int, int*));
void mutationRand1(Population *pop, int ind_parent, int *ind_mut);
void selectionMin(Population *pop);
void cec2013Fitness(Population *pop, bool mainPop, int func, int nThreads);
int getGlobalOptimum(int func);
//Main functions end

//Aux functions
double wtime();
//Aux functions end

//Clean-up fuctions
void destroyIndividual(Individual *ind);
void destroyPopulation(Population *pop);
//Clean-up functions end

//Input function
void parametersInput(Population *pop, FILE *arq);
//Input function end

//Output functions
double meanFitness(Population *pop);
double standardDeviationFitness(Population *pop, double mean);
double standardDeviationRuns(vector<vector<double>> &meanRuns, double mean, int runs, int gen);
double meanArray(double *arr, int size);
double standardDeviationArray(double *arr, int size);
void saveFitness(Population *pop, int ger, FILE *arq);
double populationDiversity(Population *pop);
void saveDiversity(Population *pop, int ger, FILE *arq);
int bestFitnessIndex(Population *pop);
//Output functions end

//Visualization Functions
void printDataIndividual(Individual *ind, int d);
void printDataPop(Population *pop);
void printDataTrial(Population *pop);
void printFitIndividual(Individual *ind);
void printFitPop(Population *pop);
//Visualization Functions end

#endif
