#include"DE.hpp"

//cec2013 stuff
double *OShift,*M,**y,**z,*x_bound;
int ini_flag=0,n_flag,func_flag;
//cec2013 stuff end

random_device rd;
mt19937_64 mt(rd());

//Main functions
Population *initializePop(FILE *input){
    Population *pop = NULL;
    Individual *ind = NULL;
    Individual *trial = NULL;
    
    pop = new Population();
    parametersInput(pop, input);
    pop->ind.reserve(pop->np);
    pop->trial.reserve(pop->np);
    
    uniform_real_distribution<double> dist(pop->li, pop->ui);
    
    for(int i=0;i<pop->np;i++){
        ind = new Individual();
        trial = new Individual();
        
        ind->fit = 10000.0;
        trial->fit = 10000.0;
        ind->data.reserve(pop->d);
        trial->data.reserve(pop->d);
        for(int j=0;j<pop->d;j++){
            ind->data.push_back(dist(mt));
            trial->data.push_back(10000.0);
        }
        
        
        pop->ind.push_back(ind);
        pop->trial.push_back(trial);
    }
    
    return pop;
}

void mainPopulationFitness(Population *pop, double (*p) (Individual*, int)){
    #pragma omp parallel for schedule(dynamic) num_threads(3)
    for(int i=0;i<pop->np;i++){
        pop->ind[i]->fit = (*p)(pop->ind[i], pop->d);
    }
}

void trialPopulationFitness(Population *pop, double (*p) (Individual*, int)){
    #pragma omp parallel for schedule(dynamic) num_threads(3)
    for(int i=0;i<pop->np;i++){
        pop->trial[i]->fit = (*p)(pop->trial[i], pop->d);
    }
}

void generateTrialPop(Population *pop){
    
    uniform_int_distribution<int> dist_mut(0, pop->np - 1);
    uniform_int_distribution<int> dist_swap(0, pop->d - 1);
    uniform_real_distribution<double> dist_cross(0.0, 1.0);
    int ind_mut[3];
    int swap;
    double swap_temp;
    float cross;
    
    
    for(int i=0;i<pop->np;i++){
        // Select individuals for mutation
        // De/Rand/1 in this case
        do{
            ind_mut[0] = dist_mut(mt);
        } while(ind_mut[0] == i);
        do{
            ind_mut[1] = dist_mut(mt);
        } while(ind_mut[1] == i || ind_mut[1] == ind_mut[0]);
        do{
            ind_mut[2] = dist_mut(mt);
        } while(ind_mut[2] == i || ind_mut[2] == ind_mut[0] || ind_mut[2] == ind_mut[1]);
        
        // Call the desired mutation strategy
        mutation(pop, i, ind_mut, mutationRand1);
        
        // Crossover Operation
        swap = dist_swap(mt);
        swap_temp = pop->trial[i]->data[swap];
        for(int j=0;j<pop->d;j++){
            cross = dist_cross(mt);
            if(cross > pop->cr){
                pop->trial[i]->data[j] = pop->ind[i]->data[j];
            }
        }
        pop->trial[i]->data[swap] = swap_temp;
    }
    
}

void mutation(Population *pop, int ind_parent, int *ind_mut, void (*p) (Population*, int, int*)){
    (*p)(pop, ind_parent, ind_mut);
}

void mutationRand1(Population *pop, int ind_parent, int *ind_mut){
    for(int i=0;i<pop->d;i++){
        //Mutation strategy DE/rand/1
        pop->trial[ind_parent]->data[i] = pop->ind[ind_mut[0]]->data[i] + 
                                          pop->f*(pop->ind[ind_mut[1]]->data[i] - pop->ind[ind_mut[2]]->data[i]);
        
        //Adjust limits of trial individual, if needed
        if(pop->trial[ind_parent]->data[i] < pop->li){
            pop->trial[ind_parent]->data[i] = (pop->li + pop->ind[ind_parent]->data[i])/2.0;
        }
        else if(pop->trial[ind_parent]->data[i] > pop->ui){
            pop->trial[ind_parent]->data[i] = (pop->ui + pop->ind[ind_parent]->data[i])/2.0;
        }
    }
}

void selectionMin(Population *pop){
    for(int i=0;i<pop->np;i++){
        if(pop->trial[i]->fit < pop->ind[i]->fit){
            for(int j=0;j<pop->d;j++){
                pop->ind[i]->data[j] = pop->trial[i]->data[j];
            }
            pop->ind[i]->fit = pop->trial[i]->fit;
        }
    }
}

void cec2013Fitness(Population *pop, bool mainPop, int func, int nThreads){
    double *f = (double*)malloc(pop->np * sizeof(double));
    vector<Individual*> *p;
    int optimum;
    
    if(mainPop){
        p = &(pop->ind);
    }
    else{
        p = &(pop->trial);
    }
    
    #pragma omp parallel for schedule(dynamic) num_threads(nThreads)
    for(int i=0;i < pop->np;i++){
        int id = omp_get_thread_num();
        test_func(&((*p)[i]->data[0]), &f[i], pop->d, 1, func, id);
    }
    
    optimum = getGlobalOptimum(func);
    for(int i = 0;i < pop->np;i++){
        (*p)[i]->fit = f[i] - optimum;
    }
    
    free(f);
}

int getGlobalOptimum(int func){
    int optimum;
    switch(func) {
        case 1 :
            optimum = -1400;
            break;
        case 2 :
            optimum = -1300;
            break;
        case 3 :
            optimum = -1200;
            break;
        case 4 :
            optimum = -1100;
            break;
        case 5 :
            optimum = -1000;
            break;
        case 6 :
            optimum = -900;
            break;
        case 7 :
            optimum = -800;
            break;
        case 8 :
            optimum = -700;
            break;
        case 9 :
            optimum = -600;
            break;
        case 10 :
            optimum = -500;
            break;
        case 11 :
            optimum = -400;
            break;
        case 12 :
            optimum = -300;
            break;
        case 13 :
            optimum = -200;
            break;
        case 14 :
            optimum = -100;
            break;
        case 15 :
            optimum = 100;
            break;
        case 16 :
            optimum = 200;
            break;
        case 17 :
            optimum = 300;
            break;
        case 18 :
            optimum = 400;
            break;
        case 19 :
            optimum = 500;
            break;
        case 20 :
            optimum = 600;
            break;
        case 21 :
            optimum = 700;
            break;
        case 22 :
            optimum = 800;
            break;
        case 23 :
            optimum = 900;
            break;
        case 24 :
            optimum = 1000;
            break;
        case 25 :
            optimum = 1100;
            break;
        case 26 :
            optimum = 1200;
            break;
        case 27 :
            optimum = 1300;
            break;
        case 28 :
            optimum = 1400;
            break;
      }
      return optimum;
}
//Main functions end

//Aux functions
double wtime(){
   struct timeval t;
   gettimeofday(&t, NULL);
   return t.tv_sec + t.tv_usec / 1000000.0;
}
//Aux functions end

//Clean-up fuctions
void destroyIndividual(Individual *ind){
    delete ind;
}

void destroyPopulation(Population *pop){
    for(int i=0;i<pop->np;i++){
        destroyIndividual(pop->ind[i]);
        destroyIndividual(pop->trial[i]);
    }
    delete pop;
}
//Clean-up functions end

//Input function
void parametersInput(Population *pop, FILE *arq){
    int np, d, aval, r;
    double li, ui, f, cr;
    int temp;
    
    temp = fscanf(arq, "np: %d\n", &np);
    temp = fscanf(arq, "d: %d\n", &d);
    temp = fscanf(arq, "aval: %d\n", &aval);
    temp = fscanf(arq, "li: %lf\n", &li);
    temp = fscanf(arq, "ui: %lf\n", &ui);
    temp = fscanf(arq, "f: %lf\n", &f);
    temp = fscanf(arq, "cr: %lf\n", &cr);
    temp = fscanf(arq, "r: %d\n", &r);
    
    rewind(arq);
    
    pop->np = np;
    pop->d = d;
    pop->aval = aval*d;
    pop->li = li;
    pop->ui = ui;
    pop->f = f;
    pop->cr = cr;
}
//Input function end

//Output functions
double meanFitness(Population *pop){
    double av = 0.0;
    
    for(int i=0;i<pop->np;i++){
        av += pop->ind[i]->fit;
    }
    
    return (av/(double)pop->np);
}

double standardDeviationFitness(Population *pop, double mean){
    double sd = 0.0;
    
    for(int i=0;i<pop->np;i++){
        sd += pow((pop->ind[i]->fit - mean), 2.0);
    }
    
    return sqrt(sd/(double)pop->np);
}

double standardDeviationRuns(vector<vector<double>> &meanRuns, double mean, int runs, int gen){
    double sd = 0.0;
    
    for(int r=0;r<runs;r++){
        sd += pow((meanRuns[r][gen] - mean), 2.0);
    }
    
    return sqrt(sd/(double)runs);
}

double meanArray(double *arr, int size){
    double mean = 0.0;
    
    for(int i=0;i<size;i++){
        mean += arr[i];
    }
    
    return (mean/(double)size);
}

double standardDeviationArray(double *arr, int size){
    double mean = meanArray(arr, size);
    double sd = 0.0;
    
    for(int i=0;i<size;i++){
        sd += pow((arr[i] - mean), 2.0);
    }
    
    return sqrt(sd/(double)size);
}

void saveFitness(Population *pop, int ger, FILE *arq){
    double m = meanFitness(pop);
    double sd = standardDeviationFitness(pop, m);
    fprintf(arq, "%8d %20.15f %20.15f %20.15f\n", ger, m, m+sd, m-sd);
}

//Population diversity as momentum of inertia, with modifications
double populationDiversity(Population *pop){
    double *cm = NULL;
    double aux;
    double div = 0.0;
    
    cm = (double*)malloc(pop->d*sizeof(double));
    
    //Compute certer of mass
    for(int i=0;i<pop->d;i++){
        aux = 0.0;
        for(int j=0;j<pop->np;j++){
            aux += pop->ind[j]->data[i];
        }
        aux = aux/(double)pop->np;
        cm[i] = aux;
    }
    
    //Compute diversity
    for(int i=0;i<pop->d;i++){
        aux = 0.0;
        for(int j=0;j<pop->np;j++){
            aux += pow((pop->ind[j]->data[i] - cm[i]), 2.0);
        }
        aux = sqrt(aux/(double)(pop->np - 1));
        div += aux;
    }
    free(cm);
    
    return (div/(double)pop->d);
}

void saveDiversity(Population *pop, int ger, FILE *arq){
    double d = populationDiversity(pop);
    fprintf(arq, "%8d %20.15f\n", ger, d);
}

int bestFitnessIndex(Population *pop){
    int index = 0;
    double best = pop->ind[0]->fit;
    
    for(int i=1;i<pop->np;i++){
        if(pop->ind[i]->fit < best){
            index = i;
            best = pop->ind[i]->fit;
        }
    }
    
    return index;
}
//Output functions end

//Visualization Functions
void printDataIndividual(Individual *ind, int d){
    printf("[");
    for(int i=0;i<d - 1;i++){
        printf("%.5f, ", ind->data[i]);
    }
    printf("%.5f]\n", ind->data[d - 1]);
}

void printDataPop(Population *pop){
    for(int i=0;i<pop->np;i++){
        printf("%d: ", i);
        printDataIndividual(pop->ind[i], pop->d);
    }
}

void printDataTrial(Population *pop){
    for(int i=0;i<pop->np;i++){
        printf("%d: ", i);
        printDataIndividual(pop->trial[i], pop->d);
    }
}

void printFitIndividual(Individual *ind){
    printf("%.20f\n", ind->fit);
}

void printFitPop(Population *pop){
    for(int i=0;i<pop->np;i++){
        printf("%d: ", i);
        printFitIndividual(pop->ind[i]);
    }
}
//Visualization Functions end


