#include"LSHADE.hpp"

//cec2013 stuff
double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;
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
    pop->aux.reserve((pop->np)*(pop->rArc)*2);
    
    pop->sf.reserve(pop->np);
    pop->scr.reserve(pop->np);
    pop->fitDiff.reserve(pop->np);
    
    pop->mcr.reserve(pop->h);
    pop->mf.reserve(pop->h);
    pop->k = 0;
    
    for(int i=0;i<pop->h;i++){
        pop->mf.push_back(0.5);
        pop->mcr.push_back(0.5);
        pop->mfcp.push_back(0.5);
    }
    
    uniform_real_distribution<double> dist(pop->li, pop->ui);
    
    
    for(int i=0;i<pop->np;i++){
        ind = new Individual();
        trial = new Individual();
        
        ind->id = i;
       
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
    uniform_int_distribution<unsigned int> dist_mut2(0, pop->np + pop->aux.size() - 1);
    uniform_int_distribution<int> dist_swap(0, pop->d - 1);
    uniform_int_distribution<int> dist_h(0, pop->h - 1);
    uniform_real_distribution<double> dist_cross(0.0, 1.0);
    int ind_mut[3];
    int swap;
    double swap_temp;
    double mut_choice;
    float cross;
    vector<int> pBest = findPBest(pop);
    int indPBest;
    int pSize;
    int ri;
    int temp;
    
    for(int i=0;i<pop->np;i++){
        ri = dist_h(mt);
        normal_distribution<double> dist_cr(pop->mcr[ri], 0.1);
        cauchy_distribution<double> dist_f(pop->mf[ri], 0.1);
        
        pop->ind[i]->cr = dist_cr(mt);
        if(pop->ind[i]->cr < 0.0){
            pop->ind[i]->cr = 0.0;
        }
        else if(pop->ind[i]->cr > 1.0){
            pop->ind[i]->cr = 1.0;
        }
        
        do{
            pop->ind[i]->f = dist_f(mt);
        } while(pop->ind[i]->f <= 0.0);
        if(pop->ind[i]->f > 1.0){
            pop->ind[i]->f = 1.0;
        }
        
        pSize = round((pop->p)*(pop->np));
        if(pSize<2){
            pSize = 2;
        }
        
        // Select individuals for mutation
        // De/current-to-pbest/1 in this case
        uniform_int_distribution<int> dist_pbest(0, pSize - 1);
        
        indPBest = dist_pbest(mt);
        ind_mut[0] = pBest[indPBest];
        do{
            ind_mut[1] = dist_mut(mt);
        } while(ind_mut[1] == i);
        do{
            ind_mut[2] = dist_mut2(mt);
        } while(ind_mut[2] == i || ind_mut[2] == ind_mut[1]);
        
        mut_choice = dist_cross(mt);
        
        // Call the desired mutation strategy
        if(mut_choice <= pop->mfcp[ri]){
            mutation(pop, i, ind_mut, mutationCurrentToPBest);
            pop->trial[i]->mut = 1;
        }
        else{
            //Sort selected individuals
            if(pop->ind[ind_mut[0]]->fit > pop->ind[ind_mut[1]]->fit){
                temp = ind_mut[0];
                ind_mut[0] = ind_mut[1];
                ind_mut[1] = temp;
            }
            if(ind_mut[2] < pop->np){
                if(pop->ind[ind_mut[1]]->fit > pop->ind[ind_mut[2]]->fit){
                    temp = ind_mut[1];
                    ind_mut[1] = ind_mut[2];
                    ind_mut[2] = temp;
                    if(pop->ind[ind_mut[0]]->fit > pop->ind[ind_mut[1]]->fit){
                        temp = ind_mut[0];
                        ind_mut[0] = ind_mut[1];
                        ind_mut[1] = temp;
                    }
                }
            }
            else{
                if(pop->ind[ind_mut[1]]->fit > pop->aux[ind_mut[2] - pop->np]->fit){
                    temp = ind_mut[1];
                    ind_mut[1] = ind_mut[2];
                    ind_mut[2] = temp;
                    if(pop->ind[ind_mut[0]]->fit > pop->aux[ind_mut[1] - pop->np]->fit){
                        temp = ind_mut[0];
                        ind_mut[0] = ind_mut[1];
                        ind_mut[1] = temp;
                    }
                }
            }
            
            mutation(pop, i, ind_mut, mutationCurrentToOrdPBest);
            pop->trial[i]->mut = 2;
        }
        
        
        // Crossover Operation
        swap = dist_swap(mt);
        swap_temp = pop->trial[i]->data[swap];
        for(int j=0;j<pop->d;j++){
            cross = dist_cross(mt);
            if(cross > pop->ind[i]->cr){
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
                                          pop->ind[ind_parent]->f*(pop->ind[ind_mut[1]]->data[i] - pop->ind[ind_mut[2]]->data[i]);
        
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
    pop->sumDiffM1 = 0.0;
    pop->sumDiffM2 = 0.0;
    for(int i=0;i<pop->np;i++){
        if(pop->trial[i]->fit == pop->ind[i]->fit){
            for(int j=0;j<pop->d;j++){
                pop->ind[i]->data[j] = pop->trial[i]->data[j];
            }
            pop->ind[i]->mut = pop->trial[i]->mut;
        }
        else if(pop->trial[i]->fit < pop->ind[i]->fit){
            //Comment the next line to remove the archive
            if(pop->trial[i]->mut == 1){
                pop->sumDiffM1 += (pop->ind[i]->fit - pop->trial[i]->fit);
            }
            else{
                pop->sumDiffM2 += (pop->ind[i]->fit - pop->trial[i]->fit);
            }
            copyIndividualToAux(pop, i);
            for(int j=0;j<pop->d;j++){
                pop->ind[i]->data[j] = pop->trial[i]->data[j];
            }
            pop->ind[i]->mut = pop->trial[i]->mut;
            pop->scr.push_back(pop->ind[i]->cr);
            pop->sf.push_back(pop->ind[i]->f);
            pop->fitDiff.push_back(fabs(pop->ind[i]->fit - pop->trial[i]->fit));
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

//JADE/SHADE functions
vector<int> findPBest(Population *pop){
    vector<Individual*> sorted;
    vector<int> ret;
    int total = (int)((pop->np)*(pop->p)) + 1;
    
    if(total <= 2){
        total = 3;
    }
    
    sorted = pop->ind;
    sort(sorted.begin(), sorted.end(), compareByFitness);
    
    for(int i=0;i<total;i++){
        ret.push_back(sorted[i]->id);
    }
    
    return ret;
}

bool compareByFitness(Individual *a, Individual *b){
    return a->fit < b->fit;
}

double wArithmeticMean(vector<double> &v, vector<double> &d, double sumDiff){
    double sum = 0.0;
    double w;
    unsigned int n = v.size();
    
    
    for(unsigned int i=0;i<n;i++){
        w = d[i]/sumDiff;
        sum += w*v[i];
    }
    
    return sum;
}

double wLehmerMean(vector<double> &v, vector<double> &d, double sumDiff){
    double sum = 0.0;
    double sumSQ = 0.0;
    double w;
    unsigned int n = v.size();
    
    for(unsigned int i=0;i<n;i++){
        w = d[i]/sumDiff;
        sum += w*v[i];
        sumSQ += w*pow(v[i], 2.0);
    }
    
    return (sumSQ/sum);
}

void updateMeans(Population *pop){
    double meanCR;
    double meanF;
    double sumDiff = 0.0;
    double maxSCR = 0.0;
    double deltaMut1;
    int dSize = pop->fitDiff.size();
    
    /////////////////Testing first
    if(pop->scr.size() > 0){
        for(int i=0;i<dSize;i++){
            sumDiff += pop->fitDiff[i];
        }
        
        
        for(int i=0;i<pop->scr.size();i++){
            if(pop->scr[i] > maxSCR){
                maxSCR = pop->scr[i];
                break;
            }
        }
        
        if((pop->mcr[pop->k] == -1.0) || (maxSCR == 0.0)){
            pop->mcr[pop->k] = -1.0;
        }
        else{
            meanCR = wLehmerMean(pop->scr, pop->fitDiff, sumDiff);
            pop->mcr[pop->k] = meanCR;
        }
        
        meanF = wLehmerMean(pop->sf, pop->fitDiff, sumDiff);
        pop->mf[pop->k] = meanF;
        
        deltaMut1 = fmax(0.2, (pop->sumDiffM1/(pop->sumDiffM1 + pop->sumDiffM2)));
        deltaMut1 = fmin(0.8, deltaMut1);
        
        pop->mfcp[pop->k] = (1-pop->c)*pop->mfcp[pop->k] + pop->c*deltaMut1;
        
        pop->k++;
        if(pop->k >= pop->h){
            pop->k = 0;
        }
    }
}

void mutationCurrentToPBest(Population *pop, int ind_parent, int *ind_mut){
    for(int i=0;i<pop->d;i++){
        pop->trial[ind_parent]->data[i] = pop->ind[ind_parent]->data[i] + 
                                          pop->ind[ind_parent]->f*(pop->ind[ind_mut[0]]->data[i] - pop->ind[ind_parent]->data[i]);
        if(ind_mut[2] < pop->np){
            pop->trial[ind_parent]->data[i] += pop->ind[ind_parent]->f*(pop->ind[ind_mut[1]]->data[i] - pop->ind[ind_mut[2]]->data[i]);
        }
        else{
            pop->trial[ind_parent]->data[i] += pop->ind[ind_parent]->f*(pop->ind[ind_mut[1]]->data[i] - pop->aux[ind_mut[2] - pop->np]->data[i]);
        }
        
        //Adjust limits of trial individual, if needed
        if(pop->trial[ind_parent]->data[i] < pop->li){
            pop->trial[ind_parent]->data[i] = (pop->li + pop->ind[ind_parent]->data[i])/2.0;
        }
        else if(pop->trial[ind_parent]->data[i] > pop->ui){
            pop->trial[ind_parent]->data[i] = (pop->ui + pop->ind[ind_parent]->data[i])/2.0;
        }
    }
}

void mutationCurrentToOrdPBest(Population *pop, int ind_parent, int *ind_mut){
    for(int i=0;i<pop->d;i++){
        pop->trial[ind_parent]->data[i] = pop->ind[ind_parent]->data[i];
        
        if(ind_mut[0] < pop->np){
            pop->trial[ind_parent]->data[i] += pop->ind[ind_parent]->f*(pop->ind[ind_mut[0]]->data[i] - pop->ind[ind_parent]->data[i]);
        }
        else{
            pop->trial[ind_parent]->data[i] += pop->ind[ind_parent]->f*(pop->aux[ind_mut[0] - pop->np]->data[i] - pop->ind[ind_parent]->data[i]);
        }
        
        if( (ind_mut[1] < pop->np) && (ind_mut[2] < pop->np )){
            pop->trial[ind_parent]->data[i] += pop->ind[ind_parent]->f*(pop->ind[ind_mut[1]]->data[i] - pop->ind[ind_mut[2]]->data[i]);
        }
        else if(ind_mut[1] >= pop->np){
            pop->trial[ind_parent]->data[i] += pop->ind[ind_parent]->f*(pop->aux[ind_mut[1] - pop->np]->data[i] - pop->ind[ind_mut[2]]->data[i]);
        }
        else if(ind_mut[2] >= pop->np){
            pop->trial[ind_parent]->data[i] += pop->ind[ind_parent]->f*(pop->ind[ind_mut[1]]->data[i] - pop->aux[ind_mut[2] - pop->np]->data[i]);
        }
        else{
            printf("Mutation ord-pbest error\n");
        }
        
        //Adjust limits of trial individual, if needed
        if(pop->trial[ind_parent]->data[i] < pop->li){
            pop->trial[ind_parent]->data[i] = (pop->li + pop->ind[ind_parent]->data[i])/2.0;
        }
        else if(pop->trial[ind_parent]->data[i] > pop->ui){
            pop->trial[ind_parent]->data[i] = (pop->ui + pop->ind[ind_parent]->data[i])/2.0;
        }
    }
}

void copyIndividualToAux(Population *pop, int index){
    Individual *ind = new Individual();
    
    ind->fit = pop->ind[index]->fit;
    ind->data.reserve(pop->d);
    for(int i=0;i<pop->d;i++){
        ind->data.push_back(pop->ind[index]->data[i]);
    }
    
    pop->aux.push_back(ind);
}

void adjustSettings(Population *pop){
    //Clears vector sf and scr
    //Remove elements form vector aux, if needed
    
    pop->sf.clear();
    pop->sf.reserve(pop->np);
    
    pop->scr.clear();
    pop->scr.reserve(pop->np);
    
    pop->fitDiff.clear();
    pop->fitDiff.reserve(pop->np);
    
    if((int)pop->aux.size() > (int)((pop->np)*(pop->rArc))){
        int n = pop->aux.size() - (int)((pop->np)*(pop->rArc));
        int index;
        for(int i=0;i<n;i++){
            uniform_int_distribution<unsigned int> dist(0, pop->aux.size() - 1);
            index = dist(mt);
            delete pop->aux[index];
            pop->aux.erase(pop->aux.begin() + index);
        }
    }
    
}

void linearPopulationReduction(Population *pop, int currentAval){
    int nnp;
    int nReduction;
    int worst;
    
    nnp = round(((pop->npf - pop->npi)/(double)pop->aval)*currentAval + pop->npi);
    
    if(nnp < pop->npf){
        nnp = pop->npf;
    }
    
    if(nnp < pop->np){
        nReduction = pop->np - nnp;
        for(int i=0;i<nReduction;i++){
            worst = 0;
            for(int j=1;j<pop->np;j++){
                if((pop->ind[j]->fit) > (pop->ind[worst]->fit)){
                    worst = j;
                }
            }
            
            delete pop->ind[worst];
            pop->ind.erase(pop->ind.begin() + worst);
            
            delete pop->trial[worst];
            pop->trial.erase(pop->trial.begin() + worst);
            
            pop->np = pop->np - 1;
            
        }
        
        for(int i=0;i<pop->np;i++){
            pop->ind[i]->id = i;
        }
        
    }
}

void alternativePopulationReduction(Population *pop, int currentAval){
    int nnp;
    int nReduction;
    int worst;
    double div = (double)pop->npf/(double)pop->npi;
    
    nnp = round(pop->npi*(pow(div, (currentAval/(double)pop->aval))));
    
    if(nnp < pop->npf){
        nnp = pop->npf;
    }
    
    if(nnp < pop->np){
        nReduction = pop->np - nnp;
        for(int i=0;i<nReduction;i++){
            worst = 0;
            for(int j=1;j<pop->np;j++){
                if((pop->ind[j]->fit) > (pop->ind[worst]->fit)){
                    worst = j;
                }
            }
            
            delete pop->ind[worst];
            pop->ind.erase(pop->ind.begin() + worst);
            
            delete pop->trial[worst];
            pop->trial.erase(pop->trial.begin() + worst);
            
            pop->np = pop->np - 1;
            
        }
        
        for(int i=0;i<pop->np;i++){
            pop->ind[i]->id = i;
        }
        
    }
}
//JADE/SHADE functions end

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
    int np, npi, npf, d, aval, h;
    double li, ui, c, p, rArc;
    int temp;
    
    temp = fscanf(arq, "npi(x d): %d\n", &npi);
    temp = fscanf(arq, "npf: %d\n", &npf);
    temp = fscanf(arq, "d: %d\n", &d);
    temp = fscanf(arq, "aval: %d\n", &aval);
    temp = fscanf(arq, "li: %lf\n", &li);
    temp = fscanf(arq, "ui: %lf\n", &ui);
    temp = fscanf(arq, "h: %d\n", &h);
    temp = fscanf(arq, "p: %lf\n", &p);
    temp = fscanf(arq, "rArc: %lf\n", &rArc);
    temp = fscanf(arq, "c: %lf\n", &c);
    
    rewind(arq);
    
    pop->npi = npi*d;
    pop->npf = npf;
    pop->np = npi*d;
    pop->d = d;
    pop->aval = aval*d;
    pop->li = li;
    pop->ui = ui;
    pop->h = h;
    pop->p = p;
    pop->rArc = rArc;
    pop->c = c;
}
//Input function end

//Auxiliary Functions
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

double wtime(){
   struct timeval t;
   gettimeofday(&t, NULL);
   return t.tv_sec + t.tv_usec / 1000000.0;
}
//Auxiliary Functions end

//Output functions

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
void printDataIndividual(Individual *ind){
    int d = ind->data.size();
    printf("[");
    for(int i=0;i<d - 1;i++){
        printf("%.5f, ", ind->data[i]);
    }
    printf("%.5f]\n", ind->data[d - 1]);
}

void printDataPop(Population *pop){
    for(int i=0;i<pop->np;i++){
        printf("%d: ", i);
        printDataIndividual(pop->ind[i]);
    }
}

void printDataTrial(Population *pop){
    for(int i=0;i<pop->np;i++){
        printf("%d: ", i);
        printDataIndividual(pop->trial[i]);
    }
}

void printFitIndividual(Individual *ind){
    printf("%.15f\n", ind->fit);
}

void printFitPop(Population *pop){
    for(int i=0;i<pop->np;i++){
        printf("%d: ", i);
        printFitIndividual(pop->ind[i]);
    }
}
//Visualization Functions end


