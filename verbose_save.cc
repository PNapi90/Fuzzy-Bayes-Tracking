#include "verbose_save.h"

using namespace std;

verbose_save::verbose_save(int N){
    this->N = N;
    integrals_save = new ofstream*[N];
    for(int i = 0;i < N;i++) integrals_save[i] = new ofstream("saved_integrals"+to_str(N-i)+".dat");

    integrals = new double**[N];
    for (int i = 0; i < N;i++){
        integrals[i] = new double*[N];
        for(int j = 0;j < N;j++){
            integrals[i][j] = new double[10000];
            for(int k = 0;k < 10000;k++) integrals[i][j][k] = 0.;
        }
    }
}

verbose_save::~verbose_save(){
    for(int i = 0;i < N;i++){
        integrals_save[i]->close();
        delete integrals_save[i];
        for(int j = 0;j < N;j++) delete integrals[i][j];
        delete[] integrals[i];
    }

    delete[] integrals;
    delete[] integrals_save;
}

void verbose_save::save_it(double* integral,int cluster_iter,int sub_cluster_id){
    for(int i = 0;i < 10000;i++) integrals[cluster_iter-1][sub_cluster_id][i] += integral[i];
}