#include <iostream>

#include "merge_class.h"

using namespace std;

merge_class::merge_class(int* lens,int am_clust,int am_perm,perm_class** A,int interactions){
    this->A = A;
    this->am_clust = am_clust;
    this->lens = lens;
    this->am_perm = am_perm;
    this->interactions = interactions;

    cluster_probs = new double[am_perm];
    cluster_tmp = new int[interactions]; 
    full_len = 0;

    ret_perms_arr = new int*[am_clust];
    for(int i = 0;i < am_clust;i++) ret_perms_arr[i] = new int[this->lens[i]];

    for(int i = 0;i < this->am_clust;i++) full_len += this->lens[i];

    ret_arr = new double[2];
}

merge_class::~merge_class(){
    
    for(int i = 0;i < am_clust;i++) delete[] ret_perms_arr[i];
    delete[] ret_perms_arr;
    delete[] cluster_probs;
    delete[] cluster_tmp;
    delete[] ret_arr;
}


void merge_class::compare_true(){
    cout << "am_clust " << am_clust << endl;
    for(int i = 0;i < am_perm;i++){
            for(int j = 0;j < am_clust;j++){
                A[j]->cout_line(i);
                cout << "| ";
            }
        cout << endl;
    }
}

void merge_class::set_probability(double prob,int pos){cluster_probs[pos] = prob;}

void merge_class::reset_probability(){for(int i = 0;i < am_perm;i++) cluster_probs[i] = 0;}

int* merge_class::return_sub_cluster_perm(int k,int j){
    for(int i = 0;i < lens[k];i++) ret_perms_arr[k][i] = A[k]->return_line(j,i);
    return ret_perms_arr[k];
}


int merge_class::get_cluster_width(int k){return A[k]->get_len_x();}

double* merge_class::get_highest_prob_cluster(){
    double max_prob = 0;
    double sum_prob = 0;
    double pos = 0;
    for(int i = 0;i < am_perm;i++){
        sum_prob += cluster_probs[i];
        if(max_prob < cluster_probs[i]){
            max_prob = cluster_probs[i];
            pos = (double) i;
        }
    }
    ret_arr[0] = max_prob;///sum_prob;
    ret_arr[1] = pos;

    return ret_arr;
}
