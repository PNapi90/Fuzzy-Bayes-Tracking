#include "perm_class.h"


#ifndef MERGE_CLASS_H
#define MERGE_CLASS_H

class merge_class{

private:

    int am_clust, full_len, am_perm;
    int interactions;
    int* lens;
    int* cluster_tmp;

    int** ret_perms_arr;

    double* cluster_probs;

    double* ret_arr;

    perm_class** A;

public:
    merge_class(int*,int,int,perm_class**,int);
    ~merge_class();
    
    void compare_true();
    void set_probability(double,int);
    void reset_probability();

    int get_cluster_width(int);

    int* return_sub_cluster_perm(int,int);
    double* get_highest_prob_cluster();
};




#endif