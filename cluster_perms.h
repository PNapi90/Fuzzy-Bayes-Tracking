#include <iostream>

#include "permutations.h"
#include "perm_class.h"
#include "merge_class.h"


#ifndef CLUST_PERM_H
#define CLUST_PERM_H

class cluster_perms{

private:
    
    int active_cluster_int;
    int active_sub_cluster_int;
    int internal_perm_counter;
	int verbosity;
	int thr_num;
	
    double prob_tmp;

    permutations* loaded_perms;

    perm_class**** A;
    merge_class*** B;

    merge_class** active_cluster;
    merge_class* active_sub_cluster;

    int* remain_len_arr;
    int*** all_possibilities;

    void set_all_possibilities(int);
    void get_perms_arr(int*,int);
    void reset_active_cluster_probs();
    void test_scope();

public:

    cluster_perms(int,int);
    ~cluster_perms();

    void set_active_cluster(int);
    void set_sub_cluster(int,bool);
    void assign_prob_to_perm(double,int);
    void reset_probabilities(int);
    void set_probability_to_sub_cluster(double,int,int);

    bool already_processed(int,int);
    double get_already_processed_prob();
    
    double* get_highest_prob_cluster();
    int* return_active_sub_cluster_perm(int,int);
    int return_active_sub_cluster_amount();
    int return_am_cluster_groups_in_active();
    int return_active_sub_cluster_width(int);


};

#endif
