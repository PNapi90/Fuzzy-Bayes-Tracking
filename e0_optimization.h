#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

#include "permutations.h"
#include "probabilities_c.h"
#include "cluster_perms.h"
#include "data_class.h"
#include "max_prob_class.h"
#include "threshold_handler.h"

#ifndef E0_OPT_H
#define E0_OPT_H

class e0_optimization{

private:

    const int max_interactions = 10;
    const double Emax_coarse = 10000.;
    const double Emin_coarse = 0.;
    const double mc2 = 511.;

	std::ofstream peak_prob_file;

    bool leading_order;
    bool last_order;
    bool solve_type_bool;
    bool peak_over_sum;
    bool robust_cluster;
    bool bad_cluster_bool;
    bool single_allowed;
    
    bool** ignore_arr;
	
	int verbosity;
	int thr_num;
    int am_clusters;
    int am_clusters_old;
    int internal_counter;
    double l2;
	double binning;
	double Edep_untorn;
	double save_sum;
	double prior_e0;
	
	double* l2_etc;
	
	double* tmp_peak;
	double** first_interaction;
	int tmp_peak_iter;
	
    int* cluster_len;
    double*** clusters;
    double* e0_coarse;
    double* integral_coarse;
    double* edep_sum_arr;
	
    double** prob_pos_of_clusters;
    double*** re_clustered_clusters;
    
    int* lens_of_clusters;

    int* already_processed;

    probabilities_c* prob;
    permutations* perm;
    cluster_perms* CL_Perms;
    data_handler* integral_data;
    data_handler* sum_data;
    max_prob_class* max_prob_obj;
    threshold_handler* thresholds;
    
    //double mean_e;
    //double sig_e;

	bool check_edep();
	
	void get_mean(double**,int);
    void reset();
    void single_point_cluster();
    void set_e0s();
    void get_clusters_total_probs(int,int***);
    
    
    void set_first_interaction(int*,int);
    void check_most_prob(int,int);

    void write_cluster(int);

    int start_int(double);
	double get_position_prior(int);
    double e0_loop(double,int);
    double e0_loop(double,int,bool);
    
	inline double prior_dist(double);
    inline double do_calculations(int*,int,int,int);
    inline double do_calculations(int*,int,int,int,bool);
    inline int factorial(int);
    inline int get_prob_n(double);
    inline double poisson_dist(int,double);
    inline double TS_poisson(double,int);

public:
    
    e0_optimization(int,int,bool);
    ~e0_optimization();

    void re_cluster(double**,int);
    void set_am_clusters(int);
    void reset_re_clustered();

	bool bad_cluster();

    int am_cluster();
    int get_cluster_len(int);
    double** get_cluster(int);
};


#endif
