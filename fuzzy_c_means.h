#include <random>

#ifndef FUZZY_MEAN_H
#define FUZZY_MEAN_H

class fuzzy_c_means{

private:

    const int max_len = 10;
    const double abort_val = 1e-3;
    const double DMAX = 80.;
	
    bool first_iter;
    int fuzzyness;
    int verbosity;
    int am_points;
    int am_clusters;
    int thr_num;    
    double dmax;
    double epsilon;
    double sum_edep;

    int* cluster_to_points;
    int* cluster_len;
    
    double** x;
    double** weights;
    double** old_weights;
    double** centroids;
    double** dist_2;

    double*** fc_clusters;

    void get_am_clusters();
    void set_dmax();
    void init_weights();
    void adjust_old_weights();
    void calc_centroids();
    void calc_dist2();
    void calc_weights();
    void calc_epsilon();
    void set_cluster_to_points();
    void print_info();
    void sort_clusters();

public:
    fuzzy_c_means(double,int,int);
    ~fuzzy_c_means();
    
    void fuzzy_clustering(double**,int,int);
    void reset();
    
    int return_am_clusters();
    int return_cluster_len(int);
	double get_edep_in_cluster();
	
    double** return_cluster(int);
};



#endif
