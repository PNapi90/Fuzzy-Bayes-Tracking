#include "fuzzy_c_means.h"
#include "RLFCM.h"
#include "e0_optimization.h"
#include "air_path.h"

#ifndef CLUSTER_H
#define CLUSTER_H

class cluster_class{
private:
    
    const int verbose_thread_num = 1;
    
    bool solve_type_bool;
    bool first_cluster;
    int am_torn;
    int am_cluster;
    int len_TS;
    int thr_num;
    int gamma_counter;
    int verbosity;

    int* cluster_lens;

    double** TS_cluster;
    double*** clusters;

    fuzzy_c_means* fc_mean;
    e0_optimization* e_zero;
    RLFCM* RL_FCM;
    air_path* air;
    
    bool check_etmp(double);
    bool check_for_air(double**,int);
    
	inline int threshold(int,int);

public:
    cluster_class(int,int,bool);
    ~cluster_class();

    void clustering(double**,int);
    void reset();

    double** get_cluster(int);
    int cluster_len(int);
    int amount();
    int get_am_torn();
};


#endif
