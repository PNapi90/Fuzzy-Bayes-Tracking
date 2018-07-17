#include "cross_sec.h"
#include "data_class.h"

#include <cmath>
#include <fstream>
#include <string>

#ifndef AIR_PATH_H
#define AIR_PATH_H

class air_path{

private:

	bool all_out;

	cross_sec* xsec;
	data_handler* data;
	
	const double r0 = 2.81794*1e-15;
    double source[3];

    const double multiplier_arr[10] = {1.33,1.33,1.33,1.1,1,0.8,0.7,0.6,0.6,0.6};
	
	double Edep_sum;
    double inner_radius_2;
	double max_dist;
	
	int verbosity;
	int thr_num;
	
	double sigma_w;
	int max_tupel[2];
	int len;
	int max_d_int;
	
    int* amount_air;
    double max_dist_to_mean;
    double* distance_to_mean;
    double* mean;
    double* gamma;
    double** dist;
    double** x;
    double*** lambda_arr;
	
	std::ofstream data_matrix;
	
	bool distance_matrix();
	bool check_angles();
	bool check_poss_diff(double**);
	bool check_num_inter(int);
	bool multiplier(double**);
	bool check_sphere();
	bool lambda_sphere();
	bool abort_statement(double,double,double);
	int air_passages(int);
	
	
	inline double line_function(double);
	
    double scalar_pr(double*);
    double get_radius_diff(double*);
    double get_lambda(double,double,int);
    double max_travel_dist(double);
    double get_factor(double**,double);
    double KN(double,double);

	double get_theta(double);
	double get_angle(double*,double*);
	
	void get_source_pos();
    void reset_lambda();
    void get_mean(double**,int);
    void get_gamma(double*,double*,int,int);

public:
    air_path(int,int);
    ~air_path();
    
    bool check_air(double**,int);
};



#endif
