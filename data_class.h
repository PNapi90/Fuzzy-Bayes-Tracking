#include <fstream>
#include <string>
#include <iostream>

#ifndef DATA_C_H
#define DATA_C_H

class data_handler{

private:

    const double dummy_val = -42.;
    const char* format = "%lf %lf %lf %lf %lf %lf";
    const double c_delim = -99999;
    const double gamma_delim = -1337;
    const int max_am_data = 10;

    int amount;
    int iter;
    int iter_c;
    int len;

    double** internal_ints;
    
    bool thrd_ctor;
    bool scnd_ctor;
    bool fth_ctor;
    bool first_comb;

    int* cluster_size_arr;
    double** input_arr;
    
	double max_integral;

    std::ifstream input;
    std::ofstream output;
    std::ofstream integral_data;
    std::ofstream cluster_data;
    std::ofstream l2_norm_file;
	std::ofstream sum_file;
	std::ofstream lambda_file;
	
    void load_input(bool);

public:
	data_handler(std::string,std::string,int,bool);
	data_handler(int);
	data_handler(int,bool);
	data_handler(int,char);
    
    ~data_handler();
    
	void save_sum_val(double,double,double);
	void save_lambda_vals(double*,int);
    void save(double**,int);
    void set_integral_rest_zero();
    void save_integral(double,double);
    void next();
    void next_iter();
    void next_gamma();
    void save_perm_total(int,int,double);
    void save_l2_norm(double,int,double,double,double);
    void next_perm();
    void next_comb();
    int get_am_clusters();

    double** get_data();
    int* get_cluster_sizes();
};

#endif
