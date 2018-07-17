#ifndef PRIOR_H
#define PRIOR_H 1

#include <fstream>
#include <iostream>
class prior_handler{

private:
	
	const int max_len = 140;
	const int max_N = 20;
	
    double** prior_list;
    double** prior_low;
    double* prior_n;

    void load_input();
    
    inline double low_E_prior(double,int);

public:
    prior_handler();
    ~prior_handler();

    double get_prior(double,int);
	double get_n_prior(int);
};


#endif
