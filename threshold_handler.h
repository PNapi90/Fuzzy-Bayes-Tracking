#ifndef THRESHOLD_HANDLER_H
#define THRESHOLD_HANDLER_H 1

#include <iostream>
#include <fstream>
#include <cmath>

class threshold_handler{
private:
	
	int maximum_inter;
	
	double* max_prob_1;
	double** thresholds;

	void set_thresholds();

public:
	threshold_handler();
	~threshold_handler();

	void print_thresholds();

	bool compare(double,double,int);
	
};

#endif
