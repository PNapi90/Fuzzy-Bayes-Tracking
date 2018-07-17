#include <fstream>
#include <iostream>

#ifndef CROSS_H
#define CROSS_H

class cross_sec{

private:
    double** comp_xsec;
    double** photo_xsec;

    double** comp_xsec_normed;
    double** photo_xsec_normed;

    double* sigma;

    void load_xsecs();
    void normalize();
    void get_sigma_from_list(double);
    double get_sigma_from_list_normed(double);

public:
    cross_sec();
    ~cross_sec();
    
    double* get_sigma(double);
    double get_sigma_normed(double);
};

#endif