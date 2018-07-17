#include <iostream>
#include <cmath>
#include <fstream>

#ifndef  G_MAT_SLIM
#define G_MAT_SLIM

class G_Matrix_slim{

private:

    const double tolerance_raw = 0.002;
    const double delta_x_arr[10] = {1,1.25,1.3,1.4,1.5,3,2.2,2.4,2.6,2.8};
    const double tolerance_arr[10] = {1,1.2,1.4,1.6,1.8,1,1,1,1,1};
    double source[3];
    const double delta_x = 0.005/2.355;
    const double mc2 = 511;
	
	bool leading_order;
    bool solve_type_bool;
	
    int interactions;
    int dim_m;
    double e0;
    double tolerance;


    double* d_x_mat;
    double* d_e_mat;
    
    double* error_ei;

    double* Edep;
    double* ei;
    double* x_scalars;
    double* x_norms;
    double* x_angles;
    double* delta_edep_arr;
    double* delta_theta_arr;
    double** x;
    double** sigma_vec;

	void get_source_pos();

    void set_derivatives();
    
    void set_sigmas();
    void set_sigmas(bool);
    
    void set_angles();
    void set_Edeps();
    void set_error_ei();
    
    void calc_d_x_mat();
    void calc_d_x_mat(int);
    
    void calc_d_e_mat();
    void calc_d_e_mat(bool);
    
    void calc_scalars();
    void calc_norms();

    double sigma_e(int);
    double x_derivatives(int,int);
    void set_tolerances();
    
    double get_delta_x(int);
    inline double delta_ei(int,double);
    inline double delta_edepi(int,double);
    inline double set_error_ei(int);
    inline double delta_theta(int,double);

public:

    G_Matrix_slim(bool);
    ~G_Matrix_slim();
    
    void set_order(bool);

    double* get_sigmas(double**,double*,int,double);
    double* get_sigmas(double**,double*,int,double,bool);
    
};


#endif
