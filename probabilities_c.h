#ifndef DER_H
#define DER_H


#include "cross_sec.h"
#include "geometry.h"
#include "G_Matrix_slim.h"
#include "prior_handler.h"

class probabilities_c{

private:
	
	const double tolerance_arr[10] = {1,2,2,2.3,1,1,1,1,1,1};
    double source[3];
    const double mc2 = 511;
    const double N = 4.415;
    const double pi = M_PI;
    const double weights[7] = {41.,216.,27.,272.,27.,216.,41.};
    const double weights_milne[5] = {7.,32.,12.,32.,7.};
    const double ti_milne[5] = {0.,0.25,0.5,0.75,1.};
    const double Z = 32;
    const double r0 = 2.81794*1e-15;
    const double mm = 1e-3;
    const double barn = 1e28;
    const bool printer = false;
    const double tolerance_raw = 0.005;
    const double d_tol = 0.0025;
    const double inner_radius = 0.235;
	
	double dist_to_agata[3];
	
	double tolerance;
	bool solve_type_bool;
	bool leading_order;
    bool exiter;
    int am_points;
    double full_integral;
    double photo_to_compton;

    double* ei;
    double* ed;
    double* Edep;
    double** x;
    double* sigma_theta_arr;
    double* sigma_edep_arr;
    double* theta;
    double* alpha;
    double* kn_normed;
    double* gauss_t;
    double* mfpath_exp;
    double* dist;
    double* d_lambda;
    double* pesc_arr;
	
	prior_handler* prior;
    cross_sec* xsec;
    geometry* geom;
    G_Matrix_slim* G_Mat_Slim;


	void get_source_pos();
    void calc_theta();
    void calc_ei_ed();
    void calc_d_lambda();
    void calc_full_integral();
    void calc_gauss_t();
    void calc_kn_normed();
    void set_tolerance(int);
    bool in_range();
    
    inline double calc_Eout(double,double);
    inline double gauss(double,double);
    inline double sigma_e(double);
    inline double Klein_Nishina(double,double);
    inline double milne_gaussian(double,double,double);
    inline double get_air_path(int);
	inline double poisson_dist(double);
    double get_kn_en();
    double kn_en(double,double,double);
    double pesc_kn(double,double);
    double sigma_theta(int);

public:
    probabilities_c(bool);
    ~probabilities_c();
    
    void set_data(double**,int);

    double set_e0(double,bool);
    double return_photo_to_compton();
    
    double get_source_to_first(double*,bool);
    
    bool check_priors(double**,int);
};


#endif
