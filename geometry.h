#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "Detector_Geometry.h"

#ifndef GEOM_H
#define GEOM_H

class geometry{
private:
	
	const double r_agata_inner = 0.235; // 17cm
	const double r_agata_outer = 0.325; // 17cm + 9cm
	const double phi_agata = 2.*M_PI;
	const double theta_agata = 37.5/180.*M_PI;
	const int bi_sec_it = 12;
	
	std::string geom_type;
	double E_in;
	
	bool* test_arr;
	double* angles;
	double* lambda_vec;
	double** mu;
	double** x_bi;
	double* esc_lens;
	
	Detector_Geometry* Agata_Setup;
	
	bool in_cin;
	Eigen::Vector3d* dir;
	Eigen::Vector3d* rot;
	Eigen::Vector3d* dir_after;
	const double pi = M_PI;
	const double weights[7] = {41.,216.,27.,272.,27.,216.,41.};
	const double nodes[7] = {0.,1.,2.,3.,4.,5.,6.};
	const double det_dim = 0.04;
	
	void set_rotation_vec();
	double run_bisek(double);
	double find_lambda_min_a(double);
	double find_lambda_min_c(double);	
	double get_lambdas(double,double);
	void test_inside();
	
public:
	geometry(std::string);
	~geometry();
	
	void escape_dist(double**,double*,double,double);
	void escape_dist(double**,double*);
	
	double* get_esc_len();
	
};

#endif
