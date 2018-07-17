#include "geometry.h"
#include <cstdlib>
#include <iostream>

using namespace std;

//ctor
geometry::geometry(string geom_type){
	this->geom_type = geom_type;
	
	test_arr = new bool[3];
	x_bi = new double*[3];
	for(int i = 0;i < 3;i++){
		x_bi[i] = new double[3];
		for(int j = 0;j < 3;j++) x_bi[i][j] = 0;
		test_arr[i] = 1;
	}	
	lambda_vec = new double[6];
	esc_lens = new double[5];
    
	angles = new double[5];
    
	rot = new Eigen::Vector3d;
	mu = new double*[2];
	dir = new Eigen::Vector3d;
	dir_after = new Eigen::Vector3d;
	rot->setZero();
	dir->setZero();
	dir_after->setZero();
	for(int i = 0;i < 2;i++) mu[i] = new double[3];
	
	Agata_Setup = new Detector_Geometry(this->geom_type);

}

//dtor
geometry::~geometry(){
	for(int i = 0;i < 2;i++){
		delete[] mu[i];
		delete[] x_bi[i];
	}
	delete[] x_bi[2];
	delete[] x_bi;
	
	delete[] lambda_vec;
	delete[] esc_lens;
	delete dir;
	delete[] angles;
	delete[] mu;
	delete rot;
	delete[] test_arr;
	delete dir_after;
	delete Agata_Setup;
}
/*
void geometry::escape_dist(double** mu_vec,double* thetas){
    for(int i = 0;i < 2;i++) for(int j = 0;j < 3;j++) mu[i][j] = mu_vec[i][j];
    //calculate direction of ingoing photon and scattering angles (and normalize vector)
    double norm_dir = 0.;
    for(int i = 0;i < 3;i++){
        (*dir)(i) = mu[1][i] - mu[0][i];
        norm_dir += pow((*dir)(i),2.);
    }
    norm_dir = sqrt(norm_dir);
    for(int i = 0;i < 3;i++) (*dir)(i) /= norm_dir;
    for(int i = 0;i < 7;i++) angles_2[i] = thetas[i];
    set_rotation_vec();  
    
    for(int i = 0;i < 7;i++){
        if(angles_2[i] != angles_2[i]) esc_lens_2[i] = -42.;
        else esc_lens_2[i] = find_lambda_min_a(angles_2[i]);
    }
}*/

void geometry::escape_dist(double** mu_vec,double* e_milne,double E,double Ec){
	for(int i = 0;i < 2;i++) for(int j = 0;j < 3;j++) mu[i][j] = mu_vec[i][j];
	E_in = E;
	
	//calculate direction of ingoing photon and scattering angles (and normalize vector)
	double norm_dir = 0.;
	for(int i = 0;i < 3;i++){
		(*dir)(i) = mu[1][i] - mu[0][i];
		norm_dir += pow((*dir)(i),2.);
	}
	norm_dir = sqrt(norm_dir);
	for(int i = 0;i < 3;i++) (*dir)(i) /= norm_dir;
    for(int i = 0;i < 5;i++) angles[i] = acos(1. - 511/e_milne[i] + 511/E);
	set_rotation_vec();
	for(int i = 0;i < 5;i++){
		if(angles[i] != angles[i]) esc_lens[i] = -42.;
		else esc_lens[i] = find_lambda_min_a(angles[i]);	
	}
	
}
/*
double* geometry::get_esc_lens(){
    return esc_lens;
}
*/

double geometry::find_lambda_min_a(double angle){
	Eigen::AngleAxis<double> aa(angle,*rot);
	(*dir_after) = aa*(*dir);
	
	//1 cm stepsize
	double stepsize = 0.01;
	int maxstep = 1000;
	double lambda = 0.;
	
	for(int i = 0;i < 3;i++){
		x_bi[0][i] = mu[1][i];
		x_bi[2][i] = x_bi[0][i];
	}	
	for(int i = 0;i < maxstep;i++){
		for(int j = 0;j < 3;j++) x_bi[2][j] += stepsize*(*dir_after)(j);
		test_inside();
		if(!test_arr[2]){
			lambda = run_bisek(stepsize);
			lambda += i*stepsize;
			break;
		}
		for(int j = 0;j < 3;j++) x_bi[0][j] = x_bi[2][j];
	}
	return lambda;
}

double geometry::run_bisek(double dist){
	double stepper = dist/2.;
	double sign = 0.;
	for(int j = 0;j < 3;j++) x_bi[1][j] = x_bi[0][j] + stepper*(*dir_after)(j);
	for(int i = 0;i < bi_sec_it;i++){
		test_inside();
		if(test_arr[1] == test_arr[0]){
			for(int j = 0;j < 3;j++) x_bi[0][j] = x_bi[1][j];
			sign = 1.;
		}
		else if(test_arr[1] == test_arr[2]){
			for(int j = 0;j < 3;j++) x_bi[2][j] = x_bi[1][j];
			sign = -1.;
		}
		stepper += sign*dist/(4.*(i+1));
		for(int j = 0;j < 3;j++) x_bi[1][j] += stepper*(*dir_after)(j);
	}
	return stepper;
}

void geometry::test_inside(){
	double r = 0.;
	double phi = 0.;
	double theta = 0.;
	for(int i = 0;i < 3;i++){
		r = 0.;
		for(int j = 0;j < 3;j++) r += pow(x_bi[i][j],2.);
		r = sqrt(r);	
		phi = atan2(x_bi[i][1],x_bi[i][0]);
		theta = acos(x_bi[i][2]/r);
		if(!Agata_Setup->check_inside(r,abs(theta),phi)) test_arr[i] = false;
		else test_arr[i] = true;
	}	
}

void geometry::set_rotation_vec(){
	//calculate a random rotation-vector perpendicular to
	//photon direction dir
	int iter;
	for(int i = 0;i < 3;i++) (*rot)(i) = (double)(rand() % 500 + 1)/500.;
	for(int i = 0;i < 3;i++){
		if((*dir)(i) != 0) iter = i;
		break; 
	}
	//set rotation-vector to allow perpendicular character
	if(iter == 0) (*rot)(0) = -((*dir)(1)*(*rot)(1)+(*dir)(2)*(*rot)(2))/(*dir)(0);
	else if(iter == 1) (*rot)(1) = -((*dir)(0)*(*rot)(0)+(*dir)(2)*(*rot)(2))/(*dir)(1);
	else if(iter == 2) (*rot)(2) = -((*dir)(0)*(*rot)(0)+(*dir)(1)*(*rot)(1))/(*dir)(2);
	
	//normalize rotation-vector
	double summe = 0.;
	for(int i = 0;i < 3;i++) summe += pow((*rot)(i),2.);
	summe = sqrt(summe);
	for(int i = 0;i < 3;i++) (*rot)(i) /= summe;
}

double* geometry::get_esc_len(){return esc_lens;}
