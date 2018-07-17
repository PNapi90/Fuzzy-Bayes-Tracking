#include <iostream>

#ifndef RL_FCM_H
#define RL_FCM_H

class RLFCM{

private:

	const double abort_val = 2e-3;

	int am_points;
	int am_clusters;
	int old_am_clusters;

	double r[3];
	double epsilon;

	double* alphas;

	double** weights;
	double** x;		
	double** centroids;
	double** old_centroids;
	double** dist_2;
	
	
	void save_old_centroids();
	void calc_new_centroids();
	void renorm_alphas_weights();
	void calc_new_am_clusters();
	void calc_new_alphas_r3();
	void calc_new_r1_r2(int);
	void calc_new_dist();
	void calc_new_weights();
	
	double get_new_epsilon();
	inline double get_min(double,double);
	
	

public:
	RLFCM();
	~RLFCM();

	void robust_clustering(double**,int);
	
	int get_am_clusters();
};


#endif
