#include "RLFCM.h"
#include <cmath>

using namespace std;

RLFCM::RLFCM(){
	
	alphas = new double[20];
	x = new double*[20];
	centroids = new double*[20];
	old_centroids = new double*[20];
	dist_2 = new double*[20];
	weights = new double*[20];
	for(int i = 0;i < 20;i++){
		weights[i] = new double[20];
		dist_2[i] = new double[20];
		x[i] = new double[3];
		centroids[i] = new double[3];
		old_centroids[i] = new double[3];
	}
	
}

RLFCM::~RLFCM(){
	for(int i = 0;i < 20;i++){
		delete[] weights[i];
		delete[] dist_2[i];
		delete[] x[i];
		delete[] centroids[i];
		delete[] old_centroids[i];
	}
	delete[] alphas;
	delete[] x;
	delete[] centroids;
	delete[] old_centroids;
	delete[] dist_2;
	delete[] weights;
}


void RLFCM::robust_clustering(double** data,int am_points){
	this->am_points = am_points;
	am_clusters = am_points;
	for(int i = 0;i < am_points;i++){
		for(int j = 0;j < 3;j++){
			x[i][j] = data[i][j+1]*1e-3;
			centroids[i][j] = x[i][j];
		}
		alphas[i] = 1./((double) am_clusters);
	}
	for(int i = 0;i < 3;i++) r[i] = 1;
	epsilon = 10;
	int iter = 1;
	int dummy = 1;
	int am_loops = 0;
	double epsilon_old = 0;
	while(epsilon > abort_val && am_loops < 5){
		calc_new_dist();
		calc_new_weights();
		calc_new_r1_r2(iter);
		calc_new_alphas_r3();
		calc_new_am_clusters();
		renorm_alphas_weights();
		calc_new_centroids();
		epsilon_old = epsilon;
		epsilon = get_new_epsilon();
		if(epsilon == epsilon_old) am_loops += 1;
		iter += 1;
	}
	
	
}

inline double RLFCM::get_min(double a,double b){
	if(a >= b) return b;
	return a;
}

double RLFCM::get_new_epsilon(){
	double dist = 0;
	double max_dist = 0;
	for(int k = 0;k < old_am_clusters;k++){
		dist = 0;
		for(int j = 0;j < 3;j++) dist += pow(centroids[k][j] - old_centroids[k][j],2);
		dist = sqrt(dist);
		if(max_dist < dist) max_dist = dist;
	}
	return max_dist;
}


void RLFCM::save_old_centroids(){
	old_am_clusters = am_clusters;
	for(int k = 0;k < am_clusters;k++){
		for(int j = 0;j < 3;j++) old_centroids[k][j] = centroids[k][j];
	}
}

void RLFCM::calc_new_centroids(){
	double weight_sum = 0;
	for(int k = 0;k < am_clusters;k++){
		for(int i = 0;i < am_points;i++){
			weight_sum += weights[i][k];
		}
		for(int j = 0;j < 3;j++){
			centroids[k][j] = 0;
			for(int i = 0;i < am_points;i++){
				centroids[k][j] += weights[i][k]*x[i][j];
			}
			centroids[k][j] /= weight_sum;
		}
	}
}

void RLFCM::renorm_alphas_weights(){
	double renorm_alpha = 0;
	double renorm_weights[am_points] = {0};
	for(int t = 0;t < am_clusters;t++){
		renorm_alpha += alphas[t];
	}
	for(int i = 0;i < am_points;i++){
		for(int t = 0;t < am_clusters;t++){
			renorm_weights[i] += weights[i][t];
		}
	}
	
	for(int k = 0;k < am_clusters;k++){
		alphas[k] /= renorm_alpha;
	}
	for(int i = 0;i < am_points;i++){
		for(int k = 0;k < am_clusters;k++){
			weights[i][k] /= renorm_weights[i];
		}
	}
	
}

void RLFCM::calc_new_am_clusters(){
	
	
	save_old_centroids();
	
	double one_over_n = 1./((double) am_points);
	
	int reducer = 0;
	for(int k = 0;k < am_clusters;k++){
		if(alphas[k] < one_over_n) reducer += 1;
	}
	am_clusters -= reducer;
}

void RLFCM::calc_new_alphas_r3(){
	
	double one_over_n = 1./((double) am_points);
	
	double old_alphas[am_clusters] = {0};
	for(int k = 0;k < am_clusters;k++) old_alphas[k] = alphas[k];
	
	double tmp_sum = 0;
	double weight_sum = 0;
	for(int t = 0;t < am_clusters;t++) tmp_sum += old_alphas[t]*log(old_alphas[t]);
	for(int k = 0;k < am_clusters;k++){
		weight_sum = 0;
		for(int i = 0;i < am_points;i++) weight_sum += weights[i][k];
		alphas[k] = one_over_n*weight_sum + r[2]/r[0]*alphas[k]*(log(alphas[k]) - tmp_sum);
	}
	double max_weight = 0.;
	double max_entropy = old_alphas[0]*tmp_sum;
	double tmp_weight_sum = 0;
	for(int k = 0;k < am_clusters;k++){
		tmp_weight_sum = 0;
		for(int i = 0;i < am_points;i++) tmp_weight_sum += weights[i][k];
		tmp_weight_sum *= one_over_n;
		
		if(max_weight < tmp_weight_sum) max_weight = tmp_weight_sum;
		if(max_entropy < old_alphas[k]*tmp_sum) max_entropy = old_alphas[k]*tmp_sum;
		
	}
	
	r[2] = get_min(1./((double) am_clusters),(1. - max_weight)/(-max_entropy));
}


void RLFCM::calc_new_r1_r2(int iter){
	r[0] = exp(-iter/10.);
	r[1] = exp(-iter/100.);
}

void RLFCM::calc_new_dist(){
	double tmp = 0;
	for(int i = 0;i < am_points;i++){
		for(int k = 0;k < am_clusters;k++){
			tmp = 0;
			for(int s = 0;s < 3;s++) tmp += pow(x[i][s] - centroids[k][s],2);
			dist_2[i][k] = tmp;
		}
	}
}


void RLFCM::calc_new_weights(){
	double tmp_sum = 0;
	for(int i = 0;i < am_points;i++){
		for(int k = 0;k < am_clusters;k++){
			tmp_sum = 0;
			for(int t = 0;t < am_clusters;t++){
				tmp_sum += exp((-dist_2[i][t] + r[0]*log(alphas[t]))/r[1]);
			}
			weights[i][k] = exp((-dist_2[i][k] + r[0]*log(alphas[k]))/r[1])/tmp_sum;
		}
	}
}

int RLFCM::get_am_clusters(){return am_clusters;}
