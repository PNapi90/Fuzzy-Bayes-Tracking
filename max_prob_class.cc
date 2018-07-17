#include "max_prob_class.h"

max_prob_class::max_prob_class(){
	probs = new double*[300];
	sum_val = new double*[300];
	l2_peak = new double*[300];
	l2_sum = new double*[300];
	ent_peak = new double*[300];
	ent_sum = new double*[300];
	return_arr = new double[6];
	
	for(int i = 0;i < 300;++i){
		probs[i] = new double[300];
		sum_val[i] = new double[300];
		l2_peak[i] = new double[300];
		l2_sum[i] = new double[300];
		ent_peak[i] = new double[300];
		ent_sum[i] = new double[300];
		for(int j = 0;j < 300;++j){
			probs[i][j] = 0;
			sum_val[i][j] = 0;
			l2_peak[i][j] = 0;
			l2_sum[i][j] = 0;
			ent_peak[i][j] = 0;
			ent_sum[i][j] = 0;
		}
	}
}

max_prob_class::~max_prob_class(){
	for(int i = 0;i < 300;++i){
		delete[] probs[i];
		delete[] sum_val[i];
		delete[] l2_peak[i];
		delete[] l2_sum[i];
		delete[] ent_peak[i];
		delete[] ent_sum[i];
	}
	delete[] probs;
	delete[] sum_val;
	delete[] return_arr;
	delete[] l2_peak;
	delete[] l2_sum;
	delete[] ent_peak;
	delete[] ent_sum;
}

void max_prob_class::store(double* vals,int c_iter,int k){
	probs[c_iter][k] = vals[0];
	sum_val[c_iter][k] = vals[1];
	l2_peak[c_iter][k] = vals[2];
	l2_sum[c_iter][k] = vals[3];
	ent_peak[c_iter][k] = vals[4];
	ent_sum[c_iter][k] = vals[5];
}

double max_prob_class::get_prob(int c_iter,int k){return probs[c_iter][k];}

double max_prob_class::get_sum(int c_iter,int k){return sum_val[c_iter][k];}

double* max_prob_class::get_all_vals(int c_iter,int k){
	return_arr[0] = probs[c_iter][k];
	return_arr[1] = sum_val[c_iter][k];
	return_arr[2] = l2_peak[c_iter][k];
	return_arr[3] = l2_sum[c_iter][k];
	return_arr[4] = ent_peak[c_iter][k];
	return_arr[5] = ent_sum[c_iter][k];
	
	return return_arr;
}
