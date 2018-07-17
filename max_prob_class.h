#ifndef MAX_PROB_H
#define MAX_PROB_H 1

class max_prob_class{

private:
	
	double* return_arr;
	
	double** l2_peak;
	double** l2_sum;
	double** ent_peak;
	double** ent_sum;
	
	double** probs;
	double** sum_val;
public:
	max_prob_class();
	~max_prob_class();

	void store(double*,int,int);
	
	double get_prob(int,int);
	double get_sum(int,int);
	
	double* get_all_vals(int,int);
};

#endif
