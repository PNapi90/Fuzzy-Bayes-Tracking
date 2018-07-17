#include "prior_handler.h"

using namespace std;

prior_handler::prior_handler(){
    prior_list = NULL;
    prior_n = NULL;
    prior_low = NULL;
    load_input();
    
}

prior_handler::~prior_handler(){
    for(int i = 0;i < max_len;++i) delete[] prior_list[i];
    for(int i = 0;i < 9;++i) delete[] prior_low[i];
    delete[] prior_low;
    delete[] prior_list;
    delete[] prior_n;
}


void prior_handler::load_input(){
    
    ifstream input_stream("inputdata/prior_normed.dat");
    
    prior_list = new double*[max_len];
    prior_n = new double[max_N];
    double full_sum = 0;
    for(int i = 0;i < max_N;++i) prior_n[i] = 0;
    for(int i = 0;i < max_len;++i){
        prior_list[i] = new double[max_N];
        for(int j = 0;j < max_N;++j){
			input_stream >> prior_list[i][j];
			prior_n[j] += prior_list[i][j];
			full_sum += prior_n[j];
		}
    }
        
    ifstream low_stream("inputdata/prior_low.dat");
    
    prior_low = new double*[9];
    
    for(int i = 0;i < 9;++i){
        prior_low[i] = new double[max_N];
        for(int j = 0;j < max_N;++j){
			low_stream >> prior_low[i][j];
			prior_n[j] += prior_low[i][j];
			full_sum += prior_n[j];
		}
    }
    for(int i = 0;i < max_N;++i) prior_n[i] /= full_sum;
    
    
    
}

double prior_handler::get_n_prior(int N){return prior_n[N-1];}

double prior_handler::get_prior(double e0,int N){
	if(e0 < 100) return low_E_prior(e0,N);
	
    double interpol = e0/100.;
    int e0_pos = (int) interpol - 1;
   // cout << "e0 -> " << e0 << " " << interpol << " " << e0_pos << endl;
    return (prior_list[e0_pos+1][N-1] - prior_list[e0_pos][N-1])/100.*(interpol - e0_pos + 1.) + prior_list[e0_pos][N-1];
}


inline double prior_handler::low_E_prior(double e0,int N){
	if(e0 < 10){
		if(N == 1) return 1;
		return 0;
	}
	
	double interpol = e0/10.;
    int e0_pos = (int) interpol - 1;
    //cout << "e0 " << e0 << " " << interpol << " " << e0_pos << endl;
    if(e0_pos < 8) return (prior_low[e0_pos+1][N-1] - prior_low[e0_pos][N-1])/10.*(interpol - e0_pos + 1.) + prior_low[e0_pos][N-1];
    
    return (prior_list[0][N-1] - prior_low[e0_pos][N-1])/10.*(interpol - e0_pos + 1.) + prior_low[e0_pos][N-1];

}
