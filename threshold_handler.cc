#include "threshold_handler.h"

using namespace std;

threshold_handler::threshold_handler(){
	thresholds = new double*[20];
	for(int i = 0;i < 20;++i) thresholds[i] = new double[2];
	maximum_inter = 5;
	
	max_prob_1 = new double[2];
	
	set_thresholds();

}

threshold_handler::~threshold_handler(){
	for(int i = 0;i < 20;++i) delete[] thresholds[i];
	delete[] thresholds;
	delete[] max_prob_1;
}


void threshold_handler::print_thresholds(){
	cout << "-------------\n";
	cout << "thresholds" << endl;
	cout << "-------------\n";
	for(int i = 0;i < maximum_inter;++i){
		for(int j = 0;j < 2;++j) cout << thresholds[i][j] << "\t";
		cout << endl;
	}
	cout << ".............\n";
	cout << max_prob_1[0] << "\t" << max_prob_1[1] << endl;
	cout << "-------------\n";
}

void threshold_handler::set_thresholds(){
	ifstream file("inputdata/thresholds.dat");
	for(int i = 0;i < maximum_inter;++i) for(int j = 0;j < 2;++j) file >> thresholds[i][j];
	
	ifstream max_file("inputdata/thresholds_max.dat");
	max_file >> max_prob_1[0];
	max_file >> max_prob_1[1];
}


bool threshold_handler::compare(double val_peak,double val_sum,int am_points){
	if(am_points > maximum_inter) return true;
	
	am_points -= 1;
	bool statements[4] = {false};
	
	if(am_points == 0){
		statements[0] = log(val_peak) > thresholds[am_points][0];
		statements[1] = log(val_sum) > thresholds[am_points][1];
		statements[2] = log(val_peak) < max_prob_1[0];
		statements[3] = log(val_sum) < max_prob_1[1];
		
		bool all_statements = true;
		for(int i = 0;i < 4;++i) all_statements = all_statements && statements[i];
		return all_statements;
	}
	else{
		statements[0] = log(val_peak) > thresholds[am_points][0];
		statements[1] = log(val_sum) > thresholds[am_points][1];
		return (statements[0] && statements[1]);
	}
	return false; 
}
