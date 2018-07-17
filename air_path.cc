#include <iostream>

#include "air_path.h"

using namespace std;

air_path::air_path(int verbosity,int thr_num){
	
	
	
	data_matrix.open("files/distances/distances_"+to_string(thr_num)+".dat");
	data = new data_handler(thr_num,'c');
	len = 0;
	this->thr_num = thr_num;
	this->verbosity = verbosity;
	xsec = new cross_sec();
	get_source_pos();
    inner_radius_2 = pow(235,2);
    gamma = new double[3];
    amount_air = new int[20];

    for(int i = 0;i < 20;i++) amount_air[i] = 0;

    lambda_arr = new double**[2];
    for(int s = 0;s < 2;s++){
        lambda_arr[s] = new double*[20];
        for(int i = 0;i < 20;i++){
            lambda_arr[s][i] = new double[20];
            for(int j = 0;j < 20;j++){
                lambda_arr[s][i][j] = 0;
            }
        }
    }
    dist = new double*[20];
    for(int i = 0;i < 20;i++){
        dist[i] = new double[20];
        for(int j = 0;j < 20;j++) dist[i][j] = 0;
    }
    mean = new double[3];
	for(int i = 0;i < 3;++i) mean[i] = 0;
    distance_to_mean = new double[20];
    for(int i = 0;i < 20;++i) distance_to_mean[i] = 0;
    sigma_w = 0;
    all_out = false;
}

air_path::~air_path(){
    
    for(int i = 0;i < 20;i++) delete[] dist[i];
    delete[] dist;

    for(int i = 0;i < 2;i++){
        for(int j = 0;j < 20;j++) delete[] lambda_arr[i][j];
        delete[] lambda_arr[i];
    }
    delete[] lambda_arr;
    delete[] mean;
    delete[] amount_air;
    delete[] gamma;
    delete[] distance_to_mean;
    delete xsec;
    delete data;
}


void air_path::get_source_pos(){
	ifstream source_data("inputdata/source_pos.dat");
	for(int i = 0;i < 3;++i){
		source_data >> source[i];
		source[i] *= 1e3;
	}
	
	if(verbosity > 0 && thr_num == 1){
		cout << "\n. . . . . . . . . . \n";
		cout << "source position: ";
		cout << source[0] << " " << source[1] << " " << source[2];
		cout << "\n. . . . . . . . . . \n";
		cout << endl;
	}
}


void air_path::get_mean(double** x,int len){
	sigma_w = 0;
	if(verbosity > 0 && thr_num == 1) cout << "mean at: ";
	for(int i = 0;i < 3;i++){
		mean[i] = 0;
		for(int j = 0;j < len;j++){
			mean[i] += x[j][0]*x[j][i+1];
		}
		mean[i] /= Edep_sum;
		if(verbosity > 0 && thr_num == 1) cout << mean[i] << " ";
	}
	if(verbosity > 0 && thr_num == 1) cout << endl;
	double tmp_norm = 0;
	for(int j = 0;j < len;j++){
		tmp_norm = 0;
		for(int i = 0;i < 3;i++){
			tmp_norm += pow(x[j][i+1] - mean[i],2);
		}
		sigma_w += tmp_norm;
	}
	sigma_w = sqrt(sigma_w);
	sigma_w /= (double) len;
	
	tmp_norm = 0;
	double sigma_mean = 0;
	double mean_dist_dist = 0;
	max_dist_to_mean = 0;
	//if(verbosity > 0) cout << "dists" << endl;
	for(int i = 0;i < len;i++){
		tmp_norm = 0;
		for(int j = 0;j < 3;j++) tmp_norm += pow(mean[j] - x[i][j+1],2);
		distance_to_mean[i] = sqrt(tmp_norm);
		//if(verbosity > 0) cout << distance_to_mean[i] << " ";
		if(max_dist_to_mean < distance_to_mean[i]){
			max_dist_to_mean = distance_to_mean[i];
			max_d_int = i;
		}
		//mean_dist_dist += distance_to_mean[i];
	}
	
	for(int i = 0;i < len;i++){
		//distance_to_mean[i] *= (Edep_sum - x[i][0])/Edep_sum;
		mean_dist_dist += distance_to_mean[i]; 
	}
	//if(verbosity > 0) cout << endl;
	mean_dist_dist /= (double) len;
	for(int i = 0;i < len;i++) sigma_mean = pow(mean_dist_dist - distance_to_mean[i],2);
	sigma_mean = sqrt(sigma_mean);
	//sigma_w = sigma_mean;
	//if(verbosity > 0) cout <<"\n mean_mean "<<mean_dist_dist << " " << sigma_mean << endl;
	
	//if(verbosity > 0) cout << "\n-> sigma: " << sigma_w << endl;
	sigma_w = sigma_mean;
}
bool air_path::check_sphere(){
		
	double radius = 0;
	double sphere_rad = pow(sigma_w,2);
	int am_in = 0;
	for(int i = 0;i < len;i++){
		for(int j = 0;j < 3;j++) radius += pow(x[i][j+1] - mean[j],2);
		if(radius < sphere_rad) am_in += 1;
		radius = 0;
	}
	if(am_in < len){
		if(am_in > 0) all_out = false;
		else all_out = true;
		return true;
	}
	return false;
}

bool air_path::abort_statement(double lam_over_max_dist,double diff_over_esum,double sig_e_over_diff){
	
	bool below_thres = false;
	if(lam_over_max_dist >= 1) return true;
	if(lam_over_max_dist < 0.38) return false;
	//if(lam_over_max_dist < 0.42) below_thres = true;	
	if(diff_over_esum < line_function(lam_over_max_dist)) return false;//below_thres = true;
	return true;
	if(below_thres && sig_e_over_diff < 0.175) return false;
	return true;
}


bool air_path::lambda_sphere(){
	bool check_it = check_sphere();
	double diff = Edep_sum - x[max_d_int][0];
	if(verbosity > 0 && thr_num == 1) cout << "diff " << diff << endl;
	if(diff < 0.5*Edep_sum) diff = x[max_d_int][0];
	double* sigma = xsec->get_sigma(diff);
	double lambda = 1./(4.415*(sigma[0] + sigma[1]))*1e3;
	if(verbosity > 0 && thr_num == 1){
		//cout << "sphere: " << check_sphere() << " " << all_out << endl;
		cout << "Energy remains: " << diff << " lambda: " << lambda << " " << sigma_w <<" "<< lambda/sigma_w << endl;
		cout << lambda/diff<<" "<< lambda/Edep_sum << endl;
		cout << diff/Edep_sum << " " << line_function(lambda/max_dist_to_mean) << endl;
	}
	double len_tmp = (double) len;
	double mean_e = Edep_sum/len_tmp;
	double sig_e = 0;
	for(int i = 0;i < len;i++){
		sig_e += pow(x[i][0] - mean_e,2);
	}
	sig_e = sqrt(sig_e)/len_tmp;
	
	double* tmp_stream = new double[6];
	tmp_stream[0] = Edep_sum;
	tmp_stream[1] = diff;
	tmp_stream[2] = lambda;
	tmp_stream[3] = sigma_w;
	tmp_stream[4] = max_dist_to_mean;
	tmp_stream[5] = sig_e;
	
	
	data->save_lambda_vals(tmp_stream,6);
	delete[] tmp_stream;
	
	return abort_statement(lambda/max_dist_to_mean,diff/Edep_sum,sig_e/diff);
	
	//return check_it; 
	//return true;
	/*
	double len_tmp = (double) len;
	double mean_e = Edep_sum/len_tmp;
	double sig_e = 0;
	for(int i = 0;i < len;i++){
		sig_e += pow(x[i][0] - mean_e,2);
	}
	sig_e = sqrt(sig_e)/len_tmp;
	
	double* tmp_stream = new double[6];
	tmp_stream[0] = Edep_sum;
	tmp_stream[1] = diff;
	tmp_stream[2] = lambda;
	tmp_stream[3] = sigma_w;
	tmp_stream[4] = max_dist_to_mean;
	tmp_stream[6] = sig_e;
	
	
	data->save_lambda_vals(tmp_stream,6);
	delete[] tmp_stream;
	
	double radius = 0;
	double sphere_rad = pow(2*lambda,2);
	int amount_in_sphere = 0;
	for(int i = 0;i < len;i++){
		for(int j = 0;j < 3;j++) radius += pow(x[i][j+1] - mean[j],2);
		if(radius < sphere_rad) amount_in_sphere += 1;
		radius = 0;
	}
	if(amount_in_sphere == len && all_out) return true;
	if(amount_in_sphere < len && all_out) return false;
	if(amount_in_sphere < len && !all_out) return true;
	
	return true;*/
}

inline double air_path::line_function(double val){
    double t = 0.86;
    double end_point = 1.45;
    double slope = (-t + 0.5)/end_point;	
    return slope*val + t;
}

bool air_path::check_air(double** x,int len){
	
	this->x = x;
	this->len = len;
	
	reset_lambda();
	max_dist = 0;
	
	Edep_sum = 0;
	for(int i = 0;i < len;i++) Edep_sum += x[i][0];
	get_mean(x,len);
    double scalar_product = 0;
    double radius_diff_2 = 0;
    for(int s = 0;s < 2;s++){
        for(int i = 0;i < len;i++){
            radius_diff_2 = get_radius_diff(x[i]);
            //cout << x[i][0] << " ";
            for(int j = i+1;j < len;j++){
                get_gamma(x[i],x[j],i,j);
                scalar_product = scalar_pr(x[i]);
                lambda_arr[s][i][j] = get_lambda(scalar_product,radius_diff_2,s);
                //if(verbosity > 0)cout << " -> lam " << lambda_arr[s][i][j] << endl;
                if(lambda_arr[s][i][j] < dist[i][j] && lambda_arr[s][i][j] > 0) amount_air[i] += 1;
            }
            //cout << endl;
        }
    }
    int air = air_passages(len);
    if(air == 1) return true;
    else if(air > 1) return false;
    return distance_matrix();    
//if(check_sphere()) 
    if(!lambda_sphere()) return false;
    if(!distance_matrix()) return false;
    return true;
	//else return false;
	
    //if(check_angles()) return true;
    return false;
    if(Edep_sum > 1400) return false;
    return true;
    
    
    if(check_poss_diff(x)) return true;
    //else if(max_dist < max_travel_dist(Edep_sum)) return true;
    return false;
}

bool air_path::distance_matrix(){
	double full_sum = 0;
	double* sigma_tmp = NULL;
	double lambda = 0;
	for(int i = 0;i < len;++i) full_sum += x[i][0];
	double diff_tmp = 0;
    double distance_matrix[len][len] = {0.};
    for(int i = 0;i < len;++i){
		diff_tmp = full_sum - x[i][0];
		sigma_tmp = xsec->get_sigma(diff_tmp);
		lambda = 1./(4.415*(sigma_tmp[0] + sigma_tmp[1]))*1e3;
        for(int j = 0;j < len;++j){
            for(int k = 0;k < 3;++k) distance_matrix[i][j] += pow(x[i][k+1] - x[j][k+1],2);
            distance_matrix[i][j] = sqrt(distance_matrix[i][j]);
            data_matrix << distance_matrix[i][j] << "\n";
            distance_matrix[i][j] /= lambda;
		}
    }
    data_matrix << len << "\n";
    data_matrix << -99999 << "\n";
    data_matrix.flush();
    
    for(int i = 0;i < len;++i){
    	for(int j = 0;j < len;++j) if(distance_matrix[i][j] >= 2 && distance_matrix[j][i] >= 2) return false;
    }
    return true;


}



bool air_path::check_angles(){
	
	double diff = Edep_sum - x[max_d_int][0];
	if(verbosity > 0 && thr_num == 1) cout << "diff " << diff << endl;
	if(diff < 0.5*Edep_sum) diff = x[max_d_int][0];
	double* sigma = xsec->get_sigma(diff);
	double lambda = 1./(4.415*(sigma[0] + sigma[1]))*1e3;
	
	double meanNorm = 0;
	double norm_furthest = 0;
	for(int i = 0;i < 3;i++){
		meanNorm += pow(mean[i] - source[i],2);
		norm_furthest += pow(x[max_d_int][i+1] - source[i],2);
	}
	
	meanNorm = sqrt(meanNorm);
	norm_furthest = sqrt(norm_furthest);
	
	double max_angle = atan2(2*lambda,meanNorm);
	
	double gamma_angle = pow(meanNorm,2) + pow(max_dist_to_mean,2) - pow(norm_furthest,2);
	gamma_angle /= 2*meanNorm*max_dist_to_mean;
	gamma_angle = acos(gamma_angle);
	
	
	double x_min = max_dist_to_mean - 5.;
	//double x_max = max_dist_to_mean + 5.;
	double alpha_min = atan2(sin(gamma_angle)*x_min,meanNorm - cos(gamma_angle)*x_min);
	
	double comp_angle = acos(meanNorm/(2*norm_furthest) + norm_furthest/(2*meanNorm)
							 - pow(max_dist_to_mean,2)/(2*norm_furthest*meanNorm));
	int i = max_tupel[0];
	int j = max_tupel[1];
	
	double angle_pos = get_angle(x[i],x[j]);
	double angle_E = get_theta(x[i][0]);
							 
	double KN_normed = KN(Edep_sum,angle_pos)/xsec->get_sigma(Edep_sum)[1];
	
	if(verbosity > 0 && thr_num == 1) cout <<"KN_Norm " << KN_normed << " pos "<<angle_pos*180/M_PI <<" en "<<angle_E*180/M_PI << endl;
	if(verbosity > 0 && thr_num == 1) cout << "max_dist " << max_dist_to_mean << " " << lambda <<  " " <<sigma_w << " " << sigma_w/lambda << " " <<sigma_w/max_dist_to_mean << endl;
	if(verbosity > 0 && thr_num == 1) cout <<"angles: " << max_angle*180/M_PI << " " << comp_angle*180/M_PI  << " min angle " << alpha_min*180/M_PI<< endl;
	if(max_angle > alpha_min) return true;
	return false;
	
	
	
	
}


bool air_path::check_num_inter(int len){
	int typical = 0;
	if(Edep_sum > 1600) typical = 4;
	else if(Edep_sum > 800) typical = 3;
	else typical = 2;
	
	if(typical == len) return true;
	return false;
	
	
}

double air_path::get_angle(double* x,double* y){
	double norm1 = 0;
	double norm2 = 0;
	double scalar = 0;
	for(int i = 0;i < 3;i++){
		//cout << x[i+1] << " " << y[i+1] << endl;
		norm1 += pow(x[i+1] - source[i],2.);
		norm2 += pow(y[i+1] - x[i+1],2.);
		scalar += (x[i+1] - source[i])*(y[i+1] - x[i+1]);
	}
	return acos(scalar/sqrt(norm1*norm2));
	
}

double air_path::get_theta(double Ex){
	return acos(1 - 511/(Edep_sum - Ex) + 511/Edep_sum);
}

double air_path::KN(double E,double theta){
	double E_over_e = 1./(1 + E/511*(1. - cos(theta)));
	return 32*pow(r0,2)*pow(E_over_e,2)*(E_over_e + 1./E_over_e - pow(sin(theta),2))*1e28;
}

bool air_path::multiplier(double** x){
	int i = max_tupel[0];
	int j = max_tupel[1];
	
	double* sigma = xsec->get_sigma(Edep_sum);
	double N = 4.415;
	double lambda = 1./(N*sigma[0] + N*sigma[1])*1e3;
	if(verbosity > 0 && thr_num == 1) cout << "lam " << lambda << " "  << sigma_w <<" " << lambda/sigma_w << " "<< max_dist_to_mean <<  endl;
		
	if(lambda/sigma_w > 4) return true;
	return false;
}

double air_path::get_factor(double** x,double mu){
	
	int i = max_tupel[0];
	int j = max_tupel[1];
	
	double angle_pos = get_angle(x[i],x[j]);
	double angle_E = get_theta(x[i][0]);
	double diff = Edep_sum-x[max_d_int][0];
	if(diff < 0.5*Edep_sum) diff = x[max_d_int][0];
	double* sigma = xsec->get_sigma(diff);
	double KN_normed = KN(Edep_sum,angle_pos)/sigma[1];
	
	//if(verbosity > 0) cout <<"KN_Norm " << KN_normed << endl;
	double N = 4.415;
	double lambda = 1./(N*sigma[0] + N*sigma[1])*1e3;
	if(verbosity > 0 && thr_num == 1) cout << "lam " << lambda << " "  << sigma_w <<" " << lambda/sigma_w << " "<< max_dist_to_mean <<  endl;
	
	if(lambda/max_dist_to_mean > 1) return 4;
	else return 1;
	
		
	if(2*lambda > max_dist_to_mean && sigma_w/max_dist_to_mean < 0.5 ) return 4;
	else return 1;
	
	/*
	if(sigma_w < 0.5*lambda) return 3;
	if(sigma_w > 0.5*lambda && sigma_w < 1.5*lambda) return 4;
	else if(sigma_w > 1.5*lambda && sigma_w < 3*lambda) return 3.;
	else return 1;
	
	
	
	if(angle_E != angle_E) return 2;
	if(verbosity > 0) cout << angle_pos*180./M_PI << " " << angle_E*180/M_PI << endl;
	if(abs(abs(angle_pos) - angle_E) < 20./180.*M_PI) return 1+log(10);
	else if(abs(abs(angle_pos) - angle_E) < 60./180.*M_PI) return 0.5+log(10);
	else return log(10);*/
}

bool air_path::check_poss_diff(double** x){
	
	return multiplier(x);
	
	
	double Ei = x[max_tupel[0]][0];
	double* sigma = xsec->get_sigma(Edep_sum-Ei);
	double mu = 4.415*(sigma[0] + sigma[1]);
	double mfp_max = get_factor(x,mu)/mu*1000;
	
	double mfp_min = 0.1/mu*1000;
	double air_tmp = 0;
	bool air_bool = true;
	for(int i = 0;i < 2;i++){
		if(lambda_arr[i][max_tupel[0]][max_tupel[1]] < 0 || lambda_arr[i][max_tupel[0]][max_tupel[1]] > dist[max_tupel[0]][max_tupel[1]]) air_bool = false;
	}
	if(air_bool) air_tmp = lambda_arr[0][max_tupel[0]][max_tupel[1]] - lambda_arr[1][max_tupel[0]][max_tupel[1]];
	max_dist -= air_tmp;
	//cout<<"-> " << Ei << " " << Edep_sum-Ei <<" atmp " <<air_tmp << " " << mfp_max << " " << mfp_min << " " << max_dist << endl;
	if(max_dist_to_mean > mfp_min && max_dist_to_mean < mfp_max) return true;
	//cout << "b" << endl;
	return false; 
}

double air_path::max_travel_dist(double E){
	double* sigma = xsec->get_sigma(E);
	double mu = 4.415*(sigma[0] + sigma[1]);
	double en_diff = 3/mu*1000;
	if(en_diff < 120) return 120;
	return en_diff;
}


int air_path::air_passages(int len){
    int count_air = 0;
    for(int i = 0;i < len;i++){
        if(amount_air[i] > 0) count_air += 1;
    }
    return count_air;
}


double air_path::get_lambda(double scalar_product,double radius_diff_2,int s){
     if(pow(scalar_product,2) - radius_diff_2 <= 0) return 0;
     return -scalar_product + pow(-1,s)*sqrt(pow(scalar_product,2) - radius_diff_2);
}

double air_path::get_radius_diff(double* x){
    double val = 0;
    for(int i = 0;i < 3;i++) val += pow(x[i+1],2);
    return val - inner_radius_2;
}


double air_path::scalar_pr(double* x){
    double val = 0;
    for(int i = 0;i < 3;i++) val += x[i+1]*gamma[i];
    return val;
}

void air_path::reset_lambda(){
    for(int s = 0;s < 2;s++){
        for(int i = 0;i < 20;i++){
            for(int j = 0;j < 20;j++) lambda_arr[s][i][j] = 0;
        }
    }
    for(int i = 0;i < 20;i++) amount_air[i] = 0;
}

void air_path::get_gamma(double* x,double* y,int i,int j){
    double norm = 0;
    //if(verbosity > 0)cout << endl;
    for(int i = 0;i < 3;i++){
        gamma[i] = -x[i+1] + y[i+1];
        norm += pow(gamma[i],2);
    }
    norm = sqrt(norm);
    dist[i][j] = norm;
    //if(verbosity > 0) cout << dist[i][j] << " ";
    if(max_dist < dist[i][j]){
		max_dist = dist[i][j];
		max_tupel[0] = i;
		max_tupel[1] = j;
	}
    for(int i = 0;i < 3;i++) gamma[i] /= norm;

}
