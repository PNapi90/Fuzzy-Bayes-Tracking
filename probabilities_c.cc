#include "probabilities_c.h"

using namespace std;

probabilities_c::probabilities_c(bool solve_type_bool){
	
	get_source_pos();
    this->solve_type_bool = solve_type_bool;
	am_points = 0;
    Edep = new double[10];
    x = new double*[10];
    for(int i = 0;i < 10;i++){
		x[i] = new double[3];
		for(int j = 0;j < 3;j++) x[i][j] = 0;
	}
    dist = new double[10];
    theta = new double[10];
    alpha = new double[10];
    ei = new double[11];
    ed = new double[10];
    gauss_t = new double[10];
    mfpath_exp = new double[10];
    d_lambda = new double[10];
    pesc_arr = new double[7];
    kn_normed = new double[10];
	
	prior = new prior_handler();
    xsec = new cross_sec();
    geom = new geometry("LNL");
	leading_order = false;
    G_Mat_Slim = new G_Matrix_slim(this->solve_type_bool);
}

probabilities_c::~probabilities_c(){
    for(int i = 0;i < 10;i++) delete[] x[i];
    delete[] x;
    delete[] Edep;
    delete[] theta;
    delete[] ei;
    delete[] ed;
    delete[] gauss_t;
    delete[] mfpath_exp;
    delete[] d_lambda;
    delete[] pesc_arr;
    delete[] dist;
    delete[] alpha;
    delete[] kn_normed;
	
	delete prior;
    delete xsec;
    delete geom;
    delete G_Mat_Slim;
}


bool probabilities_c::check_priors(double** xtmp,int length){
	return false;
        double prior_tot = 0;
	for(int i = 0;i < length;++i) prior_tot += prior->get_prior(xtmp[i][0],1);
	prior_tot /= (double) length;
	return (prior_tot > 0.2);
}

void probabilities_c::get_source_pos(){
	ifstream source_data("inputdata/source_pos.dat");
	for(int i = 0;i < 3;++i) source_data >> source[i];
}


double probabilities_c::set_e0(double E,bool leading_order){
    
    this->leading_order = leading_order;
    exiter = false;
    ei[0] = E;
    G_Mat_Slim->set_order(this->leading_order);
    
    if(solve_type_bool) sigma_edep_arr = G_Mat_Slim->get_sigmas(x,Edep,am_points,E,solve_type_bool);
    else sigma_theta_arr = G_Mat_Slim->get_sigmas(x,Edep,am_points,E);
    
    calc_theta();
    calc_ei_ed();
    if(exiter) return 0.;

    if(am_points > 1){
		calc_gauss_t();
		calc_d_lambda();
		calc_kn_normed();
	}
	calc_full_integral();
	
    return full_integral;
}
void probabilities_c::set_tolerance(int len){
	//if(leading_order) tolerance = tolerance_raw*tolerance_arr[len-1];
	//else 
	tolerance = tolerance_raw;//*tolerance_arr[len-1];
}

void probabilities_c::set_data(double** data,int len){
    am_points = len;
    set_tolerance(len);
    //cout << "---"<<endl;
    for(int i = 0;i < len;i++){
        Edep[i] = data[i][0];
        for(int j = 0;j < 3;j++) x[i][j] = data[i][j+1]*mm;
		//cout << Edep[i] << " ";
		//or(int j = 0;j < 3;++j) cout << x[i][j] << " ";
		//cout << endl;
    }
    //cout << "---"<<endl;
}

inline double probabilities_c::Klein_Nishina(double e,double cos_angle){
    double E_over_e = 1./(1. + e/mc2*(1. - cos_angle));
    return Z*pi*pow(r0,2)*pow(E_over_e,2.)*(E_over_e + 1./E_over_e - pow(cos_angle,2))*barn;
}

void probabilities_c::calc_kn_normed(){
	double* sigma_tmp = NULL;
    for(int i = 0;i < am_points-1;i++){
		sigma_tmp = xsec->get_sigma(ei[i]);
        kn_normed[i] = Klein_Nishina(ei[i],theta[i])/(sigma_tmp[0] + sigma_tmp[1]);
    }
}

void probabilities_c::calc_theta(){
    double norm1 = 0.;
    double norm2 = 0.;

    for(int i = 0;i < am_points;i++){
        theta[i] = 0.;
        norm1 = 0;
        norm2 = 0;
        if(i == 0){
            for(int j = 0;j < 3;j++){
                theta[i] += (x[i][j] - source[j])*(x[i+1][j] - x[i][j]);
                norm1 += pow(x[i][j] - source[j],2.);
                norm2 += pow(x[i+1][j] - x[i][j],2.);
            }
            dist[i] = sqrt(norm1);
            theta[i] /= dist[i]*sqrt(norm2);
        }
        else if(i < am_points-1){
            for(int j = 0;j < 3;j++){
                theta[i] += (x[i][j] - x[i-1][j])*(x[i+1][j] - x[i][j]);
                norm1 += pow(x[i][j] - x[i-1][j],2.);
                norm2 += pow(x[i+1][j] - x[i][j],2.);
            }
            dist[i] = sqrt(norm1);
            theta[i] /= dist[i]*sqrt(norm2);
        }
        else{
			for(int j = 0;j < 3;j++) norm1 += pow(x[i][j] - x[i-1][j],2.);
            dist[i] = sqrt(norm1);
		}
    }
}


double probabilities_c::sigma_theta(int k){
    //double ret_val = 0;
    //for(int i = 0;i < 2;i++) ret_val += pow(sigma_theta_arr[k][i],2.);
    return 1;//sqrt(ret_val);
}

void probabilities_c::calc_ei_ed(){
    if(solve_type_bool){
        for(int i = 1;i < am_points+1;i++){
            ei[i] = ei[i-1]/(1. + ei[i-1]/mc2*(1. - theta[i-1]));
            ed[i-1] = ei[i-1] - ei[i];
            if(ei[i] < 0 || ed[i-1] < 0 ){
                exiter = true;
                break;
            }
        }
    }
    else{
		//if(leading_order && ei[0] < 1400)cout <<"\n-------\ne0 "<<ei[0] << endl;
        for(int i = 1;i < am_points+1;i++){
            ed[i-1] = Edep[i-1];
            ei[i] = ei[i-1] - ed[i-1];
            if(i < am_points){
				alpha[i-1] = acos(1 + mc2/ei[i-1] - mc2/ei[i]);
				//if(leading_order && ei[0] < 1400) cout <<"ens "<< ei[i] << " " << ed[i-1] << " " << alpha[i-1]*180./M_PI << endl;
				if(alpha[i-1] != alpha[i-1]) {
					exiter = true;
					break;
				}
				
			}
            if(ei[i] < 0 || ed[i-1] < - 5){
                exiter = true;
                break;
            }
        }
        //if(leading_order && ei[0] < 1400)cout << "\n-------\n ";
    }
}

inline double probabilities_c::milne_gaussian(double val,double mean,double sig){
	auto gaussian = [](double t,double sig,double mu){
		return 1./(sqrt(2.*M_PI)*sig)*exp(-0.5*pow((t-mu)/sig,2));
	};
	double summed_peak = 0.;
	double point = 0.;
	double len_diff = sig/10.;
	double a = val - len_diff/2.;
	
	for(int i = 0;i < 5;i++){
		point = a + ti_milne[i]*len_diff;
		summed_peak += weights_milne[i]*gaussian(point,sig,mean);
	}
	return len_diff*summed_peak/90.;
}

void probabilities_c::calc_gauss_t(){
	double angle_temp = 0;
    if(solve_type_bool){
        for(int i = 0;i < am_points-1;i++){
			//if(sigma_edep_arr[i] !=sigma_edep_arr[i] ||sigma_edep_arr[i] == 0 ) cout << "si " << sigma_edep_arr[i] << endl;
            gauss_t[i] = milne_gaussian(ed[i],Edep[i],sigma_edep_arr[i]);
        }
    }
    else{
        for(int i = 0;i < am_points-1;i++){ 
			angle_temp = alpha[i];
			//if(angle_temp > M_PI/2) angle_temp = M_PI - angle_temp;
            gauss_t[i] = milne_gaussian(abs(acos(theta[i])),abs(angle_temp),sigma_theta_arr[i]);
            //if(leading_order && ei[0] < 1200) cout <<"ans " <<acos(theta[i])*180/M_PI << " " << angle_temp*180/M_PI << " " <<sigma_theta_arr[i]*180/M_PI<< endl;
        }
    }
}

inline double probabilities_c::sigma_e(double E){
    return tolerance/2.355*E;
}

inline double probabilities_c::get_air_path(int pos){
	
	if(pos == 0){
		double norm_gamma = 0;
		double gamma_vec[3] = {0.};
		for(int i = 0;i < 3;i++){
			gamma_vec[i] = x[pos][i] - source[i];
			norm_gamma += pow(gamma_vec[i],2);
		}
		norm_gamma = sqrt(norm_gamma);
		double outside_shell = pow(source[2],2) - pow(inner_radius,2);
		double g_times_s = gamma_vec[2]/norm_gamma*source[2];
		return -g_times_s + sqrt(pow(g_times_s,2) - outside_shell);
	}
	
	pos -= 1;
	double norm_gamma = 0;
	double norm_x2 = 0;
	double gamma_vec[3] = {0.};
	for(int i = 0;i < 3;i++){
		gamma_vec[i] = x[pos+1][i] - x[pos][i];
		norm_x2 += pow(x[pos][i],2);
		norm_gamma += pow(gamma_vec[i],2);
	}
	norm_gamma = sqrt(norm_gamma);
	for(int i = 0;i < 3;i++) gamma_vec[i] /= norm_gamma;
	
	double scalar_product = 0;
	for(int i = 0;i < 3;i++) scalar_product += gamma_vec[i]*x[pos][i];
	
	double inside_shell = norm_x2 - pow(inner_radius,2);
	
	if(pow(scalar_product,2) <= inside_shell) return 0.;
	
	double lambdas = 0.;
	for(int i = 0;i < 2;i++){
		lambdas = -scalar_product + pow(-1,i)*sqrt(pow(scalar_product,2) - inside_shell);
		if(lambdas < 0 || lambdas > norm_gamma) return 0.;
	}
	return 2*sqrt(pow(scalar_product,2) - inside_shell);
}

double probabilities_c::get_source_to_first(double* point,bool dummy){
	
	
	
	double norm_gamma = 0;
	double gamma_vec[3] = {0.};
	for(int i = 0;i < 3;i++){
		gamma_vec[i] = point[i]*mm - source[i];
		norm_gamma += pow(gamma_vec[i],2);
	}
	norm_gamma = sqrt(norm_gamma);
	double outside_shell = pow(source[2],2) - pow(inner_radius,2);
	double g_times_s = gamma_vec[2]/norm_gamma*source[2];
	return norm_gamma + g_times_s - sqrt(pow(g_times_s,2) - outside_shell);
}


void probabilities_c::calc_d_lambda(){
    double* sigma_t = NULL;
    double mu_t = 0.;
    double d_tmp = 0;
    double d_calc = 0;
    //if(leading_order)cout << "b" << endl;
    for(int i = 0;i < am_points;i++){
        sigma_t = xsec->get_sigma(ei[i]);
        mu_t = N*(sigma_t[0] + sigma_t[1]);
		d_tmp = get_air_path(i);
		//if(i > 0 && leading_order) cout << "dtmp " << d_tmp << " true dist " << dist[i] << endl;
		d_calc = dist[i] - d_tmp;
		d_lambda[i] = d_calc*mu_t;
        mfpath_exp[i] = (exp(-(d_calc-d_tol/2.)*mu_t) - exp(-(d_calc+d_tol/2.)*mu_t));
        mfpath_exp[i] *= exp(-d_calc*mu_t);
    }
    //if(leading_order)cout << "e" << endl;
    sigma_t = NULL;
}


inline double probabilities_c::calc_Eout(double Ein,double angle){
    return Ein/(1. + Ein/mc2*(1. - cos(angle)));
}


void probabilities_c::calc_full_integral(){

    full_integral = 1.;
    double etmp_sum = 0;
	//if(leading_order && ei[0] < 1200)cout << "--------" << endl;
    for(int i = 0;i < am_points-1;i++){
        full_integral *= gauss_t[i]*kn_normed[i]*mfpath_exp[i];
        etmp_sum += Edep[i];
        if(full_integral == 0) return;
        //cout << ei[0] <<" "<< i <<" "<<gauss_t[i] <<" "<< mfpath_exp[i] <<" "<<kn_normed[i] << endl;
    }
    etmp_sum += Edep[am_points-1];
    if(full_integral != full_integral) cout << "_---------------_" << endl;
    double* sigma_tmp = xsec->get_sigma(ei[am_points-1]);
    double sigma_P = sigma_tmp[0];
    double kn_en = get_kn_en();
    double sigma_norm_val = sigma_P + sigma_tmp[1];
    double photo_prob = sigma_P*gauss(ei[am_points-1],Edep[am_points-1])/sigma_norm_val;

    full_integral *= photo_prob + kn_en/sigma_norm_val;
    full_integral *= prior->get_prior(ei[0],am_points)*prior->get_n_prior(am_points);
    //cout <<"-> full " << full_integral <<" " << photo_prob + kn_en/sigma_norm_val <<" prior " << prior->get_prior(ei[0],am_points)*prior->get_n_prior(am_points) <<endl;
    //poisson_dist(etmp_sum);
    
    if(kn_en == 0 && photo_prob == 0) photo_to_compton = -42;
    else photo_to_compton = (photo_prob - kn_en)/(photo_prob + kn_en);

}

inline double probabilities_c::poisson_dist(double esum){
	auto factorial = [](int n){
		int val = 1;
		for(int i = 1;i <= n;++i) val *= i;
		return val;
	};
	
	double mu_val = 0;
	if(esum < 500) mu_val = 1.;
	else if(esum < 1000) mu_val = 1.3;
	else if(esum < 1500) mu_val = 1.7;
	else if(esum < 2000) mu_val = 2.5;
	else if(esum < 2500) mu_val = 3;
	else mu_val = 3.5;
	return pow(mu_val,am_points)/factorial(am_points)*exp(-mu_val);
}

double probabilities_c::return_photo_to_compton(){return photo_to_compton;}

inline double probabilities_c::gauss(double e,double E){
	double sigma_last = 0;
	if(am_points > 1){
		for(int i = 0;i < am_points;i++){
			sigma_last += pow(tolerance*Edep[i]/2.355,2);
		}
		sigma_last = sqrt(sigma_last);
	}
	else sigma_last = tolerance*Edep[am_points-1]/2.355;
    return milne_gaussian(e,E,sigma_last);
}

bool probabilities_c::in_range(){
    if(ei[0] < 1500 && printer) return true;
    return false;
}

double probabilities_c::kn_en(double E,double e,double Edep_l){
    double angle = acos(1 - mc2/E + mc2/e);
    //escape not possible due to energy of photon
    if(angle != angle) return 0.;
    double ret_val =  Z*pi*pow(r0,2.)*pow(E/e,2.)*(E/e + e/E - pow(sin(angle),2.0))
            *sin(angle)*abs(mc2/(sqrt(1 - pow(1 + mc2/e - mc2/E,2.))*pow(E,2.)))
            *gauss(e-E,Edep_l);

    return ret_val*barn;
}

double probabilities_c::get_kn_en(){
    double** mu = new double*[2];
    for(int i = 0;i < 2;i++){
		mu[i] = new double[3];
		for(int j = 0;j < 3;j++) mu[i][j] = 0;
	}
    if(am_points > 1){
        for(int i = 0;i < 2;i++){
            for(int j = 0;j < 3;j++) mu[i][j] = x[am_points-2+i][j];
        }
    }
    else{
        for(int j = 0;j < 3;j++){
            mu[0][j] = source[j];
            mu[1][j] = x[0][j];
        }
    }
    double pcomp = 0.;
    double Ec = ei[am_points];
    double int_energies[5];
    for(int i = 0;i < 5;i++){
        int_energies[i] = Ec*(1 - tolerance) + i*0.25*2.*tolerance*Ec;
        if(int_energies[i] < 0.) int_energies[i] = 0.;
    }
    geom->escape_dist(mu,int_energies,ei[am_points-1],Ec);
    double temp_val = 0;
    double sig_norming = xsec->get_sigma(ei[am_points-1])[1];
    double* esc_lens = geom->get_esc_len();
    
    for(int i = 0;i < 5;i++){
        temp_val = kn_en(int_energies[i],ei[am_points-1],ed[am_points-1]);
        if(temp_val != temp_val) temp_val = 0.;
        pcomp += temp_val*pesc_kn(xsec->get_sigma(int_energies[i])[1],esc_lens[i]);
    }

    pcomp *= 2.*tolerance*Ec/90.;

    for(int i = 0;i < 2;i++) delete[] mu[i];
    delete[] mu;

    if(pcomp != pcomp) pcomp = 0;

    return pcomp;
}

double probabilities_c::pesc_kn(double sigma,double esc_len){
    if(esc_len == -42) return 0.;
    return exp(-esc_len*N*sigma);
}
