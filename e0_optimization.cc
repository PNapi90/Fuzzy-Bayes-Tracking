#include "e0_optimization.h"
		
using namespace std;

e0_optimization::e0_optimization(int thr_num,int verbosity,bool solve_type_bool){
    
    binning = 1.;
    
    l2_etc = new double[4];
    
    this->solve_type_bool = solve_type_bool;
    peak_prob_file.open("files/peak_prob/peak_prob_file_Eu_n_opt_"+to_string(thr_num)+".dat");
    last_order = false;
    leading_order = false;
    this->thr_num = thr_num;
    this->verbosity = verbosity;
    if(this->verbosity > 0 && thr_num == 1){
		cout << "*** Verbosity set to " << verbosity << " ***" <<  endl;
		cout << "------------------------------------------------\n";
		cout << "using bins of width: " << binning << " keV" << endl;
		cout << "------------------------------------------------\n";
    }
    
    if(this->verbosity > 1) integral_data = new data_handler(thr_num,this->solve_type_bool);
    thresholds = new threshold_handler();
    
    if(this->verbosity > 0 && thr_num == 1) thresholds->print_thresholds();
    am_clusters_old = 0;
    sum_data = new data_handler(thr_num);
    ignore_arr = new bool*[20];
    for(int i = 0;i < 20;i++){
		ignore_arr[i] = new bool[20];
		for(int j = 0;j < 20;j++) ignore_arr[i][j] = false;
	}
	single_allowed = false;
    prob = new probabilities_c(this->solve_type_bool);    
    perm = new permutations();
    CL_Perms = new cluster_perms(this->verbosity,thr_num);
    
    if(this->verbosity > 0 && thr_num == 1 ){
		cout << "======================\n";
		cout << "cluster perms set\n";
		cout << "======================\n" << endl;
	}
    internal_counter = 0;
    
    first_interaction = new double*[100];
    for(int i = 0;i < 100;++i){
		first_interaction[i] = new double[3];
		for(int j = 0;j < 3;++j) first_interaction[i][j] = 0;
	}
    cluster_len = new int[10];
    clusters = new double**[10];
    for(int i = 0;i < 10;i++){
        clusters[i] = new double*[20];
        
        for(int j = 0;j < 20;j++) clusters[i][j] = new double[4];
    }
        
    already_processed = new int[3];
    
    peak_over_sum = false;
    
    //mean_e = 0;
    //sig_e = 0;
    prior_e0 = 0;
    robust_cluster = false;
    tmp_peak_iter = 0;
    tmp_peak = NULL;
    save_sum = 0;
    
    max_prob_obj = new max_prob_class();
    
    set_e0s();
}

e0_optimization::~e0_optimization(){
	peak_prob_file.close();
    reset_re_clustered();
    delete CL_Perms;
    //cout << " -> deleted CL_Perms object" << endl;
    delete prob;
    //cout << " -> deleted prob_cluster object" << endl;   
    delete perm;
    //cout << " -> deleted permutation object" << endl;
	if(verbosity > 1) delete integral_data;
    for(int i = 0;i < 10;i++){
        for(int j = 0;j < 20;j++){
            delete[] clusters[i][j];
        }
        delete[] clusters[i]; 
    }
    for(int i = 0;i < 100;++i) delete[] first_interaction[i];
	delete[] first_interaction;
    delete sum_data;
    for(int i = 0;i < 20;i++) delete[] ignore_arr[i];
	delete[] ignore_arr;
    delete[] clusters;
    delete[] cluster_len;
    delete[] e0_coarse;
    delete[] integral_coarse;
    delete[] already_processed;
    delete[] l2_etc;
    delete max_prob_obj;
    delete thresholds;
}

void e0_optimization::get_mean(double** x,int len){
	double len_tmp = (double) len;
	double mean_e = Edep_untorn/len_tmp;
	double sig_e = 0;
	for(int i = 0;i < len;i++){
		sig_e += pow(x[i][0] - mean_e,2);
	}
	sig_e = sqrt(sig_e)/len_tmp;
	if(verbosity > 0 && thr_num == 1) cout << "robustness analysis: " << mean_e << " " << sig_e << " " << mean_e/(2.355*sig_e) << endl;
	if(mean_e/(2.355*sig_e) > 5 || mean_e/(2.355*sig_e) < 2) robust_cluster = false;
	else robust_cluster = true;
}

bool e0_optimization::bad_cluster(){return bad_cluster_bool;}

void e0_optimization::re_cluster(double** x,int len){
	
	bad_cluster_bool = false;
	Edep_untorn = 0;
    cluster_len[internal_counter] = len;
    for(int j = 0;j < cluster_len[internal_counter];j++){
        for(int k = 0;k < 4;k++) clusters[internal_counter][j][k] = x[j][k];
        Edep_untorn += x[j][0];
    }
    
    
    if(len == 1){
		if(verbosity > 0 && thr_num == 1 ) cout << "amount of points in Fuzzy cluster: " << 1 << endl;
		single_point_cluster();
		return;
	}
	if(len == 0){
		cerr << "zero encoutered" << endl;
		am_clusters = 0;
		am_clusters_old = 0;
		return;
	}
	get_mean(x,len);
    
    int permutation_val = cluster_len[internal_counter];
    if(verbosity > 0 && thr_num == 1){
		if(robust_cluster) cout << "\n-> robust cluster <-\n" << endl;
		else cout << "\n-> non-robust cluster <-\n" << endl;
		cout << "amount of points in Fuzzy cluster: " << permutation_val << endl;
    }  
    CL_Perms->set_active_cluster(len);
    
    int* sub_cluster = NULL;
    int sub_cluster_amount = 0;
    int sub_cluster_width = 0;

    int all_possible_perms_len = factorial(len);
	tmp_peak = new double[all_possible_perms_len];
    tmp_peak_iter = 0;
    
    int cluster_iter = 1;
    int all_possible_perms = CL_Perms->return_am_cluster_groups_in_active();
    
    double prob_tmp = 0;
    double prob_total = 1;
    
    double* ret_arr = NULL; 
	
	double*** prob_tmp_arr = new double**[all_possible_perms+1];
	for(int i = 0;i < all_possible_perms + 1;i++){
		prob_tmp_arr[i] = new double*[all_possible_perms_len];
		for(int j = 0;j < all_possible_perms_len;j++){
			prob_tmp_arr[i][j] = new double[len];
			for(int k = 0;k < len;k++) prob_tmp_arr[i][j][k] = 0.;
		}
	}
	
    prob_pos_of_clusters = new double*[all_possible_perms];
    for(int i = 0;i < all_possible_perms;i++) prob_pos_of_clusters[i] = new double[2];
    double p_total_max = 0;
    int m_iter = 0;
    double pos_prior = 1;
    
    double** sum_saver = NULL;
    double*** rest_saver = NULL;
    double* store_arr = new double[6];    
    
    double edep_tmp = 0;
    while(cluster_iter <= all_possible_perms){
		if(cluster_iter == all_possible_perms) last_order = true;
		else last_order = false;
        
        CL_Perms->set_sub_cluster(cluster_iter,true);
        sub_cluster_amount = CL_Perms->return_active_sub_cluster_amount();
        
        sum_saver = new double*[sub_cluster_amount];
        rest_saver = new double**[sub_cluster_amount];
        for(int i = 0;i < sub_cluster_amount;++i){
			sum_saver[i] = new double[all_possible_perms_len];
			rest_saver[i] = new double*[all_possible_perms_len];
			for(int j = 0;j < all_possible_perms_len;++j)rest_saver[i][j] = new double[4];
		}
        if(verbosity > 0 && thr_num == 1 ){
			cout << "-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n" << endl;
			cout << "cluster iter " << cluster_iter << " all_possible_perms " << all_possible_perms << endl;
			cout << "beginning loop -> " << all_possible_perms_len << " * ";
			cout << sub_cluster_amount << " = " << all_possible_perms_len*sub_cluster_amount << endl;
		}
		
        p_total_max = 0;
        m_iter = 0;
        if(last_order) single_allowed = prob->check_priors(x,len);
        for(int j = 0;j < all_possible_perms_len;j++){
            prob_total = 1;
            for(int k = 0;k < sub_cluster_amount;k++){
                if(cluster_iter == 1) leading_order = true;
                else leading_order = false;
                
                sub_cluster = CL_Perms->return_active_sub_cluster_perm(j,k);
                sub_cluster_width = CL_Perms->return_active_sub_cluster_width(k);

                if(sub_cluster_width >= 1) prob_tmp = do_calculations(sub_cluster,sub_cluster_width,k,j);
                else prob_tmp = 0;
                sum_saver[k][j] = save_sum;
                for(int s = 0;s < 4;++s) rest_saver[k][j][s] = l2_etc[s];
                //cout << "sum_saver " << j << " " << k << " " << sum_saver[k][j] << endl;
                
                //if(cluster_iter > 1) set_first_interaction(sub_cluster,k);
        
                if(prob_tmp == 0){
					
                    //prob_total = 0;
                    //break;
                }
                prob_tmp_arr[cluster_iter][j][k] = prob_tmp;
                prob_total *= prob_tmp;
            }
            //if(cluster_iter > 1) pos_prior = get_position_prior(sub_cluster_amount);
            //prob_total *= pos_prior;
            //prob_total *= poisson_dist(sub_cluster_amount,2);
            if(p_total_max < prob_total){
					m_iter = j;
					p_total_max = prob_total;
			}
			if(verbosity > 1 && thr_num == 1) integral_data->save_perm_total(cluster_iter,j,prob_total);
            CL_Perms->assign_prob_to_perm(prob_total,j);
        }
        
        for(int k = 0;k < sub_cluster_amount;++k){
			store_arr[0] = prob_tmp_arr[cluster_iter][m_iter][k];
			store_arr[1] = sum_saver[k][m_iter];
			for(int s = 0;s < 4;++s) store_arr[s+2] = rest_saver[k][m_iter][s];
			max_prob_obj->store(store_arr,cluster_iter,k);
		}
		
        if(verbosity > 1) integral_data->next_comb();
        
        ret_arr = CL_Perms->get_highest_prob_cluster();
        /*if(ret_arr[0] == 0 && leading_order){
			for(int j = 0;j < all_possible_perms;j++) CL_Perms->assign_prob_to_perm(tmp_peak[j],j);
			ret_arr = CL_Perms->get_highest_prob_cluster();
			robust_cluster = true;
			if(verbosity > 0) cout << "-> switching to absolute values!" << endl;
		}*/
        if(verbosity > 0 && thr_num == 1){
			cout << "\n=================================\n";
			cout << "Total prob " << ret_arr[0] << " @ perm # " << ret_arr[1] << endl;
			for(int s = 0;s < sub_cluster_amount;s++) cout << prob_tmp_arr[cluster_iter][m_iter][s] << " ";  
			cout << "\n=================================\n";
			cout << "\n-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --" << endl;
        }
        if(verbosity > 1 && leading_order) check_most_prob(cluster_iter,m_iter);
        for(int s = 0;s < sub_cluster_amount;s++){
			if(prob_tmp_arr[cluster_iter][m_iter][s] < 0.01) ignore_arr[cluster_iter-1][s] = true;
		}
        for(int s = 0;s < 2;s++) prob_pos_of_clusters[cluster_iter-1][s] = ret_arr[s];
        
        for(int k = 0;k < sub_cluster_amount;++k) prob_pos_of_clusters[cluster_iter - 1][0] /= sum_saver[k][m_iter];
        prob_pos_of_clusters[cluster_iter - 1][0] *= poisson_dist(sub_cluster_amount,get_prob_n(Edep_untorn));//*TS_poisson(sub_cluster_amount,1);
        cluster_iter += 1;
        
        for(int k = 0;k < sub_cluster_amount;++k){
			delete[] sum_saver[k];
			for(int j = 0;j < all_possible_perms_len;++j) delete[] rest_saver[k][j];
			delete[] rest_saver[k];
		}
        delete[] sum_saver;
        delete[] rest_saver;
        rest_saver = NULL;
        sum_saver = NULL; 
    }
	if(verbosity > 1) integral_data->next_gamma();
    write_cluster(all_possible_perms);

    CL_Perms->reset_probabilities(len);
	
	for(int i = 0;i < all_possible_perms+1;i++){
		for(int j = 0;j < all_possible_perms_len;j++) delete[] prob_tmp_arr[i][j];
		delete[] prob_tmp_arr[i];
	}
	delete[] prob_tmp_arr;
	
    for(int i = 0;i < all_possible_perms;i++) delete[] prob_pos_of_clusters[i];
    delete[] prob_pos_of_clusters;
    peak_over_sum = false;
    
    for(int i = 0;i < 20;i++){
		for(int j = 0;j < 20;j++) ignore_arr[i][j] = false;
	}
	delete[] tmp_peak;
	delete[] store_arr;
}

inline double e0_optimization::TS_poisson(double mu,int n){return pow(mu,n)/factorial(n)*exp(-mu);}

double e0_optimization::get_position_prior(int sub_cluster_amount){
	int max_perms_am = factorial(sub_cluster_amount)/(factorial(sub_cluster_amount - 2)*2);
	int norm_iter = 0;
	
	double norms[max_perms_am] = {0};
	
	for(int i = 0;i < sub_cluster_amount;++i){
		for(int j = i+1;j < sub_cluster_amount;++j){
			for(int k = 0;k < 3;++k) norms[norm_iter] += pow(first_interaction[i][k] - first_interaction[j][k],2);
			norms[norm_iter] = sqrt(norms[norm_iter]);
			norm_iter += 1;
		}
	}
	double pos_prior_tmp = 1;
	for(int i = 0;i < max_perms_am;i++) pos_prior_tmp *= prior_dist(norms[i]);
	return pos_prior_tmp;
}

inline int e0_optimization::get_prob_n(double Edep){
	if(Edep < 1500) return 1;
	else if(Edep < 2000) return 1.5;
	else if(Edep < 3000) return 2.;
	else if(Edep < 4000) return 2.5;
	else return 3;
        //if(len <= 2) return 1;
//	else if(len <= 5) return 2;
//	else if(len <= 7) return 2.5;
//	else return 3;
}

inline double e0_optimization::prior_dist(double norm){
	return 1./((60. - sqrt(2*M_PI)))*(1. - exp(-0.5*pow(norm/(2.355*5*sqrt(2)),2)));
}


void e0_optimization::set_first_interaction(int* sub_cluster,int k){
	for(int i = 0;i < 3;++i) first_interaction[k][i] = clusters[internal_counter][sub_cluster[0]-1][i+1];
}


void e0_optimization::check_most_prob(int cluster_iter,int m_iter){
	bool SAVE = true;
	int sub_cluster_amount = CL_Perms->return_active_sub_cluster_amount();
    int j = m_iter;
    double prob_tmp = 0;
    
    int* sub_cluster = NULL;
    int sub_cluster_width = 0;
    for(int k = 0;k < sub_cluster_amount;k++){
		sub_cluster = CL_Perms->return_active_sub_cluster_perm(j,k);
        sub_cluster_width = CL_Perms->return_active_sub_cluster_width(k);
        prob_tmp = do_calculations(sub_cluster,sub_cluster_width,k,j,SAVE);
	}
	integral_data->next_iter();
}


void e0_optimization::reset_re_clustered(){
	if(re_clustered_clusters != NULL){
		if(lens_of_clusters != NULL){
			for(int i = 0;i < am_clusters_old;i++){
				if(re_clustered_clusters[i] != NULL){
					for(int j = 0;j < lens_of_clusters[i];j++){
						if(re_clustered_clusters[i][j] != NULL) delete[] re_clustered_clusters[i][j];
					}
					delete[] re_clustered_clusters[i];
				}
			}
			delete[] lens_of_clusters;
		}
		delete[] re_clustered_clusters;
	}
}

bool e0_optimization::check_edep(){
	if(abs(Edep_untorn - 1173.2) < 10) return true;
	if(abs(Edep_untorn - 1332.5) < 10) return true;
	return false;
}

inline double e0_optimization::poisson_dist(int n,double mu){
	return pow(mu,n)/((double) factorial(n))*exp(-mu);
}



void e0_optimization::write_cluster(int all_possible_perms){
           
    
    double max_prob = 0;
    double max_pos_t = 0;
    
    int best_comb = 0;
	int max_looper = all_possible_perms - 1;
	if(single_allowed) max_looper += 1;
	
    for(int i = 0;i < max_looper;i++){
        if(max_prob < prob_pos_of_clusters[i][0]){
            max_prob = prob_pos_of_clusters[i][0];
            max_pos_t = prob_pos_of_clusters[i][1];
            best_comb = i;
            if(verbosity > 0 && thr_num == 1) cout << "maxs " << max_prob << " @ perm " << i+1 << endl;
        }
    }
    if(max_prob < 0.15 && best_comb == 0) bad_cluster_bool = true;
    if(best_comb != 0 && check_edep() && verbosity > 0 && thr_num == 1) cout << "torn!" << endl;
    if(max_prob == 0){
		if(all_possible_perms == 0) all_possible_perms = 1;
        max_prob = prob_pos_of_clusters[all_possible_perms-1][0];
        max_pos_t = prob_pos_of_clusters[all_possible_perms-1][1];
        best_comb = all_possible_perms-1;
    }
    best_comb += 1;
    int max_pos = (int) max_pos_t;
    CL_Perms->set_sub_cluster(best_comb,false);
    am_clusters = CL_Perms->return_active_sub_cluster_amount();
        
    lens_of_clusters = new int[am_clusters];
    re_clustered_clusters = new double**[am_clusters];
    am_clusters_old = am_clusters;
    
    int* sub_cluster = NULL;
    int tmp_iter = 0;
    double edep_tmp = 0;
    
    bool compare_val = true;
    
    bool nulled = false;
    int null_am = 0;
    double* save_ptr_arr = NULL;
    
    double* smallest_dist_arr = new double[3];
    double smallest_dist = 0;
	
    for(int i = 0;i < am_clusters;i++){
		edep_tmp = 0;
		smallest_dist = 0;
        lens_of_clusters[tmp_iter] = CL_Perms->return_active_sub_cluster_width(i);
        
        compare_val = thresholds->compare(max_prob_obj->get_prob(best_comb,i),max_prob_obj->get_sum(best_comb,i),lens_of_clusters[tmp_iter]);
        if(!compare_val && false){
			null_am += 1;
			nulled = true;
			continue;
		}
        re_clustered_clusters[tmp_iter] = new double*[lens_of_clusters[tmp_iter]];
        sub_cluster = CL_Perms->return_active_sub_cluster_perm(max_pos,tmp_iter);
        for(int j = 0;j < lens_of_clusters[tmp_iter];j++){
            re_clustered_clusters[i][j] = new double[4];
            for(int k = 0;k < 4;k++){
                re_clustered_clusters[tmp_iter][j][k] = clusters[internal_counter][sub_cluster[j]-1][k];
            }
            if(j == 0){
				for(int k = 0;k < 3;++k) smallest_dist_arr[k] = clusters[internal_counter][sub_cluster[j]-1][k+1];
				smallest_dist = prob->get_source_to_first(smallest_dist_arr,true);
			}
            edep_tmp += re_clustered_clusters[tmp_iter][j][0];
        }
        save_ptr_arr = max_prob_obj->get_all_vals(best_comb,i);
        peak_prob_file << edep_tmp << "\t";
        for(int s = 0;s < 6;++s) peak_prob_file << save_ptr_arr[s] << "\t";
        peak_prob_file << smallest_dist << "\t" <<lens_of_clusters[tmp_iter] << "\n";
        tmp_iter += 1;
    }
    if(nulled){
		if(verbosity > 0 && thr_num == 1){
			cout << "~> deleted " << null_am << " out of " << am_clusters_old << " clusters" << endl;
		}
		for(int i = tmp_iter;i < am_clusters_old;++i) re_clustered_clusters[i] = NULL;
	}
    
    am_clusters = tmp_iter;
    for(int i = 0;i < 8;++i) peak_prob_file << -99999 << "\t";
    peak_prob_file << -99999 << "\n";
    peak_prob_file.flush();
    
    delete[] smallest_dist_arr;
}

void e0_optimization::single_point_cluster(){
	int* sub_cluster = new int[1];
	sub_cluster[0] = {1};
	
	double prob_tmp = do_calculations(sub_cluster,1,0,0);
	//prob_tmp /= save_sum;
	if(!thresholds->compare(prob_tmp,save_sum,1) && false){
		am_clusters_old = 1;
		am_clusters = 0;
		re_clustered_clusters = NULL;
		lens_of_clusters = NULL;
		if(verbosity > 0 && thr_num == 1){
			cout << "~> deleted 1 out of 1 clusters" << endl;
		}
		return;
	}
	lens_of_clusters = new int[am_clusters];
	am_clusters = 1;
	am_clusters_old = 1;
	lens_of_clusters[0] = 1;
	re_clustered_clusters = new double**[am_clusters];
	double edep_tmp = 0;
	
	double* smallest_dist_arr = new double[3];
    double smallest_dist = 0;
	
	
    for(int i = 0;i < am_clusters;i++){
    re_clustered_clusters[i] = new double*[lens_of_clusters[i]];    
		for(int j = 0;j < lens_of_clusters[i];j++){
			edep_tmp += clusters[internal_counter][j][0];
			re_clustered_clusters[i][j] = new double[4];
			for(int k = 0;k < 4;k++){
				re_clustered_clusters[i][j][k] = clusters[internal_counter][j][k];
			}
			if(j == 0){
				for(int k = 0;k < 3;++k) smallest_dist_arr[k] = clusters[internal_counter][sub_cluster[j]-1][k+1];
				smallest_dist = prob->get_source_to_first(smallest_dist_arr,true);
			}
		}
		peak_prob_file << edep_tmp << "\t" << prob_tmp << "\t" << save_sum <<"\t";
		for(int s = 0;s < 4;++s) peak_prob_file << l2_etc[s] << "\t";
		peak_prob_file << smallest_dist << "\t" << 1 << "\n";
	}	
	for(int i = 0;i < 8;++i) peak_prob_file << -99999 << "\t";
	peak_prob_file << -99999 << "\n";
	peak_prob_file.flush();
	
	delete[] sub_cluster;
	delete[] smallest_dist_arr;
}

inline double e0_optimization::do_calculations(int* sub_cluster,int sub_cluster_width,int k,int pos_in_c,bool SAVE){

    bool first_perm = false;
    double Edep_sum = 0.;
    double** data_arr = new double*[sub_cluster_width];


    for(int i = 0;i < sub_cluster_width;i++){
        data_arr[i] = new double[4];
        Edep_sum += clusters[internal_counter][sub_cluster[i]-1][0];
    }

    for(int j = 0;j < sub_cluster_width;j++){
        for(int k = 0;k < 4;k++){
            data_arr[j][k] = clusters[internal_counter][sub_cluster[j]-1][k];
        }
    }
    if(Edep_sum > 10000) return 0;
    prob->set_data(data_arr,sub_cluster_width);

    double ret_val = e0_loop(Edep_sum,sub_cluster_width,SAVE);
    for(int i = 0;i < sub_cluster_width;i++) delete[] data_arr[i];
    delete[] data_arr;

    return ret_val;
}


inline double e0_optimization::do_calculations(int* sub_cluster,int sub_cluster_width,int k,int pos_in_c){

    if(CL_Perms->already_processed(k,pos_in_c)) return CL_Perms->get_already_processed_prob();

    double Edep_sum = 0.;
    double** data_arr = new double*[sub_cluster_width];

    //bool CL_Perms->check_already_processed();

    for(int i = 0;i < sub_cluster_width;i++){
        data_arr[i] = new double[4];
        Edep_sum += clusters[internal_counter][sub_cluster[i]-1][0];
    }

    for(int j = 0;j < sub_cluster_width;j++){
        for(int k = 0;k < 4;k++){
            data_arr[j][k] = clusters[internal_counter][sub_cluster[j]-1][k];
        }
    }
    if(Edep_sum > 10000) return 0;
    prob->set_data(data_arr,sub_cluster_width);

    double ret_val = e0_loop(Edep_sum,sub_cluster_width);
    
    for(int i = 0;i < sub_cluster_width;i++) delete[] data_arr[i];
    delete[] data_arr;
    
    CL_Perms->set_probability_to_sub_cluster(ret_val,k,pos_in_c);

    return ret_val;
}

double e0_optimization::e0_loop(double Edep_sum,int sub_cluster_width,bool SAVE){
	return 1;
    double val = 0.;
	int start = start_int(Edep_sum);

    int looper = 0;
    double sum = 0;
    double peak = 0;
    double non_peak = 0.;
	
    l2 = 0;
    while(e0_coarse[looper] < 9998){
        integral_coarse[looper] = prob->set_e0(e0_coarse[looper],leading_order);
        integral_data->save_integral(e0_coarse[looper],integral_coarse[looper]);
        
        if(e0_coarse[looper] - Edep_sum > -2 && e0_coarse[looper] - Edep_sum < 38 && integral_coarse[looper] > 0){
			peak += integral_coarse[looper];
		}
	
        else if(looper - start > 40) non_peak += integral_coarse[looper];

        sum += integral_coarse[looper];
        l2 += pow(integral_coarse[looper],2);
        
        looper += 1;
    }
   
    if(leading_order && verbosity > 1){
		integral_data->save_l2_norm(Edep_sum,cluster_len[internal_counter],l2,peak,sum);
        integral_data->next();
    }
    if(sum == 0 || l2 == 0) val = 0;
    else val = peak/sum;

    reset();
    return val;
}


double e0_optimization::e0_loop(double Edep_sum,int sub_cluster_width){
    
    double peak_range[2] = {0};
    if(solve_type_bool){
		peak_range[0] = -40;
		peak_range[1] = 40;
	}
	else{
		peak_range[0] = -5;
		peak_range[1] = 85;
	}
    
    double val = 0.;
	int start = start_int(Edep_sum);
    int looper = start ;
    int max_coarse_int = 0;
    double max_coarse = 0.;
    double sum = 0;
    double peak = 0;
    double non_peak = 0.;
    int am_non_zero = 0;
    double ent = 0.;
    double ent_tot = 0;
    double l2_peak;
    l2 = 0;
    bool peaks_true = false;
	//cout << "START @ " << start << " -> " << e0_coarse[start] << " " << Edep_sum << endl;
    if(leading_order){
		if(abs(Edep_untorn - 1172.3) < 6 || abs(Edep_untorn - 1332.5) < 6) peaks_true = true;
	}
    while(e0_coarse[looper] < 9998){
        integral_coarse[looper] = prob->set_e0(e0_coarse[looper],leading_order);
        if(verbosity > 1 && thr_num == 1) integral_data->save_integral(e0_coarse[looper],integral_coarse[looper]);
        
        if(e0_coarse[looper] - Edep_sum > peak_range[0] && e0_coarse[looper] - Edep_sum < peak_range[1] && integral_coarse[looper] > 0){
			peak += integral_coarse[looper];
			ent -= integral_coarse[looper]*log(integral_coarse[looper]);
			l2_peak += integral_coarse[looper]*integral_coarse[looper];
		}
        else if(e0_coarse[looper] - Edep_sum > 20) non_peak += integral_coarse[looper];
        
        if(integral_coarse[looper] > 0){
			am_non_zero++;
			ent_tot -= integral_coarse[looper]*log(integral_coarse[looper]);
		}
        sum += integral_coarse[looper];
        l2 += pow(integral_coarse[looper],2);

        if(integral_coarse[looper] > max_coarse){
            max_coarse_int = looper;
            max_coarse = integral_coarse[looper];
        }
        looper += 1;
    }
    
    
    if(leading_order && peaks_true) sum_data->save_sum_val(Edep_untorn,sum,sqrt(l2));
    if(verbosity > 0  && thr_num == 1 && false){
		cout <<"----\n";
		cout <<"Hs peak "<<ent << " | l2 " << sqrt(l2)   << " | Hs_tot "<< ent_tot << endl;
		cout << "Hs_peak/l2 " << ent/sqrt(l2) <<" | Hs_peak/Hs_tot " << ent/ent_tot << endl;
		cout << "peak prob " << peak << " | total prob " << sum  << " | peak/total " << peak/sum << endl;
		cout <<"----\n";
	}
	bool esc = false;
	if(!leading_order && am_non_zero < 20) esc = true;
	
    if(sum == 0 || l2 == 0 || esc){
		reset();
		return 0;
	}
    //if(robust_cluster) val = peak;
    //else{
    if(leading_order){
		tmp_peak[tmp_peak_iter] = peak;
		tmp_peak_iter += 1;
	}
	if(!robust_cluster && abs(peak/sum - 1) < 1e-4) val = 0;
	//else 
	//if(!robust_cluster) 
	val = peak;
	save_sum = sum;
	
	l2_etc[0] = l2_peak;
	l2_etc[1] = l2;
	l2_etc[2] = ent;
	l2_etc[3] = ent_tot;
	///sum*prior_e0;
	//sqrt(peak)/sqrt(l2);
		//val = log(peak) + log(binning) - log(10000.) - log(sum);
	//else val = peak/sqrt(l2);
   
    reset();
    return val;
}

void e0_optimization::set_e0s(){
    integral_coarse = new double[10000];
    e0_coarse = new double[10000];
	prior_e0 = 1./(10000./binning);
    for(int i = 0;i < 10000;i++){
        e0_coarse[i] = ((double) i+1)*(Emax_coarse - Emin_coarse)/(10000./binning);
        integral_coarse[i] = 0.;
    }
}

void e0_optimization::reset(){for(int i = 0;i < 10000;i++) integral_coarse[i] = 0.;}

int e0_optimization::start_int(double E){
    double value = 0.98*E/binning;
    return (int) value-1;
}

int e0_optimization::am_cluster(){return am_clusters;}

int e0_optimization::get_cluster_len(int i){
	//if(verbosity > 0) cout << "-> length cluster " << i+1 << ": " << lens_of_clusters[i] << endl;
	return lens_of_clusters[i];
}

double** e0_optimization::get_cluster(int i){return re_clustered_clusters[i];}

void e0_optimization::set_am_clusters(int am_clusters){this->am_clusters = am_clusters;}

inline int e0_optimization::factorial(int a){
    int s = 1;
    for(int i = 1;i <= a;i++) s *= i;
    return s;
}
