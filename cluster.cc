#include "cluster.h"

using namespace std;

cluster_class::cluster_class(int thr_num,int verbosity,bool solve_type_bool){
	this->thr_num = thr_num;
	this->solve_type_bool = solve_type_bool;

	//if(this->thr_num != verbose_thread_num) this->verbosity = 0;
	//else
	this->verbosity = verbosity;
	
	air = new air_path(this->verbosity,this->thr_num);
	RL_FCM = new RLFCM();
    fc_mean = new fuzzy_c_means(1.5,this->verbosity,thr_num);
    e_zero = new e0_optimization(this->thr_num,this->verbosity,this->solve_type_bool);

    am_cluster = 0;
    cluster_lens = new int[30];

    clusters = new double**[10];
    for(int i = 0;i < 10;i++){
        clusters[i] = new double*[30];
        for(int j = 0;j < 30;j++){
            clusters[i][j] = new double[4];
            for(int k = 0;k < 4;k++) clusters[i][j][k] = 0.;
        }
    }

    TS_cluster = new double*[30];
    for(int i = 0;i < 30;i++) TS_cluster[i] = new double[4];
    first_cluster = true;
    gamma_counter = 0;
    am_torn = 0;

}

cluster_class::~cluster_class(){
    for(int i = 0;i < 10;i++){
        for(int j = 0;j < 30;j++) delete[] clusters[i][j];
        delete[] clusters[i];
    }
    delete[] clusters;
    delete[] cluster_lens;

    for(int i = 0;i < 30;i++) delete[] TS_cluster[i];
    delete[] TS_cluster;
    //cout << "************************************\n";
    delete fc_mean;
    delete RL_FCM;
    //cout << " -> deleted fuzzy_c_means objects" << endl;
    delete air;
    //cout << " -> deleted cluster decision object" << endl;
    delete e_zero;
    //cout << " -> deleted e0_optimization object" << endl;
    //cout << " -> deleted cluster object" << endl;
    //cout << "************************************\n\n";
    cout << "\n===============\n";
    cout << "thread " << thr_num << " done";
    cout << "\n===============\n" << endl;
    //cout << "* * * all clusters processed * * *\n\n";
    //cout << "=> moving to BAYES-TRACKING of Compton-escaped gammas ... \n" << endl;
}
inline int cluster_class::threshold(int c_am,int len){

	double root = sqrt((double) len);
	double remain = root - ((int) root);
	if(c_am > root){
		if(remain >= 0.5) return (int) root + 1;
		return (int) root; 
	}
	return c_am;
	
}


void cluster_class::clustering(double** x,int len){
	
	int full_cluster_iter = 0;
	
	//reset arrays in e0_optimization
    if(!first_cluster) e_zero->reset_re_clustered();
    else first_cluster = false;
	
    for(int i = 0;i < len;i++) for(int j = 0;j < 4;j++) TS_cluster[i][j] = x[i][j];
    //pre cluster using fuzzy c means
    if(verbosity > 0 && thr_num == 1) cout << "* * * * * Fuzzy pre-clustering * * * * *" << endl;
    double etmp = 0;
    for(int i = 0;i < len;i++) etmp += x[i][0];
    bool air_passage = air->check_air(x,len);
    int RL_FCM_clusters = 1;
    
    if(!air_passage && len > 2){
		RL_FCM->robust_clustering(x,len);
		RL_FCM_clusters = threshold(RL_FCM->get_am_clusters(),len);
    }
    fc_mean->fuzzy_clustering(x,len,RL_FCM_clusters);
    e_zero->set_am_clusters(len);
    
    bool tmp_b = false;
    double edep_comp = fc_mean->get_edep_in_cluster();
    if(abs(edep_comp - 1173.2) < 10) tmp_b = true;
    if(abs(edep_comp - 1332.5) < 10) tmp_b = true;
    
    int full_am_clusters = 0;
    double** dummy_cluster = NULL;
    
    double ed_temp = 0;
	bool already_counted = false;
	
    //validate and re-cluster using e0-optimization method
    for(int c = 0;c < fc_mean->return_am_clusters();c++){
		e_zero->re_cluster(fc_mean->return_cluster(c),fc_mean->return_cluster_len(c));
		
		if(e_zero->bad_cluster() && fc_mean->return_am_clusters() == 1 && false){
		    e_zero->reset_re_clustered();
		    fc_mean->reset();
		    fc_mean->fuzzy_clustering(x,len,RL_FCM_clusters+1);
            e_zero->set_am_clusters(len);
            c--;
            cout << "Bad cluster encoutered -> Refuzzy!" << endl;
            continue;
		}
		//extract clusters out of e_zero class
		am_cluster = e_zero->am_cluster();
		full_am_clusters += am_cluster;
		
		if(verbosity > 0 && thr_num == 1) cout <<"re_clustered -> amount of new clusters: "<< am_cluster << endl;
    
		if(am_cluster == 0){
			cerr << "zero in thr " << thr_num << endl;
			cerr << "at gamma # " << gamma_counter << endl;
		}
		ed_temp = 0;
		already_counted = false;
		
		for(int i = 0;i < am_cluster;i++){
			ed_temp = 0;
			cluster_lens[full_cluster_iter] = e_zero->get_cluster_len(i);
			if(verbosity > 0 && thr_num == 1){
				cout <<"cluster len "<< cluster_lens[full_cluster_iter];
				cout <<"\n= = = = = = = = = = = = = = = = = = = = = = ="<<endl;
			}
			dummy_cluster = e_zero->get_cluster(i);
			for(int j = 0;j < cluster_lens[full_cluster_iter];j++){
				for(int k = 0;k < 4;k++) clusters[full_cluster_iter][j][k] = dummy_cluster[j][k];
			}
			for(int j = 0;j < cluster_lens[full_cluster_iter];j++){
				ed_temp += dummy_cluster[j][0];
				if(verbosity > 0 && thr_num == 1){
					for(int k = 0;k < 4;k++){
						cout << dummy_cluster[j][k] << "\t";
					}
					cout << endl;
				}
			}
			if(tmp_b && check_etmp(ed_temp) && !already_counted){
				am_torn += 1;
				already_counted = true;
			}
			if(verbosity > 0 && thr_num == 1){
				cout << "----------------\n";
				cout << "Edep: " << ed_temp << " keV";
				cout << "\n----------------\n";
				cout <<"= = = = = = = = = = = = = = = = = = = = = = =\n"<<endl;
			}
			full_cluster_iter += 1;
		}
	}    
	//if(am_cluster > 1 && tmp_b) am_torn += 1; 
    am_cluster = full_am_clusters;
    gamma_counter += 1;
}

bool cluster_class::check_for_air(double** x,int len){
	double max_dist = 0;
	double** dist_mat = new double*[len];
	for(int i = 0;i < len;i++){
		dist_mat[i] = new double[len];
		for(int j = 0;j < len;j++){
			dist_mat[i][j] = 0;
			for(int k = 0;k < 3;k++) dist_mat[i][j] += pow(x[i][k+1] - x[j][k+1],2);
			dist_mat[i][j] = sqrt(dist_mat[i][j]);
			if(verbosity > 0 && thr_num == 1) cout << dist_mat[i][j] << " ";
			if(max_dist < dist_mat[i][j]) max_dist = dist_mat[i][j];
		}
		if(verbosity > 0 && thr_num == 1) cout << endl;
	}
	
	if(max_dist > 70 && max_dist < 100) return false;
	return true;
	
}

bool cluster_class::check_etmp(double e){
	if(abs(e - 1172.3) < 10) return false;
	if(abs(e - 1332.5) < 10) return false;
	return true; 
}

void cluster_class::reset(){
    for(int i = 0;i < am_cluster;i++){
        for(int j = 0;j < cluster_lens[i];j++){
            for(int k = 0;k < 4;k++) clusters[i][j][k] = 0.;
        }
        cluster_lens[i] = 0;
    }
    am_cluster = 0;
    fc_mean->reset();
}

double** cluster_class::get_cluster(int i){return clusters[i];}

int cluster_class::amount(){return am_cluster;}

int cluster_class::cluster_len(int pos){return cluster_lens[pos];}

int cluster_class::get_am_torn(){return am_torn;}
