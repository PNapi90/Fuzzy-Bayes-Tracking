#include <iostream>

#include "fuzzy_c_means.h"


using namespace std;

fuzzy_c_means::fuzzy_c_means(double fuzz,int verbosity,int thr_num){
	
	this->verbosity = verbosity;
	this->thr_num = thr_num;
    fuzzyness = fuzz;
    epsilon = 10;
    set_dmax();

    weights = new double*[max_len];
    centroids = new double*[max_len]; 
    dist_2 = new double*[max_len];
    x = new double*[max_len];
    old_weights = new double*[max_len];

    for(int i = 0;i < max_len;i++){
        weights[i] = new double[max_len];
        centroids[i] = new double[3];
        x[i] = new double[4];
        dist_2[i] = new double[max_len];
        old_weights[i] = new double[max_len];
    }

    cluster_to_points = new int[max_len];
    cluster_len = new int[max_len];
    for(int i = 0;i < max_len;i++) cluster_len[i] = 0;
    
    fc_clusters = new double**[max_len];
    for(int i = 0;i < max_len;i++){
        fc_clusters[i] = new double*[max_len];
        for(int j = 0;j < max_len;j++) fc_clusters[i][j] = new double[4];
    }
}

fuzzy_c_means::~fuzzy_c_means(){
    for(int i = 0;i < max_len;i++){
        for(int j = 0;j < max_len;j++) delete[] fc_clusters[i][j];
        delete[] fc_clusters[i];
        delete[] weights[i];
        delete[] centroids[i];
        delete[] dist_2[i];
        delete[] x[i];
        delete[] old_weights[i];
    }
    delete[] fc_clusters;
    delete[] weights;
    delete[] old_weights;
    delete[] centroids;
    delete[] dist_2;
    delete[] x;
    delete[] cluster_to_points;
    delete[] cluster_len;
}



void fuzzy_c_means::init_weights(){
    
    double norm = 0.;
    for(int i = 0;i < am_clusters;i++){
        for(int k = 0;k < am_points;k++) weights[i][k] = ((double) (rand() % 1000 + 1))/1000.;
    }
    for(int k = 0;k < am_points;k++){
		for(int i = 0;i < am_clusters;i++) norm += weights[i][k];
		for(int i = 0;i < am_clusters;i++){
			weights[i][k] /= norm;
			old_weights[i][k] = weights[i][k];
			//cout << weights[i][k] << " ";
		}
		//cout << endl;
		norm = 0;
	}
}

void fuzzy_c_means::fuzzy_clustering(double** data,int len,int am_clusters){
    this->am_clusters = am_clusters;
    sum_edep = 0;
    if(len == 1){
        am_clusters = 1;
        am_points = 1;
        cluster_len[0] = 1;
        for(int k = 0;k < am_points;k++){
            for(int j = 0;j < 4;j++) fc_clusters[0][k][j] = data[k][j];
        }
        if(verbosity > 0 && thr_num == 1) print_info();
        return;
    }
    for(int i = 0;i < len;i++){
        sum_edep += data[i][0];
        for(int j = 0;j < 4;j++) x[i][j] = data[i][j];
    }

    if(verbosity > 0 && thr_num == 1){
		cout << "\n###################################" << endl;
		for(int i = 0;i < len;i++){
			for(int j = 0;j < 4;j++) cout << x[i][j] << "\t";
			cout << endl;
		}
		cout << "###################################" << endl;
		cout << "------------------------------\n";
		cout << "Edep in Cluster: " << sum_edep << " keV";
		cout << "\n------------------------------\n" << endl;
    }
    
    am_points = len;
    //get_am_clusters();
    init_weights();
    first_iter = true;
    epsilon = 10;
    //cout << abort_val << " " << epsilon << endl;
    while(epsilon > abort_val){
        calc_centroids();
        calc_dist2();
        calc_weights();
        calc_epsilon();
        adjust_old_weights();
        //cout << abort_val << " " << epsilon << endl;
    }

    set_cluster_to_points();
}

void fuzzy_c_means::adjust_old_weights(){
    for(int i = 0;i < am_clusters;i++){
        if(!first_iter) for(int k = 0;k < am_points;k++) old_weights[i][k] = weights[i][k];
        else for(int k = 0;k < am_points;k++) old_weights[i][k] = 0.;
    }
}

void fuzzy_c_means::calc_centroids(){
    double sum_weights;
    for(int i = 0;i < am_clusters;i++){
        sum_weights = 0.;
        for(int j = 0;j < 3;j++){
            centroids[i][j] = 0.;
            for(int k = 0;k < am_points;k++){
                centroids[i][j] += pow(weights[i][k],fuzzyness)*x[k][j+1];
                if(j == 0) sum_weights += pow(weights[i][k],fuzzyness);
            }
            centroids[i][j] /= sum_weights;
        }
    }
}

void fuzzy_c_means::calc_dist2(){
    for(int i = 0;i < am_clusters;i++){
        for(int k = 0;k < am_points;k++){
            dist_2[i][k] = 0.;
            for(int j = 0;j < 3;j++) dist_2[i][k] += pow(x[k][j+1] - centroids[i][j],2.);
        }
       
    }
}

void fuzzy_c_means::calc_weights(){
    double sum_val;
    double exponent = 1./(fuzzyness - 1.);
    for(int i = 0;i < am_clusters;i++){
        for(int k = 0;k < am_points;k++){
            sum_val = 0.;
            for(int j = 0;j < am_clusters;j++){ 
				if(j != i) sum_val += pow(dist_2[i][k]/dist_2[j][k],exponent);
				else sum_val += 1;
			}
            weights[i][k] = 1./sum_val;
        }
    }
}

void fuzzy_c_means::calc_epsilon(){
    if(first_iter){
        first_iter = false;
        //return;
    }
    double norm_val = 0.;
    double max = 0;
    for(int i = 0;i < am_clusters;i++){
        for(int k = 0;k < am_points;k++){
			norm_val = abs(weights[i][k] - old_weights[i][k]);
			if(max < norm_val) max = norm_val;
			//cout << norm_val << endl;
		}
    }
    epsilon = max;
}

void fuzzy_c_means::reset(){
    epsilon = 10.;
    for(int i = 0;i < max_len;i++){
        for(int j = 0;j < 3;j++){
            centroids[i][j] = 0.;
            x[i][j] = 0.;
        }
        x[i][3] = 0.;
        cluster_to_points[i] = -3;
        cluster_len[i] = 0;
    }
}

void fuzzy_c_means::set_cluster_to_points(){
    //cout <<"do the thing " << am_points <<" "<<am_clusters<< endl;
    double max_weight = 0.;
    int max_cl_int = 0;
    
    for(int k = 0;k < am_points;k++){
        max_weight = 0.;
        max_cl_int = 0;
        for(int i = 0;i < am_clusters;i++){
            if(max_weight < weights[i][k]){
                max_weight = weights[i][k];
                max_cl_int = i;
            }
        }
        cluster_to_points[k] = max_cl_int;
        for(int j = 0;j < 4;j++){
            fc_clusters[max_cl_int][cluster_len[max_cl_int]][j] = x[k][j];
        }
        cluster_len[max_cl_int]++;
    }
    sort_clusters();
    if(verbosity > 0 && thr_num == 1) print_info();
    
}

void fuzzy_c_means::sort_clusters(){
	double*** tmp_cluster = new double**[am_clusters];
    int* tmp_len = new int[am_clusters];
    for(int i = 0;i < am_clusters;++i){
		tmp_len[i] = 0;
		tmp_cluster[i] = new double*[am_points];
		for(int j = 0;j < am_points;++j){
			tmp_cluster[i][j] = new double[4];
			for(int k = 0;k < 4;++k) tmp_cluster[i][j][k] = 0;
		}
	}
    
    int iter = 0;
    for(int i = 0;i < am_clusters;++i){
        if(cluster_len[i] > 0){
            for(int j = 0;j < cluster_len[i];++j){
                for(int k = 0;k < 4;++k) tmp_cluster[iter][j][k] = fc_clusters[i][j][k];
            }
            tmp_len[iter] = cluster_len[i];
            iter += 1;
        }
    }
    double edep_tmp = 0;
    for(int i = 0;i < iter;++i){
        edep_tmp = 0;
        for(int j = 0;j < tmp_len[i];++j){
            for(int k = 0;k < 4;++k) fc_clusters[i][j][k] = tmp_cluster[i][j][k];
            edep_tmp += fc_clusters[i][j][0];
        }
        cluster_len[i] = tmp_len[i];
        //Edep_arr_2[i] = edep_tmp;
    }
    am_clusters = iter;
    for(int i = 0;i < am_clusters;++i){
		for(int j = 0;j < am_points;++j) delete[] tmp_cluster[i][j];
		delete[] tmp_cluster[i];
	}
	delete[] tmp_cluster;
	delete[] tmp_len;
	
}

void fuzzy_c_means::print_info(){
	cout << "*** Amount of fuzzy clusters: " << am_clusters << " ***" << endl;
    cout << "= = = = = = = = = = = = = = = = = =" << endl;
    for(int i = 0;i < am_clusters;i++){
		cout << "*************************\n";
        for(int k = 0;k < cluster_len[i];k++){
            for(int j = 0;j < 4;j++){
                cout << fc_clusters[i][k][j] << " ";
            }
            cout << endl;
        }
        cout << "*************************\n";
        cout << "--------------------------\n";
        cout << "len of cluster: " << cluster_len[i];
        cout << "\n--------------------------\n";
    }
    cout << "= = = = = = = = = = = = = = = = = =" << endl;
}

void fuzzy_c_means::set_dmax(){
    dmax = DMAX;
}

void fuzzy_c_means::get_am_clusters(){
    am_clusters = 1;
    
    double** distances_tmp = new double*[am_points];
    for(int i = 0;i < am_points;i++){
		distances_tmp[i] = new double[am_points];
		for(int j = 0;j < am_points;j++) distances_tmp[i][j] = 0;
	}
	for(int i = 0;i < am_points;i++){
		cout << i+1 << " -> ";
		for(int j = 0;j < am_points;j++){
			if(i >= j||true){
				for(int k = 0;k < 3;k++){
					distances_tmp[i][j] += pow(x[i][k] - x[j][k],2);
				}
				distances_tmp[i][j] = sqrt(distances_tmp[i][j]);
				cout << distances_tmp[i][j] << "\t";
			}
			
		}
		cout << endl;
	}
	/*
    double norm_val = 0.;
    int iter = 0;
    if(am_points == 4){
		while(iter < 3){	
			for(int k = iter;k < am_points;k++){
				if(k != iter) for(int j = 0;j < 3;j++) norm_val += pow(x[iter][j+1] - x[k][j+1],2.);
			}
			iter += 1;
			norm_val = sqrt(norm_val);
			norm_val = 0;
			if(norm_val > dmax) am_clusters += 1;
		}
	}
    */
    for(int i = 0;i < am_clusters;i++) cluster_len[i] = 0;
    //am_clusters = 2;
    
    for(int i = 0;i < am_points;i++) delete[] distances_tmp[i];
	delete[] distances_tmp;
}

int fuzzy_c_means::return_am_clusters(){
	int amount = 0;
	for(int i = 0;i < am_clusters;++i) if(cluster_len[i] > 0) amount += 1;
	return amount;	
}

int fuzzy_c_means::return_cluster_len(int pos){return cluster_len[pos];}

double** fuzzy_c_means::return_cluster(int pos){return fc_clusters[pos];}

double fuzzy_c_means::get_edep_in_cluster(){return sum_edep;}
