#include "cluster_perms.h"

using namespace std;

cluster_perms::cluster_perms(int verbosity,int thr_num){
    
    this->thr_num = thr_num;
    this->verbosity = verbosity;
    internal_perm_counter = 0;

    A = new perm_class***[10];
    B = new merge_class**[10];

    remain_len_arr = new int[10];

    loaded_perms = new permutations();

    all_possibilities = new int**[10];
    for(int i = 0;i < 10;i++){
        all_possibilities[i] = new int*[100];
        for(int j = 0;j < 100;j++){
            all_possibilities[i][j] = new int[i+1];
            for(int k = 0;k < i+1;k++) all_possibilities[i][j][k] = 0;
        }
    }
    for(int i = 1;i <= 10;i++) set_all_possibilities(i);

    prob_tmp = 0;

    //test_scope();
}

cluster_perms::~cluster_perms(){
    for(int k = 0;k < 10;k++){
        for(int i = 0;i < remain_len_arr[k];i++){
            for(int j = 0;j <= k;j++) if(A[k][i][j] != NULL) delete A[k][i][j];
            delete[] A[k][i];
            delete B[k][i];
        }
        delete[] B[k];
        delete[] A[k];
    }
    delete[] A;
    delete[] B;
    for(int i = 0;i < 10;i++){
        for(int j = 0;j < 100;j++) delete[] all_possibilities[i][j];
        delete[] all_possibilities[i];
    }
    delete[] all_possibilities;
    delete loaded_perms;
    delete[] remain_len_arr;
}

void cluster_perms::test_scope(){
    int iter = 1;
    int remain_len;
    while(iter < 11){
        remain_len = loaded_perms->return_remains_len(iter);
        for(int i = 0;i < remain_len;i++){
            B[iter-1][i]->compare_true();
            cout << "================================\n";
        }
        iter += 1;
    }
}

void cluster_perms::get_perms_arr(int* arr,int iter){
    
    auto factorial = [](int a){
        int ret_val = 1;
        for(int i = 1;i <= a;i++) ret_val *= i;
        return ret_val;
    };

    for(int i = 1;i <= iter;i++) arr[i-1] = factorial(iter)/(factorial(i)*factorial(iter-i));

}

void cluster_perms::set_all_possibilities(int iter){
	int remain_len = 0;
    int** remains = loaded_perms->return_remains(iter);
    remain_len = loaded_perms->return_remains_len(iter);
    if(iter == 1) remain_len += 1;
    //cout << iter <<" -> !remain_len! " << remain_len << endl;
    remain_len_arr[iter-1] = remain_len;

    int* am_clusters = loaded_perms->return_am_clusters(iter);

    A[iter-1] = new perm_class**[remain_len];
    
    for(int i = 0;i < remain_len;i++){
        A[iter-1][i] = new perm_class*[iter];
        for(int j = 0;j < iter;j++){
            if(remains[i][j] > 0) A[iter-1][i][j] = new perm_class(remains[i][j],iter);
            else A[iter-1][i][j] = NULL;
        }
    }
    
    int** all_perms = loaded_perms->get_perms_full(iter);
    int am_perms_complete = 1;
    for(int i = 1;i <= iter;i++) am_perms_complete *= i;

    int start_pos = 0;
    for(int i = 0;i < remain_len;i++){
        for(int j = 0;j < iter;j++){
            if(remains[i][j] > 0){
                if(j == 0){
					A[iter-1][i][j]->set_id(iter,remains[i][j],0,i);
					A[iter-1][i][j]->feed_complete(all_perms,am_perms_complete,0);
				}
                else{
                    start_pos += remains[i][j-1];
                    A[iter-1][i][j]->set_id(iter,remains[i][j],start_pos,i);
                    A[iter-1][i][j]->feed_complete(all_perms,am_perms_complete,start_pos);
                }
                
            }
        }
        start_pos = 0;
    }
	//cout << "save for " << iter << " complete " << endl;
    B[iter-1] = new merge_class*[remain_len];
    for(int i = 0;i < remain_len;i++){
        B[iter-1][i] = new merge_class(remains[i],am_clusters[i],am_perms_complete,A[iter-1][i],iter);
        //B[iter-1][i]->compare_true();
        //cout << "================================\n";
    }
    //delete[] am_perms_arr;
}

//void cluster_perms::check_already_processed(){
    
//}

void cluster_perms::set_probability_to_sub_cluster(double val,int k,int j){
    return;
    A[active_cluster_int-1][active_sub_cluster_int-1][k]->set_probability(val,j);
}

bool cluster_perms::already_processed(int k,int j){
    /*if(A[active_cluster_int-1][active_sub_cluster_int-1][k]->check_recurrent(j)){
        prob_tmp = A[active_cluster_int-1][active_sub_cluster_int-1][k]->get_probability(j);
        return true; 
    }*/
    return false;
}

double cluster_perms::get_already_processed_prob(){return prob_tmp;}

void cluster_perms::set_active_cluster(int iter){
    active_cluster = B[iter-1];
    active_cluster_int = iter;
}

void cluster_perms::set_sub_cluster(int iter,bool p_val){
    active_sub_cluster = active_cluster[iter-1];
    active_sub_cluster_int = iter-1;
    
    int* remains = loaded_perms->return_remains(active_cluster_int)[iter-1];
    if(verbosity > 0 && p_val && thr_num == 1){
		cout << "\n***** active remains cluster *****" << endl;
		for(int i = 0;i < active_cluster_int;i++) cout << remains[i] << " ";
		cout << "\n**********************************\n" << endl;
	}
}

void cluster_perms::reset_probabilities(int iter){for(int i = 0;i < remain_len_arr[iter-1];i++) B[iter-1][i]->reset_probability();}

void cluster_perms::assign_prob_to_perm(double prob,int j){active_sub_cluster->set_probability(prob,j);}

void cluster_perms::reset_active_cluster_probs(){active_sub_cluster->reset_probability();}

double* cluster_perms::get_highest_prob_cluster(){return active_sub_cluster->get_highest_prob_cluster();}

int* cluster_perms::return_active_sub_cluster_perm(int j,int k){return active_sub_cluster->return_sub_cluster_perm(k,j);}

int cluster_perms::return_am_cluster_groups_in_active(){return loaded_perms->return_remains_len(active_cluster_int);}

int cluster_perms::return_active_sub_cluster_width(int k){return active_sub_cluster->get_cluster_width(k);}

int cluster_perms::return_active_sub_cluster_amount(){return loaded_perms->return_am_clusters(active_cluster_int,active_sub_cluster_int);}
