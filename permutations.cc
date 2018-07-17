#include <iostream>

#include "permutations.h"

using namespace std;

permutations::permutations(){
    load_poss_combs();
    arr = new int*[100000];
    for(int i = 0;i < 100000;i++) arr[i] = new int[20];

    remain_arr = new int**[10];
    remain_lens = new int[10];
	for(int i = 0;i < 10;i++) remain_lens[i] = 0;
    am_clusters = new int*[10];
    for(int i = 0;i < 10;i++){
        am_clusters[i] = new int[500];
        for(int j = 0;j < 500;j++) am_clusters[i][j] = 0;
    }

    //cout << "Calculating remain array for clustering ...\t";
    for(int i = 0;i < 10;i++){
        remain_arr[i] = new int*[1000];
        for(int j = 0;j < 1000;j++){
			remain_arr[i][j] = new int[20];
			for(int k = 0;k < 20;k++) remain_arr[i][j][k] = 0;
		}
        if(i == 0) continue;
        calc_remain_combs(i+1);
        save_important_remains(i+1);
    }

    full_perm_arr = new int**[10];
    int f_len = 0;
    for(int i = 0;i < 10;i++){
        f_len = factorial(i+1);
        full_perm_arr[i] = new int*[f_len];
        for(int j = 0;j < f_len;j++) full_perm_arr[i][j] = new int[i+1];
    }

    for(int i = 0;i < 10;i++) set_full_perm_arr(i+1);

    //cout << "done\n" << endl;
}

permutations::~permutations(){
    for(int i = 0;i < max_interactions;i++){
        for(int j = 0;j <= i;j++){
            for(int k = 0;k < max_perms[i][j];k++) delete[] permutation_sub_arr[i][j][k];
            delete[] permutation_sub_arr[i][j];
        }
        delete[] permutation_sub_arr[i];
    }

    for(int i = 0;i < 10;i++){
        for(int j = 0;j < 1000;j++) delete[] remain_arr[i][j];
        delete[] remain_arr[i];
        delete[] am_clusters[i];
    }
    for(int i = 0;i < 100000;i++) delete[] arr[i];
	
    int f_len = 0;
    for(int i = 0;i < 10;i++){
		f_len = factorial(i+1);
        for(int j = 0;j < f_len;j++) delete[] full_perm_arr[i][j];
        delete[] full_perm_arr[i];
    }
    delete[] full_perm_arr;
	
    for(int i = 0;i < max_interactions;i++) delete[] max_perms[i];
    delete[] max_perms;
	
    delete[] am_clusters;
    delete[] arr;
    delete[] permutation_sub_arr;
    delete[] remain_lens;
    delete[] remain_arr;
}

void permutations::set_full_perm_arr(int iter){
    int sortarray[iter];
    int sub_iter = 0;
    for(int i = 0;i < iter;i++) sortarray[i] = i+1;
    sort(sortarray,sortarray+iter);
    do{
        for(int i = 0;i < iter;i++) full_perm_arr[iter-1][sub_iter][i] = sortarray[i];
        sub_iter += 1;
    }while(next_permutation(sortarray,sortarray+iter));
}


inline int permutations::factorial(int a){
    int ret_val = 1;
    for(int i = 1;i <= a;i++) ret_val *= i;
    return ret_val;
}

void permutations::load_poss_combs(){
    get_max_perms();

    permutation_sub_arr = new int***[max_interactions];
    for(int i = 0;i < max_interactions;i++){
        permutation_sub_arr[i] = new int**[i+1];
        for(int j = 0;j <= i;j++){
            permutation_sub_arr[i][j] = new int*[max_perms[i][j]];
            for(int k = 0;k < max_perms[i][j];k++) permutation_sub_arr[i][j][k] = new int[j+1];
        }
    }

    string name = "Perms_";
    string ending = ".dat";
    string folder;
    
    ifstream*** data_arr = new ifstream**[max_interactions];
    for(int i = 0;i < max_interactions;i++){
        data_arr[i] = new ifstream*[i+1];
        folder = "perm_list/ints/"+to_string(i+1)+"_int/";
        for(int j = 0;j < i+1;j++) data_arr[i][j] = new ifstream((folder+name+to_string(j+1)+ending).c_str());
    }
    
    for(int i = 0;i < max_interactions;i++){
        for(int j = 0;j <= i;j++){
            for(int k = 0;k < max_perms[i][j];k++){
                for(int l = 0;l <= j;l++){
                    *data_arr[i][j] >> permutation_sub_arr[i][j][k][l];
                }
            }
        }
    }

    for(int i = 0;i < max_interactions;i++){
        for(int j = 0;j <= i;j++) delete data_arr[i][j];
        delete[] data_arr[i];
    }
    delete[] data_arr;


}

void permutations::get_max_perms(){
    
    max_perms = new int*[max_interactions];

    for(int i = 0;i < max_interactions;i++){
        max_perms[i] = new int[i+1];
        for(int j = 0;j <= i;j++){
            max_perms[i][j] = k_out_n(i+1,j+1);
        }
    }
}

int*** permutations::get_perms(int interaction){
    return permutation_sub_arr[interaction-1];
}
int* permutations::return_max_perms(int interaction){
    return max_perms[interaction-1];
}

int** permutations::get_perms_full(int interaction){return full_perm_arr[interaction-1];}

inline int permutations::k_out_n(int n,int k){
    return factorial(n)/(factorial(k)*factorial(n-k));
}

void permutations::cout_perms(int interaction){
    cout << "in get perms with int " << interaction <<endl;
    for(int i = 0;i < interaction;i++){
        cout <<"draw "<<i+1<<" -> max perms: " << max_perms[interaction-1][i] << endl;
        for(int j = 0;j < max_perms[interaction-1][i];j++){
            for(int k = 0;k <= i;k++) cout << permutation_sub_arr[interaction-1][i][j][k] << " ";
            cout << "\t";
        }
        cout << endl;
    }
    cout << endl;
}



void permutations::calc_remain_combs(int balls){
    int am_bars = balls-1;
    int len = am_bars+balls;
    int count = 0;
    int x[len];
    int int_iter = 0;
    int a = 0;
    for(int i = 0;i < len;i++){
        if(i < am_bars) x[i] = 1;
        else x[i] = 0;
    }
    sort(x,x+len);
    int sum = 0;
    do{
        sum = 0;
        for(int i = 0;i < len;i++){
            if(x[i] == 0) sum++;
            else{
                arr[int_iter][a] = sum;
                sum = 0;
                a += 1;
                count += 1;
            }
            if(i == len-1 && count == am_bars){
                arr[int_iter][a] = sum;
            }
        }
        int_iter += 1;
        a = 0;
        count  = 0;
    }while(std::next_permutation(x,x+len));
    leni = int_iter;
}

void permutations::save_important_remains(int iterator){
    
    int old = 0;
    bool plotter = false;
    int true_ones = 0;
    for(int i = 0;i < leni;i++){
        plotter = false;
        old = arr[i][0];
        for(int j = 1;j < iterator;j++){
            if(old >= arr[i][j] && old > 0){
                plotter = true;
            }
            else if(old < arr[i][j]){
                plotter = false;
                break;
            }
            old = arr[i][j];
        }
        for(int j = 0;j < iterator;j++){
            if(plotter){
                remain_arr[iterator-1][true_ones][j] = arr[i][j];
            }
        }
        if(plotter) true_ones++;
    }
    

    remain_lens[iterator-1] = true_ones;
    
    for(int j = 0;j < remain_lens[iterator-1];j++){
        for(int k = 0;k < iterator;k++){
            if(remain_arr[iterator-1][j][k] > 0) am_clusters[iterator-1][j] += 1;
        }
    } 
}

int** permutations::return_remains(int i){return remain_arr[i-1];}

int permutations::return_remains_len(int i){return remain_lens[i-1];}

int* permutations::return_am_clusters(int i){return am_clusters[i-1];}

int permutations::return_am_clusters(int i,int j){return am_clusters[i-1][j];}
