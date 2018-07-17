#include <iostream>

#include "merge_class.h"

using namespace std;

merge_class::merge_class(int* lens,int am_clust,int am_perm,perm_class** A){
    this->A = A;
    this->am_clust = am_clust;
    this->lens = lens;
    this->am_perm = am_perm;

    full_len = 0;
    for(int i = 0;i < this->am_clust;i++) full_len += this->lens[i]; 

    true_cluster = new int*[1000];
    for(int i = 0;i < 1000;i++) true_cluster[i] = new int[full_len];

}

merge_class::~merge_class(){
    for(int i = 0;i < 1000;i++) delete true_cluster[i];
    delete[] true_cluster;
}


void merge_class::compare_true(){
    cout << "am_clust " << am_clust << endl;
    for(int i = 0;i < am_perm;i++){
            for(int j = 0;j < am_clust;j++){
                A[j]->cout_line(i);
                cout << "| ";
            }
        cout << endl;
    }
}

/*
void merge_class::compare_true(){
    cout << "am_clust " << am_clust << endl;

    int true_count = 0;

    bool* skip = new bool[am_perm];


    for(int i = 0;i < am_perm;i++){
        for(int j = 0;j < am_clust;j++){
            if(A[j]->test_for_double()) true_count++;
        }
        if(true_count == am_clust) skip[i] = true;
        else skip[i] = false;
        true_count = 0;
    }

    int*** all_perms_arr = new int**[am_clust];
    int* all_perms_lens = new int[am_clust];
    for(int i = 0;i < am_clust;i++){
        all_perms_arr[i] = A[i]->get_arr();
        all_perms_lens[i] = A[i]->get_len_x();
    }
    true_count = 0;

    bool same = false;
    for(int i = 0;i < am_perm;i++){
        if(!skip[i]){
            for(int j = 0;j < am_clust;j++){
                /*for(int s = 0;s < am_clust;s++){
                    if(s != j && all_perms_lens[s] == all_perms_lens[j] && all_perms_lens[s] > 1){
                        
                        for(int p = 0;p <= i;p++){
                            for(int o = 0;o < all_perms_lens[j];o++){
                                if(all_perms_arr[s][i][o] == all_perms_arr[j][p][o]) true_count++;
                            }
                            if(true_count == all_perms_lens[j]){
                                same = true;
                                break;
                            }
                        }
                    }
                    if(same) break;
                }
                if(same) break;
                A[j]->cout_line(i);
                cout << "| ";
            }
            //same = false;
            //true_count = 0;

        }
        cout << endl;
    }
    delete[] skip;
    delete[] all_perms_arr;
    delete[] all_perms_lens;
}

*/