#include "perm_class.h"
#include <iostream>

using namespace std;

perm_class::perm_class(int len_x){
    this->len_x = len_x;

    arr_perms = new int*[1000];
    for(int i = 0;i < 1000;i++) arr_perms[i] = new int[this->len_x];

    arr_true = new int*[1000];
    for(int i = 0;i < 1000;i++) arr_true[i] = new int[this->len_x];
    cout << "created perm_class obj with len_x " << len_x << endl;
}


perm_class::~perm_class(){
    for(int i = 0;i < 1000;i++){
        delete arr_perms[i];
        delete arr_true[i];
    }
    delete[] arr_true;
    delete[] arr_perms;
}



void perm_class::feed_complete(int** complete_perms,int am_perm,int len_x_all,int start_pos){
    this->len_x_all = len_x_all;
    this->complete_perms = complete_perms;
    for(int i = 0;i < am_perm;i++){
        for(int j = 0;j < len_x;j++) cout << complete_perms[i][j] <<" ";
        cout << endl;
    }
    this->start_pos = start_pos;
    this->am_perm = am_perm;
    for(int i = 0;i < this->am_perm;i++){
        for(int j = 0;j < len_x;j++) arr_perms[i][j] = this->complete_perms[i][j+this->start_pos];
    }
    iter = 0;
    true_len = 0;
}


bool perm_class::test_for_double(){
    int true_count = 0;
    bool return_val = false;
    for(int i = 0;i < iter;i++){
        for(int j = 0;j < len_x;j++){
            if(arr_perms[i][j] == complete_perms[iter][j+start_pos]) true_count++;
        }
        if(true_count == len_x){
            return_val = true;
            /*for(int j = 0;j < len_x;j++){
                already[iter_a][j] = arr_perms[i][j];
            }
            iter_a++;*/
            break;
        }
        true_count = 0;
    }
    if(!return_val){
        for(int j = 0;j < len_x;j++) arr_true[true_len][j] = complete_perms[iter][j+start_pos];
        true_len++;
    }
    iter++;
    return return_val;
}

void perm_class::cout_line(){
    for(int i = 0;i < len_x;i++) cout << arr_perms[iter][i] << " ";
}


void perm_class::cout_line(int iter_2){
    for(int i = 0;i < len_x;i++) cout << arr_perms[iter_2][i] << " ";
}

int** perm_class::get_arr(){return arr_perms;}

int perm_class::get_len_x(){return len_x;}