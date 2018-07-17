#include "hold_all_perms.h"

using namespace std;

hold_all_perms::hold_all_perms(){
    
    get_all_perms_lens();
    get_all_perms();
    
}

hold_all_perms::~hold_all_perms(){
    for(int i = 0;i < 10;i++){
        for(int j = 0;j < len_all_perms[i];j++) delete hold_all_perms_arr[i][j];
        delete[] hold_all_perms_arr[i];
    }
    delete[] hold_all_perms_arr;
    delete[] len_all_perms;
}

void hold_all_perms::get_all_perms_lens(){
    auto factorial = [](int a){
        int ret_val = 1;
        for(int i = 1;i <= a;i++) ret_val *= i;
        return ret_val;
    };

    len_all_perms = new int[10];

    for(int i = 0;i < 10;i++) len_all_perms[i] = factorial(i+1);
} 

void hold_all_perms::get_all_perms(){
    hold_all_perms::hold_all_perms_arr = new int**[10];
    for(int i = 0; i < 10;i++){
        hold_all_perms_arr[i] = new int*[len_all_perms[i]];
        for(int j = 0;j < len_all_perms[i];j++) hold_all_perms_arr[i][j] = new int[i+1];
    }

    for(int i = 0;i < 10;i++) do_permutation(i+1);
}



void hold_all_perms::do_permutation(int iter){
    int sortarray[iter];
    for(int i = 0;i < iter;i++) sortarray[i] = i+1;

    int iter_2 = 0;
    sort(sortarray,sortarray+iter);
    do{
        for(int i = 0;i < iter;i++) hold_all_perms_arr[iter-1][iter_2][i] = sortarray[i];
        iter_2++;
    }while(next_permutation(sortarray,sortarray+iter));
}