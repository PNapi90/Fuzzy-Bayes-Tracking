#include "perm_class.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

perm_class::perm_class(int len_x,int len_x_all){
    this->len_x = len_x;
    this->len_x_all = len_x_all;

    generate_arr_perms();
}


perm_class::~perm_class(){
    if(save_mode){
		delete[] recurrent_arr;
		delete[] prob_arr;
    }
    //for(int i = 0;i < all_perms_len;i++){
        //delete arr_perms[i];
        //delete arr_true[i];
    //}
    //delete[] arr_true;
    //delete[] arr_perms;
}

void perm_class::set_id(int interactions,int sub_cl_len,int position_x,int position_y){
	this->interactions = interactions;
	this->sub_cl_len = sub_cl_len;
	this->position_x = position_x;
	this->position_y = position_y;
	//cout << "ID: " << interactions << " " << sub_cl_len << " " << position_x << " " << position_y << endl;
}

void perm_class::generate_arr_perms(){
    auto factorial = [](int a){
        int ret_val = 1;
        for(int i = 1;i <= a;i++) ret_val *= i;
        return ret_val;
    };

    all_perms_len = factorial(len_x_all);
    
    //arr_perms = new int*[all_perms_len];
    //for(int i = 0;i < all_perms_len;i++) arr_perms[i] = new int[len_x];

    //arr_true = new int*[all_perms_len];
    //for(int i = 0;i < all_perms_len;i++) arr_true[i] = new int[len_x];
    //cout << "created perm_class obj with len_x " << len_x << endl;

}


void perm_class::feed_complete(int** complete_perms,int am_perm,int start_pos){
    
    this->complete_perms = complete_perms;
    //cout << "= = = = = = = = =\n";
    //cout << am_perm << " " << len_x << endl;
    //cout << "- - - - - - - -\n";
    
    this->start_pos = start_pos;
    this->am_perm = am_perm;
    //for(int i = 0;i < this->am_perm;i++){
    //    for(int j = 0;j < len_x;j++) arr_perms[i][j] = this->complete_perms[i][j+this->start_pos];
    //}
    iter = 0;
    true_len = 0;
    //cout << "am_perm" << am_perm <<endl;
    
    if(save_mode){
		recurrent_arr = new int[am_perm];
		prob_arr = new double[am_perm];
		get_recurrent_arr();
		save_recurrent_arr();
	}
	//else load_recurrent_arr();
}


void perm_class::save_recurrent_arr(){
	
	string data_name = "files/get_recurrent/";
	data_name += to_string(interactions) + "int/";
	data_name += to_string(sub_cl_len) + "_" + to_string(position_x);
	data_name += "_" + to_string(position_y) + ".dat";
	
	ofstream data_stream(data_name);
	for(int i = 0;i < am_perm;i++){
		data_stream << recurrent_arr[i] << "\n";
	}
	
	data_stream.close();

}



inline bool perm_class::check_all(int i,int j){
    int num_trues = 0;
    for(int k = 0;k < len_x;k++){
        if(complete_perms[i][k+start_pos] == complete_perms[j][k+start_pos]) num_trues += 1;
    }
    if(num_trues == len_x) return true;
    return false;
}

void perm_class::get_recurrent_arr(){
	bool recurrent = false;
    for(int i = 0;i < am_perm;i++){
		recurrent = false;
        for(int j = 0;j < i;j++){
			if(len_x >= len_x_all - 1) break;
            if(check_all(i,j)){
				recurrent_arr[i] = j;
				recurrent = true;
				break;
			}
		}
        if(!recurrent) recurrent_arr[i] = i;
        if(am_perm == 24 && false) cout << i << " -> " << recurrent_arr[i] <<"\n=========== "<<endl;   
    }
}

bool perm_class::check_recurrent(int j){
    if(recurrent_arr[j] != j) return true;
    return false;
}

double perm_class::get_probability(int j){return prob_arr[recurrent_arr[j]];}

void perm_class::set_probability(double val,int j){
    prob_arr[j] = val;
}


bool perm_class::test_for_double(){
    int true_count = 0;
    bool return_val = false;
    for(int i = 0;i < iter;i++){
        for(int j = 0;j < len_x;j++){
            if(complete_perms[i][j+start_pos] == complete_perms[iter][j+start_pos]) true_count += 1;
        }
        if(true_count == len_x){
            return_val = true;
            break;
        }
        true_count = 0;
    }
    if(!return_val){
        //for(int j = 0;j < len_x;j++) arr_true[true_len][j] = complete_perms[iter][j+start_pos];
        true_len += 1;
    }
    iter += 1;
    return return_val;
}

void perm_class::cout_line(){
    for(int i = 0;i < len_x;i++) cout << arr_perms[iter][i] << " ";
}


void perm_class::cout_line(int iter_2){
    for(int i = 0;i < len_x;i++) cout << complete_perms[iter_2][i+start_pos] << " ";
}

//int** perm_class::get_arr(){return arr_perms;}

int perm_class::get_len_x(){return len_x;}

int perm_class::return_line(int k,int j){return complete_perms[k][j+start_pos];}
