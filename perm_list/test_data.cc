#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;


int k_out_n(int n,int k){
    auto factorial = [](int a){
        int ret_val = 1;
        if(a > 0) for(int i = 1;i <= a;i++) ret_val *= i;
        return ret_val;
    };
    return factorial(n)/(factorial(k)*factorial(n-k));
}
/*
void load_poss_combs(){
    get_max_perms();
    
    string name = "Perms_";
    string ending = ".dat";
    string folder;
    
    ifstream*** data_arr = new ifstream**[max_interactions];
    for(int i = 0;i < max_interactions;i++){
        //len[i] = get_perms(i+1);
        data_arr[i] = new ifstream*[i+1];
        folder = "ints/"+to_string(i+1)+"_int/";
        for(int j = 0;j < i+1;j++) data_arr[i][j] = new ifstream((folder+name+to_string(j+1)+ending).c_str());
    }
    
    for(int i = 0;i < max_interactions;i++){
        for(int j = 0;j <= i;j++){
            for(int k = 0;k <= j){
                for(int l = 0;l < max_perms[i][j];j++){
                    *data_arr[i][j] >> permutation_sub_arr[i][j][k][l];
                }
            }
        }
    }
}

void get_max_perms(){
    for(int i = 0;i < max_interactions;i++){
        for(int j = 0;j <= i;j++){
            max_perms[i][j] = k_out_n(i+1,j+1);
        }
    }
}
*/
int main(){
    
    int max_interactions = 4;

    int** max_perms = new int*[max_interactions];
    for(int i = 0;i < max_interactions;i++) max_perms[i] = new int[i+1];

    for(int i = 0;i < max_interactions;i++){
        cout << "------\n" << i+1 << endl;
        for(int j = 0;j <= i;j++){
            max_perms[i][j] = k_out_n(i+1,j+1);
            cout << i+1 <<" " << j+1 <<"    "<<max_perms[i][j] << "\n";
        }
        cout <<"o o o o o"<<endl;
    }
    cout << endl;
    string name = "Perms_";
    string ending = ".dat";
    string folder;
    
    ifstream*** data_arr = new ifstream**[max_interactions];
    for(int i = 0;i < max_interactions;i++){
        //len[i] = get_perms(i+1);
        data_arr[i] = new ifstream*[i+1];
        folder = "ints/"+to_string(i+1)+"_int/";
        for(int j = 0;j < i+1;j++) data_arr[i][j] = new ifstream((folder+name+to_string(j+1)+ending).c_str());
    }
    int dummy_val = 0;
    for(int i = 0;i < max_interactions;i++){
        for(int j = 0;j <= i;j++){
            for(int l = 0;l < max_perms[i][j];l++){
                for(int k = 0;k <= j;k++){
                    *data_arr[i][j] >> dummy_val;
                    cout << dummy_val << "\t";
                }
                cout << endl;
            }
            cout << " - - - - - - - - - - - - " << endl;
        }
        cout << " # # # # # # # # # # # # # # # " << endl;
    }




    return 0;
}