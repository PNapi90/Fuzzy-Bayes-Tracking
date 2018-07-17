#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

void comb(int N,int K,ofstream* data,ofstream* miss){
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's
    int iter = 0;
    int iter_intern = 0;
    do {
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i]) *data << i+1 <<"\t";
            else *miss << i+1<<"\t";
        }
        *data << "\n";
        *miss << "\n";
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

int main(){
    int max_interactions = 10;
    string name = "Perms_";
    string name_miss = "Perms_miss_";
    string ending = ".dat";
    string folder;
    ofstream*** data_arr = new ofstream**[max_interactions];
    for(int i = 0;i < max_interactions;i++){
        //len[i] = get_perms(i+1);
        data_arr[i] = new ofstream*[i+1];
        folder = "ints/"+to_string(i+1)+"_int/";
        for(int j = 0;j < i+1;j++) data_arr[i][j] = new ofstream((folder+name+to_string(j+1)+ending).c_str());
    }

    ofstream*** miss_arr = new ofstream**[max_interactions];
    for(int i = 0;i < max_interactions;i++){
        //len[i] = get_perms(i+1);
        miss_arr[i] = new ofstream*[i+1];
        folder = "ints/"+to_string(i+1)+"_int/";
        for(int j = 0;j < i+1;j++) miss_arr[i][j] = new ofstream((folder+name_miss+to_string(j+1)+ending).c_str());
    }
    
    
    for(int i = 0;i < max_interactions;i++){
        for(int j = 0;j <= i;j++){
            comb(i+1,j+1,data_arr[i][j],miss_arr[i][j]);
        }
    }
    for(int i = 0;i < max_interactions;i++){
        for(int j = 0;j <= i;j++){
            data_arr[i][j]->close();
            miss_arr[i][j]->close();
        }
    }
    return 0;
}