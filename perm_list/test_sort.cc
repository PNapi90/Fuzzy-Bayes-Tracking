#include <iostream>
#include <algorithm>
#include <fstream>


using namespace std;


bool test_for_double(int** doubles,int* all,int len){
    for(int i = 0;i < len;i++){
        //cout << doubles[i][0]<< "  " << doubles[i][1] << "\t" << all[0] << " " << all[1] << endl;
        if(doubles[i][0] == all[0] && doubles[i][1] == all[1]) return false;
    }
    return true;
}


int main(){

    //ofstream data("testi.txt");
    int len = 4;
    int sortarray[len];
    int** all_perms = new int*[24];
    for(int i = 0;i < 24;i++) all_perms[i] = new int[4];
    for(int i = 0;i < len;i++) sortarray[i] = i+1;
    sort(sortarray,sortarray+len);
    int p = 0;
    do{
        for(int i = 0;i < len;i++) all_perms[p][i] = sortarray[i];
        p++;
    }while(next_permutation(sortarray,sortarray+len));

    //data.close();
    //cout << "hello" << endl;
    int perms_all = 1;
    int** twos = new int*[30];
    int** ones = new int*[30];
    for(int i = 0;i < 30;i++){
        twos[i] = new int[2];
        ones[i] = new int[2];
    }
    
    for(int i = 1;i <= len;i++) perms_all *= i;
    
    int iter = 0;
    //cout << perms_all << endl;
    for(int i = 0;i < perms_all;i++){
        if(test_for_double(twos,all_perms[i],iter)){
            for(int j = 0;j < 2;j++){
                twos[iter][j] = all_perms[i][j];
                ones[iter][j] = all_perms[i][j+2];
                //cout << twos[iter][j] << " ";
            }
            //cout << endl;
            iter++;
        }
    }
    for(int i = 0;i < iter;i++){
        for(int j = 0;j < 2;j++) cout << twos[i][j] << " ";
        for(int j = 0;j < 2;j++) cout << ones[i][j] << " ";
        cout << endl;
    }


    return 0;
}