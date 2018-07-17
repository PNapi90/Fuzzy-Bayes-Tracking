#include <iostream>
#include <algorithm>

using namespace std;

void comb(int N,int K,int** arr){
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's
    int iter = 0;
    int iter_2 = 0;
    //cout << "=======================\n" << N << " " << K << "\n------------------------" << endl;
    do {
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i]){
                arr[iter][iter_2] = i+1;
                iter_2++;
            }
        }
        iter++;
        iter_2 = 0;
        //cout << "\n";
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    //cout << "=======================\n";
}

int main(){
    int*** arr = new int**[5];
    for(int i = 0;i < 5;i++){
        arr[i] = new int*[200];
        for(int j = 0;j < 200;j++) arr[i][j] = new int[5];
    }
    int o = 0;
    for(int i = 1;i < 5;i++) comb(5,i,arr[i-1]);
    for(int i = 0;i < 10;i++){
        //nt i = 0;
        for(int j = 0;j < 3;j++) cout << arr[2][i][j] <<" ";
        cout << arr[0][4-o][0] << " " << arr[0][4-(2-o) ][0] << " ";
        cout << endl;
        o++;
        if(o == 5) o = 0;
    }
    return 0;
}