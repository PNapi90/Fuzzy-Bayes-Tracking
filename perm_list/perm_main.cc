#include <iostream>
#include <algorithm>

#include "perm_class.h"
#include "merge_class.h"

using namespace std;



int main(){

    int len = 4;

    int sortarray[len];
    for(int i = 0;i < len;i++) sortarray[i] = i+1;

    sort(sortarray,sortarray + len);
    
    int am_perms = 1;
    for(int i = 1;i <= len;i++) am_perms *= i;

    int** all_perms = new int*[am_perms];
    for(int i = 0;i < am_perms;i++) all_perms[i] = new int[len];

    int iter = 0;
    do{
        for(int i = 0;i < len;i++) all_perms[iter][i] = sortarray[i];
        iter++;
    }while(next_permutation(sortarray,sortarray + len));

    int** all_poss = new int*[5];
    for(int i = 0;i < 5;i++) all_poss[i] = new int[4];
    int am_clust[5] = {1,2,2,3,4};
    for(int i = 0;i < 5;i++) for(int j = 0;j < 4;j++) all_poss[i][j] = 0;
    all_poss[0][0] = 4;
    all_poss[1][0] = 3;
    all_poss[1][1] = 1;
    all_poss[2][0] = 2;
    all_poss[2][1] = 2;
    all_poss[3][0] = 2;
    all_poss[3][1] = 1;
    all_poss[3][2] = 1;
    for(int i = 0;i < 4;i++) all_poss[4][i] = 1;

    perm_class*** A = new perm_class**[5];
    
    for(int i = 0;i < 5;i++){
        A[i] = new perm_class*[4];
        for(int j = 0;j < 4;j++) if(all_poss[i][j] > 0) A[i][j] = new perm_class(all_poss[i][j]);
    }
    
    int start_pos = 0;
    for(int i = 0;i < 5;i++){
        for(int j = 0;j < 4;j++){
            if(all_poss[i][j] > 0){
                
                if(j == 0) A[i][j]->feed_complete(all_perms,am_perms,len,0);
                else{
                    start_pos += all_poss[i][j-1];
                    A[i][j]->feed_complete(all_perms,am_perms,len,start_pos);
                }
                
            }
        }
        start_pos = 0;
    }

    merge_class** B = new merge_class*[5];
    for(int i = 0;i < 5;i++){
        B[i] = new merge_class(all_poss[i],am_clust[i],am_perms,A[i]);
        B[i]->compare_true();
        cout << "================================\n";
    } 

    /*
    int true_count = 0;
    int poss = 0;
    for(int i = 0;i < 5;i++){
        for(int j = 0;j < 4;j++){
            if(all_poss[i][j] > 0){
                poss++;
                if(A[i][j]->test_for_double()) true_count++;
            }
        }
        if(true_count < poss){
            for(int j = 0;j < 4;j++){
                if(all_poss[i][j] > 0){
                    A[i][j]->cout_line();
                    cout << "| ";
                }
            }
            cout << endl;
        }
        poss = 0;
        true_count = 0;
    }
    */
    for(int i = 0;i < 5;i++){
        for(int j = 0;j < 4;j++) if(all_poss[i][j] > 0) delete A[i][j];
        delete[] A[i];
    }
    delete[] A;
    
    for(int i = 0;i < am_perms;i++) delete all_perms[i];
    delete[] all_perms;

    return 0;
}