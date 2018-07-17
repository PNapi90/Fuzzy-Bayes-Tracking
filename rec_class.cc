#include "rec_class.h"


using namespace std;


rec_class::rec_class(int a){
    
    am = a;
    arr = new int*[100000];
    for(int i = 0;i < 100000;i++) arr[i] = new int[200];
    //out = new ofstream*[10];
    //for(int i = 1;i <= 10;i++) out[i-1] = new ofstream("out_of/test_"+to_string(i)+".dat");

}


rec_class::~rec_class(){

    for(int i = 0;i < 200;i++) delete arr[i];
    delete[] arr;
    //for(int i = 0;i < 10;i++) out[i]->close();
}


void rec_class::bars(int balls){
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
                //cout << int_iter << endl;
                arr[int_iter][a] = sum;
                cout << sum << " | ";
                sum = 0;
                a++;
                count++;
            }
            if(i == len-1 && count == am_bars){
                cout << sum;
                arr[int_iter][a] = sum;
            }
        }
        int_iter++;
        a = 0;
        cout <<" || "<< int_iter << endl;
        count  = 0;
    }while(std::next_permutation(x,x+len));
    leni = int_iter;
}

void rec_class::cout_it(){
    //int plot_arr[20] = {0,0,0,0,0};
    int old = 0;
    bool plotter = false;
    int true_ones = 0;
    for(int i = 0;i < leni;i++){
        plotter = false;
        old = arr[i][0];
        for(int j = 1;j < am;j++){
            if(old >= arr[i][j] && old > 0){
                plotter = true;
                //cout <<old <<" "<<arr[i][j] << " / ";
            }
            else if(old < arr[i][j]){
                plotter = false;
                break;
            }
            old = arr[i][j];
            //else plotter = false;
        }
        //cout << endl;
        for(int j = 0;j < am;j++){
            if(plotter){
                cout << arr[i][j] << " ";
            }
        }
        if(plotter){
            true_ones++;
            cout << endl;
        }
    }
    cout << "\ntrue ones " << true_ones << endl;
}


