/*#include <iostream>

#include "rec_class.h"

using namespace std;


int main(){

    rec_class r;

    r.rec(6,1,0);
    r.cout_it();

    cout << endl;

    return 0;
}*/


#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

/*
void bars(int am_bars,int balls){
    int len = am_bars+balls;
    int count = 0;
    int x[len];
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
                cout << sum << " | ";
                sum = 0;
                count++;
            }
            if(i == len-1 && count == am_bars) cout << sum;
        }
        cout << endl;
        count  = 0;
    }while(std::next_permutation(x,x+len));
    
}*/

#include "rec_class.h"


int main(){
    int stuff = 10;
    rec_class r(stuff);
    r.bars(stuff);
    r.cout_it();

    return 0;
}

