#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

void save_in_cluster(double***,double**,int,int,int*);

int main(){

    ifstream data("gamma_file_target_no_doubles_Eu.dat");

    double* line = new double[6];

    double** tmp_arr = new double*[30];
    for(int i = 0;i < 30;++i) tmp_arr[i] = new double[6];

    double*** write_arr = new double**[40];
    for(int i = 0;i < 40;++i){
        write_arr[i] = new double*[100000];
        for(int j = 0;j < 100000;++j) write_arr[i][j] = new double[6];
    }
    
    bool abort_t = false;
    double edep = 0;
    int iter = 0;

    int am_gam[40] = {0};
    int gam_iter = 0;
    int* counter_arr = new int[40];
    for(int i = 0;i < 40;++i) counter_arr[i] = 0;
    
    double* edep_tot = new double[40*1000];
    for(int i = 0;i < 40000;++i) edep_tot[i] = 0;
    int tot_iter = 0;
    
    while(true){

        if(gam_iter == 40) break;

        for(int i = 0;i < 6;++i) data >> tmp_arr[iter][i];
        if(tmp_arr[iter][0] == -99999){
            iter += 1;
            if(abs(edep - 1332.5) < 6 || true){
                edep_tot[tot_iter] = edep;
                tot_iter += 1;
                save_in_cluster(write_arr,tmp_arr,iter,gam_iter,counter_arr);
                am_gam[gam_iter] += 1;
                if(am_gam[gam_iter] % 10 == 0) cout << "thr " << gam_iter << " -> " << (double) am_gam[gam_iter]/10. << " %" << endl;
                if(am_gam[gam_iter] == 1000){
                    gam_iter += 1;
                }
            }
            edep = 0;
            iter = 0;
            for(int k = 0;k < 30;++k) for(int j = 0;j < 6;++j) tmp_arr[k][j] = 0;
            continue;
        }
        edep += tmp_arr[iter][0];
        iter += 1;
    }
    string fol = "files/input/Eu/";
    ofstream** data_out = new ofstream*[40];
    for(int i = 0;i < 40;++i) data_out[i] = new ofstream(fol+"gamma_file_target_Eu_p_small_"+to_string(i)+".dat");

    for(int i = 0;i < 40;++i){
        for(int j = 0;j < counter_arr[i];++j){
            for(int k = 0;k < 5;++k) *data_out[i] << write_arr[i][j][k] << "\t";
            *data_out[i] << write_arr[i][j][5] << "\n";
        }
        data_out[i]->close();
    }
    ofstream tot_data("edep_tot.dat");
    for(int i = 0;i < tot_iter;++i) tot_data << edep_tot[i] << "\n";
    tot_data.close();

    return 0;
}


void save_in_cluster(double*** writer,double** tmp,int len,int pos,int* counter_arr){
    for(int i = 0;i < len;++i){
        for(int j = 0;j < 6;++j) writer[pos][counter_arr[pos]][j] = tmp[i][j];
        counter_arr[pos] += 1;
    }
}
