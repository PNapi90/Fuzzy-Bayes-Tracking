#include "data_class.h"
#include <string>
using namespace std;


//ctor -> opens in- and output files and load input data
data_handler::data_handler(string input_s,string output_s,int len,bool doit){
    scnd_ctor = false;
    this->len = len;
    iter = 0;
    iter_c = 0;
    input.open(input_s.c_str());
    output.open(output_s.c_str());
	thrd_ctor = false;
	fth_ctor = false;
    load_input(false);
}

data_handler::data_handler(int thr_num,bool solve_type_bool){
	if(thr_num == 1) cout << "2nd ctor of data_handler called" << endl;
	if(solve_type_bool) integral_data.open("temp_int_data_en.dat");
	else integral_data.open("temp_int_data_an.dat");
    
    cluster_data.open("temp_cluster_data.dat");
    l2_norm_file.open("files/l2_norms/l2norms_all_"+to_string(thr_num)+".dat");
	scnd_ctor = true;
    first_comb = true;
    internal_ints = new double*[10000];
    for (int i = 0; i < 10000;i++){
        internal_ints[i] = new double[2];
        internal_ints[i][0] = (double) i + 1.;
        internal_ints[i][1] = 0.;
    }
    iter = 0;
    thrd_ctor = false;
    fth_ctor = false;
}


//ctor default
data_handler::data_handler(int thr_num){
	//3rd ctor
	sum_file.open("files/sum_files/sum_file_"+to_string(thr_num)+".dat");
	thrd_ctor = true;
	scnd_ctor = false;
	fth_ctor = false;
}

data_handler::data_handler(int thr_num,char c){
	//4th ctor
	lambda_file.open("files/lambda_files/lambda_file_"+to_string(thr_num)+".dat");
	thrd_ctor = false;
	scnd_ctor = false;
	fth_ctor = true;
}

//dtor
data_handler::~data_handler(){
	if(scnd_ctor){
        integral_data.close();
        cluster_data.close();
        for (int i = 0; i < 10000; ++i){
            delete[] internal_ints[i];
        }
        delete[] internal_ints;
        
    }
    if(thrd_ctor) sum_file.close(); 
    if(!thrd_ctor && !scnd_ctor && !fth_ctor){
		output.close();
		for(int i = 0;i < len;i++){
			delete[] input_arr[i];
		}
		delete[] input_arr;
		delete[] cluster_size_arr;
	}
	if(fth_ctor) lambda_file.close();
}

void data_handler::save_lambda_vals(double* data_stream,int len){
	for(int i = 0;i < len;i++) lambda_file << data_stream[i] << "\t";
	lambda_file << "\n";
	lambda_file.flush();
}

void data_handler::save_sum_val(double E,double sum,double l2){
	sum_file << E << "\t" << sum << "\t"<< l2 << "\n";
	sum_file.flush();
}

void data_handler::save_l2_norm(double start,int interactions,double l2,double peak,double sum){
	l2_norm_file << start << "\t" << interactions << "\t" << l2 << "\t" << peak << "\t" << sum  << "\n";
	l2_norm_file.flush();
}

void data_handler::save_perm_total(int cluster_iter,int j,double p){
    cluster_data << cluster_iter << "\t" << j << "\t" << p << "\n";
}

//loads input file to input_arr
void data_handler::load_input(bool doit){
    cluster_size_arr = new int[len];
    input_arr = new double*[len];
    for(int i = 0;i < len;i++){
        input_arr[i] = new double[4];
        cluster_size_arr[i] = 0;
    }
    string line;
    double dummy_arr[4];
    double drop_it[2];

    while(input.good()){
        getline(input,line,'\n');
        if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&dummy_arr[0],&dummy_arr[1],&dummy_arr[2],&dummy_arr[3],&drop_it[0],&drop_it[1]);
        if(dummy_arr[0] != c_delim){
            for(int i = 0;i < 4;i++){
                input_arr[iter][i] = dummy_arr[i];
                if(i == 1 || i == 3) input_arr[iter][i] *= -1;
                if(doit) cout << input_arr[iter][i] << " ";
            }
            iter += 1;
            cluster_size_arr[iter_c]++;
            if(doit) cout << cluster_size_arr[iter_c] << endl;
        }
        else{
			iter_c += 1;
			if(doit) cout << "-------------------------------\n";
		}
    }
}

void data_handler::next_perm(){iter = 0;}

void data_handler::next_comb(){
    for(int i = 0;i < 10000;i++) internal_ints[i][1] = 0.;
    first_comb = false;
    integral_data << -3232 << "\t" << -3232 << "\n";
}

//save integral data
void data_handler::save_integral(double e,double integral){
	integral_data << e << "\t" << integral << "\n";
}

void data_handler::next_gamma(){integral_data << gamma_delim << "\t" << gamma_delim << "\n";}

void data_handler::next(){integral_data << -99999 << "\t" << -99999 <<"\t "<< -99999 <<"\t"<< -99999 <<"\t"<< -99999 << "\n";}

void data_handler::next_iter(){integral_data << -42 << "\t" << -42 << "\n";}

//saves cluster x in output file (with delimiter)
void data_handler::save(double** x,int cluster_size){
    for(int i = 0;i < cluster_size;i++){
        for(int j = 0;j < 4;j++) output << x[i][j] << "\t";
        output << cluster_size << "\n";
    }
    for(int i = 0;i < 4;i++) output << dummy_val << "\t";
    output << dummy_val << "\n";
    output.flush();
}

//passes input_arr to cluster alg.
double** data_handler::get_data(){return input_arr;}

//passes (TStamp) cluster_size_arr to cluster alg.
int* data_handler::get_cluster_sizes(){return cluster_size_arr;}

int data_handler::get_am_clusters(){return iter_c;}
