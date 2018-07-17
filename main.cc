#include <iostream>
#include <string>
#include <thread>

#include "thread_handle.h"

using namespace std;

int main(int argc,const char* argv[]){
	
	
	bool verbose_flag = false;
    bool solve_type_bool = false;
    bool type_flag = false;
    bool accept_all = false;
	int verbosity = 0;
	for(int i = 0;i < argc;i++){
 		if(string(argv[i]) == "-v"){
			verbose_flag = true;
			continue;
		}
		if(verbose_flag){
			if(string(argv[i]) == "0") verbosity = 0;
			else if(string(argv[i]) == "1") verbosity = 1;
			else if(string(argv[i]) == "2") verbosity = 2;
			else{
				cout << "\nVerbosity has to be either 0, 1 or 2";
				cout << "\n=> Verbosity set to 0\n" << endl;
				verbosity = 0;
			}
            verbose_flag = false;
		}
        if(string(argv[i]) == "-s"){
            type_flag = true;
            continue;
        }
        if(type_flag){
            if(string(argv[i]) == "e") solve_type_bool = true;
            else if(string(argv[i]) == "a") solve_type_bool = false;
            type_flag = false;
        }
        if(string(argv[i]) == "-all"){
			accept_all = true;
			continue;
		}
	}
	
    int am_thr = 40;
    int offset = 0;
    
    cout.precision(4);

    //string dataname_dummy = "files/input/gamma_file_target_p_small_";
    string dataname_dummy = "files/input/Eu/gamma_file_target_Eu_p_small_";
    string dat_end = ".dat";

    string dataname_tmp;
    //string outputname = "files/output/gamma_file_target_out_p_small_21_";
    string outputname = "files/output/Eu/gamma_file_target_out_Eu_p_small_54321_";
    string outputname_tmp;
    

    thread_handle** run = new thread_handle*[am_thr];
     
    for(int i = 0+offset;i < am_thr+offset;i++){
        dataname_tmp = dataname_dummy + to_string(i) + dat_end;
        outputname_tmp = outputname + to_string(i) + dat_end;
        if(am_thr == 1 && false){
			dataname_tmp = "test_stuff/testitagain.dat";
			outputname_tmp = "test_perf_energies"+dat_end;
		}
        run[i-offset] = new thread_handle(dataname_tmp,outputname_tmp,i+1,verbosity,solve_type_bool);
    }   
    run[0]->call_title(am_thr);
    
    //init threads
    thread t[am_thr];
    for(int i = 0;i < am_thr;i++) t[i] = run[i]->threading();
    for(int i = 0;i < am_thr;i++) t[i].join();
    
    for(int i = 0;i < am_thr;i++) delete run[i];
    delete[] run;

	return 0;

}
