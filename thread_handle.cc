#include "thread_handle.h"

using namespace std;


thread_handle::thread_handle(string input,string output,int thread_num,int verbosity,bool solve_type_bool){
    this->input = input;
    this->output = output;
    this->thread_num = thread_num;
	this->verbosity = verbosity;
    this->solve_type_bool = solve_type_bool;

	bool doit = false;
	if(this->thread_num == 1) doit = true; 
    data = new data_handler(input,output,amount,doit);
}

thread_handle::~thread_handle(){
    delete data;
}

void thread_handle::do_calculations(){

    //unclustered data, but seperated by TStamps
    double** TS_clustered = data->get_data();
    int* cluster_sizes = data->get_cluster_sizes();

    int cluster_pos = 0;
    int cluster_iter = 0;
    int amount_clusters = data->get_am_clusters();

    double** cluster_dummy = new double*[30];
    for(int i = 0;i < 30;i++) cluster_dummy[i] = new double[4];


    //cluster data for each TStamp
    cluster_class cluster(thread_num,verbosity,solve_type_bool);
    int internal_iter = 0;
    double etmp_sum = 0;
    while(cluster_iter < amount_clusters){
		etmp_sum = 0;
        if(abort_statement(cluster_sizes[cluster_iter])){
            for(int i = 0;i < cluster_sizes[cluster_iter];i++){
                for(int j = 0;j < 4;j++) cluster_dummy[i][j] = TS_clustered[i+cluster_pos][j];
                etmp_sum += cluster_dummy[i][0];
            }
			if(abs(etmp_sum - 1173.2 - 1332.5) < 6 || true){ 
				//do clustering
				cluster.clustering(cluster_dummy,cluster_sizes[cluster_iter]);
        
				//save new clusters in output
				for(int i = 0;i < cluster.amount();i++){
					data->save(cluster.get_cluster(i),cluster.cluster_len(i));
				}
			}
        }
        
        //increment to next cluster
        cluster_pos += cluster_sizes[cluster_iter];
        cluster_iter += 1;

        cluster.reset();
        internal_iter += 1;
    }

    for(int i = 0;i < 30;i++) delete[] cluster_dummy[i];
    delete[] cluster_dummy;
	//cout << "am_torn in thr " << thread_num << " -> " << cluster.get_am_torn() << endl;
    //cout << "\n!done!\n" << endl;
}


thread thread_handle::threading(){return thread([=] {do_calculations();} );}


void thread_handle::call_title(int am_threads){
	string singular = " thread";
	string plural = " threads";
	
	string pl_string;
	if(am_threads == 1) pl_string = singular;
	else pl_string = plural;
	
    cout << "\nWelcome to FUZZY BAYES-TRACKING (FBT)" << endl;
    cout << "\nstarting threaded mode with "<< am_threads << pl_string <<"\n\n"; 
}


bool thread_handle::abort_statement(int i){return (i <= 5 && i > 0);}
