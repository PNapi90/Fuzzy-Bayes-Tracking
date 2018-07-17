#include <iostream>
#include <string>

#include "data_class.h"
#include "cluster.h"

using namespace std;

int main(){

    cout << "\nWelcome to FUZZY BAYES-TRACKING (FBT)" << endl;
    cout << "\nstarting FUZZY BAYES-CLUSTERING ... \n\n"; 

    cout.precision(4);

    int amount = 10000;

    string input = "test2.dat";
    string output = "clustered_test2.dat";
    data_handler data(input,output,amount);
    
    //unclustered data, but seperated by TStamps
    double** TS_clustered = data.get_data();
    int* cluster_sizes = data.get_cluster_sizes();

    int cluster_pos = 0;
    int cluster_iter = 0;
    int amount_clusters = data.get_am_clusters();

    double** cluster_dummy = new double*[30];
    for(int i = 0;i < 30;i++) cluster_dummy[i] = new double[4];


    //cluster data for each TStamp
    cluster_class cluster;
	int internal_iter = 0;
    while(cluster_iter < amount_clusters-1){
		if(cluster_sizes[cluster_iter] < 6){
            for(int i = 0;i < cluster_sizes[cluster_iter];i++){
                for(int j = 0;j < 4;j++) cluster_dummy[i][j] = TS_clustered[i+cluster_pos][j];
            }

            //do clustering
		    cluster.clustering(cluster_dummy,cluster_sizes[cluster_iter]);
        
            //save new clusters in output
            for(int i = 0;i < cluster.amount();i++){
                data.save(cluster.get_cluster(i),cluster.cluster_len(i));
            }
        }
        
        //increment to next cluster
        cluster_pos += cluster_sizes[cluster_iter];
        cluster_iter++;

        cluster.reset();
        internal_iter++;
    }

    for(int i = 0;i < 30;i++) delete cluster_dummy[i];
    delete[] cluster_dummy;

    cout << "\n!done!\n" << endl;

    return 0;
}
