#include "cross_sec.h"

using namespace std;

cross_sec::cross_sec(){

    sigma = new double[2];
    load_xsecs();

}

cross_sec::~cross_sec(){

    for(int i = 0;i < 64;i++) delete[] comp_xsec[i];
    for(int i = 0;i < 91;i++) delete[] photo_xsec[i];
    for(int i = 0;i < 91;i++) delete[] photo_xsec_normed[i];
    for(int i = 0;i < 64;i++) delete[] comp_xsec_normed[i];
    delete[] comp_xsec;
    delete[] comp_xsec_normed;
    delete[] photo_xsec;
    delete[] photo_xsec_normed;
    delete[] sigma;
}

void cross_sec::load_xsecs(){
    //init arrays
    comp_xsec = new double*[64];
    photo_xsec = new double*[91];
    comp_xsec_normed = new double*[64];
    photo_xsec_normed = new double*[91];
    
    //load xsecs from NIST files
    // - Compton -  
    ifstream comptonfile("inputdata/compton_sigma.dat");
    for(int i = 0;i < 64;i++){
        comp_xsec[i] = new double[2];
        comp_xsec_normed[i] = new double[2];
        if(!comptonfile){
            cout<<"could not find Compton cross secs data"<<endl;
            break;
        }
        for(int j = 0;j < 2;j++){
            comptonfile >> comp_xsec[i][j];
        }
    }
    // - photo -    
    ifstream photofile("inputdata/photo_sigma.dat");
    for(int i = 0;i < 91; ++i){
        photo_xsec[i] = new double[2];
        photo_xsec_normed[i] = new double[2];
        if(!photofile){
            cout<<"could not find photo cross secs data"<<endl;
            break;
        }
        for(int j = 0;j < 2;j++){
            photofile >> photo_xsec[i][j];
        }
    }
    normalize();
}

void cross_sec::normalize(){
    double sum = 0;
    for(int i = 0;i < 91;i++) sum += photo_xsec[i][1];
    for(int i = 0;i < 91;i++){
        photo_xsec_normed[i][0] = photo_xsec[i][0];
        photo_xsec_normed[i][1] /= sum;
    }

    sum = 0;
    for(int i = 0;i < 64;i++) sum += comp_xsec[i][1];
    for(int i = 0;i < 64;i++){
        comp_xsec_normed[i][0] = comp_xsec[i][0];
        comp_xsec_normed[i][1] /= sum;
    }
    
}

void cross_sec::get_sigma_from_list(double E){
    for(int j = 0;j < 64;j++){
        if(E >= comp_xsec[j][0] && E < comp_xsec[j+1][0]){
            sigma[1] = comp_xsec[j][1]*(comp_xsec[j+1][0]-E)/(comp_xsec[j+1][0]-comp_xsec[j][0])
                    +comp_xsec[j+1][1]*(E - comp_xsec[j][0])/(comp_xsec[j+1][0]-comp_xsec[j][0]);
            break;
        }
    }   
    
    for(int j = 0;j < 91;j++){
        if(E >= photo_xsec[j][0] && E < photo_xsec[j+1][0]){
            sigma[0] = photo_xsec[j][1]*(photo_xsec[j+1][0]-E)/(photo_xsec[j+1][0]-photo_xsec[j][0])
                    +photo_xsec[j+1][1]*(E - photo_xsec[j][0])/(photo_xsec[j+1][0]-photo_xsec[j][0]);
            break;
        }
    }
}

double cross_sec::get_sigma_from_list_normed(double E){  
    double sig_ret = 0;
    for(int j = 0;j < 91;j++){
        if(E >= photo_xsec_normed[j][0] && E < photo_xsec_normed[j+1][0]){
            sig_ret = photo_xsec_normed[j][1]*(photo_xsec_normed[j+1][0]-E)/(photo_xsec_normed[j+1][0]-photo_xsec_normed[j][0])
                    +photo_xsec_normed[j+1][1]*(E - photo_xsec_normed[j][0])/(photo_xsec_normed[j+1][0]-photo_xsec_normed[j][0]);
            break;
        }
    }
    return sig_ret;
}


double* cross_sec::get_sigma(double E){
    get_sigma_from_list(E*1e-3);
    return sigma;
}

double cross_sec::get_sigma_normed(double E){
    return get_sigma_from_list_normed(E*1e-3);
}
