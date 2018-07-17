#include "G_Matrix_slim.h"

using namespace std;

G_Matrix_slim::G_Matrix_slim(bool solve_type_bool){

    this->solve_type_bool = solve_type_bool;
	get_source_pos();
    x_scalars = new double[10];
    x_norms = new double[10];
    x_angles = new double[10];

	delta_edep_arr = new double[10];
	delta_theta_arr = new double[10];

    d_x_mat = new double[30];
    d_e_mat = new double[10];
    Edep = new double[10];
    error_ei = new double[10];
    ei = new double[10];
    x = new double*[10+3];
    for(int i = 0;i < 10+3;i++) x[i] = new double[3];

    sigma_vec = new double*[10];
    for(int i = 0;i < 10;i++) sigma_vec[i] = new double[2];

}

G_Matrix_slim::~G_Matrix_slim(){
    
    for (int i = 0; i < 10;i++){
        delete[] x[i];
        delete[] sigma_vec[i];
    }
    for(int i = 10;i < 13;i++) delete[] x[i];
    delete[] d_e_mat;
    delete[] d_x_mat;
    delete[] ei;
    delete[] Edep;
    delete[] sigma_vec;
    delete[] x;
    delete[] x_scalars;
    delete[] x_norms;
    delete[] x_angles;
    
    delete[] delta_edep_arr;
    delete[] delta_theta_arr;
    delete[] error_ei;


}

void G_Matrix_slim::get_source_pos(){
	ifstream source_data("inputdata/source_pos.dat");
	for(int i = 0;i < 3;++i) source_data >> source[i];
}


void G_Matrix_slim::set_order(bool leading_order){this->leading_order = leading_order;}


void G_Matrix_slim::set_tolerances(){	
	//if(leading_order) tolerance = tolerance_raw*tolerance_arr[interactions-1];
	//else
	tolerance = tolerance_raw;//*tolerance_arr[interactions-1];
}


double* G_Matrix_slim::get_sigmas(double** data_x,double* ed,int interactions,double e){
    ei[0] = e;
    this->interactions = interactions;
    dim_m = interactions*3 + 3 + 1;
    for(int i = 0;i < 3;i++) x[0][i] = source[i];
    for(int i = 0;i < interactions;i++){
        Edep[i] = ed[i];
        ei[i+1] = ei[i] - Edep[i];
        for(int j = 0;j < 3;j++) x[i+1][j] = data_x[i][j];
    }
   	
	set_tolerances();
	calc_norms();
	calc_scalars();
	set_angles();
	
	
    set_derivatives();
    calc_d_e_mat();
    set_sigmas();

    return delta_theta_arr;
}

double* G_Matrix_slim::get_sigmas(double** data_x,double* ed,int interactions,double e,bool a){
    ei[0] = e;
    this->interactions = interactions;
    dim_m = interactions*3 + 3 + 1;
    for(int i = 0;i < 3;i++) x[0][i] = source[i];
    for(int i = 0;i < interactions;i++){
        for(int j = 0;j < 3;j++) x[i+1][j] = data_x[i][j];
    }

    set_derivatives();
    set_angles();
    set_Edeps();
    set_sigmas(a);

    return delta_edep_arr;
}

double G_Matrix_slim::x_derivatives(int k,int i){
	double val = 0;
	if(i == 0){
		for(int j = 0;j < 3;j++){
			val += (x[k][j] - x[k+1][j])/(x_norms[k-1]*x_norms[k]) - x_angles[k-1]*(x[k-1][j] - x[k][j])/pow(x_norms[k-1],2);
		}
	}
	else if(i == 1){
		for(int j = 0;j < 3;j++){
			val += (-2*x[k][j] + x[k+1][j] - x[k-1][j])/(x_norms[k-1]*x_norms[k]);
			val += -x_angles[k-1]*((x[k][j] - x[k-1][j])/pow(x_norms[k-1],2) + (x[k][j] - x[k+1][j])/pow(x_norms[k],2)) ;
		}
	}
	else{
		for(int j = 0;j < 3;j++){
			val += (x[k][j] - x[k-1][j])/(x_norms[k-1]*x_norms[k]) - x_angles[k-1]*(x[k+1][j] - x[k][j])/pow(x_norms[k],2);
		}
	}
	return val;
}


void G_Matrix_slim::calc_d_x_mat(){
	int i = 0;
	double dxtmp = 0;
	if(interactions > 1){
		for(int k = 1;k <= interactions-1;k++){
			dxtmp = 0;
			for(int i = 0;i < 3;i++) dxtmp += x_derivatives(k,i);
			d_x_mat[k-1] = dxtmp;
			//if(leading_order && ei[0] < 1200) cout <<"dxmat "<< k-1 << " " <<d_x_mat[k-1] << endl;
		}
	}
	else{
		dxtmp = 0;
		//dxtmp += x_derivatives(1,0);
		//dxtmp += x_derivatives(1,2);
		d_x_mat[0] = 0;//dxtmp;
	}
}

void G_Matrix_slim::calc_d_x_mat(int k){
    int iter = 0;
    int j = 1;
    int j_iter = 0;
    if(interactions > 1){
        for(int i = 0;i < interactions*3;i++){
            if(i < (k+3)*3){
                if(i < 3){
                    d_x_mat[i] = (x[j][iter] - x[j+1][iter])/(x_norms[j-1]*x_norms[j])
                                + ((x[j][iter] - x[j-1][iter])*x_scalars[j-1])/(pow(x_norms[j-1],3.)*x_norms[j]);
                    iter += 1;
                }
                else if(i < (k+3)*3 - 3){
                    d_x_mat[i] = (x[j+1][iter] - x[j-1][iter])/(x_norms[j-1]*x_norms[j])
                                - x_scalars[j-1]/pow(x_norms[j-1]*x_norms[j],3.)
                                *(-pow(x_norms[j-1],2.)*(x[j+1][iter] - x[j][iter]) 
                                + pow(x_norms[j],2.)*(x[j][iter] - x[j-1][iter]));
                    iter += 1;                     
                }
                if((k+3)*3 < dim_m - 3 && i >= (k+3)*3 - 3){
                        d_x_mat[i] = (x[j][iter] - x[j-1][iter])/(x_norms[j-1]*x_norms[j])
                                    + ((x[j+1][iter] - x[j][iter])*x_scalars[j-1])/(pow(x_norms[j],3.)*x_norms[j-1]);
                        iter += 1;
                        j_iter += 1;
                }
                
                if(iter == 3){
                    iter = 0;
                    if(j_iter == 3){
                        j_iter = 0;
                        j += 1;
                    }
                }
            }
            else d_x_mat[i] = 0.;
            }
    }
    else{
        for(int i = 0;i < interactions*3;i++){
            if(i < (k+3)*3){
                d_x_mat[i] = (x[0][iter] - x[1][iter])/(x_norms[0])
                            + ((x[0][iter] - x[1][iter])*x_scalars[0])/pow(x_norms[0],3.);
                iter += 1;
                if(i % 3 == 0 && i > 0) iter = 0;
            }
            else d_x_mat[i] = 0.;
        }
    }
}


void G_Matrix_slim::set_Edeps(){
	for(int i = 0;i < interactions-1;i++){
		Edep[i]= ei[i]*(1. - 1./(1. + ei[i]/mc2*(1. - x_angles[i])));
		ei[i+1] = ei[i] - Edep[i];
	}
}



void G_Matrix_slim::set_angles(){
	for(int i = 0;i < interactions-1;i++){
		x_angles[i] = x_scalars[i]/(x_norms[i]*x_norms[i+1]);
	}
}

void G_Matrix_slim::calc_scalars(){
    for(int i = 1;i < interactions;i++){
        x_scalars[i-1] = 0.;
        for(int j = 0;j < 3;j++){
            x_scalars[i-1] += (x[i+1][j] - x[i][j])*(x[i][j] - x[i-1][j]);
        }
    }   
}
void G_Matrix_slim::calc_norms(){
    for(int i = 1;i < interactions+1;i++){
        x_norms[i-1] = 0.;
        for(int j = 0;j < 3;j++){
            x_norms[i-1] += pow(x[i][j] - x[i-1][j],2.);
        }
        x_norms[i-1] = sqrt(x_norms[i-1]);
    }
    for(int i = 0;i < interactions;i++){
        if(x_norms[i] != x_norms[i]) cout << "NaN in norms" << endl;
    }
}


void G_Matrix_slim::set_derivatives(){
    calc_scalars();
    calc_norms();
    //for(int i = 0;i < interactions-1;i++) calc_d_x_mat(i);
	calc_d_x_mat();
}

double G_Matrix_slim::get_delta_x(int i){
	//return delta_x;
	i += 1;
	double return_val = delta_x;
	//if(leading_order && ei[0] < 1200) cout << delta_x << " " <<  delta_x_arr[interactions-1]  << " " << delta_x*delta_x_arr[interactions-1]<< endl;
	//return_val = delta_x*delta_x_arr[interactions-1];
	//else return delta_x;
	
	double distance_to_nominal = 0;
	for(int j = 0;j < 3;j++){
		distance_to_nominal += pow(x[i][j],2);
	}
	distance_to_nominal = sqrt(distance_to_nominal);
	if(distance_to_nominal < 0.235 + 0.02) return_val *= 3;
	return return_val;
}

void G_Matrix_slim::set_sigmas(){
	int iter = 0;
	//errors for angles from x[i]'s
	double angle_err_x[10] = {0.};
	double sin_angle = 0;
	for(int i = 0;i < interactions-1;i++) angle_err_x[i] = 0;
	int i = 0;
	double temp_dist = 0;
	for(int i = 0;i < interactions-1;i++){
		//cout << "dmax " << d_x_mat[iter] << " dx " << delta_x << endl;
		//angle_err_x[i] = pow(1/sqrt(1. - pow(x_angles[i],2))*d_x_mat[iter]*get_delta_x(i),2.);
        sin_angle = sin(acos(x_angles[i]));
        temp_dist = get_delta_x(i);
        angle_err_x[i] = pow(d_x_mat[i]/sin_angle*temp_dist,2);
        //if(leading_order && ei[0] < 1200)cout <<"an_err "<<sqrt(angle_err_x[i]) << " angle " << acos(x_angles[i])*180./M_PI << " sin " << sin_angle <<" dx "<<temp_dist<< endl;
	}
	//errors for angles from Edep[i]'s
	double angle_err_e[10] = {0.};
	for(int i = 0;i < interactions-1;i++){
		if(i > 0) angle_err_e[i] = pow(d_e_mat[i],2); //angle_err_e[i-1] 
		else angle_err_e[i] = pow(d_e_mat[i],2);
		//if(leading_order && ei[0] < 1200)cout << "derr " << angle_err_e[i] << " demat " <<  d_e_mat[i] << endl;
	}
	//if(leading_order && ei[0] < 1400) cout << "errs" << endl;
	for(int i = 0;i < interactions-1;i++){
		delta_theta_arr[i] = sqrt(angle_err_e[i] + angle_err_x[i]);
		//if(leading_order && ei[0] < 1400) cout << angle_err_e[i]*180/M_PI << " " << angle_err_x[i]*180/M_PI << endl;
	}
	//if(leading_order && ei[0] < 1400)cout << "..." << endl;
}
/*
void G_Matrix_slim::set_sigmas(){
    
    int iter = 0;
    for(int i = 0;i < interactions-1;i++){
        
        if(i == 0) delta_theta_arr[i] = 0.;
        else delta_theta_arr[i] = delta_theta_arr[i-1];
        
        for(int j = 0;j < 3;j++){
            delta_theta_arr[i] += pow(d_x_mat[iter]*delta_x,2.);
            iter += 1;
        }
        delta_theta_arr[i] += pow(d_e_mat[i],2);
    }   
    for(int i = 0;i < interactions-1;i++) delta_theta_arr[i] = sqrt(delta_theta_arr[i]);
}
*/
void G_Matrix_slim::set_sigmas(bool a){
    
    int iter = 0;
    for(int i = 0;i < interactions-1;i++){
        if(i == 0){
            sigma_vec[i][0] = 0.;
            sigma_vec[i][1] = 0.;
        }
        else{
            sigma_vec[i][0] = sigma_vec[i-1][0];
            sigma_vec[i][1] = sigma_vec[i-1][1];
        }
        for(int j = 0;j < 3;j++){
            sigma_vec[i][0] += pow(d_x_mat[iter]*delta_x,2.);
            iter += 1;
        }
    }
    calc_d_e_mat(true);
}

inline double G_Matrix_slim::delta_ei(int i,double c_factor){
	double a = 1. - c_factor + ei[i-1]/mc2*(1. - x_angles[i])*pow(c_factor,2);
	double b = Edep[i-1]/mc2*(1. - x_angles[i])*pow(c_factor,2);
	
	return a - b;
}

inline double G_Matrix_slim::delta_edepi(int i,double c_factor){
	double a = -ei[i-1]/mc2*(1. - x_angles[i])*pow(c_factor,2);
	double b = 1 - c_factor - Edep[i-1]/mc2*(1. - x_angles[i])*pow(c_factor,2);
	
	return a - b;
}

inline double G_Matrix_slim::delta_theta(int i,double c_factor){
	return pow(ei[i-1] - Edep[i-1],2)/mc2*pow(c_factor,2);
}

inline double G_Matrix_slim::set_error_ei(int i){
	double error = 0;
	
	if(i == 0) return 0.;
	
	for(int j = 0;j < i;j++){
		error += pow(delta_edep_arr[j],2.);
	}
	return sqrt(error);
}

void G_Matrix_slim::calc_d_e_mat(bool a){
	double c_factor = 0;
	double delta_edep_old = 0;
	
	for(int i = 0;i < interactions - 1;i++){
		if(i == 0){
			c_factor = pow(1. + ei[i]/mc2*(1. - x_angles[i]),-1);
			delta_edep_old = abs(ei[i]*ei[i]/mc2*pow(c_factor,2));
		}
		else{
			c_factor = pow(1. + (ei[i-1] - Edep[i-1])/mc2*(1. - x_angles[i]),-1);
			delta_edep_old = pow(delta_theta(i,c_factor)*sigma_vec[i][0],2) + pow(delta_edepi(i,c_factor)*delta_edep_old,2.)
							+ pow(delta_ei(i,c_factor)*set_error_ei(i),2);
		}
		delta_edep_old = sqrt(delta_edep_old);
		delta_edep_arr[i] = delta_edep_old;
	}
}


void G_Matrix_slim::calc_d_e_mat(){
    double c_angle, delta_cos;
    double ei_err = 0;
    double eip1_err = 0;
    double sin_angle = 0;
    for(int i = 0;i < interactions-1;i++){
        c_angle = 1. - mc2/ei[i+1] + mc2/ei[i];
        ei_err = mc2*(pow(ei[i] - Edep[i],-2) - pow(ei[i],-2));
        eip1_err = mc2*(pow(ei[i+1],-2) - pow(ei[i+1] + Edep[i],-2));
        delta_cos = pow(ei_err*sigma_e(i),2) + pow(eip1_err*sigma_e(i+1),2);
        //delta_cos = pow(mc2,2)*(pow(sigma_e(i+1),2)/pow(ei[i+1],4) + pow(sigma_e(i),2)/pow(ei[i],4));
        delta_cos = sqrt(delta_cos);
        sin_angle = sin(acos(c_angle));
        d_e_mat[i] = delta_cos/sin_angle;
        //d_e_mat[i] = abs(delta_cos/sqrt(1. - pow(c_angle,2.)));
    }

}

double G_Matrix_slim::sigma_e(int i){
    double val = 0.;
    for(int k = 0;k < i;k++){
        val += pow(tolerance/2.355*Edep[k],2.);
    }
    return sqrt(val);
}
