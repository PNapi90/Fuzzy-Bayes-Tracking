#ifndef PERMCLASS_H
#define PERMCLASS_H

class perm_class{

private:
	
	const bool save_mode = false;
	
    int len_x,len_x_all,am_perm,start_pos,iter,true_len;
    int all_perms_len;
    
    int interactions;
	int sub_cl_len;
	int position_x;
	int position_y;

    int* recurrent_arr;
    double* prob_arr;
    
    int** arr_perms;
    int** complete_perms;
    int** arr_true;

    void generate_arr_perms();
    void get_recurrent_arr();
	
	void save_recurrent_arr();
	
    inline bool check_all(int,int);

public:
    perm_class(int,int);
    ~perm_class();
    
    void feed_complete(int**,int,int);
    void set_probability(double,int);

    bool test_for_double();
    bool check_recurrent(int);

    void cout_line();
    void cout_line(int);
    void set_id(int,int,int,int);

    double get_probability(int);
    int return_line(int,int);

    int** get_arr();
    int get_len_x();
};




#endif
