#ifndef PERMCLASS_H
#define PERMCLASS_H

class perm_class{

private:

    int len_x,len_x_all,am_perm,start_pos,iter,true_len;
    
    int** arr_perms;
    int** complete_perms;
    int** arr_true;

public:
    perm_class(int);
    ~perm_class();
    
    void feed_complete(int**,int,int,int);
    bool test_for_double();
    void cout_line();
    void cout_line(int);

    int** get_arr();
    int get_len_x();
};




#endif