#include <algorithm>

#ifndef HOLD_ALL_PERMS_H
#define HOLD_ALL_PERMS_H

class hold_all_perms{

private:

    int* len_all_perms;
    int*** hold_all_perms_arr;

    void do_permutation(int);
    void get_all_perms_lens();
    void get_all_perms();

public:
    hold_all_perms();
    ~hold_all_perms();

    int** return_perm_arr(int);
    
};

#endif