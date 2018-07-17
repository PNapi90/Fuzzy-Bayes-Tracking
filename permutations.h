#include <fstream>
#include <algorithm>

#ifndef PERMUTE_H
#define PERMUTE_H

class permutations{

private:

    const int max_interactions = 10;

    int leni;
    int* remain_lens;
    int** arr;
    int** am_clusters;
    int*** remain_arr;
    int*** full_perm_arr;
    int** max_perms;
    int**** permutation_sub_arr;

    void load_poss_combs();
    void get_max_perms();
    void calc_remain_combs(int);
    void save_important_remains(int);
    void set_full_perm_arr(int);

    inline int k_out_n(int,int);
    inline int factorial(int a);


public:
    permutations();
    ~permutations();

    int** return_remains(int);
    int*** get_perms(int);
    int* return_max_perms(int);
    int* return_am_clusters(int);
    int return_am_clusters(int,int);
    int** get_perms_full(int);

    int return_remains_len(int);

    void cout_perms(int);
};

#endif
