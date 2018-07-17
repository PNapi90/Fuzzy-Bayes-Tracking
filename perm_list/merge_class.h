#include "perm_class.h"


#ifndef MERGE_CLASS_H
#define MERGE_CLASS_H

class merge_class{

private:

    int am_clust, full_len, am_perm;
    int* lens;

    int** true_cluster;

    perm_class** A;

public:
    merge_class(int*,int,int,perm_class**);
    ~merge_class();
    
    void compare_true();
};




#endif