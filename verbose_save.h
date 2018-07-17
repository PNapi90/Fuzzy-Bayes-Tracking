#include <fstream>
#include <string>

#ifndef  VERB_SAVE_H
#define VERB_SAVE_H

class verbose_save{

private:

    int N;
    int iter;
    double*** integrals;

    std::ofstream*** integrals_save;



public:
    verbose_save(int);
    ~verbose_save();
    
};





#endif