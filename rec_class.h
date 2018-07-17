#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

#ifndef R_CLASS
#define R_CLASS

class rec_class{

private:
    int am;
    int** arr;
    int leni;
    //std::ofstream** out;

public:
    rec_class(int);
    ~rec_class();
    
    void bars(int);
    void cout_it();

};



#endif