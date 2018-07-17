#include <iostream>
#include <string>
#include <thread>

#include "data_class.h"
#include "cluster.h"

#ifndef THREAD_HANDLE
#define THREAD_HANDLE


class thread_handle{

private:

    const int amount = 100000;
    int thread_num;
    int verbosity;
    bool solve_type_bool;

    std::string input, output;

    data_handler* data;
	
	bool abort_statement(int);
	
    void do_calculations();

public:
    thread_handle(std::string,std::string,int,int,bool);
    ~thread_handle();
    
    void call_title(int);

    std::thread threading();
};

#endif
