CC      = g++
CFLAGS  = -g -pthread -O3 -std=c++11 -Wall
LIBFLAGS = -I ~/eigen_dir/

OBJECTS = Detector_Geometry.o threshold_handler.o max_prob_class.o prior_handler.o RLFCM.o air_path.o thread_handle.o G_Matrix_slim.o merge_class.o cluster_perms.o perm_class.o cluster.o data_class.o geometry.o permutations.o e0_optimization.o probabilities_c.o fuzzy_c_means.o cross_sec.o

all: Program

%.o: %.cc
	$(CC) -c $< $(CFLAGS) $(LIBFLAGS)
	
Program: $(OBJECTS) main.o
	$(CC) $(OBJECTS) main.o -pthread -o fbt_agata_all

clean: 
	rm -f *.o $(PROGRAM)
