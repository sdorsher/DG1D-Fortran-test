MODULES = tools.o kinds.o structures.o
#LIBS = /usr/local/lib/liblapack.a
F90 = gfortran
#F90 = gfortran

F90FLAGS = -O3 -g

all:	driver.x schw_wave.x conf_wave.x

driver.x:	driver.o $(MODULES)
	$(F90) $(F90FLAGS) -L /usr/local/lib/ -lblas -llapack driver.o $(MODULES) -o driver.x

driver.o:	driver.f90 $(MODULES)
	$(F90) $(F90FLAGS) -c driver.f90

schw_wave.x:	schw_wave.o $(MODULES) schw_physics.o
	$(F90) $(F90FLAGS) -L /usr/local/lib/ -lblas -llapack schw_wave.o $(MODULES) schw_physics.o $(LIBS) -o schw_wave.x

schw_wave.o:	schw_wave.f90 $(MODULES) schw_physics.o
	$(F90) $(F90FLAGS) -c schw_wave.f90

conf_wave.x:	conf_wave.o $(MODULES) schw_physics.o
	$(F90) $(F90FLAGS) -L /usr/local/lib/ -lblas -llapack conf_wave.o $(MODULES) schw_physics.o $(LIBS) -o conf_wave.x

conf_wave.o:	conf_wave.f90 $(MODULES) schw_physics.o
	$(F90) $(F90FLAGS) -c conf_wave.f90

kinds.o:	kinds.f90
	$(F90) $(F90FLAGS) -c kinds.f90

tools.o:	tools.f90 structures.o kinds.o
	$(F90) $(F90FLAGS) -c tools.f90

structures.o:	structures.f90 kinds.o
	$(F90) $(F90FLAGS) -c structures.f90

schw_physics.o:	schw_physics.f90 kinds.o
	$(F90) $(F90FLAGS) -c schw_physics.f90

clean:
	rm -rf *.o *.mod *.x
