COMPILER = gfortran
FLAGS = -O3

OBJS = constants.o mdSubroutines.o analysisSubroutines.o

equilSim: equilWater.f90 testMod.o $(OBJS)
	$(COMPILER) $(FLAGS) -o equilSim equilWater.f90 testMod.o $(OBJS)

analSim: analSim.f90 $(OBJS)
	$(COMPILER) $(FLAGS) -o analSim analSim.f90 $(OBJS)

waterBox: waterBox.f90
	$(COMPILER) $(FLAGS) -o waterBox waterBox.f90

drawVelDist: drawVelDist.f90
	$(COMPILER) $(FLAGS) -o drawVelDist drawVelDist.f90

constants.o: constants.f90
	$(COMPILER) $(FLAGS) -c constants.f90

mdSubroutines.o: mdSubroutines.f90 constants.o
	$(COMPILER) $(FLAGS) -c mdSubroutines.f90

analysisSubroutines.o: analysisSubroutines.f90 constants.o
	$(COMPILER) $(FLAGS) -c analysisSubroutines.f90

testMod.o: testMod.f90 constants.o
	$(COMPILER) $(FLAGS) -c testMod.f90
