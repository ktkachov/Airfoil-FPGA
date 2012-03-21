ifdef MAXCOMPILERDIR
include $(MAXCOMPILERDIR)/lib/Makefile.include
endif
CC = gcc
CCFLAGS = -Wall -m64 -pedantic -std=c99 -g 
METIS_INC = -I/homes/kt208/metis/metis-install/include
METIS_LDINC = -L/homes/kt208/metis/metis-install/lib
LDFLAGS = -lm -lmetis 
LD = gcc

mesh: mesh_part.c airfoil_utils.h airfoil_kernels.h 
	$(CC) -o $@ $(CCFLAGS) $(METIS_INC) $(METIS_LDINC) $< $(LDFLAGS)


# The mesh_fpga targets should be compiled on a maxstation
mesh_fpga.o: mesh_part.c airfoil_utils.h airfoil_kernels.h
	$(CC) -c -o $@ -DRUN_FPGA $(CCFLAGS) $(MAXCOMPILER_INC) $(METIS_INC) $< 

mesh_fpga: mesh_fpga.o ResCalc.o
	$(LD) -o $@ $^ -lc $(LDFLAGS) $(METIS_LDINC) $(MAXCOMPILER_LIBS)

ResCalc.o: ResCalcSim.max
	$(MAXFILECOMPILE) $^ $@ ResCalc

airfoil: airfoil.cpp
	g++ -o $@ -Wall -Werror -O3 -lm $<

%.svg: %.dot
	./dot2svg $<
	cp $< graphs/

%.dot: mesh
	./$<

graph: meshColoured.svg meshSchedule.svg
	cp $^ graphs/

clean:
	rm -f mesh_fpga
	rm -f mesh
	rm -f airfoil
	rm -f *.svg
	rm -f *.dot
	rm -f *.o
