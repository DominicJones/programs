SRC = ../src


LIB = ../lib
LADIR = $(LIB)/lapack
LA = $(LADIR)/lapack95.a \
     $(LADIR)/lapack.a \
     $(LADIR)/blas.a \
     $(LADIR)/tmglib.a
LAMOD = $(LADIR)/lapack95_modules

FTN = gfortran
EXT = f90
OPT = -O3
EXE = spmom

PREC = precision.o
SOLV = solv_data.o solv_res.o solv_precond.o solv_direct.o solv_bicg.o solv_cg.o
MOM = mom_data.o mom_math.o mom_closure.o mom_solve.o
CFD1 = cfd_data.o
GRD = grd_data.o grd_solv.o grd_math.o grd_constr.o grd_read.o grd_conn.o grd_geom.o
CFD2 = cfd_constr.o cfd_outp.o
SPR = spr_data.o spr_inject.o spr_impact.o spr_hydro.o spr_conv.o
CFD3 = cfd_define.o cfd_setup.o cfd_solve.o
MAIN = main.o

$(EXE): $(PREC) $(SOLV) $(MOM) $(CFD1) $(GRD) $(CFD2) $(SPR) $(CFD3) $(MAIN)
	$(FTN) $(OPT) -o $@ $^  $(LA)

%.o: $(SRC)/%.$(EXT)
	$(FTN) $(OPT) -c $<  -I $(LAMOD)

clean:
	rm -f $(SRC)/*~ *~ *.o *.mod $(EXE)
