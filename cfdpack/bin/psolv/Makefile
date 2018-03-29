include ../../make.inc

INC+=$(PARMS_INC) $(OPENMPI_INC) $(LAPACK_INC) $(PARMETIS_INC) $(UTILS_INC)
LIB+=$(PARMS_LIB) $(OPENMPI_LIB) $(LAPACK_LIB) $(PARMETIS_LIB) $(UTILS_LIB)

LIB+=-lgfortran -lm
FLG+=-D"USE_MPI"

srcdir:=
srcdir+=$(shell find $(SRC_ROOT)/base -type d)
srcdir+=$(shell find $(SRC_ROOT)/psolv -type d)

incdir:=
incdir+=$(srcdir)

incdir:=$(subst $(space),$(space)-I,$(incdir))
INC+=-I$(incdir)

vpath
vpath %.cpp $(subst $(space),:,$(srcdir))
vpath %.c   $(subst $(space),:,$(srcdir))
vpath %.f90 $(subst $(space),:,$(srcdir))


SRC:=

SRC+=prec_m.f90
SRC+=constants_m.f90
SRC+=alloc_m.f90
SRC+=vector_utils_m.f90
SRC+=matrix_utils_m.f90
SRC+=graph_m.f90
SRC+=hash_m.f90
SRC+=solv_m.f90
SRC+=geometry_m.f90
SRC+=transport_m.f90
SRC+=psolv_cfd.f90

SRC+=parms_intf.cpp

SRC+=utils.cpp
SRC+=graph.cpp
SRC+=mshf.cpp
SRC+=mesh.cpp
SRC+=xgrph.cpp
SRC+=psolv.cpp

OBJ:=$(addsuffix .o,$(basename $(SRC)))
BIN:=psolv

FC_MPI_LINK:=$(shell $(FC_MPI) -showme:link)


$(BIN): $(OBJ)
	$(CXX_MPI) $(CXX_FLG) $(FLG) -o $@ $^ $(LIB) $(FC_MPI_LINK)

%.o: %.cpp
	$(CXX_MPI) $(CXX_FLG) $(FLG) -c $< $(INC)

%.o: %.c
	$(CC_MPI) $(FC_FLG) $(FLG) -c $< $(INC)

%.o: %.f90
	$(FC_MPI) $(FC_FLG) $(FLG) -c $< $(INC)

PHONY: clean
clean:
	@rm -f *.o *.mod *~ gmon.out

PHONY: distclean
distclean:
	$(MAKE) clean
	@rm -f $(BIN)
	@find ../ -name \*~ -exec rm {} \;
