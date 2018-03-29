include ../../make.inc

INC:=$(PARMETIS_INC) $(METIS_INC) $(UTILS_INC)
LIB:=$(PARMETIS_LIB) $(METIS_LIB) $(UTILS_LIB)
FLG:=-D"USE_MPI"

srcdir:=
srcdir+=$(shell find $(SRC_ROOT)/base -type d)
srcdir+=$(shell find $(SRC_ROOT)/pgmsh -type d)

incdir:=
incdir+=$(srcdir)

incdir:=$(subst $(space),$(space)-I,$(incdir))
INC+=-I$(incdir)

vpath
vpath %.cpp $(subst $(space),:,$(srcdir))

SRC:=
SRC+=utils.cpp
SRC+=graph.cpp
SRC+=mshf.cpp
SRC+=mesh.cpp
SRC+=xgrph.cpp
SRC+=ghst_comm.cpp
SRC+=pgmsh.cpp

OBJ:=$(addsuffix .o,$(basename $(SRC)))
BIN:=pgmsh


$(BIN): $(OBJ)
	$(CXX_MPI) $(CXX_FLG) $(FLG) -o $@ $^ $(LIB)

%.o: %.cpp
	$(CXX_MPI) $(CXX_FLG) $(FLG) -c $< $(INC)


PHONY: clean
clean:
	@rm -f *.o *~ gmon.out

PHONY: distclean
distclean:
	$(MAKE) clean
	@rm -f $(BIN)
	@find ../ -name \*~ -exec rm {} \;
