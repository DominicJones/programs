include ../make.inc


LIB:=
LIB+=-lgfortran -lm


srcdir:=
srcdir+=$(shell find $(SRC_ROOT) -type d)


vpath
vpath %.c $(subst $(space),:,$(srcdir))
vpath %.f90 $(subst $(space),:,$(srcdir))


SRC1:=
SRC1+=const_m.f90  alloc_m.f90  sparse_m.f90 sort_m.f90  hash_m.f90  list_m.f90
SRC1+=mesh_data_m.f90  geom_data_m.f90  pde_data_m.f90
SRC1+=mesh_format_m.f90  read_mesh_m.f90  connectivity_m.f90  write_field_m.f90  write_mesh_m.f90
SRC1+=vector_utils_m.f90  matrix_utils_m.f90  misc_utils_m.f90
SRC1+=geometry_m.f90  gradient_m.f90  limiters_m.f90
SRC1+=pde_utils_m.f90  lin_eqns_solvers_m.f90  turbulence_m.f90  spalart_allmaras_m.f90
SRC1+=transport_m.f90  pdes_m.f90  algorithms_m.f90
SRC1+=gpde_main.f90

SRC2:=
SRC2+=main.c

SRC:=
SRC+=$(SRC1) $(SRC2)


OBJ:=$(addsuffix .o,$(basename $(SRC)))

BIN:=gpde


$(BIN): $(OBJ)
	$(CC) $(LD_FLG) $(FLG) -o $@ $^ $(LIB)


%.o: %.f90
	$(FC) $(FC_FLG) $(FLG) -c $< $(INC)


%.o: %.c
	$(CC)  $(CC_FLG) $(FLG) -c $< $(INC)


.PHONY: clean
clean:
	@rm -f *.o *.mod *~ gmon.out


.PHONY: distclean
distclean:
	$(MAKE) clean
	@rm -f $(BIN)
	@find ../ -name \*~ -exec rm {} \;
