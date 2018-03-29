include ../etc/make/Makefile.variables

target_path:=../bin/
target:=flowlab.exe

lib_path:=../lib
include ../lib/Makefile.lib

f90_flags+=-O2 -cpp
cpp_flags+=-O2

ifeq ($(strip $(mpi)),yes)
  f90_compiler:=$(f90_mpi_compiler)
  f90_flags+=-D"USE_MPI"

  cpp_compiler:=$(cpp_mpi_compiler)
  cpp_flags+=-D"USE_MPI"

  libraries+=$(shell $(f90_mpi_compiler) -showme:link)
endif

libraries+=-lgfortran

f90_sources:=
f90_sources+=prec_m.f90
f90_sources+=constants_m.f90
f90_sources+=alloc_m.f90
f90_sources+=vector_utils_m.f90
f90_sources+=matrix_utils_m.f90
f90_sources+=graph_m.f90
f90_sources+=hash_m.f90
f90_sources+=solv_m.f90
f90_sources+=geometry_m.f90
f90_sources+=transport_m.f90
f90_sources+=navierstokes.f90

# cpp_sources:=$(wildcard *.cpp)
cpp_sources:=
cpp_sources+=map.cpp
cpp_sources+=ghst_comm.cpp
cpp_sources+=parms_intf.cpp
cpp_sources+=utils.cpp
cpp_sources+=graph.cpp
cpp_sources+=mshf.cpp
cpp_sources+=mesh.cpp
cpp_sources+=xgrph.cpp
cpp_sources+=partition.cpp
cpp_sources+=main.cpp

sources:=$(f90_sources) $(cpp_sources)

objects:=$(addsuffix .o,$(basename $(sources)))
obj_dep_rules:=$(cpp_sources:.cpp=.mk)

objects:=$(addprefix $(target_path), $(objects))
obj_dep_rules:=$(addprefix $(target_path), $(obj_dep_rules))
target:=$(addprefix $(target_path), $(target))

make_target_path:=mkdir -p $(target_path)

.PHONY: clean default

default: $(target)

$(target): $(objects)
	@$(make_target_path)
	$(cpp_compiler) $(cpp_flags) $^ -o $@ $(libraries)

$(target_path)%.o: $(target_path)%.mk
	@$(make_target_path)
	$(cpp_compiler) $(cpp_flags) -c $(basename $*).cpp -o $@ $(includes)

$(target_path)%.mk: %.cpp
	@$(make_target_path)
	$(cpp_compiler) $(cpp_flags) -MM $< >> $@ $(includes)

$(target_path)%.o: %.f90
	@$(make_target_path)
	$(f90_compiler) $(f90_flags) -c $< -o $@ $(includes) -J$(target_path)

clean:
	@rm -rf $(target_path)

ifneq ($(MAKECMDGOALS), clean)
-include $(obj_dep_rules)
endif
