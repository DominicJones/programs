include ../etc/make/Makefile.variables

target_path:=../bin/
target:=$(target_name).exe

lib_path:=../lib
include ../lib/Makefile.lib

includes+=$(JSONCPP_INC)
includes+=$(JSONF90_INC)

libraries+=$(JSONCPP_LIB)
libraries+=$(JSONF90_LIB)

f90_flags+=-cpp
cpp_flags+=

ifeq ($(strip $(mpi)),yes)
  f90_compiler:=$(f90_mpi_compiler)
  f90_flags+=-D"USE_MPI"

  cpp_compiler:=$(cpp_mpi_compiler)
  cpp_flags+=-D"USE_MPI"

  libraries+=$(shell $(f90_mpi_compiler) -showme:link)
endif


build:=

ifeq ($(strip $(build)),debug)
  cpp_flags+=-g -p -O0
  f90_flags+=-g -p -O0
  f90_flags+=-fbacktrace
  f90_flags+=-fbounds-check
  f90_flags+=-frange-check
  f90_flags+=-fcheck-array-temporaries
  f90_flags+=-ffpe-trap=invalid,zero,underflow,overflow # ,denormal,precision
  f90_flags+=-finit-real=zero # (zero|inf|-inf|nan)
  f90_flags+=-finit-local-zero
  f90_flags+=-finit-integer=0
endif

ifeq ($(strip $(build)),)
  cpp_flags+=-O3
  f90_flags+=-O3
endif





libraries+=-lgfortran
libraries+=-static-libstdc++

# cpp_compiler:=$(GCC_HOME)/bin/g++ -std=c++14 -static-libstdc++
# libraries+=-L$(dir $(shell gcc --print-file-name libgfortran.so))
# libraries+=-L$(dir $(shell $(GCC_HOME)/bin/g++ --print-file-name libgfortran.so))
# libraries+=-L$(dir $(shell gcc --print-file-name libc.so))
# libraries+=-dynamiclib -lgfortran
# libraries+=-static -lgfortran
# libraries+=-static-libgfortran
# libraries+=-static-libgcc

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
cpp_sources+=machine.cpp
cpp_sources+=parms_intf.cpp
cpp_sources+=graph.cpp
cpp_sources+=mshf.cpp
cpp_sources+=mesh.cpp
cpp_sources+=geometry.cpp
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
	@rm -f $(target_path)/*

ifneq ($(MAKECMDGOALS), clean)
-include $(obj_dep_rules)
endif
