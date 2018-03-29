src_path:=src/
lib_path:=lib/


LIB:=
LIB+=json-cpp
LIB+=json-f90
ifeq ($(strip $(mpi)),yes)
  LIB+=ParMetis-3.1.1
  LIB+=pARMS_3.3
endif

LIB:=$(addprefix $(lib_path), $(LIB))
LIB_MAKE:=$(foreach item, $(LIB), $(MAKE) -C $(item);)
LIB_CLEAN:=$(foreach item, $(LIB), $(MAKE) -C $(item) clean;)

.PHONY: all
all: lib-make flowlab

flowlab:
	$(MAKE) -C $(src_path)

flowlab-clean:
	$(MAKE) -C $(src_path) clean

lib-make:
	$(LIB_MAKE)

lib-clean:
	$(LIB_CLEAN)

rm-junk:
	@find ./ -name \*~ -exec rm {} \;
