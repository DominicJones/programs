include ../../etc/make/Makefile.variables

target:=jsonf90

sources:=
sources+=fson_string_m.f90
sources+=fson_value_m.f90
sources+=fson_path_m.f90
sources+=fson.f90

objects:=$(addsuffix .o,$(basename $(sources)))

include ../../etc/make/Makefile.staticlib

include ../../etc/make/Makefile.compile
