SHELL:=/bin/bash
export PATH:=$(PATH):$(PWD)

empty:=
space:=$(empty) $(empty)

ROOT?=../..
BIN_ROOT:=$(ROOT)/bin
SRC_ROOT:=$(ROOT)/src
INC_ROOT:=$(ROOT)/inc
LIB_ROOT:=$(ROOT)/lib


M4:=m4
M4_FLG:=-P

CPP:=cpp
CPP_FLG:=-traditional-cpp -x c
CPP_POST:=sed "/^\#/d"

CC:=gcc
CXX:=g++
DC:=dmd
FC:=gfortran

CC_MPI:=mpicc
CXX_MPI:=mpic++
FC_MPI:=mpif90


LAPACK_PATH:=
LAPACK_INC:=
LAPACK_LIB:=-llapack -lblas

OPENMPI_PATH:=
OPENMPI_INC:=-I/usr/include/openmpi
OPENMPI_LIB:=

PARMETIS_PATH:=$(LIB_ROOT)/ParMetis-3.1.1 $(HOME)/bin/ParMetis-3.1.1
PARMETIS_INC:=$(foreach p, $(PARMETIS_PATH), -I$(p))
PARMETIS_LIB:=$(foreach p, $(PARMETIS_PATH), -L$(p)) -lparmetis -lmetis

UTILS_PATH:=$(LIB_ROOT)/utils
UTILS_INC:=-I$(UTILS_PATH)
UTILS_LIB:=-L$(UTILS_PATH) -lutils

PARMS_PATH:=$(LIB_ROOT)/pARMS_3.3 $(HOME)/bin/pARMS_3.3
PARMS_INC:=$(foreach p, $(PARMS_PATH), -I$(p)/include)
PARMS_LIB:=$(foreach p, $(PARMS_PATH), -L$(p)/lib) -lparms
PARMS_FLG:=-D"USE_MPI" -D"REAL=double" -D"DBL" -D"HAS_BLAS" -D"FORTRAN_UNDERSCORE" -D"VOID_POINTER_SIZE_8"


BUILD:=
CC_FLG:=
CXX_FLG:=
FC_FLG:=
LD_FLG:=

ifeq ($(strip $(BUILD)),DEBUG)
  CC_FLG+=-g -p -O1
  CXX_FLG+=-g -p -O1
  FC_FLG+=-g -p -O1
  FC_FLG+=-fbacktrace
  FC_FLG+=-fbounds-check
  FC_FLG+=-frange-check
  FC_FLG+=-fcheck-array-temporaries
  FC_FLG+=-ffpe-trap=invalid,zero,underflow,overflow # ,denormal,precision
  LD_FLG+=-g -p -O1
endif

ifeq ($(strip $(BUILD)),DEBUG_ZERO)
  CC_FLG+=-g -p -O1
  FC_FLG+=-g -p -O1
  FC_FLG+=-fbacktrace
  FC_FLG+=-fbounds-check
  FC_FLG+=-frange-check
  FC_FLG+=-fcheck-array-temporaries
  FC_FLG+=-ffpe-trap=invalid,zero,underflow,overflow # ,denormal,precision
  FC_FLG+=-finit-real=zero # (zero|inf|-inf|nan)
  FC_FLG+=-finit-local-zero
  FC_FLG+=-finit-integer=0
  LD_FLG+=-g -p -O1
endif

ifeq ($(strip $(BUILD)),PROFILE)
  CC_FLG+=-g -p -O1
  CXX_FLG+=-g -p -O1
  FC_FLG+=-g -p -O1
  LD_FLG+=-g -p -O1
endif

ifeq ($(strip $(BUILD)),RELEASE)
  CC_FLG+=-O3
  CXX_FLG+=-O3
  FC_FLG+=-O3
  LD_FLG+=-O3
endif

ifeq ($(strip $(BUILD)),)
  CC_FLG+=-O2
  CXX_FLG+=-O2
  FC_FLG+=-O2
  LD_FLG+=-O2
endif
