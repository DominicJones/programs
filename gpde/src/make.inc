SHELL:=/bin/bash
export PATH:=$(PATH):$(PWD)


empty:=
space:=$(empty) $(empty)


# ROOT is set with respect to Makefile in bin/
ROOT:=..


SRC_ROOT:=$(ROOT)/src
INC_ROOT:=$(ROOT)/inc
LIB_ROOT:=$(ROOT)/lib


CC:=gcc
FC:=gfortran


LAPACK_PATH:=
LAPACK_INC:=
LAPACK_LIB:=-llapack -lblas


BUILD:=
CC_FLG:=
FC_FLG:=
LD_FLG:=

ifeq ($(strip $(BUILD)),DEBUG)
  CC_FLG+=-g -p -O1
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
  FC_FLG+=-g -p -O1
  LD_FLG+=-g -p -O1
endif

ifeq ($(strip $(BUILD)),RELEASE)
  CC_FLG+=-O3
  FC_FLG+=-O3
  LD_FLG+=-O3
endif

ifeq ($(strip $(BUILD)),)
  CC_FLG+=-O2
  FC_FLG+=-O2
  LD_FLG+=-O2
endif
