# --------------------------------------------------------------
# GNUmakefile for examples module.
# --------------------------------------------------------------
name := tauola_decayer
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all

ifdef TAUOLALOCATION
all: lib bin
  include $(G4INSTALL)/config/binmake.gmk

  INCFLAGS  += -I$(TAUOLALOCATION)/include
  LDLIBS    += -L$(TAUOLALOCATION)/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaHEPEVT
else
all:
	@echo 'ERROR - Tauola not defined in the environment !'
endif
