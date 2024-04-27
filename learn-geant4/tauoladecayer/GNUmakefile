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

ifdef Tauola++_DIR
all: lib bin
  include $(G4INSTALL)/config/binmake.gmk

  INCFLAGS  += -I$(Tauola++_DIR)/include
  LDLIBS    += -L$(Tauola++_DIR)/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaHEPEVT
else
all:
	@echo 'ERROR - Tauola not defined in the environment !'
endif
