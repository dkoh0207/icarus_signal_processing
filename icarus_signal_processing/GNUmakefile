#
# This is an example GNUmakefile for my packages
#

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

# include options for this package
INCFLAGS  = -I. \
			-I/sdf/group/neutrino/koh0207/DSP/icarus_signal_processing/icarus_signal_processing/Filters \
			-I/sdf/group/neutrino/koh0207/DSP/icarus_signal_processing/icarus_signal_processing/Detection
CXXFLAGS=-Werror

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(MYSW_DIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(MYSW_DIR)/Makefile/GNUmakefile.CORE
