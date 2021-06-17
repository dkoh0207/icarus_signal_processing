#
# This is an example GNUmakefile for my packages
#

# specific names for this package
SOURCES = $(wildcard *.cxx) $(wildcard */*.cxx)
# FMWK_HEADERS = LinkDef.h ./Filters/LinkDef.h ./Detection/LinkDef.h
# HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
HEADERS = $(wildcard *.h) $(wildcard */*.h)

# include options for this package
INCLUDE := -I. 
INCFLAGS+=$(INCLUDE)
CXXFLAGS+=$(INCLUDE)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(MYSW_DIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(MYSW_DIR)/Makefile/GNUmakefile.CORE
