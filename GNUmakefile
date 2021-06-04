
ifndef MYSW_DIR
ERROR_MESSAGE := $(error MYSW_DIR is not set... run configure.sh!)
endif

OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

INCFLAGS := -I.

include $(MYSW_DIR)/Makefile/Makefile.${OSNAME}

SUBDIRS := ./icarus_signal_processing ./icarus_signal_processing/Filters ./icarus_signal_processing/Detection
CXXFLAGS+=$(INCFLAGS)

.phony: all clean

all:
	@echo "Start building..."
	@for i in $(SUBDIRS); do ( echo "" && echo "Compiling $$i..." && $(MAKE) --directory=$$i) || exit $$?; done
	@echo "Done!"
clean:
	@echo "Cleaning..."
	@for i in $(SUBDIRS); do ( echo "" && echo "Cleaning $$i..." && $(MAKE) clean --directory=$$i) || exit $$?; done
	@echo "Done!"

