######
# Automatic building of the yamsvc programs
######

D_COMPILER = ldmd2
D_FLAGS =--compiler=$(D_COMPILER) -Ithirdparty/BioD -O2 -debug -inline -cov -Ithirdparty/sambamba -Isource

THIS_MAKEFILE := $(lastword $(MAKEFILE_LIST))
MAKEFILE_DIR := $(realpath $(dir $(realpath $(THIS_MAKEFILE))))

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
uname_O := $(shell sh -c 'uname -o 2>/dev/null || echo not')
uname_R := $(shell sh -c 'uname -r 2>/dev/null || echo not')
uname_P := $(shell sh -c 'uname -p 2>/dev/null || echo not')
uname_V := $(shell sh -c 'uname -v 2>/dev/null || echo not')

PWD := $(MAKEFILE_DIR)
BUILD_DIR := $(abspath $(MAKEFILE_DIR)/../build/$(uname_S)_$(uname_M))
OUT_DIR := $(abspath $(MAKEFILE_DIR)/../bin)

INSTALL := install

PROGRAMS := yamsvc-regionbed fastq-seqstat yamsvc-caller
TARGET_BINS := $(addprefix $(BUILD_DIR)/, $(PROGRAMS))

INCL_OBJS := source/yamsvc/utils.d


all: $(TARGET_BINS)
install: all
	$(INSTALL) -d -m 755 '$(OUT_DIR)'
	$(INSTALL) $(TARGET_BINS) '$(OUT_DIR)'

$(BUILD_DIR): 
	mkdir -p $@

$(BUILD_DIR)/yamsvc-regionbed: source/yamsvc_regionbed.d $(INCL_OBJS) | $(BUILD_DIR)
	rdmd --force $(D_FLAGS) -of$@ $<

$(BUILD_DIR)/yamsvc-caller: source/yamsvc_caller.d $(INCL_OBJS) | $(BUILD_DIR)
	rdmd --force $(D_FLAGS) -of$@ $<

$(BUILD_DIR)/fastq-seqstat: source/fastq/seqstat.d $(INCL_OBJS) | $(BUILD_DIR)
	rdmd --force $(D_FLAGS) -of$@ $<

.PHONY: clean

clean:
	rm -rf $(BUILD_DIR)