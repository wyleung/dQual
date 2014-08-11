######
# Automatic building of the yamsvc programs
######


uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
uname_O := $(shell sh -c 'uname -o 2>/dev/null || echo not')
uname_R := $(shell sh -c 'uname -r 2>/dev/null || echo not')
uname_P := $(shell sh -c 'uname -p 2>/dev/null || echo not')
uname_V := $(shell sh -c 'uname -v 2>/dev/null || echo not')

PWD := $(shell pwd)
BUILD_DIR := $(PWD)/builds/$(uname_S)_$(uname_M)
OUT_DIR := $(PWD)/bin

INSTALL := install

PROGRAMS := yamsvc-regionbed fastq-seqstat yamsvc-caller
TARGET_BINS := $(addprefix $(BUILD_DIR)/, $(PROGRAMS))

all: $(TARGET_BINS)
install: all
	$(INSTALL) -d -m 755 '$(OUT_DIR)'
	$(INSTALL) $(TARGET_BINS) '$(OUT_DIR)'

$(BUILD_DIR): 
	mkdir -p $@

$(BUILD_DIR)/yamsvc-regionbed: yamsvc_regionbed.d yamsvc/utils.d | $(BUILD_DIR)
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD -of$@ $<

$(BUILD_DIR)/yamsvc-caller: yamsvc_caller.d yamsvc/utils.d | $(BUILD_DIR)
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD -of$@ $<

$(BUILD_DIR)/fastq-seqstat: fastq/seqstat.d yamsvc/utils.d | $(BUILD_DIR)
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD -of$@ $<

