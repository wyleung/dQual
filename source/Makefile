
PWD := $(shell pwd)
OUT_DIR := $(PWD)/bin

TARGET := yamsvc-regionbed fastq-seqstat yamsvc-caller

all: $(addprefix $(OUT_DIR)/, $(TARGET))

$(OUT_DIR)/yamsvc-regionbed: yamsvc_regionbed.d
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD -of$@ $<

$(OUT_DIR)/yamsvc-caller: yamsvc_caller.d
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD -of$@ $<

$(OUT_DIR)/fastq-seqstat: fastq/seqstat.d
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD -of$@ $<

