
all: yamsvc-regionbed fastq-qc yamsvc-caller

yamsvc-regionbed: yamsvc_regionbed.d
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD -of$@ $<

yamsvc-caller: yamsvc_caller.d
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD -of$@ $<

fastq-qc: fastq/test.d
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD -of$@ $<

