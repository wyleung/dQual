
all: yamsvc-regionbed

main: main.d
	rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/BioD $<
    # rdmd --force --compiler=ldmd2  --build-only  -O2 -release  -cov -I`pwd`/sambamba -I`pwd`/sambamba/BioD hello.d

yamsvc-regionbed: main
	cp $< $@

