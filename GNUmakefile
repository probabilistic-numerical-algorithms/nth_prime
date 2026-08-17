.DELETE_ON_ERROR:
.DEFAULT_GOAL := default

MKDIR_P = mkdir -p

CC = gcc
CFLAGS =

DC = gdc
DFLAGS =

FC = gfortran
FFLAGS =

CFLAGS += -std=gnu23
CFLAGS += -O3
CFLAGS += -march=native -mtune=native
CFLAGS += -fopenmp
CFLAGS += -funroll-loops
CFLAGS += -flto=auto

DFLAGS += -frelease
DFLAGS += -O3
DFLAGS += -march=native -mtune=native
DFLAGS += -funroll-loops
DFLAGS += -flto=auto

FFLAGS += -std=f2023
FFLAGS += -O3
FFLAGS += -march=native -mtune=native
FFLAGS += -fopenmp
FFLAGS += -funroll-loops
FFLAGS += -flto=auto

.PHONY: default all clean
default: all

nth-prime-64-c: c23/nth_prime_64.c
	$(CC) $(CFLAGS) -DSTANDALONE=1 $(<) -lm -o $(@)
all:: nth-prime-64-c
clean:: ; -rm -f nth-prime-64-c

nth-prime-64-d: dlang/nth_prime_64.d
	$(DC) $(DFLAGS) -fversion=standalone $(<) -lm -o $(@)
all:: nth-prime-64-d
clean:: ; -rm -f nth-prime-64-d

nth-prime-64-f: $(patsubst %,f2023/%.f90, nth_prime_64 types_mod main)
	$(FC) $(FFLAGS) $(^) -o $(@)
all:: nth-prime-64-f
clean:: ; -rm -f nth-prime-64-f
clean:: ; -rm -f nth_prime_64_mod.mod types_mod.mod
