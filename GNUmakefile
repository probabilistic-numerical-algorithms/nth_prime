.DELETE_ON_ERROR:
.DEFAULT_GOAL := default

PKGCONF := pkgconf
GMP_LIBS := $(shell $(PKGCONF) --libs gmp 2>/dev/null || echo "")

MKDIR_P = mkdir -p

DC = gdc
DFLAGS =

CC = gcc
CFLAGS =

CFLAGS += -std=gnu23
CFLAGS += -O3
CFLAGS += -march=native -mtune=native
CFLAGS += -fopenmp
CFLAGS += -funroll-loops
CFLAGS += -flto=auto

DFLAGS += -frelease
DFLAGS += -march=native -mtune=native
DFLAGS += -O3
DFLAGS += -funroll-loops
DFLAGS += -flto=auto

LIMBS_SIZE = LIMBS_256

.PHONY: default all clean
default: all

nth-prime-c: c23/nth_prime.c
	$(CC) $(CFLAGS) -DSTANDALONE=1 -D$(LIMBS_SIZE) \
	    -o $(@) $(<) -lm $(GMP_LIBS)
all:: nth-prime-c
clean::
	-rm -f nth-prime-c

nth-prime-64-c: c23/nth_prime_64.c
	$(CC) $(CFLAGS) -DSTANDALONE=1 \
	    -o $(@) $(<) -lm $(GMP_LIBS)
all:: nth-prime-64-c
clean::
	-rm -f nth-prime-64-c

nth-prime-d: dlang/nth_prime.d
	$(DC) $(DFLAGS) -fversion=standalone \
	    -fversion=$(LIMBS_SIZE) \
	    -o $(@) $(<) -lm $(GMP_LIBS)
all:: nth-prime-d
clean::
	-rm -f nth-prime-d

nth-prime-64-d: dlang/nth_prime_64.d
	$(DC) $(DFLAGS) -fversion=standalone \
	    -o $(@) $(<) -lm $(GMP_LIBS)
all:: nth-prime-64-d
clean::
	-rm -f nth-prime-64-d
