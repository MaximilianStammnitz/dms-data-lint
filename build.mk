# dms_data_lint - integrity and format checks of paired-end FASTQ files (nf-core/deepmutscan).
# Needs ISA-L (headers and library; tested with 2.31.0 and 2.32.1) and POSIX threads:
#   make -f build.mk                          # ISA-L installed under /usr/local
#   make -f build.mk ISAL_PREFIX=/usr         # e.g. Debian/Ubuntu package libisal-dev
CC          ?= cc
CFLAGS      ?= -O3 -Wall -Wextra
ISAL_PREFIX ?= /usr/local

dms_data_lint: dms_data_lint.c
	$(CC) $(CFLAGS) -I$(ISAL_PREFIX)/include -o $@ $< -L$(ISAL_PREFIX)/lib -lisal -lpthread $(LDFLAGS)

clean:
	rm -f dms_data_lint
