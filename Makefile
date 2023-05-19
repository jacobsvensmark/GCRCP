FC=gfortran-7

src_path=for/

DEBUG_FLAGS=-g -O0 -fbacktrace -Wall -fcheck=all
CFLAGS=-O3
reversal_src=tji95_reversal.for
reversal_out=tji95_reversal.out
present_B_src=tji95_present_B.for
present_B_out=tji95_present_B.out

all: tji95_reversal.out  tji95_present_B.out

tji95_reversal.out: $(src_path)$(reversal_src)
	$(FC) $(CFLAGS) -o $(reversal_out) $(src_path)$(reversal_src)

tji95_present_B.out: $(src_path)$(present_B_src)
	$(FC) $(CFLAGS) -o $(present_B_out) $(src_path)$(present_B_src)

debug: debug_reversal debug_present_B

debug_reversal: tji95_reversal.out
	$(FC) $(DEBUG_FLAGS) $(CFLAGS) -o $(reversal_out)  $(src_path)$(reversal_src)
debug_present_B: tji95_present_B.out
	$(FC) $(DEBUG_FLAGS) $(CFLAGS) -o $(present_B_out) $(src_path)$(present_B_src)
