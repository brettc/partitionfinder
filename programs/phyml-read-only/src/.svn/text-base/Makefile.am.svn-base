# Uncomment (i.e. remove the `#' character at the begining of) 
# the two lines below if you want to use MPI.

# CC=mpicc
# DEFS=-DUNIX -D$(PROG) -DDEBUG -DMPI

# Comment the line below if you want to use MPI.
# Uncomment the line below if you don't want to use MPI.

DEFS = -DUNIX -D$(PROG) -DDEBUG $(REVISION)

# REVISION = $(shell if test -e .svn; then echo -DREV=\"`LC_ALL=C svn info | sed -n '/^Rev/p'| sed -e 's/^Revision:\ //'`\"; fi;)


if WANT_PHYTIME
bin_PROGRAMS = phytime
PROG = PHYTIME
else
if WANT_PHYCONT
bin_PROGRAMS = phycont
PROG = PHYCONT
else
if WANT_PART
bin_PROGRAMS = part
PROG = PART
else
if WANT_RWRAP
PROG = RWRAP
else
if WANT_TIPORDER
bin_PROGRAMS = tiporder
PROG = TIPORDER
else 
if WANT_M4
bin_PROGRAMS = m4
PROG = M4
else 
# if WANT_PHYML
bin_PROGRAMS = phyml
PROG = PHYML
# endif
endif
endif
endif
endif
endif
endif

if WANT_PHYTIME
phytime_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
spr.c spr.h\
times.c times.h\
m4.c m4.h\
draw.c draw.h\
rates.c rates.h\
mcmc.c mcmc.h\
stats.c stats.h\
mg.c mg.h\
tiporder.c tiporder.h
phytime_LDADD = -lm
else 
if WANT_PHYCONT
phycont_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
spr.c spr.h\
times.c times.h\
m4.c m4.h\
draw.c draw.h\
rates.c rates.h\
mcmc.c mcmc.h\
stats.c stats.h\
mg.c mg.h\
tiporder.c tiporder.h
# continuous.c continuous.h
phycont_LDADD = -lm
else 
# if WANT_RWRAP
# lib_LTLIBRARIES = librwrap.la
# librwrap_la_SOURCES = main.c \
# utilities.c utilities.h\
# optimiz.c optimiz.h\
# lk.c lk.h\
# bionj.c bionj.h\
# models.c models.h\
# free.c free.h\
# help.c help.h\
# simu.c simu.h\
# eigen.c eigen.h\
# pars.c pars.h\
# alrt.c alrt.h\
# interface.c interface.h\
# cl.c cl.h\
# mg.c mg.h\
# spr.c spr.h\
# m4.c m4.h\
# draw.c draw.h\
# mcmc.c mcmc.h\
# stats.c stats.h\
# times.c times.h\
# tiporder.c tiporder.h\
# rates.c rates.h\
# rwrapper.c rwrapper.h
# librwrap_la_LIBADD = -lm
# librwrap_la_LDFLAGS = -I/Library/Frameworks/R.framework/Resources/include -I/Library/Frameworks/R.framework/Resources/include/i386 -I/usr/local/include -flat_namespace -undefined suppress -shared -module
# librwrap_la_CFLAGS=-std=gnu99 -fPIC
# else 
if WANT_PART
part_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
mg.c mg.h\
spr.c spr.h\
m4.c m4.h\
draw.c draw.h\
stats.c stats.h\
tiporder.c tiporder.h
part_LDADD = -lm
else
if WANT_TIPORDER
tiporder_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
mg.c mg.h\
times.c times.h\
mcmc.c mcmc.h\
rates.c rates.h\
spr.c spr.h\
m4.c m4.h\
draw.c draw.h\
stats.c stats.h\
tiporder.c tiporder.h
tiporder_LDADD = -lm
else
if WANT_M4
m4_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
mg.c mg.h\
times.c times.h\
mcmc.c mcmc.h\
rates.c rates.h\
spr.c spr.h\
draw.c draw.h\
stats.c stats.h\
tiporder.c tiporder.h\
m4.c m4.h
m4_LDADD = -lm
else
# if WANT_PHYML
phyml_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
mpi_boot.c  mpi_boot.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h
phyml_LDADD = -lm
# endif
# endif
endif
endif
endif
endif
endif

all-am:	intro $(bin_PROGRAMS)$(EXEEXT)
	@echo ""
	@echo "Done."

intro:	
	@echo ""
	@echo ""
	@echo ":: Building $(PROG). Version $(VERSION) ::"
	@echo ""
	@echo ""

