############################################################################
# Makefile of FAST program.
#
# Created:                 Evandro Parente Junior
#                                                            
# Modified:  24-Apr-2015   Evandro Parente Junior
#            Added _DEBUG_ in DEFINE.
#
# Modified:  29-Apr-2015   Elias Saraiva Barroso
#            Added IGA folder.
#
# Modified:  13-Mar-2018   Elias Saraiva Barroso
#            Moved IGA shape and load files and rename iga folder to geom.
############################################################################

# Aplication constants

DIR = .
CTRL	= $(DIR)/ctrl
ELC	= $(DIR)/elc
GEOM	= $(DIR)/geom
LOAD	= $(DIR)/load
MAIN	= $(DIR)/main
MAT	= $(DIR)/mat
MLIB	= $(DIR)/mlib
NODE	= $(DIR)/node
OBJ	= $(DIR)/obj$(if $(filter 1,$(USE_PETSC)),-petsc,)
SEC	= $(DIR)/sec
SHP	= $(DIR)/shp
SPRING  = $(DIR)/spring
HEAT    = $(DIR)/heat

# Directory constants

DIRS =	$(CTRL) $(ELC) $(GEOM) $(LOAD) $(MAIN) $(MAT) $(MLIB) \
        $(NODE) $(OBJ) $(SHP) $(SEC)  $(SPRING)  $(HEAT)

# Compilation parameters

CC	= g++
POLICE	= -Wall
#DEBUG	= 
DEBUG	= -g 
#DEBUG	= -g -pg
#OPTMIZ	=
OPTMIZ	= -O3
#DEFINE	= -D_UNIX_
DEFINE	= -D_UNIX_ -D_DEBUG_
#DEFINE	= -D_UNIX_ -D_DEBUG_ -D_OMP_ 
#DEFINE	= -D_UNIX_ -D_DEBUG_ -D_OMP_ -D_EXPR_
INCLUDE	= -I$(CTRL) -I$(ELC) -I$(GEOM) -I$(LOAD) -I$(MAIN) -I$(MAT) -I$(MLIB) \
          -I$(NODE) -I$(SHP) -I$(SEC) -I$(SPRING) -I$(HEAT)
CFLAGS	= $(INCLUDE) $(POLICE) $(DEBUG) $(OPTMIZ) $(DEFINE)
#CFLAGS = -fopenmp $(INCLUDE) $(POLICE) $(DEBUG) $(OPTMIZ) $(DEFINE)
SYSLIBS	= -lm

# Optional PETSc KSP/PC backend; native solvers remain the default.
# PETSC_CFLAGS/PETSC_LIBS can be supplied for installations without pkg-config.
ifeq ($(USE_PETSC),1)
PETSC_CFLAGS ?= $(shell pkg-config --cflags PETSc 2>/dev/null)
PETSC_LIBS ?= $(shell pkg-config --libs PETSc 2>/dev/null)
PETSC_LIBDIR ?= $(shell pkg-config --variable=libdir PETSc 2>/dev/null)
ifeq ($(strip $(PETSC_LIBS)),)
$(error PETSc not found: set PKG_CONFIG_PATH or PETSC_CFLAGS and PETSC_LIBS)
endif
CFLAGS += -DFAST_WITH_PETSC $(PETSC_CFLAGS)
SYSLIBS += $(PETSC_LIBS)
ifneq ($(strip $(PETSC_LIBDIR)),)
SYSLIBS += -Wl,-rpath,$(PETSC_LIBDIR)
endif
endif


# Aplication modules

CTRLMOD	=		\
	ctrl		\
	ctrlarclen	\
	ctrldsp		\
	ctrlinc		\
	ctrlgenalpha	\
	ctrllins	\
	ctrllinstab	\
	ctrlmodal	\
	ctrlnewmark	\
	ctrlnlstab	\
	ctrlnr		\
	ctrlhht 	\
	ctrlpath	\
	ctrlqstc	\
	ctrlsec		\
	ctrlheat	\
	field

ELCMOD	=		\
	anmodel		\
	element		\
	elmfrm3d	\
	elmfrm3dcr	\
        elmgrid         \
	elminterf	\
	elmparam	\
	elmparamC1	\
	elmdsh		\
	elmplfrm	\
	elmpltrs	\
	elmtrs3d	\
	intpoint

GEOMOD	=		\
	bernbasis	\
	bspbasis	\
	knotvec		\
	cpdata		\
	ctrlpnt		\
	patch		\
	subpat   	

LOADMOD	=		\
	load		\
	loadiga		\
	prescdof	\
        timefunc

MAINMOD	=		\
	input		\
	main		\
	utl

MATMOD	=		\
	cmodel		\
	cmodel1d	\
	material

MLIBMOD	=		\
	linsolver	\
	mat		\
	matvec		\
	sysmat		\
	eig		\
	vec	

NODEMOD	=		\
	node

SECMOD	=		\
	section		\
	secbar  	\
	secanalysis	

SHPMOD	=		\
	shape		\
	shpinf		\
	shpinterf	\
	shpinterfiga	\
	shpline		\
	shpplane	\
	shpsolid	\
	shpsurf		\
	shpshell	\
	shpiga   	\
	shpbsp   	\
	shpbez   	

SPRINGMOD  =            \
        spring		\
        sprprop

HEATMOD  =		\
        heat		


# Object modules

OBJS	= $(CTRLMOD:%=$(OBJ)/%.o)	\
	  $(ELCMOD:%=$(OBJ)/%.o)	\
	  $(GEOMOD:%=$(OBJ)/%.o) 	\
	  $(LOADMOD:%=$(OBJ)/%.o)	\
	  $(MAINMOD:%=$(OBJ)/%.o)	\
	  $(MATMOD:%=$(OBJ)/%.o)	\
	  $(MLIBMOD:%=$(OBJ)/%.o)	\
	  $(NODEMOD:%=$(OBJ)/%.o)	\
	  $(SECMOD:%=$(OBJ)/%.o)	\
	  $(SHPMOD:%=$(OBJ)/%.o)	\
	  $(SPRINGMOD:%=$(OBJ)/%.o)	\
	  $(HEATMOD:%=$(OBJ)/%.o)

# Executable (To use the profile, use '-static -pg' )

FIRST: $(DIRS)
	$(MAKE) $(DIR)/fast

$(OBJ):
	mkdir -p $@

$(DIR)/fast     : $(DIRS) $(OBJS) FORCE
		$(CC) -o $@ $(OBJS) $(SYSLIBS)
#		$(CC) -fopenmp -o $@ $(OBJS) $(SYSLIBS)
#                 $(CC) -o $@ $(OBJS) -static -pg $(SYSLIBS)
#		$(CC) -fopenmp -o $@ $(OBJS) -static -pg $(SYSLIBS)

# Object compilation

$(OBJ)/%.o	: $(CTRL)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(ELC)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(GEOM)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(LOAD)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(MAIN)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(MAIN)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(MAT)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(MLIB)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(NODE)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(SEC)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(SHP)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(SPRING)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

$(OBJ)/%.o	: $(HEAT)/%.cpp
		  $(CC) -c $(CFLAGS) -o $@ $<

# Dependencies generation

depend:	FORCE
	@echo ""
	@echo "Generating CTRL dependencies..."
	-@g++ $(INCLUDE) -MM $(CTRL)/*.cpp > fastdep.tmp
	@echo "Generating ELC dependencies..."
	-@g++ $(INCLUDE) -MM $(ELC)/*.cpp >> fastdep.tmp
	@echo "Generating GEOM dependencies..."
	-@g++ $(INCLUDE) -MM $(GEOM)/*.cpp >> fastdep.tmp
	@echo "Generating LOAD dependencies..."
	-@g++ $(INCLUDE) -MM $(LOAD)/*.cpp >> fastdep.tmp
	@echo "Generating MAT dependencies..."
	-@g++ $(INCLUDE) -MM $(MAT)/*.cpp >> fastdep.tmp
	@echo "Generating MLIB dependencies..."
	-@g++ $(INCLUDE) -MM $(MLIB)/*.cpp >> fastdep.tmp
	@echo "Generating NODE dependencies..."
	-@g++ $(INCLUDE) -MM $(NODE)/*.cpp >> fastdep.tmp
	@echo "Generating SEC dependencies..."
	-@g++ $(INCLUDE) -MM $(SEC)/*.cpp >> fastdep.tmp
	@echo "Generating SHP dependencies..."
	-@g++ $(INCLUDE) -MM $(SHP)/*.cpp >> fastdep.tmp
	@echo "Generating SPRING dependencies..."
	-@g++ $(INCLUDE) -MM $(SPRING)/*.cpp >> fastdep.tmp
	@echo "Generating HEAT dependencies..."
	-@g++ $(INCLUDE) -MM $(HEAT)/*.cpp >> fastdep.tmp
	@sed -e '1,$$s/^\([^ ]\)/$$(OBJ)\/\1/' < fastdep.tmp > fastdep && \
          rm -f fastdep.tmp
	@echo "Generation completed."

FORCE:

# Existing model regressions and matrix-level solver contracts.
TESTOBJS = $(addprefix $(OBJ)/,linsolver.o sysmat.o vec.o mat.o matvec.o)
$(OBJ)/linear-solver-test: $(DIR)/tests/linear_solver.cpp $(TESTOBJS) $(DIR)/makefile
	$(CC) $(CFLAGS) -o $@ $< $(TESTOBJS) $(SYSLIBS)

.PHONY: test
test: $(DIR)/fast $(OBJ)/linear-solver-test
	$(OBJ)/linear-solver-test native
	python3 $(DIR)/tests/solvers.py $(DIR)/fast $(if $(filter 1,$(USE_PETSC)),--petsc,)
ifeq ($(USE_PETSC),1)
	$(OBJ)/linear-solver-test petsc -ksp_type gmres -pc_type ilu
	$(OBJ)/linear-solver-test petsc-spd -ksp_type cg -pc_type jacobi
	$(OBJ)/linear-solver-test petsc -ksp_type preonly -pc_type lu
endif

# Track header changes for both build variants.
$(OBJS): $(wildcard $(CTRL)/*.h $(ELC)/*.h $(GEOM)/*.h $(LOAD)/*.h $(MAIN)/*.h $(MAT)/*.h $(MLIB)/*.h $(NODE)/*.h $(SEC)/*.h $(SHP)/*.h $(SPRING)/*.h $(HEAT)/*.h)
$(OBJS): | $(OBJ)

# dependencies


###################################################### End of file #########
