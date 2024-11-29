CC=gcc
CFLAGS = -fPIC -Wall -O2 -march=native
#CFLAGS = -fPIC -Wall -O0
LIBS = -lm
LD_FLAGS = -shared -fPIC

SRC =  src/cov_eigen_v.c \
	src/eigen_v.c \
	src/covariance.c \
	src/cov_eigen_v_mesh.c \
	src/kdtree.c \
	src/compute_medial_axis.c \
	src/dc_patches.c

HDR = 	src/ft_constants.h src/faulttools.h src/SortPCtx.h src/kdtree.h

INC = src

# Replace .c with .o
_OBJ = $(SRC:.c=.o)

# Replace src/ with obj/
OBJ = $(patsubst src/%, obj/%, $(_OBJ))

VT_LIBNAME = faulttools.so
VT_LIB = lib/$(VT_LIBNAME)

all: $(OBJ)
	@$(MAKE) $(VT_LIB)

lib/%.so: $(OBJ)
	@/bin/echo -n "[creating shared object] "
	@mkdir -p $(@D)
	$(CC) $(LD_FLAGS) -o $(VT_LIB) $(OBJ) $(LIBS)

build_lib:
	@/bin/echo -n "[creating shared object] "
	$(CC) $(LD_FLAGS) -o lib/$(VT_LIBNAME) $(OBJ) $(LIBS)

$(OBJ): obj/%.o : src/%.c $(HDR)
	@/bin/echo -n "[compiling source] "
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@ -I./$(INC)

clean:
	@rm -f lib/$(VT_LIBNAME)
	@rm -f $(OBJ)
