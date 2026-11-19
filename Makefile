#=======================================================================
# PRIMS-F Makefile (Fortran 2018; audited with GNU Fortran 14.2)
# - Estrutura: ./scr (fontes), ./tmp (obj/mod), ./bin (exe)
# - Compilação paralela segura: dependências explícitas de módulos
# - Targets úteis:
#     make            -> release (bin/prims_f)
#     make check_all  -> debug + runtime checks (bin/prims_f_check_all)
#     make asan       -> AddressSanitizer (bin/prims_f_asan)
#=======================================================================
.DEFAULT_GOAL := all

# -------------------- diretórios --------------------
DIR_SRC := scr
DIR_OBJ := tmp/release
DIR_BIN := bin
DIR_OUT := output
DIR_DOC := doc

# -------------------- compilador --------------------
FC := gfortran

# Production code no longer passes internal procedures as callbacks.
# No compiler-specific executable/heap trampoline option is required.
CALLBACK_FLAGS ?=

# -------------------- flags base --------------------
FFLAGS_BASE := -std=f2018 -Wall -Wextra \
  -Wno-unused-parameter -Wno-unused-dummy-argument \
  -fimplicit-none -ffree-line-length-none \
  -J$(DIR_OBJ) -I$(DIR_OBJ) -I$(DIR_SRC) \
  -fopenmp $(CALLBACK_FLAGS) -Wa,--noexecstack

LDFLAGS_BASE := -fopenmp -Wl,-z,noexecstack
LIBS         ?= -llapack -lblas

# Release (padrão)
FFLAGS_RELEASE := $(FFLAGS_BASE) -O2
LDFLAGS_RELEASE := $(LDFLAGS_BASE)

# Debug pesado (linha exata + bounds + NaN traps)
FFLAGS_CHECKALL := $(FFLAGS_BASE) -O0 -g3 -fbacktrace -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-999 \
  -fno-omit-frame-pointer
LDFLAGS_CHECKALL := $(LDFLAGS_BASE)

# Sanitizers (ótimo p/ SIGSEGV silencioso)
FFLAGS_ASAN := $(FFLAGS_BASE) -O0 -g3 -fbacktrace -fcheck=all \
  -fsanitize=address,undefined -fno-omit-frame-pointer
LDFLAGS_ASAN := $(LDFLAGS_BASE) -fsanitize=address,undefined

# Permite sobrescrever pela linha de comando: make FFLAGS=...
FFLAGS ?= $(FFLAGS_RELEASE)
LDFLAGS ?= $(LDFLAGS_RELEASE)

# -------------------- alvo --------------------------
EXE := $(DIR_BIN)/prims_f

# -------------------- fontes/objs -------------------
# IMPORTANTE:
#   - NÃO use "gfortran scr/*.f90 ..." para o executável principal:
#       (1) a ordem dos módulos importa
#       (2) arquivos "program" de debug (ex.: debug_*.f90) causam múltiplos mains
#   - Aqui a ordem é explícita e o make sabe compilar com -j.
SRC := \
  $(DIR_SRC)/mod_kinds.f90 \
  $(DIR_SRC)/mod_rootfinding.f90 \
  $(DIR_SRC)/mod_void_helpers.f90 \
  $(DIR_SRC)/module_flow_regime.f90 \
  $(DIR_SRC)/module_config.f90 \
  $(DIR_SRC)/module_lapack_small.f90 \
  $(DIR_SRC)/module_chexal_1992.f90 \
  $(DIR_SRC)/module_closure_law.f90 \
  $(DIR_SRC)/module_find.f90 \
  $(DIR_SRC)/module_catenary.f90 \
  $(DIR_SRC)/module_steady_state.f90 \
  $(DIR_SRC)/module_jacobian.f90 \
  $(DIR_SRC)/module_matrix_roe.f90 \
  $(DIR_SRC)/module_roe_solver.f90 \
  $(DIR_SRC)/module_feed_contracts.f90 \
  $(DIR_SRC)/module_boundary_conditions.f90 \
  $(DIR_SRC)/module_solver_utils.f90 \
  $(DIR_SRC)/module_time_step.f90 \
  $(DIR_SRC)/module_active_set.f90 \
  $(DIR_SRC)/module_simulate_pipe_riser.f90 \
  $(DIR_SRC)/main.f90

OBJ := $(patsubst $(DIR_SRC)/%.f90,$(DIR_OBJ)/%.o,$(SRC))

.PHONY: all clean dirs check_all asan release debug help test

help:
	@echo "Targets:"
	@echo "  make / make release   -> bin/prims_f (release)"
	@echo "  make check_all        -> bin/prims_f_check_all (-fcheck=all, -ffpe-trap, etc.)"
	@echo "  make asan             -> bin/prims_f_asan (Address/UB sanitizers)"
	@echo "  make clean            -> remove tmp + bin executables"

release: all

test:
	FC='$(FC)' LIBS='$(LIBS)' CALLBACK_FLAGS='$(CALLBACK_FLAGS)' python3 scripts/run_unit_tests.py

all: $(EXE)

# -------------------- diretórios como alvos (seguro com -j) -----------
$(DIR_BIN) $(DIR_OBJ) $(DIR_OUT) $(DIR_DOC):
	@mkdir -p $@

dirs: $(DIR_BIN) $(DIR_OBJ) $(DIR_OUT) $(DIR_DOC)

# -------------------- link --------------------------
$(EXE): $(OBJ) | $(DIR_BIN)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

# -------------------- compile -----------------------
$(DIR_OBJ)/%.o: $(DIR_SRC)/%.f90 | $(DIR_OBJ)
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(DIR_OBJ) $(DIR_BIN)/prims_f $(DIR_BIN)/prims_f_check_all $(DIR_BIN)/prims_f_asan

# -------------------- builds alternativos ------------------------------
# Separate object directories avoid races between clean/all and mixing flags.
check_all:
	$(MAKE) all DIR_OBJ=tmp/check EXE=bin/prims_f_check_all \
	  FFLAGS='-std=f2018 -O0 -g3 -Wall -Wextra -fimplicit-none -ffree-line-length-none -fopenmp $(CALLBACK_FLAGS) -Wa,--noexecstack -Jtmp/check -Itmp/check -Iscr -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-999' \
	  LDFLAGS='$(LDFLAGS_CHECKALL)'

asan:
	$(MAKE) all DIR_OBJ=tmp/asan EXE=bin/prims_f_asan \
	  FFLAGS='-std=f2018 -O0 -g3 -Wall -Wextra -fimplicit-none -ffree-line-length-none -fopenmp $(CALLBACK_FLAGS) -Wa,--noexecstack -Jtmp/asan -Itmp/asan -Iscr -fcheck=all -fbacktrace -fsanitize=address,undefined -fno-omit-frame-pointer' \
	  LDFLAGS='$(LDFLAGS_ASAN)'

# -----------------------------------------------------------------------------
# Dependências explícitas (ordem segura para make -j com módulos)
# -----------------------------------------------------------------------------
$(DIR_OBJ)/mod_rootfinding.o:         $(DIR_OBJ)/mod_kinds.o
$(DIR_OBJ)/mod_void_helpers.o:        $(DIR_OBJ)/mod_kinds.o $(DIR_OBJ)/mod_rootfinding.o
$(DIR_OBJ)/module_config.o:           $(DIR_OBJ)/mod_kinds.o
$(DIR_OBJ)/module_lapack_small.o:     $(DIR_OBJ)/mod_kinds.o
$(DIR_OBJ)/module_flow_regime.o:      $(DIR_OBJ)/mod_void_helpers.o

# catenary usa mod_rootfinding
$(DIR_OBJ)/module_catenary.o:         $(DIR_OBJ)/mod_kinds.o $(DIR_OBJ)/mod_rootfinding.o

$(DIR_OBJ)/module_solver_utils.o:     $(DIR_OBJ)/mod_void_helpers.o
$(DIR_OBJ)/module_chexal_1992.o:      $(DIR_OBJ)/mod_kinds.o
$(DIR_OBJ)/module_closure_law.o:      $(DIR_OBJ)/mod_void_helpers.o $(DIR_OBJ)/module_chexal_1992.o
$(DIR_OBJ)/module_find.o:             $(DIR_OBJ)/module_closure_law.o $(DIR_OBJ)/mod_rootfinding.o $(DIR_OBJ)/module_flow_regime.o

$(DIR_OBJ)/module_steady_state.o:     $(DIR_OBJ)/module_find.o $(DIR_OBJ)/module_catenary.o $(DIR_OBJ)/module_solver_utils.o
$(DIR_OBJ)/module_jacobian.o:         $(DIR_OBJ)/module_closure_law.o $(DIR_OBJ)/module_find.o
$(DIR_OBJ)/module_matrix_roe.o:       $(DIR_OBJ)/module_find.o $(DIR_OBJ)/module_jacobian.o

# roe_solver precisa da matriz de Roe e lapack_small
$(DIR_OBJ)/module_roe_solver.o:       $(DIR_OBJ)/module_matrix_roe.o $(DIR_OBJ)/module_lapack_small.o

# boundary_conditions usa module_roe_solver
$(DIR_OBJ)/module_feed_contracts.o: $(DIR_OBJ)/module_config.o

$(DIR_OBJ)/module_boundary_conditions.o: $(DIR_OBJ)/module_feed_contracts.o \
                                      $(DIR_OBJ)/module_steady_state.o \
                                      $(DIR_OBJ)/module_roe_solver.o $(DIR_OBJ)/module_jacobian.o \
                                      $(DIR_OBJ)/module_find.o $(DIR_OBJ)/module_closure_law.o \
                                      $(DIR_OBJ)/module_config.o

$(DIR_OBJ)/module_time_step.o:        $(DIR_OBJ)/module_lapack_small.o

$(DIR_OBJ)/module_simulate_pipe_riser.o: \
                                      $(DIR_OBJ)/module_active_set.o \
                                      $(DIR_OBJ)/module_roe_solver.o $(DIR_OBJ)/module_boundary_conditions.o \
                                      $(DIR_OBJ)/module_time_step.o $(DIR_OBJ)/module_catenary.o $(DIR_OBJ)/module_solver_utils.o \
                                      $(DIR_OBJ)/module_flow_regime.o $(DIR_OBJ)/module_config.o

$(DIR_OBJ)/main.o:                    $(DIR_OBJ)/module_simulate_pipe_riser.o $(DIR_OBJ)/module_steady_state.o

$(DIR_OBJ)/module_active_set.o: $(DIR_OBJ)/module_time_step.o
