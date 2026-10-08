# FAST

FAST is the Finite element AnalisyS Tool. See
`FASTv2.4.4_UsersManual.pdf` for model definitions and analysis algorithms.

## Build and run

The default build needs a C++11 or newer compiler and GNU Make. Eigen is bundled.
It uses FAST's existing solvers and has no PETSc dependency:

```sh
make -j4
./fast /path/to/model           # reads model.dat, writes model.pos
make test                      # requires Python 3
```

Copy example `.dat` files from `dat/` into a working directory before running
them. FAST writes results next to its input. `-silent` suppresses progress output.

## Optional PETSc linear solvers

Install a **real, double-precision PETSc** development package, for example
`libpetsc-real-dev` on Debian/Ubuntu, or build PETSc from
[the official installation instructions](https://petsc.org/release/install/).
PETSc is an external, optional dependency; no PETSc source or binaries need to
be added to this repository. This integration was tested with PETSc 3.23.7.

For a source installation, expose its pkg-config file:

```sh
export PKG_CONFIG_PATH="$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig:$PKG_CONFIG_PATH"
make -j4 USE_PETSC=1
make USE_PETSC=1 test
```

Use the compiler compatible with your PETSc installation, e.g.
`make USE_PETSC=1 CC=mpicxx` for an MPI build. If pkg-config is unavailable,
provide `PETSC_CFLAGS` and `PETSC_LIBS` explicitly. The pkg-config library directory is embedded as a
runtime search path; for manual flags, set `PETSC_LIBDIR` if the shared library
is outside the system search path. Native and PETSc builds use
separate `obj/` and `obj-petsc/` directories; running `make` again restores the
native executable without cleaning. Header edits automatically trigger rebuilds.

### Choose the backend in a model

Existing inputs keep their original solver and matrix storage. To opt in, add:

```text
%HEADER.ANALYSIS.LINEAR.SOLVER
'petsc'

%HEADER.ANALYSIS.PETSC.OPTIONS
'-ksp_type gmres -pc_type ilu -ksp_rtol 1e-10 -ksp_max_it 1000'
```

Use `'native'` or omit `LINEAR.SOLVER` for FAST's original solvers. This selector
is independent of `%HEADER.ANALYSIS.MATRIX.TYPE`: all existing matrix formats
work with PETSc, including analysis controllers that allocate skyline matrices
directly. Without `USE_PETSC=1`, requesting PETSc produces a clear error.
Options strings use single quotes, like other FAST labels, and can span lines.

FAST passes PETSc options through to `KSPSetFromOptions`; there is no FAST-side
list of allowed methods. Any KSP method and PC supported by your installed PETSc
can be selected, subject to their matrix and algorithm requirements. Examples:

| Method | Options |
| --- | --- |
| GMRES with incomplete LU | `-ksp_type gmres -pc_type ilu` |
| Conjugate gradient for positive definite systems | `-ksp_type cg -pc_type jacobi` |
| BiCGStab | `-ksp_type bcgs -pc_type jacobi` |
| Direct LU | `-ksp_type preonly -pc_type lu` |
| MUMPS direct solve, if PETSc was built with MUMPS | `-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps` |

GMRES/ILU is the default. FAST's matrix tolerance and maximum iteration settings
supply defaults where the analysis sets them; PETSc's `-ksp_rtol`, `-ksp_atol`,
and `-ksp_max_it` override them. For nonlinear and eigenvalue analyses, explicit
PETSc tolerances keep the inner linear solve independent of the outer iteration.

### Reusable options files and command-line overrides

Keep solver tuning in a separate file:

```text
%HEADER.ANALYSIS.LINEAR.SOLVER
'petsc'

%HEADER.ANALYSIS.PETSC.OPTIONS.FILE
'solver.opts'
```

`solver.opts` is resolved relative to the input model's directory and contains
ordinary PETSc options (no surrounding FAST quotes):

```text
-ksp_type gmres
-pc_type ilu
-ksp_rtol 1e-10
-ksp_converged_reason
```

Use inline options to override that file, or command-line options to override
both. PETSc's standard `PETSC_OPTIONS` environment variable is also supported;
precedence is environment, model options file, inline options, then command line.
For example:

```sh
./fast /path/to/model -ksp_type preonly -pc_type lu -silent
./fast /path/to/model -options_file /path/to/solver.opts -ksp_view
```

### Behavior and scope

PETSc replaces the equation solves, including repeated right-hand sides, MPC
linear solves, dynamic analyses, and inner solves of FAST's existing eigenvalue
algorithms. Nonconvergence, unavailable solvers, missing options files, and PETSc
errors stop the run with a nonzero exit status instead of printing completed
results. Native solver behavior is retained.

Assembly remains serial. Run one process; distributed assembly and SLEPc
eigenvalue algorithms are outside this integration. PETSc uses a sequential AIJ
matrix populated from FAST's stored entries in linear time, expanding symmetric
halves without a dense intermediate. Each solve rebuilds the PETSc matrix and
KSP to honor reassembly; factorization caching is a future optimization. FAST's
critical-point stability diagnostics still obtain Crout pivots by factoring a
separate native copy, leaving the PETSc matrix intact.

`make test` checks every original storage/solver choice and reference results
for frame, continuum, and thermal models. With `USE_PETSC=1`, it also checks
multiple KSP/PC choices, nonsymmetric matrices, reassembly, aliased right-hand
sides, stability pivots, options precedence, dynamic/nonlinear regressions, and
configuration/nonconvergence errors. Generated test files use temporary
directories and do not overwrite the example data.
