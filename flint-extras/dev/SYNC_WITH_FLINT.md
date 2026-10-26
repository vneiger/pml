# Parts of the build system shared with FLINT

PML's build system (`configure.ac`, `acinclude.m4`, `Makefile.in`, `config/`)
is adapted from FLINT's. Some parts are kept as **verbatim copies** of FLINT's,
so that PML detects the host and chooses its compiler flags (in particular
`-march`) exactly as FLINT does: PML must be compiled for the same instruction
set as FLINT (see the section "instruction set consistency with FLINT" of
`configure.ac`).

Last synced with FLINT commit `aa406e15`.

## Verbatim copies

| PML file | part | FLINT counterpart |
|---|---|---|
| `config/config.guess`, `config.sub`, `configfsf.guess`, `configfsf.sub` | whole files | same paths |
| `acinclude.m4` | the host patterns `define(..._PATTERN, ...)` | same |
| `configure.ac` | body of the section "architecture specifics" (apart from the issue URL) | same section |
| `configure.ac` | body of the section "CFLAGS" | same section |

In the section "architecture specifics", the variables `asm_path` and
`param_path` are unused by PML; they are kept so that the copy stays verbatim.

(`config/install-sh` and `config/ltmain.sh` are not tracked, neither in PML nor
in FLINT: `bootstrap.sh` generates them.)

## Checking and updating

From `flint-extras/`, with a FLINT checkout in `FLINT_DIR`:

    dev/check_flint_sync.sh FLINT_DIR       # one line per part: ok / DRIFT
    dev/check_flint_sync.sh -v FLINT_DIR    # also show the diffs

Run it before a release, and after changes to FLINT's configure. To update:

- copy the drifted files of `config/` from FLINT;
- replace the drifted pattern definitions of `acinclude.m4`, or the body of
  the drifted section of `configure.ac`, by FLINT's (in "architecture
  specifics", keep PML's issue URL);
- check whether the change needs something on PML's side, for example a new
  `have_avx2` / `have_avx512` setting or a new configure option;
- update the FLINT commit above, and run `./bootstrap.sh && ./configure`.

## Parts adapted from FLINT, not tracked

The rest of `configure.ac` and `acinclude.m4` (options, checks of GMP, MPFR,
FLINT, shared library names, substitutions) and `Makefile.in` are adapted to
PML and differ from FLINT's on purpose. FLINT options without effect on PML
(GC, NTL, BLAS, assembly, reentrant, thread-safe, GMP internals) were removed.
