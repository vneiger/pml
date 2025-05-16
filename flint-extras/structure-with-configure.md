See structure.md 

Directory flint-extras

---

#### Installation  

1. run `./autogen` with usual configure options, e.g. `--with-gmp=/opt/homebrew --with-flint=/usr/local`  or `configure`, directly 
2. make 
3. make install 
4. make check

---

#### Directories, handled by configure 

`fmpz_extra`
`nmod_extra, nmod_mat_extra, nmod_vec_extra`

The tree structure of flint-extras is kept, with extra directories and files 

- .h in respective directories, symbolic links to `include/` with configure 

- .c in respective src directories (no makefile), objects built in `lib/` by top make 

- local test: flint-extras original directories unchanged, makefiles (modified for macOS also) should work 

- top directory tests: new, for building objects using make check  

- local directories tests: slightly modified .c files from original test directories, accessed via top make check 


---

#### Management 

- `configure.ac`:  modify `AC_CONFIG_LINKS` for include files 
- `include/Makefile.am`: add include files for installation  
- `lib/Makefile.am`: add src directories and files for building objects and the library 
- `tests/Makefile.am`: add tests directories and files for make check 

---

#### Directories, not yet handled by configure 

Should work as usually in flint-extras (possibly slighlty modified makefiles for macOS also)

---

#### 
