DIRS = nmod_extra nmod_vec_extra nmod_poly_extra fmpz_extra nmod_mat_extra
	#  nmod_mat_poly_extra nmod_poly_mat_extra sagemath_extra

all:    clean
	$(foreach dir, $(DIRS), if [ -d $(dir) ]; then cd $(dir)/src ; make ; cd ../.. ; fi ;)

clean:
	rm -f lib/*
	rm -f include/*.h
	$(foreach dir, $(DIRS), if [ -d $(dir)/src ]; then cd $(dir)/src ; make clean ; cd ../.. ; fi ;)
	$(foreach dir, $(DIRS), if [ -d $(dir)/test ]; then cd $(dir)/test ; make clean ; cd ../.. ; fi ;)
	$(foreach dir, $(DIRS), if [ -d $(dir)/timings ]; then cd $(dir)/timings ; make clean ; cd ../.. ; fi ;)

.PHONY: doc
doc:
	doxygen doxygen.conf
