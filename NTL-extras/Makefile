DIRS = util sage_output magma_output lzz_p_extra vec_lzz_p_extra lzz_pX_extra lzz_pX_middle_product lzz_pX_CRT mat_lzz_p_extra structured_lzz_p mat_lzz_pX_extra lzz_pXY structured_lzz_pX lzz_pX_seq

all:	clean
	$(foreach dir, $(DIRS), if [ -d $(dir) ]; then cd $(dir)/src ; make ; cd ../.. ; fi ;)

clean:
	rm -f lib/*
	rm -f include/*.h
	$(foreach dir, $(DIRS), if [ -d $(dir)/src ]; then cd $(dir)/src ; make clean ; cd ../.. ; fi ;)
	$(foreach dir, $(DIRS), if [ -d $(dir)/test ]; then cd $(dir)/test ; make clean ; cd ../.. ; fi ;)
	$(foreach dir, $(DIRS), if [ -d $(dir)/timings ]; then cd $(dir)/timings ; make clean ; cd ../.. ; fi ;)
	$(foreach dir, $(DIRS), if [ -d $(dir)/tune ]; then cd $(dir)/tune ; make clean ; cd ../.. ; fi ;)

doc:
	cd include/ ; doxygen doxygen.conf ; cd ../
