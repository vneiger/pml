


AUTOMAKE_OPTIONS = gnu

EXTRA_DIST=  perfpublisher.sh 

AM_CPPFLAGS = -I../../include 

LDADD =  ../../lib/libpml.la

PERFPUBLISHERFILE=tests-report.xml

TESTS = test_pluq 

test_pluq_SOURCES=test_pluq.c 


check_PROGRAMS = $(TESTS)

EXTRA_PROGRAMS = 


# Perfpublisher script interaction - AC 2016/06/09
perfpublisher:
	+chmod 777 ./perfpublisher.sh
	+./perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(TESTS)" "$(CXX)"


