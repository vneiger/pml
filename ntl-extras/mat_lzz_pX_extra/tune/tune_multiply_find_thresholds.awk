BEGIN {
    print "#ifndef __THRESHOLDS_WASKMAN_EVALUATE__H";
    print "#define __THRESHOLDS_WASKMAN_EVALUATE__H";
}

{
    print $0;
}

END {
    print "#endif";

    print "// Local Variables:";
    print "// mode: C++";
    print "// tab-width: 4";
    print "// indent-tabs-mode: nil";
    print "// c-basic-offset: 4";
    print "// End:";
    print "// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\\:0,t0,+0,=s";
}
