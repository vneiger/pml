BEGIN {
    print "#ifndef __THRESHOLDS_MP_NAIVE_EVALUATE__H";
    print "#define __THRESHOLDS_MP_NAIVE_EVALUATE__H";
}

{
    print $0;
}

END {
    print "#endif";
}
