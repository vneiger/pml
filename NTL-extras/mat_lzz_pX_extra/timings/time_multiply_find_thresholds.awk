BEGIN {
    print "#ifndef __THRESHOLDS_WASKMAN_EVALUATE__H";
    print "#define __THRESHOLDS_WASKMAN_EVALUATE__H";
}

{
    print $0;
}

END {
    print "#endif";
}
