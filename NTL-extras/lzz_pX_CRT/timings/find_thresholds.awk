BEGIN {
    OK = 0;
    DONE = 0;
}

{
    if (OK == 1 && DONE == 0){
        printf("#define MIN_GEOMETRIC_FFT_%s %i\n", type, $1);
        printf("#define MAX_GEOMETRIC_FFT_%s %i\n", type, $2);
        DONE = 1;
    }
    if ($1 == "prime" && $2 == p)
        OK = 1;
}


END {
}
