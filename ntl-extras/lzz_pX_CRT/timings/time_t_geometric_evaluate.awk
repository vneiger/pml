BEGIN {
    ind = 0;
    dmin = 99999;
    dmax = 0;
}

{
    p = $1;
    d = $2;
    t_fft_direct = $3;
    t_fft_full = $4;
    t_direct_full = $5;
    t_fft_half = $6;
    t_direct_half = $7;

    if (d < dmin)
        dmin = d;
    if (d > dmax)
        dmax = d;
    
    index_p = 0;

    for (i = 1; i <= ind; i++)
        if (P[i] == p)
            index_p = i;

    if (index_p == 0) {
        ind++;
        index_p = ind;
        P[ind] = p;
    }

    FF[index_p][d] = t_fft_full;
    DF[index_p][d] = t_direct_full;
    FH[index_p][d] = t_fft_half;
    DH[index_p][d] = t_direct_half;
}

END {
    mult = 0.95;

    print "full";
    for (index_p in P){
        start = dmin;
        print "prime", P[index_p];
        for (d = dmin+1; d <= dmax; d++){
            if ((DF[index_p][d-1] <= mult*FF[index_p][d-1]) && (DF[index_p][d] > mult*FF[index_p][d])){
                print start, d;
            }
            if ((DF[index_p][d-1] > mult*FF[index_p][d-1]) && (DF[index_p][d] <= mult*FF[index_p][d])){
                start = d;
            }
        }
    }

    print "half";
    for (index_p in P){
        start = dmin;
        print "prime", P[index_p];
        for (d = dmin+1; d <= dmax; d++){
            if ((DH[index_p][d-1] <= mult*FH[index_p][d-1]) && (DH[index_p][d] > mult*FH[index_p][d])){
                print start, d-1;
            }
            if ((DH[index_p][d-1] > mult*FH[index_p][d-1]) && (DH[index_p][d] <= mult*FH[index_p][d])){
                start = d;
            }
        }
    }
}
