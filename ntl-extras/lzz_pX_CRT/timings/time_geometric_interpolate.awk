BEGIN {
    ind = 0;
    dmin = 99999;
    dmax = 0;
}

{
    p = $1;
    d = $2;
    t_fft = $3;
    t_direct = $4;

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

    F[index_p][d] = t_fft;
    D[index_p][d] = t_direct;
}

END {
    mult = 0.95;

    print "full";
    for (index_p in P){
        start = dmin;
        print "prime", P[index_p];
        for (d = dmin+1; d <= dmax; d++){
            if ((D[index_p][d-1] <= mult*F[index_p][d-1]) && (D[index_p][d] > mult*F[index_p][d])){
                print start, d;
            }
            if ((D[index_p][d-1] > mult*F[index_p][d-1]) && (D[index_p][d] <= mult*F[index_p][d])){
                start = d;
            }
        }
    }
}
