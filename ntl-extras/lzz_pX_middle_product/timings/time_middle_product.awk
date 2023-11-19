BEGIN {
    print "# middle product timings";
    print "# prime | average ratio T_mid/T_naive | average ratio T_mid/T_direct";
    print "# (prime = 0 for an FFT prime)"
}

{
    p = $1;
    i = $2;
    j = $3;
    tmid = $4;
    tnaive = $5;
    tdirect = $6;

    MN[p] += tmid/tnaive;  
    MD[p] += tmid/tdirect; 
    m[p] += 1;
}

END {
    for (p in MN) {
	print p, MN[p]/m[p], MD[p]/m[p];
    }
}
