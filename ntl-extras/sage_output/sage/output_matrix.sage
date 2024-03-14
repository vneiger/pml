def output_matrix(mat, fname):
    f = open(fname, 'w')
    f.write(str(mat.nrows()))
    f.write('\n')
    f.write(str(mat.ncols()))
    for i in range(mat.nrows()):
        for j in range(mat.ncols()):
            f.write('\n')
            t = mat[i][j].coefficients(false)
            for q in t:
                f.write(str(q))
                f.write(' ');
    f.write('\n')
    f.close()
