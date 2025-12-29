void IO_CP_import(char *path, int *n, double (**L)[2][2], char **LA);

void IO_FOLD_import(char *path,
    int *Vn, double (**V)[2],
    int *En, int (**EV)[2], char **EA,
    int *Fn, int ***FV, int **FVn
);

void IO_FOLD_export(char *path, char *title,
    int Vn, double (*V)[2],
    int En, int (*EV)[2], char *EA,
    int Fn, int **FV, int *FVn
);
