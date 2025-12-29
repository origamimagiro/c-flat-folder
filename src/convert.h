void X_V_EV_2_L(double (*V)[2], int En, int (*EV)[2], double (*L)[2][2]);

void X_LA_EL_2_EA(char *LA, int En, int **EL, int *ELn, char *EA);

double X_L_2_V_EV_EL(
    const int Ln, const double (*L)[2][2],
    int *Vn, double (**V)[2],
    int *En, int (**EV)[2], int ***EL, int **ELn
);

void X_L_eps_2_V_EV_EL(
    const int Ln, const double (*L)[2][2], const double eps,
    int *Vn, double (**V)[2],
    int *En, int (**EV)[2], int ***EL, int **ELn
);

void X_V_EV_2_VV_FV(
    int Vn, double (*V)[2],
    int En, int (*EV)[2],
    int ***VV, int **VVn,
    int *Fn, int ***FV, int **FVn
);

void X_V_EV_EA_FV_2_Vf_Ff(
    int Vn, double (*V)[2],
    int En, int (*EV)[2], char *EA,
    int Fn, int **FV, int *FVn,
    double (**Vf)[2], int **Ff
);
