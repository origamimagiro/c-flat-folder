#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <time.h>

#include "io.h"
#include "svg.h"
#include "utility.h"
#include "convert.h"

struct FoldData {
    char *title;
    double eps, feps;
    int Vn, En, Fn, Pn, Sn, Cn;
    double time;
};

double lap(clock_t *t, int v) {
    clock_t t0 = clock();
    double out = ((double) t0 - *t)*1000/CLOCKS_PER_SEC;
    if (v) { printf("Time: %lf\n", out); }
    *t = t0;
    return out;
}

void process_file(
    char *path, char *file, char *ext,
    char *fold_out, char *svg_out,
    struct FoldData *data
) {
    printf("Processing: %s\n", path);

    int Vn = 0, En = 0, Fn = 0;
    double eps = 1;

    double (*V)[2] = NULL;
    int (*EV)[2] = NULL, **FV = NULL, *FVn = NULL;
    char *EA = NULL;

    if (strcmp(ext, "cp") == 0) {
        int Ln = 0;
        char *LA = NULL;
        double (*L)[2][2] = NULL;
        IO_CP_import(path, &Ln, &L, &LA);
        assert(Ln > 0, "Unable to parse CP file");

        int **EL = NULL, *ELn = NULL;
        eps = X_L_2_V_EV_EL(Ln, L, &Vn, &V, &En, &EV, &EL, &ELn);

        int **VV = NULL, *VVn = NULL;
        X_V_EV_2_VV_FV(Vn, V, En, EV, &VV, &VVn, &Fn, &FV, &FVn);
        free(L); free(VV); free(VVn);

        EA = malloc(En*sizeof(char));
        X_LA_EL_2_EA(LA, En, EL, ELn, EA);
        free(LA); free(EL); free(ELn);

    } else if (strcmp(ext, "fold") == 0) {
        IO_FOLD_import(path, &Vn, &V, &En, &EV, &EA, &Fn, &FV, &FVn);
        assert(Vn > 0, "Unable to parse FOLD file");

        int Vn_ = 0, En_ = 0, (*EV_)[2] = NULL, **EL = NULL, *ELn = NULL;
        double (*V_)[2] = NULL, (*L)[2][2] = malloc(En*2*2*sizeof(double));
        X_V_EV_2_L(V, En, EV, L);
        eps = X_L_2_V_EV_EL(En, L, &Vn_, &V_, &En_, &EV_, &EL, &ELn);

        int has_extra_vertices = (Vn != Vn_);
        if (has_extra_vertices) {
            free(V); free(EV); free(FV); free(FVn);
            Vn = Vn_; En = En_; V = V_; EV = EV_;
            V_ = NULL; EV_ = NULL; FV = NULL; FVn = NULL;
            char *EA_ = EA;

            int **VV = NULL, *VVn = NULL;
            X_V_EV_2_VV_FV(Vn, V, En, EV, &VV, &VVn, &Fn, &FV, &FVn);
            free(VV); free(VVn);

            EA = malloc(En*sizeof(char));
            X_LA_EL_2_EA(EA_, En, EL, ELn, EA);
            free(EA_);
        } else if (Fn == 0) {
            int **VV = NULL, *VVn = NULL;
            X_V_EV_2_VV_FV(Vn, V, En, EV, &VV, &VVn, &Fn, &FV, &FVn);
            free(VV); free(VVn);
        }

        free(L); free(V_); free(EV_); free(EL); free(ELn);

    } else {
        assert(ext != NULL, "File extension unrecognized (expects CP or FOLD)");
    }

    if (fold_out != NULL) {
        char filename[256] = {0};
        sprintf(filename, "./%s%s", fold_out, file);
        int off = 1 + strlen(fold_out) + strlen(file) - strlen(ext);
        sprintf(filename + off, ".fold");
        printf(" - Writing FOLD: %s\n", filename);
        IO_FOLD_export(filename, file, Vn, V, En, EV, EA, Fn, FV, FVn);
    }

    int Pn = 0, Sn = 0, Cn = 0;
    double cell_eps;

    int *Ff = NULL, (*SP)[2] = NULL, **CP = NULL, *CPn = NULL;
    double (*Vf)[2] = NULL, (*P)[2] = NULL;

    {
        X_V_EV_EA_FV_2_Vf_Ff(Vn, V, En, EV, EA, Fn, FV, FVn, &Vf, &Ff);

        double (*L)[2][2] = malloc(En*2*2*sizeof(double));
        X_V_EV_2_L(Vf, En, EV, L);

        int **SE = NULL, *SEn = NULL;
        cell_eps = X_L_2_V_EV_EL(En, L, &Pn, &P, &Sn, &SP, &SE, &SEn);
        free(SE); free(SEn); free(L);

        int **PP = NULL, *PPn = NULL;
        X_V_EV_2_VV_FV(Pn, P, Sn, SP, &PP, &PPn, &Cn, &CP, &CPn);
        free(PP); free(PPn);
    }

    if (svg_out != NULL) {
        char filename[256] = {0};
        sprintf(filename, "./%s%s", svg_out, file);

        int off = 1 + strlen(svg_out) + strlen(file) - strlen(ext);

        sprintf(filename + off, ".svg");
        printf(" - Writing SVG: %s\n", filename);
        SVG_V_EV_EA_FV_2_path(Vn, V, En, EV, EA, Fn, FV, FVn, eps, filename);

        sprintf(filename + off, "_folded.svg");
        printf(" - Writing SVG: %s\n", filename);
        SVG_V_EV_EA_FV_2_path(Vn, Vf, En, EV, EA, Fn, FV, FVn, eps, filename);

        sprintf(filename + off, "_overlap.svg");
        printf(" - Writing SVG: %s\n", filename);
        SVG_V_EV_EA_FV_2_path(Pn, P, Sn, SP, NULL, Cn, CP, CPn, cell_eps, filename);
    }

    free(V); free(EV); free(EA); free(FV); free(FVn);
    free(Vf); free(Ff); free(P); free(SP); free(CP); free(CPn);

    *data = (struct FoldData) {
        .title = file, .eps = eps, .feps = cell_eps,
        .Vn = Vn, .En = En, .Fn = Fn, .Pn = Pn, .Sn = Sn, .Cn = Cn
    };
}

void get_ext(char *path, char *ext) {
    ext[0] = 0;
    char *ext_p = NULL;
    for (char *c = path; *c != 0; ++c) {
        if (*c == '.') { ext_p = c + 1; }
    }
    if (ext_p == NULL) { return; }
    strncpy(ext, ext_p, 256);
    for (char *c = ext; *c != 0; ++c) { // make lowercase
        *c |= ((*c >= 'A') && (*c <= 'Z')) ? 0x60 : 0x0;
    }
}

void add_data(FILE *file, struct FoldData *data) {
    if (file == NULL) { return; }
    fprintf(file, "%s,%i,%i,%i,%.17lf,%i,%i,%i,%.17lf,%lf",
        data->title, data->Vn, data->En, data->Fn, data->eps,
        data->Pn, data->Sn, data->Cn, data->feps, data->time
    );
    fprintf(file, "\n");
    printf("    - Found cp epsilon: %.17g\n", data->eps);
    printf("    - Found %i Vertices\n", data->Vn);
    printf("    - Found %i Edges\n",    data->En);
    printf("    - Found %i Faces\n",    data->Fn);
    printf("    - Found cell eps: %.17g\n", data->feps);
    printf("    - Found %i Points\n",   data->Pn);
    printf("    - Found %i Segments\n", data->Sn);
    printf("    - Found %i Cells\n",    data->Cn);
    printf("    - Time: %lf\n",    data->time);
}

int main(int argc, char **argv) {
    clock_t t0 = clock();
    assert(argc >= 2, "Requires at least 1 command line argument: a path");

    FILE *v_out    = NULL;
    char *svg_out  = NULL;
    char *fold_out = NULL;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-v"   ) == 0) { v_out = fopen(argv[++i], "w");  continue; }
        if (strcmp(argv[i], "-svg" ) == 0) { svg_out = argv[++i];   continue; }
        if (strcmp(argv[i], "-fold") == 0) { fold_out = argv[++i];  continue; }
        assert(0, "Found command line option that is not -svg or -fold");
    }

    char *path = argv[1];

    char ext[256] = {0};
    get_ext(path, ext);

    struct FoldData data = {};
    if (v_out != NULL) {
        fprintf(v_out, "title,vertices,edges,faces,eps,points,segments,cells,feps,time\n");
    }

    clock_t t1 = clock();
    if (ext[0] == 0) {
        printf("No extension found, processing path as folder\n");

        struct dirent **dr;
        int n = scandir(path, &dr, NULL, alphasort);

        assert(n != -1, "Could not find folder at input path");

        for (int i = 0; i < n; ++i) {
            struct dirent *de = dr[i];
            if (de == NULL) { continue; }
            char *title = de->d_name;
            get_ext(title, ext);
            if ((ext[0] != 0) &&
                ((strcmp(ext, "cp") == 0) || (strcmp(ext, "fold") == 0))
            ) {
                char filepath[256];
                sprintf(filepath, "%s%s", path, title);
                process_file(filepath, title, ext, fold_out, svg_out, &data);
                data.time = lap(&t1, 0);
                add_data(v_out, &data);
            }
            free(de);
        }
        free(dr);

    } else {
        printf("%s Extension found, processing path as file\n", ext);
        assert((strcmp(ext, "cp") == 0) || (strcmp(ext, "fold") == 0),
            "Input file must have extension CP or FOLD"
        );

        char *title;
        for (char *c = path; *c != 0; ++c) {
            if (*c == '/') { title = c + 1; }
        }

        process_file(path, title, ext, fold_out, svg_out, &data);
        data.time = lap(&t1, 0);
        add_data(v_out, &data);
    }

    if (v_out != NULL) { fclose(v_out); }

    lap(&t0, 1);
    return 0;
}
