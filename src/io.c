#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "assert.h"
#include "da.h"

static
int seek_past(FILE *f, char *s) {
    char c;
    while ((c = fgetc(f)) != EOF) {
        for (char *a = s; *a != 0; ++a) {
            if (c == *a) { return c; }
        }
    }
    return 0;
}

void IO_CP_import(char *path, int *Ln, double **L, char **LA) {
    assert((*L == NULL) && (*LA == NULL), "Outputs should start NULL");
    *Ln = 0;

    FILE *file = fopen(path, "r");
    assert(file != NULL, "Unable to open file at path");

    struct DA L_ = {sizeof(double)}, LA_ = {sizeof(char)};
    {
        char T[] = "?BMVFU";
        unsigned a, a_lim = strlen(T);
        while (fscanf(file, "%i", &a) == 1) {
            if (a >= a_lim) { a = 0; }
            double c[4];
            for (int i = 0; i < 4; ++i) {
                if (fscanf(file, "%lf", &(c[i]))) { continue; }
                fclose(file);
                DA_empty(&L_);
                DA_empty(&LA_);
                return;
            }
            DA_push(&LA_, T + a);
            for (int i = 0; i < 4; ++i) { DA_push(&L_, c + i); }
        }
    }

    fclose(file);
    *Ln = LA_.n;
    *L = DA_freeze(&L_);
    *LA = DA_freeze(&LA_);
}

void IO_FOLD_import(char *path,
    int *Vn,  double (**V)[2],
    int *En,  int (**EV)[2], char **EA,
    int *Fn,  int ***FV, int **FVn
) {
    assert((*V == NULL) && (*EV == NULL) && (*EA == NULL) &&
        (*FV == NULL) && (*FVn == NULL), "Outputs should start NULL");
    int EAn = *Vn = *En = *Fn = 0;

    FILE *file = fopen(path, "r");
    if (file == NULL) { return; }
    struct { char *V, *EV, *EA, *FV; } S = {
        "vertices_coords", "edges_vertices",
        "edges_assignment", "faces_vertices",
    };
    char buff[32];
    while (seek_past(file, "\"")) {
        fscanf(file, "%31s", buff);

        if ((*Vn == 0) && !strncmp(buff, S.V, strlen(S.V))) {

            if (!seek_past(file, "[")) { continue; }
            struct DA P = {sizeof(double)};
            double x, y;
            while (fscanf(file, " [ %lf , %lf ] , ", &x, &y) == 2) {
                DA_push(&P, &x);
                DA_push(&P, &y);
            }
            *Vn = P.n/2;
            *V = (double (*)[2]) DA_freeze(&P);

        } else if ((*En == 0) && !strncmp(buff, S.EV, strlen(S.EV))) {

            if (!seek_past(file, "[")) { continue; }
            struct DA EV_ = {sizeof(int)};
            int u, v;
            while (fscanf(file, " [ %i , %i ] , ", &u, &v) == 2) {
                DA_push(&EV_, &u);
                DA_push(&EV_, &v);
            }
            assert((EAn == 0) || (2*EAn == EV_.n), "EV data not aligned with EA");
            *En = EV_.n/2;
            *EV = (int (*)[2]) DA_freeze(&EV_);

        } else if ((EAn == 0) && !strncmp(buff, S.EA, strlen(S.EA))) {

            if (!seek_past(file, "[")) { continue; }
            struct DA EA_ = {sizeof(char)};
            char c;
            while (fscanf(file, " \"%c\" , ", &c) == 1) {
                DA_push(&EA_, &c);
            }
            assert((*En == 0) || (*En == EA_.n), "EA data not aligned with EV");
            EAn = EA_.n;
            *EA = (char *) DA_freeze(&EA_);

        } else if ((*Fn == 0) && !strncmp(buff, S.FV, strlen(S.FV))) {

            if (!seek_past(file, "[")) { continue; }
            struct DA FV_ = {sizeof(struct DA)};
            int sum = 0;
            while (seek_past(file, "[]") != ']') {
                struct DA V = {sizeof(int)};
                int f;
                while (fscanf(file, " %i , ", &f) == 1) { DA_push(&V, &f); }
                sum += V.n;
                DA_push(&FV_, &V);
                seek_past(file, "]");
            }
            *Fn = FV_.n;
            *FVn = malloc((*Fn)*sizeof(int));
            *FV = malloc((*Fn)*sizeof(int*) + sum*sizeof(int));
            int *FVD = (int*) (*FV + *Fn);

            for (int i = 0; i < *Fn; ++i) {
                struct DA *fV = DA_getp(&FV_, i);
                int fn = (*FVn)[i] = fV->n;
                int *fV_ = (*FV)[i] = FVD;
                FVD += fn;
                
                for (int j = 0; j < fn; ++j) {
                    DA_get(fV, j, &(fV_[j]));
                }
                DA_empty(fV);
            }
            DA_empty(&FV_);

        }
    }
    assert(EAn == *En, "EA and EV data not aligned");
    fclose(file);
}

void IO_FOLD_export(char *path, char *title,
    int Vn,  double (*V)[2],
    int En,  int (*EV)[2], char *EA,
    int Fn,  int **FV, int *FVn
) {
    FILE *file = fopen(path, "w");
    assert(file != NULL, "Unable to open file at path");

    fprintf(file,
        "{\n"
        "  \"file_spec\": 1.2,\n"
        "  \"file_creator\": \"FFc\",\n"
        "  \"file_classes\": [\"singleModel\"],\n"
        "  \"file_title\": \"%s\",\n", title
    );

    fprintf(file, "  \"vertices_coords\": [\n");
    for (int i = 0; i < Vn; ++i) {
        fprintf(file, "    [%.17g, %.17g]%s\n",
            V[i][0], V[i][1], (i == Vn - 1) ? "" : ",");
    }
    fprintf(file, "  ],\n");

    fprintf(file, "  \"edges_vertices\": [\n");
    for (int i = 0; i < En; ++i) {
        fprintf(file, "    [%i, %i]%s\n",
            EV[i][0], EV[i][1], (i == En - 1) ? "" : ",");
    }
    fprintf(file, "  ],\n");

    fprintf(file, "  \"edges_assignment\": [\n");
    for (int i = 0; i < En; ++i) {
        fprintf(file, "    \"%c\"%s\n", EA[i], (i == En - 1) ? "" : ",");
    }
    fprintf(file, "  ],\n");

    fprintf(file, "  \"faces_vertices\": [\n");
    for (int i = 0; i < Fn; ++i) {
        fprintf(file, "    [");
        int fn = FVn[i];
        int *F = FV[i];
        for (int j = 0; j < fn; ++j) {
            if (j > 0) { fprintf(file, ", "); }
            fprintf(file, "%i", F[j]);
        }
        fprintf(file, "]%s\n", (i == Fn - 1) ? "" : ",");
    }
    fprintf(file, "  ]\n");

    fprintf(file, "}");
    fclose(file);
}
