#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "point.h"
#include "utility.h"

void SVG_add_svg_open(
    FILE *file,
    double  x, double  y, double  w, double  h,
    double vx, double vy, double vw, double vh
) {
    fprintf(file, "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        "x=\"%.2lf\" y=\"%.2lf\" width=\"%.2lf\" height=\"%.2lf\" "
        "viewBox=\"%.2lf %.2lf %.2lf %.2lf\">\n", x, y, w, h, vx, vy, vw, vh);
}

void SVG_add_svg_close(
    FILE *file
) {
    fprintf(file, "</svg>\n");
}

void SVG_add_line(
    FILE *file,
    double x1, double y1,
    double x2, double y2,
    double width, char *color
) {
    fprintf(file, "<line x1=\"%.2lf\" y1=\"%.2lf\" x2=\"%.2lf\" y2=\"%.2lf\" "
        "stroke=\"%s\" stroke-width=\"%.2lf\" stroke-line-cap=\"round\"></line>\n", 
        x1, y1, x2, y2, color, width
    );
}

void SVG_add_polygon(
    FILE *file,
    double (*P)[2], int Pn,
    char *color
) {
    fprintf(file, "<polygon points=\"");
    for (int i = 0; i < Pn; ++i) {
        if (i > 0) { fprintf(file, " "); }
        fprintf(file, "%.2lf,%.2lf", P[i][0], P[i][1]);
    }
    fprintf(file, "\" fill=\"%s\"></polygon>\n", color);
}

void SVG_add_disk(
    FILE *file,
    double cx, double cy, double r,
    char *fill, char *stroke
) {
    fprintf(file, "<circle cx=\"%.2lf\" cy=\"%.2lf\" r=\"%.2lf\" "
        "fill=\"%s\" stroke=\"%s\"></circle>\n", cx, cy, r, fill, stroke);
}

void SVG_add_text(
    FILE *file,
    double x, double y,
    char *color, char *text
) {
    fprintf(file, "<text x=\"%.2lf\" y=\"%.2lf\" "
                  "fill=\"%s\">%s</text>\n", x, y, color, text);
}

void SVG_V_EV_EA_FV_2_path(
    int Vn, double (*V_)[2],
    int En, int (*EV)[2], char *EA,
    int Fn, int **FV, int *FVn,
    double eps, char *path
) {
    FILE *file = fopen(path, "w");
    assert(file != NULL, "Unable to open file at path");

    struct Point *V = (struct Point *) V_;

    struct Point min = {{{ INFINITY,  INFINITY}}};
    struct Point max = {{{-INFINITY, -INFINITY}}};
    for (int i = 0; i < Vn; ++i) {
        struct Point v = V[i];
        for (int j = 0; j < 2; ++j) {
            min.C[j] = (v.C[j] < min.C[j]) ? v.C[j] : min.C[j];
            max.C[j] = (v.C[j] > max.C[j]) ? v.C[j] : max.C[j];
        }
    }
    struct Point w = P_sub(max, min);
    double s = (w.x < w.y) ? w.y : w.x;
    double t = 1000;
    double b = 100;
    SVG_add_svg_open(file,
        -t/2, -t/2, t, t,
        -(b + t)/2, -(b + t)/2, b + t, b + t
    );
    for (int i = 0; i < Fn; ++i) {
        int fn = FVn[i];
        int *fV = FV[i];
        struct Point c = P_scale(P_polygon_centroid(fV, fn, V), t, s, min, w);
        struct Point *P = malloc(fn*sizeof(struct Point));
        for (int j = 0; j < fn; ++j) {
            struct Point v = V[fV[j]];
            P[j] = P_scale(v, t, s, min, w);
            P[j] = P_add(P_mul(P_sub(P[j], c), 0.5), c);
        }
        SVG_add_polygon(file, (double (*)[2]) P, fn, "#DDD");
        free(P);
    }
    char C[][5] = {"#000", "#F00", "#00F", "#AAA", "#FF0"};
    for (int i = 0; i < En; ++i) {
        struct Point p[2];
        for (int j = 0; j < 2; ++j) {
            struct Point v = V[EV[i][j]];
            p[j] = P_scale(v, t, s, min, w);
        }
        char *c = "#000";
        if (EA != NULL) {
            switch (EA[i]) {
                case 'B': c = C[0]; break;
                case 'M': c = C[1]; break;
                case 'V': c = C[2]; break;
                case 'F': c = C[3]; break;
                case 'U': c = C[4]; break;
            }
        }
        SVG_add_line(file, p[0].x, p[0].y, p[1].x, p[1].y, 1, c);
    }
    char buff[256];
    for (int i = 0; i < Fn; ++i) {
        int fn = FVn[i];
        int *fV = FV[i];
        struct Point c = P_scale(P_polygon_centroid(fV, fn, V), t, s, min, w);
        SVG_add_disk(file, c.x, c.y, 1, "#000", "none");
        sprintf(buff, "%i", i);
        c.x += 3;
        c.y -= 3;
        SVG_add_text(file, c.x, c.y, "#000", buff);
    }
    for (int i = 0; i < Vn; ++i) {
        struct Point v = V[i];
        struct Point c = P_scale(v, t, s, min, w);
        SVG_add_disk(file, c.x, c.y, eps*t/s, "none", "#070");
        c.x += 3;
        c.y -= 3;
        sprintf(buff, "%i", i);
        SVG_add_text(file, c.x, c.y, "#070", buff);
    }
    SVG_add_svg_close(file);
    fclose(file);
}
