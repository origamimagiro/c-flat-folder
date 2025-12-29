void SVG_add_svg_open(FILE *file, double  x, double  y, double  w, double  h,
                                  double vx, double vy, double vw, double vh);
void SVG_add_svg_close(FILE *file);
void SVG_add_line(FILE *file, double x1, double y1,
                              double x2, double y2, double width, char *color);
void SVG_add_polygon(FILE *file, double (*P)[2], int Pn, char *color);
void SVG_add_disk(FILE *file, double x, double y, double r, char *color);
void SVG_add_text(FILE *file, double x, double y, char *color, char *text);

void SVG_V_EV_EA_FV_2_path(int Vn, double (*V)[2],
    int En, int (*EV)[2], char *EA,
    int Fn, int **FV, int *FVn,
    double eps, char *path
);
