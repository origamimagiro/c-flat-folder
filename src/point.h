struct Point { union { double C[2]; struct { double x, y; }; }; };

struct Point P_sub(struct Point p, struct Point q);
struct Point P_add(struct Point p, struct Point q);
struct Point P_mul(struct Point p, double s);
struct Point P_div(struct Point p, double s);
double P_dot(struct Point p, struct Point q);
double P_magsq(struct Point p);
double P_mag(struct Point p);
struct Point P_unit(struct Point p);
struct Point P_perp(struct Point p);
double P_ang(struct Point p);
double P_distsq(struct Point p, struct Point q);
double P_dist(struct Point p, struct Point q);
double P_cross(struct Point p, struct Point q);

double P_scale_vals(struct Point *P, int Pn, struct Point *min, struct Point *diag);

struct Point P_scale(struct Point p, double t, double s, struct Point min, struct Point diag);

int P_point_comp(struct Point p1, struct Point p2, double eps);
int P_line_intersect(
    struct Point p1, struct Point p2,
    struct Point p3, struct Point p4, double eps, struct Point *out
);

double P_polygon_area2(int *P, int Pn, struct Point *V);
struct Point P_polygon_centroid(int *P, int Pn, struct Point *V);
