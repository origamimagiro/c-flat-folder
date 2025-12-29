#include <math.h>
#include <stdio.h>

struct Point {
    union {
        struct { double x, y; };
        double C[2];
    };
};

struct Point P_add(struct Point p, struct Point q) {
    return (struct Point) {.x = p.x + q.x, .y = p.y + q.y};
}

struct Point P_sub(struct Point p, struct Point q) {
    return (struct Point) {.x = p.x - q.x, .y = p.y - q.y};
}

struct Point P_mul(struct Point p, double s) {
    return (struct Point) {.x = (p.x)*s, .y = (p.y)*s};
}

struct Point P_div(struct Point p, double s) {
    return (struct Point) {.x = (p.x)/s, .y = (p.y)/s};
}

double P_dot(struct Point p, struct Point q) {
    return (p.x)*(q.x) + (p.y)*(q.y);
}

double P_magsq(struct Point p) {
    return P_dot(p, p);
}

double P_mag(struct Point p) {
    return sqrt(P_magsq(p));
}

struct Point P_unit(struct Point p) {
    return P_div(p, P_mag(p));
}

struct Point P_perp(struct Point p) {
    return (struct Point) {.x = p.y, .y = -p.x};
}

double P_ang(struct Point p) {
    double a = atan2(p.y, p.x);
    return a + ((a < 0) ? 2*M_PI : 0);
}

double P_distsq(struct Point p, struct Point q) {
    return P_magsq(P_sub(q, p));
}

double P_dist(struct Point p, struct Point q) {
    return sqrt(P_distsq(p, q));
}

double P_cross(struct Point p, struct Point q) {
    return (p.x)*(q.y) - (q.x)*(p.y);
}

double P_scale_vals(struct Point *P, int Pn,
    struct Point *min, struct Point *diag
) {
    struct Point *max = diag;
    min->x = min->y = INFINITY; 
    max->x = max->y = -INFINITY; 
    for (int i = 0; i < Pn; ++i) {
        struct Point *p = P + i;
        for (int j = 0; j < 2; ++j) {
            min->C[j] = (p->C[j] < min->C[j]) ? p->C[j] : min->C[j];
            max->C[j] = (p->C[j] > max->C[j]) ? p->C[j] : max->C[j];
        }
    }
    *max = P_sub(*max, *min);
    return diag->C[(diag->x < diag->y) ? 1 : 0];
}

struct Point P_scale(struct Point p, double t, double s,
    struct Point min, struct Point diag
) {
    return P_mul(P_sub(P_sub(p, min), P_div(diag, 2)), t/s);
}

int P_point_comp(struct Point p1, struct Point p2, double eps) {
    struct Point d = P_sub(p1, p2);
    return ((fabs(d.y) > eps) ? ((d.y < 0) ? -1 : 1) :
           ((fabs(d.x) > eps) ? ((d.x < 0) ? -1 : 1) : 0));
}

int P_line_intersect(
    struct Point p1, struct Point p2,
    struct Point p3, struct Point p4,
    double eps, struct Point *out
) {
    struct Point d12 = P_sub(p1, p2);
    struct Point d34 = P_sub(p3, p4);
    double denom = P_cross(d12, d34);
    if (denom < eps*eps) { return 0; }
    double c12 = P_cross(p1, p2);
    double c34 = P_cross(p3, p4);
    *out = P_div(P_sub(P_mul(d34, c12), P_mul(d12, c34)), denom);
    return 1;
}

double P_polygon_area2(int *P, int Pn, struct Point *V) {
    double area = 0; 
    struct Point p = V[P[Pn - 1]];
    for (int i = 0; i < Pn; ++i) {
        struct Point q = V[P[i]];
        area += (p.x + q.x)*(q.y - p.y);
        p = q;
    }
    return area;
}

struct Point P_polygon_centroid(int *P, int Pn, struct Point *V) {
    struct Point p = {0};
    for (int i = 0; i < Pn; ++i) {
        p = P_add(p, V[P[i]]);
    }
    return P_div(p, Pn);  
}
