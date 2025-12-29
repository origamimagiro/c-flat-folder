#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "assert.h"
#include "svg.h"
#include "da.h"
#include "pq.h"
#include "compare.h"
#include "avl.h"
#include "point.h"
#include "hash.h"

#include "svg.h"

struct Segment {
    int P[2];
    struct Point U;
    double a, d;
};

static
double point_seg_dist(struct Point *p, struct Segment *s) {

    return s->d - P_dot(P_perp(s->U), *p);
}

static
int on_line(struct Point *p, struct Segment *s, double eps) {

    return fabs(point_seg_dist(p, s)) < eps;
}

int seg_comp(int i, int j, void *ctx) {

    struct Segment *S = ctx;
    int d = S[i].P[0] - S[j].P[0];

    return (d == 0) ? (S[i].P[1] - S[j].P[1]) : d;
}

struct Sweep_Context {
    struct Segment *L;
    struct DA *P, *S;
    double eps;
    int curr;
};

static
int AVL_point_comp(int i, int j, void *ctx) {

    struct Sweep_Context *C = ctx;
    struct Point pi; DA_get(C->P, i, &pi);
    struct Point pj; DA_get(C->P, j, &pj);
    return P_point_comp(pi, pj, C->eps);
}

static
int AVL_segment_comp(int i, int j, void *ctx) {

    struct Sweep_Context *C = ctx;
    struct Segment *si = DA_getp(C->S, i);  // this is always S0
    struct Segment *sj = DA_getp(C->S, j);
    struct Point *p = DA_getp(C->P, C->curr);

    double d = point_seg_dist(p, sj);
    int out = (-d > 0) ? 1 : -1;

    if (fabs(d) < C->eps) {     // curr on segment j
        if (si->P[1] != -1) {   // si is a newly inserted segment
            struct Point *q = DA_getp(C->P, si->P[1]);
            if (on_line(q, sj, C->eps)) { return 0; }   // segment exists
        }
        out = ((sj->a - si->a) > 0) ? 1 : -1;
    }
    return out;
}

int line_comp(int i, int j, void *ctx) {

    struct Sweep_Context *C = ctx;
    struct Point p;  DA_get(C->P, C->curr, &p);
    struct Point pi; DA_get(C->P, (C->L)[i].P[1], &pi);
    struct Point pj; DA_get(C->P, (C->L)[j].P[1], &pj);
    double dj = P_distsq(p, pj);
    double di = P_distsq(p, pi);

    return dj - di;
}

void X_V_EV_2_L(double (*V)[2], int En, int (*EV)[2], double (*L)[2][2]) {
    for (int i = 0; i < En; ++i) {
        L[i][0][0] = V[EV[i][0]][0];
        L[i][0][1] = V[EV[i][0]][1];
        L[i][1][0] = V[EV[i][1]][0];
        L[i][1][1] = V[EV[i][1]][1];
    }
}

void X_LA_EL_2_EA(char *LA, int En, int **EL, int *ELn, char *EA) {
    for (int i = 0; i < En; ++i) {
        sort(EL[i], ELn[i], decreasing, NULL);
        EA[i] = LA[EL[i][0]];
    }
}

void sweep_svg(int f, int Pn, struct DA *P, struct DA *V2P,
    struct AVL *Q, struct DA *S);

void X_L_eps_2_V_EV_EL(
    const int Ln, const double (*L)[2][2], const double eps,
    int *Vn, double (**V)[2],
    int *En, int (**EV)[2], int ***EL, int **ELn
) {
    assert((*V == NULL) && (*EV == NULL) &&
        (*EL == NULL) && (*ELn == NULL), "Outputs should start NULL");

    struct Segment *LL = malloc(Ln*sizeof(struct Segment));

    struct DA P = {sizeof(struct Point)};   // struct Point[]
    struct DA PL  = {sizeof(struct DA)};    // int[][]
    struct DA S = {sizeof(struct Segment)}; // struct Segment[]
    struct DA SL  = {sizeof(struct DA)};    // int[][]

    struct Sweep_Context X = {.L = LL, .P = &P, .S = &S, .eps = eps};
    struct AVL Q = {.comp = AVL_point_comp, .ctx = &X};

    for (int li = 0; li < Ln; ++li) {   // compress close input points
        struct Segment *l = LL + li;
        struct Point p[2];
        p[0] = *((struct Point*) &(L[li][0]));
        p[1] = *((struct Point*) &(L[li][1]));

        if (P_point_comp(p[0], p[1], eps) > 0) {  // order ends
            struct Point t = p[0]; p[0] = p[1]; p[1] = t;
        }

        for (int i = 0; i < 2; ++i) {
            l->P[i] = P.n;
            DA_push(&P, &(p[i]));             // add point
            int pi = AVL_insert(&Q, l->P[i]);
            if (pi != l->P[i]) {       // point exists
                l->P[i] = pi;
                DA_pop(&P, NULL);
            }
        }

        struct DA pL = {sizeof(int)};
        while (P.n > PL.n) { DA_push(&PL, &pL); }

        int i0 = l->P[0];
        int i1 = l->P[1];
        struct Point p0; DA_get(&P, i0, &p0);
        struct Point p1; DA_get(&P, i1, &p1);
        struct Point d = P_sub(p1, p0);

        l->U = P_unit(d);
        l->a = (d.y < eps) ? 0 : P_ang(d);
        l->a = (l->a < 0) ? 0 : l->a;
        l->d = P_dot(P_perp(l->U), p0);

        DA_push(DA_getp(&PL, i0), &li);
    }

    double SCALE;
    {
        struct Point p; DA_get(&P, 0, &p);
        double x_min = p.x, y_min = p.y, x_max = p.x, y_max = p.y;
        for (int i = 1; i < P.n; ++i) {
            DA_get(&P, i, &p);
            x_min = (p.x < x_min) ? p.x : x_min;
            y_min = (p.y < y_min) ? p.y : y_min;
            x_max = (p.x > x_max) ? p.x : x_max;
            y_max = (p.y > y_max) ? p.y : y_max;
        }
        double h = y_max - y_min;
        double w = x_max - x_min;
        SCALE = (h > w) ? h : w;
    }

    {   // initialize sentinal segment
        struct Segment s0 = {};
        s0.U.x = -1;
        s0.a = INFINITY;
        s0.P[1] = -1;
        DA_push(&S, &s0);
        struct DA sL0 = {sizeof(int)};
        DA_push(&SL, &sL0);
    }

    struct DA S1  = {sizeof(int)}; // int[]
    struct DA V2P = {sizeof(int)}; // int[]
    struct AVL T = {.comp = AVL_segment_comp, .ctx = &X};

    // int Pn = P.n, cc = 0;
    // char filepath[256];
    // sprintf(filepath, "./out/%.15lf.csv", 100000*eps);
    // printf("%s\n", filepath);
    // FILE *out = fopen(filepath, "w");

    while (AVL_size(&Q) > 0) {    // line sweep
        // AVL_print(&Q);
        // sweep_svg(cc++, Pn, &P, &V2P, &Q, &S);

        int *ppi = AVL_first(&Q);
        int pi = *ppi;
        AVL_remove(&Q, ppi);
        X.curr = pi;
        struct Point p; DA_get(&P, pi, &p);
        struct DA *L = DA_getp(&PL, pi);

        DA_empty(&S1);

        {   // initialize query segment
            struct Segment *s0 = DA_getp(&S, 0);
            s0->P[0] = pi;
            s0->d = P_dot(P_perp(s0->U), p);
        }

        {   // horizontal entering segment
            int *prev = AVL_prev_past(&T, 0);
            if ((prev != NULL) &&
                (((struct Segment*) DA_getp(&S, *prev))->a == 0)
            ) {
                DA_push(&S1, prev);
                AVL_remove(&T, prev);
            }
        }

        while (1) {     // non-horizontal entering segment(s)
            int *next = AVL_next_past(&T, 0);
            if (!next || !on_line(&p, DA_getp(&S, *next), eps)) { break; }
            DA_push(&S1, next);
            AVL_remove(&T, next);
        }

        if ((L->n == 0) && (S1.n < 2)) {  // nothing starting
            if (S1.n == 0) { continue; }    // nothing entering

            int si; DA_get(&S1, 0, &si);  // one segment entering
            struct DA *L1 = DA_getp(&SL, si);
            int ends = 0;

            for (int i = 0; i < L1->n; ++i) {
                int li; DA_get(L1, i, &li);
                struct Point q; DA_get(&P, LL[li].P[1], &q);
                if (P_point_comp(q, p, eps) <= 0) {
                    ends = 1;
                    break;
                }
            }
            if (!ends) {
                AVL_insert(&T, si); // passes through
                continue;           // point unused
            }
        }

        if (S1.n == 1) { // L->n > 0
            int si; DA_get(&S1, 0, &si);
            struct Segment *s = DA_getp(&S, si);
            int all_parallel = 1;
            for (int i = 0; i < L->n; ++i) {
                int li; DA_get(L, i, &li);
                struct Point pL; DA_get(&P, LL[li].P[1], &pL);
                if (!on_line(&pL, s, eps)) {
                    all_parallel = 0;
                }
            }
            if (all_parallel) { // vertex unused (one segment passing)
                AVL_insert(&T, si);
                for (int i = 0; i < L->n; ++i) {
                    int li; DA_get(L, i, &li);
                    DA_push(DA_getp(&SL, si), &li); 
                }
                continue;
            }
        }

        DA_push(&V2P, &pi); // point used as vertex

        for (int i = 0; i < S1.n; ++i) {    // close entering segments
            int si; DA_get(&S1, i, &si);
            struct DA *sL = DA_getp(&SL, si);
            ((struct Segment*) DA_getp(&S, si))->P[1] = pi;

            for (int j = 0; j < sL->n; ++j) {   // add lines that pass through
                int li; DA_get(sL, j, &li);
                struct Point q; DA_get(&P, LL[li].P[1], &q);
                int c = P_point_comp(q, p, eps);
                if (c <= 0) { continue; }   // line ended (possibly before)
                DA_push(L, &li);            // line continues
            }
        }

        sort((int*) L->A, L->n, line_comp, &X);

        for (int i = 0; i < L->n; ++i) {    // add lines to new segments
            int li; DA_get(L, i, &li);
            struct Segment s = LL[li];
            s.P[0] = pi;

            int si = S.n;
            struct DA sL = {sizeof(int)};
            DA_push(&S, &s);
            DA_push(&SL, &sL);
            int sj = AVL_insert(&T, si);    // add segment

            if (sj != si) {                 // existing segment
                DA_pop(&S, &s);
                DA_pop(&SL, NULL);
                DA_push(DA_getp(&SL, sj), &li);
            } else {                        // new segment
                ((struct Segment*) DA_getp(&S, si))->P[1] = -1;
                DA_push(DA_getp(&SL, si), &li);
            }
        }

        for (int i = 0; i < 2; ++i) {       // check for new intersections
            ((struct Segment*) DA_getp(&S, 0))->a = (i == 0) ? -1 : INFINITY;
            int *l = AVL_prev_past(&T, 0);
            if (l == NULL) { continue; }    // no segment to left
            int *r = AVL_next(&T, l);
            if (r == NULL) { continue; }    // no segment to right

            struct Segment *sl = DA_getp(&S, *l);
            struct Segment *sr = DA_getp(&S, *r);
            struct Point pl; DA_get(&P, sl->P[0], &pl);
            struct Point pr; DA_get(&P, sr->P[0], &pr);
            struct Point x = {};

            int c;
            if (!P_line_intersect(  // no intersection
                    pl, P_add(pl, P_mul(sl->U, SCALE)),
                    pr, P_add(pr, P_mul(sr->U, SCALE)), eps, &x)
                || ((c = P_point_comp(x, p, eps)) == 0) // x is p
                || ((c < 0) && (x.y <= p.y))            // x is behind p
            ) { continue; }

            int xi = P.n;
            DA_push(&P, &x);
            int pj = AVL_insert(&Q, xi);    // add intersection
            if (pj == xi) {        // intersection is new
                struct DA pL = {sizeof(int)};
                DA_push(&PL, &pL);
            } else {                // intersection exists
                DA_pop(&P, NULL);
            }
        }
    }
    // printf("(%i)\n", cc);
    // cc = eps*10000;
    // sweep_svg(cc, Pn, &P, &V2P, &Q, &S);
    // printf("%i, %i\n", AVL_size(&Q), AVL_size(&T));
    // fclose(out);

    {   // check termination
        if (AVL_size(&T) != 0) {
            *Vn = *En = 0;
            free(LL);

            AVL_empty(&Q);
            AVL_empty(&T);

            DA_empty(&P);
            DA_empty(&S);
            DA_empty(&V2P);
            DA_empty(&S1);

            for (int i = 0; i < PL.n; ++i) { DA_empty(DA_getp(&PL, i)); }
            DA_empty(&PL);

            for (int i = 0; i < SL.n; ++i) { DA_empty(DA_getp(&SL, i)); }
            DA_empty(&SL);

            return;
        }
        assert(AVL_size(&Q) == 0, "Q not empty\n");
        assert(AVL_size(&T) == 0, "T not empty\n");

        for (int si = 1; si < S.n; ++si) {
            struct Segment *s = DA_getp(&S, si);
            assert(s->P[1] != -1, "Segment not terminated");
        }
    }

    {   // map segment indices from points to vertices
        int *P2V = malloc(P.n*sizeof(int));
        for (int i = 0; i < P.n; ++i) { P2V[i] = -1; }
        for (int i = 0, ii; i < V2P.n; ++i) {
            DA_get(&V2P, i, &ii);
            P2V[ii] = i;
        }

        ((struct Segment*) DA_getp(&S, 0))->P[0] = -1;
        for (int i = 1; i < S.n; ++i) {
            int *sP = ((struct Segment*) DA_getp(&S, i))->P;
            int a = P2V[sP[0]];
            int b = P2V[sP[1]];
            sP[0] = (a < b) ? a : b;
            sP[1] = (a < b) ? b : a;
        }
        free(P2V);
    }

    *Vn = V2P.n;
    *En = S.n - 1;

    *V   = malloc(2*(*Vn)*sizeof(double));
    *EV  = malloc(2*(*En)*sizeof(int));
    *ELn = malloc((*En)*sizeof(int));

    {   // assign vertex data
        for (int i = 0, ii; i < *Vn; ++i) {
            DA_get(&V2P, i, &ii);
            struct Point p; DA_get(&P, ii, &p);
            (*V)[i][0] = p.x;
            (*V)[i][1] = p.y;
        }
    }

    {   // assign edge data
        int *SI = malloc(S.n*sizeof(int));  // sort segments
        for (int i = 0; i < S.n; ++i) { SI[i] = i; }
        sort(SI, S.n, seg_comp, S.A);

        int sum = 0;
        for (int i = 0; i < *En; ++i) {
            int si = SI[i + 1];
            struct Segment *s = DA_getp(&S, si);
            sum += ((*ELn)[i] = ((struct DA*) DA_getp(&SL, si))->n);
            (*EV)[i][0] = s->P[0];
            (*EV)[i][1] = s->P[1];
        }

        *EL = malloc((*En)*sizeof(int*) + sum*sizeof(int));
        int *ELD = (int *) (*EL + *En);

        for (int i = 0; i < *En; ++i) {
            int Ln = (*ELn)[i];
            int *L = (int*) ((struct DA*) DA_getp(&SL, SI[i + 1]))->A;
            int *L_ = (*EL)[i] = ELD;
            ELD += Ln;
            for (int j = 0; j < Ln; ++j) {
                L_[j] = L[j];
            }
        }
        free(SI);
    }

    free(LL);

    AVL_empty(&Q);
    AVL_empty(&T);

    DA_empty(&P);
    DA_empty(&S);
    DA_empty(&V2P);
    DA_empty(&S1);

    for (int i = 0; i < PL.n; ++i) { DA_empty(DA_getp(&PL, i)); }
    DA_empty(&PL);

    for (int i = 0; i < SL.n; ++i) { DA_empty(DA_getp(&SL, i)); }
    DA_empty(&SL);
}

void sweep_svg(int f, int Pn, struct DA *P, struct DA *V2P,
    struct AVL *Q, struct DA *S
) {
    double t = 1000;
    double b = 100;
    struct Point min, diag;
    double ts = P_scale_vals((struct Point*) P->A, Pn, &min, &diag);

    char buff_[256] = {};
    char *buff = buff_;

    sprintf(buff, "./out/output_%i.svg", f);
    FILE *file = fopen(buff, "w");
    SVG_add_svg_open(file, -t/2, -t/2, t, t, -(b + t)/2, -(b + t)/2, b + t, b + t);
    for (int i = 1; i < S->n; ++i) {
        sprintf(buff, "%i", i);
        struct Segment *s = DA_getp(S, i);
        struct Point p; DA_get(P, s->P[0], &p);
        p = P_scale(p, t, ts, min, diag);
        struct Point q = P_add(p, P_mul(s->U, t/30));
        struct Point tc = q;
        char *c = "#F00";
        if (s->P[1] != -1) {
            c = "#00F";
            DA_get(P, s->P[1], &q);
            q = P_scale(q, t, ts, min, diag);
            tc = P_div(P_add(p, q), 2);
        }
        SVG_add_line(file, p.x, p.y, q.x, q.y, 1, c);
        SVG_add_text(file, tc.x, tc.y, c, buff);
    }
    for (int i = 0; i < V2P->n; ++i) {
        int pj; DA_get(V2P, i, &pj);
        sprintf(buff, "%i", pj);
        struct Point p; DA_get(P, pj, &p);
        p = P_scale(p, t, ts, min, diag);
        char *c = (pj < Pn) ? "#F00" : "#000";
        SVG_add_disk(file, p.x, p.y, 2, c);
        SVG_add_text(file, p.x, p.y, "#000", buff);
    }
    for (int *i = AVL_first(Q); i != NULL; i = AVL_next(Q, i)) {
        sprintf(buff, "%i", *i);
        struct Point p; DA_get(P, *i, &p);
        p = P_scale(p, t, ts, min, diag);
        char *c = (*i < Pn) ? "#F00" : "#000";
        SVG_add_disk(file, p.x, p.y, 2, c);
        SVG_add_text(file, p.x, p.y, c, buff);
    }
    SVG_add_svg_close(file);
    fclose(file);
}

double X_L_2_V_EV_EL(
    const int Ln, const double (*L)[2][2],
    int *Vn, double (**V)[2],
    int *En, int (**EV)[2], int ***EL, int **ELn
) {
    assert((*V == NULL) && (*EV == NULL) &&
        (*EL == NULL) && (*ELn == NULL), "Outputs should start NULL");

    double d = INFINITY;
    for (int i = 0; i < Ln; ++i) {  // min line length
        struct Point u = *((struct Point*) &(L[i][0]));
        struct Point v = *((struct Point*) &(L[i][1]));
        double d_ = P_dist(u, v);
        if (d_ < d) { d = d_; }
    }

    const int N = 50, k = 3;
    int nV = 0, nE = 0, count = 0, k_ = 0, i_ = 3;
    int Vn_, En_, (*EV_)[2] = NULL, **EL_ = NULL, *ELn_ = NULL;
    double (*V_)[2] = NULL;
    for (int i = 0; i < N; ++i) {
        free(V_); free(EV_); free(EL_); free(ELn_);
        V_ = NULL; EV_ = NULL; EL_ = NULL; ELn_ = NULL;

        double eps = d/(1 << (i + 3));

        X_L_eps_2_V_EV_EL(Ln, L, eps, &Vn_, &V_, &En_, &EV_, &EL_, &ELn_);

        if (Vn_ == 0) { nV = nE = count = 0; continue; }

        count = ((Vn_ == nV) && (En_ == nE)) ? (count + 1) : 1;
        nV = Vn_; nE = En_;

        if (count <= k_) { continue; }

        free(*V); free(*EV); free(*EL); free(*ELn);
        *Vn = Vn_; *V = V_; *En = En_; *EV = EV_; *EL = EL_, *ELn = ELn_; 
        Vn_ = 0; En_ = 0; V_ = NULL; EV_ = NULL; EL_ = NULL; ELn_ = NULL; 
        k_ = count; i_ = i + 3;

        if (k_ == k) { break; }
    }
    int eps_i = (1 << (i_ - k_));
    double eps = d/eps_i;
    return eps;
}

int by_angle(int i, int j, void *ctx) {
    struct { double (*V)[2]; int curr; } *C = ctx;
    struct Point p  = *((struct Point*) &(C->V[C->curr]));
    struct Point pi = *((struct Point*) &(C->V[i]));
    struct Point pj = *((struct Point*) &(C->V[j]));
    double ai = P_ang(P_sub(pi, p));
    double aj = P_ang(P_sub(pj, p));
    // assert(ai != aj, "Should not have equal angles around vertex");
    return (ai < aj) ? -1 : 1;
}

int by_area(int i, int j, void *ctx) {
    struct { struct DA *F; struct Point *V; } *C = ctx;
    struct DA *Vi = DA_getp(C->F, i);
    struct DA *Vj = DA_getp(C->F, j);
    double ai = P_polygon_area2((int *) Vi->A, Vi->n, C->V);
    double aj = P_polygon_area2((int *) Vj->A, Vj->n, C->V);
    return (ai == aj) ? 0 : ((aj < ai) ? -1 : 1);
}

static
unsigned long long encode(unsigned i, unsigned j) {
    return ((unsigned long long) i << (sizeof(int)*8)) | j;
}

void X_V_EV_2_VV_FV(
    int Vn, double (*V)[2],
    int En, int (*EV)[2],
    int ***VV, int **VVn,
    int *Fn, int ***FV, int **FVn
) {
    assert((*VV == NULL) && (*VVn == NULL) &&
        (*FV == NULL) && (*FVn == NULL), "Outputs should start NULL");

    {
        struct DA *A = calloc(Vn, sizeof(struct DA));
        for (int i = 0; i < Vn; ++i) { A[i].s = sizeof(int); }
        for (int i = 0; i < En; ++i) {
            int a = EV[i][0];
            int b = EV[i][1];
            DA_push(&(A[a]), &b);
            DA_push(&(A[b]), &a);
        }
        *VVn = malloc(Vn*sizeof(int));
        *VV = malloc(Vn*sizeof(int*) + 2*En*sizeof(int));
        int *VVD = (int*) (*VV + Vn);
        for (int i = 0; i < Vn; ++i) {
            struct DA *vA = &(A[i]);
            int vn = (*VVn)[i] = vA->n;
            int *vV = (*VV)[i] = VVD;
            VVD += vn;
            for (int j = 0; j < vn; ++j) {
                DA_get(vA, j, &(vV[j]));
            }
            struct { double (*V)[2]; int curr; } ctx = {V, i};
            sort(vV, vn, by_angle, &ctx);
            DA_empty(vA);
        }
        free(A);
    }
    {
        struct DA F = {sizeof(struct DA)};
        struct HM seen = {sizeof(long long), 0};

        for (unsigned vi = 0; vi < Vn; ++vi) {

            int *Ai = (*VV)[vi];
            for (int j = 0; j < (*VVn)[vi]; ++j) {

                unsigned vj = Ai[j];
                unsigned long long key = encode(vi, vj);

                if (HM_get(&seen, &key, 0)) { continue; }

                HM_set(&seen, &key, 0);
                struct DA fV = {sizeof(int)};
                DA_push(&fV, &vi);
                unsigned a = vi, b = vj;

                while (b != vi) {


                    DA_push(&fV, &b);

                    unsigned t = a; a = b;
                    int *A = (*VV)[b];
                    int An = (*VVn)[b];

                    for (int k = 0; k < An; ++k) {

                        if (A[k] != t) { continue; }

                        b = (k == 0) ? A[An - 1] : A[k - 1];
                        break;
                    }

                    key = encode(a, b);
                    HM_set(&seen, &key, 0);
                }
                assert(fV.n > 2, "Found face with only 2 sides");
                DA_push(&F, &fV);
            }
        }
        HM_empty(&seen);

        int *I = malloc(F.n*sizeof(int));
        int sum = 0;

        for (int i = 0; i < F.n; ++i) {
            I[i] = i;
            sum += ((struct DA*) DA_getp(&F, I[i]))->n;
        }

        struct { struct DA *F; double (*V)[2]; } ctx = {&F, V};
        sort(I, F.n, by_area, &ctx);

        *Fn = F.n - 1;      // exclude outer face with negative area
        *FVn = malloc((*Fn)*sizeof(int));
        *FV = malloc((*Fn)*sizeof(int*) + sum*sizeof(int));
        int *FVD = (int*) (*FV + *Fn);

        DA_empty(DA_getp(&F, I[*Fn])); // empty face
        for (int i = 0; i < *Fn; ++i) {

            struct DA *fV = DA_getp(&F, I[i]);
            int fn = (*FVn)[i] = fV->n;
            int *fV_ = (*FV)[i] = FVD;
            FVD += fn;

            for (int j = 0; j < fn; ++j) {
                DA_get(fV, j, &(fV_[j]));
            }
            DA_empty(fV);
        }
        DA_empty(&F);
        free(I);
    }
}

struct Edata {
    int u, v, f, a, p;
    double d;
};

static
int E_comp(int i, int j, void *ctx) {
    struct Edata *E = ctx;
    return (E[i].d < E[j].d) ? -1 : 1;
}

void X_V_EV_EA_FV_2_Vf_Ff(
    int Vn, double (*V)[2],
    int En, int (*EV)[2], char *EA,
    int Fn, int **FV, int *FVn,
    double (**Vf)[2], int **Ff
) {
    assert((*Vf == NULL) && (*Ff == NULL), "Outputs should start NULL");

    struct Point *P = (struct Point*) V;

    struct Edata *E = malloc(2*En*sizeof(struct Edata));

    // map from directed vertex pair to edge
    struct HM E_map = {sizeof(long long), sizeof(int)};
    for (int ei = 0; ei < En; ++ei) {
        for (int i = 0; i < 2; ++i) {
            unsigned en = 2*ei + i;
            struct Edata *e = E + en;
            e->u = EV[ei][i];
            e->v = EV[ei][(i + 1) % 2];
            e->a = (EA[ei] != 'F');
            e->p = 0;
            e->f = -1;
            struct Point pu = P[e->u];
            struct Point pv = P[e->v];
            e->d = P_dist(pu, pv);
            unsigned long long k = encode(e->u, e->v);
            HM_set(&E_map, &k, &en);
        }
    }

    for (int fi = 0; fi < Fn; ++fi) {
        int fn = FVn[fi];
        int *F = FV[fi];
        int u = F[fn - 1];
        for (int i = 0; i < fn; ++i) {
            int v = F[i];
            unsigned long long k = encode(v, u);
            unsigned i = 0;
            HM_get(&E_map, &k, &i);
            E[i].f = fi;
            u = v;
        }
    }

    *Vf = malloc(2*Vn*sizeof(double));
    *Ff = malloc(Fn*sizeof(int));

    struct Point *Pf = (struct Point*) *Vf;

    for (int i = 0; i < Vn; ++i) { Pf[i].x = INFINITY; }
    for (int i = 0; i < Fn; ++i) { (*Ff)[i] = -1; }

    struct PQ Q = {.comp = E_comp, .ctx = E};

    {
        int v0 = FV[0][0];
        int v1 = FV[0][1];
        Pf[v0] = P[v0];
        Pf[v1] = P[v1];
        unsigned long long k = encode(v1, v0);
        unsigned i = 0;
        HM_get(&E_map, &k, &i);
        PQ_insert(&Q, i);
    }

    while (Q.n > 0) {
        int ei = PQ_extract(&Q);
        struct Edata *e = &(E[ei]);

        assert(e->f >= 0, "Edge went outside paper boundary");
        if ((*Ff)[e->f] != -1) { continue; }
        (*Ff)[e->f] = e->p;

        struct Point x = P_unit(P_sub(P[e->v], P[e->u]));
        struct Point y = P_perp(x);
        struct Point xf = P_unit(P_sub(Pf[e->v], Pf[e->u]));
        struct Point yf = P_perp(xf);

        int fn = FVn[e->f];
        int *F = FV[e->f];
        int u = F[fn - 1];
        for (int i = 0; i < fn; ++i) {
            int v = F[i];
            if (Pf[v].x == INFINITY) {
                struct Point p = P_sub(P[v], P[e->u]);
                struct Point dx = P_mul(xf, P_dot(p, x));
                struct Point dy = P_mul(yf, P_dot(p, y)*(e->p ? -1 : 1));
                Pf[v] = P_add(P_add(dx, dy), Pf[e->u]);
            }
            unsigned long long k = encode(u, v);
            unsigned ej = 0;
            HM_get(&E_map, &k, &ej);
            struct Edata *e_ = &(E[ej]);
            if ((e_->f != -1) && ((*Ff)[e_->f] == -1)) {
                PQ_insert(&Q, ej);
                e_->p = (e_->a) ? !(e->p) : e->p;
            }
            u = v;
        }
    }

    free(E);

    PQ_empty(&Q);
    HM_empty(&E_map);
}

