#include <stdlib.h>
#include "assert.h"
#include "da.h"

static
void heapify_down(int i, int *A, int n, int (*comp)(int, int, void*), void *ctx) {
    for (int p = i, c = 2*p + 1, t; c < n; p = c, c = 2*p + 1) {
        if ((c + 1 < n) && (comp(A[c], A[c + 1], ctx) < 0)) { ++c; }
        if (comp(A[c], A[p], ctx) <= 0) { break; }
        t = A[p]; A[p] = A[c]; A[c] = t;
    }
}

static
void build_heap(int *A, int n, int (*comp)(int, int, void*), void *ctx) {
    for (int i = n - 1; i >= 0; --i) {
        heapify_down(i, A, n, comp, ctx);
    }
}

void sort(int *A, int n, int (*comp)(int, int, void*), void *ctx) {
    build_heap(A, n, comp, ctx);
    for (int i = n - 1, t; i >= 0; --i) {
        t = A[0]; A[0] = A[i]; A[i] = t;;
        heapify_down(0, A, i, comp, ctx);
    }
}

static
void heapify_up(int i, int *A, int (*comp)(int, int, void*), void *ctx) {
    for (int c = i, p, t; c > 0; c = p) {
        p = (c - 1)/2;
        if (comp(A[c], A[p], ctx) <= 0) { break; }
        t = A[p]; A[p] = A[c]; A[c] = t;
    }
}

struct PQ {
    int (*comp)(int, int, void*);
    void *ctx;
    union {
        struct DA D;
        struct { int s, n, m, *A; };
    };
};

static
int min(int a, int b, void *ctx) { return b - a; }

void PQ_empty(struct PQ *Q) { DA_empty(&(Q->D)); }

void PQ_build(struct PQ *Q, int *A, int n) {
    PQ_empty(Q);
    Q->s = sizeof(int); Q->n = Q->m = n; Q->A = A;
    if (Q->comp == 0) { Q->comp = min; }
    build_heap(A, n, Q->comp, Q->ctx);
}

void PQ_insert(struct PQ *Q, int x) {
    if (Q->comp == 0) { Q->comp = min; }
    if (Q->s == 0) { Q->s = sizeof(int); }
    DA_push(&(Q->D), &x);
    heapify_up(Q->n - 1, Q->A, Q->comp, Q->ctx);
}

int PQ_best(struct PQ *Q) { return Q->A[0]; }

#include <stdio.h>

int PQ_extract(struct PQ *Q) {
    int out; DA_get(&(Q->D), 0, &out);
    int a; DA_pop(&(Q->D), &a);
    if (Q->n > 0) {
        *((int*) DA_getp(&(Q->D), 0)) = a;
        heapify_down(0, Q->A, Q->n, Q->comp, Q->ctx);
    }
    return out;
}

#include <stdio.h>

void PQ_print(struct PQ *Q) {
    printf("[");
    for (int i = 0; i < Q->n; ++i) {
        if (i > 0) { printf(","); }
        printf("%i", Q->A[i]);
    }
    printf("]\n");
}

static char BUFF[256];
static char *B = BUFF;

static
void print_op(struct PQ *Q, char *B) { printf("%-17s | ", B); PQ_print(Q); }

static
struct PQ
test_init() {
    struct PQ Q = {0}; sprintf(B, "Q = {0}"); print_op(&Q, B);
    return Q;
}

static
void test_insert(struct PQ *Q, int i) {
    PQ_insert(Q, i); sprintf(B, "insert(%2i)", i); print_op(Q, B);
}

static
void test_extract(struct PQ *Q) {
    int i = PQ_extract(Q); sprintf(B, "extract() => %-2i", i); print_op(Q, B);
}

static
void test_empty(struct PQ *Q) {
    PQ_empty(Q); sprintf(B, "empty()"); print_op(Q, B);
}

void PQ_test() {
    #include "compare.h"
    int n = 20;
    int *A = malloc(n*sizeof(int));
    for (int i = 0; i < n; ++i) { A[i] = i; }
    for (int i = 0; i < n; ++i) {
        printf("%s%i", (i > 0) ? "," : "[", A[i]);
    } printf("]\n");
    printf("sort decreasing\n");
    sort(A, n, decreasing, 0);
    for (int i = 0; i < n; ++i) {
        printf("%s%i", (i > 0) ? "," : "[", A[i]);
    } printf("]\n");
    struct PQ Q = test_init();
    for (int i = n - 1; i >= 0; --i) { test_insert(&Q, i); }
    for (int i = 0; i < (n/2); ++i) { test_extract(&Q); }
    test_empty(&Q);
    printf("change to max PQ\n");
    Q.comp = max;
    for (int i = 0; i < n; ++i) { test_insert(&Q, i); }
    for (int i = 0; i < n; ++i) { test_extract(&Q); }
}
