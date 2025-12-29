#include <stdlib.h>
#include <stdio.h>
#include "da.h"

struct Node {
    int x;
    unsigned p, h;
    union {
        unsigned c[2];
        struct { unsigned l, r; };
    };
};

struct AVL {
    int (*comp)(int, int, void*);
    void *ctx;
    struct DA A;
};

/////////////////////////////////////////////////////

static
struct Node* node(struct AVL *T, unsigned i) { return DA_getp(&(T->A), i); }

static
unsigned idx(struct AVL *T, struct Node *X) {
    return ((int) ((char *) X - T->A.A))/sizeof(struct Node);
}

static
int ht(struct AVL *T, int r, struct Node *X) {
    return (X->c[r] == 0) ? 0 : node(T, X->c[r])->h;
}

static
int skew(struct AVL *T, struct Node *X) {
    return ht(T, 1, X) - ht(T, 0, X);
}

static
void update(struct AVL *T, struct Node *X) {
    int hL = ht(T, 0, X), hR = ht(T, 1, X);
    X->h = 1 + ((hL < hR) ? hR : hL);
}

static
struct Node* find(struct AVL *T, int x) {
    struct Node *X = node(T, 0);
    while (1) {
        int c = T->comp(x, X->x, T->ctx);
        if (c == 0) { break; }
        unsigned j = (c < 0) ? X->l : X->r;
        if (j == 0) { break; }
        X = node(T, j);
    }
    return X;
}

static
void rotate(struct AVL *T, int r, struct Node *D) {
    unsigned d = idx(T, D), b = D->c[!r], e = D->c[r], Dx = D->x;
    struct Node *B = node(T, b);
    unsigned a = B->c[!r], c = B->c[r], Bx = B->x;
    B->c[!r] = c; B->c[r] = e; B->x = Dx;
    D->c[!r] = a; D->c[r] = b; D->x = Bx; //  __D_      _B__
    if (a != 0) { node(T, a)->p = d; }    // _B_ E  =>  A _D_
    if (e != 0) { node(T, e)->p = b; }    // A C          C E
    update(T, B);
    update(T, D);
}

static
void maintain(struct AVL *T, struct Node *X) {
    while (1) {
        update(T, X);
        int s = skew(T, X);
        int r = !(s >= 2);
        if ((!r) || (s <= -2)) {
            struct Node *Y = node(T, X->c[!r]);
            if (skew(T, Y) == ((2*r) - 1)) { rotate(T, !r, Y); }
            rotate(T, r, X);
        }
        if (idx(T, X) == 0) { break; }
        X = node(T, X->p);
    }
}

static
struct Node* end(struct AVL *T, int r, struct Node *X) {
    while (X->c[r] != 0) { X = node(T, X->c[r]); }
    return X;
}

static
struct Node* adj(struct AVL *T, int r, struct Node *X) {
    if (T->A.n == 0) { return NULL; }
    if ((X == NULL) || (X->c[r] != 0)) {
        return end(T, !r, node(T, (X == NULL) ? 0 : X->c[r]));
    }
    struct Node *Y = X;
    int has = 0;
    while ((idx(T, Y) != 0) && !(has = (node(T, node(T, Y->p)->c[!r]) == Y))) {
        Y = node(T, Y->p);
    }
    return has ? node(T, Y->p) : NULL;
}

static
void remove_key(struct AVL *T, struct Node *X) {
    int r = (X->c[0] == 0);
    unsigned c = X->c[r];
    while (c != 0) {
        while (node(T, c)->c[!r] != 0) {
            c = node(T, c)->c[!r];
        }
        X->x = node(T, c)->x;
        X = node(T, c);
        c = X->c[r];
    }
    c = X->p;
    unsigned i = idx(T, X);
    r = (node(T, c)->r == i);
    struct Node Y;
    DA_pop(&(T->A), &Y);
    unsigned j = T->A.n;
    if (T->A.n == 0) { return; }
    if (i != j) {
        node(T, Y.p)->c[node(T, Y.p)->r == j] = i;
        if (Y.l) { node(T, Y.l)->p = i; }
        if (Y.r) { node(T, Y.r)->p = i; }
        *node(T, i) = Y;
        c = (c == j) ? i : c;
    }
    node(T, c)->c[r] = 0;
    maintain(T, node(T, c));
}

static
struct Node* neighbor(struct AVL *T, int r, int x) {
    if (T->A.n == 0) { return NULL; }
    struct Node *X = find(T, x);
    int c = T->comp(x, X->x, T->ctx);
    if ((c == 0) || (r != (c < 0))) { X = adj(T, r, X); }
    return X;
}

/////////////////////////////////////////////////////

unsigned AVL_size(struct AVL *T) { return T->A.n; }
void AVL_empty(struct AVL *T) { DA_empty(&(T->A)); *T = (struct AVL) {}; }

int AVL_insert(struct AVL *T, int x) {
    if (T->A.s == 0) { T->A.s = sizeof(struct Node); }
    unsigned i; int c;
    if (T->A.n > 0) {
        i = idx(T, find(T, x));
        c = T->comp(x, node(T, i)->x, T->ctx);
        if (c == 0) { return node(T, i)->x; }
    }
    unsigned j = T->A.n;
    DA_push(&(T->A), &((struct Node) {}));
    struct Node *Y = node(T, j);
    if (T->A.n > 1) {
        node(T, i)->c[c > 0] = j;
        Y->p = i;
    }
    Y->x = x;
    maintain(T, Y);
    return x;
}

void AVL_remove(struct AVL *T, int *x) { remove_key(T, (struct Node*) x); }

int* AVL_find(struct AVL *T, int x) {
    struct Node *X = find(T, x);
    return T->comp(x, X->x, T->ctx) ? NULL : ((int *) X);
}

int* AVL_next_past(struct AVL *T, int x) { return (int*) neighbor(T, 1, x); }
int* AVL_prev_past(struct AVL *T, int x) { return (int*) neighbor(T, 0, x); }

int* AVL_first(struct AVL *T) { return T->A.n ? &(end(T, 0, node(T, 0))->x) : NULL; }
int* AVL_last (struct AVL *T) { return T->A.n ? &(end(T, 1, node(T, 0))->x) : NULL; }

int* AVL_next(struct AVL *T, int *x) { return (int *) adj(T, 1, (struct Node*) x); }
int* AVL_prev(struct AVL *T, int *x) { return (int *) adj(T, 0, (struct Node*) x); }

/////////////////////////////////////////////////////

static
void print_node(struct Node *X) {
    if (X == NULL) { return; }
    printf("(x: %i, h: %i, p: %u, l: %u, r: %u)", X->x, X->h, X->p, X->l, X->r);
}

void AVL_print_struct(struct AVL *T) {
    printf("AVL tree: %p\n", T);
    printf("  comp: %p, ctx: %p, n: %i\n", T->comp, T->ctx, T->A.n);
    for (int i = 0; i < T->A.n; ++i) {
        printf("%i: ", i); print_node(node(T, i)); printf("\n");
    }
}

static
int node_sprint(struct Node *X, char *buff) {
    return snprintf(buff, buff ? 16 : 0, "%i", X->x);
}

static
int str_width(struct AVL *T, int i) {
    struct Node *X = node(T, i);
    int out = node_sprint(X, NULL);
    out += X->l ? str_width(T, X->l) : 0;
    out += X->r ? str_width(T, X->r) : 0;
    return out;
}

static
int str_gen(struct AVL *T, unsigned i, int *x, int y, int l, int w, char *S) {
    struct Node *X = node(T, i);
    char *row = S + y*w;
    int x0 = X->l ? str_gen(T, X->l, x, y + 1, 1, w, S) : *x;
    int x1 = *x;
    for (int j = x0; j < x1; ++j) { row[j] = '_'; }
    char c = *(row + (*x) + node_sprint(X, NULL));
    int x2 = (*x += node_sprint(X, row + (*x)));
    row[x2] = c;
    int x3 = X->r ? str_gen(T, X->r, x, y + 1, 0, w, S) : *x;
    for (int j = x2; j < x3; ++j) { row[j] = '_'; }
    return l ? x1 : x2;
}

void AVL_print(struct AVL *T) {
    if (T->A.n == 0) { return; }
    int w = str_width(T, 0) + 1;
    int h = node(T, 0)->h;
    char *S = malloc((h*w + 1)*sizeof(char));
    for (int i = 0; i < h*w; ++i) { S[i] = ((i + 1) % w) ? ' ' : '\n'; }
    S[h*w] = 0;
    int x = 0;
    str_gen(T, 0, &x, 0, 0, w, S);
    printf("%s\n", S);
    free(S);
}

int AVL_verbose = 1;

void AVL_test() {
    #include "compare.h"
    int n = 50;
    AVL_verbose = 1;
    struct AVL T = {.comp = decreasing};
    printf("Testing AVL_insert\n");
    for (int x = 0; x < n; ++x) {
        AVL_insert(&T, x);
        printf("\ninserted %i:\n", x);
        AVL_print(&T);
    }
    printf("Testing AVL_empty\n");
    AVL_empty(&T);
    T.comp = increasing;
    for (int x = 0; x < n; ++x) {
        AVL_insert(&T, x);
        printf("\ninserted %i:\n", x);
        AVL_print(&T);
    }
    printf("Testing forward traversal:\n");
    for (int *x = AVL_first(&T); x != NULL; x = AVL_next(&T, x)) {
        printf("%i,", *x);
    }
    printf("\n");
    printf("Testing backward traversal:\n");
    for (int *x = AVL_last(&T); x != NULL; x = AVL_prev(&T, x)) {
        printf("%i,", *x);
    }
    printf("\n");
}

// int main() { AVL_test(); }
