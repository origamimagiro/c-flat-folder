#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>

struct DA {
    int s, n, m;
    char *A;
};

void DA_push(struct DA *A, void *x) {
    if (A->n == A->m) { A->A = realloc(A->A, (A->m = A->m ? A->m*2 : 2)*A->s); }
    memcpy(A->A + A->n*A->s, (char*) x, A->s);
    ++(A->n);
}

void DA_pop(struct DA *A, void *x) {
    assert(A->n > 0);
    if (4*A->n == A->m) { A->A = realloc(A->A, (A->m /= 2)*A->s); }
    --(A->n);
    if (x != NULL) { memcpy((char*) x, A->A + A->n*A->s, A->s); }
}

void DA_get(struct DA *A, unsigned i, void *x) {
    assert((i >= 0) && (i < A->n));
    memcpy(x, A->A + i*A->s, A->s);
}

void DA_set(struct DA *A, unsigned i, void *x) {
    assert((i >= 0) && (i < A->n));
    memcpy(A->A + i*A->s, x, A->s);
}

void* DA_getp(struct DA *A, unsigned i) {
    assert((i >= 0) && (i < A->n));
    return A->A + i*A->s;
}

void DA_empty(struct DA *A) {
    A->n = A->m = 0;
    memset(A->A, 0, A->m*A->s);
    free(A->A);
    A->A = NULL;
}

void* DA_freeze(struct DA *A) {
    char *out = realloc(A->A, A->n*A->s);
    A->n = A->m = 0;
    return out;
}

//////////////////////////////////////

#include <stdio.h>

void DA_print(struct DA *A, void (*print_item)(void *)) {
    printf("{s: %i, n: %i, m: %i, A:%p}: [", A->s, A->n, A->m, A->A);
    for (int i = 0; i < A->n; ++i) {
        if (i) { printf(","); }
        print_item(A->A + i*A->s);
    }
    printf("]\n");
}

void print_x(void *x) {
    printf("%i", *((int*) x));
}

//////////////////////////////////////

int DA_test() {
    struct DA A = {sizeof(int)};
    DA_print(&A, print_x);
    int n = 10;
    for (int i = 0; i < n; ++i) {
        DA_push(&A, &i);
        DA_print(&A, print_x);
    }
    int x = 0;
    for (int i = 0; i < n; ++i) {
        DA_get(&A, i, &x);
        printf("%i: %i\n", i, x);
    }
    for (int i = 0; i < n + 1; ++i) {
        DA_pop(&A, &x);
        printf("popped: %i\n", x);
        DA_print(&A, print_x);
    }
    DA_empty(&A);
    return 0;
}
