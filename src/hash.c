#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "assert.h"

const unsigned SI = sizeof(int);        // Requires = 4
const unsigned SL = sizeof(long long);  // Requires = 8

struct HM {        // 48 bytes
    unsigned k, v;  // number of bytes in key and value resepctively
    int (*eq)(unsigned, const void*, const void*);      // key equality
    unsigned long long (*hash)(unsigned, const void*);  // key to hash
    unsigned n;     // # stored items
    unsigned m;     // size of table
    unsigned f;     // index of next free
    unsigned *A;    // next data (size 2*m)
    char *K, *V;    // key, val data (each size 2*m)
};

unsigned long long HM_HASH_INT(unsigned kn, const void *k) {
    unsigned long long h = 0; memcpy(&h, k, (kn <= SL) ? kn : SL);
    return ((11400714785074694791ULL)*(h + ((kn <= SL) ? 0 :
        HM_HASH_INT(kn - SL, ((char *) k) + SL)))) % ((unsigned long) -59);
}
int HM_EQ_INT(unsigned kn, const void *a, const void *b) {
    return !memcmp(a, b, kn);
}
unsigned long long HM_HASH_STR(unsigned kn, const void *k) {
    return HM_HASH_INT(strlen(*((char **) k)), *((char **) k));
}
int HM_EQ_STR(unsigned kn, const void *a, const void *b) {
    return !strcmp(*((char **) a), *((char **) b));
}

int HM_set(struct HM *M, void *k, void *v);
static
void rebuild(struct HM *M, unsigned m) {
    if (M->eq == 0) {
        M->eq   = HM_EQ_INT;
        M->hash = HM_HASH_INT;
    }
    unsigned *A = M->A;
    char *K = M->K, *V = M->V;
    M->A = calloc(2*m, sizeof(int));
    M->K = calloc(2*m, M->k);
    M->V = calloc(2*m, M->v);
    M->f = m; M->n = 0;
    unsigned m_ = M->m; M->m = m;
    char *k = K, *v = V;
    for (int i = 0; i < 2*m_; ++i) {
        if (A[i]) { HM_set(M, k, v); }
        k += M->k;
        v += M->v;
    }
    free(A); free(K); free(V);
}

static
int find(struct HM *M, void *k, unsigned *i) {
    unsigned a = M->A[*i = (M->hash(M->k, k) % M->m)];
    if (a == 0) { return 0; }
    while (!M->eq(M->k, k, M->K + (*i)*M->k)) {
        if (a == 1) { return 0; }
        a = M->A[*i = a];
    }
    return 1;
}

int HM_get(struct HM *M, void *k, void *v) {
    unsigned i;
    if ((M->m == 0) || !find(M, k, &i)) { return 0; }
    memcpy(v, M->V + i*M->v, M->v);
    return 1;
}

int HM_set(struct HM *M, void *k, void *v) {
    if (M->n == M->m) { rebuild(M, M->m ? 2*M->m : 2); }
    unsigned i;
    if (find(M, k, &i)) { memcpy(M->V + i*M->v, v, M->v); return 1; }
    if (M->A[i]) {
        i = M->A[i] = M->f;
        do { M->f = ((M->f + 1 - M->m) % M->m) + M->m; } while (M->A[M->f]);
    }
    M->A[i] = 1;
    memcpy(M->K + i*M->k, k, M->k);
    memcpy(M->V + i*M->v, v, M->v);
    ++(M->n);
    return 0;
}

int HM_del(struct HM *M, void *k, void *v) {
    if (M->m == 0) { return 0; }
    if ((2*M->n == M->m) && (M->m > 2)) { rebuild(M, M->m/2); }
    unsigned i;
    if (!find(M, k, &i)) { return 0; }
    memcpy(v, M->V + i*M->v, M->v);
    unsigned a = M->A[i];
    if (a != 1) {
        memcpy(M->K + i*M->k, M->K + a*M->k, M->k);
        memcpy(M->V + i*M->v, M->V + a*M->v, M->v);
        M->A[i] = M->A[a];
        i = a;
    }
    memset(M->K + i*M->k, 0, M->k);
    memset(M->V + i*M->v, 0, M->v);
    M->A[i] = 0;
    --(M->n);
    return 1;
}

void HM_empty(struct HM *M) {
    M->n = M->m = M->f = 0;
    free(M->A); free(M->K); free(M->V);
    M->A = NULL; M->K = NULL; M->V = NULL;
}

/////////////////////////////////////////////////

void HM_print(struct HM *M, int verbose) {
    printf("{k: %u, v: %u, m: %u, n: %u, f: %u}\n",
            M->k, M->v, M->m, M->n, M->f);
    if (!verbose) { return; }
    for (unsigned i = 0; i < 2*M->m; ++i) {
        if (i == M->m) { printf("-----\n"); }
        printf("   %u: {k: ", i);
        char *Mk = M->K + i*M->k;
        if (M->eq == HM_EQ_STR) {
            printf("%s", *((char **) Mk));
        } else if ((M->eq == HM_EQ_INT) && (M->k == SL)) {
            printf("%llu", *((unsigned long long *) Mk));
        } else {
            for (int j = M->k - 1; j >= 0; --j) {
                printf("%.2X", Mk[j] & 0xFF);
            }
        }
        printf(", v: ");
        char *Mv = M->V + i*M->v;
        if (M->v == SL) {
            printf("%llu", *((unsigned long long *) Mv));
        } else {
            for (int j = M->v - 1; j >= 0; --j) {
                printf("%.2X", Mv[j] & 0xFF);
            }
        }
        printf(", next: %u}\n", M->A[i]);
    }
}

static
unsigned long long rand_ul() {
    return (((long long unsigned) rand()) << sizeof(int)) + rand();
}

int HM_test() {
    int n = 10;
    srand(0);
    unsigned long long *K = malloc(n*sizeof(long long));
    unsigned long long *V = malloc(n*sizeof(long long));
    for (int i = 0; i < n; ++i) {
        K[i] = rand_ul();
        V[i] = rand_ul();
    }
    {
        struct HM M = {sizeof(long long), sizeof(long long)};
        HM_print(&M, 1);
        for (int i = 0; i < n; ++i) {
            printf("Inserting: %llu, %llu\n", K[i], V[i]);
            assert(!HM_set(&M, K + i, V + i), "key already exists");
            HM_print(&M, 1);
        }
        for (int i = 0; i < n; ++i) {
            unsigned long long v;
            printf("Getting: %llu\n", K[i]);
            assert(HM_get(&M, K + i, &v), "key not found");
            printf("Found: %llu\n", v);
            assert(v == V[i], "value does not match");
        }
        for (int i = 0; i < n; ++i) {
            unsigned long long v;
            printf("Deleting: %llu\n", K[i]);
            assert(HM_del(&M, K + i, &v), "key not found");
            assert(v == V[i], "value does not match");
        }
        HM_print(&M, 1);
        HM_empty(&M);
        struct HM S = {sizeof(char*), sizeof(long long), HM_EQ_STR, HM_HASH_STR};
        char (*C)[16] = calloc(1, n*16);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 15; ++j) {
                C[i][j] = 'a' + (rand() % 26);
            }
            char *k = C[i];
            printf("\n**** setting (%s, %llu) ****\n", k, V[i]);
            assert(!HM_set(&S, &k, V + i), "key already exists");
            HM_print(&S, 1);
        }
        HM_print(&S, 1);
        for (int i = 0; i < n; ++i) {
            unsigned long long v;
            char *k = C[i];
            assert(HM_get(&S, &k, &v), "key not found");
            assert(v == V[i], "value does not match");
        }
        for (int i = 0; i < n; ++i) {
            unsigned long long v;
            char *k = C[i];
            assert(HM_del(&S, &k, &v), "key not found");
            assert(v == V[i], "value does not match");
            printf("\n**** deleting (%s, %llu) ****\n", k, v);
        }
        HM_print(&S, 1);
        HM_empty(&S);
    }
    {
        struct HM M = {sizeof(int), sizeof(int)};
        HM_print(&M, 1);
        for (int i = 0; i < n; ++i) {
            unsigned k = (unsigned) K[i], v = (unsigned) V[i];
            printf("Inserting: %X, %X\n", k, v);
            assert(!HM_set(&M, &k, &v), "key already exists");
            HM_print(&M, 1);
        }
        for (int i = 0; i < n; ++i) {
            unsigned k = (unsigned) K[i], v = (unsigned) V[i];
            unsigned v_ = 0;
            printf("Getting: %X\n", k);
            assert(HM_get(&M, &k, &v_), "key not found");
            assert(v == v_, "value does not match");
        }
        for (int i = 0; i < n; ++i) {
            unsigned k = (unsigned) K[i], v = (unsigned) V[i];
            unsigned v_ = 0;
            printf("Deleting: %X\n", k);
            assert(HM_del(&M, &k, &v_), "key not found");
            assert(v == v_, "value does not match");
        }
        HM_print(&M, 1);
        HM_empty(&M);
    }
    {
        struct HM M = {sizeof(long long)};
        HM_print(&M, 1);
        for (int i = 0; i < n; ++i) {
            printf("Adding: %llu\n", K[i]);
            assert(!HM_set(&M, K + i, 0), "key already exists");
            HM_print(&M, 1);
        }
        for (int i = 0; i < n; ++i) {
            printf("Has: %llu\n", K[i]);
            assert(HM_get(&M, K + i, 0), "key not found");
        }
        for (int i = 0; i < n; ++i) {
            printf("Remove: %llu\n", K[i]);
            assert(HM_del(&M, K + i, 0), "key not found");
        }
        HM_print(&M, 1);
        HM_empty(&M);
    }
    return 0;
}

// int main() { HM_test(); }
