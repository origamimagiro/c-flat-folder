#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

struct HM {        // 48 bytes
    unsigned k, v;  // number of bytes in key and value resepctively
    int (*eq)(unsigned, const void*, const void*);    // key equality
    unsigned long (*hash)(unsigned, const void*);     // key to hash
    unsigned s;     // k + v + sizeof(unsigned)
    unsigned n;     // # stored items
    unsigned m;     // size of table
    unsigned f;     // index of next free
    char *A;        // data (2*m*s bytes allocated)
};

unsigned long HM_HASH_INT(unsigned kn, const void *k) {
    unsigned long h = 0; memcpy(&h, k, (kn <= 8) ? kn : 8);
    return ((11400714785074694791UL)*(h + ((kn <= 8) ? 0 :
        HM_HASH_INT(kn - 8, ((char *) k) + 8)))) % ((unsigned long) -59);
}
int HM_EQ_INT(unsigned kn, const void *a, const void *b) {
    return !memcmp(a, b, kn);
}
unsigned long HM_HASH_STR(unsigned kn, const void *k) {
    return HM_HASH_INT(strlen(*((char **) k)), *((char **) k));
}
int HM_EQ_STR(unsigned kn, const void *a, const void *b) {
    return !strcmp(*((char **) a), *((char **) b));
}

const unsigned SU = sizeof(unsigned);
int HM_set(struct HM *M, void *k, void *v);
static
void rebuild(struct HM *M, unsigned m) {
    if (M->m == 0) {
        assert((M->k > 0) && !(M->k % SU) && !(M->v % SU));
        M->s = SU + M->k + M->v;
        if (M->eq == 0) {
            M->eq   = HM_EQ_INT;
            M->hash = HM_HASH_INT;
        }
    }
    char *A = M->A;
    M->A = calloc(2*m, M->s);
    M->f = m; M->n = 0;
    unsigned m_ = M->m; M->m = m;
    char *c = A;
    for (int i = 0; i < 2*m_; ++i, c += M->s) {
        if (!(*((unsigned *) c))) { continue; }
        HM_set(M, c + SU, c + SU + M->k);
    }
    free(A);
}

static
int find(struct HM *M, void *k, char **c) {
    *c = M->A + (M->hash(M->k, k) % M->m)*M->s;
    unsigned a = *((unsigned *) *c);
    if (a == 0) { return 0; }
    while (!M->eq(M->k, k, *c + SU)) {
        if (a == 1) { return 0; }
        a = *((unsigned *) (*c = M->A + a*M->s));
    }
    return 1;
}

int HM_get(struct HM *M, void *k, void *v) {
    char *c;
    if ((M->m == 0) || !find(M, k, &c)) { return 0; }
    memcpy(v, c + SU + M->k, M->v);
    return 1;
}

int HM_set(struct HM *M, void *k, void *v) {
    if (M->n == M->m) { rebuild(M, M->m ? 2*M->m : 2); }
    char *c;
    if (find(M, k, &c)) { memcpy(c + SU + M->k, v, M->v); return 1; }
    if (*((unsigned *) c)) {
        c = M->A + (*((unsigned *) c) = M->f)*M->s;
        do { M->f = ((M->f + 1 - M->m) % M->m) + M->m; }
        while (*((unsigned *) (M->A + M->f*M->s)));
    }
    *((unsigned *) c) = 1;
    memcpy(c + SU, k, M->k);
    memcpy(c + SU + M->k, v, M->v);
    ++(M->n);
    return 0;
}

int HM_del(struct HM *M, void *k, void *v) {
    if (M->m == 0) { return 0; }
    if ((2*M->n == M->m) && (M->m > 2)) { rebuild(M, M->m/2); }
    char *c;
    if (!find(M, k, &c)) { return 0; }
    memcpy(v, c + SU + M->k, M->v);
    if (*((unsigned *) c) != 1) {
        char *c_ = M->A + *((unsigned *) c)*M->s;
        memcpy(c, c_, M->s);
        c = c_;
    }
    memset(c, 0, M->s);
    --(M->n);
    return 1;
}

void HM_empty(struct HM *M) { M->n = M->m = M->f = 0; free(M->A); }

int HM_next(struct HM *M, void *kp, void *vp) {
    char *k = *((char **) kp);
    char *a = (k == NULL) ? M->A : (k + M->k + M->v);
    char *lim = M->A + 2*M->m*M->s;
    while ((a < lim) && (*((unsigned *) a) == 0)) { a += M->s; }
    if (a == lim) {
        *((char **) kp) = *((char **) vp) = NULL;
        return 0;
    }
    *((char **) vp) = (*((char **) kp) = a + SU) + M->k;
    return 1;
}

/////////////////////////////////////////////////

void HM_print(struct HM *M, int verbose) {
    printf("{k: %u, v: %u, m: %u, n: %u, f: %u}\n",
            M->k, M->v, M->m, M->n, M->f);
    if (!verbose) { return; }
    for (unsigned i = 0; i < 2*M->m; ++i) {
        if (i == M->m) { printf("-----\n"); }
        printf("   %u: {k: ", i);
        char *c = M->A + i*M->s;
        char *Mk = c + SU;
        if (M->eq == HM_EQ_STR) {
            printf("%s", *((char **) Mk));
        } else if ((M->eq == HM_EQ_INT) && (M->k == 8)) {
            printf("%lu", *((unsigned long *) Mk));
        } else {
            for (int j = M->k - 1; j >= 0; --j) {
                printf("%.2X", Mk[j] & 0xFF);
            }
        }
        printf(", v: ");
        char *Mv = c + SU + M->k;
        if (M->v == 8) {
            printf("%lu", *((unsigned long *) Mv));
        } else {
            for (int j = M->v - 1; j >= 0; --j) {
                printf("%.2X", Mv[j] & 0xFF);
            }
        }
        printf(", next: %u}\n", *((unsigned *) c));
    }
}

static
unsigned long rand_ul() {
    return (((long unsigned) rand()) << 32) + rand();
}

int HM_test() {
    int n = 10;
    srand(0);
    unsigned long *K = malloc(n*sizeof(unsigned long));
    unsigned long *V = malloc(n*sizeof(unsigned long));
    for (int i = 0; i < n; ++i) {
        K[i] = rand_ul();
        V[i] = rand_ul();
    }
    {
        struct HM M = {sizeof(long), sizeof(long)};
        HM_print(&M, 1);
        for (int i = 0; i < n; ++i) {
            printf("Inserting: %lu, %lu\n", K[i], V[i]);
            assert(!HM_set(&M, K + i, V + i));
            HM_print(&M, 1);
        }
        for (int i = 0; i < n; ++i) {
            unsigned long v;
            printf("Getting: %lu\n", K[i]);
            assert(HM_get(&M, K + i, &v));
            printf("Found: %lu\n", v);
            assert(v == V[i]);
        }
        unsigned long *k = 0, *v = 0, i = 0;
        while (HM_next(&M, &k, &v)) {
            printf("%lu: [%lu, %li]\n", i++, *k, *v);
        }
        for (int i = 0; i < n; ++i) {
            unsigned long v;
            printf("Deleting: %lu\n", K[i]);
            assert(HM_del(&M, K + i, &v));
            assert(v == V[i]);
        }
        HM_print(&M, 1);
        HM_empty(&M);
        struct HM S = {sizeof(long), sizeof(long), HM_EQ_STR, HM_HASH_STR};
        char (*C)[16] = calloc(1, n*16);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 15; ++j) {
                C[i][j] = 'a' + (rand() % 26);
            }
            char *k = C[i];
            printf("\n**** setting (%s, %lu) ****\n", k, V[i]);
            assert(!HM_set(&S, &k, V + i));
            HM_print(&S, 1);
        }
        HM_print(&S, 1);
        for (int i = 0; i < n; ++i) {
            unsigned long v;
            char *k = C[i];
            assert(HM_get(&S, &k, &v));
            assert(v == V[i]);
        }
        for (int i = 0; i < n; ++i) {
            unsigned long v;
            char *k = C[i];
            assert(HM_del(&S, &k, &v));
            assert(v == V[i]);
            printf("\n**** deleting (%s, %lu) ****\n", k, v);
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
            assert(!HM_set(&M, &k, &v));
            HM_print(&M, 1);
        }
        for (int i = 0; i < n; ++i) {
            unsigned k = (unsigned) K[i], v = (unsigned) V[i];
            unsigned v_ = 0;
            printf("Getting: %X\n", k);
            assert(HM_get(&M, &k, &v_));
            assert(v == v_);
        }
        for (int i = 0; i < n; ++i) {
            unsigned k = (unsigned) K[i], v = (unsigned) V[i];
            unsigned long v_ = 0;
            printf("Deleting: %X\n", k);
            assert(HM_del(&M, &k, &v_));
            assert(v == v_);
        }
        HM_print(&M, 1);
        HM_empty(&M);
    }
    {
        struct HM M = {sizeof(long)};
        HM_print(&M, 1);
        for (int i = 0; i < n; ++i) {
            printf("Adding: %lX\n", K[i]);
            assert(!HM_set(&M, K + i, 0));
            HM_print(&M, 1);
        }
        for (int i = 0; i < n; ++i) {
            printf("Has: %lX\n", K[i]);
            assert(HM_get(&M, K + i, 0));
        }
        for (int i = 0; i < n; ++i) {
            printf("Remove: %lX\n", K[i]);
            assert(HM_del(&M, K + i, 0));
        }
        HM_print(&M, 1);
        HM_empty(&M);
    }
    return 0;
}

// int main() { HM_test(); }
