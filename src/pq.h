void sort(int *A, int n, int (*comp)(int, int, void*), void *ctx);

struct PQ {
    int (*comp)(int, int, void*);
    void *ctx;
    struct { int s, n, m, *A; };
};

void PQ_build(struct PQ *Q, int *A, int n);
int  PQ_best(struct PQ *Q);
void PQ_insert(struct PQ *Q, int x);
int  PQ_extract(struct PQ *Q);
void PQ_empty(struct PQ *Q);
void PQ_test();
