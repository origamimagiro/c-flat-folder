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

int HM_set(struct HM *M, void *kp, void *vp);
int HM_get(struct HM *M, void *kp, void *vp);
int HM_del(struct HM *M, void *kp, void *vp);

int HM_next(struct HM *M, void *kp, void *vp);

void HM_empty(struct HM *M);

unsigned long long HM_HASH_INT(unsigned k, const void *kp);
unsigned long long HM_HASH_STR(unsigned k, const void *kp);

int HM_EQ_INT(unsigned k, const void *a, const void *b);
int HM_EQ_STR(unsigned k, const void *a, const void *b);

void HM_print(struct HM *M, int verbose);
