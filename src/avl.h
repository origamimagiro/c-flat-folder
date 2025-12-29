struct AVL {
    int (*comp)(int, int, void*);
    void *ctx;
    struct { int s, n, m, *A; };
};

int  AVL_size(struct AVL *T);
void AVL_empty(struct AVL *T);

int  AVL_insert(struct AVL *T, int x);
void AVL_remove(struct AVL *T, int *x);

int* AVL_find(struct AVL *T, int x);
int* AVL_next_past(struct AVL *T, int x);
int* AVL_prev_past(struct AVL *T, int x);

int* AVL_first(struct AVL *T);
int* AVL_last(struct AVL *T);

int* AVL_next(struct AVL *T, int *x);
int* AVL_prev(struct AVL *T, int *x);

void AVL_print(struct AVL *T);
void AVL_print_struct(struct AVL *T);

void AVL_test();
