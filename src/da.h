struct DA {
    int s, n, m;
    char *A;
};

void DA_push(struct DA *A, void *x);
void DA_pop(struct DA *A, void *x);
void DA_get(struct DA *A, unsigned i, void *x);
void DA_set(struct DA *A, unsigned i, void *x);
void* DA_getp(struct DA *A, unsigned i);

void DA_empty(struct DA *A);
void* DA_freeze(struct DA *A);

void DA_print(struct DA *A, void (*print_item)(void *));
