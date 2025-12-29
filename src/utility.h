#define print_A(A, n, F) do {\
    printf("[");\
    for (int i = 0; i < (n); ++i) {\
        if (i > 0) { printf(","); }\
        printf(F, (A)[i]);\
    }\
    printf("]\n");\
} while (0)

void assert(int condition, char *message);
