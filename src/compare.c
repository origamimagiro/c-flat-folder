int max(int a, int b, void *ctx) { return (a < b) ? -1 : (a > b); }
int min(int a, int b, void *ctx) { return (a > b) ? -1 : (a < b); }

int maxA(int a, int b, void *A) { return min(((int*) A)[a], ((int*) A)[b], 0); }
int minA(int a, int b, void *A) { return max(((int*) A)[a], ((int*) A)[b], 0); }

int increasing(int a, int b, void *ctx) { return max(a, b, ctx); }
int decreasing(int a, int b, void *ctx) { return min(a, b, ctx); }

int increasingA(int a, int b, void *A) { return maxA(a, b, A); }
int decreasingA(int a, int b, void *A) { return minA(a, b, A); }

