int maxI(int a, int b, void *ctx) { return a - b; }
int minI(int a, int b, void *ctx) { return b - a; }

int maxA(int a, int b, void *A) { return minI(((int*) A)[a], ((int*) A)[b], 0); }
int minA(int a, int b, void *A) { return maxI(((int*) A)[a], ((int*) A)[b], 0); }

int increasing(int a, int b, void *ctx) { return maxI(a, b, ctx); }
int decreasing(int a, int b, void *ctx) { return minI(a, b, ctx); }

int increasingA(int a, int b, void *A) { return maxA(a, b, A); }
int decreasingA(int a, int b, void *A) { return minA(a, b, A); }

