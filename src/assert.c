#include <stdlib.h>
#include <stdio.h>

void assert(int condition, char *message) {
    if (condition) { return; }
    printf("ERROR: %s\n", message);
    abort();
}
