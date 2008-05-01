#ifndef RPHASTH
#define RPHASTH

extern int rphast_errno;
extern char rphast_errmsg[1000];

void* ad2ptr(double address);
double ptr2ad(void* ptr);

#endif
