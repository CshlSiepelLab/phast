/****************************************************
rph_util.c
RPHAST low-level utility functions.
These functions should not be needed by the user.

Alexandra Denby
Last updated 4/26/08
****************************************************/
#include <rph_util.h>
#include <stdio.h>
#include <stdlib.h>

int rphast_errno;
char rphast_errmsg[1000];

void init(){
  rphast_errno=0;
  strcpy(rphast_errmsg,"");
}

void* ad2ptr(double address){
  return((void*)(unsigned long)address);
}

double ptr2ad(void* ptr){
  return((double)(unsigned long)ptr);
}

