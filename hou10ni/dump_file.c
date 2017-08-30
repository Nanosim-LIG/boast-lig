/*
  File used to record all input parameters in binary format for BOAST 
*/

#include <stdlib.h>
#include <stdio.h>

void dump_file_(char* filename, void* mem, size_t sz){
        FILE *f = fopen(filename, "wb");
        fwrite(mem,1,sz,f);
        fclose(f);
}

int main(void){
        int cost;
        size_t sz=sizeof(cost);
        cost=0;
        dump_file_("cost.in",&cost,sz);
        return 0;
}
