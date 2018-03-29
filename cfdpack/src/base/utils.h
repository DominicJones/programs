#ifndef utils_h
#define utils_h

int icmp (const void *a, const void *b);
void uniq (int m0, int* lst0, int* m1, int** lst1);
void setcp (int* lst, int n_lst, int iset);
void setuq (int* lst, int n_lst, int iset);
void print_list (char* name, int m, int* lst);
void split_list (int m, int* lst, int* flg, int* m0, int** lst0, int* m1, int** lst1);
void merge_lists (int m0, int* lst0, int m1, int* lst1, int* m, int** lst);
int offset_lb (int m, int* blk, int val);

#endif
