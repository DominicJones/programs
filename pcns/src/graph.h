// -*- C++ -*-
#pragma once

struct mshf_t;
struct mesh_t;
struct ghst_t;

struct grph_t {
  int mblk;
  int *nblk, *blk;
  int mrow;
  int *nrow, *row;
  int *col;
};

int connectivity (struct mshf_t* mshf, struct mesh_t* mesh);

int match_faces (int n_fce, int* fce_vrt_ptr, int* fce_vrt_lst, int* fce_elms, int* ghst_fce, int** fce_mtch);

void print_range (char const * name, int mrow, int* row);
void print_graph (char const * name, int mrow, int* row, int* col);
void print_igraph (char const * name, int mrow, int* idx, int* row, int* col);
void offset_to_count (int m, int* row, int** nrow);
void count_to_offset (int m, int* nrow, int** row);
void add_pole_entry (int os, int m, int** row, int** col);
void purge_graph (int m, int* row, int** col);

void split_graph (
  int m_grph, int* grph_ia, int* grph_ja, int* flg,
  int* m_prt0, int** prt0_ia, int** prt0_ja,
  int* m_prt1, int** prt1_ia, int** prt1_ja);

void merge_graphs (
  int m_prt0, int* prt0_ia, int* prt0_ja,
  int m_prt1, int* prt1_ia, int* prt1_ja,
  int* m_grph, int** grph_ia, int** grph_ja);

void reorder_graph (int m_grph, int* grph_ia, int* grph_ja, int* ord);
void tuple_to_graph (int m_ord, int m_tup, int* tuple, int* m_grph, int** grph_ia, int** grph_ja);
void graph_to_matrix (int m_grph, int* grph_ia, int* grph_ja, int* m_mat, int** mat_ia, int** mat_ja);
void graph_to_dual (int m_grph, int* grph_ia, int* grph_ja, int* m_dual, int** dual_ia, int** dual_ja);
int shortest_dual_row (int n_row, int* row, int m_dual, int* dual_ia, int* dual_ja);

void mesh_to_dual (int m_msh, int* msh_ia, int* msh_ja, int os, int min_mtch, int** dual_ia, int** dual_ja);
void sort_column_entries(int m_grph, int* grph_ia, int* grph_ja, int** grph_ja_s);


void find_ghost_elements (
  struct ghst_t* ghst,
  int n_elm, int* elm_blk, int* elm_neig_ptr, int* elm_neig_lst);

int communicate_ghost_elements (
  struct ghst_t* ghst, struct ghst_t* locl,
  int* elm_blk,
  int* elm_tag_ptr, int* elm_tag_lst,
  int* elm_vrt_ptr, int* elm_vrt_lst);


int icmp (const void *a, const void *b);
void uniq (int m0, int* lst0, int* m1, int** lst1);
void setcp (int* lst, int n_lst, int iset);
void setuq (int* lst, int n_lst, int iset);
void set_last_cp (int* lst, int n_lst, int iset);
void print_list (char const * name, int m, int* lst);
void split_list (int m, int* lst, int* flg, int* m0, int** lst0, int* m1, int** lst1);
void merge_lists (int m0, int* lst0, int m1, int* lst1, int* m, int** lst);
int offset_lb (int m, int* blk, int val);
