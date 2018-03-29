#ifndef graph_h
#define graph_h

struct grph_t {
  int mblk;
  int *nblk, *blk;
  int mrow;
  int *nrow, *row;
  int *col;
};

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

#endif
