#pragma once
#include "flag.hpp"

string filename_prefix(string type_of_file = "F", string directory = "");

vector<flag> &get_unlabeled_flags(int Kn);
vector<vector<flag>> &get_flags();
vector<flag> &get_forbidden_subflags();

#ifdef G_COLORED_VERTICES
void add_color(int color);
#endif

void generate_unlabeled_flags_of_size(
    int Kn, bool force_generate_flags = false,
    bool remove_duplicates_while_loading = false,
    bool remove_forbidden_wile_loading = false, int verbose_output = 0,
    bool dump_unlabeled_while_generating = false);

// Creates all subflags of flags of size sizeKn
void generate_all_unlabeled_subflags_from_size(int sizeKn, int verbose_output);

void generate_labeled_flags(int flag_size, int type_size, int verbose_output,
                            bool types_from_file = false);

void generate_all_types_containing_one_subtype(int big_type_size,
                                               const flag &type,
                                               vector<flag> &flag_list);

void get_labeled_flags_of_one_type(int flag_size, const flag &type,
                                   vector<flag> &flag_list);

bool dump_unlabeled_flags(int sizeKn);

bool dump_labeled_flags(int sizeKn);

bool load_flags_and_coefficients_from_file(string filename,
                                           vector<flag_coeff> &flag_list,
                                           int verbose_output = 0);

bool load_labeled_flags_from_file(int sizeKn, int verbose_output = 0);

#ifdef G_COLORED_VERTICES
void count_exact_number_of_colored_vertices(int Kn);
#endif

void load_forbidden();
void hydrate_flag_product_cache(int Kn, int idx_in_list, int typesize,
                                ofstream &outfile);

void hydrate_projection_cache(flag type, int vertices, int new_vertices,
                              int remaining_labeled, ofstream &outfile);

void hydrate_multiplication_cache(flag type, int size_left, int size_right,
                                  ofstream &outfile);
