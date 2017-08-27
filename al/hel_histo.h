#ifndef HEL_HISTO_H
#define HEL_HISTO_H

#include "hel_init.h"

LL find_real_key_bin(hel_histo_t* histo, hel_data_t* data, hel_preprocessing_t* preprocessing);

void compute_single_histogram(hel_algo_mode_t algo_mode, ZZX* hist, LL subkey_number, hel_data_t* data, hel_histo_t* histo, hel_enum_input_t* enum_input, hel_preprocessing_t* preprocess);

void swap_column(LL** target, LL index_col_1, LL index_col_2);

int cmp_array( const void* a, const void *b);

void search_convo_order(LL** index, LL* convo_order, LL updated_nb_subkey);

ZZX* compute_histograms_procedure(hel_algo_mode_t algo_mode, hel_preprocessing_t* preprocessing, hel_data_t* data, hel_histo_t* histo, hel_enum_input_t* enum_input);

LL hel_compute_histogram(hel_algo_mode_t algo_mode, hel_preprocessing_t* preprocessing, hel_histo_t* histo, hel_data_t* data, hel_enum_input_t* enum_input);

void get_real_key_info(hel_real_key_info_t* real_key_info, hel_preprocessing_t* preprocessing, hel_histo_t* histo, hel_data_t* data );

#endif
