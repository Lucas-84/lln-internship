#ifndef HEL_ENUM_H
#define HEL_ENUM_H

#include "hel_histo.h"
#include "aes.h"


LL*** convert_list_into_array( hel_int_list_t*** key_list, LL nb_bins, LL updated_nb_subkey);

LL** get_binary_hist_from_big_hist(ZZX* hists, LL* big_hist_size,LL* convolution_order, LL updated_nb_subkey);

LL get_index_bin_bound(ZZX* hist, ZZ bound,ZZ* nb_total_elem,LL boolean_end);

void get_enum_info( hel_enum_info_t* enum_info, hel_enum_input_t* enum_input , hel_preprocessing_t* preprocessing, hel_histo_t* histo,hel_real_key_info_t* real_key_info);

LL set_enum_input(hel_enum_input_t* enum_input,  hel_enum_info_t* enum_info , hel_preprocessing_t* hel_preprocessing, hel_histo_t* histo);

LL test_key_equality(unsigned char* ct1, unsigned char* ct2, LL size);

LL print_file_key_original(LL** current_key_facto,LL test_key_boolean, LL updated_nb_subkey, LL merge_value, unsigned char** pt_ct);

LL print_file_key_up_to_key(LL** current_key_facto,LL enumerate_to_real_key_rank, LL updated_nb_subkey, LL merge_value, unsigned char** pt_ct, LL* real_key);

void decompose_bin(LL current_small_hist, LL current_index, hel_preprocessing_t* preprocessing, hel_enum_input_t* enum_input, LL* found_key);

void start_recursive_enumeration( hel_preprocessing_t* preprocessing, hel_histo_t* histo, hel_enum_input_t* enum_input, hel_enum_info_t* enum_info);

#endif
