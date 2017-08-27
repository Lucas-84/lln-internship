#ifndef HEL_INIT_H
#define HEL_INIT_H

#include "hel_util.h"



//preprocessing
hel_preprocessing_t* hel_alloc_preprocessing();
LL hel_init_preprocessing( hel_preprocessing_t* preprocessing, LL merge_value);
void hel_free_preprocessing( hel_preprocessing_t* preprocessing );


//histogram
hel_histo_t* hel_alloc_histo();
LL hel_init_histo(hel_histo_t* histo, LL nb_bins);
void hel_free_histo(hel_histo_t* histo);


//data
hel_data_t* hel_alloc_data();
LL hel_init_data(hel_data_t* data, hel_preprocessing_t* preprocessing, double** score_mat_init, LL* key);
void hel_free_data( hel_data_t* data, hel_preprocessing_t* preprocessing);


//enum param
hel_enum_input_t* hel_alloc_enum_input();
LL hel_init_enum_input(hel_enum_input_t* enum_input, hel_preprocessing_t* preprocessing, hel_data_t* data, hel_histo_t* histo , ZZ bound_start, ZZ bound_end, LL test_key_boolean ,unsigned char** pt_ct, LL enumerate_to_real_key_rank, LL up_to_bound_bool);
void hel_free_enum_input(hel_enum_input_t* enum_input, hel_histo_t* histo, hel_preprocessing_t* preprocessing);




//param
hel_param_t* hel_alloc_param();
LL hel_init_param(hel_param_t* param, hel_algo_mode_t algo_mode, LL merge_value, LL nb_bins, double** score_mat_init, LL* key_init, ZZ bound_start, ZZ bound_end , LL test_key_boolean, unsigned char** pt_ct,LL enumerate_to_real_key_rank, LL up_to_bound_bool);
void hel_free_param(hel_param_t* param);


//result
hel_real_key_info_t* hel_alloc_real_key_info();
LL hel_init_real_key_info(hel_real_key_info_t* real_key_info);
void hel_free_real_key_info( hel_real_key_info_t* real_key_info);

hel_enum_info_t* hel_alloc_enum_info();
LL hel_init_enum_info(hel_enum_info_t* enum_info);
void hel_free_enum_info( hel_enum_info_t* enum_info);

hel_result_t* hel_alloc_result();
LL hel_init_result(hel_result_t* result);
void hel_free_result(hel_result_t* result);


//accessor
ZZ hel_result_get_estimation_rank(hel_result_t* result);
ZZ hel_result_get_estimation_rank_min(hel_result_t* result);
ZZ hel_result_get_estimation_rank_max(hel_result_t* result);
double hel_result_get_estimation_time(hel_result_t* result);
ZZ hel_result_get_enum_rank(hel_result_t* result);
ZZ hel_result_get_enum_rank_min(hel_result_t* result);
ZZ hel_result_get_enum_rank_max(hel_result_t* result);
double hel_result_get_enum_time_preprocessing(hel_result_t* result);
double hel_result_get_enum_time(hel_result_t* result);
LL hel_result_is_key_found(hel_result_t* result);


#endif
