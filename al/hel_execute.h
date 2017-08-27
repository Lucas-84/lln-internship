#ifndef HEL_EXECUTE_H
#define HEL_EXECUTE_H

#include "hel_enum.h"

LL hel_execute_procedure(hel_param_t* param, hel_result_t* result);

hel_result_t* hel_execute_rank(LL merge_value, LL nb_bins, double** score_mat_init,LL* key_init);

hel_result_t* hel_execute_enum(LL merge_value, LL nb_bins, double** score_mat_init, LL* key_init, ZZ bound_start, ZZ bound_end,LL test_key_bool, unsigned char** pt_ct, LL enumerate_to_real_key_rank_bool, LL up_to_bound_bool);

#endif
