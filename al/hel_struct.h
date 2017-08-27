#ifndef HEL_STRUCT_H
#define HEL_STRUCT_H

#include <NTL/ZZXFactoring.h>
#include <NTL/RR.h>
#include "hel_util.h"
using namespace std;
using namespace NTL;
typedef long long LL;

struct hel_int_list_t{

	LL val;
	struct hel_int_list_t* next;

};

typedef struct hel_int_list_t hel_int_list_t;



typedef enum{

	RANK,
	ENUM

}hel_algo_mode_t;

typedef struct{

	LL updated_nb_subkey; //number of subkeys lists after merging
	LL* updated_nb_key_value; //number of subkeys value per lists after merging
 	LL merge_value;

}hel_preprocessing_t;



typedef struct{

	LL nb_bins;
	RR* width;
	ZZX* hists;

	LL hists_alloc_bool;

}hel_histo_t;

typedef struct{

	double** log_probas;
	LL* real_key;
	double* log_probas_real_key;
	double shift;

	LL log_probas_alloc_bool;
	LL real_key_alloc_bool;
	LL log_probas_real_key_alloc_bool;

}hel_data_t;

typedef struct{

	LL** binary_hists;
	LL* binary_hists_size;
	LL*** key_list_bin;
	hel_int_list_t*** key_list_bin2;
	LL** index_list;
	LL* convolution_order;
	LL** key_factorization;
	LL bin_to_start;
	LL bin_to_end;
	ZZ* bound_start;
	ZZ* bound_end;
	LL test_key_boolean;
	unsigned char** pt_ct;
	LL enumerate_to_real_key_rank; //if set to 1, enumeration is launched only if this rank is less than the one specified in "bound"
	LL* real_key;
	LL up_to_bound_bool;

	LL binary_hists_alloc_bool;
	LL binary_hists_size_alloc_bool;
	LL key_list_bin_alloc_bool;
	LL key_list_bin2_alloc_bool;
	LL index_list_alloc_bool;
	LL convolution_order_alloc_bool;
	LL key_factorization_alloc_bool;
	LL pt_ct_alloc_bool;
	LL real_key_alloc_bool;

}hel_enum_input_t;


typedef struct{

	LL bin_real_key;  //bin where the real key is found
	ZZ* bound_real_key; //associated value of key number

	LL bin_bound_min; //bin of the min bound of the real key
	ZZ* bound_min; //associated value of key number

	LL bin_bound_max; //bin of the max bound of the real key
	ZZ* bound_max; //associated value of key number

	double rank_estimation_time;

}hel_real_key_info_t;


typedef struct{

	LL bin_enum_start; //bin where the enumeration starts according to the user's provided bound
	ZZ* nb_key_enum_start; //associated sumed number of key from the last bin to this bin

	LL bin_enum_end; //bin where the enumeration ends according to the user's provided bound
	ZZ* nb_key_enum_end; //associated value of key number

	LL bin_found_key; //bin where the real key is found (if found)
	ZZ* nb_key_enum_found_key; //associated sumed number of key from the last bin to this bin

	LL bin_found_key_bound_min; //bin associated to the min bound of the real key (if found)
	ZZ* nb_key_enum_found_key_bound_min; //associated sumed number of key from the last bin to this bin

	LL bin_found_key_bound_max; //bin associated to the max bound of the real key (if found)
	ZZ* nb_key_enum_found_key_bound_max; //associated sumed number of key from the last bin to this bin

	LL bound_bin_enum_start; //bin associated to the min bound where the enumeration starts according to the user's provided bound
	ZZ* bound_nb_key_enum_start; //associated sumed number of key from the last bin to this bin

	LL bound_bin_enum_end; //bin associated to the max bound where the enumeration ends according to the user's provided bound
	ZZ* bound_nb_key_enum_end; //associated sumed number of key from the last bin to this bin

	double preprocessing_time;
	double enum_time;

    LL found_key_boolean;


}hel_enum_info_t;

typedef struct{

	hel_enum_info_t* enum_info;
	hel_real_key_info_t* real_key_info;

}hel_result_t;


typedef struct{

	hel_algo_mode_t algo_mode;
	hel_preprocessing_t* preprocessing;
	hel_histo_t* histo;
	hel_data_t* data;
	hel_enum_input_t* enum_input;

}hel_param_t;

#endif
