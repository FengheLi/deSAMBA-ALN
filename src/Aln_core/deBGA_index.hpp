#ifndef DEBGA_INDEX_HPP_
#define DEBGA_INDEX_HPP_

#include <stdint.h>
#include <stdlib.h>
#include <string>
extern "C"
{
#include "../clib/desc.h"
#include "../clib/utils.h"
}

#include <vector>

#define LEN_KMER 25 //kmer length used for search
#define START_POS_REF 0

#define UNITIG_MAGIC_NUM 2048

//MEM
struct MEM_Within_UNITIG
{
	uint64_t uid;
	uint32_t seed_id;
	uint32_t read_pos;
	uint32_t uni_pos_off;
	uint32_t length;
	uint32_t pos_n;

	static int cmp(const void * a, const void * b)
	{
		MEM_Within_UNITIG* sm1 = (MEM_Within_UNITIG *)a;
		MEM_Within_UNITIG* sm2 = (MEM_Within_UNITIG *)b;

	    if(sm1->uid > sm2->uid)
	        return 1;
	    if(sm1->uid < sm2->uid)
	        return -1;
	    else
	    {
	        //according to read_pos and uni_pos_off detected if there is an inversion
	        if(sm1->read_pos > sm2->read_pos)
	            return 1;
	        if(sm1->read_pos < sm2->read_pos)
	            return -1;
	        else    return 0;
	    }
	}

	void print(){ fprintf(stderr, "vertex_MEM: uid:[%ld], seed_id:[%d], read_pos:[%d], uni_pos_off:[%d], length:[%d], pos_n:[%d] \n",
				uid, seed_id, read_pos, uni_pos_off, length, pos_n);}

};

//MEM joint region in uinpath
struct Chain_Within_Unitig
{

	uint64_t uid;
	uint32_t read_offset;
	uint32_t unitig_offset;
	uint32_t length_in_read; //length in read
	uint32_t length_in_ref;// length in reference
	uint32_t ref_pos_n;
	uint32_t cov;

	uint32_t read_end_id;
	uint32_t direction_mode;

	void print(){
		fprintf(stderr, "Chain_Within_Unitig: read_end_id[%d], direction_mode[%d], uid:[%ld], read_pos:[%d], uni_pos_off:[%d],"
				" length_in_read:[%d], length_in_ref:[%d], reference_position_num:[%d], cov:[%d] \n",
				read_end_id, direction_mode, uid, read_offset, unitig_offset,
				length_in_read, length_in_ref, ref_pos_n, cov);}
};

struct deBGA_INDEX{

	uint64_t* buffer_ref_seq = NULL;//original reference
	uint64_t* buffer_seq = NULL;// UNITIG sequence

	uint64_t* buffer_seqf = NULL;//start offset of [UINTIG] in [buffer_seq]
	uint64_t* buffer_off_g = NULL;//start offset of [KMER] in [buffer_seq]
	uint64_t* buffer_p = NULL;//start offset of [UINTIG] in [buffer_ref_seq]
	uint64_t* buffer_pp = NULL;//start index of [UINTIG] in [buffer_p], used to search offset of [UINTIG] in [buffer_ref_seq]
	uint64_t* buffer_hash_g = NULL;//index for the first 14 base pair, used to search kmer in the hash table
	uint32_t* buffer_kmer_g = NULL;//used to store other part of a kmer(except the first 14 byte)

	uint32_t chr_search_index_size = 0;
	uint32_t * chr_search_index = NULL;

	uint64_t reference_len = 0;
	std::vector<uint64_t>chr_end_position;
	std::vector<std::string>chr_names;

	uint64_t result_ref_seq = 0;
	uint64_t result_seq = 0;
	uint64_t result_seqf = 0;
	uint64_t result_p = 0;
	uint64_t result_pp = 0;
	uint64_t result_pu = 0;
	uint64_t result_hash_g = 0;
	uint64_t result_kmer_g = 0;
	uint64_t result_off_g = 0;
	uint64_t result_ref_g = 0;

	int load_index_file(char *index_dir);

	//search kmer in the index, the search result "index of kmer" stored in range, range[0] is the start index, and range[1] is the end index
	//return -1 when search failed, otherwise return 0
	bool search_kmer(int len_k, uint64_t kmer, int64_t range[2], int8_t seed_offset);
	//function: given an index of kmer, search the MEM within a UNITIG
	//return 0 when successfully stored data, -1 when failed
	int MEM_extern_within_UNITIG_for_KMER(uint64_t kmer_index, std::vector<MEM_Within_UNITIG>& vertexm_v, uint64_t *read_bit, uint32_t read_off, uint32_t read_length, int len_k, uint32_t & max_right_i);
	//function: merge all MEMs within a unipath, the result stored in vertexu_v
	void merge_seed_in_unipath(std::vector<MEM_Within_UNITIG>& vertexm_v, std::vector<Chain_Within_Unitig>& vertexu_v);

	void get_refseq(uint8_t *ref, uint32_t len, uint32_t start);
	void get_refseq_char_bin(std::vector<uint8_t> &ref, uint32_t start, uint32_t len);
	void get_refseq_char(std::vector<char> &ref, uint32_t start, uint32_t len);
	void get_ref_jump(uint8_t *ref, uint32_t start1, uint32_t len1, uint32_t start2, uint32_t len2);
	void get_ref_onebyone(uint8_t *ref, uint32_t start, uint32_t len, uint32_t pre_pos);

	void building_chr_index();
	//get chr ID and begin
	int get_chromosome_ID(uint64_t position);
	void global_position_to_chrID_and_position_in_chr(uint64_t position_global, int & chrID, int &position_in_chr);
	void chrID_and_position_in_chr_to_global_position(uint64_t &position_global, int chrID, int position_in_chr);

	void building_bam_header();

	void free_memory();

};
uint64_t load_data_from_file(const char * path_name, const char * fn, void ** data, bool useCalloc, int additional_byte);
int build_deBGA_index(int argc, char *argv[]);

#endif /* DEBGA_INDEX_HPP_ */
