
#ifndef READ_REALIGNMENT_HPP_
#define READ_REALIGNMENT_HPP_

#include <stdint.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <cstring>
#include <ctype.h>

#include <algorithm>

#include <vector>
#include <map>

#include "../Aln_core/deBGA_index.hpp"
#include "../cpp_lib/get_option_cpp.hpp"
extern "C"
{
#include "../clib/kalloc.h"
#include "../clib/ksw2.h"

#include "../clib/desc.h"
}

#define SEED_OFFSET (K_T - LEN_KMER)
#define LEN_KMER_LEFT (32 - LEN_KMER) //LEN_KMER_LEFT = 32 - LEN_KMER
#define SEED_STEP 5
#define UNI_POS_N_MAX 32

#define SEED_STEP 5
#define MIN_CHAIN_SCORE 20

#define GAP_OPEN 24
#define GAP_EXT 5
#define GAP_OPEN2 42
#define GAP_EXT2 0
#define MATCH_SCORE 2
#define MISMATCH_SCORE 36

//original
//#define GAP_OPEN 16
//#define GAP_EXT 1
//#define GAP_OPEN2 32
//#define GAP_EXT2 0
//#define MATCH_SCORE 2
//#define MISMATCH_SCORE 12


#define ZDROP_SCORE 400
#define BANDWIDTH 500

#define OUTPUT "./aln.sam"

struct MAP_PARA{
	//basic option
	int thread_n;

	//score option
	int match_D;
	int mismatch_D;
	int gap_open_D;
	int gap_ex_D;
	int gap_open2_D;
	int gap_ex2_D;
	int bw;
	int zdrop_D; //for DNA zdrop = 400, 200 for RNA

	int max_use_read;
	//alignment option

	//input option
	char * indexDir;
	char * read_fastq1;
	char * read_fastq2;

	int get_option(int argc, char *argv[]){
		options_list l;
		l.show_command(stderr, argc + 1, argv - 1);
	    l.add_title_string("\n");
	    l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  fc_aln  [Options] <IndexDir> [Read.1.fa][Read.2.fa]\n");
	    l.add_title_string("  Basic:   \n");
	    l.add_title_string("    <IndexDir>      FOLDER   the directory contains index\n");
	    l.add_title_string("    [ReadFiles.fa]  FILES    reads files, FASTQ(A)(or fa.gz/fq.gz) format. nly one file accepted, \n");
	    l.add_title_string("                             for pair end NGS read, read 1 and 2 in a pair should store together.\n");
	    l.add_title_string("                             Using [signal] command to generate this type of file\n");
	    l.add_title_string("    [ori_header.sam]  FILES  Header file of original BAM/CRAM file, using [signal] command or [samtools view -H] to generate it\n");

		//thread number
		l.add_option("thread", 			't', "Number of threads", true, 1); l.set_arg_pointer_back((void *)&thread_n);
		//score
		l.add_option("gap-open1", 		'O', "Gap open penalty 1", true, GAP_OPEN); l.set_arg_pointer_back((void *)&gap_open_D);
			l.l.back().add_help_msg("a k-long gap costs min{O+k*E,P+k*F}..");
		l.add_option("gap-open2", 		'P', "Gap open penalty 2.", true, GAP_OPEN2); l.set_arg_pointer_back((void *)&gap_open2_D);
		l.add_option("gap-extension1", 	'E', "Gap extension penalty 1.", true, GAP_EXT);l.set_arg_pointer_back((void *)&gap_ex_D);
		l.add_option("gap-extension2", 	'F', "Gap extension penalty 2.", true, GAP_EXT2);l.set_arg_pointer_back((void *)&gap_ex2_D);
		l.add_option("match-score", 	'M', "Match score for SW-alignment.", true, MATCH_SCORE);l.set_arg_pointer_back((void *)&match_D);
		l.add_option("mis-score", 		'm', "Mismatch score for SW-alignment.", true, MISMATCH_SCORE);l.set_arg_pointer_back((void *)&mismatch_D);
		l.add_option("zdrop", 			'z', "Z-drop score for splice/non-splice alignment.", true, ZDROP_SCORE);l.set_arg_pointer_back((void *)&zdrop_D);
		l.add_option("band-width", 		'w', "Bandwidth used in chaining and DP-based alignment.", true, BANDWIDTH);l.set_arg_pointer_back((void *)&bw);
		//output
		l.add_option("max_use_read", 	'R', "Max number of read to alignment, used for debug or test PG.", true, MAX_int32t); l.set_arg_pointer_back((void *)&max_use_read);
		if(l.default_option_handler(argc, argv)) {
			l.show_c_value(stderr);
			return 1;
		}
		l.show_c_value(stderr);

		if (argc - optind < 3)
			return l.output_usage();
		xassert((thread_n >= 1) && (thread_n <= 48), "Input error: thread_n cannot be less than 1 or more than 48\n");

		indexDir = strdup(argv[optind]);
		if (indexDir[strlen(indexDir) - 1] != '/') strcat(indexDir, "/");
		read_fastq1 = strdup(argv[optind + 1]);
		read_fastq2 = strdup(argv[optind + 2]);
		return 0;
	}
};

struct SEED_IN_REF
{
	int32_t seed_id;
	int32_t read_begin;
	int32_t read_end;
	int32_t chr_ID;
	int32_t ref_begin;
	int32_t ref_end;
	int32_t cov;

	int32_t read_end_id;
	int32_t direction_mode;
	int32_t used_in_clustering;

	void set(int32_t seed_id, int32_t read_begin,	int32_t read_end,
			int32_t chr_ID, int32_t ref_begin, int32_t ref_end, int32_t cov, int32_t read_end_id,
			int32_t direction_mode){
		this->seed_id = seed_id;
		this->read_begin = read_begin;
		this->read_end = read_end;

		this->chr_ID = chr_ID;
		this->ref_begin = ref_begin;
		this->ref_end = ref_end;
		this->cov = cov;

		this->read_end_id = read_end_id;
		this->direction_mode = direction_mode;
		this->used_in_clustering = false;
	}

	void print(){
		fprintf(stderr,
				"SEED_IN_REF: "
				"read_end_id:[%d], direction_mode:[%d], "
				"read_region:[%d-%d], seed_id:[%d], "
				"ref_region [%d:%d-%d], cov:[%d], diff:[%d]\n",
				read_end_id, direction_mode,
				read_begin, read_end, seed_id,
				chr_ID, ref_begin, ref_end, cov, ref_begin - read_begin);}

	static inline int cmp_by_ref_position(const SEED_IN_REF &a, const SEED_IN_REF &b){
		if(a.read_end_id != b.read_end_id)
			return a.read_end_id < b.read_end_id;
		if(a.direction_mode != b.direction_mode)
			return a.direction_mode < b.direction_mode;
		if(a.chr_ID != b.chr_ID)
			return a.chr_ID < b.chr_ID;
		return a.ref_end < b.ref_end;
	}

};

#define PARAM 2
#define INTRON_PENALTY 2

struct Cluster_Node{

	void init(SEED_IN_REF *sir){
		this->pre_node = NULL;
		this->sir = sir;
		this->sdp_score = sir->cov;;
	}

	void setNewScore(int score, Cluster_Node *pre_node){
		if(score > this->sdp_score){
			this->sdp_score = score;
			this->pre_node = pre_node;
		}
	}
	Cluster_Node *pre_node;
	SEED_IN_REF *sir;
    int sdp_score;

    void show(){
    	fprintf(stderr, "-->[%d] ", sdp_score);
    	sir->print();
    }
};


void print_cigar_list(uint32_t* cigar, int cigar_l);
int cigar_adjust(uint32_t *cigar_l_, uint32_t *cigar, bool add_blank, int delete_small_tail);
void print_cigar_list(uint32_t* cigar, int cigar_l);
void print_cigar_list(std::vector<uint32_t>& c);
void print_X_E_sequence(uint32_t* bam_cigar,int cigar_n, uint8_t* qseq, uint8_t *tseq);

struct Cluster{
	std::vector<Cluster_Node> node_list;
	std::vector<Cluster_Node*> max_chain;

	int max_score_chaining_score;
	//chaining information
	int chrID;
	int ref_position_bg;
	int ref_position_ed;

	int read_end_id;
	int read_position_bg;
	int read_position_ed;
	int direction;
	//ksw-alignmnet
	int alignment_score;
	std::vector<uint32_t> kswCigar;
	uint64_t ksw_final_position_global;
	uint32_t ksw_final_position;


	void clear(){
		node_list.clear();
		max_chain.clear();
		max_score_chaining_score = 0;
		chrID = 0;
		read_position_bg = MAX_int32t;
		read_position_ed = 0;

		read_end_id = -1;
		ref_position_bg = MAX_int32t;
		ref_position_ed = 0;
		direction = 0;
	}

	void generateMAXChain(){
		Cluster_Node *max_socreNode = NULL;
		for(Cluster_Node&n : node_list){
			if(max_score_chaining_score < n.sdp_score){
				max_score_chaining_score = n.sdp_score;
				max_socreNode = &n;
			}
		}
		if(max_socreNode != NULL){
			//set chrID and direction of the chain
			chrID = max_socreNode->sir->chr_ID;
			direction = max_socreNode->sir->direction_mode;
			read_end_id = max_socreNode->sir->read_end_id;
			//store
			max_chain.emplace(max_chain.begin(), max_socreNode);
			while(max_socreNode->pre_node != NULL){
				max_chain.emplace(max_chain.begin(), max_socreNode->pre_node);
				max_socreNode = max_socreNode->pre_node;
			}
		}
		//set the region of the chain
		for(Cluster_Node* cp: max_chain){
			read_position_bg = MIN(read_position_bg, cp->sir->read_begin);
			read_position_ed = MAX(read_position_ed, cp->sir->read_end);

			ref_position_bg = MIN(ref_position_bg, cp->sir->ref_begin);
			ref_position_ed = MAX(ref_position_ed, cp->sir->ref_end);
		}
	}

	static inline int cmp_by_max_score(const Cluster &a, const Cluster &b){
		return a.max_score_chaining_score > b.max_score_chaining_score;
	}

	void show(char endl){
		fprintf(stderr, "Cluster: read_end_id[%d]\t direction [%d]\t, read region:[%d-%d]\t", read_end_id, direction, read_position_bg, read_position_ed);
		fprintf(stderr, "ref region:[%d:%d-%d]\t",  chrID, ref_position_bg, ref_position_ed);
		fprintf(stderr, "max_chaining score:[%d]\t",  max_score_chaining_score);
		if(false){
			fprintf(stderr, "Node list is:{ ");
			for(Cluster_Node & n :node_list)
				n.show();
			fprintf(stderr, "}");
		}
		fprintf(stderr, "%c", endl);
	}

	void showKSW_align_rst(uint8_t *read_bin[2][2], int read_l, deBGA_INDEX * idx){
		fprintf(stderr, "Read_end_id[%d]\t direction [%d]\t, chaining read region:[%d-%d]\t", read_end_id, direction, read_position_bg, read_position_ed);
		fprintf(stderr, "chaining ref region:[%d:%d-%d]\t",  chrID, ref_position_bg, ref_position_ed);
		fprintf(stderr, "alignment_score is %d\t KSW position is [%d;%d] CIGAR is ", alignment_score, chrID, ksw_final_position);
		print_cigar_list(kswCigar);
		std::vector<uint8_t > ref;
		idx->get_refseq_char_bin(ref, ksw_final_position_global, 2*read_l);
		print_X_E_sequence(&(kswCigar[0]), kswCigar.size(), read_bin[read_end_id][direction], &(ref[0]));
	}
};
//
//#define OUT_BAM_CIGAR_STR   "MIDNSHP=XB"
//struct CIGAR_PATH{
//	CIGAR_PATH(char type_, uint16_t size_){
//		switch(type_){
//		case 'M': type = 0; break;
//		case 'I': type = 1; break;
//		case 'D': type = 2; break;
//		case 'N': type = 3; break;
//		case 'S': type = 4; break;
//		case 'H': type = 5; break;
//		case 'P': type = 6; break;
//		case '=': type = 7; break;
//		case 'X': type = 8; break;
//		case 'B': type = 9; break;
//		default : type = 0; break;
//		}
//		size = size_;
//	}
//	CIGAR_PATH(uint32_t bin_cigar){
//		type = (int)((bin_cigar & 0xf));
//		size = (bin_cigar >> 0x4);
//	}
//	CIGAR_PATH(const CIGAR_PATH &cp){
//		type = cp.type;
//		size = cp.size;
//	}
//
//	//return true when merge to the previous node; false when not merge and need to new a node
//	bool try_merge(const CIGAR_PATH &cp){
//		if(cp.size < 0){
//			xassert(cp.type == 2, ""); //deletion
//			if(type == 0){
//				size += cp.size;
//				xassert(size > 0, "");
//				return true;
//			}else if(type == 2){
//				size -= cp.size;
//				xassert(size > 0, "");
//				return true;
//			}else{
//				xassert(0, "");
//			}
//		}else if(type == cp.type || cp.size == 0){
//			size += cp.size;
//			return true;
//		}
//		return false;
//	}
//
//	uint8_t type;
//	int16_t size;
//};

struct KSW_ALN_handler{

#define KSW_ALN_left_extend 0
#define KSW_ALN_right_extend 1
#define KSW_ALN_end_to_end 2
//#define KSW_ALN_left_extend_forward_cigar 3

	//read info
	uint8_t *read_str;
	//uint32_t read_cigar[100];
	//uint32_t n_cigar;
	uint32_t total_q_len;

	//sequence buff
	uint8_t *qseq;
	uint8_t *tseq;
	uint8_t *qseq_rev;
	uint32_t qlen;
	uint32_t tlen;
	uint8_t type;

	//result
	bool is_simple_aln;
	uint32_t simple_NM;
	ksw_extz_t ez;

	//kmalloc
	void *km;

	//mapping options
	int8_t mata_D[25];

	int8_t match_D;
	int8_t mismatch_D;
	int8_t gap_open_D;
	int8_t gap_ex_D;
	int8_t gap_open2_D;
	int8_t gap_ex2_D;
	uint16_t zdrop_D; //for DNA zdrop = 400, 200 for RNA
	int bandwith;

	int flag;
	std::vector<uint32_t> cigar;

	void showCigar(){
		fprintf(stderr, "Cigar is ");
		print_cigar_list(&(cigar[0]), cigar.size());
		fprintf(stderr, "\n");
	}
	int adjustCIGAR(){
		uint32_t n_cigar = cigar.size();
		int adj_size = cigar_adjust(&n_cigar, &(cigar[0]), false, 15);
		cigar.resize(n_cigar);
		return adj_size;
	}

	deBGA_INDEX * idx;

	void copy_option(MAP_PARA 	*o);
	void ksw_gen_mat_D();
	void init(MAP_PARA 	*o, deBGA_INDEX * idx_);
	void free_memory();
	void setRead(uint8_t *read_str_);
	void alignment(bool printLog, int ref_chrID,  int read_st_, int read_ed_, int ref_st_, int ref_ed_, int type_);
	void align_non_splice();
	int get_misMatch(int read_st_, int read_ed_, int ref_st_, int ref_ed_);

};

//used to store temporary alignment results
struct MAX_IDX_OUTPUT{
	//alignment scores
	uint32_t align_score;
	uint32_t chain_score;

	//position in read
	uint32_t max_index;
	uint32_t read_bg;

	//mapq
	uint8_t mapq;

	//mate information
	bool has_mate;
	uint32_t mate_chrID;
	uint32_t mate_ref_bg;
	//position in reference
	uint32_t chrID;
	uint32_t ref_bg;
	int direction;

	//alignment cigar
	std::vector<uint32_t> cigar;

	//index in results
	int rst_idx;

	static inline int cmp_chain_score(const MAX_IDX_OUTPUT &a,const MAX_IDX_OUTPUT&b)	{
		if(a.chain_score != b.chain_score)	return a.chain_score < b.chain_score;
		else									return a.max_index > b.max_index;
	}

	static inline int cmp_align_score(const MAX_IDX_OUTPUT &a,const MAX_IDX_OUTPUT &b)	{
		if(a.align_score != b.align_score)	return a.align_score < b.align_score;
		else									return a.max_index > b.max_index;
	}

	void print(FILE * f, char endl){ fprintf(f, "chain_score:[%d], align_score:[%d], read_bg:[%d], chrID:[%d], ref_bg:[%d], %c", chain_score, align_score, read_bg, chrID, ref_bg, endl); }

};

//used for alignment of single end read
#define MAX_READ_LEN 1600
#define MAX_OUTPUT_NUMBER 6
struct Single_Read_handler{

	uint random_seed = 0;
	//result
	std::vector<MAX_IDX_OUTPUT> result;

	MAX_IDX_OUTPUT * primary_result;
	MAX_IDX_OUTPUT * secondary_result;

	// step 1:
	void read_register(kseq_t * c_read_, int read_end_id){
		c_read = c_read_;
		read_l = c_read->seq.l;
		read_seq_char = c_read->seq.s;
		binary_read_2_bit();
		this->read_end_id = read_end_id;
	}

	// step2 :aligning read into reference and store all result in [result]，total result number stored in [result_num], and
	// primary ID and second ID stored in 	[int primary_index;] and [int secondary_index;]
	//private function and variable:
	uint64_t read_l;
	char * read_seq_char;
	kseq_t * c_read;
	uint8_t bin_read[2][MAX_READ_LEN];
	int read_end_id;

	//buff:
	std::vector<MEM_Within_UNITIG> MEM_within_unitig_rst_l;

	//final result
	int postion;
	int chrID;
	int MAPQ;

	//output buff
	char sam_string[2048];
	kstring_t s;

	//ori_rst[0] is the original result of read 1 and ori_rst[1] is the second read
	void search_and_chain_seed_within_UNITIG(int direction_mode, std::vector<Chain_Within_Unitig> &CWU_l, deBGA_INDEX *idx);
	void binary_read_2_bit();
};

struct PE_score{
	//results part
	int max_same;
	int max_score;
	//properly mated
	bool read_pair_is_proper_mated;
	int cur_isize; // current isize

	//bool: this value indicated whether the pan-genome gain better results compared with original results
	//bool: this value indicated whether read pair nearly full-match-aligned to the original genome or pan-genome
	//bool read_pair_full_match; // when all reads in a pair has 5 mis-matches or other mutations, this value will be set to true.

	MAX_IDX_OUTPUT *max_1;
	MAX_IDX_OUTPUT *max_2;

	//bam/cram status
	int max_isize; //the max ISZIE the that can be accepted, 99% of reads have ISIZEs shorter than that value
	int min_isize; //the min ISZIE the that can be accepted, 99% of reads have ISIZEs longer than that value
	int normal_read_len;
	int min_filter_score;

	void init(int max_isize_, int min_isize_, int normal_read_len_, int min_filter_score_){
		max_isize = max_isize_ + 200;
		min_isize = min_isize_ - 200;
		min_isize = MAX(0, min_isize);
		normal_read_len = normal_read_len_;
		min_filter_score = min_filter_score_;
		clear();
	}

	void clear(){
		max_same = 1;
		max_score = 0;
		max_1 = NULL;
		max_2 = NULL;
		cur_isize = 0;
		read_pair_is_proper_mated = false;
		//read_pair_full_match = false;
	}

	void read_get_best_pairing_results(Single_Read_handler &SE_h, uint32_t *random_seed_p){
		clear();
		//for both original
		int se0_rst_num = SE_h.result.size(); int se1_rst_num = SE_h.result.size();
		//when using both end
		for(int i = 0; i < se0_rst_num; i++){
			for(int j = 0; j < se1_rst_num; j++)
				store_pair_end_score(&SE_h.result[i], &SE_h.result[j], random_seed_p);
		}
	}

	void set_primary_secondary_mate(Single_Read_handler &SE_h){
		for(int i = 0; i < 2; i++){
			MAX_IDX_OUTPUT *c_rst = get_cur_output(i);
			if(c_rst == NULL)
				continue;
			//set primary results
			SE_h.primary_result = c_rst;
			//set secondary
			SE_h.secondary_result = NULL;
			if(SE_h.result.size() > 1){
				SE_h.secondary_result = (c_rst->rst_idx == 0)?&(SE_h.result[1]):&(SE_h.result[0]);
			}
			//set mate
			MAX_IDX_OUTPUT *mate_rst = get_cur_output(1 - i);
			if(mate_rst != NULL){
				c_rst->has_mate = true;
				c_rst->mate_chrID = mate_rst->chrID;
				c_rst->mate_ref_bg = mate_rst->ref_bg;
				//when one end is original and mate read is new alignment, set the original has same SV info as new alignment
			}
			else{
				c_rst->has_mate = false;
			}
		}
	}

private:

	MAX_IDX_OUTPUT *get_cur_output(int read_ID){ return (read_ID == 0)?max_1:max_2; }

	//a sub-function used in "store_pair_end_score"
	void set_score(int new_score, MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2, int ISIZE, uint32_t *random_seed_p){
		bool store_new_best = true;
		if(new_score > max_score) max_same = 1;
		else if(new_score == max_score){
			max_same ++;
			if(rand_r(random_seed_p) % max_same != 0) store_new_best = false;
		}
		if(store_new_best){
			max_1 = se1; max_2 = se2; max_score = new_score; cur_isize = ISIZE; read_pair_is_proper_mated = (cur_isize > 0);
		}
	}

	//return true isize (> 0); when not proper paired, return 0;
	int get_isize(int p1, int p2, int d1, int d2){
		if(d1 == d2) return 0;
		int isize = normal_read_len + ((d1 == FORWARD)?(p2 - p1):(p1 - p2));
		return (isize < max_isize && isize > min_isize)?isize: 0;
	}
	//it is a sub-function for "store_pair_end_score"
	//return isize (>0) when reads are properly mapped; otherwise return 0
	int read_is_proper_mated(MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2){
		if(se1 == NULL || se2 == NULL || se1->chrID != se2->chrID){ return 0; }
		return 0;
	}

	//calculate the final score of single end results part se1 and se2
	int store_pair_end_score(MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2, uint32_t *random_seed_p){
		//detect whether reads is nearby
		int ISIZE =	read_is_proper_mated(se1, se2);
		//detect whether one of reads is new alignment
		int basic_score = ((se1 == NULL)?0:se1->align_score) + ((se2 == NULL)?0:se2->align_score); //basic scores
		int final_score = basic_score +
				+ ((ISIZE > 0)?(0):(-60)); //not nearby penalty

		//if(read_pair_full_match == false && final_score > FINAL_SCORE_CUTOFF && ISIZE > 0){ read_pair_full_match = true; }

		if(final_score >= max_score)
			set_score(final_score, se1, se2, ISIZE, random_seed_p);
		return 0;
	}
};

#define K_T 25 		//kmer length used for building deBGA index

//data used for each ALIGN-PROCESSIONG thread
struct Classify_buff_pool
{
	//pair end handler:
	PE_score ps;
	//single end handler
	Single_Read_handler SEH[2];


	std::vector<Cluster> cluster_l[2];//graph node
	//data used for both end
	deBGA_INDEX *idx;// =  SE_h[0].idx;
	std::vector<Chain_Within_Unitig> CWU_l;
	std::vector<SEED_IN_REF> SIR_l;

	random_data rand_buff;
	uint random_seed = 0;
	kstring_t s;
	KSW_ALN_handler kswh;

	void init(MAP_PARA 	*o, deBGA_INDEX * idx_){
		char * c_statebuf = (char *)xcalloc(128, 1);
		initstate_r(rand_r(&random_seed), c_statebuf, 128, &rand_buff);
		kswh.init(o, idx_);
		s.s = (char *)xcalloc(0xfffff,1);//2^20=1M
		s.m = 0xfffff;
		idx = idx_;
	}

};

//data used for PIPELINE each threads
struct CLASSIFY_THREAD_DATA
{
	kseq_t	*	seqs1; // first  read in a pair
	kseq_t	*	seqs2; // second read in a pair
	//OUTPUT_RST * rst; 	//mapping result for read1 and read2
	int 		readNum;
	void 	*	share_data_pointer;// register for the shared data
	std::vector<std::string> result1;
	std::vector<std::string> result2;
};

//data SHARED by all threads
struct CLASSIFY_SHARE_DATA
{
	//shared
	deBGA_INDEX * idx = NULL;//deBGA index
	kstream_t 	*_fp1 = NULL; //fastq files 1
	kstream_t 	*_fp2 = NULL; //fastq files 1

	//deBUG:
	FILE * debug_1 = NULL;
	FILE * debug_2 = NULL;

	MAP_PARA 	*o; 		//mapping option
	Classify_buff_pool 	*buff = NULL;//data used for each classify thread
	CLASSIFY_THREAD_DATA *data = NULL;//for each pipeline thread
};

struct deCOY_CLASSIFY_MAIN{
	CLASSIFY_SHARE_DATA *share = NULL;
	void init_run(int argc, char *argv[]);
private:
	static void *classify_pipeline(void *shared, int step, int tid, void *_data);									//pipeline
	static void *classify_pipeline_single_thread(void *shared);
	static int load_reads(kstream_t *_fp1, kstream_t *_fp2, kseq_t *_seqs1, kseq_t *_seqs2, int n_needed, MAP_PARA *o);
	static void inline worker_for(void *_data, long data_index, int thread_index); 									//function in step 2
	static void output_results( int readNum, std::vector<std::string> &result1, std::vector<std::string> &result2);
};

#endif /* READ_REALIGNMENT_HPP_ */
