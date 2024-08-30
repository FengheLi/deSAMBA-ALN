
#include <string.h>

#include <getopt.h>
#include <algorithm>
#include <stdio.h>
#include <sys/time.h>

#include "zlib.h"
#include "align.hpp"
#include "../cpp_lib/cpp_utils.hpp"

extern "C"{
#include "../clib/kthread.h"
#include "../clib/utils.h"
#include "../clib/desc.h"
}

#include <math.h>

//**************************************class: CLASSIFY_MAIN********************************/
#define PIPELINE_T_NUM 3//one for reading; one for classifying; one for writing
#define STEP_NUM PIPELINE_T_NUM
#define N_NEEDED 2000000 //2M read pair per time
extern void align_read_pair(kseq_t * read1,  kseq_t * read2, Classify_buff_pool * buff, std::string & rst_1,std::string & rst_2, FILE * debug_1, FILE * debug_2);

void outputSAM_header(FILE * output, int argc, char *argv[],
		std::vector<uint64_t> &chr_end_position,
		std::vector<std::string> &chr_names){
	fprintf(stdout, "@HD\tVN:1.6\tSO:unsorted\n");
	if(!chr_end_position.empty())
		fprintf(stdout, "@SQ\tSN:%s\tLN:%u\n", chr_names[0].c_str(), (uint32_t )(chr_end_position[0]));
	for(uint i = 1; i < chr_end_position.size(); i++)
		fprintf(stdout, "@SQ\tSN:%s\tLN:%u\n", chr_names[i].c_str(), (uint32_t )(chr_end_position[i] - chr_end_position[i - 1]));
	fprintf(stdout, "@PG\tID:%s\tPN:%s\tVN:%s\tCL: %s ", PACKAGE_NAME, PACKAGE_NAME, PACKAGE_VERSION, PACKAGE_NAME);
	for(int i = 0; i < argc; i++)
		fprintf(stdout, "%s ", argv[i]);
	fprintf(stdout, "\n");
}

void deCOY_CLASSIFY_MAIN::init_run(int argc, char *argv[]){

	//share data
	share = (CLASSIFY_SHARE_DATA*)xcalloc(1, sizeof(CLASSIFY_SHARE_DATA));

	//loading parameters
	share->o = (MAP_PARA *)xcalloc(1, sizeof(MAP_PARA));//mapping option
	if(share->o->get_option(argc, argv) != 0)
		return;
	//loading index
	share->idx = (deBGA_INDEX *)xcalloc(1, sizeof(deBGA_INDEX));
	//original header:
	share->idx->load_index_file(share->o->indexDir);

	//get pipeline threads data
	share->data = (CLASSIFY_THREAD_DATA *)xcalloc(PIPELINE_T_NUM, sizeof(CLASSIFY_THREAD_DATA));
	for(int i = 0; i < PIPELINE_T_NUM; i++)
	{
		share->data[i].seqs1	= (kseq_t *)xcalloc(N_NEEDED, sizeof(kseq_t));
		share->data[i].seqs2	= (kseq_t *)xcalloc(N_NEEDED, sizeof(kseq_t));

		share->data[i].result1.resize(N_NEEDED);
		share->data[i].result2.resize(N_NEEDED);

		share->data[i].share_data_pointer = share;//anti-register for share data
	}

	//get alignment thread data
	share->buff = (Classify_buff_pool*)xcalloc(share->o->thread_n, sizeof(Classify_buff_pool));
	for(int i = 0; i < share->o->thread_n; i++){
		share->buff[i].init(share->o, share->idx);
		share->buff[i].ps.init(0,0,0, 0);
	}

	//running alignment
	fprintf(stderr,"Start classify\n");
	double cpu_time = cputime();
	struct timeval start;
	gettimeofday(&start, NULL);
	//open SAM/BAM files
	//fprintf(stderr, "Open SAM/BAM file [%s]\n", share->o->read_fastq1);
	gzFile fp1 = xzopen(share->o->read_fastq1, "rb");
	gzFile fp2 = xzopen(share->o->read_fastq2, "rb");
	//output bam file
	share->_fp1 = ks_init(fp1);
	share->_fp2 = ks_init(fp2);
	//
	share->debug_1 = fopen("/media/fenghe/Data/window_linux_share/eclipse-workspace/panSVR_aln/Release/demo_r1.log", "w");
	share->debug_2 = fopen("/media/fenghe/Data/window_linux_share/eclipse-workspace/panSVR_aln/Release/demo_r2.log", "w");

	//fprintf(stderr, "Output sam header\n");
	outputSAM_header(stdout, argc, argv, share->idx->chr_end_position, share->idx->chr_names);
	fprintf(stderr, "Processing file: [%s].\n", share->o->read_fastq1);

	bool run_in_single_thread_mode = false;
	if(run_in_single_thread_mode){
		classify_pipeline_single_thread(share);
	}else{
		kt_pipeline(PIPELINE_T_NUM, classify_pipeline, share, STEP_NUM);
	}
	//close FASTA files
	gzclose(fp1);
	gzclose(fp2);

	if(share->debug_1 != NULL) fclose(share->debug_1);
	if(share->debug_2 != NULL) fclose(share->debug_2);

	fprintf(stderr, "Classify CPU: %.3f sec\n", cputime() - cpu_time);
}

#define MAX_read_size 100000000//100M
void *deCOY_CLASSIFY_MAIN::classify_pipeline(void *shared, int step, int tid, void *_data) {
	CLASSIFY_SHARE_DATA * s = (CLASSIFY_SHARE_DATA*) shared;
	CLASSIFY_THREAD_DATA *d = (CLASSIFY_THREAD_DATA *) s->data;
	//step0: read read data from files; step1: process; step2: output result
	if 		(step == 0)	{ 	if((s->data[tid].readNum = load_reads(s->_fp1, s->_fp2, s->data[tid].seqs1, s->data[tid].seqs2, N_NEEDED, s->o))) 	return (void *)1; }
	else if (step == 1)	{	kt_for(s->o->thread_n, worker_for, s->data + tid, s->data[tid].readNum);  							return (void *)1; }
	else if (step == 2)	{	output_results(	s->data[tid].readNum, d->result1, d->result2);    	return (void *)1;}
	return 0;
}

#define MAX_read_size 100000000//1000M
#define N_NEEDED_simple 30000000//30M reads
void *deCOY_CLASSIFY_MAIN::classify_pipeline_single_thread(void *shared) {
//	CLASSIFY_SHARE_DATA * s = (CLASSIFY_SHARE_DATA*) shared;
//	//step0: read read data from files; step1: process; step2: output result
//	s->data[0].readNum = load_reads(s->_fp1, s->_fp2, s->data[0].seqs1, s->data[0].seqs2, N_NEEDED_simple, s->o);
//	//kt_for(s->o->thread_n, worker_for, s->data + 0, s->data[0].readNum);
//	CLASSIFY_THREAD_DATA *d = (CLASSIFY_THREAD_DATA *) s->data;
//	//align_read_pair(d->seqs1 + data_index, d->seqs2 + data_index, s->buff + thread_index, d->b1 + data_index, d->b2 + data_index, d->ori_b1 + data_index, d->ori_b2 + data_index);
//	for(int data_index = 0; data_index < s->data[0].readNum;data_index++ ){
//		align_read_pair(d->seqs1 + data_index, d->seqs2 + data_index, s->buff + 0, d->result1[data_index], d->result2[data_index]);
//	}
//	output_results(s->data[0].readNum, d->result1, d->result2);

	return 0;
}

//function in step 1
static int load_read_number = 0;
int deCOY_CLASSIFY_MAIN::load_reads(kstream_t *_fp1, kstream_t *_fp2, kseq_t *_seqs1, kseq_t *_seqs2, int n_needed, MAP_PARA *o)
{
	kseq_t *temp1 = _seqs1;
	kseq_t *temp2 = _seqs2;
	int i, rst1 = 0, rst2 = 0, total_length = 0;
	for( i = 0; i < n_needed &&  total_length < MAX_read_size && load_read_number < o->max_use_read; ++i)
	{
		temp1[i].f = _fp1; //load two reads from same file
		temp2[i].f = _fp2;
		rst1 = kseq_read(temp1 + i);
		rst2 = kseq_read(temp2 + i);
		if(rst1 < 0 || rst2 < 0)
			break;
		total_length += (temp1[i].seq.l + temp2[i].seq.l);
		load_read_number++;
	}
	return i;
}

//function in step 2
void inline deCOY_CLASSIFY_MAIN::worker_for(void *_data, long data_index, int thread_index)
{ // kt_for() callback
	CLASSIFY_THREAD_DATA *d = (CLASSIFY_THREAD_DATA *) _data;
	CLASSIFY_SHARE_DATA  *s = (CLASSIFY_SHARE_DATA  *) (d->share_data_pointer);

	align_read_pair(d->seqs1 + data_index, d->seqs2 + data_index, s->buff + thread_index, d->result1[data_index], d->result2[data_index], s->debug_1,s->debug_2);//, d->b1 + data_index, d->b2 + data_index, d->ori_b1 + data_index, d->ori_b2 + data_index);
}

//GLOBAL VALUABLE
static int read_block_count = 0;
void deCOY_CLASSIFY_MAIN::output_results( int readNum, std::vector<std::string> &result1, std::vector<std::string> &result2){
		fprintf(stderr, "Processing %d reads, at block ID %d\n", readNum, read_block_count++);
		for(int i = 0; i < readNum; i++){
			//output the SAM results
			if(!result1[i].empty()){
				fprintf(stdout, "%s\n", result1[i].c_str());
			}
			if(!result2[i].empty()){
				fprintf(stdout, "%s\n", result2[i].c_str());
			}
		}

}

//**************************************class: single_end_handler********************************/

uint8_t charToDna5n[] =
{
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,//4 the last but one
    /*   		 A     C           G                    N */
    /*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//'Z'
    /*          	      T */
    /*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
    /*   		 a     c           g                    n */
    /* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*          	      t */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

uint64_t inline getKmer(uint32_t read_off, uint64_t *read_bit){
	uint32_t index_word = read_off >> 5;
	uint32_t index_in_word = (read_off & 0X1f);
	uint64_t full_kmer = (((read_bit[index_word]) << ((index_in_word) << 1)) | ((index_in_word == 0)?0:(read_bit[index_word+ 1] >> ((32 - index_in_word) << 1))));
	return full_kmer >> (LEN_KMER_LEFT << 1);
	//return (((read_bit[index_word] & kmerMask[index_in_word]) << ((index_in_word - LEN_KMER_LEFT) << 1)) | (read_bit[index_word+ 1] >> ((64 - LEN_KMER - index_in_word) << 1)));
}
//
//static int sort_output(Graph_handler &g, deBGA_INDEX * idx, std::vector<MAX_IDX_OUTPUT> &rst, int direction, uint32_t *random_seed){
//	//output the best
//	if(g.seed_size == 0) return 0;
//
//	//get max score
//	g.max_index = MAX_uint32_t; g.max_distance = 0;
//	g.same_top_max_distance_id_list.clear();
//	g.same_top_max_distance_id_list.emplace_back(g.max_index);
//	int i = g.seed_size - 1;
//   	for(; i >= 0; i--){
//   		if(g.dist_path[i].already_used)
//   			continue;
//   		float c_dist = g.dist_path[i].dist;
//		if (g.max_distance < c_dist){
//			g.max_distance = c_dist;
//			g.max_index = i;
//			g.same_top_max_distance_id_list.clear();
//			g.same_top_max_distance_id_list.emplace_back(i);
//		}
//		else if(g.max_distance == c_dist){
//			g.same_top_max_distance_id_list.emplace_back(i);
//		}
//   	}
//
////   	if(g.max_distance == 127){
////   		printf(" ");
////   	}
//   	if(g.max_index == MAX_uint32_t)//without reult
//   		return 0;
//
//   	int node_num_already_used = 0;
//   	int node_num_not_used = 0;
//   	//get chain for that score
//   	uint32_t same_size = g.same_top_max_distance_id_list.size();
//   	if(same_size > 1){ //random select a new chain
//   		g.max_index = g.same_top_max_distance_id_list[rand_r(random_seed) % same_size];
//   	}
//
//	int first_node = g.max_index;
//	int original_first_node = first_node;
//
//	//debug code::
//	//fprintf(stderr, "\n\ng.max_index: [%d], g.max_distance:[%f]\n", g.max_index, g.max_distance);
//	for(; first_node != -1; ){
//		//debug code::
//		//g.vertexArr_[0][first_node].print();
//		if(g.dist_path[first_node].already_used == true)
//			node_num_already_used++;
//		else
//			node_num_not_used++;
//		g.dist_path[first_node].already_used = true;
//		int next_node = g.dist_path[first_node].pre_node;
//		if(next_node == -1)
//			break;
//		first_node = next_node;
//	}
//	int original_final_node = first_node;
//
//	//total node number within region is far bigger then all used node in a chain: treat as STR or VNTR
//	if(original_first_node - original_final_node > ((node_num_not_used + node_num_already_used + 5) << 1)){
//		for(int node_ID = original_final_node; node_ID < original_first_node; node_ID ++)
//			g.dist_path[node_ID].already_used = true;
//	}
//
////	if(node_num_already_used >= node_num_not_used){
////		return sort_output(g, idx, rst, direction, random_seed);
////	}
//
//	//debug code::
//	//fprintf(stderr, "read_begin: [%d], ref_begin:[%d]\n", g.vertexArr_[0][first_node].read_begin, g.vertexArr_[0][first_node].ref_begin);
//	//set result
//	int ref_begin = g.SIR_l_BC[0][first_node].ref_begin;
//	int chr_ID = idx->get_chromosome_ID(ref_begin);
//	rst.emplace_back();
//	rst.back().direction = direction;
//	rst.back().max_index = g.max_index;
//	rst.back().chain_score = g.max_distance;
//   	rst.back().read_bg = g.SIR_l_BC[0][first_node].read_begin;
//   	rst.back().chrID = chr_ID;
//   	if(chr_ID == 0)
//   		rst.back().ref_bg = ref_begin;//end of last chr is the begin of this chr
//   	else
//   		rst.back().ref_bg = ref_begin - idx->chr_end_position[chr_ID - 1];//end of last chr is the begin of this chr
//  	rst.back().ref_bg = ref_begin - idx->chr_end_position[chr_ID - 1];//end of last chr is the begin of this chr
//	return 1;
//}

static void binary_read_64_bit(int read_len, uint8_t * read_seq, uint64_t* read_bit64_FOR){
	unsigned long int read_bit_char = (((uint16_t )((read_len >> 5) + 1)) << 3);
	memset(read_bit64_FOR, 0, read_bit_char);
	for (int r_i = 0;r_i < read_len; r_i++)
		read_bit64_FOR[r_i >> 5] |= (((uint64_t )read_seq[r_i]) << ((31 - (r_i & 0X1f)) << 1));
}

#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

inline char segment_type_to_cigar_code(const int id)
{
	switch (id)
	{
	case BAM_CMATCH:
		return 'M';
	case BAM_CINS:
		return 'I';
	case BAM_CDEL:
		return 'D';
	case BAM_CREF_SKIP:
		return 'N';
	case BAM_CSOFT_CLIP:
		return 'S';
	case BAM_CHARD_CLIP:
		return 'H';
	case BAM_CPAD:
		return 'P';
	case BAM_CEQUAL:
		return '=';
	case BAM_CDIFF:
		return 'X';
	case BAM_CBACK:
		return 'B';
	default:
		return 'X';
	}
}

void print_cigar_list(uint32_t* cigar, int cigar_l)
{
	for (int i = 0; i < cigar_l; ++i)
	{
		int type = (int)((cigar[i] & BAM_CIGAR_MASK));
		char c_type = segment_type_to_cigar_code(type);
		int length = (cigar[i] >> BAM_CIGAR_SHIFT);
		fprintf(stderr, "%d%c",length,  c_type);
	}
	fprintf(stderr, " ");
}

void print_cigar_list(std::vector<uint32_t>& c){
	print_cigar_list(&(c[0]), c.size());
}

int cigar_adjust(uint32_t *cigar_l_, uint32_t *cigar, bool add_blank, int delete_small_tail){
	FILE * output = stderr;
	if(*cigar_l_ == 0) return 0;
	int cigar_l = *cigar_l_;
	if(false){
		fprintf(output, "\noriginal cigar: \t");
		print_cigar_list(cigar, cigar_l);
		fprintf(output, "\n");
	}

	//merge adjacent SAME type CIGAR
	{
		int store_idx = 0;
		for(int cigar_ID = 1; cigar_ID < (int)cigar_l; cigar_ID++){
			int type 		= cigar[cigar_ID] & BAM_CIGAR_MASK;
			int type_store 	= cigar[store_idx] & BAM_CIGAR_MASK;
			if(type == type_store){
				cigar[store_idx] += ((cigar[cigar_ID] >> BAM_CIGAR_SHIFT) << BAM_CIGAR_SHIFT);
			}else{
				cigar[++store_idx] = cigar[cigar_ID];
			}
		}
		cigar_l = store_idx + 1;
	}

	int M_len = 0; int stable_cigar_ID = 0;
	//search to find a stable long match:
	for(int cigar_ID = 0;cigar_ID < (int)cigar_l; cigar_ID++){
	  int cigar_len = cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
	  int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
	  if(type == 0){//M
		if(cigar_len > delete_small_tail){ stable_cigar_ID = cigar_ID; break; }
		M_len += cigar_len;
	  }
	}
	int position_adjust = 0;
	int insertion_len = 0;
	if(stable_cigar_ID != 0){//need to be adjust at begin of alignment
		position_adjust += M_len;
		insertion_len += M_len;
		for(int cigar_ID = 0;cigar_ID < stable_cigar_ID; cigar_ID++){
			int cigar_len = cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
			if(type == 1) insertion_len += cigar_len;//INS
			else if(type == 2) position_adjust += cigar_len; //DEL
		}
		//store new cigar list:
		int new_cigar_len = 0;
		if(insertion_len != 0)
		cigar[new_cigar_len++] = (insertion_len << BAM_CIGAR_SHIFT) + 1;
		for(int cigar_ID = stable_cigar_ID; cigar_ID < (int)cigar_l; cigar_ID++)
			cigar[new_cigar_len++] = cigar[cigar_ID];
		cigar_l = new_cigar_len;
	}

	//adjust from tail of a cigar
	M_len = 0; stable_cigar_ID = 0;
	//search to find a stable long match:
	for(int cigar_ID = cigar_l - 1; cigar_ID >= 0; cigar_ID--){
		int cigar_len = cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
		int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
		if(type == 0){//M
			if(M_len + cigar_len > delete_small_tail){
				stable_cigar_ID = cigar_ID;
				break;
			}
			M_len += cigar_len;
		}
	}
	insertion_len = 0;
	if(stable_cigar_ID != cigar_l - 1){//need to be adjust at end of alignment
		insertion_len += M_len;
		for(int cigar_ID = cigar_l - 1;cigar_ID > stable_cigar_ID; cigar_ID--){
			int cigar_len = cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
			if(type == 1) insertion_len += cigar_len;//INS
		}
		//store new cigar list:
		if(insertion_len != 0){
			cigar[stable_cigar_ID + 1] = (insertion_len << BAM_CIGAR_SHIFT) + 1;
			cigar_l = stable_cigar_ID + 2;
		}
		else
			cigar_l = stable_cigar_ID + 1;
	}

	if(false){
		fprintf(output, "AfterAdj cigar: \t");
		print_cigar_list(cigar, cigar_l);
		fprintf(output, "\t new pos %d\n", position_adjust);
	}
 	xassert((int)(*cigar_l_) >= cigar_l, "");
	if(add_blank){
		for(int i = cigar_l; i < (int)(*cigar_l_); i++)
			cigar[i] = 0;
	}
	else
		*cigar_l_ = cigar_l;
	return position_adjust;
}

int getgapScore(KSW_ALN_handler &kswh, int length){
	int s1 = kswh.gap_open_D + (length*kswh.gap_ex_D);
	int s2 = kswh.gap_open2_D + (length*kswh.gap_ex2_D);
	return MIN(s1, s2);
}

int getAlignmentScore(KSW_ALN_handler &kswh, uint32_t* bam_cigar,int cigar_n, uint8_t* qseq, uint8_t *tseq){
	int number_of_match = 0;
	int number_of_mismatch = 0;
	int gapScore = 0;
	int tseq_i = 0;
	int qseq_i = 0;
	for(int cigar_ID = 0;cigar_ID < cigar_n; cigar_ID++){
	  int cigar_len = bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
	  int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
	  switch(type){
	  case 0:
		for(int i = 0; i < cigar_len; i++, qseq_i++){
		  if(qseq[qseq_i] != tseq[tseq_i])	number_of_mismatch++;//1 == M; 2 == X; 0 == -;
		  else 								number_of_match++;
		  tseq_i++;
		}
		break;//M
	  case 1: qseq_i += cigar_len; gapScore += getgapScore(kswh, cigar_len);break;//I, print nothing
	  case 2: tseq_i += cigar_len; gapScore += getgapScore(kswh, cigar_len);break;//D, print '-'
	  case 3: qseq_i += cigar_len; tseq_i += cigar_len; break;//N, print N
	  case 4: qseq_i += cigar_len; break;//S, print -
	  default: fprintf(stderr, "ERROR CIGAR  %d %d ", type, cigar_len);
	  }
	}
	int score = kswh.match_D * number_of_match - kswh.mismatch_D * number_of_mismatch - gapScore;
	if(score < 0) score = 0;
	return score;
}

void print_X_E_sequence(uint32_t* bam_cigar,int cigar_n, uint8_t* qseq, uint8_t *tseq){
	FILE * log_output = stderr;
	int output_index = 0;
	int seq_i = 0;
	//for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
	for(int cigar_ID = 0;cigar_ID < cigar_n; cigar_ID++){
	  int cigar_len = bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
	  int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
	  switch(type){
	  case 0:
		for(int i = 0; i < cigar_len; i++, seq_i++){
		  if(qseq[seq_i] == tseq[output_index])	fprintf(log_output, "=");//1 == M; 2 == X; 0 == -;
		  else									fprintf(log_output, "X");//1 == M; 2 == X; 0 == -;
		  output_index++;
		}
		break;//M
	  case 1: seq_i += cigar_len; break;//I, print nothing
	  case 2: for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;} break;//D, print '-'
	  case 3: for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "N"); output_index++;} break;//N, print N
	  case 4: seq_i += cigar_len; break;//S, print nothing
	  default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
	  }
	}
	fprintf(log_output, "\n");
}

#define MAX_BIN_READ_LEN_64 50 //1600/32 = 50
#define MIN_STR_REPEAT_COUNT 4 //when a kmer count >= MIN_STR_REPEAT_COUNT, that kmer will be treated as an STR kmer
#define MIN_STR_DETECT_LEN 15 //when a read has more than 15 STR-kmer, it will be treat as an STR read

//#define KMER_COUNT_ANALYSIS 1
int global_read_N = 0;
int get_ksw_score(bool printLog, uint8_t *read_bin[2][2], std::vector<Cluster*> &cluster_l,
		int first_node, int read_l, KSW_ALN_handler &kswh){
	printLog = false;
	int total_alignment_score = 0;
	for(Cluster* cluster_c : cluster_l){
		if(printLog)
			cluster_c->show(' ');
		//clear and set reads
		kswh.cigar.clear();
		kswh.setRead(read_bin[cluster_c->read_end_id][cluster_c->direction]);
		//debug code::
		int ref_chrID 	 = cluster_c->chrID;
		int MEM_ref_beg = cluster_c->ref_position_bg;
		int MEM_ref_end = cluster_c->ref_position_ed + 1;

		int MEM_read_beg = cluster_c->read_position_bg;
		int MEM_read_end = cluster_c->read_position_ed + 1;

		int EXT_left_len  = MEM_read_beg + 30;
		int EXT_right_len = read_l - MEM_read_end + 30;
		bool extern_left = (MEM_read_beg > 0);

		if(extern_left)	{
			kswh.alignment(printLog, ref_chrID, 0, 			  MEM_read_beg, MEM_ref_beg - EXT_left_len, MEM_ref_beg,	KSW_ALN_left_extend);
			if(printLog) print_cigar_list(kswh.ez.cigar, kswh.ez.n_cigar);
		}
		if(MEM_read_beg < MEM_read_end)	{
			kswh.alignment(printLog, ref_chrID, MEM_read_beg, MEM_read_end, MEM_ref_beg, 			    MEM_ref_end,	KSW_ALN_end_to_end);
			if(printLog) print_cigar_list(kswh.ez.cigar, kswh.ez.n_cigar);
		}
		if(MEM_read_end < read_l){
			kswh.alignment(printLog, ref_chrID, MEM_read_end, read_l, 		MEM_ref_end, MEM_ref_end + EXT_right_len, 	KSW_ALN_right_extend);
			if(printLog) print_cigar_list(kswh.ez.cigar, kswh.ez.n_cigar);
		}
		cluster_c->ksw_final_position = MEM_ref_beg;
		if(extern_left)
			cluster_c->ksw_final_position -= EXT_left_len;
		int adj_size = kswh.adjustCIGAR();
		cluster_c->ksw_final_position += adj_size;

		//store the cigar
		std::swap(cluster_c->kswCigar, kswh.cigar);
		//edgeIns to Clip
		if(!cluster_c->kswCigar.empty()){
			if((cluster_c->kswCigar[0] & BAM_CIGAR_MASK) == BAM_CINS){
				cluster_c->kswCigar[0] += 3;
			}
			if((cluster_c->kswCigar.back() & BAM_CIGAR_MASK) == BAM_CINS){
				cluster_c->kswCigar.back() += 3;
			}
		}

		//
		//generate the alignment score:
		std::vector<uint8_t > ref;
		kswh.idx->chrID_and_position_in_chr_to_global_position(cluster_c->ksw_final_position_global, ref_chrID, cluster_c->ksw_final_position);
		kswh.idx->get_refseq_char_bin(ref, cluster_c->ksw_final_position_global, 2*read_l);
		uint8_t * read_bin_p = read_bin[cluster_c->read_end_id][cluster_c->direction];
		cluster_c->alignment_score = getAlignmentScore(kswh, &(cluster_c->kswCigar[0]), cluster_c->kswCigar.size(), read_bin_p, &(ref[0]));
		total_alignment_score += cluster_c->alignment_score;
		if(printLog) {
			if(true){
				fprintf(stderr, "alignment_score is %d\t", cluster_c->alignment_score);
				print_cigar_list(cluster_c->kswCigar);
				print_X_E_sequence(&(cluster_c->kswCigar[0]), cluster_c->kswCigar.size(), read_bin_p, &(ref[0]));
			}

			if(true){
				cluster_c->show(' ');

				fprintf(stderr, "\n");
				fprintf(stderr, "read is ");
				uint8_t * read_bin_p = read_bin[cluster_c->read_end_id][cluster_c->direction];
				for(int i = 0 ; i < read_l; i++)
					fprintf(stderr, "%c", "ACGTN"[read_bin_p[i]]);
				fprintf(stderr, "\nref  is ");
				for(uint i = 0 ; i < ref.size(); i++)
					fprintf(stderr, "%c", "ACGTN"[ref[i]]);
				fprintf(stderr, "\n");
			}
		}
	}
	return total_alignment_score;
}

#define MAX_CHAIN_SOCRE_DIFF 30 //when chain score results are less the [best chain score - MIN_CHAIN_SOCRE_DIFF], the search will be break;
#define MIN_CHAIN_SOCRE 30 //when chain score results are less the [best chain score - MIN_CHAIN_SOCRE_DIFF], the search will be break;
#define MIN_ALN_SOCRE 40

#define  BAM_PAIRED         0x001
#define  BAM_PROPER_PAIR    0x002
#define  BAM_UNMAPPED       0x004
#define  BAM_MATE_UNMAPPED  0x008
#define  BAM_STRAND  		0x010
#define  BAM_MATE_STRAND    0x020
#define  BAM_FIRST_READ     0x040
#define  BAM_SECOND_READ    0x080
#define  BAM_SECONDARY      0x100
#define  BAM_FILTER         0x200
#define  BAM_DUPLICATE      0x400
#define  BAM_SUPPLEMENTARY  0x800

char getReverseChar(char c){
	switch(c){
	case 'A': case 'a': return 'T';
	case 'C': case 'c': return 'G';
	case 'G': case 'g': return 'C';
	case 'T': case 't': return 'A';
	}
	return 'N';
}

void getReverseStr_char(char * s, int len){
	int len_half = (len >> 1);
	for(int i = 0; i < len_half; i++){
		char tmp = s[i];
		s[i] = getReverseChar(s[len - 1 - i]);
		s[len - 1 - i] = getReverseChar(tmp);
	}
	if(len & 0x1){
		s[len_half] = getReverseChar(s[len_half]);
	}
}

void getReverseStr_qual_char(char * q, int len){
	int len_half = (len >> 1);
	for(int i = 0; i < len_half + 1; i++){
		int ri = len - 1 - i;
		char tmp = q[i];
		q[i ] = q[ri];
		q[ri] = tmp;
	}
}

void Single_Read_handler::binary_read_2_bit(){
	for (int r_i = 0; read_seq_char[r_i]; r_i++){
		char tmp_char = read_seq_char[r_i];
		if (tmp_char == 'N') tmp_char = "ACGT"[rand_r(&random_seed)%4]; //random
		uint8_t c_tmp = charToDna5n[(uint8_t)tmp_char];
		bin_read[0][r_i] = c_tmp;
		bin_read[1][(read_l - r_i - 1)] = (c_tmp ^ 0X3);
	}
}

bool bam_has_clip_or_unmapped_new(MAX_IDX_OUTPUT *c_max, int min_clip_len)
{
	if(c_max == NULL || c_max->cigar.empty()) return true;//unmapped
	int total_cigar_len = 0;
	for(uint32_t c: c_max->cigar){
		if((c & BAM_CIGAR_MASK) == 0){
			total_cigar_len += (c >> BAM_CIGAR_SHIFT);
		}
	}
	return (total_cigar_len >= min_clip_len);
}

void search_and_chain_seed_within_UNITIG(int direction_mode, int read_l, int read_end_id,
		std::vector<Chain_Within_Unitig> &CWU_l, deBGA_INDEX *idx, std::vector<MEM_Within_UNITIG> &MEM_within_unitig_rst_l,
		uint8_t *read_bin[2][2]){
	//step 1: binary reads
	uint64_t read_bit_64[MAX_BIN_READ_LEN_64];
	binary_read_64_bit(read_l, read_bin[read_end_id][direction_mode], read_bit_64);
	MEM_within_unitig_rst_l.clear();
	uint32_t kmer_number = read_l - LEN_KMER + 1;
	//test get repeat kmer, used for seed step 5.
	//seed list store whether a kmer can be selected as seed, 0:NO; >0: YES
	uint32_t max_search_right = 0;
	for(uint32_t read_off = 0; read_off < kmer_number; read_off += SEED_STEP) //every seed_l[tid] we choose a seed
	{
		if(read_off + LEN_KMER - 1 <= max_search_right) continue;
		uint64_t kmer_bit = getKmer(read_off, read_bit_64);
		int64_t range[2];
		bool search_rst = idx->search_kmer(LEN_KMER, kmer_bit, range, SEED_OFFSET);//search kmer in the index
		if(!search_rst) continue;

		uint64_t hit_bg = range[0];
		uint64_t hit_ed = range[1];

		if ((hit_ed - hit_bg + 1) > UNI_POS_N_MAX)
			continue;
		//for each kmer in the index:
		uint32_t max_right_i = 1;
		for (uint64_t hit_i = hit_bg; hit_i <= hit_ed; ++hit_i)
			idx->MEM_extern_within_UNITIG_for_KMER(hit_i, MEM_within_unitig_rst_l,
					read_bit_64, read_off, read_l, LEN_KMER, max_right_i);
		//skip seed when MEM covered
		max_search_right = read_off + LEN_KMER + max_right_i - 1;
	}
	int CWU_l_old_size = CWU_l.size();
	//merge seed in the unipath:
	idx->merge_seed_in_unipath(MEM_within_unitig_rst_l, CWU_l);

	//extend seed to reference, but the MEM is not done
	for(uint i = CWU_l_old_size; i < CWU_l.size(); i++){
		CWU_l[i].read_end_id = read_end_id;
		CWU_l[i].direction_mode = direction_mode;
	}
}


void Single_Read_handler::search_and_chain_seed_within_UNITIG(int direction_mode,
		std::vector<Chain_Within_Unitig> &CWU_l, deBGA_INDEX *idx){
	//step 1: binary reads
	uint64_t read_bit_64[MAX_BIN_READ_LEN_64];
	binary_read_64_bit(read_l, bin_read[direction_mode], read_bit_64);
	MEM_within_unitig_rst_l.clear();
	uint32_t kmer_number = read_l - LEN_KMER + 1;
	//test get repeat kmer, used for seed step 5.
	//seed list store whether a kmer can be selected as seed, 0:NO; >0: YES
	uint32_t max_search_right = 0;
	for(uint32_t read_off = 0; read_off < kmer_number; read_off += SEED_STEP) //every seed_l[tid] we choose a seed
	{
		if(read_off + LEN_KMER - 1 <= max_search_right) continue;
		uint64_t kmer_bit = getKmer(read_off, read_bit_64);
		int64_t range[2];
		bool search_rst = idx->search_kmer(LEN_KMER, kmer_bit, range, SEED_OFFSET);//search kmer in the index
		if(!search_rst) continue;

		uint64_t hit_bg = range[0];
		uint64_t hit_ed = range[1];

		if ((hit_ed - hit_bg + 1) > UNI_POS_N_MAX)
			continue;
		//for each kmer in the index:
		uint32_t max_right_i = 1;
		for (uint64_t hit_i = hit_bg; hit_i <= hit_ed; ++hit_i)
			idx->MEM_extern_within_UNITIG_for_KMER(hit_i, MEM_within_unitig_rst_l,
					read_bit_64, read_off, read_l, LEN_KMER, max_right_i);
		//skip seed when MEM covered
		max_search_right = read_off + LEN_KMER + max_right_i - 1;
	}
	int CWU_l_old_size = CWU_l.size();
	//merge seed in the unipath:
	idx->merge_seed_in_unipath(MEM_within_unitig_rst_l, CWU_l);

	//extend seed to reference, but the MEM is not done
	for(uint i = CWU_l_old_size; i < CWU_l.size(); i++){
		CWU_l[i].read_end_id = read_end_id;
		CWU_l[i].direction_mode = direction_mode;
	}
}

#define MAX_REF_DIS 300
#define MAX_READ_DIS 300
#define MAX_SEARCH_STEP 100 // = max read length / search step length  = 400 / 5 = 80
#define MAX_ABS_GAP 25

void sdp_chainging_seed_in_ref(std::vector<SEED_IN_REF>& SIR_l, std::vector<Cluster> &cluster_l, bool show_log){
	uint32_t seed_size = SIR_l.size();
	if(SIR_l.empty())
		return;
	std::sort(SIR_l.begin(), SIR_l.end(), SEED_IN_REF::cmp_by_ref_position);
	cluster_l.clear();

	//greedy clustering all nearby seeds to clusters
	for (uint32_t i = 0; i < seed_size; ++i){
		//skip already used node
		if(SIR_l[i].used_in_clustering)
			continue;
		//init the cluster
		Cluster & cur_cluster = cluster_l.emplace_back();
		cur_cluster.clear();
		cur_cluster.node_list.emplace_back();
		cur_cluster.node_list.back().init(&(SIR_l[i]));
		SIR_l[i].used_in_clustering = true;

		//search forward
		for (uint32_t j = i + 1; j < seed_size && j < i + MAX_SEARCH_STEP; ++j)
		{
			 //two mem generated from the same seed will not be connected.
			if(SIR_l[j].read_end_id != SIR_l[i].read_end_id)		break;
			if(SIR_l[j].chr_ID != SIR_l[i].chr_ID)					break;
			if(SIR_l[j].direction_mode != SIR_l[i].direction_mode)	break;
			if(SIR_l[j].seed_id == SIR_l[i].seed_id)				continue;

			int32_t dis_ref = (SIR_l[j].ref_begin - SIR_l[i].ref_end);
			int32_t dis_read = (SIR_l[j].read_begin - SIR_l[i].read_end);
			uint32_t abs_gap = ABS_U(dis_read, dis_ref);

			if(dis_ref > MAX_REF_DIS)  	break;
			if(dis_read > MAX_READ_DIS) continue;
			if(abs_gap > MAX_ABS_GAP) 	continue;

			cur_cluster.node_list.emplace_back();
			cur_cluster.node_list.back().init(&(SIR_l[j]));
			SIR_l[j].used_in_clustering = true;
		}
	}

	//dynamic programming chaining all nearby seeds in a cluster
	for (Cluster & clu: cluster_l){
		std::vector<Cluster_Node> &cur_node_list = clu.node_list;
		//DP for all nodes in the cluster： from the first node to last node
		for(uint nodeID = 0; nodeID < cur_node_list.size(); nodeID++){
			//search all the previous nodes to find the previous node
			for(uint nodeID_pre = 0; nodeID_pre < nodeID; nodeID_pre++){
				int32_t dis_ref  = (cur_node_list[nodeID].sir->ref_begin - cur_node_list[nodeID_pre].sir->ref_begin);
				int32_t dis_read = (cur_node_list[nodeID].sir->read_begin - cur_node_list[nodeID_pre].sir->read_begin);
				uint32_t abs_gap = ABS_U(dis_read, dis_ref);

				int overlap_read = cur_node_list[nodeID_pre].sir->read_end - cur_node_list[nodeID].sir->read_begin + 1;
				if(overlap_read < 0)
					overlap_read = 0;

				int penalty = (abs_gap == 0)?0: ((ABS(abs_gap)) + 3);// ABS(gap) * 1 + 3; open:3; ext:1
				uint32_t cur_sdp_score = cur_node_list[nodeID_pre].sdp_score +
						cur_node_list[nodeID].sir->cov - overlap_read - penalty;
				cur_node_list[nodeID].setNewScore(cur_sdp_score, &(cur_node_list[nodeID_pre]));
			}
		}
		//generate the max chain for the cluster
		clu.generateMAXChain();
	}
	return;
}

#define POS_N_MAX_T1 30
#define POS_N_MAX_T2 300

struct AlignmnetResult{

	int read1_pos_chaining = 0;
	int read2_pos_chaining = 0;

	int read1_pos_ksw = 0;
	int read2_pos_ksw = 0;

	int mapQ1;
	int mapQ2;

	int ISIZE_ksw;

	int ISIZE_chain;
	int chaining_score;
	std::vector<Cluster *> used_Cluster_l;
	Cluster * read1_cluster;
	Cluster * read2_cluster;

	//ksw alignment score
	int base_level_score;


	void storeCluster(Cluster * c){
		used_Cluster_l.emplace_back(c);
	}

	void setChaining_score_and_ISIZE(int read1_pos_chaining, int read2_pos_chaining, int chaining_score, int ISIZE_chain){
		this->read1_pos_chaining = read1_pos_chaining;
		this->read2_pos_chaining = read2_pos_chaining;
		this->chaining_score = chaining_score;
		this->ISIZE_chain = ISIZE_chain;
	}

	static inline int cmp_by_chaining_score(const AlignmnetResult &a, const AlignmnetResult &b){
		return a.chaining_score > b.chaining_score;
	}

	static inline int cmp_by_base_level_score(const AlignmnetResult &a, const AlignmnetResult &b){
		return a.base_level_score > b.base_level_score;
	}

	void showBasic(uint8_t *read_bin[2][2], int read_l, deBGA_INDEX * idx){
		fprintf(stderr, "Base_level_score %d, ISIZE %d \n", base_level_score, ISIZE_chain);
		for(Cluster* cluster_c : used_Cluster_l)
			cluster_c->showKSW_align_rst(read_bin, read_l, idx);
	}

};

void store_seed_in_ref_NEARBY_POS(deBGA_INDEX *idx, std::vector<Chain_Within_Unitig> &chain_with_in_UNI_l,
		int &seed_in_ref_id, std::vector<SEED_IN_REF> & seed_in_ref_l_nearby, bool show_log, int mate_read_end_id, int mate_ref_beg_global){
	for(Chain_Within_Unitig & cwu: chain_with_in_UNI_l){
		if(cwu.read_end_id == mate_read_end_id)
			continue;
		for(uint32_t m = 0; m < cwu.ref_pos_n; m++){
			uint64_t ref_beg_global = idx->buffer_p[m + idx->buffer_pp[cwu.uid]] + cwu.unitig_offset - 1;
			if(mate_ref_beg_global > ref_beg_global + 2000 || mate_ref_beg_global < ref_beg_global - 2000 )
				continue;
			int chrID; int ref_bg_in_chr;
			idx->global_position_to_chrID_and_position_in_chr(ref_beg_global, chrID, ref_bg_in_chr);
			int ref_end_in_chr = ref_bg_in_chr + cwu.length_in_ref - 1;
			uint32_t read_end = cwu.read_offset + cwu.length_in_read - 1;
			//global to ref-position in chromosomo
			seed_in_ref_l_nearby.emplace_back();
			seed_in_ref_l_nearby.back().set(
					seed_in_ref_id++, cwu.read_offset, read_end,
					chrID, ref_bg_in_chr, ref_end_in_chr, cwu.cov,
					cwu.read_end_id, cwu.direction_mode);
			if(show_log && false){
				seed_in_ref_l_nearby.back().print();
			}
		}
	}
}



void store_seed_in_ref(int complex_mode, deBGA_INDEX *idx, std::vector<Chain_Within_Unitig> &chain_with_in_UNI_l,
		int &seed_in_ref_id,
		std::vector<SEED_IN_REF> & seed_in_ref_l, bool show_log
	){
	for(Chain_Within_Unitig & cwu: chain_with_in_UNI_l){
		if(complex_mode == 0 && cwu.ref_pos_n > POS_N_MAX_T1)
			continue;
		if(complex_mode == 1 && (cwu.ref_pos_n <= POS_N_MAX_T1 || cwu.ref_pos_n > POS_N_MAX_T2))
			continue;

		for(uint32_t m = 0; m < cwu.ref_pos_n; m++){
			uint64_t ref_beg_global = idx->buffer_p[m + idx->buffer_pp[cwu.uid]] + cwu.unitig_offset - 1;
			int chrID; int ref_bg_in_chr;
			idx->global_position_to_chrID_and_position_in_chr(ref_beg_global, chrID, ref_bg_in_chr);
			int ref_end_in_chr = ref_bg_in_chr + cwu.length_in_ref - 1;
			uint32_t read_end = cwu.read_offset + cwu.length_in_read - 1;
			//global to ref-position in chromosomo
			seed_in_ref_l.emplace_back();
			seed_in_ref_l.back().set(
					seed_in_ref_id++, cwu.read_offset, read_end,
					chrID, ref_bg_in_chr, ref_end_in_chr, cwu.cov,
					cwu.read_end_id, cwu.direction_mode);
			if(show_log && false){
				seed_in_ref_l.back().print();
			}
		}
	}
	if(show_log && false){
		fprintf(stderr, "Chaining begin in complex_mode:[%d]\n", complex_mode);
		for(SEED_IN_REF & s: seed_in_ref_l)
			if(false) s.print();
	}
}

void cluster_pair_end(std::vector<Cluster> &cluster_l_curmode, std::vector<AlignmnetResult> &alignmnetResult){
	//pairing / chaining the two end of reads
	for(uint i = 0; i < cluster_l_curmode.size(); i++){
		//i == j
		alignmnetResult.emplace_back();
		alignmnetResult.back().storeCluster(&( cluster_l_curmode[i]));
		int read_pos_chaining = cluster_l_curmode[i].ref_position_bg;
		int ISIZE = 0;
		int final_score = cluster_l_curmode[i].max_score_chaining_score - 80;
		if(cluster_l_curmode[i].read_end_id == 0)
			alignmnetResult.back().setChaining_score_and_ISIZE(read_pos_chaining, -1, final_score, ISIZE);
		else
			alignmnetResult.back().setChaining_score_and_ISIZE(-1 ,read_pos_chaining, final_score, ISIZE);

		//i < j
		for(uint j = i + 1; j < cluster_l_curmode.size(); j++){
			if(cluster_l_curmode[i].read_end_id == cluster_l_curmode[j].read_end_id)
				continue;
			//detect whether reads is nearby
			int ISIZE_chaining =	0;
			if(cluster_l_curmode[i].chrID != cluster_l_curmode[j].chrID){
				ISIZE_chaining = MAX_int32t;
			}else if(cluster_l_curmode[i].direction == 0){
				ISIZE_chaining = cluster_l_curmode[j].ref_position_ed - cluster_l_curmode[i].ref_position_bg;
			}else
				ISIZE_chaining = cluster_l_curmode[i].ref_position_ed - cluster_l_curmode[j].ref_position_bg;

			//detect whether one of reads is new alignment
			int final_chaining_score = cluster_l_curmode[i].max_score_chaining_score + cluster_l_curmode[j].max_score_chaining_score +
					+ ((ISIZE_chaining < 1000 && ISIZE_chaining > -1000)?(0):(-80)); //not nearby penalty
			if(final_chaining_score > 60){//at least 60
				int read1_pos_chaining = cluster_l_curmode[i].ref_position_bg;
				int read2_pos_chaining = cluster_l_curmode[j].ref_position_bg;
				if(cluster_l_curmode[i].read_end_id == 1)
					std::swap(read1_pos_chaining, read2_pos_chaining);
				alignmnetResult.emplace_back();
				alignmnetResult.back().storeCluster(&( cluster_l_curmode[i]));
				alignmnetResult.back().storeCluster(&( cluster_l_curmode[j]));
				alignmnetResult.back().setChaining_score_and_ISIZE(read1_pos_chaining, read2_pos_chaining, final_chaining_score, ISIZE_chaining);
			}
		}
	}
	std::sort(alignmnetResult.begin(), alignmnetResult.end(), AlignmnetResult::cmp_by_chaining_score);
}

int base_level_alignment(std::vector<AlignmnetResult> &alignmnetResult, bool show_log, uint8_t *read_bin[2][2], int read_len, KSW_ALN_handler &kswh){
	//base level alignment
	//ksw-base-level alignment
	int cur_ksw_idx = 0;
	int true_ksw_count = 0;
	for(AlignmnetResult & ar: alignmnetResult){
		if(cur_ksw_idx ++ >= 8)
			break;
		if(ar.chaining_score < alignmnetResult[0].chaining_score - 100)
			break;
		true_ksw_count++;
		int base_level_score = get_ksw_score(show_log & false, read_bin, ar.used_Cluster_l, 0, read_len, kswh);
		base_level_score += ((ar.ISIZE_chain < -1000 || ar.ISIZE_chain > 1000)?(-150):0);
		ar.base_level_score = base_level_score;
		//set the ksw score/ISIZE/positions
		ar.read1_cluster = NULL;
		ar.read2_cluster = NULL;

		for(Cluster* cluster_c : ar.used_Cluster_l){
			if(cluster_c->read_end_id == 0) ar.read1_cluster = cluster_c;
			if(cluster_c->read_end_id == 1) ar.read2_cluster = cluster_c;
		}

		ar.read1_pos_ksw = ((ar.read1_cluster == NULL)?(-1):(ar.read1_cluster->ksw_final_position));
		ar.read2_pos_ksw = ((ar.read2_cluster == NULL)?(-1):(ar.read2_cluster->ksw_final_position));
	}
	std::sort(alignmnetResult.begin(), alignmnetResult.begin() + true_ksw_count, AlignmnetResult::cmp_by_base_level_score);
	return true_ksw_count;
}

struct BENCH_data{
	bool readOverlapSV;
	bool run_aln;
	int SV_pos[2] = {0};
	int SV_end[2] = {0};
	int SVLen[2] = {0};
	int read_IN_SV[2] = {0};
	int mapping_IN_SV[2] = {0};
	int diff_read_POS[2] = {0};
	int diff_read_pos2[2] = {0};
	int read1_pos_true = 0;
	int read2_pos_true = 0;
	int middle_in_INS[2] = {0};

	void scan_comment(char * read1_comment, char * read2_comment){
		readOverlapSV = false;
		run_aln = false;
		std::vector<std::string> item_value;
		split_string(item_value, read1_comment , "_");
		read1_pos_true = atoi(item_value[1].c_str());

		if(item_value.size() > 3){
			readOverlapSV = true;
			SV_pos[0] = atoi(item_value[3].c_str());
			SV_end[0] = atoi(item_value[4].c_str());
			SVLen[0] =  atoi(item_value[6].c_str()) - atoi(item_value[5].c_str()) - (SV_end[0] - SV_pos[0]);
			read_IN_SV[0] = read1_pos_true - SV_pos[0];
		}

		split_string(item_value, read2_comment , "_");
		read2_pos_true = atoi(item_value[1].c_str());

		if(item_value.size() > 3){
			readOverlapSV = true;
			SV_pos[1] = atoi(item_value[3].c_str());
			SV_end[1] = atoi(item_value[4].c_str());
			SVLen[1] =  atoi(item_value[6].c_str()) - atoi(item_value[5].c_str()) - (SV_end[1] - SV_pos[1]);
			read_IN_SV[1] = read2_pos_true - SV_pos[1];
		}

		int DEBUG_POS = 41050311;
		bool run_aln = false;
		if((read1_pos_true > DEBUG_POS - 200  && read1_pos_true < DEBUG_POS + 200) ||
				(read2_pos_true > DEBUG_POS - 200  && read2_pos_true < DEBUG_POS + 200) ||
				(SV_pos[0] == DEBUG_POS || SV_pos[1] == DEBUG_POS)){
			run_aln = true;
		}
		if(SVLen[0] > 0 && read_IN_SV[0] > -60 && read_IN_SV[0] < SVLen[0] - 90)
			middle_in_INS[0] = true;
		if(SVLen[1] > 0 && read_IN_SV[1] > -60 && read_IN_SV[1] < SVLen[1] - 90)
			middle_in_INS[1] = true;

		//if(run_aln == false) return;
		if(readOverlapSV == false)  return;
		if(middle_in_INS[0] && middle_in_INS[1]) return;
	}

	void bench(int read1_pos_final, int read2_pos_final, int *mapq, kseq_t * read1, kseq_t * read2){
		//benchmarks
		diff_read_POS[0] = read1_pos_final - read1_pos_true;
		diff_read_POS[1] = read2_pos_final - read2_pos_true;

		diff_read_pos2[0] = read1_pos_final - (read1_pos_true - SVLen[0]);
		diff_read_pos2[1] = read2_pos_final - (read2_pos_true - SVLen[1]);

		int wrong_aln[2] = {0};
		//check:
		int ERROR_region = 2000;
		if((diff_read_POS[0] > ERROR_region || diff_read_POS[0] < - ERROR_region) &&
			(diff_read_pos2[0] > ERROR_region || diff_read_pos2[0] < - ERROR_region)){
			wrong_aln[0] = 1;
		}
		if((diff_read_POS[1] > ERROR_region || diff_read_POS[1] < - ERROR_region) &&
			(diff_read_pos2[1] > ERROR_region || diff_read_pos2[1] < - ERROR_region)){
			wrong_aln[1] = 1;
		}

		mapping_IN_SV[0] = read2_pos_final - SV_pos[0];
		mapping_IN_SV[1] = read2_pos_final - SV_pos[1];

		if(readOverlapSV){
			if(true){
				fprintf(stderr, "SV_search\n");
				fprintf(stderr, "SV_read1 %s:wrong_aln:%s;mid_ins:%s;POS=%d;END=%d;LENGTH=%d;SVTYPE=%s;READ_IN_SV=%d;mapping_IN_SV=%d;\n",
						read1->name.s, (wrong_aln[0]==0)?"false":"true",(middle_in_INS[0]==0)?"false":"true", SV_pos[0],SV_end[0],SVLen[0],
								(SVLen[0] > 0)?"INS":"DEL",read_IN_SV[0], mapping_IN_SV[0]);
				fprintf(stderr, "SV_read2 %s:wrong_aln:%s;mid_ins:%s;POS=%d;END=%d;LENGTH=%d;SVTYPE=%s;READ_IN_SV=%d;mapping_IN_SV=%d;\n",
						read2->name.s, (wrong_aln[1]==0)?"false":"true",(middle_in_INS[1]==0)?"false":"true", SV_pos[1],SV_end[1],SVLen[1],
								(SVLen[1] > 0)?"INS":"DEL",read_IN_SV[1], mapping_IN_SV[1]);
				fprintf(stderr, "READ1( %d :mapq %d ) bench CMP: %d, %d, %d, %d, %d \n", global_read_N,
						mapq[0], read1_pos_final, (read1_pos_true - SVLen[0]), read1_pos_true, diff_read_POS[0],
						diff_read_pos2[0]);
				fprintf(stderr, "READ2( %d :mapq %d ) bench CMP: %d, %d, %d, %d, %d \n", global_read_N,
						mapq[1], read2_pos_final, (read2_pos_true - SVLen[1]), read2_pos_true, diff_read_POS[1], diff_read_pos2[1]);
				//debug show reads
				if(false){
					fprintf(stderr, "[read1]\n"); kseq_show(stderr, read1);
					fprintf(stderr, "[read2]\n"); kseq_show(stderr, read2);
				}
			}
		}
	}
};

void setMAPQ_and_finalPos(std::vector<AlignmnetResult> &alignmnetResult, Cluster * read_alignment_CIGAR[2],
		int *pos_final, int *mapq){
	pos_final[0] = pos_final[1] = -1;
	mapq[0] = mapq[1] = 0;
	read_alignment_CIGAR[0] = read_alignment_CIGAR[1] = NULL;

	if(alignmnetResult.size() >= 1){
		pos_final[0] = alignmnetResult[0].read1_pos_ksw;
		pos_final[1] = alignmnetResult[0].read2_pos_ksw;
		read_alignment_CIGAR[0] =  alignmnetResult[0].read1_cluster;
		read_alignment_CIGAR[1] =  alignmnetResult[0].read2_cluster;
		if(alignmnetResult.size() == 1){
			mapq[0] = mapq[1] = 60;
		}else{
			int mapq_ALL = (alignmnetResult[0].base_level_score - alignmnetResult[1].base_level_score);
			mapq_ALL = MIN(60,mapq_ALL);
			if(pos_final[0] == alignmnetResult[1].read1_pos_ksw)
				mapq[0] = 60;
			else
				mapq[0] = mapq_ALL;

			if(pos_final[1] == alignmnetResult[1].read2_pos_ksw)
				mapq[1] = 60;
			else
				mapq[1] = mapq_ALL;
		}
	}
	if(pos_final[0] == -1 ) {pos_final[0] = pos_final[1]; mapq[0] = mapq[1];}
	if(pos_final[1] == -1 ) {pos_final[1] = pos_final[0]; mapq[1] = mapq[0];}
}

uint8_t setFlag(int read_end_id, int direction, int isPairEndRead, int isProper_paired,
		int isUnmapped, int isMateUnmapped){
	uint8_t flag =
			((read_end_id == 0)?BAM_FIRST_READ:BAM_SECOND_READ) +
			((direction == REVERSE)?BAM_MATE_STRAND:BAM_STRAND) +
			((isPairEndRead)?BAM_PAIRED:0) +
			((isProper_paired)?BAM_PROPER_PAIR:0) +
			((isUnmapped)?BAM_UNMAPPED:0) +
			((isMateUnmapped)?BAM_MATE_UNMAPPED:0);
	return flag;
}
void generateSamStringUnmapped(int read_end_id, std::string & rst, kseq_t * c_read, kstring_t &s_buff){
	uint8_t flag = setFlag(read_end_id, FORWARD, 1, 0, 1, 1);

	s_buff.l = 0;
	sprintf(s_buff.s + s_buff.l,	"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s",
			c_read->name.s, flag, c_read->seq.s, c_read->qual.s); s_buff.l += strlen(s_buff.s + s_buff.l);//read name//flag//chr_name //mapq, set to blank at first//ref begin
	//final
	s_buff.s[s_buff.l] = 0;
	//output result:
	rst = s_buff.s;
}

void generateSamStringSingleUnmapped(std::string & rst, kseq_t * c_read, int mapq, Cluster *mate_align, int ISIZE, deBGA_INDEX *idx, kstring_t &s){
	//skip original alignment on decoy reference:
	//init:
	int direction = 1-mate_align->direction;
	int read_end_id = mate_align->read_end_id;
	uint8_t flag = setFlag(read_end_id, direction, 1, 0, 1, 0);

	s.l = 0;
	//BASIC PART:
	int read_l = c_read->seq.l;
	sprintf(s.s + s.l,	"%s\t%d\t%s\t%d\t%d\t", c_read->name.s, flag, idx->chr_names[mate_align->chrID].c_str(),
			mate_align->ksw_final_position + 1 , mapq); s.l += strlen(s.s + s.l);//read name//flag//chr_name //mapq, set to blank at first//ref begin
	//CIGAR
	sprintf(s.s + s.l, "*\t"	); s.l += strlen(s.s + s.l);
	// mate pos and ISIZE: isize is simple 500 for mated read
	if(mate_align != NULL)
		sprintf(s.s + s.l, "%s\t%d\t%d\t", idx->chr_names[mate_align->chrID].c_str(), mate_align->ksw_final_position + 1, 0);
	else
		sprintf(s.s + s.l,	"*\t0\t0\t");
	s.l += strlen(s.s + s.l);
	//seq//qual //align score
	sprintf(s.s + s.l,  "%s\t%s\t" , c_read->seq.s, c_read->qual.s); s.l += strlen(s.s + s.l);
	sprintf(s.s + s.l, "RC:Z:%s", c_read->comment.s); s.l += strlen(s.s + s.l);
	//final
	s.s[s.l] = 0;
	//output result:
	rst = s.s;
}

//set tid to be -1 to discard the output of BAM
void generateSamString(std::string & rst, Cluster *read_align,  kseq_t * c_read, int mapq,
		Cluster *mate_align, int ISIZE, deBGA_INDEX *idx, kstring_t &s)
{	//get primary alignment and mapQ
	//skip original alignment on decoy reference:
	//init:
	int direction = read_align->direction;
	int read_end_id = read_align->read_end_id;
	//todo::TOO simple
	int	isProper_paired = (ISIZE < -4000 || ISIZE > 4000 )?0:1;
	uint8_t flag = setFlag(read_end_id, direction, 1, isProper_paired, 0, 0);
	s.l = 0;
	//BASIC PART:
	int read_l = c_read->seq.l;
	sprintf(s.s + s.l,	"%s\t%d\t%s\t%d\t%d\t", c_read->name.s, flag, idx->chr_names[read_align->chrID].c_str(),
			read_align->ksw_final_position + 1, mapq); s.l += strlen(s.s + s.l);//read name//flag//chr_name //mapq, set to blank at first//ref begin
	//CIGAR
	for(uint32_t c : read_align->kswCigar)
	{ sprintf(s.s + s.l, "%d%c", (c >> BAM_CIGAR_SHIFT), BAM_CIGAR_STR[c & BAM_CIGAR_MASK]); s.l += strlen(s.s + s.l); }
	sprintf(s.s + s.l, "\t"	); s.l ++;
	// mate pos and ISIZE: isize is simple 500 for mated read
	int ABS_isize = ABS(ISIZE);
	int isize = (direction == FORWARD)?(ABS_isize):(-ABS_isize);
	if(mate_align != NULL)
		sprintf(s.s + s.l, "%s\t%d\t%d\t", idx->chr_names[mate_align->chrID].c_str(), mate_align->ksw_final_position + 1, isize);
	else
		sprintf(s.s + s.l, "%s\t%d\t%d\t", idx->chr_names[read_align->chrID].c_str(), read_align->ksw_final_position + 1, isize);
	s.l += strlen(s.s + s.l);
	//seq//qual //align score
	if(direction == REVERSE){
		getReverseStr_char(c_read->seq.s, read_l);
		getReverseStr_qual_char(c_read->qual.s, read_l);
	}
	sprintf(s.s + s.l,  "%s\t%s\t" "AS:i:%d\t", c_read->seq.s, c_read->qual.s, read_align->alignment_score); s.l += strlen(s.s + s.l);
	if(direction == REVERSE){ //restore to original after used
		getReverseStr_char(c_read->seq.s, read_l);
		getReverseStr_qual_char(c_read->qual.s, read_l);
	}
	//chaining score ,when it is new alignment results
	sprintf(s.s + s.l, "CS:i:%d\t", read_align->max_score_chaining_score); s.l += strlen(s.s + s.l);

//	//secondary alignment
//	if(secondary_result != NULL){//???? only one??
//		sprintf(s.s + s.l, "XA:Z:%d,%d,%d,", secondary_result->chrID, secondary_result->ref_bg, secondary_result->read_bg); s.l += strlen(s.s + s.l);
//		sprintf(s.s + s.l, "%d,%c\t", secondary_result->align_score, (secondary_result->direction == FORWARD)?'F':'R');	s.l += strlen(s.s + s.l);
//	}
	//todo::
	sprintf(s.s + s.l, "RC:Z:%s", c_read->comment.s); s.l += strlen(s.s + s.l);

	//final
	s.s[s.l] = 0;
	//output result:
	rst = s.s;
}
void binary_read_2_bit(uint8_t *read_bin[2], int read_l, char * read_seq_char,uint *random_seed){
	for (int r_i = 0; read_seq_char[r_i]; r_i++){
		char tmp_char = read_seq_char[r_i];
		if (tmp_char == 'N') tmp_char = "ACGT"[rand_r(random_seed)%4]; //random
		uint8_t c_tmp = charToDna5n[(uint8_t)tmp_char];
		read_bin[0][r_i] = c_tmp;
		read_bin[1][(read_l - r_i - 1)] = (c_tmp ^ 0X3);
	}
}

void align_read_pair(kseq_t * read1, kseq_t * read2, Classify_buff_pool * buff, std::string & rst_1,std::string & rst_2,
		FILE * debug_1, FILE * debug_2){
	bool show_log = true;
	global_read_N++;
	//fprintf(stderr, "\n\nglobal_read_N is %d\n", global_read_N);
	int read_l = read1->seq.l;
	//get comments
	BENCH_data benchData;
	benchData.scan_comment(read1->comment.s, read2->comment.s);

	//restore reads to FASTQ
	if(false){
		fprintf(debug_1, "@%s_%s!%s %s\n%s\n+\n%s\n", read1->name.s	, read1->comment.s, read2->comment.s, read1->comment.s, read1->seq.s, read1->qual.s);
		fprintf(debug_2, "@%s_%s!%s %s\n%s\n+\n%s\n", read2->name.s, read1->comment.s, read2->comment.s, read2->comment.s, read2->seq.s, read2->qual.s);
		return;
	}
	//if(strcmp(read1->comment.s, "@chr22_41052938_1_41050311_41050311_41097776_41100475") != 0) return;

	//step 1: generate all unitig-seed for both single end reads
	//step 1: align both single end reads separately
	deBGA_INDEX *idx = buff->idx;
	std::vector<Chain_Within_Unitig> &chain_with_in_UNI_l = buff->CWU_l;
	chain_with_in_UNI_l.clear();

	if(show_log && false)
		fprintf(stderr, "Begin searching seed in unitig\n\n\n");
	//store the binary reads
	uint8_t read_bin_BUFF[2][2][MAX_READ_LEN];
	uint8_t *read_bin[2][2];
	for(int read_end_id = 0; read_end_id < 2; read_end_id++){
		read_bin[read_end_id][0] = read_bin_BUFF[read_end_id][0];
		read_bin[read_end_id][1] = read_bin_BUFF[read_end_id][1];
		binary_read_2_bit(read_bin[read_end_id], read_l, (read_end_id == 0)?read1->seq.s:read2->seq.s, &(buff->random_seed));
	}

	//todo:: BUFF
	std::vector<MEM_Within_UNITIG> MEM_within_unitig_rst_l;
	//search seeds in index
	for(int read_end_id = 0; read_end_id < 2; read_end_id++)
		for(int direction_mode = 0; direction_mode < 2; direction_mode ++)
			search_and_chain_seed_within_UNITIG(direction_mode, read_l, read_end_id, chain_with_in_UNI_l, idx, MEM_within_unitig_rst_l, read_bin);

	if(show_log &&false){
		for(Chain_Within_Unitig & cwu: chain_with_in_UNI_l)
			 cwu.print();
		fprintf(stderr, "End searching seed in unitig\n\n\n");
	}
	//chaining the read pair using the simple "Chain_Within_Unitig"
	//additional chaining the read pair using the complexz "Chain_Within_Unitig"
	std::vector<SEED_IN_REF> & seed_in_ref_l = buff->SIR_l;


	int seed_in_ref_id = 0;
	//simple result
	int pos_final[2];
	int mapq[2];
	Cluster * read_align_rst[2];

	std::vector<AlignmnetResult> alignmnetResult;
	alignmnetResult.clear();
	int true_ksw_count = 0;
	for(int complex_mode = 0; complex_mode < 1; complex_mode ++){
		if(show_log && false)fprintf(stderr, "In complex_mode [%d]\n", complex_mode);

		std::vector<Cluster> &cluster_l_curmode = buff->cluster_l[complex_mode];//graph node
		seed_in_ref_l.clear();
		store_seed_in_ref(complex_mode, idx, chain_with_in_UNI_l, seed_in_ref_id, seed_in_ref_l, show_log);
		sdp_chainging_seed_in_ref(seed_in_ref_l, cluster_l_curmode, show_log);//chaining all result in the reference
		cluster_pair_end(cluster_l_curmode, alignmnetResult);//pairing / chaining the two end of reads
		true_ksw_count = base_level_alignment(alignmnetResult, show_log, read_bin, read1->seq.l, buff->kswh);//ksw-base-level alignment
		//store the final results
		setMAPQ_and_finalPos(alignmnetResult, read_align_rst, pos_final, mapq);
	}
	if(mapq[0] < 0) mapq[0] = 0;
	if(mapq[1] < 0) mapq[1] = 0;

	//nearby read recover
	{
		int mate_read_id = -1;
		uint64_t mate_position_global = -1;
		if(read_align_rst[0] == NULL && read_align_rst[1] != NULL){
			mate_read_id = 1;
			mate_position_global = read_align_rst[1]->ksw_final_position_global;
		}
		if(read_align_rst[1] == NULL && read_align_rst[0] != NULL){
			mate_read_id = 0;
			mate_position_global = read_align_rst[0]->ksw_final_position_global;
		}
		if(mate_read_id != -1){
			fprintf(stderr, "Nearby read recover begin!\n\n");

			std::vector<SEED_IN_REF> seed_in_ref_l_nearby;
			store_seed_in_ref_NEARBY_POS(idx, chain_with_in_UNI_l, seed_in_ref_id, seed_in_ref_l_nearby,
					show_log, mate_read_id, mate_position_global);
			if(!seed_in_ref_l_nearby.empty()){
				std::vector<Cluster> cluster_l_curmode;
				std::vector<AlignmnetResult> alignmnetResult;
				sdp_chainging_seed_in_ref(seed_in_ref_l_nearby, cluster_l_curmode, show_log);//chaining all result in the reference
				cluster_pair_end(cluster_l_curmode, alignmnetResult);//pairing / chaining the two end of reads
				true_ksw_count = base_level_alignment(alignmnetResult, show_log, read_bin, read1->seq.l, buff->kswh);//ksw-base-level alignment
				int cur_ksw_idx = 0;
				for(AlignmnetResult & ar: alignmnetResult){
					if(cur_ksw_idx++ >= true_ksw_count)
						break;
					if(show_log && false) fprintf(stderr, "Begin @[%d]\n", cur_ksw_idx);
					ar.showBasic(read_bin, read1->seq.l, idx);
					if(false){
						fprintf(stderr, "Test Begin @[%d]\n", cur_ksw_idx);
						get_ksw_score(show_log, read_bin, ar.used_Cluster_l, 0, read1->seq.l, buff->kswh);
					}
				}
				fprintf(stderr, "\n\n\n\n");
			}
		}
	}

	//generate SAM string:
	//
	int ISIZE = 0;
	kstring_t s_buff;
	char sam_string[2048];
	s_buff.s = sam_string; s_buff.m = 2048;
	if(read_align_rst[0] == NULL && read_align_rst[1] == NULL){
		generateSamStringUnmapped(0, rst_1, read1, s_buff);
		generateSamStringUnmapped(1, rst_2, read2, s_buff);
	}else if(read_align_rst[0] == NULL){
		generateSamStringSingleUnmapped(rst_1, read1, mapq[0], read_align_rst[1], ISIZE, idx, s_buff);
		generateSamString(rst_2, read_align_rst[1], read2, mapq[1], read_align_rst[0], ISIZE, idx, s_buff);
	}else if(read_align_rst[1] == NULL){
		generateSamString(rst_1, read_align_rst[0], read1, mapq[0], read_align_rst[1], ISIZE, idx, s_buff);
		generateSamStringSingleUnmapped(rst_2, read2, mapq[1], read_align_rst[0], ISIZE, idx, s_buff);
	}
	else{
		generateSamString(rst_1, read_align_rst[0], read1, mapq[0], read_align_rst[1], ISIZE, idx, s_buff);
		generateSamString(rst_2, read_align_rst[1], read2, mapq[1], read_align_rst[0], ISIZE, idx, s_buff);
	}

	benchData.bench(pos_final[0], pos_final[1], mapq, read1, read2);

	if(false){
		std::vector<Cluster> cc = buff->cluster_l[0];
		std::sort(cc.begin(), cc.end(), Cluster::cmp_by_max_score);
		int max_show = 30;
		for (Cluster & clu: cc){
			if(max_show-- <= 0) break;
			clu.show('\n');
		}
	}

	if(false){
		fprintf(stderr, "T3");
		int cur_ksw_idx = 0;
		for(AlignmnetResult & ar: alignmnetResult){
			if(cur_ksw_idx++ >= true_ksw_count)
				break;
			if(show_log && false) fprintf(stderr, "Begin @[%d]\n", cur_ksw_idx);
			ar.showBasic(read_bin, read1->seq.l, idx);
			if(false){
				fprintf(stderr, "Test Begin @[%d]\n", cur_ksw_idx);
				get_ksw_score(show_log, read_bin, ar.used_Cluster_l, 0, read1->seq.l, buff->kswh);
			}
		}
		fprintf(stderr, "\n\n\n\n");
	}

}

void getReverseStr_uint8_t(uint8_t * q, int len){
	int len_half = (len >> 1);
	for(int i = 0; i < len_half; i++){
		int ri = len - 1 - i;
		uint8_t tmp = q[i ];
		q[i ] = q[ri];
		q[ri] = tmp;
	}
}

//**************************************class: KSW_ALN_handler********************************/

void KSW_ALN_handler::copy_option(MAP_PARA 	*o){
	match_D= o->match_D ;
	mismatch_D= o->mismatch_D ;
	gap_open_D= o->gap_open_D ;
	gap_ex_D= o->gap_ex_D ;
	gap_open2_D= o->gap_open2_D ;
	gap_ex2_D= o->gap_ex2_D ;
	zdrop_D= o->zdrop_D ; //for DNA zdrop = 400, 200 for RNA
	bandwith = 200;
	flag = 0;
}

void KSW_ALN_handler::ksw_gen_mat_D(){
	int8_t l,k,m;
	for (l = k = 0; l < 4; ++l) {
		for (m = 0; m < 4; ++m) {
			mata_D[k] = l == m ? match_D : -(mismatch_D);	/* weight_match : -weight_mismatch */
			k++;
		}
		mata_D[k] = 0; // ambiguous base
		k++;
	}
	for (m = 0; m < 5; ++m) {
		mata_D[k] = 0;
		k++;
	}
}

void KSW_ALN_handler::init(MAP_PARA 	*o, deBGA_INDEX * idx_){
	tseq = (uint8_t *)xmalloc(1600);
	qseq_rev = (uint8_t *)xmalloc(1600);
	km = km_init();
	memset(&ez, 0, sizeof(ksw_extz_t));
	//mapping options
	copy_option(o);
	ksw_gen_mat_D();
	//index
	idx = idx_;
}

void KSW_ALN_handler::free_memory(){
	if(ez.cigar)   kfree(km, (void *)(ez.cigar));
	km_destroy(km);
	free(tseq);
	free(qseq_rev);
}

void KSW_ALN_handler::setRead(uint8_t *read_str_){
	read_str = read_str_;
	//n_cigar = 0;
	total_q_len = 0;
}

void KSW_ALN_handler::align_non_splice()
{
	if ((int64_t)tlen * qlen > 1000000)
	{
		ksw_reset_extz(&ez);
		if(ez.m_cigar < 2)
		{
			ez.n_cigar = 2;
			ez.m_cigar = (ez.n_cigar)<<2;
			ez.cigar = (uint32_t *)krealloc(km, (void *)ez.cigar, (ez.m_cigar)<<4);
		}
		ez.n_cigar = 2;
		ez.cigar[0] = qlen<<4 | 1;
		ez.cigar[1] = tlen<<4 | 3;
		ez.score = 0;
	}
	else{
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}
}

int KSW_ALN_handler::get_misMatch(int read_st_, int read_ed_, int ref_st_, int ref_ed_){
	qlen = read_ed_ - read_st_;
	qseq = read_str + read_st_;
	tlen = ref_ed_ - ref_st_;
	if(ref_ed_ < ref_st_){//when ref_len < 0
		tlen = 0;
		qlen += (ref_st_ - ref_ed_);
	}

	xassert(tlen < 1600, "");
	idx->get_refseq(tseq, tlen, ref_st_);
	int c_simple_NM = 0;
	for(uint32_t i = 0; i < qlen; c_simple_NM += (qseq[i] != tseq[i])?1:0, i++);
	if(c_simple_NM > 3)	c_simple_NM = 3;
	return c_simple_NM;
}

void KSW_ALN_handler::alignment(bool printLog, int ref_chrID, int read_st_, int read_ed_, int ref_st_, int ref_ed_, int type_){
	qlen = read_ed_ - read_st_;
	qseq = read_str + read_st_;
	tlen = ref_ed_ - ref_st_;
	if(ref_ed_ < ref_st_){//when ref_len < 0
		tlen = 0;
		qlen += (ref_st_ - ref_ed_);
	}
	ez.n_cigar = 0;
	xassert(tlen < 1600, "");
	uint64_t position_global;
	idx->chrID_and_position_in_chr_to_global_position(position_global, ref_chrID, ref_st_);
	idx->get_refseq(tseq, tlen, position_global);
	type = type_;

	if(printLog){
		fprintf(stderr, "KSW_ALN_handler::alignment: show read\n");
		for(uint i = 0; i < qlen; i++){
			fprintf(stderr, "%c", "ACGTN"[qseq[i]]);
		}
		fprintf(stderr, "\nKSW_ALN_handler::alignment: show ref\n");
		for(uint i = 0; i < tlen; i++){
			fprintf(stderr, "%c", "ACGTN"[tseq[i]]);
		}
		fprintf(stderr, "\nKSW_ALN_handler::alignment: end show\n");
	}

	if(type == KSW_ALN_left_extend){
		getReverseStr_uint8_t(tseq, tlen);
		memcpy(qseq_rev, qseq, qlen);
		getReverseStr_uint8_t(qseq_rev, qlen);
		qseq = qseq_rev;
	}

	total_q_len += qlen;

	//try simple compare
	is_simple_aln = false;
	simple_NM = 0;
	if(qlen == 0 || tlen == 0){
		is_simple_aln = true;
		simple_NM = qlen + tlen;
	}else if(qlen == tlen || (type != KSW_ALN_end_to_end && qlen <= tlen)){//when [global search and length is equal] or [extension search]
		for(uint32_t i = 0; i < qlen && simple_NM < 6; simple_NM += (qseq[i] != tseq[i])?1:0, i++);
		if(simple_NM == 1 || (simple_NM < 6 && ((simple_NM << 3) < qlen)))
			is_simple_aln = true;
	}

	if(!is_simple_aln)
		 align_non_splice();

	//store score and cigar:
	//set read score
	if(is_simple_aln){
		if(type == KSW_ALN_left_extend){
			if(tlen > qlen)
				cigar.emplace_back(((tlen - qlen) << BAM_CIGAR_SHIFT) + BAM_CDEL);
			cigar.emplace_back(((qlen) << BAM_CIGAR_SHIFT) + BAM_CMATCH);
		}
		else if(type == KSW_ALN_end_to_end){
			cigar.emplace_back(((qlen) << BAM_CIGAR_SHIFT) + BAM_CMATCH);
		}else if(type == KSW_ALN_right_extend){
			cigar.emplace_back(((qlen) << BAM_CIGAR_SHIFT) + BAM_CMATCH);
			if(tlen > qlen)
				cigar.emplace_back(((tlen - qlen) << BAM_CIGAR_SHIFT) + BAM_CDEL);
		}
	}
	else{
		if(type == KSW_ALN_end_to_end){ //global
			for(int i = 0; i < ez.n_cigar; i++)
				cigar.emplace_back(ez.cigar[i]);
		}
		else if(type == KSW_ALN_left_extend){ // left extension, soft clip at end
			for(int i = ez.n_cigar - 1; i >= 0; i--)
				cigar.emplace_back(ez.cigar[i]);
		}else if(type == KSW_ALN_right_extend){ //right extension, soft clip at begin
			for(int i = 0; i < ez.n_cigar; i++)
				cigar.emplace_back(ez.cigar[i]);
		}
	}
}

