#include "../Aln_core/deBGA_index.hpp"

#include <cstring>
extern "C"
{
#include "../clib/utils.h"
}


#define POS_N_MAX 500
#define POS_N_MAX_LEVEL2 8000
#define RANDOM_NUM 500

#define KMER_LEN_FIRST_LEVEL 14
#define UPPER_PART_SHIFT ((LEN_KMER - KMER_LEN_FIRST_LEVEL) << 1) //（20 - 14）X 2 = 12
#define LOWER_PART_MASK 0X3fffff //k = 22

int binsearch_range(uint64_t key, uint32_t *v, int64_t n,  int64_t *range, int8_t k_off)
{
    int64_t l=0, r=n-1, m;
    uint32_t tmp = 0;
    range[0] = range[1] = -1;
    if (k_off == 0){
        while(l <= r){
            m = (l + r)/2;
            if (key < v[m])       	r = m - 1;
            else if (key > v[m])  	l = m + 1;
            else                 	{ range[0] = range[1] = m; return 1;}
        }
    }
    else {
        while (l <= r)  {
            m = (l+r)/2;
            tmp = v[m] >> k_off;
            if (tmp == key){
                range[0] = range[1] = m;
                //run low bound
                int64_t sl=l, sr=m-1, sm;
                while (sl <= sr)
                {
                    sm = (sl+sr)/2;
                    tmp = v[sm] >> k_off;
                    if (tmp == key) 		{ range[0] = sm; sr = sm-1; }
                    else if (tmp > key) 	sr = sm - 1;
                    else    				sl = sm + 1;
                }

                //run upper bound
                sl = m+1; sr = r;
                while (sl <= sr)
                {
                    sm = (sl+sr)/2;
                    tmp = v[sm] >> k_off;
                    if (tmp == key) 		{ range[1] = sm; sl = sm+1; }
                    else if (tmp > key) 	sr = sm - 1;
                    else    				sl = sm + 1;
                }
                return 1;
            }
            else if (tmp > key) r = m - 1;
            else l = m + 1;
        }
    }
    return -1;
}

int64_t binsearch_interval_unipath64(uint64_t x, uint64_t v[], uint64_t n)
{
    int64_t low = 0, high = n - 1, mid = 0;
    while ( low <= high ) {
        mid = ((int64_t )(low + high)) >> 1;
        if(x < v[mid])            high = mid - 1;
        else if(x > v[mid])       low = mid + 1;
        else  /*found match*/     return mid;
    }
    return high;
}

#define ROUTE_LENGTH_MAX 1024
//load from file, return the true load data size(in byte)
uint64_t load_data_from_file(const char * path_name, const char * fn, void ** data, bool useCalloc, int additional_byte){
	char full_fn[ROUTE_LENGTH_MAX] = {0};
	strcpy(full_fn, path_name);
	if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
	strcat(full_fn, fn);
	FILE *fp_us_b = xopen (full_fn, "rb" );

	fseek(fp_us_b, 0, SEEK_END);// non-portable
	uint64_t file_size = ftell(fp_us_b);
	rewind(fp_us_b);

	if(useCalloc)
		(*data) = (uint64_t* ) xcalloc (file_size + additional_byte, 1);
	else
		(*data) = (uint64_t* ) xmalloc (file_size + additional_byte);

	xread ((*data), 1, file_size, fp_us_b);
	fclose(fp_us_b);
	return file_size;
}

#define MAX_CHR_NAME_LENGTH 200
int deBGA_INDEX::load_index_file(char *index_dir)
{
    //        buffer_ref_seq:store reference seq;load from dm file ref.seq
    result_ref_seq = (load_data_from_file(index_dir, "ref.seq", ((void ** )&buffer_ref_seq), true, 536) >> 3);
    //        buffer_seq:store unipath seq(concatenate all unipaths's seq);load from dm file unipath.seqb
    result_seq = (load_data_from_file(index_dir, "unipath.seqb", ((void ** )&buffer_seq), false, 0) >> 3);
//        buffer_seqf:store unipath's offset on unipath seq;load from dm file unipath.seqfb
    result_seqf = (load_data_from_file(index_dir, "unipath.seqfb", ((void ** )&buffer_seqf), false, 0) >> 3);
    //buffer_p: store every unipath's set of positions on reference;load from dm file unipath.pos
    result_p = (load_data_from_file(index_dir,  "unipath.pos", ((void ** )&buffer_p), false, 0) >> 3);
    //buffer_pp: unipath's pointer to array buffer_p;load from dm file unipath.posp
    result_pp = (load_data_from_file(index_dir,   "unipath.posp", ((void ** )&buffer_pp), false, 0) >> 3);
    //buffer_hash_g: store kmer's hash part; load from dm file unipath_g.hash
    result_hash_g = (load_data_from_file(index_dir,  "unipath_g.hash", ((void ** )&buffer_hash_g), false, 0) >> 3);
// buffer_kmer_g:store kmer's kmer part;load from dm file unipath_g.kmer
    result_kmer_g = (load_data_from_file(index_dir,  "unipath_g.kmer", ((void ** )&buffer_kmer_g), false, 0) >> 2);
//buffer_off_g: store kmer's offset on unipath seq; load from dm file unipath_g.offset
    result_off_g = (load_data_from_file(index_dir,  "unipath_g.offset", ((void ** )&buffer_off_g), false, 0) >> 2);

    //********************************************************************************************
    //read chr names and length from unipath.chr
    char unichr[ROUTE_LENGTH_MAX] = {0};
    strcpy(unichr, index_dir);
    strcat(unichr, "unipath.chr");

    FILE *fp_chr = xopen (unichr, "r" );
    uint32_t chr_line_n = 0;
    char chr_line_content[MAX_CHR_NAME_LENGTH];
    fscanf(fp_chr,"%s",chr_line_content);
    while(!feof(fp_chr))
    {
        if ((chr_line_n & 0X1) == 0)	{
        	chr_names.emplace_back(chr_line_content);
        }
        else{
        	uint32_t cur_chr_end_position = 0;
        	sscanf(chr_line_content, "%u", &cur_chr_end_position);
        	chr_end_position.emplace_back(cur_chr_end_position);
        }

        fflush(stdout);
        chr_line_n++;
        fscanf(fp_chr,"%s",chr_line_content);
    }
    chr_names.emplace_back("*");
    reference_len = chr_end_position.back();
    fclose(fp_chr);

    //building search index for CHR
    building_chr_index();
    building_bam_header();

    return 0;
}

//search kmer in the index, the search result "index of kmer" stored in range, range[0] is the start index, and range[1] is the end index
//return false when search failed, otherwise return 0
bool deBGA_INDEX::search_kmer(int len_k, uint64_t kmer, int64_t range[2], int8_t seed_offset){
	if (len_k == KMER_LEN_FIRST_LEVEL){ //k = 14
		range[0] = buffer_hash_g[kmer];
		range[1] = buffer_hash_g[kmer + 1] - 1;
		if (range[1] < range[0])
			return false;
	}
	else{
		uint64_t seed_kmer = (kmer & LOWER_PART_MASK);
		uint64_t seed_hash = (kmer >> UPPER_PART_SHIFT);

		int result = binsearch_range(seed_kmer, buffer_kmer_g + buffer_hash_g[seed_hash],
				buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], range, seed_offset<<1);
		if (result == -1)	return false;
		range[0] += buffer_hash_g[seed_hash];
		range[1] += buffer_hash_g[seed_hash];
	}
	return true;
}

//function: given an index of kmer, search the MEM within a UNITIG
//return 0 when successfully stored data, -1 when failed
int deBGA_INDEX::MEM_extern_within_UNITIG_for_KMER(uint64_t kmer_index, std::vector<MEM_Within_UNITIG>& mem_search_rst,
		uint64_t *read_bit, uint32_t read_off, uint32_t read_length, int len_k, uint32_t & max_right_i){

	uint64_t kmer_pos_uni = buffer_off_g[kmer_index];//this kmer's offset on unipath seq
	//find the UID of this kmer
	int64_t seed_id_r = binsearch_interval_unipath64(kmer_pos_uni, buffer_seqf, result_seqf);
	uint64_t ref_pos_n = buffer_pp[seed_id_r + 1] - buffer_pp[seed_id_r];

	uint32_t uni_offset_s_l = kmer_pos_uni - buffer_seqf[seed_id_r];
	uint32_t uni_offset_s_r = buffer_seqf[seed_id_r + 1] - (kmer_pos_uni + len_k);
	//extend the kmer to a exact match
	uint32_t left_i, right_i;
	//extend left
	for(left_i = 1; (left_i <= uni_offset_s_l) && (left_i <= read_off); left_i++){
		uint8_t ref_base = ((buffer_seq[(kmer_pos_uni - left_i) >> 5] >> ((31 - ((kmer_pos_uni - left_i) & 0X1f)) << 1)) & 0X3);
		uint8_t read_base = ((read_bit[(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3);
		if(ref_base != read_base)
			break;
	}
	//extend right
	for(right_i = 1; (right_i <= uni_offset_s_r) && (right_i <= read_length - read_off - len_k); right_i++)
	{
		uint8_t ref_base = ((buffer_seq[(kmer_pos_uni + len_k - 1 + right_i) >> 5] >> ((31 - ((kmer_pos_uni + len_k - 1 + right_i) & 0X1f)) << 1)) & 0X3);
		uint8_t read_base = ((read_bit[(read_off + len_k - 1 + right_i) >> 5] >> ((31 - ((read_off + len_k - 1 + right_i) & 0X1f)) << 1)) & 0X3);
		if(ref_base != read_base)
			break;
	}

	uint32_t read_pos = read_off + 1 - left_i;
	uint32_t mem_length = len_k + left_i + right_i - 2;

	mem_search_rst.emplace_back();
	MEM_Within_UNITIG & vertexm = mem_search_rst.back();
	vertexm.uid = seed_id_r;
	vertexm.seed_id = mem_search_rst.size() - 1;
	vertexm.read_pos = read_pos;
	vertexm.uni_pos_off = uni_offset_s_l + 1 - left_i;
	vertexm.length = mem_length;
	vertexm.pos_n = ref_pos_n;

	if (right_i > max_right_i)
		max_right_i = right_i;

	return 0;
}

#define waitingLen 3 //at most 2 bases in gap, disable multy-SNP
#define	Eindel 1// disable INDLE in SEED merging
//function: merge all MEMs within a unipath, the result stored in vertexu_v
void deBGA_INDEX::merge_seed_in_unipath(std::vector<MEM_Within_UNITIG>& MWU_l, std::vector<Chain_Within_Unitig>& CWU_l){
	uint32_t mem_i = MWU_l.size();
	//merge seeds in the same unipath
	uint32_t s1, e1;
	if(mem_i == 0){} // do nothing
	else if (mem_i == 1)
	{
		CWU_l.emplace_back();
		Chain_Within_Unitig & vertexu = CWU_l.back();
		MEM_Within_UNITIG & vertexm = MWU_l.back();
		vertexu.uid = vertexm.uid;
		vertexu.read_offset = vertexm.read_pos;
		vertexu.unitig_offset = vertexm.uni_pos_off; //whether to set uni_pos_off to the leftest position
		vertexu.ref_pos_n = vertexm.pos_n;
		vertexu.length_in_read = vertexm.length;
		vertexu.length_in_ref = vertexm.length;
		vertexu.cov = vertexm.length;
	}else{
		qsort(&(MWU_l[0]), mem_i, sizeof(MEM_Within_UNITIG), MEM_Within_UNITIG::cmp);
		uint64_t uni_id_temp = MWU_l[0].uid;
		uint32_t j = 0;
		uint32_t cov = 0;
		while (j < mem_i)
		{
			s1 = j;
			cov = MWU_l[s1].length;
			j++;
			while ((uni_id_temp == MWU_l[j].uid) && (MWU_l[j].uni_pos_off > MWU_l[j-1].uni_pos_off) && (j < mem_i))
			{
				int diff = (int)(MWU_l[j].read_pos - MWU_l[j-1].read_pos - MWU_l[j-1].length);
				if (diff > waitingLen)
					break;
				int c_eindel = (MWU_l[j].uni_pos_off - MWU_l[j-1].uni_pos_off) - (MWU_l[j].read_pos - MWU_l[j-1].read_pos);
				if (std::abs(c_eindel) < Eindel)
				{
					cov += (diff > 0)? MWU_l[j].length : (diff + MWU_l[j].length);
					++j;
				}
				else
					break;
			}
			e1 = j - 1;
			//store vertexu_v
			CWU_l.emplace_back();
			Chain_Within_Unitig & vertexu = CWU_l.back();
			vertexu.uid = MWU_l[s1].uid;
			vertexu.read_offset = MWU_l[s1].read_pos;
			vertexu.unitig_offset = MWU_l[s1].uni_pos_off; //whether to set uni_pos_off to the leftest position
			vertexu.ref_pos_n = MWU_l[s1].pos_n;
			vertexu.cov = cov;
			cov = 0;

			//cal length
			if (s1 == e1)
			{
				vertexu.length_in_read = MWU_l[s1].length;
				vertexu.length_in_ref = MWU_l[s1].length;
			}else{
				vertexu.length_in_read = MWU_l[e1].read_pos + MWU_l[e1].length - MWU_l[s1].read_pos;
				vertexu.length_in_ref = MWU_l[e1].uni_pos_off + MWU_l[e1].length - MWU_l[s1].uni_pos_off;
			}

			uni_id_temp = MWU_l[j].uid;
		}
	}
}

void deBGA_INDEX::get_refseq(uint8_t *ref, uint32_t len, uint32_t start)
{
	uint32_t m;
    for (m = 0; m < len; ++m)
    {
        ref[m] = (buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3;
    }
}

void deBGA_INDEX::get_refseq_char_bin(std::vector<uint8_t> &ref, uint32_t start, uint32_t len){
    for (uint32_t m = 0; m < len; ++m)
    	ref.emplace_back((buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3);
}

void deBGA_INDEX::get_refseq_char(std::vector<char> &ref, uint32_t start, uint32_t len){
    for (uint32_t m = 0; m < len; ++m)
    	ref.emplace_back("ACGTN"[(buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3]);
    ref.emplace_back(0);
}

void deBGA_INDEX::get_ref_jump(uint8_t *ref, uint32_t start1, uint32_t len1, uint32_t start2, uint32_t len2)
{
	uint32_t m;
	uint32_t k = 0;

    for (m = 0; m < len1; ++m)
    {
        ref[k++] = (buffer_ref_seq[(m + start1) >> 5] >> ((31 - ((m + start1) & 0X1f)) << 1)) & 0X3;
    }

    for (m = 0; m < len2; ++m)
    {
        ref[k++] = (buffer_ref_seq[(m + start2) >> 5] >> ((31 - ((m + start2) & 0X1f)) << 1)) & 0X3;
    }
}

void deBGA_INDEX::get_ref_onebyone(uint8_t *ref, uint32_t start, uint32_t len, uint32_t pre_pos)
{
	uint32_t m;
	uint32_t k = pre_pos;
    for (m = 0; m < len; ++m)
    {
        ref[k++] = (buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3;
    }
}

void deBGA_INDEX::free_memory(){
	if(buffer_ref_seq 	!= NULL){ free(buffer_ref_seq); buffer_ref_seq = NULL; 	} //original reference
	if(buffer_seq 		!= NULL){ free(buffer_seq); 	buffer_seq = NULL; 		}// UNITIG sequence
	if(buffer_seqf != NULL)		{ free(buffer_seqf); 	buffer_seqf = NULL; }//start offset of [UINTIG] in [buffer_seq]
	if(buffer_off_g != NULL)	{ free(buffer_off_g); 	buffer_off_g = NULL; }//start offset of [KMER] in [buffer_seq]
	if(buffer_p != NULL)		{ free(buffer_p); 		buffer_p = NULL; }//start offset of [UINTIG] in [buffer_ref_seq]
	if(buffer_pp != NULL)		{ free(buffer_pp); 		buffer_pp = NULL; }//start index of [UINTIG] in [buffer_p], used to search offset of [UINTIG] in [buffer_ref_seq]
	if(buffer_hash_g != NULL)	{ free(buffer_hash_g); 	buffer_hash_g = NULL; }//index for the first 14 base pair, used to search kmer in the hash table
	if(buffer_kmer_g != NULL)	{ free(buffer_kmer_g); 	buffer_kmer_g = NULL; }//used to store other part of a kmer(except the first 14 byte)
}

#define CHR_SEARCH_IDX_STEP 0x4000 //16384 = 2 ^ 14
void deBGA_INDEX::building_chr_index(){
	chr_search_index_size = (reference_len >> 14) + 2;
	chr_search_index = (uint32_t *)xcalloc(chr_search_index_size, sizeof(uint32_t));
	uint32_t pos_index_size = 0;
	for(uint32_t i = 0; i < chr_end_position.size(); i++){
		uint32_t pos_index = chr_end_position[i] / CHR_SEARCH_IDX_STEP;
		while(pos_index >= pos_index_size){
			chr_search_index[pos_index_size++] = i;
		}
	}
	chr_search_index[pos_index_size] = chr_end_position.size();// store final block
}

//get the chr_ID that a 'position' belong to.
int deBGA_INDEX::get_chromosome_ID(uint64_t position_global)
{
	int chr_id_rst = 0;
	uint64_t pos_index = position_global / CHR_SEARCH_IDX_STEP;
	int low = chr_search_index[pos_index];
	int high = chr_search_index[pos_index + 1];
	int mid;
	uint64_t pos = position_global + 1;

	while ( low <= high )
	{
		mid = (low + high) >> 1;
		if(pos < (chr_end_position[mid] - 1))
		{
			high = mid - 1;
		}
		else if(pos > (chr_end_position[mid] - 1))
		{
			low = mid + 1;
		}
		else
		{
			return mid;
		}
		chr_id_rst = low;
	}
	return chr_id_rst;
}

void deBGA_INDEX::global_position_to_chrID_and_position_in_chr(uint64_t position_global, int & chrID, int &position_in_chr){
	chrID = get_chromosome_ID(position_global);
	if(chrID == 0)
		position_in_chr = position_global - UNITIG_MAGIC_NUM;
	else
		position_in_chr = position_global - chr_end_position[chrID - 1] - UNITIG_MAGIC_NUM;
}

void deBGA_INDEX::chrID_and_position_in_chr_to_global_position(uint64_t &position_global, int chrID, int position_in_chr){
	if(chrID == 0)
		position_global = position_in_chr + 0 + UNITIG_MAGIC_NUM;
	else
		position_global = position_in_chr + chr_end_position[chrID - 1] + UNITIG_MAGIC_NUM;
}

void deBGA_INDEX::building_bam_header(){
	//todo::void

}

int deBGA_index_usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:	%s fc_index <ref.fa> <index_route>\n", PACKAGE_NAME);
	fprintf(stderr, "		build deBGA index file with default 22-kmer. You can get more deBGA information from https://github.com/HongzheGuo/deBGA");
	fprintf(stderr, "\n");
	return 1;
}

int get_bin_dir(char *bin, char *dir)
{
	char *end = strrchr(bin, '/');
	if (end == NULL)
		return 1;

	bin[end-bin] = '\0';
	strcpy(dir, bin);
	return 0;
}

int build_deBGA_index_core(char *dir, char *ref_fa, char *index_route)
{
	char cmd[1024];
	sprintf(cmd, "%sdeBGA index %s %s", dir, ref_fa, index_route);
	fprintf(stderr, "[deBGA_index] Executing deBGA index ...\n");
	if (system(cmd) != 0)
	{
		fprintf(stderr, "\n[deBGA_index] Indexing undoing, deBGA index exit abnormally. \n");
		exit(1);
	}
	fprintf(stderr, "[deBGA_index] Done!\n");

    return 1;
}

int build_deBGA_index(int argc, char *argv[])
{
	// clock_t t = clock();
	char *ref_fa = 0;
	char *index_route = 0;

    char bin_dir[1024];
    char panSVR_path[1024];
    int r;

	if(optind + 2 > argc)	return deBGA_index_usage();

	ref_fa = strdup(argv[optind]);
	index_route = strdup(argv[optind + 1]);

    r = readlink("/proc/self/exe", bin_dir, 2048);
    if (r < 0 || r >= 2048)
    {
        fprintf(stderr, "Failed, could not find the program of deBGA!\n");
    }
    bin_dir[r] = '\0';

    if (!get_bin_dir(bin_dir, panSVR_path))
        strcat(panSVR_path, "/");

    char path_deBGA[1024];
    strcpy(path_deBGA, panSVR_path);
    strcat(path_deBGA, "deBGA");

    if((access(path_deBGA, F_OK)) == -1)
    {
        fprintf(stderr, "[Wrong!] %s not exist, please check!\n", path_deBGA);
        exit(1);
    }

	char cmd[1024];
	sprintf(cmd, "%sdeBGA index %s %s", panSVR_path, ref_fa, index_route);
	fprintf(stderr, "[deBGA_index] Executing deBGA index ...\n");
	if (system(cmd) != 0)
	{
		fprintf(stderr, "\n[deBGA_index] Indexing undoing, deBGA index exit abnormally. \n");
		exit(1);
	}
	fprintf(stderr, "[deBGA_index] Done!\n");

	return 0;
}

