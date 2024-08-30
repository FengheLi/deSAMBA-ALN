#include<stdio.h>
#include<string.h>
#include <time.h>

#include "Aln_core/align.hpp"
#include"cpp_lib/get_option_cpp.hpp"
extern "C"
{
#include "clib/utils.h"
#include "clib/desc.h"
}

int classify_main(int argc, char *argv[])
{
	fprintf(stderr, "\n\n V1.21\n\n");
	//get option
	deCOY_CLASSIFY_MAIN cm;
	cm.init_run(argc, argv);
	return 0;
}

int main(int argc, char *argv[])
{
	fprintf(stderr, "\n\n main version V2.00\n\n");
	COMMAND_HANDLER ch;
	ch.add_function("fc_index", "The 2ed step of force calling (S2)", build_deBGA_index);
	ch.add_help_msg_back("Getting signal reads from bam files");
	ch.add_function("fc_aln", "The 4th step of force calling (S4)", classify_main);
	ch.add_help_msg_back("Aligning all signal reads into anchor reference");

	return ch.run(argc, argv);
}
