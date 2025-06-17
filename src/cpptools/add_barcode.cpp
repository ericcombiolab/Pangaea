#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <htslib/sam.h>
#include "cmdline.h"
#include <dirent.h>


std::string parse_records(const int& min_l, int& read_paired, int& has_barcode, const std::unordered_map<std::string, std::string>& mapping, std::vector<bam1_t*>& bamAln_vec, bam_hdr_t* bamHeader)
{
	// std::string read1seq, read1qual, read2seq, read2qual;
	std::string readname;
	const int supplementaryAlignment = 0x00000800, secondaryAlignment = 0x00000100, isunmapped = 0x00000004;
	const int Read1 = 0x00000040, Read2 = 0x00000080, reverse = 0x00000010;
	// map query reads name to barcode

	std::vector<std::string> satisfied_barcodes;
	// find all possible barcodes
	for (auto&& bamAln : bamAln_vec)
	{
		int flag = bamAln->core.flag;
		int read_unmapped = 1;
		bool isSupplementary = flag & supplementaryAlignment, isSecondary = flag & secondaryAlignment;
		bool isRead1 = flag & Read1, isRead2 = flag & Read2, record_unmapped = flag & isunmapped, isreverse = flag & reverse;
		int mappingQuality = bamAln->core.qual;
		// or long read length
		int contigLength = bamHeader->target_len[bamAln->core.tid];
		int readLength = bamAln->core.l_qseq;
		int pos = bamAln->core.pos + 1, endPos = bam_endpos(bamAln) + 1;
		int reflen = sam_hdr_tid2len(bamHeader, bamAln->core.tid);
		readname = bam_get_qname(bamAln);
		

		// if (record_unmapped || pos <= 500 || endPos >= contigLength - 500)
		if (record_unmapped)
			continue;

		std::string contigname;
		if (sam_hdr_tid2name(bamHeader, bamAln->core.tid))
			contigname = bamHeader->target_name[bamAln->core.tid];
		if (contigname.empty())
			continue;
		if (mapping.find(contigname) == mapping.end()) {
			std::cout << "Find a long read without assigned barcode: " << contigname << std::endl;
			continue;
		}

		int alignmentColumns = 0;
		uint32_t* cigarPointer = bam_get_cigar(bamAln);
		for (int i = 0; i < bamAln->core.n_cigar; ++i)
		{
			char cigarOperator = bam_cigar_opchr(cigarPointer[i]);
			if (cigarOperator == 'M' || cigarOperator == 'I' || cigarOperator == 'D')
				alignmentColumns += bam_cigar_oplen(cigarPointer[i]);
		}
		if (alignmentColumns < min_l)
			continue;

		satisfied_barcodes.emplace_back(mapping.at(contigname));

		read_unmapped = 0;
		bam_destroy1(bamAln);
	}

	std::string ret;
	read_paired = 0;
	std::string barcode;
	has_barcode = 0;
	if (satisfied_barcodes.size() >= 1)
		barcode = satisfied_barcodes[rand() % satisfied_barcodes.size()];
	if (!readname.empty())
	{
		read_paired = 1;
		if (!barcode.empty()) {
			has_barcode = 1;
			ret = "paired" + readname + " " + barcode + "\n";
			// std::cout<<ret<<std::endl;
			// ret = "@" + readname + "\t" + barcode + "\n" + read1seq + "\n+\n" + read1qual + "\nPaired@" + readname + "\t" + barcode + "\n" + read2seq + "\n+\n" + read2qual + "\n";
		}
		else
			{
				ret = "paired" + readname + "\n";
				// std::cout<<ret<<std::endl;
			}
			// ret = "@" + readname + "\n" + read1seq + "\n+\n" + read1qual + "\nPaired@" + readname + "\n" + read2seq + "\n+\n" + read2qual + "\n";
	}
	// else if (!read1seq.empty())
	// 	ret = "@" + readname + "/1\n" + read1seq + "\n+\n" + read1qual + "\n";
	// else if (!read2seq.empty())
	// 	ret = "@" + readname + "/2\n" + read2seq + "\n+\n" + read2qual + "\n";
	return ret;
}

int main(int argc, char* argv[])
{
	cmdline::parser argParser;
	argParser.add<std::string>("bam", 'b', "path to name sorted (samtools sort -n) bam file.", true);
	argParser.add<std::string>("map", 'm', "path to long read-barcode mapping file.", true);
	argParser.add<int>("aligned_length", 'l', "the threshold of aligned bases", false, 60);
	argParser.add<std::string>("output", 'o', "prefix to output.", true);
	argParser.parse_check(argc, argv);

	srand(2023);
	std::string output = argParser.get<std::string>("output");
	std::string mapping = argParser.get<std::string>("map");
	std::string bam = argParser.get<std::string>("bam");
	int min_l = argParser.get<int>("aligned_length");

	std::cout << "Start loading long read to barcode mapping." << std::endl;
	std::unordered_map<std::string, std::string> lr_to_barcode;
	std::ifstream mapping_f(mapping.c_str());
	std::string line;
	while (getline(mapping_f, line))
	{
		std::string name, barcode;
		double coverage = 0;
		size_t pos = 0;
		pos = line.find(' ');
		name = line.substr(0, pos);
		barcode = line.substr(pos + 1);
		if (lr_to_barcode.find(name) != lr_to_barcode.end())
			std::cout << "Barcode for long read " << name << " is already in the list." << std::endl;
		else
			lr_to_barcode[name] = barcode;
	}
	std::cout << "Print the first 10 samples for checking." << std::endl;
	int cnt = 0;
	for (auto&& i : lr_to_barcode) {
		std::cout << i.first << " --> " << i.second << std::endl;
		if (++cnt == 10)
			break;
	}

	std::cout << "Start reading bam file (must be name sorted)." << std::endl;
	unsigned long long cnt_records = 0;
	samFile* bamFile = hts_open(bam.c_str(), "r");
	bam_hdr_t* bamHeader = sam_hdr_read(bamFile);
	bam1_t* bamAln = bam_init1();
	std::string last_readname;
	std::ofstream out_barcoded_map(output + "_map.txt");

	// initialize openmp and bam reader
	const int openmp_size = 10000000;
	int openmp_index = 0;
	int cnt_unpaired = 0, cnt_barcode = 0, cnt_no_barcode = 0;
	std::vector<bam1_t*> bamAln_tmp;
	std::vector<std::vector<bam1_t*>> bamAln_vec(openmp_size);
	std::vector<int> read_paired_vec(openmp_size);
	std::vector<int> has_barcode_vec(openmp_size);
	std::vector<std::string> rets(openmp_size);

	// start reading bam line by line
	while (sam_read1(bamFile, bamHeader, bamAln) >= 0)
	{
		// logging lines processed
		if (cnt_records % 10000000 == 0)
			std::cout << "\rProcessed " << cnt_records << " alignment records. " << std::flush;
		cnt_records += 1;

		std::string current_readname;
		if (bam_get_qname(bamAln))
			current_readname = bam_get_qname(bamAln);
		if (current_readname.empty())
			continue;

		// only process the records when all the alignments for the read are loaded in the memory
		if (current_readname.compare(last_readname))
		{
			// store alignment records and process them together in parallel by openmp
			bamAln_vec[openmp_index++] = std::move(bamAln_tmp);
			bamAln_tmp.clear();
			if (openmp_index == openmp_size)
			{
#pragma omp parallel for
				for (int i = 0; i < openmp_index; ++i)
				{
					rets[i] = std::move(parse_records(min_l, read_paired_vec[i], has_barcode_vec[i], lr_to_barcode, bamAln_vec[i], bamHeader));
				}
				for (int i = 0; i < openmp_index; ++i)
				{
					if (read_paired_vec[i]) {
						std::string readmap;
						readmap = rets[i].substr(rets[i].find("paired") + 6);
						out_barcoded_map << readmap;
						if (has_barcode_vec[i])
							cnt_barcode += 1;
						else cnt_no_barcode += 1;
					}
					else cnt_unpaired += 1;
				}
				openmp_index = 0;
				bamAln_vec.clear();
				bamAln_vec.resize(openmp_size);
				read_paired_vec.clear();
				read_paired_vec.resize(openmp_size);
				has_barcode_vec.clear();
				has_barcode_vec.resize(openmp_size);
				rets.clear();
				rets.resize(openmp_size);
			}
			last_readname = current_readname;
		}
		bamAln_tmp.emplace_back(bamAln);
		bamAln = bam_init1();
	}

	bamAln_vec[openmp_index++] = std::move(bamAln_tmp);
	bamAln_tmp.clear();
#pragma omp parallel for
	for (int i = 0; i < openmp_index; ++i)
	{
		rets[i] = std::move(parse_records(min_l, read_paired_vec[i], has_barcode_vec[i], lr_to_barcode, bamAln_vec[i], bamHeader));
	}
	for (int i = 0; i < openmp_index; ++i)
	{
		if (read_paired_vec[i]) {
			std::string readmap;
			readmap = rets[i].substr(rets[i].find("paired") + 6);
			out_barcoded_map << readmap;
			if (has_barcode_vec[i])
				cnt_barcode += 1;
			else cnt_no_barcode += 1;
		}
		else cnt_unpaired += 1;
	}
	std::cout << "\rProcessed " << cnt_records << " alignment records. " << std::endl;

	std::cout << "With barcode: " << cnt_barcode << "; without barcode: " << cnt_no_barcode << "; unpaired (not outputted): " << cnt_unpaired << std::endl;
	out_barcoded_map.close();
	bam_hdr_destroy(bamHeader);
	bam_destroy1(bamAln);
	sam_close(bamFile);

	return 0;
}
