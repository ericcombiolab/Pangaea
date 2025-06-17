#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <htslib/sam.h>
#include "cmdline.h"
#include <dirent.h>


std::pair<std::string, std::string> get_seq_qual(bam1_t* bamAln, int readLength, bool isreverse)
{
	std::string readseq, readqual;
	uint8_t* seq = bam_get_seq(bamAln);
	uint8_t* qual = bam_get_qual(bamAln);
	if (!isreverse)
	{
		for (int i = 0; i < readLength; i++) {
			readseq += seq_nt16_str[bam_seqi(seq, i)];
		}
		for (int i = 0; i < readLength; i++) {
			readqual += (qual[i] + 33);
		}
	}
	else
	{
		for (int i = readLength - 1; i >= 0; i--) {
			char base = seq_nt16_str[bam_seqi(seq, i)];
			switch (base)
			{
			case 'A':
				readseq += 'T';
				break;
			case 'T':
				readseq += 'A';
				break;
			case 'C':
				readseq += 'G';
				break;
			case 'G':
				readseq += 'C';
				break;
			default:
				readseq += base;
				break;
			}
		}
		for (int i = readLength - 1; i >= 0; i--) {
			readqual += (qual[i] + 33);
		}
	}
	return std::make_pair(readseq, readqual);
}

void parse_reads(bam1_t* bamAln, bool& isRead1, bool& isRead2, int& readLength, std::string& read1seq, std::string& read1qual, std::string& read2seq, std::string& read2qual, bool isreverse)
{
	if (isRead1 && read1seq.empty())
	{
		auto read = std::move(get_seq_qual(bamAln, readLength, isreverse));
		read1seq = read.first;
		read1qual = read.second;
	}
	else if (isRead2 && read2seq.empty())
	{
		auto read = std::move(get_seq_qual(bamAln, readLength, isreverse));
		read2seq = read.first;
		read2qual = read.second;
	}
	return;
}

std::string parse_records(const double& idt, const int& min_l, int& read_unmapped, int& read_paired, std::vector<bam1_t*>& bamAln_vec, const std::unordered_set<std::string>& contigs, bam_hdr_t* bamHeader)
{
	std::string read1seq, read1qual, read2seq, read2qual;
	read_unmapped = 1;
	std::string readname;
	const int supplementaryAlignment = 0x00000800, secondaryAlignment = 0x00000100, isunmapped = 0x00000004;
	const int Read1 = 0x00000040, Read2 = 0x00000080, reverse = 0x00000010;
	for (auto&& bamAln : bamAln_vec)
	{
		int flag = bamAln->core.flag;
		bool isSupplementary = flag & supplementaryAlignment, isSecondary = flag & secondaryAlignment;
		bool isRead1 = flag & Read1, isRead2 = flag & Read2, record_unmapped = flag & isunmapped, isreverse = flag & reverse;
		int mappingQuality = bamAln->core.qual;
		int contigLength = bamHeader->target_len[bamAln->core.tid];
		int readLength = bamAln->core.l_qseq;
		int pos = bamAln->core.pos + 1, endPos = bam_endpos(bamAln) + 1;
		int reflen = sam_hdr_tid2len(bamHeader, bamAln->core.tid);
		readname = bam_get_qname(bamAln);
		parse_reads(bamAln, isRead1, isRead2, readLength, read1seq, read1qual, read2seq, read2qual, isreverse);
		// if (record_unmapped || pos <= 500 || endPos >= contigLength - 500)
		if (record_unmapped)
			continue;

		std::string contigname;
		if (sam_hdr_tid2name(bamHeader, bamAln->core.tid))
			contigname = bamHeader->target_name[bamAln->core.tid];
		if (contigname.empty() || contigs.find(contigname) == contigs.end())
			continue;

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

		int NM = 0;
		uint8_t* tmpNM = bam_aux_get(bamAln, "NM");
		if (tmpNM != NULL)
			NM = bam_aux2i(tmpNM);
		double blastIdentity = 1.0 * (alignmentColumns - NM) / alignmentColumns;
		if (blastIdentity < idt)
			continue;
		read_unmapped = 0;
		bam_destroy1(bamAln);
	}

	std::string ret;
	read_paired = 0;
	if (!read1seq.empty() && !read2seq.empty())
	{
		read_paired = 1;
		ret = "@" + readname + "/1\n" + read1seq + "\n+\n" + read1qual + "\n" + "@" + readname + "/2\n" + read2seq + "\n+\n" + read2qual + "\n";
	}
	else if (!read1seq.empty())
		ret = "@" + readname + "/1\n" + read1seq + "\n+\n" + read1qual + "\n";
	else if (!read2seq.empty())
		ret = "@" + readname + "/2\n" + read2seq + "\n+\n" + read2qual + "\n";
	return ret;
}

int main(int argc, char* argv[])
{
	cmdline::parser argParser;
	argParser.add<std::string>("bam", 'b', "Path to name sorted bam file.", true);
	argParser.add<std::string>("cov", 'c', "Path to contig coverage.", true);
	argParser.add<int>("int", 'f', "contig coverage cutoff", false, 100);
	argParser.add<double>("identity", 'i', "the threshold of identity", false, 0.95);
	argParser.add<int>("aligned_length", 'l', "the threshold of aligned bases", false, 60);
	argParser.add<std::string>("output", 'o', "Prefix to output.", true);
	argParser.parse_check(argc, argv);

	std::string output = argParser.get<std::string>("output");
	std::string cov = argParser.get<std::string>("cov");
	int cov_cutoff = argParser.get<int>("int");
	std::string bam = argParser.get<std::string>("bam");
	double idt = argParser.get<double>("identity");
	int min_l = argParser.get<int>("aligned_length");

	std::unordered_set <std::string> contigs;
	std::ifstream contig_cov(cov);
	std::string line;
	std::ofstream out_list(output + ".list");
	while (getline(contig_cov, line))
	{
		std::string name;
		double coverage = 0;
		size_t pos1, pos2 = 0;
		pos1 = line.find('\t');
		name = line.substr(0, pos1);
		if (!name.compare("contigName")) continue;
		pos2 = line.find('\t', pos1 + 1);
		pos1 = line.find('\t', pos2 + 1);
		coverage = std::stod(line.substr(pos2 + 1, pos1 - pos2 - 1));
		if (coverage >= cov_cutoff)
		{
			out_list << name << '\n';
			contigs.emplace(name);
		}
	}
	out_list.close();

	unsigned long long cnt_records = 0;
	samFile* bamFile = hts_open(bam.c_str(), "r");
	bam_hdr_t* bamHeader = sam_hdr_read(bamFile);
	bam1_t* bamAln = bam_init1();
	std::string last_readname;
	std::ofstream out_low_abd(output + ".low_abd.fq");

	const int openmp_size = 10000000;
	int openmp_index = 0;
	std::vector<bam1_t*> bamAln_tmp;
	std::vector<std::vector<bam1_t*>> bamAln_vec(openmp_size);
	std::vector<int> read_unmapped_vec(openmp_size);
	std::vector<int> read_paired_vec(openmp_size);
	std::vector<std::string> rets(openmp_size);

	while (sam_read1(bamFile, bamHeader, bamAln) >= 0)
	{
		if (cnt_records % 10000000 == 0)
			std::cout << "\rProcessed " << cnt_records << " alignment records. " << std::flush;
		cnt_records += 1;

		std::string current_readname;
		if (bam_get_qname(bamAln))
			current_readname = bam_get_qname(bamAln);
		if (current_readname.empty())
			continue;

		if (current_readname.compare(last_readname))
		{
			bamAln_vec[openmp_index++] = std::move(bamAln_tmp);
			bamAln_tmp.clear();
			if (openmp_index == openmp_size)
			{
#pragma omp parallel for
				for (int i = 0; i < openmp_index; ++i)
				{
					rets[i] = std::move(parse_records(idt, min_l, read_unmapped_vec[i], read_paired_vec[i], bamAln_vec[i], contigs, bamHeader));
				}
				for (int i = 0; i < openmp_index; ++i)
				{
					if (read_unmapped_vec[i])
					{
						if (read_paired_vec[i])
							out_low_abd << rets[i];
					}
				}
				openmp_index = 0;
				bamAln_vec.clear();
				bamAln_vec.resize(openmp_size);
				read_unmapped_vec.clear();
				read_unmapped_vec.resize(openmp_size);
				read_paired_vec.clear();
				read_paired_vec.resize(openmp_size);
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
		rets[i] = std::move(parse_records(idt, min_l, read_unmapped_vec[i], read_paired_vec[i], bamAln_vec[i], contigs, bamHeader));
	}
	for (int i = 0; i < openmp_index; ++i)
	{
		if (read_unmapped_vec[i])
		{
			if (read_paired_vec[i])
				out_low_abd << rets[i];
		}
	}
	std::cout << "\rProcessed " << cnt_records << " alignment records. " << std::endl;

	out_low_abd.close();
	bam_hdr_destroy(bamHeader);
	bam_destroy1(bamAln);
	sam_close(bamFile);

	return 0;
}
