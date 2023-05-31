#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>
#include "cmdline.h"
#include "gzstream.h"
#include "util.h"

int main(int argc, char* argv[])
{
    cmdline::parser argParser;
    argParser.add<std::string>("reads1", '1', "Path to reads1 (gzipped and barcode sorted).", false);
    argParser.add<std::string>("reads2", '2', "Path to reads2 (gzipped and barcode sorted).", false);
    argParser.add<std::string>("interleaved", 'i', "Path to interleaved reads (gzipped and barcode sorted).", false);
    argParser.add<std::string>("long_reads", 'r', "Path to longreads (gzipped).", false);

    argParser.add<std::string>("clusters", 'c', "Path to clusters tsv.", true, "");
    argParser.add<std::string>("output", 'o', "Prefix to output.", true, "");
    argParser.parse_check(argc, argv);

    std::string reads1 = argParser.get<std::string>("reads1");
    std::string reads2 = argParser.get<std::string>("reads2");
    std::string interleaved = argParser.get<std::string>("interleaved");
    std::string long_reads = argParser.get<std::string>("long_reads");
    std::string clusters = argParser.get<std::string>("clusters");
    std::string output = argParser.get<std::string>("output");

    std::unordered_map<std::string, unsigned int> barcode2cluster;
    std::fstream clustersf(clusters.c_str(), std::fstream::in);
    // std::vector<std::fstream> output_reads1, output_reads2, output_barcode, output_fq;
    std::vector<std::fstream> output_barcode, output_fq;
    std::string line;
    unsigned int cluster_id = 0;
    while (getline(clustersf, line))
    {
        std::size_t pos = line.find_first_of('\t');
        if (line.substr(0, pos).compare("-1") == 0)
            continue;
        // std::fstream f1((output + "_bin" + line.substr(0, pos) + "_1.fq").c_str(), std::fstream::out);
        // std::fstream f2((output + "_bin" + line.substr(0, pos) + "_2.fq").c_str(), std::fstream::out);
        std::fstream f3((output + "_bin" + line.substr(0, pos) + ".barcode").c_str(), std::fstream::out);
        std::fstream f4((output + "_bin" + line.substr(0, pos) + ".fq").c_str(), std::fstream::out);
        // output_reads1.push_back(std::move(f1));
        // output_reads2.push_back(std::move(f2));
        output_barcode.push_back(std::move(f3));
        output_fq.push_back(std::move(f4));
        std::string barcode;
        while (pos != line.size())
        {
            while (++pos < line.size() && line.at(pos) != ',')
                barcode += line.at(pos);
            barcode2cluster[barcode] = cluster_id;
            barcode.clear();
        }
        cluster_id += 1;
    }
    clustersf.close();

    if (!interleaved.empty())
    {
        std::cout << "interleaved file: "<<interleaved.c_str()<<std::endl;
        igzstream readsf(interleaved.c_str());
        unsigned long long cnt_line = 0;
        bool isSatisfied = false;
        std::string read_line, reads_seq;
        std::pair<std::string, std::string> p;

        while (getline(readsf, read_line))
        {
            switch (++cnt_line % 8)
            {
            case 1:
                p = std::move(getBarcode(read_line));
                if (barcode2cluster.find(p.second) != barcode2cluster.end()) {
                    isSatisfied = true;
                    reads_seq += (p.first + "\tBX:Z:" + p.second + "-1\n");
                } 
                else
                    isSatisfied = false;
                break;
            case 0:
                if (cnt_line % 8000000 == 0)
                    std::cout << "\rProcessed " << cnt_line / 8 << " read pairs. " << std::flush;
                if (isSatisfied)
                {
                    reads_seq += (read_line + '\n');
                    output_barcode.at(barcode2cluster.at(p.second)) << p.second << '\n';
                    output_fq.at(barcode2cluster.at(p.second)) << reads_seq;
                    reads_seq.clear();
                }
                break;
            default:
                if (isSatisfied)
                {
                    reads_seq += (read_line + '\n');
                    break;
                }
            }
        }
    }
    else if (!long_reads.empty()) {
        igzstream readsf(long_reads.c_str());
        unsigned long long cnt_line = 0;
        bool isSatisfied = false;
        std::string read_line, reads_seq;
        std::pair<std::string, std::string> p;

        while (getline(readsf, read_line))
        {
            switch (++cnt_line % 4)
            {
            case 1:
                p = std::move(getBarcode(read_line));
                if (barcode2cluster.find(p.second) != barcode2cluster.end()) {
                    isSatisfied = true;
                    reads_seq += (p.first + "\tBX:Z:" + p.second + "\n");
                } 
                else
                    isSatisfied = false;
                break;
            case 0:
                if (cnt_line % 4000000 == 0)
                    std::cout << "\rProcessed " << cnt_line / 4 << " read pairs. " << std::flush;
                if (isSatisfied)
                {
                    reads_seq += (read_line + '\n');
                    output_barcode.at(barcode2cluster.at(p.second)) << p.second << '\n';
                    output_fq.at(barcode2cluster.at(p.second)) << reads_seq;
                    reads_seq.clear();
                }
                break;
            default:
                if (isSatisfied)
                {
                    reads_seq += (read_line + '\n');
                    break;
                }
            }
        }
    }
    else if (!reads1.empty() && !reads2.empty()) {
        igzstream reads1f(reads1.c_str()), reads2f(reads2.c_str());
        unsigned long long cnt_line = 0;
        bool isSatisfied = false;
        std::string line1, line2, reads_seq1, reads_seq2;
        std::pair<std::string, std::string> p1, p2;

        while (getline(reads1f, line1))
        {
            getline(reads2f, line2);
            switch (++cnt_line % 4)
            {
            case 1:
                p1 = std::move(getBarcode(line1)), p2 = std::move(getBarcode(line2));
                if (barcode2cluster.find(p1.second) != barcode2cluster.end() && !p1.first.compare(p2.first) && !p1.second.compare(p2.second))
                {
                    isSatisfied = true;
                    reads_seq1 += (p1.first + "\tBX:Z:" + p1.second + "-1\n");
                    reads_seq2 += (p2.first + "\tBX:Z:" + p2.second + "-1\n");
                    // reads_seq1 += (p1.first + "/1\n");
                    // reads_seq2 += (p2.first + "/2\n");
                }
                else
                    isSatisfied = false;
                break;
            case 0:
                if (cnt_line % 4000000 == 0)
                    std::cout << "\rProcessed " << cnt_line / 4 << " read pairs. " << std::flush;
                if (isSatisfied)
                {
                    reads_seq1 += (line1 + '\n');
                    reads_seq2 += (line2 + '\n');
                    // output_reads1.at(barcode2cluster.at(p1.second)) << reads_seq1;
                    // output_reads2.at(barcode2cluster.at(p1.second)) << reads_seq2;
                    output_barcode.at(barcode2cluster.at(p1.second)) << p1.second << '\n';
                    output_fq.at(barcode2cluster.at(p1.second)) << reads_seq1 << reads_seq2;
                    reads_seq1.clear();
                    reads_seq2.clear();
                }
                break;
            default:
                if (isSatisfied)
                {
                    reads_seq1 += (line1 + '\n');
                    reads_seq2 += (line2 + '\n');
                    break;
                }
            }
        }
        reads1f.close();
        reads2f.close();
        if (cnt_line % 4000000 == 0)
            std::cout << std::endl;
        else
            std::cout << "\rProcessed " << cnt_line / 4 << " read pairs. " << std::endl;
    }

    // for (auto &&f : output_reads1)
    //     f.close();
    // for (auto &&f : output_reads2)
    //     f.close();
    for (auto&& f : output_barcode)
        f.close();
    for (auto&& f : output_fq)
        f.close();
    return 0;
}