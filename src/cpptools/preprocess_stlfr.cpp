#include "cmdline.h"
#include "gzstream.h"
#include <unordered_map>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>

int main(int argc, char* argv[])
{
    cmdline::parser argParser;
    argParser.add<std::string>("reads1", '1', "Path to reads1 (gzipped).", true);
    argParser.add<std::string>("reads2", '2', "Path to reads2 (gzipped).", true);
    argParser.add("number", 'n', "Use barcode index numbers.");
    argParser.add("library", 'l', "Add library after barcode.");
    argParser.add<std::string>("whitelist", 'w', "Barcode whitelist (needed if -n is not set).", false, "/path/to/whitelist");
    argParser.add<std::string>("output", 'o', "Prefix to output.", true);
    argParser.parse_check(argc, argv);

    std::string reads1 = argParser.get<std::string>("reads1");
    std::string reads2 = argParser.get<std::string>("reads2");
    std::string output = argParser.get<std::string>("output");
    std::string wl = argParser.get<std::string>("whitelist");
    bool number = argParser.exist("number");
    bool library = argParser.exist("library");
    std::vector<std::string> barcodes;

    std::string line;
    std::vector<std::string> barcodeList;
    if (!number)
    {
        std::ifstream barcodeListf(wl);
        while (getline(barcodeListf, line))
            barcodeList.emplace_back(line.substr(0, line.find_first_of('\t')));
        barcodeListf.close();
    }

    std::string barcode, identifier, line1, line2;
    unsigned long long cnt_line = 0;

    // Support both compressed (.gz) and uncompressed input
    std::ifstream reads1_plain, reads2_plain;
    igzstream reads1_gz, reads2_gz;
    std::istream* reads1f_ptr = nullptr;
    std::istream* reads2f_ptr = nullptr;
    bool reads1_is_gz = (reads1.size() > 3 && reads1.substr(reads1.size() - 3) == ".gz");
    bool reads2_is_gz = (reads2.size() > 3 && reads2.substr(reads2.size() - 3) == ".gz");
    if (reads1_is_gz) {
        reads1_gz.open(reads1.c_str());
        reads1f_ptr = &reads1_gz;
    }
    else {
        reads1_plain.open(reads1.c_str());
        reads1f_ptr = &reads1_plain;
    }
    if (reads2_is_gz) {
        reads2_gz.open(reads2.c_str());
        reads2f_ptr = &reads2_gz;
    }
    else {
        reads2_plain.open(reads2.c_str());
        reads2f_ptr = &reads2_plain;
    }
    std::ofstream reads1o((output + "_1.fq").c_str());
    std::ofstream reads2o((output + "_2.fq").c_str());
    while (getline(*reads1f_ptr, line1))
    {
        getline(*reads2f_ptr, line2);
        if (cnt_line % 4 == 0)
        {
            if (cnt_line % 40000000 == 0)
                std::cout << "\rProcessed " << cnt_line / 4 << " reads ...\n" << std::flush;

            std::size_t pos1 = line1.find_first_of('#');
            std::size_t pos2 = line1.find_first_of('/', pos1 + 1);
            barcode = line1.substr(pos1 + 1, pos2 - pos1 - 1);

            std::string barcode_trans;
            std::size_t i = 0;
            std::string bc1, bc2, bc3;
            while (barcode.at(i) != '_')
                bc1 += barcode.at(i++);
            while (barcode.at(++i) != '_')
                bc2 += barcode.at(i);
            while (++i < barcode.size())
                bc3 += barcode.at(i);
            if (bc1.compare("0") && bc2.compare("0") && bc1.compare("0"))
            {
                if (!number)
                    barcode_trans = barcodeList[atoi(bc1.c_str()) - 1] + barcodeList[atoi(bc2.c_str()) - 1] + barcodeList[atoi(bc3.c_str()) - 1];
                else
                    barcode_trans = barcode;
            }
            if (barcode_trans.empty())
                identifier = line1.substr(0, pos1);
            else
            {
                if (library)
                    barcode_trans += "-1";
                identifier = line1.replace(pos1, std::string::npos, "\tBX:Z:" + barcode_trans);
            }
            reads1o << identifier << '\n';
            reads2o << identifier << '\n';
        }
        else
        {
            reads1o << line1 << '\n';
            reads2o << line2 << '\n';
        }
        cnt_line += 1;
    }
    std::cout << "\rProcessed " << cnt_line / 4 << " reads ...\n" << std::flush;
    // Only close file streams if they were opened
    if (!reads1_is_gz && reads1_plain.is_open()) reads1_plain.close();
    if (!reads2_is_gz && reads2_plain.is_open()) reads2_plain.close();
    if (reads1_is_gz) reads1_gz.close();
    if (reads2_is_gz) reads2_gz.close();
    reads1o.close();
    reads2o.close();
    return 0;
}
