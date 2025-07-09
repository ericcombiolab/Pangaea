#include <gzstream.h>
#include <cmdline.h>
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>

struct readPair
{
    std::string header;
    std::string read1, read2;
};

int main(int argc, char* argv[])
{
    cmdline::parser argParser;
    argParser.add<std::string>("reads1", '1', "Path to reads1.", true, "");
    argParser.add<std::string>("reads2", '2', "Path to reads2.", true, "");
    argParser.add<std::string>("l1", 'l', "Path to l1.", true, "");
    argParser.add<std::string>("output", 'o', "Prefix to output.", true, "");
    argParser.parse_check(argc, argv);
    std::string reads1 = argParser.get<std::string>("reads1");
    std::string reads2 = argParser.get<std::string>("reads2");
    std::string l1 = argParser.get<std::string>("l1");
    std::string output = argParser.get<std::string>("output");

    // Support both compressed (.gz) and uncompressed input
    std::ifstream reads1_plain, reads2_plain, bcf_plain;
    igzstream reads1_gz, reads2_gz, bcf_gz;
    std::istream* reads1f_ptr = nullptr;
    std::istream* reads2f_ptr = nullptr;
    std::istream* bcf_ptr = nullptr;
    bool reads1_is_gz = (reads1.size() > 3 && reads1.substr(reads1.size() - 3) == ".gz");
    bool reads2_is_gz = (reads2.size() > 3 && reads2.substr(reads2.size() - 3) == ".gz");
    bool bcf_is_gz = (l1.size() > 3 && l1.substr(l1.size() - 3) == ".gz");
    if (reads1_is_gz) { reads1_gz.open(reads1.c_str()); reads1f_ptr = &reads1_gz; }
    else { reads1_plain.open(reads1.c_str()); reads1f_ptr = &reads1_plain; }
    if (reads2_is_gz) { reads2_gz.open(reads2.c_str()); reads2f_ptr = &reads2_gz; }
    else { reads2_plain.open(reads2.c_str()); reads2f_ptr = &reads2_plain; }
    if (bcf_is_gz) { bcf_gz.open(l1.c_str()); bcf_ptr = &bcf_gz; }
    else { bcf_plain.open(l1.c_str()); bcf_ptr = &bcf_plain; }
    unsigned long long cnt_line, out_cnt = 0;
    readPair read;
    std::string line1, line2, lineBC, barcode;
    std::unordered_map<std::string, std::vector<readPair>> barcode2reads;
    std::vector<std::string> barcodes;

    std::ofstream reads1o((output + "_1.fq").c_str()), reads2o((output + "_2.fq").c_str()), barcodewl((output + ".wl").c_str());

    std::cout << "Start scanning reads ..." << std::endl;
    while (getline(*reads1f_ptr, line1))
    {
        cnt_line += 1;
        getline(*reads2f_ptr, line2);
        getline(*bcf_ptr, lineBC);
        int flag = cnt_line % 4;
        if (flag == 1)
            read.header = line1;
        else if (flag == 2)
        {
            barcode = lineBC;
            read.read1 = "\n" + line1 + "\n+\n";
            read.read2 = "\n" + line2 + "\n+\n";
            //read.header = read.header.substr(0, read.header.find_first_of(' ')) + ":" + barcode;
            read.header = read.header.substr(0, read.header.find_first_of(' ')) + "\tBX:Z:" + barcode + "-1";
        }
        else if (flag == 0)
        {
            read.read1 = read.read1 + line1 + "\n";
            read.read2 = read.read2 + line2 + "\n";
            if (barcode.size() != 18)
            {
                std::cout << "Wrong barcode length." << std::endl;
                continue;
            }
            barcodewl << barcode << "\n";
            reads1o << read.header << read.read1;
            reads2o << read.header << read.read2;
            ++out_cnt;
            if (out_cnt % 10000000 == 0)
                std::cout << "Outputed " << out_cnt << " reads ..." << std::endl;
        }
    }
    std::cout << "Outputed " << out_cnt << " reads ...\nDone." << std::endl;
    barcodewl.close();
    reads1o.close();
    reads2o.close();
    return 0;
}