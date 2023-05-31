#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include "cmdline.h"
#include <dirent.h>
#include "gzstream.h"


std::string getReadname(const std::string& line) 
{
    std::string read_name = "";
    std::string read_type = "";

    // infer reads type
    if (!read_type.compare("")){
        if (line.find("BX:Z") != std::string::npos)
            read_type = "10x";
        else if (line.find_first_of('#') != std::string::npos)
            read_type = "stLFR";
    }
    if (read_type.compare("stLFR") == 0){
        std::size_t pos1 = line.find_first_of('@');
        std::size_t pos2 = line.find_first_of('/', pos1 + 1);
        read_name = line.substr(pos1+1, pos2-1);
    }
    else{
        if (line.find_first_of("@") != std::string::npos){
            read_name = line.substr(1,line.find_first_of(" \r\t\n") - 1);
        }
    }
    return read_name;
}

int main(int argc, char* argv[])
{
	cmdline::parser argParser;
	argParser.add<std::string>("map", 'm', "path to short read-barcode mapping file.", true);
    argParser.add<std::string>("fastq1", '1', "path to fastq1 file.", true);
    argParser.add<std::string>("fastq2", '2', "path to fastq2 file.", true);
	argParser.add<std::string>("output", 'o', "outputfile", true);
	argParser.parse_check(argc, argv);
    std::cout << "Start loading long read to barcode mapping." << std::endl;
	std::string output_path = argParser.get<std::string>("output");
    std::string mapping = argParser.get<std::string>("map");
    std::string fastq1_path = argParser.get<std::string>("fastq1");
    std::string fastq2_path = argParser.get<std::string>("fastq2");

	std::unordered_map<std::string, std::string> read_to_barcode;
	std::ifstream mapping_f(mapping.c_str());
	std::string line;
	while (getline(mapping_f, line))
	{
		std::string name, barcode;
		size_t pos = 0;
		pos = line.find("BX:Z:");
        if (pos == std::string::npos)
        {
            barcode = " ";
        } else {
		    barcode = line.substr(pos);
        }
		name = line.substr(0, pos);
        //remove the begining and end space
        name.erase(0, name.find_first_not_of(" \t\r\n"));
        name.erase(name.find_last_not_of(" \t\r\n") + 1);
		if (read_to_barcode.find(name) != read_to_barcode.end())
			std::cout << "Barcode for short read " << name << " is already in the list." << std::endl;
		else
			read_to_barcode[name] = barcode;
	}
	std::cout << "Print the first 10 samples for checking." << std::endl;
	int cnt = 0;
	for (auto&& i : read_to_barcode) {
		std::cout << i.first << " --> " << i.second << std::endl;
		if (++cnt == 10)
			break;
	}
    std::cout<<read_to_barcode.size()<<std::endl;
    int map_size = read_to_barcode.size();

    // start to read the fastq file use openmp to parallel
    std::cout << "Start reading fastq file." << std::endl;

    igzstream fastq1(fastq1_path.c_str());
    igzstream fastq2(fastq2_path.c_str());
    std::string line1, line2;
    std::ofstream output(output_path.c_str());

    int read_cnt = 0;
    int openmp_size = 10000000;
	int openmp_index = 0;
    bool stop = false;
    std::vector<std::string> reads1(openmp_size);
	std::vector<std::string> reads2(openmp_size);
    std::vector<std::string> interleaved_result(openmp_size);
    std::vector<std::string> barcode(openmp_size);
    std::vector<std::string> read_name(openmp_size);
    std::vector<std::vector<std::pair<std::string, std::string>>> results(openmp_size);
    // for debug
    int cnt_name = 0;
    std::cout<<" print the first 10 names for checking:"<<std::endl;
    do
    {
        if (getline(fastq1, line1) && getline(fastq2, line2)){
            read_cnt++;
            read_name[openmp_index] = std::move(getReadname(line1));
            if ( cnt_name++ < 10) {
                std::cout<<cnt_name<<" name:"<<read_name[openmp_index]<<"->"<<read_to_barcode[read_name[openmp_index]]<<std::endl;
            }
            // break;
            getline(fastq1, line1);
            getline(fastq2, line2);
            // read three lines of fastq1 at a time
            reads1[openmp_index] = std::move(line1);
            reads2[openmp_index] = std::move(line2);
            for (int i = 0; i < 2; ++i)
            {
                reads1[openmp_index] += '\n';
                getline(fastq1, line1);
            
                reads1[openmp_index] += (i==0)?"+":line1;

                reads2[openmp_index] += '\n';
                getline(fastq2, line2);
                reads2[openmp_index] += (i==0)?"+":line2;
            }
            // std::cout<<reads1[openmp_index]<<std::endl;
            openmp_index++;
        } else {
            openmp_size = openmp_index;
            stop = true;
            std::cout<<"total reads :<< "<<read_cnt<<std::endl;
        }
        if (openmp_index >= openmp_size) {
            std::cout<<"openmp_index:"<<openmp_index<<std::endl;
            #pragma omp parallel for
            for (int i = 0; i < openmp_index; ++i)
            {
                if (read_to_barcode.find(read_name[i]) != read_to_barcode.end() && read_to_barcode[read_name[i]] != " ") {
                    barcode[i] = read_to_barcode[read_name[i]];
                    interleaved_result[i] = "@" + read_name[i] + " " + barcode[i] + "-1\n" + reads1[i] + "\n" + "@" + read_name[i] + " " + barcode[i] + "-1\n" + reads2[i] + "\n";
                }
                else {
                    barcode[i] = "";
                    interleaved_result[i] = "@" + read_name[i] + "\n" + reads1[i] + "\n" + "@" + read_name[i] + "\n" + reads2[i] + "\n";

                }
            }
            for (int i = 0; i < openmp_index; ++i)
            {
                // std::cout<< interleaved_result[i];
                output << interleaved_result[i];
            }
            openmp_index = 0;
        }
        if (read_cnt % 10000000 == 0){ 
            std::cout << read_cnt << " reads have been processed." << std::endl;
            // break;
        }
    } while (!stop);
    std::cout << "Total " << read_cnt << " reads have been processed." << std::endl;
    output.close();
    return 0;
}