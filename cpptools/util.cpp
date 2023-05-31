#include "util.h"

// stLFR or 10x (TELLSEQ)
std::string read_type = "";
std::pair<std::string, std::string> getBarcode(std::string& line)
{
    // infer reads type
    if (!read_type.compare("")){
        if (line.find("BX:Z") != std::string::npos)
            read_type = "10x";
        else if (line.find_first_of('#') != std::string::npos)
            read_type = "stLFR";
    }
    // extract barcode
    std::string read_name, barcode;
    if (!read_type.compare("stLFR")){
        std::size_t pos1 = line.find_first_of('#');
        std::size_t pos2 = line.find_first_of('/', pos1 + 1);
        read_name = line.substr(0, pos1);
        barcode = line.substr(pos1 + 1, pos2 - pos1 - 1);
        if (!barcode.compare("0_0_0"))
            barcode.clear();
    }
    else{
        read_name = line.substr(0, line.find_first_of(" \r\t\n"));
        std::size_t pos1 = line.find("BX:Z");
        if (pos1 != std::string::npos){
            std::size_t pos2 = line.find_first_of('-', pos1 + 5);
            barcode = line.substr(pos1 + 5, pos2 - pos1 - 5);
        }
    }
    return std::make_pair(read_name, barcode);
}

std::string getReadname(const std::string& line) 
{
    if (!read_type.compare("")){
        if (line.find_first_of("@") != std::string::npos)
            read_type = "fq";
        else if (line.find_first_of('>') != std::string::npos)
            read_type = "fa";
    }
    std::string read_name;
    if (!read_type.compare("fq") || !read_type.compare("fa")) {
        read_name = line.substr(0,line.find_first_of(" \r\t\n"));
    }
    return read_name;
}


