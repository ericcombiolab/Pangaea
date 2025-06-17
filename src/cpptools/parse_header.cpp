#include <iostream>
#include <fstream>
int main(int argc, char* argv[]) {
    std::ifstream fin(argv[1]);
    std::string line;
    unsigned long cnt = 0;
    while (getline(fin, line)) {
        if (line.at(0) == '>'){
            cnt += 1;
            std::cout << ">contig_" << cnt << "\n";
        }
        else
            std::cout << line << "\n";
    }
    fin.close();
    return 0;
}
