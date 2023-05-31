#include <iostream>
#include <string>
#include <fstream>

std::pair<std::string, std::string> getBarcode(std::string& line);
std::string getReadname(const std::string& line);