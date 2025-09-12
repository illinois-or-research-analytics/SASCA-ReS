#ifndef UTILS_H
#define UTILS_H

#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

class Utils {
    public:
        static inline char GetDelimiter(std::string filepath) {
            std::ifstream edgelist(filepath);
            std::string line;
            getline(edgelist, line);
            if (line.find(',') != std::string::npos) {
                return ',';
            } else if (line.find('\t') != std::string::npos) {
                return '\t';
            } else if (line.find(' ') != std::string::npos) {
                return ' ';
            }
            throw std::invalid_argument("Could not detect filetype for " + filepath);
        }

        static inline std::unordered_map<int, std::string> GetIndexToHeaderMap(char delimiter, std::string filepath) {
            std::unordered_map<int, std::string> index_to_header_map;
            std::ifstream input_nodelist(filepath);
            std::string line;
            std::getline(input_nodelist, line);
            std::stringstream ss(line);
            std::string current_value;
            int index = 0;
            while(std::getline(ss, current_value, delimiter)) {
                index_to_header_map[index] = current_value;
                index ++;
            }
            return index_to_header_map;
        }

        static inline std::unordered_map<std::string, int> GetHeaderToIndexMap(char delimiter, std::string filepath) {
            std::unordered_map<std::string, int> header_to_index_map;
            std::ifstream input_nodelist(filepath);
            std::string line;
            std::getline(input_nodelist, line);
            std::stringstream ss(line);
            std::string current_value;
            int index = 0;
            while(std::getline(ss, current_value, delimiter)) {
                header_to_index_map[current_value] = index;
                index ++;
            }
            return header_to_index_map;
        }
};
#endif
