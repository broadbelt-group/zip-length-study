#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <numeric>
#include <sys/stat.h>

using chain_struc = std::array<int, 1>;
const chain_struc broken_chain = {-1};    // Represents a monomer radical, used in many functions

void save_state(const std::string& filename, std::ofstream& output_file, std::vector<chain_struc>& intact_chains, std::vector<chain_struc>& radical_chains,
    long& C_monomer, long& C_rad_RR, long& C_rad_BS, long& C_intact_chains, long& C_bonds, 
    double& time_s, double& total_mass, double& repeat_unit_mass, double& end_group_mass, 
    long double& Mn, long double& Mw) {
    // Write to output file
    double time_min = time_s/60.0;

    long double polymer_wt_frac = 0;
    for (const chain_struc& chain : intact_chains){
        if (chain != broken_chain){
            polymer_wt_frac += 2*end_group_mass + repeat_unit_mass*(std::accumulate(chain.begin(), chain.end(), 0)+1);
        }
    }
    polymer_wt_frac /= total_mass;

    double monomer_yield = (double) C_monomer * repeat_unit_mass / total_mass;
    
    // Write to output file
    output_file 
    << time_min << "," << monomer_yield << "," << polymer_wt_frac << "," << Mn << "," << Mw << "," 
    << C_monomer << "," << C_rad_RR << "," << C_rad_BS << "," << C_intact_chains  << "," << C_bonds
    << std::endl;
    output_file.flush();
}

// Helper for read_params
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \n\r\t");
    size_t last = str.find_last_not_of(" \n\r\t");
    return (first == std::string::npos || last == std::string::npos) ? "" : str.substr(first, last - first + 1);
}

// Function to read parameters from a file into a map
std::unordered_map<std::string, std::string> read_params(const std::string& filename) {
    std::unordered_map<std::string, std::string> params;
    std::ifstream infile(filename);
    std::string line;

    while (std::getline(infile, line)) {
        // Remove any comments (everything after '//')
        size_t comment_pos = line.find("//");
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // Trim leading and trailing whitespace
        line = trim(line);

        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            key = trim(key);
            value = trim(value);

            params[key] = value;
        }
    }
    return params;
}

bool file_exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}