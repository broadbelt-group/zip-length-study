#ifndef SAVE_VARIABLES_H
#define SAVE_VARIABLES_H

using chain_struc = std::array<int, 1>;

// For saving variables
void save_state(const std::string& filename, std::ofstream& output_file, std::vector<chain_struc>& intact_chains, std::vector<chain_struc>& radical_chains,
    long& C_monomer, long& C_rad_RR, long& C_rad_BS, long& C_intact_chains, long& C_bonds, long& C_side_product,
    double& time_s, double& total_mass, double& repeat_unit_mass, double& end_group_mass, long double& Mn, long double& Mw);

// For input and output
std::unordered_map<std::string, std::string> read_params(const std::string& filename);
bool file_exists(const std::string& name);

#endif // SAVE_VARIABLES_H
