 .using chain_struc = std::array<int, 1>;         // Array holding structure of a chain. Assuming 5 types of bonds in a linear chain. 

// Generate initial configurations
void generate_distribution(std::vector<chain_struc>& intact_chains, long& total_chains, double& ru_mass, long& Mn, long& Mw, double& total_mass, double& end_group_mass);

// Rate calculation and reaction selection functions
double calc_total_propensity(std::unordered_map<std::string, std::vector<double>>& reaction_rates, long& C_rad_RR, long& C_rad_BS, long& C_bonds, long& C_monomer);
std::string select_reaction(double rand2, double total_rate, std::unordered_map<std::string, std::vector<double>>& reaction_rates);

// Reaction mechanism implementations
//  - unimolecular reactions
void initiation_random_scission(std::vector<chain_struc>& intact_chains, std::vector<chain_struc>& radical_chains, std::vector<long long>& cumulative_weights, std::vector<size_t>& removed_chain_idxs, long& C_rad_BS, long& C_rad_RR, long& C_bonds, long& C_intact_chains);
void depropagation_end_beta_scission(std::vector<chain_struc>& radical_chains, long& C_monomer, long& C_rad_BS);
void termination_recombination(std::vector<chain_struc>& intact_chains, std::vector<chain_struc>& radical_chains, std::vector<long long>& cumulative_weights, std::vector<size_t>& removed_chain_idxs, long& C_rad_BS, long& C_rad_RR, long &C_bonds, long& C_intact_chains);
void side_reaction(long& C_monomer, long& C_side_product);

// Molecular weight output functions
long double calc_mol_wt(const std::vector<chain_struc>& intact_chains, std::string mol_wt_type, double repeat_unit_mass, double end_group_mass);
void save_bond_distribution(const std::vector<chain_struc>& intact_chains, double conversion, double repeat_unit_mass, double end_group_mass);