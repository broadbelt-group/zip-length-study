#include <vector>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <cmath> 
#include <algorithm> 
#include "reaction_functions.h"
#include <cassert>
#include <random>
#include <sstream>
#include <iostream>
#include <fstream> 
#include <array>
#include <numeric>

using chain_struc = std::array<int, 1>;             // Array holding structure of a chain. Assuming 1 type of bond in a linear chain. 
const chain_struc monomer_radical = {0};            // Represents a monomer radical, used in many functions
const chain_struc broken_chain = {-1};    // Represents a monomer radical, used in many functions

// Declare external global RNG (defined in main.cpp)
extern std::mt19937 global_rng;
extern std::uniform_real_distribution<double> uniform_dist;

void generate_distribution(std::vector<chain_struc>& intact_chains, long& total_chains, double& ru_mass, long& Mn, long& Mw, double& total_mass, double& end_group_mass){
    intact_chains.clear();
    intact_chains.reserve(total_chains);

    // Compute shape (k) and scale (theta) for Gamma distribution
    double k = 1.0 / ((static_cast<double>(Mw) / Mn) - 1.0);
    double theta = static_cast<double>(Mn) / ru_mass / k;

    std::gamma_distribution<double> gamma_dist(k, theta);

    // Generate chain length vector with the correct DPn
    std::vector<double> chain_lengths(total_chains);
    for (auto& L : chain_lengths) {
        L = gamma_dist(global_rng);
        if (L < 1.0) L = 1.0;  // avoid zeros
    }

    // Rescale if there are any mismatches with the actual Mn
    double Mn_actual = std::accumulate(chain_lengths.begin(), chain_lengths.end(), 0.0) / total_chains * ru_mass;
    double scale_factor = Mn / Mn_actual;
    for (auto& L : chain_lengths) L *= scale_factor;

    // Fill intact_chains
    total_mass = 0;
    for (size_t i = 0; i < static_cast<size_t>(total_chains); ++i) {
        double weight = chain_lengths[i] * ru_mass;
        int bonds_per_chain = static_cast<int>(std::round(weight / ru_mass)) - 1;
        if (bonds_per_chain < 0) bonds_per_chain = 0;

        chain_struc chain;
        chain = {bonds_per_chain};
        intact_chains.push_back(chain);

        // Tabulate total mass of system and max number of monomers/repeat units
        total_mass += 2*end_group_mass + (bonds_per_chain + 1)*ru_mass;
    }
}

// Reaction rate calculations
double calc_total_propensity(std::unordered_map<std::string, std::vector<double>>& reaction_rates, long& C_rad_RR, long& C_rad_BS, long& C_bonds, long& C_monomer) {
    // Calculate rates for unimolecular reactions
    reaction_rates["Initiation"][1] = reaction_rates["Initiation"][0] * C_bonds;
    reaction_rates["Depropagation"][1] = reaction_rates["Depropagation"][0] * C_rad_BS;

    // Rates for bimolecular reactions
    reaction_rates["Termination"][1] = reaction_rates["Termination"][0] * C_rad_RR * C_rad_RR;
    reaction_rates["SideRxn"][1] = reaction_rates["SideRxn"][0] * C_monomer * C_monomer;

    // Sum across the reaction channels
    double total_rate = 0;                                      // units of bonds/s
    for (const auto& [key, value] : reaction_rates) {
        total_rate += value[1];
    }

    return total_rate;
}

std::string select_reaction(double rand2, double total_rate, std::unordered_map<std::string, std::vector<double>>& reaction_rates) {
        double id = rand2 * total_rate;

        double tmp = 0;
        for (const auto& rxn : reaction_rates) {
            const std::string& rxn_name = rxn.first;
            const std::vector<double>& rxn_rate = rxn.second;
            
            if (id < rxn_rate[1] + tmp) {
                return rxn_name;
            } else {
                tmp += rxn_rate[1];
            }
        }

        return "ERROR";
    } 

// Bond fission helper function: weighted chain selection
int weighted_random_selection(const std::vector<chain_struc>& intact_chains, int target_index, std::vector<long long>& cumulative_weights, long long sum_weights) {
    // Generate random number
    std::uniform_int_distribution<long long> dis(0, sum_weights - 1);
    long long random_weight = dis(global_rng);
    
    // Binary search
    auto it = std::upper_bound(cumulative_weights.begin(), cumulative_weights.end(), random_weight);
    return std::distance(cumulative_weights.begin(), it);
}

// Unimolecular reactions
void initiation_random_scission(std::vector<chain_struc>& intact_chains, std::vector<chain_struc>& radical_chains, std::vector<long long>& cumulative_weights, std::vector<size_t>& removed_chain_idxs,
    long& C_rad_BS, long& C_rad_RR, long& C_bonds, long& C_intact_chains) {
    // Figure out which index needs to be targeted
    int target_index = 0;

    // Need to make sure there is a chain to break
    assert(intact_chains.size() >= 1);

    // Randomly select an intact chain to break (weighted by molecular weight)
    size_t chain_index = weighted_random_selection(intact_chains, target_index, cumulative_weights, cumulative_weights.back());
    chain_struc chain_to_break = intact_chains[chain_index];

    // Need to make sure it has bonds of that type available
    while (chain_to_break[target_index] < 1) {
        // Find a different chain
        std::uniform_int_distribution<size_t> dis(0, intact_chains.size() - 1);
        chain_index = dis(global_rng);
        chain_to_break = intact_chains[chain_index];
    }

    // If it does, save iterator
    std::vector<chain_struc>::iterator it = intact_chains.begin();
    std::advance(it, chain_index);

    /* TEMP NEW: update cumulative weights and sum to reflect broken chain */
    for (int i = chain_index; i < cumulative_weights.size(); ++i){
        cumulative_weights[i] -= chain_to_break[0];
    }
    // set the weight of the broken chain to 0
    // cumulative_weights[chain_index] = 0;
    /* END TEMP BLOCK */

    // Before removing a bond, remove all bonds from C_bonds for propensity tracking
    for (int i = 0; i < 1; ++i){
        C_bonds -= chain_to_break[i];
    }

    // Remove a bond
    chain_to_break[target_index] = chain_to_break[target_index] -  1;

    // Distribute the rest of the bonds on the chain randomly
    int tot_1 = chain_to_break[0];
    std::uniform_int_distribution<int> dis_frag(0, tot_1);
    int r_1 = dis_frag(global_rng);
    chain_struc chain_frag_1 = {r_1};
    chain_struc chain_frag_2 = {tot_1 - r_1};

    // Add the broken chains to radical container
    radical_chains.push_back(chain_frag_1);
    radical_chains.push_back(chain_frag_2);

    // "Delete" the broken chain from the main array (no matter what)
    // intact_chains.erase(it);

    /* TEMP NEW: instead of erasing chain, just set all bonds to -1? so it should never be used if chosen */
    intact_chains[chain_index] = {-1};
    removed_chain_idxs.push_back(chain_index);
    // This way intact_chains and cumulative_weights have corresponding units
    // and the deleted index is saved in another list (can be used in recombination)
    /* END TEMP BLOCK */

    C_rad_RR = radical_chains.size();
    C_rad_BS = C_rad_RR;
    C_intact_chains -= 1;

    for (const chain_struc chain : radical_chains){
        if (chain == monomer_radical) {
            C_rad_BS -= 1;
        }
    }

    assert(C_rad_BS <= C_rad_RR);
}

void depropagation_end_beta_scission(std::vector<chain_struc>& radical_chains, long& C_monomer, long& C_rad_BS) {
    // Need to make sure there is a chain to break
    assert(radical_chains.size() >= 1);

    // Select a chain, make sure it is not a monomer radical
    std::uniform_int_distribution<size_t> dis(0, radical_chains.size() - 1);
    size_t chain_index = dis(global_rng);
    chain_struc chain_to_break = radical_chains[chain_index];

    if (chain_to_break == monomer_radical){
        // Check to make sure that there are multiple possible radical chains to choose from
        assert(radical_chains.size() > 1);

        // If the chain chosen is a monomer radical, pick again
        while (chain_to_break == monomer_radical) {
            chain_index = dis(global_rng);
            chain_to_break = radical_chains[chain_index];
        }
    }

    // Remove a bond based on weights
    int total_bonds = 0;
    for (int bond_type : chain_to_break){
        total_bonds += bond_type;
    }

    // Select based on weight
    std::uniform_int_distribution<int> dis_bond(0, total_bonds - 1);
    int randomNum = dis_bond(global_rng);
    int cumulativeWeight = 0;
    int bond_pos = -1;
    for (int i = 0; i < chain_to_break.size(); ++i) {
        cumulativeWeight += chain_to_break[i];
        if (randomNum < cumulativeWeight) {
            bond_pos = i;
            break;
        }
    }

    // Alter the chain after making sure bond position is updated correctly
    assert(bond_pos >= 0);
    radical_chains[chain_index][bond_pos] -= 1;

    // Update monomer concentration for propensity calculation
    C_monomer++;
    // If a monomer radical is generated, then C_radical_BS should be decremented since it cannot participate in beta scission anymore
    if (radical_chains[chain_index] == monomer_radical) {
        C_rad_BS -= 1;
    }
}

void termination_recombination(std::vector<chain_struc>& intact_chains, std::vector<chain_struc>& radical_chains, std::vector<long long>& cumulative_weights, std::vector<size_t>& removed_chain_idxs,
    long& C_rad_BS, long& C_rad_RR, long &C_bonds, long& C_intact_chains) {
    // Need to make sure there are at least two chains to combine
    assert(radical_chains.size() >= 2);
    
    // Initialize bond counts (assuming just one type now)
    int type0 = 0;

    // Select one chain, make sure that there are bonds available, add the bonds to the counts
    std::uniform_int_distribution<size_t> dis_1(0, radical_chains.size() - 1);
    size_t chain_index_1 = dis_1(global_rng);
    chain_struc chain_to_break = radical_chains[chain_index_1];
    type0 += chain_to_break[0];
    
    // Remove first chain
    std::vector<chain_struc>::iterator it_1 = radical_chains.begin();
    std::advance(it_1, chain_index_1);
    radical_chains.erase(it_1);
    
    // Select another chain, make sure that there are bonds available
    std::uniform_int_distribution<size_t> dis_2(0, radical_chains.size() - 1);
    size_t chain_index_2 = dis_2(global_rng);
    chain_struc chain_to_break_2 = radical_chains[chain_index_2];
    type0 += chain_to_break_2[0];

    std::vector<chain_struc>::iterator it_2 = radical_chains.begin();
    std::advance(it_2, chain_index_2);
    radical_chains.erase(it_2);

    // Combine chains
    chain_struc new_chain = {type0};

    // Add extra bond formed
    new_chain[0]++;

    // Append fragments to the intact array
    // intact_chains.push_back(new_chain);
    /* TEMP NEW */
    // Put new chain into an intact_chains spot with a removed chain
    size_t new_ic_idx = removed_chain_idxs.back();
    assert(intact_chains[new_ic_idx] == broken_chain);
    removed_chain_idxs.pop_back();
    intact_chains[new_ic_idx] = new_chain;

    // Then also update the cumulative weights
    for (size_t i = new_ic_idx; i < cumulative_weights.size(); ++i){
        cumulative_weights[i] += new_chain[0];
    }
    // Also for the one that was 0, need to add the previous weight (extra?)
    // cumulative_weights[new_ic_idx] += cumulative_weights[new_ic_idx-1];

    // Add back to C_intact_chain count
    /* END TEMP NEW */

    // Update concentrations for future propensity calculations
    for (int i = 0; i < 1; ++i){
        C_bonds += new_chain[i];
    }

    C_rad_RR = radical_chains.size();
    C_rad_BS = C_rad_RR;
    C_intact_chains += 1;

    for (const chain_struc chain : radical_chains){
        if (chain == monomer_radical) {
            C_rad_BS -= 1;
        }
    }

    assert(C_rad_BS <= C_rad_RR);
}

void side_reaction(long& C_monomer, long& C_side_product) {
    C_monomer -= 2;
    C_side_product += 1;
}

// Molecular weight output functions
long double calc_mol_wt(const std::vector<chain_struc>& intact_chains, std::string mol_wt_type, double repeat_unit_mass, double end_group_mass){
    long double sum_Ni_Mi = 0;
    long double sum_Ni_Mi2 = 0;
    long double sum_Ni = intact_chains.size();

    // Loop through chains and add to sum
    for (const chain_struc& chain : intact_chains){
        if (chain != broken_chain){
            int chain_monomers = std::accumulate(chain.begin(), chain.end(), 0) +  1;
            double weight_chain = 2*end_group_mass + chain_monomers*repeat_unit_mass;
            sum_Ni_Mi2 += weight_chain*weight_chain;
            sum_Ni_Mi += weight_chain;
        }
    }

    // Return appropriate value
    if (mol_wt_type == "Mn"){
        return sum_Ni_Mi/sum_Ni;
    } else if (mol_wt_type == "Mw"){
        return sum_Ni_Mi2/sum_Ni_Mi;
    } else {
        return 0;
    }
}

void save_bond_distribution(const std::vector<chain_struc>& intact_chains, double conversion, double repeat_unit_mass, double end_group_mass){
    // Convert fractional conversion to integer percent
    int conversion_pct = (int)(conversion * 100);

    // Create output filestream with specific filename
    std::string filename = "chain_weights_" + std::to_string(conversion_pct) + "pct_conversion.csv";
    std::ofstream output_file;
    output_file.open(filename);

    // Create frequency vector and write to output file
    if (output_file.is_open()) {
        output_file << "chain_weight_g_mol\n";
        for (const chain_struc& chain : intact_chains){
            if (chain != broken_chain){
                int chain_monomers = std::accumulate(chain.begin(), chain.end(), 0) + 1;
                double weight_chain = 2*end_group_mass + chain_monomers*repeat_unit_mass;
                output_file << weight_chain << "\n";
            }
        }
        output_file.close();
    } else {
        std::cout << "Unable to write bond distribution to file" << std::endl;
    }
}