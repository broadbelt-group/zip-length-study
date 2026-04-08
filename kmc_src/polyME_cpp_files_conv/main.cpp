#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <array>
#include <limits>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream> 
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <chrono>
#include <unistd.h>
#include <random>
#include "reaction_functions.h"
#include "save_variables.h"

// Zip length study code
// author: Shivani Kozarekar

using chain_struc = std::array<int, 1>;             // Array holding structure of a chain
const chain_struc monomer_radical = {0};            // Represents a monomer radical, used in many functions

// Global RNG - must be initialized before main loop
std::mt19937 global_rng;
std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

// Main function
int main() {
    // Keep time
    auto start_time = std::chrono::high_resolution_clock::now();

    // Seed the random number generator robustly
    // Use multiple entropy sources to avoid collision on HPC systems
    std::random_device rd;
    std::seed_seq seed{rd(), static_cast<unsigned int>(std::time(nullptr)), 
                      static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count())};
    global_rng.seed(seed);
    
    // Constants
    const double Na = 6.022E23;             // Avogadro's number (chains per mol)
    const double R = 1.987204258640E-3;     // ideal gas constant (kcal/mol/K)
    const int BOND_TYPES = 5;               // Number of bond types

    // Read parameters from params.inp into map
    std::unordered_map<std::string, std::string> params = read_params("params.inp");
    
    // Begin reading variables out of map
    int max_simulation_time = std::stoi(params["time"]) * 60;                   // max simulation time (seconds)
    int T = std::stoi(params["T"]) + 273;                                       // temperature (K)
    long total_chains = static_cast<long>(std::stod(params["total_chains"]));   // number of initial chains to simulate
    long Mn_0 = std::stol(params["Mn"]);                                        // number-average molecular weight (g/mol)
    long Mw_0 = std::stol(params["Mw"]);                                        // weight-average molecular weight (g/mol)
    double repeat_unit_mass = std::stod(params["repeat_unit_mass"]);            // molar mass of repeat unit (g/mol)
    double end_group_mass = std::stod(params["end_group_mass"]);                // molar mass of one end group (g/mol)
    double density = std::stod(params["density"]);                              // density of the polymer (g/mL or g/cm3)
    double V_chain = Mn_0/density/Na*0.001;                                     // reaction volume (L/chain)                           
    double V = total_chains*V_chain;                                            // scaled reaction volume (total in L)

    // Kinetic parameters
    double A_I = std::stod(params["A_I"]);
    double Ea_I = std::stod(params["Ea_I"]);

    double A_D = std::stod(params["A_D"]);
    double Ea_D = std::stod(params["Ea_D"]);

    double A_T = std::stod(params["A_T"]);
    double Ea_T = std::stod(params["Ea_T"]);

    double A_SR = std::stod(params["A_SR"]);
    double Ea_SR = std::stod(params["Ea_SR"]);

    // Reaction rate constant map
    // first element: name of reaction | second element: tuple of {rate constant k_KMC from DFT/lit, total rate initialized to 0}
    std::unordered_map<std::string, std::vector<double>> reaction_rates_and_params = {
        // First order (k_micro same as k_macro)
        {"Initiation", {A_I*exp(-1*Ea_I/R/T), 0}},
        {"Depropagation", {A_D*exp(-1*Ea_D/R/T), 0}},

        // Second order (need to scale)
        {"Termination", {A_T*exp(-1*Ea_T/R/T)/V/Na, 0}}, // No 2 because reactants are distinguishable
        {"SideRxn", {2*A_SR*exp(-1*Ea_SR/R/T)/V/Na, 0}}, // Reactants indistingushable
    };

    // Output filestreams
    std::string filename = "time_based_data.csv";
    std::ofstream output_file;
    output_file.open(filename, std::ios::app);

    std::string zip_length_filename = "zip_length.csv";
    std::ofstream zip_length_file;
    zip_length_file.open(zip_length_filename, std::ios::app);

    // Variables to load reaction time data into
    std::vector<chain_struc> intact_chains;
    std::vector<chain_struc> radical_chains;
    double monomer_yield;
    double time_s = 0;
    double total_mass = 0;

    // Concentration variables (in units of # of ____)
    long C_monomer = 0;
    long C_rad_RR = 0;
    long C_rad_BS = 0;
    long C_intact_chains = 0;
    long C_bonds = 0;
    long C_side_product = 0;

    /* SPEEDUP */
    std::vector<long long> cumulative_weights;
    std::vector<size_t> removed_chain_idxs;
    long long sum_weights = 0;
    /* END BLOCK */

    // Start from scratch
    output_file << "Time (min),Monomer Conversion,Polymer Wt Fraction,Side Product Wt Fraction,Mn (g/mol),Mw (g/mol),"
    << "# Monomers,# Radical Chains RR,# Radical Chains BS,# Intact Chains,# Bonds,# Side Prods"
    << std::endl;
    output_file.flush();

    zip_length_file << "Zip Length (#),Zip End Time (s)" << std::endl;
    zip_length_file.flush();

    // Create intact chain initial config from Mns and Mws (log normal distribution), also updates total_mass
    generate_distribution(intact_chains, total_chains, repeat_unit_mass, Mn_0, Mw_0, total_mass, end_group_mass);

    // Write initial bond distribution (0% conversion)
    save_bond_distribution(intact_chains, 0.0, repeat_unit_mass, end_group_mass);

    for (const chain_struc& chain : intact_chains){
        for (size_t i = 0; i < 1; ++i){
                C_bonds += chain[i];
        }
        C_intact_chains = intact_chains.size();
    }
    std::cout << "Total Mass: " << total_mass << std::endl;
    std::cout << "Boyd ZL0: " << (reaction_rates_and_params["Depropagation"][0]*std::sqrt(repeat_unit_mass/density/1000/reaction_rates_and_params["Initiation"][0]/(reaction_rates_and_params["Termination"][0]*V*Na))) << std::endl;
    std::cout << "kI: " << reaction_rates_and_params["Initiation"][0] << std::endl;
    std::cout << "kD: " << reaction_rates_and_params["Depropagation"][0] << std::endl;
    std::cout << "kT: " << reaction_rates_and_params["Termination"][0]*V*Na << std::endl;

    /* SPEEDUP */
    cumulative_weights.reserve(intact_chains.size());
    size_t target_index = 0;
    for (const chain_struc& chain : intact_chains) {
        sum_weights += chain[target_index];
        cumulative_weights.push_back(sum_weights);
    }
    /* END BLOCK */ 

    // Main KMC loop
    double conversion_print_threshold = 0.10;  // Save MWD at every 0.10 (10%) monomer conversion
    int csv_print_time = 0;                 // seconds
    std::string rxn_name;
    double tau, rand1, rand2;    

    bool in_block = false;
    int zip_count = 0;

    while (time_s < max_simulation_time) {
        double total_propensity = calc_total_propensity(reaction_rates_and_params, C_rad_RR, C_rad_BS, C_bonds, C_monomer);

        if (total_propensity == 0){break;}

        // Full KMC iteration
        rand1 = uniform_dist(global_rng);
        rand2 = uniform_dist(global_rng);
        
        // Ensure rand1 is not exactly 0 or 1 to prevent log(0) or log(inf)
        rand1 = std::max(1e-16, std::min(1.0 - 1e-16, rand1));

        // Calculate tau using rand1
        // Added protection against numerical issues
        tau = (1.0 / total_propensity) * std::log(1.0 / rand1);
        
        // Safety check: if tau becomes inf or nan, something went wrong
        if (!std::isfinite(tau)) {
            std::cerr << "ERROR: tau is not finite! total_propensity = " << total_propensity << ", rand1 = " << rand1 << std::endl;
            std::cerr << "C_rad_RR = " << C_rad_RR << ", C_rad_BS = " << C_rad_BS << ", C_bonds = " << C_bonds << ", C_monomer = " << C_monomer << std::endl;
            break;
        }
        
        // Determine which reaction occurs using rand2
        rxn_name = select_reaction(rand2, total_propensity, reaction_rates_and_params);
        
        // Execute reaction event
        if (rxn_name.find("Initiation") != std::string::npos) {
            initiation_random_scission(intact_chains, radical_chains, cumulative_weights, removed_chain_idxs, C_rad_BS, C_rad_RR, C_bonds, C_intact_chains);

            if (!in_block) {
                in_block = true;
            }
        } else if (rxn_name == "Depropagation"){
            depropagation_end_beta_scission(radical_chains, C_monomer, C_rad_BS);

            zip_count += 1;
        } else if (rxn_name == "Termination") {
            termination_recombination(intact_chains, radical_chains, cumulative_weights, removed_chain_idxs, C_rad_BS, C_rad_RR, C_bonds, C_intact_chains);

            zip_length_file << zip_count << "," << time_s+tau << std::endl;
            zip_length_file.flush();
            zip_count = 0;
            in_block = false;
        } else if (rxn_name == "SideRxn"){
            side_reaction(C_monomer, C_side_product);
        }

        // Output to file
        if (time_s > csv_print_time) {
            long double Mn = calc_mol_wt(intact_chains, "Mn", repeat_unit_mass, end_group_mass);
            long double Mw = calc_mol_wt(intact_chains, "Mw", repeat_unit_mass, end_group_mass);

            double monomer_conversion = static_cast<double>(C_monomer) * repeat_unit_mass / total_mass;
            double polymer_wt_fraction = 0;
            for (const auto& chain : intact_chains) {
                if (chain[0] > 0) {
                    polymer_wt_fraction += (chain[0] + 1) * repeat_unit_mass + 2 * end_group_mass;
                }
            }
            polymer_wt_fraction /= total_mass;
            double sp_wt_frac = static_cast<double>(C_side_product) * 396.18 / total_mass;

            output_file << std::fixed << std::setprecision(6) << time_s/60.0 << ","
                        << monomer_conversion << "," << polymer_wt_fraction << "," << sp_wt_frac << ","
                        << Mn << "," << Mw << ","
                        << C_monomer << "," << C_rad_RR << "," << C_rad_BS << ","
                        << C_intact_chains << "," << C_bonds << "," << C_side_product << std::endl;
            output_file.flush();

            csv_print_time += 10;                   // print every 10 seconds
        }

        // Calculate current monomer conversion for MWD saving
        double current_conversion = static_cast<double>(C_monomer) * repeat_unit_mass / total_mass;

        // If monomer conversion exceeds the next threshold, save weight of each chain
        if (current_conversion >= conversion_print_threshold) {
            save_bond_distribution(intact_chains, current_conversion, repeat_unit_mass, end_group_mass);
            conversion_print_threshold += 0.10; // increment by 0.10 (10%) conversion
        }

        // Increment time
        time_s += tau;
    }

    // Save final state
    long double Mn_final = calc_mol_wt(intact_chains, "Mn", repeat_unit_mass, end_group_mass);
    long double Mw_final = calc_mol_wt(intact_chains, "Mw", repeat_unit_mass, end_group_mass);

    double monomer_conversion_final = static_cast<double>(C_monomer) * repeat_unit_mass / total_mass;
    double polymer_wt_fraction_final = 0;
    for (const auto& chain : intact_chains) {
        if (chain[0] > 0) {
            polymer_wt_fraction_final += (chain[0] + 1) * repeat_unit_mass + 2 * end_group_mass;
        }
    }
    polymer_wt_fraction_final /= total_mass;
    double sp_wt_frac_final = static_cast<double>(C_side_product) * 396.18 / total_mass;

    output_file << std::fixed << std::setprecision(6) << time_s/60.0 << ","
                << monomer_conversion_final << "," << polymer_wt_fraction_final << "," << sp_wt_frac_final << ","
                << Mn_final << "," << Mw_final << ","
                << C_monomer << "," << C_rad_RR << "," << C_rad_BS << ","
                << C_intact_chains << "," << C_bonds << "," << C_side_product << std::endl;
    output_file.flush();

    // Close output files
    output_file.close();
    zip_length_file.close();

    // Keep time
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = end_time - start_time;
    std::cout << "Simulation time: " << std::round(ms_double.count() / 1000 / 60 * 10.0) / 10.0 << " min" << std::endl;
    
    return 0;
}
