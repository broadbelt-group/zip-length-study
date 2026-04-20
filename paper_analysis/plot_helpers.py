import os
import re 
import pandas as pd
import numpy as np

def read_data(wdir):
    data = {}
    chain_weight_files = []

    for file in os.listdir(wdir):
        if "time_based" in file:
            data['all_data_df'] = pd.read_csv(os.path.join(wdir, file))
        
        if "zip_length" in file:
            data['zip_length_df'] = pd.read_csv(os.path.join(wdir, file))

        if "chain_weights" in file:
            chain_weight_files.append(os.path.join(wdir, file))

        if "params.inp" == file:
            variables = {}
            with open(os.path.join(wdir, file), 'r') as file:
                for line in file:
                    # Remove comments
                    line = re.sub(r'//.*', '', line).strip()
                    if line:
                        match = re.match(r'(\w+)\s*=\s*([Ee0-9.\-+]+)', line)
                        if match:
                            key, value = match.groups()
                            try:
                                if '.' in value or 'e' in value.lower():
                                    variables[key] = float(value)
                                else:
                                    variables[key] = int(value)
                            except ValueError:
                                variables[key] = value 

            R = 1.987204258640E-3
            Na = 6.022E23
            data['max_time'] = variables['time']

            T = variables['T']+273
            kD = variables['A_D']*np.exp(-1*variables['Ea_D']/R/T)
            kT = variables['A_T']*np.exp(-1*variables['Ea_T']/R/T)
            try:
                kI = variables['A_I']*np.exp(-1*variables['Ea_I']/R/T)
            except:
                kI = variables['A_I']*np.exp(-1*variables['Ea_I_a']/R/T)
            try:
                kSR = variables['A_SR']*np.exp(-1*variables['Ea_SR']/R/T)
            except:
                kSR = 0

            data['m0'] = variables['repeat_unit_mass']
            data['kmc_vol'] = variables['total_chains']*variables['Mn']/variables['density']/Na*0.001
            data['num_chains'] = variables['total_chains']
            data['kD'] = kD
            data['kT'] = kT
            data['kI'] = kI
            data['kSR'] = kSR

            data['kDKMC'] = kD
            data['kTKMC'] = kT/data['kmc_vol']/Na
            data['kIKMC'] = kI
            data['kSRKMC'] = kSR/data['kmc_vol']/Na

    data['chain_weight_paths'] = chain_weight_files
           
    return data

def process_data(data):
    df = data['all_data_df'].copy()
    Na = 6.02E23

    if "# Bonds" not in df.columns:
        df['# Bonds'] = df['# Bonds a']
    
    # Calculate ratios and concentrations
    df["Rp/Rt KMC"] = (data['kDKMC']*df['# Radical Chains BS'])/(data['kTKMC']*df['# Radical Chains RR']*df['# Radical Chains RR'])
    df["Rp/Ri KMC"] = (data['kDKMC']*df['# Radical Chains BS'])/(data['kIKMC']*df['# Bonds'])
    
    df["Radical Conc BS"] = (df['# Radical Chains BS'])/data['kmc_vol']/Na
    df["Radical Conc RR"] = (df['# Radical Chains RR'])/data['kmc_vol']/Na
    df["Bond Conc"] = (df['# Bonds'])/data['kmc_vol']/Na
    
    df['d0(t) (g/L)'] = df['Mn (g/mol)']*data['num_chains']/data['kmc_vol']/Na
    df["Radical Conc [Eq 10b] Mod"] = np.sqrt(df['d0(t) (g/L)']*data['kI']/data['kT']/data['m0'])
    df["Radical Conc [Eq 10b] As-Is"] = np.sqrt(2*df['d0(t) (g/L)']*data['kI']/data['kT']/data['m0'])
    
    df["Mn/Mn0 (DP)"] = df['Mn (g/mol)']/data['m0']
    
    df["Corrected ZL"] = data['kD']*np.sqrt(data['m0']/(data['kT']*data['kI']*df['d0(t) (g/L)']))
            
    kD = data['kD']
    kT = data['kT']
    kI = data['kI']

    df['2a'] = df['Rp/Rt KMC']
    df['2b'] = df['Rp/Ri KMC']
    df['5a BS'] = (kD*df['Radical Conc BS'])/(kT*df['Radical Conc RR']*df['Radical Conc RR'])
    df['5a BS'] = (kD*df['Radical Conc BS'])/(kT*df['Radical Conc RR']*df['Radical Conc RR'])
    df['5b BS'] = kD*df['Radical Conc BS'] / (kI*df['Bond Conc'])
    df['5a RR'] = (kD*df['Radical Conc RR'])/(kT*df['Radical Conc RR']*df['Radical Conc RR'])
    df['5b RR'] = kD*df['Radical Conc RR'] / (kI*df['Bond Conc'])
    
    return df

def process_data_eci(data):
    df = data['all_data_df'].copy()
    Na = 6.02E23
    
    # Calculate ratios and concentrations
    df["Rp/Rt KMC"] = (data['kDKMC']*df['# Radical Chains BS'])/(data['kTKMC']*df['# Radical Chains RR']*df['# Radical Chains RR'])
    df["Mn/Mn0 (DP)"] = df['Mn (g/mol)']/data['m0']
    
    df["Radical Conc BS"] = (df['# Radical Chains BS'])/data['kmc_vol']/Na
    df["Radical Conc RR"] = (df['# Radical Chains RR'])/data['kmc_vol']/Na
    
    df['d0(t) (g/L)'] = df['Mn (g/mol)']*data['num_chains']/data['kmc_vol']/Na
    df["Radical Conc [Eq 10b] Mod"] = np.sqrt(df['d0(t) (g/L)'].iloc[0]*data['kI']/data['kT']/data['m0']/df["Mn/Mn0 (DP)"])
    df["Radical Conc [Eq 10b] Mod d(t)"] = np.sqrt(df['d0(t) (g/L)']*data['kI']/data['kT']/data['m0']/df["Mn/Mn0 (DP)"])

    
    # for ECI, use varying x but constant density
    df["Corrected ZL"] = data['kD']*np.sqrt(data['m0']*df["Mn/Mn0 (DP)"]/(data['kT']*data['kI']*df['d0(t) (g/L)'].iloc[0]))
            
    kD = data['kD']
    kT = data['kT']
    kI = data['kI']

    df['2a'] = df['Rp/Rt KMC']
    
    return df

def process_data_mom(data):
    df = data['all_data_df'].copy()
    Na = 6.02E23
    
    # Calculate ratios and concentrations
    df["Rp/Rt"] = df["netrate_BS"] / df["netrate_R"]
    df["Rp/Ri"] = df["netrate_BS"] / df["netrate_BF"]
    
    df["Radical Conc"] = df['Rhat']
    
    df['d0(t) (g/L)'] = df['d0']
    df["Radical Conc [Eq 10b]"] = np.sqrt(df['d0(t) (g/L)']*data['kI']/data['kT']/data['m0'])
    
    df["Mn/Mn0 (DP)"] = df['Mn_poly_1']/data['m0']
    
    df['Corrected ZL'] = data['kD']*np.sqrt(data['m0']/(data['kT']*data['kI']*df['d0(t) (g/L)']))
    df["MoM ZL"] = df['ZL']

    df['Time (min)'] = df['time']/60
    
    return df

def read_data_mom(wdir):
    """Read all necessary data for a given run"""    
    # Read simulation data
    all_data_df = read_graph_data(wdir)
    all_data_df = all_data_df.apply(pd.to_numeric, errors='coerce')
    all_data_df = all_data_df.iloc[1:]
    
    # Read and process parameters
    param_df = read_param_data(wdir)
    
    # Read misc params from ddat.in
    with open(os.path.join(wdir, "ddat.in"), 'r') as f:
        for line_num, line in enumerate(f, 1):  # Start enumeration from 1 for line numbers
            if line_num == 25:
                Mn0 = float(line.strip())
            if line_num == 26:
                Mw0 = float(line.strip())
            if line_num == 46:
                temp_kelvin = float(line.strip())
                break
    
    # Calculate rate constants
    rate_constants = calculate_rate_constants(param_df, temp_kelvin)
    
    # Combine all data into a dictionary
    data = {
        "all_data_df": all_data_df,
        "kI": rate_constants['kI'],
        "kD": rate_constants['kD'],
        "kT": rate_constants['kT'],
        "m0": all_data_df['Monomer_MW_1'].mode()[0],
        "T": temp_kelvin,
        "Mn0": Mn0,
        "Mw0": Mw0,
    }
    
    return data

def read_graph_data(wdir):
    """Read and process the graph.out file containing simulation data"""
    df = pd.read_csv(os.path.join(wdir, "graph.out"), sep='\t', header=0, index_col=False, engine="python")
    # Clean up column names by removing everything after ':'
    df.columns = df.columns.str.replace(r":.*", "", regex=True)
    return df

def calculate_rate_constants(param_df, T):
    """Calculate reaction rate constants from parameters"""
    R = 1.987204258640E-3
    rate_constants = {}
    
    # Mapping of reaction names to rate constant keys
    reaction_mapping = {
        'chainf_epep': 'kI',
        'recomb_epep': 'kT',
        'end_beta_ppp': 'kD',
        'siderxn_D_1s2pp': 'kSR'
    }
    
    for _, row in param_df.iterrows():
        if row["Name"] in reaction_mapping:
            k = row['A'] * np.exp(-1 * row['Ea'] / R / T)
            rate_constants[reaction_mapping[row["Name"]]] = k
            
    return rate_constants

def read_param_data(wdir):
    """Read and process the param.inp file containing reaction parameters"""
    param_df = pd.read_csv(os.path.join(wdir, "param.inp"), 
                          sep=":", header=None, 
                          index_col=0, 
                          names=["A", "Unk1", "Unk2", "Ea", "Irrelevant", "Name"])

    # Clean up whitespace and missing name entries
    param_df["Name"] = param_df["Name"].replace("--", None).str.strip()

    # Convert numeric columns from strings to floats
    numeric_cols = ["A", "Ea"]
    param_df[numeric_cols] = param_df[numeric_cols].apply(pd.to_numeric, errors="coerce")

    return param_df.dropna()