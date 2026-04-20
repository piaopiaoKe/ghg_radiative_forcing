#%% impoart package and make evertything in local scale
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Aearth=5.1e14 # m2
RE_co2_ppb, RE_n2o_ppb, RE_ch4_ppb=1.33e-5, 38.8e-5, 32e-5 # IPCC AR6 file:///Z:/NO%20flux.Data/IPCC_AR6_WGI_Chapter07_SM.pdf
# Matm=5.135e18 kg, mair_avearge=28.965e-3 kg/mol
co2_kg_ppb, n2o_kg_ppb, ch4_kg_ppb=7.80e9, 2.84e9, 7.80e9 #kg/ppb
RE_co2_kg, RE_n2o_kg, RE_ch4_kg=RE_co2_ppb/co2_kg_ppb, RE_n2o_ppb/n2o_kg_ppb, RE_ch4_ppb/ch4_kg_ppb # W/m2/kg
# scale to local by distriuting total eneryg to local
RE_co2_site, RE_n2o_site, RE_ch4_site=RE_co2_kg*Aearth, RE_n2o_kg*Aearth, RE_ch4_kg*Aearth # W/kg 

def irf_co2(years):
    # IRF parameters: using AR6-like multi-exponential (example from dynamic_characterization)
    alpha0, alpha1, alpha2, alpha3 = 0.2173, 0.2240, 0.2824, 0.2763
    tau1, tau2, tau3 = 394.4, 36.54, 4.304  # years
    """Impulse response function for CO2 fraction remaining in atmosphere."""
    # years: numpy array of time since perturbation in years
    term0 = alpha0
    exp1 = alpha1 * np.exp(-years / tau1)
    exp2 = alpha2 *  np.exp(-years / tau2)
    exp3 = alpha3 * np.exp(-years / tau3)
    return term0 + exp1 + exp2 + exp3


def irf_n2o(years):    
    tau1=109  # years
    """Impulse response function for N2O fraction remaining in atmosphere."""
    # years: numpy array of time since perturbation in years
    term = np.exp(-years/tau1)
    return term

def irf_ch4(years):    
    tau1=11.8  # years
    """Impulse response function for CH4 fraction remaining in atmosphere."""
    # years: numpy array of time since perturbation in years
    term = np.exp(-years/tau1)
    return term


# CO2 RF
" using grassland as reference, compare the radiative forcing in "
"growing forest at this place in area of A"
age = np.array([0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80])
nep_forest = -np.array([0, 0, 20, 150, 300, 420, 450, 350, 220, 140, 80, 40]) # net ecocystem productivity, gC/m2/year
nep_grass = -np.array([66] * len(nep_forest))

age_new = np.arange(0, 81, 1)
nep_forest_new = np.interp(age_new, age, nep_forest)
nep_grass_new = np.interp(age_new, age, nep_grass)  # stays 66
dnep=nep_forest_new-nep_grass_new
cstock = pd.DataFrame({
    "age": age_new,
    "nep_forest": nep_forest_new,
    "nep_grass": nep_grass_new,
    "dnep": dnep
})

# A=4e7*100 # total area of grassland
# radiative_efficiency = 1.75e-15  # W m^-2 per kg CO2, example from literature
# Constants / parameters
M_CO2 = 44.01   # molecular weight CO2 (g/mol)
M_C = 12.01     # molecular weight Carbon (g/mol)
kg_per_g = 1e-3
max_years = 80
time = age_new 
irf_vals = irf_co2(time)  # fraction (or effective residual) per time since pulse

forest_gC = cstock["nep_forest"].values
forest_kgCO2 = forest_gC * (M_CO2 / M_C)* kg_per_g #kg co2/m2 (forest)
rf_forest = np.zeros(len(age_new))
for t in range(len(age_new)):
    for i in range(t + 1):
        rf_forest[t] += forest_kgCO2[i] * irf_vals[t - i]*RE_co2_site
        
grass_gC = cstock["nep_grass"].values
grass_kgCO2 = grass_gC * (M_CO2 / M_C)* kg_per_g
rf_grass = np.zeros(len(age_new))
for t in range(len(age_new)):
    for i in range(t + 1):
        rf_grass[t] += grass_kgCO2[i] * irf_vals[t - i]*RE_co2_site
rf_net = rf_forest - rf_grass    

cstock["rf_forest"] = rf_forest
cstock["rf_grass"] = rf_grass
cstock["rf_net"] = rf_net

data=cstock.copy()
plt.figure()
plt.plot(data['age'],data['rf_net'])
plt.ylabel("Radiative forcing (W/m²)")
plt.title("Estimated CO₂ Radiative Forcing from NEP-based fluxes")
plt.show()

#% N2O emission RF
age = np.arange(0,81,1)
n2o_forest = np.array([9.4/1e3] * len(age))/1e3 # kgN2O/m2/year
n2o_grass = np.array([34/1e3] * len(age))/1e3 # kg/m2/year
# assuem pulse emission

nstock = pd.DataFrame({
    "age": age_new,
    "n2o_forest": n2o_forest,
    "n2o_grass": n2o_grass
})
time = age
irf_vals = irf_n2o(time)  # fraction (or effective residual) per time since pulse
forest_kgN2O = nstock['n2o_forest'] #kg co2/m2 (forest)
rf_forest = np.zeros(len(age))
for t in range(len(age)):
    for i in range(t + 1):
        rf_forest[t] += forest_kgN2O[i] * irf_vals[t - i]*RE_n2o_site
  
grass_kgN2O = nstock['n2o_grass'] #kg co2/m2 (forest)
rf_grass = np.zeros(len(age))
for t in range(len(age)):
    for i in range(t + 1):
        rf_grass[t] += grass_kgN2O[i] * irf_vals[t - i]*RE_n2o_site

rf_net = rf_forest - rf_grass    
nstock["rf_forest"] = rf_forest
nstock["rf_grass"] = rf_grass
nstock["rf_net"] = rf_net

data=nstock.copy()
plt.figure()
plt.plot(data['age'],data['rf_net'])
plt.ylabel("Radiative forcing (W/m²)")
plt.title("Estimated N2O Radiative Forcing from annnual fluxes")
plt.show()

#% CH4 emission RF
age = np.arange(0,81)
ch4_forest = np.array([-127.7/1e3] * len(age))/1e3 # kgN2O/m2/year,negatiev mean uptake
ch4_grass = np.array([-584.0/1e3] * len(age))/1e3 # kg/m2/year
# assuem pulse emission

ch4stock = pd.DataFrame({
    "age": age_new,
    "ch4_forest": ch4_forest,
    "ch4_grass": ch4_grass
})
time = age
irf_vals = irf_ch4(time)  # fraction (or effective residual) per time since pulse
forest_kgch4 = ch4stock['ch4_forest'] #kg co2/m2 (forest)
rf_forest = np.zeros(len(age))
for t in range(len(age)):
    for i in range(t + 1):
        rf_forest[t] += forest_kgch4[i] * irf_vals[t - i]*RE_ch4_site
  
grass_kgch4 = ch4stock['ch4_grass'] #kg co2/m2 (forest)
rf_grass = np.zeros(len(age))
for t in range(len(age)):
    for i in range(t + 1):
        rf_grass[t] += grass_kgch4[i] * irf_vals[t - i]*RE_ch4_site

rf_net = rf_forest - rf_grass    
ch4stock["rf_forest"] = rf_forest
ch4stock["rf_grass"] = rf_grass
ch4stock["rf_net"] = rf_net

data=ch4stock.copy()
plt.figure()
plt.plot(data['age'],data['rf_net'])
plt.ylabel("Radiative forcing (W/m²)")
plt.title("Estimated CH4 Radiative Forcing from annual fluxes")
plt.show()

idx = cstock['age'].isin([1, 25, 50, 80])
print(f"CO2-induced RF is {cstock.loc[idx, ['age', 'rf_net']]}")
print(f"N2O-induced RF is {nstock.loc[idx, ['age', 'rf_net']]}")
print(f"CH4-induced RF is {ch4stock.loc[idx, ['age', 'rf_net']]}")
