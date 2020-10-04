from math import pi
import pandas as pnd
import numpy as np
import COPV_tool_functions as functions

# INPUTS
# Read data from Excel:
input_data_all = pnd.read_excel(io='Tank and composite data.xlsx', sheet_name="Sheet1", index_col=0, header=None)
# Propellant and tank related inputs
prop_names = input_data_all.loc["Oxidizer/fuel (name)", 1]
mixture_fractions = input_data_all.loc["Mixture fractions (-)", 1]
m_prop = input_data_all.loc["Fuel/oxidizer mass (kg)", 1]
p_prop = input_data_all.loc["Fuel/oxidizer pressure (bar)", 1]
T_prop = input_data_all.loc["Fuel/oxidizer temperature (K)", 1]
v_ullage = input_data_all.loc["Ullage volume (-)", 1]
proof_factor = input_data_all.loc["Proof factor (-)", 1]
safety_factor = input_data_all.loc["Safety factor (-)", 1]
LD_ratio = input_data_all.loc["L/D (-)", 1]
tank_r = input_data_all.loc["Radius (m)", 1]
# Composite related inputs
weave = input_data_all.loc["Weave pattern (degrees)", 1]
layers_no = input_data_all.loc["Numbers of layers (-)", 1]
t_ply = input_data_all.loc["Ply thickness (mm)", 1]
sigma_fibres = input_data_all.loc["Fibres tensile strength (MPa)", 1]
E_fibres = input_data_all.loc["Fibres E-modulus (GPa)", 1]
rho_fibre = input_data_all.loc["Fibres density (kg/m^3)", 1]
rho_matrix = input_data_all.loc["Matrix density (kg/m^3)", 1]
translational_eff = input_data_all.loc["Translational efficiency (-)", 1]
fibre_v = input_data_all.loc["Fibres volume fraction (-)", 1]
# Liner inputs
t_lin = input_data_all.loc["Liner thickness (mm)", 1]
rho_lin = input_data_all.loc["Liner density (kg/m^3)", 1]
sigma_lin = input_data_all.loc["Liner yield strength (MPa)", 1]
E_lin = input_data_all.loc["Liner E modulus (GPa)", 1]
poisson_ratio = input_data_all.loc["Liner Poisson ratio (-)", 1]

# CHANGE UNITS AND TRANSFORM STRINGS
weave_matrix = np.array([float(teta) for teta in weave.split(',')]) * pi / 180
layers_matrix = np.array([float(ratio) for ratio in layers_no.split(',')])
p_prop *= 10**5
sigma_fibres *= 10**6
E_fibres *= 10**9
t_ply /= 10**3
t_lin /= 10**3
E_lin *= 10**9

# CALCULATE COMPOSITE PROPERTIES
rho_comp = (1-fibre_v)*rho_matrix + fibre_v*rho_fibre
E_comp = fibre_v*E_fibres
sigma_comp = translational_eff*fibre_v*sigma_fibres

# FUEL/OXIDIZER PROPERTIES
rho_prop, tank_vol = functions.properties_calc(prop_names, mixture_fractions, T_prop, m_prop, v_ullage)
print("Propellant volume is {tank_vol} m^3 \nPropellant density is {rho_prop} kg/m^3".format(tank_vol=tank_vol,
                                                                                             rho_prop=rho_prop))

# TANK GEOMETRY
tank_r, tank_l, S_cylinder, S_endcaps, LD_ratio = functions.geometry_gen(tank_vol, LD_ratio, tank_r)
print("Tank radius is {tank_r} m \nTank length is {tank_l} m".format(tank_r=tank_r, tank_l=tank_l))

# NETTING THEORY
stress_comp_endcaps = \
    functions.netting_analysis(weave_matrix, layers_matrix, t_ply, E_comp, p_prop, tank_r, safety_factor,
                               "end caps")
stress_comp_cylinder = \
    functions.netting_analysis(weave_matrix, layers_matrix, t_ply, E_comp, p_prop, tank_r, safety_factor,
                               "cylinder")
print("Composite stress in end caps is {stress_endcaps} MPa \n"
      "Composite stress in cylinder is {stress_cylinder} MPa"
      .format(stress_endcaps=stress_comp_endcaps/10**6, stress_cylinder=stress_comp_cylinder/10**6))

functions.graphical_analysis(stress_comp_cylinder, stress_comp_endcaps, sigma_comp, p_prop, safety_factor)


# TANK MASS CALCULATION
tank_mass = functions.mass_calc(layers_matrix, t_ply, t_lin, S_cylinder, S_endcaps, rho_comp, rho_lin)
print("Tank mass is {tank_mass} kg".format(tank_mass=tank_mass))
