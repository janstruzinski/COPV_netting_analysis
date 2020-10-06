from math import pi
import thermo as th
import numpy as np
import matplotlib.pyplot as plt


def properties_calc(components_names, components_fractions, fuel_temperature, fuel_mass, fuel_ullage):
    if components_fractions == "-":
        prop = th.chemical.Chemical(components_names, T=fuel_temperature)  # compute parameters
    else:
        mix_names = components_names.split(',')  # extract components of mixture
        mix_f = [float(f) for f in components_fractions.split(',')]  # extract mass fractions of components
        prop = th.mixture.Mixture(mix_names, ws=mix_f, T=fuel_temperature)  # compute mixture
    fuel_density = prop.rho  # calculate density
    fuel_volume = (fuel_mass / fuel_density) * (1 + fuel_ullage)  # calculate tank volume
    return fuel_density, fuel_volume


def geometry_gen(fuel_volume, ld_ratio, tank_radius):
    if tank_radius == "-":
        tank_diameter = (fuel_volume / (pi * (ld_ratio / 4 - 1 / 12))) ** (1 / 3)  # diameter
        tank_radius = tank_diameter / 2  # radius
        tank_length = ld_ratio * tank_diameter  # length
        surface_cylinder = (tank_length - 2 * tank_radius) * 2 * pi * tank_radius  # cylindrical surface
        surface_endcaps = 4 * pi * tank_radius ** 2  # spherical end caps surface
    else:
        tank_length = (fuel_volume - 4 / 3 * pi * tank_radius ** 3) / (pi * tank_radius ** 2)
        ld_ratio = tank_length / (2 * tank_radius)
        surface_cylinder = tank_length * 2 * pi * tank_radius
        surface_endcaps = 4 * pi * tank_radius ** 2
    return tank_radius, tank_length, surface_cylinder, surface_endcaps, ld_ratio


def netting_analysis(winding_angles, winding_layers, ply_thickness, cylinder_thickness, endcaps_thickness,
                     composite_e_modulus, burst_pressure, tank_radius, tank_section):
    if tank_section == "end caps":
        weave = winding_angles[0:2]
        layers = winding_layers[0:2]
    else:
        weave = winding_angles
        layers = winding_layers
    layer_thickness_matrix = 2 * layers * ply_thickness
    k_1 = np.sum(layer_thickness_matrix * (np.sin(weave)) ** 2 * (np.cos(weave)) ** 2)
    k_2 = np.sum(layer_thickness_matrix * (np.cos(weave)) ** 4)
    k_3 = np.sum(layer_thickness_matrix * (np.sin(weave)) ** 4)
    k_matrix = np.zeros((2, 2))
    k_matrix[0, 0] = k_3
    k_matrix[0, 1] = k_1
    k_matrix[1, 0] = k_1
    k_matrix[1, 1] = k_2
    k_inv = np.linalg.inv(k_matrix)
    if tank_section == "end caps":
        composite_strain_matrix = np.dot(k_inv, np.array([0.5, 0.5])) * burst_pressure * tank_radius / \
                                  composite_e_modulus
        apparent_e_modulus_matrix = tank_radius * burst_pressure * np.array([0.5, 0.5]) / \
                                    (composite_strain_matrix*endcaps_thickness)
    else:
        composite_strain_matrix = np.dot(k_inv, np.array([1, 0.5])) * burst_pressure * tank_radius / composite_e_modulus
        apparent_e_modulus_matrix = tank_radius * burst_pressure * np.array([1, 0.5]) \
                                    / (composite_strain_matrix*cylinder_thickness)
    print(f"Strain matrix in {tank_section} is {composite_strain_matrix}")
    composite_stress_matrix = composite_e_modulus * np.add(composite_strain_matrix[0] * np.sin(weave) ** 2,
                                                           composite_strain_matrix[1] * np.cos(weave) ** 2)
    print(f"Apparent E modulus matrix in {tank_section} is {apparent_e_modulus_matrix/10**9} GPa")
    return composite_stress_matrix, apparent_e_modulus_matrix


def residual_stress_analysis(apparent_composite_e_modulus_matrix, composite_e_modulus, composite_thickness,
                             liner_thickness, liner_e_modulus, liner_yield_strength, liner_poission_ratio, tank_radius,
                             pressure_proof, pressure_burst, winding_angles, tank_section):

    if tank_section == "cylinder":
        comp_stiffness_matrix = apparent_composite_e_modulus_matrix * composite_thickness * np.array([1, 2])/tank_radius
        liner_stiffness_matrix = liner_e_modulus * liner_thickness * \
                                 np.array([1/(1-liner_poission_ratio/2), 2/(-liner_poission_ratio+0.5)])/tank_radius
        total_stiffness_matrix = comp_stiffness_matrix+liner_stiffness_matrix
        yield_pressure = (liner_yield_strength * liner_thickness)/(tank_radius*(1-liner_poission_ratio/2))
        weave = winding_angles
    else:
        comp_stiffness_matrix = apparent_composite_e_modulus_matrix * composite_thickness * np.array([2, 2])/tank_radius
        liner_stiffness_matrix = liner_e_modulus * liner_thickness * np.array([2/(1-liner_poission_ratio),
                                                                               2/(1-liner_poission_ratio)])/tank_radius
        total_stiffness_matrix = comp_stiffness_matrix+liner_stiffness_matrix
        yield_pressure = (2 * liner_yield_strength * liner_thickness)/(tank_radius*(1-liner_poission_ratio))
        weave = winding_angles[0:2]

    beta_matrix = comp_stiffness_matrix/total_stiffness_matrix

    res_strain_matrix = (pressure_proof - yield_pressure)/comp_stiffness_matrix \
                             - pressure_proof/total_stiffness_matrix
    max_strain_matrix = (beta_matrix * pressure_burst)/comp_stiffness_matrix + res_strain_matrix
    res_stress_comp_matrix = composite_e_modulus * np.add(res_strain_matrix[0] * np.sin(weave) ** 2,
                                                          res_strain_matrix[1] * np.cos(weave) ** 2)
    max_stress_comp_matrix = composite_e_modulus * np.add(max_strain_matrix[0] * np.sin(weave) ** 2,
                                                          max_strain_matrix[1] * np.cos(weave) ** 2)

    residual_pressure_matrix = comp_stiffness_matrix * res_strain_matrix
    if tank_section == "cylinder":
        res_stress_lin_matrix = - residual_pressure_matrix * tank_radius/liner_thickness * np.array([1, 0.5])
        max_stress_lin_matrix = (1-beta_matrix)*pressure_burst*tank_radius/liner_thickness * np.array([1, 0.5]) \
                                    + res_stress_lin_matrix
    else:
        res_stress_lin_matrix = - residual_pressure_matrix * tank_radius / liner_thickness * np.array([0.5, 0.5])
        max_stress_lin_matrix = (1 - beta_matrix) * pressure_burst * tank_radius / liner_thickness \
                                * np.array([0.5, 0.5]) + res_stress_lin_matrix

    print(f"Residual strain matrix in {tank_section} is {res_strain_matrix}")
    print(f"Maximum strain matrix in {tank_section} is {max_strain_matrix}")
    print(f"Residual stress composite matrix in {tank_section} is {res_stress_comp_matrix / 10 ** 6} MPa")
    print(f"Maximum stress composite matrix in {tank_section} is {max_stress_comp_matrix / 10 ** 6} MPa")
    print(f"Residual stress liner matrix in {tank_section} is {res_stress_lin_matrix/10**6} MPa")
    print(f"Maximum stress liner matrix in {tank_section} is {max_stress_lin_matrix/10**6} MPa")
    return res_stress_comp_matrix, max_stress_comp_matrix, res_stress_lin_matrix, max_stress_lin_matrix


def graphical_analysis(composite_stress_cylinder, composite_stress_endcaps, composite_sigma_yield,
                       residual_stress_composite_cylinder, maximum_stress_composite_cylinder,
                       residual_stress_liner_cylinder, maximum_stress_liner_cylinder,
                       residual_stress_composite_endcaps, maximum_stress_composite_endcaps,
                       residual_stress_liner_endcaps, maximum_stress_liner_endcaps,
                       liner_sigma_yield, tank_pressure):
    tank_pressure /= 10 ** 5
    composite_sigma_yield /= 10 ** 6
    composite_stress_cylinder /= 10 ** 6
    composite_stress_endcaps /= 10 ** 6
    residual_stress_composite_cylinder /= 10 ** 6
    maximum_stress_composite_cylinder /= 10 ** 6
    residual_stress_liner_cylinder /= 10 ** 6
    maximum_stress_liner_cylinder /= 10 ** 6
    residual_stress_composite_endcaps /= 10 ** 6
    maximum_stress_composite_endcaps /= 10 ** 6
    residual_stress_liner_endcaps /= 10 ** 6
    maximum_stress_liner_endcaps /= 10 ** 6
    liner_sigma_yield /= 10 ** 6

    netting_copv_analysis = plt.figure()
    netting_copv_analysis.suptitle("Netting analysis of COPV winding")

    graph_netting_cylinder = netting_copv_analysis.add_subplot(121, title="Composite stresses in the cylindrical part",
                                                               xlabel="Tank pressure (bar)", ylabel="Stress (MPa)")
    xtab = [0, tank_pressure]
    ytab_hellical_1 = [0, composite_stress_cylinder[0]]
    ytab_hellical_2 = [0, composite_stress_cylinder[1]]
    ytab_hoop = [0, composite_stress_cylinder[2]]
    ytab_comp_yield = [composite_sigma_yield, composite_sigma_yield]
    graph_netting_cylinder.plot(xtab, ytab_hellical_1, label="Stress in 1st hellical layer")
    graph_netting_cylinder.plot(xtab, ytab_hellical_2, label="Stress in 2nd hellical layer")
    graph_netting_cylinder.plot(xtab, ytab_hoop, label="Stress in hoop layer")
    graph_netting_cylinder.plot(xtab, ytab_comp_yield, label="Ultimate stress of the composite")
    graph_netting_cylinder.legend()

    graph_netting_endcaps = netting_copv_analysis.add_subplot(122, title="Composite stresses in the end caps",
                                                              xlabel="Tank pressure (bar)", ylabel="Stress (MPa)")
    ytab_hellical_1 = [0, composite_stress_endcaps[0]]
    ytab_hellical_2 = [0, composite_stress_endcaps[1]]
    graph_netting_endcaps.plot(xtab, ytab_hellical_1, label="Stress in 1st hellical layer")
    graph_netting_endcaps.plot(xtab, ytab_hellical_2, label="Stress in 2nd hellical layer")
    graph_netting_endcaps.plot(xtab, ytab_comp_yield, label="Ultimate stress of the composite")
    graph_netting_endcaps.legend()

    residual_analysis = plt.figure()
    residual_analysis.suptitle("Analysis of residual and maximal stresses in COPV")

    graph_residual_cylinder = residual_analysis.add_subplot(121, title="Stresses in the cylindrical part",
                                                            xlabel="Tank pressure (bar)", ylabel="Stress (MPa)")
    ytab_hellical_1 = [residual_stress_composite_cylinder[0], maximum_stress_composite_cylinder[0]]
    ytab_hellical_2 = [residual_stress_composite_cylinder[1], maximum_stress_composite_cylinder[1]]
    ytab_comp_hoop = [residual_stress_composite_cylinder[2], maximum_stress_composite_cylinder[2]]
    ytab_liner_hoop = [residual_stress_liner_cylinder[0], maximum_stress_liner_cylinder[0]]
    ytab_liner_longitidunal = [residual_stress_liner_cylinder[1], maximum_stress_liner_cylinder[1]]
    ytab_liner_yield = [liner_sigma_yield, liner_sigma_yield]
    graph_residual_cylinder.plot(xtab, ytab_hellical_1, label="Stress in 1st hellical layer")
    graph_residual_cylinder.plot(xtab, ytab_hellical_2, label="Stress in 2nd hellical layer")
    graph_residual_cylinder.plot(xtab, ytab_comp_hoop, label="Stress in hoop layer")
    graph_residual_cylinder.plot(xtab, ytab_liner_hoop, label="Stress in liner, hoop direction")
    graph_residual_cylinder.plot(xtab, ytab_liner_longitidunal, label="Stress in liner, longitudinal direction")
    graph_residual_cylinder.plot(xtab, ytab_comp_yield, label="Ultimate stress of the composite")
    graph_residual_cylinder.plot(xtab, ytab_liner_yield, label="Yield stress of liner")
    graph_residual_cylinder.legend()

    graph_residual_endcap = residual_analysis.add_subplot(122, title="Stresses in the end caps",
                                                          xlabel="Tank pressure (bar)", ylabel="Stress (MPa)")
    ytab_hellical_1 = [residual_stress_composite_endcaps[0], maximum_stress_composite_endcaps[0]]
    ytab_hellical_2 = [residual_stress_composite_endcaps[1], maximum_stress_composite_endcaps[1]]
    ytab_liner_hoop = [residual_stress_liner_endcaps[0], maximum_stress_liner_endcaps[0]]
    ytab_liner_longitidunal = [residual_stress_liner_endcaps[1], maximum_stress_liner_endcaps[1]]
    graph_residual_endcap.plot(xtab, ytab_hellical_1, label="Stress in 1st hellical layer")
    graph_residual_endcap.plot(xtab, ytab_hellical_2, label="Stress in 2nd hellical layer")
    graph_residual_endcap.plot(xtab, ytab_liner_hoop, label="Stress in liner, hoop direction")
    graph_residual_endcap.plot(xtab, ytab_liner_longitidunal, label="Stress in liner, longitudinal direction")
    graph_residual_endcap.plot(xtab, ytab_comp_yield, label="Ultimate stress of the composite")
    graph_residual_endcap.plot(xtab, ytab_liner_yield, label="Yield stress of liner")
    graph_residual_endcap.legend()

    plt.show()


def mass_calc(thickness_cylinder, thickness_endcaps, thickness_liner, surface_cylinder, surface_endcaps,
              composite_density, liner_density):
    tank_mass = surface_cylinder * (thickness_cylinder * composite_density + thickness_liner * liner_density) \
                + surface_endcaps * (thickness_endcaps * composite_density + thickness_liner * liner_density)
    return tank_mass


