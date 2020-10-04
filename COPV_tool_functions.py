from math import pi
import thermo as th
import numpy as np
import matplotlib.pyplot as plt


def properties_calc(names, fractions, temperature, mass, ullage):
    if fractions == "-":
        prop = th.chemical.Chemical(names, T=temperature)  # compute parameters
    else:
        mix_names = names.split(',')  # extract components of mixture
        mix_f = [float(f) for f in fractions.split(',')]  # extract mass fractions of components
        prop = th.mixture.Mixture(mix_names, ws=mix_f, T=temperature)  # compute mixture
    density = prop.rho  # calculate density
    volume = (mass/density) * (1 + ullage)  # calculate tank volume
    return density, volume


def geometry_gen(volume, ld, radius):
    if radius == "-":
        d = (volume / (pi * (ld / 4 - 1 / 12))) ** (1 / 3)  # diameter
        r = d / 2  # radius
        length = ld * d  # length
        s_cylinder = (length - 2 * r) * 2 * pi * r  # cylindrical surface
        s_endcaps = 4 * pi * r ** 2  # spherical end caps surface
    else:
        r = radius
        length = (volume - 4/3 * pi * r**3)/(pi*r**2)
        ld = length/(2*r)
        s_cylinder = length * 2*pi*r
        s_endcaps = 4*pi*r**2
    return r, length, s_cylinder, s_endcaps, ld


def netting_analysis(weave_matrix, layers_matrix, ply_thickness, e_modulus, pressure, radius, sf, section):
    if section == "end caps":
        weave = weave_matrix[0:2]
        layers = layers_matrix[0:2]
    else:
        weave = weave_matrix
        layers = layers_matrix
    a_matrix = 2*layers*ply_thickness
    k_1 = np.sum(a_matrix * (np.sin(weave)) ** 2 * (np.cos(weave)) ** 2)
    k_2 = np.sum(a_matrix * (np.cos(weave)) ** 4)
    k_3 = np.sum(a_matrix * (np.sin(weave)) ** 4)
    k_matrix = np.zeros((2, 2))
    k_matrix[0, 0] = k_3
    k_matrix[0, 1] = k_1
    k_matrix[1, 0] = k_1
    k_matrix[1, 1] = k_2
    k_inv = np.linalg.inv(k_matrix)
    if section == "end caps":
        strain_matrix = np.dot(k_inv, np.array([0.5, 0.5])) * sf * pressure * radius / e_modulus
    else:
        strain_matrix = np.dot(k_inv, np.array([1, 0.5]))*sf*pressure*radius/e_modulus
    print(f"Strain matrix in {section} is {strain_matrix}")
    stress_comp = e_modulus*np.add(strain_matrix[0]*np.sin(weave)**2, strain_matrix[1]*np.cos(weave)**2)
    return stress_comp


def graphical_analysis(stress_cylinder, stress_endcaps, sigma, pressure, sf):
    pressure /= 10**5
    sigma /= 10**6
    stress_cylinder /= 10**6
    stress_endcaps /= 10**6

    analysis = plt.figure()

    graph_cylinder = analysis.add_subplot(121, title="Composite stresses in the cylindrical part",
                                          xlabel="Tank pressure (bar)", ylabel="Stress (MPa)")
    xtab = [0, pressure * sf]
    ytab_hellical_1 = [0, stress_cylinder[0]]
    ytab_hellical_2 = [0, stress_cylinder[1]]
    ytab_hoop = [0, stress_cylinder[2]]
    ytab_sigma = [sigma, sigma]
    graph_cylinder.plot(xtab, ytab_hellical_1, label="Stress in 1st hellical layer")
    graph_cylinder.plot(xtab, ytab_hellical_2, label="Stress in 2nd hellical layer")
    graph_cylinder.plot(xtab, ytab_hoop, label="Stress in hoop layer")
    graph_cylinder.plot(xtab, ytab_sigma, label="Ultimate stress of the composite")
    graph_cylinder.legend()

    graph_endcaps = analysis.add_subplot(122, title="Composite stresses in the end caps", xlabel="Tank pressure (bar)",
                                         ylabel="Stress (MPa)")
    ytab_hellical_1 = [0, stress_endcaps[0]]
    ytab_hellical_2 = [0, stress_endcaps[1]]
    graph_endcaps.plot(xtab, ytab_hellical_1, label="Stress in 1st hellical layer")
    graph_endcaps.plot(xtab, ytab_hellical_2, label="Stress in 2nd hellical layer")
    graph_endcaps.plot(xtab, ytab_sigma, label="Ultimate stress of the composite")
    graph_endcaps.legend()

    plt.show()


def mass_calc(layers, ply_thickness, t_lin, s_cylinder, s_endcaps, rho_comp, rho_lin):
    t_cylinder = 2*np.sum(layers)*ply_thickness
    t_endcaps = 2*np.sum(layers[0:2])*ply_thickness
    tank_mass = s_cylinder*(t_cylinder*rho_comp + t_lin*rho_lin) + s_endcaps*(t_endcaps*rho_comp + t_lin*rho_lin)
    return tank_mass
