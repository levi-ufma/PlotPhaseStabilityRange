# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 13:27:20 2018
@author: jherfson

Updated on Mon Nov 26 21:30:00 2018
@author: Rodolpho
"""

import math

"""
This class calculates entropies and enthalpies of oxygen gas, based on fits of experimental data, 
and then uses them to estimate iteratively the temperature corresponding to a given oxygen chemical 
potential at partial pressure of 0.21 atm. 
"""


class Thermodynamics:
    def __init__(self, *args, **kwargs):
        # Shomate equation coefficients A, B, C, D, E, F, G for 0.1 MPa, fitted from standard experimental data.
        # Source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1
        self.c1 = (31.32234, -20.23531, 57.86644, -36.50624, -
                   0.007374, -8.903471, 246.7945)  # Range: 100 K - 700 K
        self.c2 = (30.03235, 8.772972, -3.988133, 0.788313, -
                   0.741599, -11.32468, 236.1663)  # Range: 700 K - 2000 K
        self.c3 = (20.91111, 10.72071, -2.020498, 0.146449, 9.245722,
                   5.337651, 237.6185)    # Range: 2000 K - 6000 K
        # Conversion factors
        self.kJmol_to_eV = 0.01036427230133138
        self.JmolK_to_eVK = 0.00001036427230133138
        self.atm_to_MPa = float(101325.0e-6)
        # Boltzmann's constantin eV
        self.kB = float(8.6173303e-5)
        # Energy per atom at 0 K in eV for oxygen, taken from the Materials Project website (ID: mp-12957). Calculated via DFT
        self.Ho_eV = -4.93552791875*2.0
        # Difference between the enthalpy at 298.15 K and the one at 0 K, in kJ/mol
        # Source: Malcolm W. Chase Jr., NIST-JANAF Themochemical Tables, Fourth Edition, J. Phys. Chem. Ref. Data, Monograph 9, 1998, pp. 1744
        self.HrefMinusHo_eV = float(8.683)*self.kJmol_to_eV
        # Enthalpy of oxygen gas at 298.15 K in eV, assuming the Materials Project value for the enthalpy at 0 K.
        self.Href = self.HrefMinusHo_eV + self.Ho_eV
        # Pressure in atm. Change it to the desired value.
        self.pressure = 0.21

    # This function returns the entropy of oxygen gas at 0.1 MPa and given temperature in eV
    def entropy(self, temperature):

        if temperature < 100.0:
            temperature = 100.0
        if temperature > 6000.0:
            temperature = 6000.0

        # Use the appropriate Shomate equation coefficients for the given temperature
        if temperature >= 100 and temperature <= 700:
            A = self.c1[0]
            B = self.c1[1]
            C = self.c1[2]
            D = self.c1[3]
            E = self.c1[4]
            G = self.c1[6]
        elif temperature > 700 and temperature <= 2000:
            A = self.c2[0]
            B = self.c2[1]
            C = self.c2[2]
            D = self.c2[3]
            E = self.c2[4]
            G = self.c2[6]
        elif temperature > 2000 and temperature <= 6000:
            A = self.c3[0]
            B = self.c3[1]
            C = self.c3[2]
            D = self.c3[3]
            E = self.c3[4]
            G = self.c3[6]

        t = temperature/1000.0
        S_JmolK = A*math.log(t) + B*t + C*t**2/2.0 + \
            D*t**3/3.0 - E/2.0/t**2 + G
        S_eVK = S_JmolK*self.JmolK_to_eVK

        # print(S_eVK/self.kJmol_to_eV)

        return S_eVK

    # This function returns the enthalpy of oxygen gas at given temperature minus the one at 298.15 K, in eV
    def enthalpy_minusRef(self, temperature):

        if temperature < 100.0:
            temperature = 100.0
        if temperature > 6000.0:
            temperature = 6000.0

        # Use the appropriate Shomate equation coefficients for the given temperature
        if temperature >= 100 and temperature <= 700:
            A = self.c1[0]
            B = self.c1[1]
            C = self.c1[2]
            D = self.c1[3]
            E = self.c1[4]
            F = self.c1[5]
        elif temperature > 700 and temperature <= 2000:
            A = self.c2[0]
            B = self.c2[1]
            C = self.c2[2]
            D = self.c2[3]
            E = self.c2[4]
            F = self.c2[5]
        elif temperature > 2000 and temperature <= 6000:
            A = self.c3[0]
            B = self.c3[1]
            C = self.c3[2]
            D = self.c3[3]
            E = self.c3[4]
            F = self.c3[5]

        t = temperature/1000.0
        HMinusHref_kJmolK = A*t + B*t**2/2.0 + C*t**3/3.0 + D*t**4/4.0 - E/t + F
        HMinusHref_eV = HMinusHref_kJmolK*self.kJmol_to_eV

        # print(HMinusHref_eV/self.kJmol_to_eV)
        return HMinusHref_eV

    # This function calculates iteratively the temperature corresponding to a given oxygen chemical potential
    def mu_to_temperature(self, mu):
        # We need an initial guess for temperature... Why not RT? =)
        temperature = 298.0
        # Just a dummy variable for the iterations
        temperature_old = 0.0

        # Do this until the temperature converges within 0.001 K
        while abs(temperature-temperature_old) > 1.0e-3:

            temperature_old = temperature
            # The following expression was obtained by solving the equation below for T, with p = 0.21 atm and p0 = 0.1 MPa
            # u(T,p)=[h(T)-h(Tref)]+h(Tref)-T*s(T,p0)+kB*T*ln(p/p0)
            temperature = (2*mu - self.enthalpy_minusRef(temperature_old) - self.Href)/(
                self.kB*math.log(self.pressure*self.atm_to_MPa/0.1) - self.entropy(temperature_old))

            # No negative temperatures!
            if temperature < 0.0:
                temperature = 0.0

        return temperature

    def temperature_to_mu(self, temperature):

        oxygen_enthalpy = self.Ho_eV + self.HrefMinusHo_eV + \
            self.enthalpy_minusRef(temperature)

        mu = oxygen_enthalpy - temperature * \
            self.entropy(temperature) + temperature*self.kB * \
            math.log(self.pressure*self.atm_to_MPa/0.1)

        return mu

    def print_temperature_corresponding_to_mu_equals(self, mu):

        return {
            'ChemPot_eV': round(mu, 3),
            'T_Celsius': round(self.mu_to_temperature(mu)-273.0, 2),
            'T_Kelvin': round(self.mu_to_temperature(mu), 2)
        }

    def print_mu_corresponding_to_temperature_equals(self, temperature):
        return {
            'ChemPot_eV': round(self.temperature_to_mu(temperature), 3),
            'T_Celsius': round(temperature-273.0, 2),
            'T_Kelvin': round(temperature, 3)
        }


# if __name__ == "__main__":
#     mu = Thermodynamics()
#     print('temperature_to_mu: ', mu.print_mu_corresponding_to_temperature_equals(1623))
#     print('mu_to_temperature: ', mu.print_temperature_corresponding_to_mu_equals(-5.323))
