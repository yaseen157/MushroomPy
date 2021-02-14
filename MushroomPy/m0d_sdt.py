# -*- coding: utf-8 -*-

"""
MushroomPy Module

Generate von-Neumann and Chapman-Jouget pressure data for zero-dimensional analysis.

Theory, numerical methods and applications are described in the following report:

    Numerical Solution Methods for Shock and Detonation Jump Conditions, S.
    Browne, J. Ziegler, and J. E. Shepherd, GALCIT Report FM2006.006 - R3,
    California Institute of Technology Revised September, 2018

Updated December 2020
Tested with:
    Python 3.8, Windows 10

Author Email: yr3g17@soton.ac.uk
"""

__author__ = "Yaseen Reza"

from sdtoolbox.postshock import CJspeed, PostShock_fr, PostShock_eq


class ShockDetonationToolbox:
    """**Zero-dimensional** detonation analysis tool, using methods from the **Shock & Detonation Toolbox**."""

    def __init__(self, p_initial, t_initial, q_initial):
        self.P1 = p_initial
        self.T1 = t_initial
        self.q = q_initial
        self.mech = 'gri30_highT.cti'

        self.cjspeed = False

    def calc_cjspeed(self):

        if self.cjspeed is False:
            # Calculation for the Chapman-Jouget speeds
            self.cjspeed = CJspeed(P1=self.P1, T1=self.T1, q=self.q, mech=self.mech, fullOutput=False)

        return self.cjspeed

    def calc_cjpressure(self):
        # Calculation for the Chapman-Jouget speeds
        cj_speed_mps = self.calc_cjspeed()

        # Post shock (equilibrium) pressure to be determined
        gascj = PostShock_eq(U1=cj_speed_mps, P1=self.P1, T1=self.T1, q=self.q, mech=self.mech)
        cj_pressure_pa = gascj.P

        return cj_pressure_pa

    def calc_vnpressure(self):
        # Calculation for the Chapman-Jouget speeds
        cj_speed_mps = self.calc_cjspeed()

        # Post shock (frozen) pressure to be determined
        gasvn = PostShock_fr(U1=cj_speed_mps, P1=self.P1, T1=self.T1, q=self.q, mech=self.mech)
        vn_pressure_pa = gasvn.P

        return vn_pressure_pa

class SDTTools:

    def __init__():
        return

    def plot_cjcontour():
        return

if __name__ == "__main__":
    test0d = ShockDetonationToolbox(p_initial=100000, t_initial=298, q_initial="C2H4:1 O2:3")

    print(test0d.calc_cjpressure())

    pass
