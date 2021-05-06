"""This class defines two functions for calculating the concentration of a pollutant being released from a
smokestack. It uses the Gaussian Dispersion Model, which is the most common model used to estimate and predict pollution
concentrations at a given height and distance downstream of the source.

The final concentrations for each are expressed in micrograms/cubic meter (ug/m3).
All distances and heights are assumed to be in meters.
Emission rate (Q) assumes units of ug/s.
Wind speed (u) is the average wind speed of stack height (m/s).
x and y are measured along the the center plume line.
Effective stack height (h) is the sum of the physical height of the smoke stack and the plume height at a given time.
sy and sz are obtained from logarithmic graphs plotting distance downstream (x) against standard deviations of y and z,
respectively."""


import math


class pollutionconcentration:

    def __init__(self, emissionrate, windspeed, distancedownwind, horizontaldistance, verticalheight, stackheight,
                 plumerise, horizontaldispersion, verticaldispersion):
        self.Q = emissionrate
        self.u = windspeed
        self.x = distancedownwind
        self.y = horizontaldistance
        self.z = verticalheight
        self.h = stackheight + plumerise
        self.sy = horizontaldispersion
        self.sz = verticaldispersion

    def steadystateconcentration(self, Q, u, y, z, h, sy, sz):
        concentration = (Q / (2 * math.pi * u * sy * sz)) * math.exp((-1 / 2) * ((y * y) / (sy * sy))) * (
                math.exp((-0.5) * (((z - h) * (z - h)) / (sz * sz))) +
                math.exp((-0.5) * ((z + h) * (z + h)) / (sz * sz)))
        return concentration

    def maxgroundlevelconcentration(self, Q, u, sy, sz, h):
        concentration = (Q / (math.pi * u * sy * sz)) * math.exp(-0.5 * ((h * h) / (sz * sz)))
        return concentration
