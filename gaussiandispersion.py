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


class SmokeStack:

    def __init__(self, emission_rate, wind_speed, distance_downwind, horizontal_distance, vertical_height, stack_height,
                 plume_rise, horizontal_dispersion, vertical_dispersion):
        self.q = emission_rate
        self.u = wind_speed
        self.x = distance_downwind
        self.y = horizontal_distance
        self.z = vertical_height
        self.h = stack_height + plume_rise
        self.sy = horizontal_dispersion
        self.sz = vertical_dispersion

    @classmethod
    def new_stack(cls, q, u, y, z, h, sy, sz):
        return cls(Q, u, y, z, h, sy, sz)

    @staticmethod
    def steady_state_concentration(q, u, y, z, h, sy, sz):
        concentration = (q / (2 * math.pi * u * sy * sz)) * math.exp((-0.5) * ((y * y) / (sy * sy))) * (
                math.exp((-0.5) * (((z - h) * (z - h)) / (sz * sz))) +
                math.exp((-0.5) * ((z + h) * (z + h)) / (sz * sz)))
        return concentration

    @staticmethod
    def max_ground_level_concentration(q, u, sy, sz, h):
        concentration = (q / (math.pi * u * sy * sz)) * math.exp(-0.5 * ((h * h) / (sz * sz)))
        return concentration





