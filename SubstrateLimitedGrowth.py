"""Monod Kinetics is a popular model for predicting the growth of bacterial cells based on substrate concentration,
maximum substrate utilization and the decay or death rate of the bacteria.
It is most often used in wastewater treatment to monitor food to mass ratio in the reaction basin and maintain
bacterial growth in the exponential phase where the most amount of substrate will be consumed.
The limiting substrate will depend on the strain of bacteria used in the treatment process, the chemical makeup of the
influent, and the treatment goals.
The values given here for x0, t and kd are example values to demonstrate solving the ordinary differential equation
for non-steady state growth in a system with continuous influx of substrate used in a flow through system."""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


class MonodKinetics:

    def __init__(self, specific_growth_rate, max_specific_growth_rate, microbial_decay_rate, maximum_growth_rate,
                 saturation_constant, substrate_concentration_in_solution, yield_coefficient, dissolved_oxygen_deficit,
                 initial_cell_concentration, cell_concentration):
        self.mu = specific_growth_rate
        self.mu_max = max_specific_growth_rate
        self.kd = microbial_decay_rate
        self.km = maximum_growth_rate
        self.ks = saturation_constant
        self.s = substrate_concentration_in_solution
        self.y = yield_coefficient
        self.d = dissolved_oxygen_deficit
        self.x0 = initial_cell_concentration
        self.x = cell_concentration

    @staticmethod
    def max_specific_growth_rate(y, km):
        mu_max = y*km
        return mu_max

    @staticmethod
    def single_substrate_limited_growth(mu_max, s, ks, kd):
        specific_growth_rate = ((mu_max*s)/(ks + s)) - kd
        return specific_growth_rate

    @staticmethod
    def non_steady_state_continuous_flow(d, mu, x):
        dxdt = (d*x0) + (mu-kd-d)*x
        return dxdt


x0 = (0, 50, 100)
kd = 0.05
t = np.linspace(0, 200, 500)
result = odeint(MonodKinetics.non_steady_state_continuous_flow, x0, t, args=(kd,))
fig, graphical_output = plt.subplots()
graphical_output.plot(t, result[:, 0], label="x0 = 0")
graphical_output.plot(t, result[:, 1], label="x0 = 50")
graphical_output.plot(t, result[:, 1], label="x0 = 100")
graphical_output.legend()
graphical_output.set_xlabel("t")
graphical_output.set_ylabel("x0")