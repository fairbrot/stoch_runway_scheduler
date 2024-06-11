import random
import math
import numpy as np
import scipy.special as ss

def sample_pretac_delay(alpha: float, beta: float, a: float, h: float, x_bar: float) -> float:
    """
    Sample pre-tactical delay for given flight.

    Pretactical delay is Y - (a - h)
    where a is scheduled arrival time, h is time tactical uncertainty begins
    and Y is a Gamma distribution with shape and rate parameters alpha and beta.
    If alpha or beta parameters are not valid i.e. <= 0, then pre-tactical
    delay is a given constant x_bar (mean delay).

    See section 4 of 
    "A New Simheuristic Approach fo Stochastic Runway Scheduling" (2022) by Shone et al
    for more details.
    """
    if alpha > 0 and beta > 0:
        # Note that numpy is parameterised by shape=alpha and scale = 1/beta
        pretac_delay = np.random.gamma(alpha, 1/beta) - (a - h)
    else: # in this case there isn't a well defined Gamma distribution
        # In this case the pre-tactical delay is set equal to the average lateness rather than being sampled randomly.
        pretac_delay = x_bar
    return pretac_delay

def sample_cond_gamma(t: float, alpha: float, beta: float) -> float:
    """
    Sample from a Gamma distribution with shape `alpha` and rate `beta`, conditional on being
    above the value t.
    """
    q = ss.gdtr(1/beta, alpha, t) # gdtr is a fast function in scipy for evaluating cdf of gamma dist
    z = q + (1-q) * random.random() # Unform RN between q and 1
    return ss.gdtrix(1/beta, alpha, z) # gdtrix is fast function for quantile


def gamma_cond_exp(t, alpha, beta):

    # Site: https://stats.stackexchange.com/questions/338378/closed-form-conditional-expectation-of-gamma-distributed-variable

    # This is assuming mean is alpha/beta, variance is alpha/(beta^2)
    # t is the amount of time that the service has already been in progress

    denom=0
    for j in range(alpha):
        denom+=((beta*t)**j)/math.factorial(j)
    num=((beta*t)**alpha)/math.factorial(alpha)
    x=(alpha/beta)*(1+num/denom)

    return x
