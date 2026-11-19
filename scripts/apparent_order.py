#!/usr/bin/env python3
"""Observed orders from three unequal refinement levels, without extrapolation."""
import math
import sys


def difference_ratio(order, h):
    """(h0**p-h1**p)/(h1**p-h2**p), evaluated without large common factors."""
    left = math.log(h[0] / h[1])
    right = math.log(h[1] / h[2])
    if abs(order) < 1e-10:
        return (left / right) * (1 + .5 * order * (left + right))
    return math.expm1(order * left) / -math.expm1(-order * right)


def apparent_order(h, q, absolute_floor=0.0):
    """Require measurable consecutive differences; retain reversal/nonpositive flags."""
    if len(h) != 3 or len(q) != 3 or not all(math.isfinite(x) for x in (*h, *q)):
        raise ValueError('Three finite grid sizes and observables are required')
    if not h[0] > h[1] > h[2] > 0 or not math.isfinite(absolute_floor) or absolute_floor < 0:
        raise ValueError('Positive decreasing refinement scales and nonnegative floor required')
    delta = [q[0] - q[1], q[1] - q[2]]
    if not all(math.isfinite(x) for x in delta):
        raise ValueError('Consecutive differences must be representable')
    floor = max(absolute_floor, 64 * sys.float_info.epsilon * max(map(abs, q)))
    result = dict(h=h, observable=q, consecutive_differences=delta,
                  refinement_ratios=[h[0] / h[1], h[1] / h[2]], absolute_floor=floor,
                  apparent_order=None, magnitude_order_candidate=None,
                  interpretation='Observed order only; no resolved-other-axis or asymptotic proof.')
    if min(map(abs, delta)) <= floor:
        return result | dict(status='DIFFERENCE_BELOW_DECLARED_FLOOR')
    ratio = abs(delta[0] / delta[1]);result['difference_magnitude_ratio'] = ratio
    low, high = -32., 32.
    if not difference_ratio(low, h) <= ratio <= difference_ratio(high, h):
        return result | dict(status='ORDER_OUTSIDE_REPORTED_SEARCH_INTERVAL', search_interval=[low, high])
    for _ in range(100):
        middle = .5 * (low + high)
        if difference_ratio(middle, h) < ratio:
            low = middle
        else:
            high = middle
    candidate = .5 * (low + high);result['magnitude_order_candidate'] = candidate
    if (delta[0] < 0) != (delta[1] < 0):
        return result | dict(status='NONMONOTONIC_DIFFERENCES')
    result['apparent_order'] = candidate
    return result | dict(status='POSITIVE_OBSERVED_ORDER' if candidate > 0 else 'NONPOSITIVE_OBSERVED_ORDER')
