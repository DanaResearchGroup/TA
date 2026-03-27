"""Multi-segment temperature programs for thermal-analysis experiments.

A :class:`TemperatureProgram` describes the furnace temperature as a
piecewise-linear function of time, built from a sequence of
:class:`TemperatureSegment` objects (ramps, holds, and cooling segments).
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass


@dataclass(frozen=True)
class TemperatureSegment:
    """One segment of a temperature program.

    Parameters
    ----------
    rate_C_per_min : float
        Heating rate [deg C / min].  Positive = heating, negative =
        cooling, zero = isothermal hold.
    T_target_C : float
        Target temperature [deg C] for ramp segments.  For isothermal
        holds this should equal the current temperature (it is stored
        for documentation but not used to compute the profile).
    hold_min : float
        Duration of an isothermal hold [min].  Only used when
        *rate_C_per_min* is zero.
    """

    rate_C_per_min: float
    T_target_C: float
    hold_min: float = 0.0


class TemperatureProgram:
    """Piecewise-linear furnace temperature profile.

    The program starts at *T_initial_C* and executes each segment in
    order.  Ramp segments heat or cool at a constant rate until the
    target temperature is reached; hold segments keep the temperature
    constant for a specified duration.

    Parameters
    ----------
    T_initial_C : float
        Starting temperature [deg C].
    segments : list[TemperatureSegment]
        Ordered list of program segments.

    Examples
    --------
    Ramp from 50 to 200 at 10 C/min, hold 30 min, ramp to 600 at 5 C/min:

    >>> prog = TemperatureProgram(50.0, [
    ...     TemperatureSegment(10.0, 200.0),
    ...     TemperatureSegment(0.0, 200.0, hold_min=30.0),
    ...     TemperatureSegment(5.0, 600.0),
    ... ])
    """

    def __init__(
        self,
        T_initial_C: float,
        segments: list[TemperatureSegment],
    ) -> None:
        if not segments:
            raise ValueError("At least one segment is required")

        self._T_initial_C = T_initial_C
        self._segments = list(segments)

        # Pre-compute boundary arrays for fast lookup.
        # _t_bounds[i]   = cumulative time [s] at the START of segment i
        # _T_bounds_K[i] = temperature [K] at the START of segment i
        # _rates_K_per_s[i] = rate [K/s] during segment i (0 for holds)
        t_bounds: list[float] = [0.0]
        T_bounds_K: list[float] = [T_initial_C + 273.15]
        rates_K_per_s: list[float] = []

        T_cur_K = T_initial_C + 273.15
        t_cur = 0.0

        for i, seg in enumerate(segments):
            if seg.rate_C_per_min == 0.0:
                # Isothermal hold
                if seg.hold_min <= 0:
                    raise ValueError(
                        f"Segment {i}: isothermal hold requires hold_min > 0"
                    )
                duration_s = seg.hold_min * 60.0
                rates_K_per_s.append(0.0)
                t_cur += duration_s
                # Temperature unchanged
            else:
                # Ramp (heating or cooling)
                T_target_K = seg.T_target_C + 273.15
                direction = T_target_K - T_cur_K

                if abs(direction) < 1e-10:
                    raise ValueError(
                        f"Segment {i}: T_target equals current temperature; "
                        f"use rate_C_per_min=0 for an isothermal hold"
                    )
                if direction > 0 and seg.rate_C_per_min < 0:
                    raise ValueError(
                        f"Segment {i}: rate is negative but "
                        f"T_target ({seg.T_target_C}) > T_current "
                        f"({T_cur_K - 273.15:.2f})"
                    )
                if direction < 0 and seg.rate_C_per_min > 0:
                    raise ValueError(
                        f"Segment {i}: rate is positive but "
                        f"T_target ({seg.T_target_C}) < T_current "
                        f"({T_cur_K - 273.15:.2f})"
                    )

                rate_K_per_s = seg.rate_C_per_min / 60.0
                duration_s = abs(direction) / abs(rate_K_per_s)
                rates_K_per_s.append(rate_K_per_s)
                t_cur += duration_s
                T_cur_K = T_target_K

            t_bounds.append(t_cur)
            T_bounds_K.append(T_cur_K)

        self._t_bounds = t_bounds
        self._T_bounds_K = T_bounds_K
        self._rates_K_per_s = rates_K_per_s
        self._total_time_s = t_cur

    # ------------------------------------------------------------------
    # Factories
    # ------------------------------------------------------------------
    @classmethod
    def single_ramp(
        cls,
        T_initial_C: float,
        T_final_C: float,
        rate_C_per_min: float,
    ) -> TemperatureProgram:
        """Create a program consisting of a single linear ramp.

        This reproduces the legacy ``TASimulator`` behaviour.
        """
        return cls(T_initial_C, [TemperatureSegment(rate_C_per_min, T_final_C)])

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    @property
    def T_initial_C(self) -> float:
        """Starting temperature [deg C]."""
        return self._T_initial_C

    @property
    def T_final_C(self) -> float:
        """Temperature at the end of the last segment [deg C]."""
        return self._T_bounds_K[-1] - 273.15

    @property
    def total_time_s(self) -> float:
        """Total duration of the program [s]."""
        return self._total_time_s

    @property
    def segments(self) -> list[TemperatureSegment]:
        """Copy of the segment list."""
        return list(self._segments)

    @property
    def segment_boundaries(self) -> list[tuple[float, float]]:
        """``(time_s, temperature_C)`` at each segment boundary.

        Includes the initial point and the end of every segment,
        so the length is ``len(segments) + 1``.
        """
        return [
            (t, T_K - 273.15)
            for t, T_K in zip(self._t_bounds, self._T_bounds_K)
        ]

    # ------------------------------------------------------------------
    # Temperature lookup
    # ------------------------------------------------------------------
    def T_furnace_K(self, t: float) -> float:
        """Furnace temperature [K] at time *t* [s]."""
        if t <= 0.0:
            return self._T_bounds_K[0]
        if t >= self._total_time_s:
            return self._T_bounds_K[-1]

        idx = bisect.bisect_right(self._t_bounds, t) - 1
        idx = max(0, min(idx, len(self._segments) - 1))

        t_seg = t - self._t_bounds[idx]
        return self._T_bounds_K[idx] + self._rates_K_per_s[idx] * t_seg

    def T_furnace_C(self, t: float) -> float:
        """Furnace temperature [deg C] at time *t* [s]."""
        return self.T_furnace_K(t) - 273.15
