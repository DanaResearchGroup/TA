"""Tests for the TemperatureProgram class."""

import pytest

from ta.temperature_program import TemperatureProgram, TemperatureSegment


# ------------------------------------------------------------------
# Construction & properties
# ------------------------------------------------------------------
class TestConstruction:

    def test_single_ramp(self):
        prog = TemperatureProgram.single_ramp(50.0, 600.0, 10.0)

        assert prog.T_initial_C == pytest.approx(50.0)
        assert prog.T_final_C == pytest.approx(600.0)
        # 550 deg / 10 deg/min = 55 min = 3300 s
        assert prog.total_time_s == pytest.approx(3300.0)
        assert len(prog.segments) == 1

    def test_ramp_hold_ramp(self):
        prog = TemperatureProgram(50.0, [
            TemperatureSegment(10.0, 200.0),            # 15 min
            TemperatureSegment(0.0, 200.0, hold_min=30),  # 30 min hold
            TemperatureSegment(5.0, 600.0),             # 80 min
        ])

        assert prog.T_initial_C == pytest.approx(50.0)
        assert prog.T_final_C == pytest.approx(600.0)
        expected_time = (15 + 30 + 80) * 60.0
        assert prog.total_time_s == pytest.approx(expected_time)
        assert len(prog.segments) == 3

    def test_cooling_segment(self):
        prog = TemperatureProgram(50.0, [
            TemperatureSegment(10.0, 200.0),   # heat to 200
            TemperatureSegment(-5.0, 100.0),   # cool to 100
        ])

        assert prog.T_initial_C == pytest.approx(50.0)
        assert prog.T_final_C == pytest.approx(100.0)
        # 15 min ramp + 20 min cool = 35 min
        assert prog.total_time_s == pytest.approx(35 * 60.0)

    def test_segment_boundaries(self):
        prog = TemperatureProgram(50.0, [
            TemperatureSegment(10.0, 200.0),
            TemperatureSegment(0.0, 200.0, hold_min=10),
            TemperatureSegment(5.0, 300.0),
        ])

        bounds = prog.segment_boundaries
        assert len(bounds) == 4  # initial + 3 segments

        # Initial
        assert bounds[0] == pytest.approx((0.0, 50.0))
        # After ramp: 15 min = 900 s
        assert bounds[1] == pytest.approx((900.0, 200.0))
        # After hold: 900 + 600 = 1500 s
        assert bounds[2] == pytest.approx((1500.0, 200.0))
        # After second ramp: 1500 + 20*60 = 2700 s
        assert bounds[3] == pytest.approx((2700.0, 300.0))


# ------------------------------------------------------------------
# Temperature lookup
# ------------------------------------------------------------------
class TestTemperatureLookup:

    def test_single_ramp_temperatures(self):
        prog = TemperatureProgram.single_ramp(50.0, 200.0, 10.0)

        assert prog.T_furnace_C(0.0) == pytest.approx(50.0)
        # At 5 min = 300 s: 50 + 10*5 = 100
        assert prog.T_furnace_C(300.0) == pytest.approx(100.0)
        # At 15 min = 900 s: 50 + 10*15 = 200
        assert prog.T_furnace_C(900.0) == pytest.approx(200.0)

    def test_hold_segment_flat(self):
        prog = TemperatureProgram(100.0, [
            TemperatureSegment(10.0, 200.0),             # 0-600 s
            TemperatureSegment(0.0, 200.0, hold_min=10),  # 600-1200 s
        ])

        # During hold: T stays at 200
        assert prog.T_furnace_C(700.0) == pytest.approx(200.0)
        assert prog.T_furnace_C(900.0) == pytest.approx(200.0)
        assert prog.T_furnace_C(1200.0) == pytest.approx(200.0)

    def test_cooling_segment_decreasing(self):
        prog = TemperatureProgram(200.0, [
            TemperatureSegment(-10.0, 100.0),  # 10 min = 600 s
        ])

        assert prog.T_furnace_C(0.0) == pytest.approx(200.0)
        # At 5 min = 300 s: 200 - 10*5 = 150
        assert prog.T_furnace_C(300.0) == pytest.approx(150.0)
        assert prog.T_furnace_C(600.0) == pytest.approx(100.0)

    def test_past_end_returns_final_temperature(self):
        prog = TemperatureProgram.single_ramp(50.0, 200.0, 10.0)

        # Well past the end
        assert prog.T_furnace_C(99999.0) == pytest.approx(200.0, abs=0.1)

    def test_negative_time_returns_initial(self):
        prog = TemperatureProgram.single_ramp(50.0, 200.0, 10.0)

        assert prog.T_furnace_C(-1.0) == pytest.approx(50.0)

    def test_T_furnace_K(self):
        prog = TemperatureProgram.single_ramp(50.0, 200.0, 10.0)

        assert prog.T_furnace_K(0.0) == pytest.approx(50.0 + 273.15)

    def test_multi_segment_continuity(self):
        """Temperature should be continuous across segment boundaries."""
        prog = TemperatureProgram(50.0, [
            TemperatureSegment(10.0, 200.0),
            TemperatureSegment(0.0, 200.0, hold_min=10),
            TemperatureSegment(5.0, 300.0),
        ])

        bounds = prog.segment_boundaries
        for t_s, T_C in bounds:
            assert prog.T_furnace_C(t_s) == pytest.approx(T_C, abs=0.01)


# ------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------
class TestValidation:

    def test_empty_segments_raises(self):
        with pytest.raises(ValueError, match="At least one segment"):
            TemperatureProgram(50.0, [])

    def test_hold_with_zero_duration_raises(self):
        with pytest.raises(ValueError, match="hold_min > 0"):
            TemperatureProgram(50.0, [
                TemperatureSegment(0.0, 50.0, hold_min=0),
            ])

    def test_positive_rate_cooling_raises(self):
        with pytest.raises(ValueError, match="rate is positive"):
            TemperatureProgram(200.0, [
                TemperatureSegment(10.0, 100.0),  # positive rate, cooling direction
            ])

    def test_negative_rate_heating_raises(self):
        with pytest.raises(ValueError, match="rate is negative"):
            TemperatureProgram(50.0, [
                TemperatureSegment(-10.0, 200.0),  # negative rate, heating direction
            ])

    def test_ramp_to_same_temp_raises(self):
        with pytest.raises(ValueError, match="equals current temperature"):
            TemperatureProgram(200.0, [
                TemperatureSegment(10.0, 200.0),  # rate != 0 but no temperature change
            ])
