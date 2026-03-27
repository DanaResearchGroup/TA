"""Cantera reactor-network logic for TGA/DTA/DTG/MS simulation.

Uses a :class:`~cantera.ConstPressureReactor` (sample pan) connected to a
:class:`~cantera.Reservoir` (furnace) via a :class:`~cantera.Wall` with a
finite heat-transfer coefficient.  The furnace follows a prescribed
temperature program (one or more ramp / hold / cooling segments); the
sample heats through the wall, producing realistic thermal lag and
enabling proper DTA measurement.

The furnace profile is implemented via a time-dependent heat-flux
function (:class:`~cantera.Func1`) on the wall, so that the reservoir
temperature stays fixed at the initial value while the *effective*
furnace temperature follows the program.  This avoids the Cantera
limitation that a :class:`~cantera.Reservoir` state cannot be updated
mid-integration.

This is a gas-phase approximation using **phase filtering**: species are
partitioned into *condensed* (tracked for TGA) and *volatile* (tracked for
MS), but all live in the same Cantera ``Solution``.  For full multiphase
polymer modelling a Cantera ``Interface`` / bulk-phase mechanism is required.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import cantera as ct
import numpy as np

from ta.temperature_program import TemperatureProgram, TemperatureSegment


@dataclass
class SimulationResult:
    """Container for raw simulation output.

    All arrays share the same length (number of recorded time steps).

    Attributes
    ----------
    time_s : np.ndarray
        Time in seconds.
    furnace_temperature_C : np.ndarray
        Prescribed furnace (reservoir) temperature in degree C.
    reactor_temperature_C : np.ndarray
        Actual sample (reactor) temperature in degree C.
    reactor_mass : np.ndarray
        Total gas mass inside the reactor [kg] at each step.
    mass_fractions : dict[str, np.ndarray]
        Per-species mass-fraction history Y_k(t).
    mole_fractions : dict[str, np.ndarray]
        Per-species mole-fraction history X_k(t).
    wall_heat_flux : np.ndarray
        Wall heat-transfer rate [W] (positive = into reactor).
    density : np.ndarray
        Gas density [kg/m^3].
    species_names : list[str]
        All species in the mechanism.
    condensed_species : list[str]
        Species classified as condensed (tracked for TGA).
    sample_mass_mg : float
        User-specified initial sample mass [mg].
    bath_gas : str
        Bath gas species name.
    bath_gas_flow_sccm : float
        Bath gas volumetric flow rate [standard cm^3/min].
    pressure_Pa : float
        System pressure [Pa].
    temperature_program : TemperatureProgram or None
        The temperature program used for this simulation, if available.
    """

    time_s: np.ndarray
    furnace_temperature_C: np.ndarray
    reactor_temperature_C: np.ndarray
    reactor_mass: np.ndarray
    mass_fractions: dict[str, np.ndarray]
    mole_fractions: dict[str, np.ndarray]
    wall_heat_flux: np.ndarray
    density: np.ndarray
    species_names: list[str]
    condensed_species: list[str]
    sample_mass_mg: float
    bath_gas: str
    bath_gas_flow_sccm: float
    pressure_Pa: float
    temperature_program: TemperatureProgram | None = None


class TASimulator:
    """Simulate simultaneous TGA/DTA/DTG/MS using Cantera.

    The experiment consists of a sample crucible (``ConstPressureReactor``)
    inside a furnace (``Reservoir``) connected by a ``Wall`` with heat-transfer
    coefficient *U* and area *A*.  The furnace follows a temperature program
    (one or more ramp / hold / cooling segments) while the sample heats
    through the wall, producing thermal lag.

    Implementation
    --------------
    The furnace profile is realised as a time-dependent wall heat-flux
    ``Func1`` so that the single ``ReactorNet.advance()`` loop handles
    both chemistry and heat transfer without re-creating the network.

    At each recording time-step the simulator:

    1. Advances the ``ReactorNet`` to the current time (Cantera solves
       the coupled energy + species ODEs internally).
    2. Records the reactor state and wall heat-transfer rate.

    The energy equation is **enabled** on the reactor so that the DTA signal
    (temperature difference driven by endo/exothermic reactions and wall heat
    transfer) emerges naturally.

    Parameters
    ----------
    mechanism : str
        Path to a Cantera YAML mechanism file (or a built-in name
        such as ``"gri30.yaml"``).
    sample_composition : str or dict
        Initial composition, e.g. ``"C2H6:1.0"`` or ``{"C2H6": 1.0}``.
    sample_mass_mg : float
        Nominal sample mass [mg] (default 10).  Used to set the reactor
        volume so that the simulated mass matches the physical sample.
    bath_gas : str
        Bath-gas species name (default ``"AR"``).
    bath_gas_flow_sccm : float
        Bath-gas volumetric flow rate at STP [cm^3/min].
    T_initial_C : float
        Starting temperature [deg C] (default 50).  Ignored when
        *temperature_program* is provided.
    T_final_C : float
        Target temperature [deg C] (default 600).  Ignored when
        *temperature_program* is provided.
    heating_rate_C_per_min : float
        Linear heating rate [deg C / min] (default 5).  Ignored when
        *temperature_program* is provided.
    pressure_atm : float
        System pressure [atm] (default 1).
    condensed_species : list[str] or None
        Species that remain in the crucible.  If *None*, every species
        with a non-zero initial mass fraction (excluding the bath gas)
        is treated as condensed.
    dt : float
        Recording time-step [s] (default 1).
    wall_area : float
        Wall area between furnace and sample [m^2] (default 1e-4).
    htc : float
        Wall heat-transfer coefficient [W/m^2/K] (default 100).
    temperature_program : TemperatureProgram or None
        Multi-segment temperature program.  When provided, overrides
        *T_initial_C*, *T_final_C*, and *heating_rate_C_per_min*.
    """

    def __init__(
        self,
        mechanism: str = "gri30.yaml",
        sample_composition: str | dict[str, float] = "CH4:1",
        sample_mass_mg: float = 10.0,
        bath_gas: str = "AR",
        bath_gas_flow_sccm: float = 50.0,
        T_initial_C: float = 50.0,
        T_final_C: float = 600.0,
        heating_rate_C_per_min: float = 5.0,
        pressure_atm: float = 1.0,
        condensed_species: list[str] | None = None,
        dt: float = 1.0,
        wall_area: float = 1e-4,
        htc: float = 100.0,
        temperature_program: TemperatureProgram | None = None,
    ) -> None:
        self.mechanism = mechanism
        self.sample_composition = sample_composition
        self.sample_mass_mg = sample_mass_mg
        self.bath_gas = bath_gas
        self.bath_gas_flow_sccm = bath_gas_flow_sccm
        self.dt = dt
        self.wall_area = wall_area
        self.htc = htc

        # --- temperature program ---
        if temperature_program is not None:
            self.program = temperature_program
        else:
            self.program = TemperatureProgram.single_ramp(
                T_initial_C, T_final_C, heating_rate_C_per_min,
            )

        # --- derived quantities ---
        self.T_initial_K = self.program.T_initial_C + 273.15
        self.pressure_Pa = pressure_atm * ct.one_atm
        self.t_final_s = self.program.total_time_s
        self.n_steps = int(np.ceil(self.t_final_s / self.dt))

        # --- probe mechanism for species list ---
        gas = ct.Solution(mechanism)
        gas.TPX = self.T_initial_K, self.pressure_Pa, sample_composition
        self.species_names: list[str] = list(gas.species_names)

        # validate bath gas
        if bath_gas not in self.species_names:
            raise ValueError(
                f"Bath gas '{bath_gas}' not found in mechanism "
                f"'{mechanism}'.  Available (first 10): "
                f"{self.species_names[:10]}"
            )

        # --- condensed-species list ---
        if condensed_species is not None:
            missing = [
                s for s in condensed_species if s not in self.species_names
            ]
            if missing:
                raise ValueError(
                    f"Condensed species not in mechanism: {missing}"
                )
            self.condensed_species = list(condensed_species)
        else:
            # default: every species with Y > 0 that is not the bath gas
            self.condensed_species = [
                name
                for name, y in zip(gas.species_names, gas.Y)
                if y > 0 and name != bath_gas
            ]

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def run(self) -> SimulationResult:
        """Execute the thermal-analysis simulation.

        Returns
        -------
        SimulationResult
            Dataclass holding all recorded arrays.
        """
        # --- sample gas (reactor contents) ---
        gas = ct.Solution(self.mechanism)
        gas.TPX = self.T_initial_K, self.pressure_Pa, self.sample_composition
        n_sp = gas.n_species

        # --- furnace gas (reservoir contents, fixed at T_initial) ---
        furnace_gas = ct.Solution(self.mechanism)
        furnace_gas.TP = self.T_initial_K, self.pressure_Pa

        # --- reactor network ---
        reactor = ct.ConstPressureReactor(gas)
        # Scale reactor volume so that reactor mass matches sample_mass_mg
        reactor.volume = (self.sample_mass_mg * 1e-6) / gas.density

        furnace = ct.Reservoir(furnace_gas)

        # The furnace profile is realised as a time-dependent wall heat
        # flux.  Total wall heat rate [W]:
        #   Q_wall = U·A·(T_reservoir − T_reactor) + A·Q(t)
        # With T_reservoir fixed at T_initial and
        #   Q(t) = U·(T_program(t) − T_initial)  [W/m²]:
        #   Q_wall = U·A·(T_program(t) − T_reactor)
        program = self.program
        htc = self.htc
        T_init_K = self.T_initial_K
        Q_ramp = ct.Func1(lambda t: htc * (program.T_furnace_K(t) - T_init_K))

        wall = ct.Wall(
            furnace, reactor,
            A=self.wall_area, U=self.htc, Q=Q_ramp,
        )
        net = ct.ReactorNet([reactor])

        # --- pre-allocate (trimmed later if we break early) ---
        n = self.n_steps + 1  # +1 for initial state at t=0
        times = np.empty(n)
        T_furnace = np.empty(n)
        T_reactor = np.empty(n)
        mass_reactor = np.empty(n)
        Y_hist = np.empty((n, n_sp))
        X_hist = np.empty((n, n_sp))
        q_wall = np.empty(n)
        rho_hist = np.empty(n)

        # --- record initial state (t = 0) ---
        times[0] = 0.0
        T_furnace[0] = self.T_initial_K - 273.15
        T_reactor[0] = reactor.thermo.T - 273.15
        mass_reactor[0] = reactor.mass
        Y_hist[0] = reactor.thermo.Y
        X_hist[0] = reactor.thermo.X
        q_wall[0] = wall.heat_rate
        rho_hist[0] = reactor.thermo.density

        # --- heating loop ---
        n_recorded = 1
        for step in range(1, self.n_steps + 1):
            t = step * self.dt
            if t > self.t_final_s + 1e-10:
                break

            # Advance the reactor network (Cantera handles energy +
            # chemistry via the coupled ODE system; the wall Func1
            # provides the furnace profile continuously).
            net.advance(t)

            # Record state
            times[n_recorded] = t
            T_furnace[n_recorded] = program.T_furnace_C(t)
            T_reactor[n_recorded] = reactor.thermo.T - 273.15
            mass_reactor[n_recorded] = reactor.mass
            Y_hist[n_recorded] = reactor.thermo.Y
            X_hist[n_recorded] = reactor.thermo.X
            q_wall[n_recorded] = wall.heat_rate
            rho_hist[n_recorded] = reactor.thermo.density
            n_recorded += 1

        # trim to actual length
        times = times[:n_recorded]
        T_furnace = T_furnace[:n_recorded]
        T_reactor = T_reactor[:n_recorded]
        mass_reactor = mass_reactor[:n_recorded]
        Y_hist = Y_hist[:n_recorded]
        X_hist = X_hist[:n_recorded]
        q_wall = q_wall[:n_recorded]
        rho_hist = rho_hist[:n_recorded]

        # per-species dictionaries
        mass_fracs = {
            sp: Y_hist[:, i] for i, sp in enumerate(self.species_names)
        }
        mole_fracs = {
            sp: X_hist[:, i] for i, sp in enumerate(self.species_names)
        }

        return SimulationResult(
            time_s=times,
            furnace_temperature_C=T_furnace,
            reactor_temperature_C=T_reactor,
            reactor_mass=mass_reactor,
            mass_fractions=mass_fracs,
            mole_fractions=mole_fracs,
            wall_heat_flux=q_wall,
            density=rho_hist,
            species_names=self.species_names,
            condensed_species=self.condensed_species,
            sample_mass_mg=self.sample_mass_mg,
            bath_gas=self.bath_gas,
            bath_gas_flow_sccm=self.bath_gas_flow_sccm,
            pressure_Pa=self.pressure_Pa,
            temperature_program=self.program,
        )
