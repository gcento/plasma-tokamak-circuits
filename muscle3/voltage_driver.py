"""Single-channel voltage driver for the standalone MUSCLE3 test of the linear-model actor.

Stands in for the PDS controller (actor_standalone.ymmsl): sends the F_INIT
equilibrium and pf_active, then answers every actor output with the 14 power-supply
voltages, all zero except one channel that follows a waveform. Waveform and sampling are
those of utils/supply_sequence_response.m with one supply, so the two results can be
overlaid and must coincide (utils/run_supply_sequence.m, parameter muscle3_results).

F_INIT is all zero: coil currents, ip, geometric and magnetic axis, X-point and gaps
1..3. With the actor setting operating_point: f_init every output is then a deviation
from the operating point of the model (an absolute model, e.g. plasmaless_model, starts
from zero currents), as in supply_sequence_response.m.

Ports: O_I equilibrium_out, pf_active_out (to the actor's F_INIT), voltages_out (to
pf_active_s); S equilibrium_in, pf_active_in, pf_passive_in (actor's O_I).
The F_INIT data goes out before the O_I/S loop, without voltages_out: libmuscle warns
once that the first receive does not adhere to the MMSF (harmless here).

Settings: dt, t_end (shared with the actor); md_uri, results_file, model_dir (absolute
paths set by run_actor_standalone.sh; model_dir is only recorded in the results);
channel (power supply), waveform (step | pulse | bipolar), amplitude [V],
t_on [s], width [s].
Waveform on the sample index k (t_k = k*dt, u_k applied over [t_k, t_k + dt]), with
k_on = round(t_on/dt), n_w = round(width/dt) and n_h = round(n_w/2), rounding half up:
  step     u_k = A                                 for k >= k_on
  pulse    u_k = A                                 for k_on <= k < k_on + n_w
  bipolar  u_k = +A for k_on <= k < k_on + n_h, -A for k_on + n_h <= k < k_on + n_w
The actor's first step applies the operating-point voltages (du = 0) whatever the driver
sends, so t_on must be >= dt.
Results (results_file, MATLAB .mat): time, u (channel voltage), Ia, Ie, Y (observables
read back from equilibrium_in, NaN where absent), coil_names, loop_names, y_names, and
the settings channel, waveform, amplitude, t_on, width, dt, model_dir.
"""

import math

import imas
import numpy as np
import scipy.io
from libmuscle import Instance, Message
from ymmsl import Operator

TSTART = 0.0
# power supplies in the order of pf_active_s (coil{k}.voltage = supply k), the inputs of
# every model of linear_models/
SUPPLIES = ["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L", "PF1", "PF2", "PF3", "PF4",
            "PF5", "VS1", "PF6", "VS3"]
WAVEFORMS = ("step", "pulse", "bipolar")
# observables vocabulary of plant/plant_observables_to_equilibrium.m
Y_NAMES = ["rmag", "zmag", "rX", "zX", "Ip", "gap1", "gap2", "gap3", "Dist2ndSep"]
N_GAPS = 3


def _value(node):
    """Float of a filled IDS scalar, NaN when it is not filled."""
    return float(node) if node.has_value else np.nan


def observables(ts):
    """Row of Y_NAMES values of one equilibrium time slice (DD4 fields written by
    plant_observables_to_equilibrium.m): magnetic axis, the first contour_tree
    node with critical_type 1 (X-point), ip, up to 3 gaps; Dist2ndSep has no DD4 field
    and is always NaN."""
    row = np.full(len(Y_NAMES), np.nan)
    row[0] = _value(ts.global_quantities.magnetic_axis.r)
    row[1] = _value(ts.global_quantities.magnetic_axis.z)
    for node in ts.contour_tree.node:
        if node.critical_type.has_value and int(node.critical_type) == 1:
            row[2] = _value(node.r)
            row[3] = _value(node.z)
            break
    row[4] = _value(ts.global_quantities.ip)
    for k, gap in enumerate(ts.boundary.gap):
        if k >= N_GAPS:
            break
        row[5 + k] = _value(gap.value)
    return row


def round_half_up(x):
    """MATLAB's round for x >= 0 (Python's round goes half to even)."""
    return math.floor(x + 0.5)


def channel_voltage(k, wave):
    """Voltage of the driven channel at sample k (module docstring)."""
    k_on = round_half_up(wave["t_on"] / wave["dt"])
    n_w = round_half_up(wave["width"] / wave["dt"])
    a = wave["amplitude"]
    if wave["waveform"] == "step":
        return a if k >= k_on else 0.0
    if wave["waveform"] == "pulse":
        return a if k_on <= k < k_on + n_w else 0.0
    n_h = round_half_up(n_w / 2)  # bipolar
    if k_on <= k < k_on + n_h:
        return a
    return -a if k_on + n_h <= k < k_on + n_w else 0.0


def zero_equilibrium(factory, times):
    """F_INIT equilibrium: ip, geometric and magnetic axis, X-point and gaps all 0 at the
    given times (the reference the actor passes through, and the operating point of the
    observables with operating_point: f_init)."""
    eq = factory.equilibrium()
    eq.ids_properties.homogeneous_time = 1
    eq.time = np.array(times)
    eq.time_slice.resize(len(times))
    for ts, t in zip(eq.time_slice, times):
        ts.time = t
        ts.global_quantities.ip = 0.0
        ts.global_quantities.magnetic_axis.r = 0.0
        ts.global_quantities.magnetic_axis.z = 0.0
        ts.boundary.geometric_axis.r = 0.0
        ts.boundary.geometric_axis.z = 0.0
        ts.contour_tree.node.resize(1)
        ts.contour_tree.node[0].critical_type = 1
        ts.contour_tree.node[0].r = 0.0
        ts.contour_tree.node[0].z = 0.0
        ts.boundary.gap.resize(N_GAPS)
        for g, gap in enumerate(ts.boundary.gap):
            gap.name = f"gap{g + 1}"
            gap.value = 0.0
    return eq


inst = Instance(
    {
        Operator.O_I: ["equilibrium_out", "pf_active_out", "voltages_out"],
        Operator.S: ["equilibrium_in", "pf_active_in", "pf_passive_in"],
    }
)
while inst.reuse_instance():
    dt = float(inst.get_setting("dt"))
    t_end = float(inst.get_setting("t_end"))
    md_uri = inst.get_setting("md_uri", "str")
    results_file = inst.get_setting("results_file", "str")
    model_dir = inst.get_setting("model_dir", "str", default="")
    wave = {
        "channel": inst.get_setting("channel", "str"),
        "waveform": inst.get_setting("waveform", "str"),
        "amplitude": float(inst.get_setting("amplitude")),
        "t_on": float(inst.get_setting("t_on")),
        "width": float(inst.get_setting("width")),
        "dt": dt,
    }
    if wave["channel"] not in SUPPLIES:
        raise ValueError(f"channel {wave['channel']!r} is not one of {SUPPLIES}")
    if wave["waveform"] not in WAVEFORMS:
        raise ValueError(f"waveform {wave['waveform']!r} is not one of {WAVEFORMS}")
    if wave["t_on"] < dt or wave["width"] <= 0:
        raise ValueError("need t_on >= dt (the actor's first step applies the "
                         "operating-point voltages) and width > 0")
    j = SUPPLIES.index(wave["channel"])

    with imas.DBEntry(md_uri, "r") as entry:
        pfa_md = entry.get("pf_active")
        factory = entry.factory

    # F_INIT equilibrium: gives t0 (message timestamp), the reference that the actor
    # passes through and the zero operating point of the observables. It ends at
    # t_end + dt because the actor's accumulated time can exceed t_end by round-off and
    # the actor interpolates the reference at that time.
    eq = zero_equilibrium(factory, [TSTART, t_end + dt])
    inst.send("equilibrium_out", Message(TSTART, None, eq.serialize()))

    # F_INIT pf_active: coils of the machine description (same names and order as the
    # model), currents and voltages 0
    pfa = factory.pf_active()
    pfa.ids_properties.homogeneous_time = 1
    pfa.time = np.array([TSTART])
    pfa.coil.resize(len(pfa_md.coil))
    for c, c_md in zip(pfa.coil, pfa_md.coil):
        c.name = c_md.name
        c.current.data = np.array([0.0])
        c.voltage.data = np.array([0.0])
    inst.send("pf_active_out", Message(TSTART, None, pfa.serialize()))

    time, u_hist = [TSTART], [channel_voltage(0, wave)]
    ia, ie = [np.zeros(len(pfa_md.coil))], None
    y = [observables(eq.time_slice[0])]  # F_INIT values at t0
    while True:
        msg_eq = inst.receive("equilibrium_in")
        msg_pfa = inst.receive("pf_active_in")
        msg_pfp = inst.receive("pf_passive_in")
        eq_out = factory.equilibrium()
        eq_out.deserialize(msg_eq.data)
        y.append(observables(eq_out.time_slice[0]))
        pfa_out = factory.pf_active()
        pfa_out.deserialize(msg_pfa.data)
        pfp_out = factory.pf_passive()
        pfp_out.deserialize(msg_pfp.data)
        if ie is None:  # passive currents 0 at t0 (actor assumption)
            ie = [np.zeros(len(pfp_out.loop))]
            coil_names = [str(c.name) for c in pfa_out.coil]
            loop_names = [str(lp.name) for lp in pfp_out.loop]
        t = msg_pfa.timestamp
        k = round_half_up((t - TSTART) / dt)
        v = channel_voltage(k, wave)
        time.append(t)
        u_hist.append(v)
        ia.append([float(c.current.data[0]) for c in pfa_out.coil])
        ie.append([float(lp.current[0]) for lp in pfp_out.loop])
        if msg_pfa.next_timestamp is None:
            break
        # voltages applied over [t, t+dt]: v on the channel, 0 elsewhere; no coil names
        # (as the PDS controller)
        cmd = factory.pf_active()
        cmd.ids_properties.homogeneous_time = 1
        cmd.time = np.array([t])
        cmd.coil.resize(len(SUPPLIES))
        for i, c in enumerate(cmd.coil):
            c.voltage.data = np.array([v if i == j else 0.0])
        inst.send("voltages_out", Message(t, None, cmd.serialize()))

    scipy.io.savemat(
        results_file,
        {
            "time": np.array(time),
            "u": np.array(u_hist),
            "Ia": np.array(ia),
            "Ie": np.array(ie),
            "Y": np.array(y),
            "coil_names": np.array(coil_names, dtype=object).reshape(1, -1),
            "loop_names": np.array(loop_names, dtype=object).reshape(1, -1),
            "y_names": np.array(Y_NAMES, dtype=object).reshape(1, -1),
            "channel": wave["channel"],
            "waveform": wave["waveform"],
            "amplitude": wave["amplitude"],
            "t_on": wave["t_on"],
            "width": wave["width"],
            "dt": dt,
            "model_dir": model_dir,
        },
        oned_as="column",
    )
    print(
        f"voltage_driver: {len(time) - 1} steps, t = {time[0]} .. {time[-1]}, channel "
        f"{wave['channel']} {wave['waveform']} {wave['amplitude']} V, results saved to "
        f"{results_file}"
    )
