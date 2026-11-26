import re
import time
import math
import matplotlib.pyplot as plt

# ==========================================================
# USER CONFIGURATION — CHANGE SETTINGS HERE
# ==========================================================

LOG_FILE = "log.simpleFoam"
REFRESH_INTERVAL = 0.5

# ----------------------------------------------------------
# PLOTS: True = enabled, False = disabled
# ----------------------------------------------------------

PLOT_INITIAL_RESIDUALS = True
PLOT_FINAL_RESIDUALS = False
PLOT_COEFFICIENTS = True
PLOT_CONTINUITY = True

# ----------------------------------------------------------
# RESIDUAL FIELDS: enable/disable individually
# ----------------------------------------------------------

PLOT_UX = True
PLOT_UY = True
PLOT_UZ = True
PLOT_P = True
PLOT_K = True
PLOT_OMEGA = True

# ----------------------------------------------------------
# PRESSURE RESIDUAL HANDLING
#
# "first" = first pressure solve per SIMPLE iteration
# "last"  = last pressure solve per SIMPLE iteration
# "both"  = plot first and last pressure solves
# ----------------------------------------------------------

PRESSURE_MODE = "first"

# ----------------------------------------------------------
# LIVE MONITORING
# ----------------------------------------------------------

SHOW_CURRENT_ITERATION = True
# True: show the current iteration before the next
#       "Time = ..." line appears.
# False: show completed iterations only.

# ----------------------------------------------------------
# RESIDUAL PLOT SETTINGS
# ----------------------------------------------------------

RESIDUAL_LOG_SCALE = True

# ----------------------------------------------------------
# CONTINUITY PLOT SETTINGS
#
# "symlog" = handles positive, negative and near-zero values
# "linear" = ordinary linear axis
# ----------------------------------------------------------

CONTINUITY_SCALE = "symlog"

# ==========================================================
# CONFIGURATION VALIDATION
# ==========================================================

if PRESSURE_MODE not in ("first", "last", "both"):
    raise ValueError(
        "PRESSURE_MODE must be 'first', 'last', or 'both'"
    )

if CONTINUITY_SCALE not in ("linear", "symlog"):
    raise ValueError(
        "CONTINUITY_SCALE must be 'linear' or 'symlog'"
    )

FIELD_SETTINGS = {
    "Ux": PLOT_UX,
    "Uy": PLOT_UY,
    "Uz": PLOT_UZ,
    "p": PLOT_P,
    "k": PLOT_K,
    "omega": PLOT_OMEGA,
}

FIELDS = [
    field for field, enabled in FIELD_SETTINGS.items()
    if enabled
]

PLOT_RESIDUALS = (
    PLOT_INITIAL_RESIDUALS or PLOT_FINAL_RESIDUALS
)

PLOT_ANYTHING = (
    PLOT_RESIDUALS
    or PLOT_COEFFICIENTS
    or PLOT_CONTINUITY
)

if not PLOT_ANYTHING:
    raise ValueError("Enable at least one plot.")

# ==========================================================
# REGEX PATTERNS
# ==========================================================

re_time = re.compile(
    r"^\s*Time\s*=\s*([\d.eE+-]+)"
)

re_residual = re.compile(
    r"Solving for (\w+), Initial residual = "
    r"([\d.eE+-]+), Final residual = "
    r"([\d.eE+-]+)"
)

re_cd = re.compile(r"\bCd:\s*([\d.eE+-]+)")
re_cl = re.compile(r"\bCl:\s*([\d.eE+-]+)")

re_continuity = re.compile(
    r"sum local\s*=\s*([\d.eE+-]+)"
    r"\s*,\s*global\s*=\s*([\d.eE+-]+)"
    r"\s*,\s*cumulative\s*=\s*([\d.eE+-]+)"
)

# ==========================================================
# DATA STORAGE
# ==========================================================

history = []

current_iter = None
current_data = None


def new_iteration(iteration):
    return {
        "iteration": iteration,
        "initial_first": {},
        "initial_last": {},
        "final_first": {},
        "final_last": {},
        "coefficients": {},
        "continuity": {},
    }


def finalize_iteration():
    """Save the current iteration without modifying its data."""
    if current_data is not None:
        history.append(current_data)


# ==========================================================
# LOG PARSING
# ==========================================================

def process_line(line):
    global current_iter, current_data

    # ------------------------------------------------------
    # New SIMPLE iteration
    # ------------------------------------------------------

    match = re_time.search(line)

    if match:
        finalize_iteration()

        current_iter = float(match.group(1))
        current_data = new_iteration(current_iter)
        return

    if current_data is None:
        return

    # ------------------------------------------------------
    # Solver residuals
    # ------------------------------------------------------

    match = re_residual.search(line)

    if match:
        field = match.group(1)

        if field in FIELDS:
            initial = float(match.group(2))
            final = float(match.group(3))

            data = current_data

            # Preserve the first solve and continually update
            # the last solve for each field.
            if field not in data["initial_first"]:
                data["initial_first"][field] = initial
                data["final_first"][field] = final

            data["initial_last"][field] = initial
            data["final_last"][field] = final

    # ------------------------------------------------------
    # Aerodynamic coefficients
    # ------------------------------------------------------

    match = re_cd.search(line)

    if match and "Cd" not in current_data["coefficients"]:
        current_data["coefficients"]["Cd"] = float(
            match.group(1)
        )

    match = re_cl.search(line)

    if match and "Cl" not in current_data["coefficients"]:
        current_data["coefficients"]["Cl"] = float(
            match.group(1)
        )

    # ------------------------------------------------------
    # Continuity errors
    # ------------------------------------------------------

    match = re_continuity.search(line)

    if match and not current_data["continuity"]:
        current_data["continuity"] = {
            "Local": float(match.group(1)),
            "Global": float(match.group(2)),
            "Cumulative": float(match.group(3)),
        }


# ==========================================================
# PRESSURE SERIES CONFIGURATION
# ==========================================================

def get_residual_series(residual_type):
    """
    Return the enabled residual series and their labels.
    Non-pressure fields use their first solve per iteration.
    """

    series = []

    for field in FIELDS:

        if field == "p":
            if PRESSURE_MODE in ("first", "both"):
                series.append((
                    "p (first)",
                    "initial_first" if residual_type == "initial"
                    else "final_first",
                    field,
                ))

            if PRESSURE_MODE in ("last", "both"):
                series.append((
                    "p (last)",
                    "initial_last" if residual_type == "initial"
                    else "final_last",
                    field,
                ))

        else:
            series.append((
                field,
                "initial_first" if residual_type == "initial"
                else "final_first",
                field,
            ))

    return series


# ==========================================================
# PLOT SETUP
# ==========================================================

plt.ion()

plot_types = []

if PLOT_INITIAL_RESIDUALS:
    plot_types.append("initial")

if PLOT_FINAL_RESIDUALS:
    plot_types.append("final")

if PLOT_COEFFICIENTS:
    plot_types.append("coefficients")

if PLOT_CONTINUITY:
    plot_types.append("continuity")

fig, axes_array = plt.subplots(
    len(plot_types),
    1,
    figsize=(12, 4 * len(plot_types)),
    squeeze=False,
)

axes = dict(zip(plot_types, axes_array[:, 0]))
plot_lines = {}

# ----------------------------------------------------------
# Initial and final residual plots
# ----------------------------------------------------------

for residual_type in ("initial", "final"):

    if residual_type not in axes:
        continue

    ax = axes[residual_type]
    plot_lines[residual_type] = {}

    for label, _, _ in get_residual_series(residual_type):

        if RESIDUAL_LOG_SCALE:
            line, = ax.semilogy([], [], label=label)
        else:
            line, = ax.plot([], [], label=label)

        plot_lines[residual_type][label] = line

    title = (
        "Initial Solver Residuals"
        if residual_type == "initial"
        else "Final Solver Residuals"
    )

    ax.set_title(title)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Residual")
    ax.grid(True, which="both", linestyle="--")
    ax.legend(ncol=3)

# ----------------------------------------------------------
# Coefficients
# ----------------------------------------------------------

if "coefficients" in axes:

    ax = axes["coefficients"]

    cl_line, = ax.plot([], [], label=r"$C_l$")
    cd_line, = ax.plot([], [], label=r"$C_d$")

    plot_lines["coefficients"] = {
        "Cl": cl_line,
        "Cd": cd_line,
    }

    ax.set_title("Aerodynamic Coefficients")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Coefficient")
    ax.grid(True, which="both", linestyle="--")
    ax.legend()

# ----------------------------------------------------------
# Continuity errors
# ----------------------------------------------------------

if "continuity" in axes:

    ax = axes["continuity"]

    if CONTINUITY_SCALE == "symlog":
        ax.set_yscale("symlog", linthresh=1e-10)

    plot_lines["continuity"] = {}

    for field in ("Local", "Global", "Cumulative"):
        line, = ax.plot([], [], label=field)
        plot_lines["continuity"][field] = line

    ax.set_title("Continuity Errors")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Continuity Error")
    ax.grid(True, which="both", linestyle="--")
    ax.legend()

fig.tight_layout()


# ==========================================================
# BUILD DISPLAY DATA
# ==========================================================

def get_display_history():
    """
    Include the current unfinished iteration when requested,
    without adding duplicate entries to permanent history.
    """

    data = list(history)

    if SHOW_CURRENT_ITERATION and current_data is not None:
        data.append(current_data)

    return data


def values_for(data, section, field):
    values = []

    for item in data:
        value = item.get(section, {}).get(field, float("nan"))
        values.append(value)

    return values


# ==========================================================
# UPDATE PLOTS
# ==========================================================

def update_plot():
    data = get_display_history()

    x = [item["iteration"] for item in data]

    for residual_type in ("initial", "final"):

        if residual_type not in plot_lines:
            continue

        for label, section, field in get_residual_series(
            residual_type
        ):
            y = values_for(data, section, field)

            if label.endswith("(first)"):
                present = any(
                    math.isfinite(v)
                    for item in data
                    for v in [item.get(section, {}).get(field, float("nan"))]
                )
            else:
                present = any(math.isfinite(v) for v in y)

            # Only show series that actually appear in the log.
            line = plot_lines[residual_type][label]
            line.set_visible(present)
            line.set_data(x, y)

    if "coefficients" in plot_lines:
        for field, line in plot_lines["coefficients"].items():
            line.set_data(
                x,
                values_for(data, "coefficients", field),
            )

    if "continuity" in plot_lines:
        for field, line in plot_lines["continuity"].items():
            line.set_data(
                x,
                values_for(data, "continuity", field),
            )

    for ax in axes.values():
        ax.relim()
        ax.autoscale_view()

    fig.canvas.draw_idle()
    fig.canvas.flush_events()


# ==========================================================
# LIVE MONITOR
# ==========================================================

print(f"Monitoring: {LOG_FILE}")
print(f"Refresh interval: {REFRESH_INTERVAL} s")
print(f"Pressure mode: {PRESSURE_MODE}")
print("Close the plot window to stop.")

try:
    with open(LOG_FILE, "r") as log:

        while plt.fignum_exists(fig.number):

            new_data = False

            # Read every currently available line
            while True:
                line = log.readline()

                if not line:
                    break

                process_line(line)
                new_data = True

            # Redraw when data arrives, or periodically if
            # the current unfinished iteration is displayed.
            if new_data or SHOW_CURRENT_ITERATION:
                update_plot()

            plt.pause(REFRESH_INTERVAL)

except KeyboardInterrupt:
    print("Monitoring stopped.")

finally:
    plt.ioff()

    if plt.fignum_exists(fig.number):
        plt.close(fig)
