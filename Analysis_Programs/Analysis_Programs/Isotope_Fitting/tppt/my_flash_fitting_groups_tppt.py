# Refactored and methodized version of the original code with modular functions and comments

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from iminuit import cost, Minuit
import sys, os, string, inspect
from scipy.ndimage import gaussian_filter1d
from collections import defaultdict
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import medfilt

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, argrelextrema

from FLASH_Fitting_Command_Line_Args import parse_arguments
from channel_indices import indices
from fit_params import Isotopes_Lifetimes_Dict
from SpillTimeFinder import SpillTime

ln2 = np.log(2)
e = np.exp(1)

# ----------------------------
# Utility Functions
# ----------------------------

def generate_text_y(num, max_val, subplots=False):
    """Generate y-axis text positions for annotations."""
    if max_val > 50 and not subplots:
        return [max_val - 0.04 * (i + 1) * max_val for i in range(num + 1)]
    elif max_val <= 50 and not subplots:
        return [max_val - 0.06 * (i + 1) * max_val for i in range(num + 1)]
    else:
        return [0.96 * max_val - 0.08 * (i + 1) * max_val for i in range(num + 1)]

def to_geo(petsys_id):
    """Convert PETSys ID to geometric ID."""
    return 8 * indices.get(petsys_id)[0] + indices.get(petsys_id)[1]

def to_geo_channel_id(abs_channel_id, geo_channels):
    """Convert PETSys absolute channel ID to geometric ID."""
    slaveID = abs_channel_id // 4096
    chipID = (abs_channel_id - slaveID * 4096) // 64
    channelID = abs_channel_id % 64
    PCB_ChanID = 64 * (chipID % 2) + channelID
    AbsPCB_ChanID = geo_channels[geo_channels[:, 0] == PCB_ChanID][0][1]
    return 10**4 * slaveID + 10**2 * chipID + AbsPCB_ChanID % 64

def constant_background(x, A):
    """Constant background function (used for pre-spill noise fit)."""
    return A + (10**-7) * x

def create_sum_function(fit_functions):
    """
    Create a dynamic sum-of-exponentials function using selected isotopes.
    Adds custom parameter signature for minimizer compatibility.
    """
    def sum_function(t, *amps):
        return sum(a * e ** (-t / Isotopes_Lifetimes_Dict[f]) for f, a in zip(fit_functions, amps))

    param_names = [inspect.Parameter("t", inspect.Parameter.POSITIONAL_OR_KEYWORD)] + [
        inspect.Parameter(p, inspect.Parameter.POSITIONAL_OR_KEYWORD)
        for p in string.ascii_uppercase[:len(fit_functions)]
    ]
    sum_function.__signature__ = inspect.Signature(param_names)
    sum_function.__name__ = "sum_of_fits"
    return sum_function

def prepare_geo_channels():
    """Generate lookup table from PETSys ID to geometric ID."""
    return np.array([[i, to_geo(i)] for i in range(128)])

def ensure_output_dir(path):
    """Ensure that output directory exists or create it."""
    if not os.path.exists(path) and path != "":
        print("Creating output directory...")
        os.makedirs(path)

# ----------------------------
# Valley Manipulation Functions
# ----------------------------

def fill_valleys(signal, kernel_size=41, threshold=0.6):
    """
    Fill valleys in a 1D signal via log-space interpolation.

    Parameters:
        signal (np.ndarray): 1D signal array
        kernel_size (int): Size of median filter (must be odd)
        threshold (float): Fraction below local median to mark valley

    Returns:
        np.ndarray: Modified signal with valleys interpolated
    """
    if kernel_size % 2 == 0:
        raise ValueError("kernel_size must be odd")

    signal = np.asarray(signal)
    local_median = medfilt(signal, kernel_size=kernel_size)
    is_valley = signal < (local_median * threshold)

    filled_signal = signal.copy()
    valley_indices = np.where(is_valley)[0]

    if len(valley_indices) == 0:
        return filled_signal

    # Group adjacent valley indices into continuous segments
    groups = np.split(valley_indices, np.where(np.diff(valley_indices) > 1)[0] + 1)

    for group in groups:
        start = max(group[0] - kernel_size // 2, 0)
        end = min(group[-1] + kernel_size // 2 + 1, len(signal))
        window_x = np.arange(start, end)
        window_y = signal[start:end]
        window_mask = ~is_valley[start:end]  # non-valley points

        if np.sum(window_mask) < 2:
            continue

        # Interpolate using cubic spline in log-space
        log_y = np.log(window_y[window_mask])
        interp_fn = interp1d(window_x[window_mask], log_y, kind='cubic', fill_value="extrapolate")

        # Replace valley with exponentiated interpolated values
        filled_signal[group] = e ** (interp_fn(group))

    return filled_signal

def delete_valleys(values, bin_centers, kernel_size=41, threshold=0.6):
    """
    Deletes valleys from values based on a local median filter, and shifts
    the next-right bin into the deleted valley's position.

    Args:
        values (array-like): Input values (e.g., histogram heights).
        bin_centers (array-like): Bin center positions.
        kernel_size (int): Median filter size (must be odd).
        threshold (float): Valley threshold relative to local median.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Modified values and bin_centers
    """
    if kernel_size % 2 == 0:
        raise ValueError("kernel_size must be odd")

    values = np.asarray(values).copy()
    bin_centers = np.asarray(bin_centers).copy()

    local_median = medfilt(values, kernel_size=kernel_size)
    is_valley = values < (local_median * threshold)

    # Make a copy to modify
    values_mod = values.copy()
    bin_centers_mod = bin_centers.copy()
    bin_center = []

    j = 0
    for i in range(len(values) ):  # last bin can't be shifted
        if not is_valley[i]:
            bin_center.append(bin_centers[j])
            j += 1

    return values_mod[~is_valley], np.array(bin_center)

def delete_valley_points(values, bin_centers, kernel_size=41, threshold=0.6):
    """
    Deletes values and bin centers associated with detected valleys.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Cleaned values and bin centers
    """
    if kernel_size % 2 == 0:
        raise ValueError("kernel_size must be odd")

    values = np.asarray(values)
    bin_centers = np.asarray(bin_centers)

    local_median = medfilt(values, kernel_size=kernel_size)
    is_valley = values < (local_median * threshold)

    return values[~is_valley], bin_centers[~is_valley]

def stalin_sort(values, bin_centers):
    """
    Stalin sort: remove any value that is not strictly lower than the previous.

    Parameters:
        values (np.ndarray): 1D histogram values
        bin_centers (np.ndarray): Corresponding bin centers

    Returns:
        Tuple[np.ndarray, np.ndarray]: Stalin-sorted values and bin centers
    """
    sorted_vals = []
    sorted_bins = []
    last_val = float('inf')
    i=0
    for v, b in zip(values, bin_centers):
        if v < last_val:
            sorted_vals.append(v)
            sorted_bins.append(b)
            last_val = v

    return np.array(sorted_vals), np.array(sorted_bins)

# ----------------------------
# Model Fit and Plot Functions
# ----------------------------

def create_and_fit_model(args, bin_centers, values, fit_functions):
    """
    Create and fit sum-of-exponentials model to histogram data using iminuit.
    """
    # Select or generate the fit function
    fit_function = create_sum_function(fit_functions) if args.fit_function == "" else globals()[args.fit_function]

    # Create a least-squares cost function
    print(len(values))
    values = values 
    c = cost.LeastSquares(bin_centers, values, np.sqrt(values), fit_function)

    # Prepare initial parameter estimates
    if args.initial_fit_params != [""]:
        initial_fit_params = args.initial_fit_params
    else:
        # Distribute initial amplitude guesses evenly
        initial_fit_params = [
            f"{string.ascii_uppercase[i]}={max(values) / len(args.fit_isotopes)}"
            for i in range(len(args.fit_isotopes))
        ]

    # Convert list of param strings to keyword arguments for Minuit
    fit_kwargs = ",".join(initial_fit_params)
    fitter = eval(f"Minuit(c, {fit_kwargs})")

    # Enforce positive-only parameter constraints
    fitter.limits = [(0, None)] * len(fitter.values)

    # Run optimization
    fitter.migrad()
    fitter.hesse()

    # Calculate degrees of freedom
    dof = len(bin_centers) - len(fitter.values)
    return fitter, fit_function, dof

import os

def handle_batched_fitting_and_plotting(args, dfs, bins_lengths):
    group_index = 0
    batch_fig = None
    axs = None
    fitters, group_indices, dofs = [], [], []
    Depths1 = ['0mm - 12.6mm', '12.8mm - 25.4mm', '26.2mm - 38.9mm', '39.1mm - 51.7mm',]
    Depths2 = ['52.5mm - 65.1mm', '65.3mm - 78.0mm', '78.9mm - 91.5mm', '91.7mm - 104.3mm']
    for i in range(0, len(dfs)):
       
        group = pd.concat(dfs[i:i+1], ignore_index=True)
        group_index += 1
        print(group_index)
        if i < 4:
            depth = Depths2[(group_index ) % 4 * -1]
            chosen = Depths2
            print((group_index ) % 4 * -1)
        else:  
            depth = Depths1[(group_index ) % 4 * -1]
            chosen = Depths1
        print(depth)

        print(f"Group {group_index}, size: {len(group)}")

        t = group["TimeL"][group["TimeL"] >= 1] - min(group["TimeL"][group["TimeL"] >= 1])
        y_shift = 50

        for b in bins_lengths:
            binwidth = args.run_length / args.num_bins
            new_num_bins = int((max(t) - min(t)) / b)

            values, bins = np.histogram(t, bins=new_num_bins)
            bin_centers = 0.5 * (bins[1:] + bins[:-1])

            interpolated_values = fill_valleys(values, kernel_size=19, threshold=0.9)
            striped_values, striped_bins = delete_valley_points(values, bin_centers, kernel_size=19, threshold=.9)
            stalin_values, stalin_bins = stalin_sort(striped_values, bin_centers)

            fit_functions = [f for f in args.fit_isotopes if f in Isotopes_Lifetimes_Dict]

            pos_bins = [striped_bins]
            pos_values = [ striped_values]
            affidatove = ['regular', 'stripped']

            for index in range(len(pos_bins)):
                bin_centers = pos_bins[index]
                values = pos_values[index]
                fitter, fit_function, dof = create_and_fit_model(args, bin_centers, values, fit_functions)

                # Prepare 2x2 subplot figure every 4 groups
                if group_index % 4 == 1:
                    batch_fig, axs = plt.subplots(2, 2, figsize=(7,5))
                    axs = axs.flatten()

                ax = axs[(group_index) % 4 * -1]
                plot_histogram_without_residuals(ax, fitter, fit_function, args, bin_centers, values, y_shift,depth)

                fitters.append(fitter)
                group_indices.append(group_index)
                dofs.append(dof)

                # If 4 plots accumulated or last group, save batch figure and results
                if group_index % 4 == 0 or i == len(dfs) - 1:
                    fig_index = (group_index - 1) // 4 + 1
                    plot_path = os.path.join(args.output_dir, f"{args.title}_groups_{fig_index}_no_log.png")
                    plot_path_pdf = os.path.join(args.output_dir, f"{args.title}_groups_{fig_index}_no_log.pdf")
                    batch_fig.tight_layout()
                    batch_fig.savefig(plot_path)
                    batch_fig.savefig(plot_path_pdf)
                    plt.close(batch_fig)

                    latex_path = os.path.join(args.output_dir, f"{args.title}_groups_{fig_index}_results.tex")
                    export_fit_results_latex_table(fitters, args, latex_path, group_indices, dofs,Depths=chosen)

                    # Reset for next batch
                    fitters, group_indices, dofs = [], [], []



import matplotlib.pyplot as plt
import numpy as np

def plot_histogram_without_residuals(ax, fitter, fit_function, args, bin_centers, values, y_shift,depth):
    """
    Plot histogram with sum-of-exponentials fit and isotope components on a given axis.
    """
    ax.bar(
        bin_centers,
        values + y_shift,
        width=[5 for _ in range(len(bin_centers))],
        edgecolor='none',
        facecolor='blue',
        align='center',
        alpha=0.2
    )

    plot_values = np.linspace(min(bin_centers), max(bin_centers), num=25 * len(bin_centers))
    ax.set_ylim(0, 10^5)
    ax.set_yscale('log')
    flag = False
    total = sum(fitter.values)
    for i, isotope in enumerate(args.fit_isotopes):
        if isotope in Isotopes_Lifetimes_Dict and not math.isnan(fitter.errors[i]):
            decay_curve = fitter.values[i] * np.exp(-plot_values / Isotopes_Lifetimes_Dict[isotope])
            ax.plot(
                plot_values,
                decay_curve + y_shift,
                linewidth=1,
                linestyle='dotted' if i % 2 == 0 else 'dashdot',
                label=f"{isotope} ({fitter.values[i]:.1f} | {fitter.values[i]/total:.2f})"
            )
            print(f"{isotope} = {fitter.values[i]:.1f} ± {fitter.errors[i]:.1f}")
            print(type(fitter.errors[i]))
            flag = True

    sum_fit_line = fit_function(plot_values, *fitter.values)
    existing_sum_line = next((line for line in ax.get_lines() if line.get_label() == "Sum"), None)

    # Check if a plot with the label "Sum" already exists
    existing_sum_line = next((line for line in ax.get_lines() if line.get_label() == "Sum"), None)

    if existing_sum_line:
        # Remove the existing "Sum" line
        existing_sum_line.remove()

# Plot a new "Sum" line
    if flag:
        sum_fit_line = fit_function(plot_values, *fitter.values)
        ax.plot(plot_values, sum_fit_line + y_shift, linewidth=1, color="r", alpha=0.7, label="Sum")
    ax.set_title(f'{depth} Fit (reduced $\chi^2$ = {fitter.fval / (len(bin_centers) - len(fitter.values)):.2f})', fontsize=30/ 4)
    ax.legend(fontsize=6.5,ncol=2, loc='upper right')
    ax.set_xlabel('Time [s]',fontsize=42/ 4)
    ax.set_ylabel(r'PET Event Rate [s$^{-1}$]',fontsize=42/ 4)
    ax.tick_params(axis='both', which='major', labelsize=40/ 4)
    return ax



def annotate_plot(ax, fitter, args, dof):
    """
    Annotate the plot with fitted parameters, isotope half-lives, and reduced chi-squared value.
    """
    try:
        redChiSq = np.round(fitter.fval / dof,3)  # Reduced chi-squared calculation
    except ZeroDivisionError:
        redChiSq = 'inf'
    text_y_offset = 0.12
    y_start = 0.95

    for i, isotope in enumerate(args.fit_isotopes):
        y_pos = y_start - i * text_y_offset
        print(f"{chr(65 + i)} = {fitter.values[i]:.1f} ± {fitter.errors[i]:.1f}")

        if isotope in Isotopes_Lifetimes_Dict:
            half_life = Isotopes_Lifetimes_Dict[isotope] * ln2
            print(fr"     T$_{{1/2}}$ (assumed) = {half_life:.3f} s"),

    # Add reduced chi-squared
    ax.text(0.05,.01, fr"reduced $\chi^{{2}}$ = {redChiSq}", fontsize=20, transform=ax.transAxes)

def save_plot(fig, bins, args, index, title):
    """
    Save the figure in both PNG and PDF format with a consistent naming scheme.
    """
    fit_name = "_".join(args.fit_isotopes)
    if args.fit_function != "":
        fit_name = args.fit_function

    config_filename = args.run_config.replace(" ", "_")
    out_prefix = f'{args.output_dir}/PES_Activity_Fit_{title}_{index}_no_log'

    fig.savefig(out_prefix + ".png")
    fig.savefig(out_prefix + ".pdf")
    plt.close(fig)

    print(f"Saved plot to {out_prefix}.png/.pdf")


def load_prepare_and_group_by_coordinate_chunked_grouping(args,in_dir,data_file, geo_channels, coord_csv_path, coord_axis, group_size=8, chunksize=100000):
    """
    Load PET data in chunks, preprocess and group by spatial coordinates for memory-efficient handling.

    Parameters:
        args: Command-line arguments
        geo_channels: Array mapping channels
        coord_csv_path: Path to geometry CSV
        coord_axis: Axis to group by ('x', 'y', 'z')
        group_size: How many coordinates per group
        chunksize: Rows per chunk for read_csv

    Returns:
        list: List of grouped DataFrames
    """
    def to_reco_channel_id(AbsChannelID):
        portID = AbsChannelID // 131072
        slaveID = (AbsChannelID - 131072 * portID) // 4096
        chipID = (AbsChannelID - 4096 * slaveID - 131072 * portID) // 64
        channelID = AbsChannelID % 64
        return 3072 * portID + 1024 * slaveID + 64 * chipID + channelID

    # Read coordinate CSV
    mapping_df = pd.read_csv(
        coord_csv_path,
        usecols=[0, 1, 2],
        header=None,
        names=["position x", "position y", "position z"]
    )

    # Validate axis
    coord_axis = coord_axis.lower()
    if coord_axis not in ['x', 'y', 'z']:
        raise ValueError("coord_axis must be one of 'x', 'y', or 'z'.")

    # Precompute sorted coordinate groupings
    full_coords = mapping_df[f"position {coord_axis}"].dropna().unique()
    unique_coords_desc = sorted(full_coords, reverse=True)
    unique_coords_desc = unique_coords_desc[1:]
    coord_groups = [
        unique_coords_desc[i * group_size:(i + 1) * group_size]
        for i in range((len(unique_coords_desc) + group_size - 1) // group_size)
    ]

    grouped_dict = defaultdict(list)  # Temporary storage for grouped chunks
    dfs = pd.DataFrame()
    # print(coord_groups)
    # print(group_size)
    # print(len(coord_groups))

    # Read the PET dataset in chunks
    for chunk in pd.read_csv(os.path.join(in_dir, data_file), delimiter="\t", chunksize=chunksize):
        chunk.columns = ["TimeL", "ChargeL", "ChannelIDL", "TimeR", "ChargeR", "ChannelIDR"]

        # Time normalization
        chunk["TimeL"] /= 1e12
        chunk["TimeR"] /= 1e12

        # Spill time subtraction
        spill_time_end = args.spill_time_end
        if spill_time_end < 0.0:
            spill_time_end = SpillTime(chunk["TimeL"], args.spill_time_finder_window)
            if not args.dont_write and args.file:
                with open(args.file[0], 'a') as file:
                    file.write(f"\n--spill_time_end {spill_time_end}")
        
        chunk["TimeL"] -= spill_time_end

        chunk = chunk[chunk["TimeL"] > args.offset]

        # Compute reco channel IDs
        chunk["RecoChannelIDL"] = chunk["ChannelIDL"].apply(to_reco_channel_id)
        chunk["RecoChannelIDR"] = chunk["ChannelIDR"].apply(to_reco_channel_id)

        # Map to coordinates
        chunk["_coordL"] = chunk["RecoChannelIDL"].map(mapping_df[f"position {coord_axis}"])
        chunk["_coordR"] = chunk["RecoChannelIDR"].map(mapping_df[f"position {coord_axis}"])
        chunk = chunk.dropna(subset=["_coordL", "_coordR"])

        if chunk.empty:
            continue

        # Assign chunk to appropriate coordinate groups
        if group_size == -1:
            dfs = pd.concat([dfs, chunk], ignore_index=True)
        else:
            for group_coords in coord_groups:
                group_mask = chunk["_coordL"].isin(group_coords) & chunk["_coordR"].isin(group_coords)
                group = chunk[group_mask]
                if not group.empty:
                    group_label = round(np.mean(group_coords), 3)
                    grouped_dict[group_label].append(group)

    # Final assembly of grouped data
    if group_size == -1:
        return [dfs]
    final_grouped_dfs = [pd.concat(grouped_dict[label], ignore_index=True) for label in sorted(grouped_dict)]
    return final_grouped_dfs
def export_fit_results_latex_table(fitters, args, output_path, group_indices, dofs, Depths):
    """
    Export a LaTeX-formatted table for multiple fit groups.

    Parameters:
        fitters (list): List of Minuit fitted objects
        args: Parsed command-line arguments (must contain 'title' and 'fit_isotopes')
        output_path (str): Full file path (with .tex) to save output
        group_indices (list): List of group indices for labeling
        dofs (list): List of degrees of freedom for each fit
    """
    # Assume args.fit_isotopes contains formatted names like "C10 (19.3 s)"
    isotope_labels = args.fit_isotopes
    n_isotopes = len(isotope_labels)

    header_cells = ["\\textbf{Depth (mm)}"]
    header_cells += [f"\\makecell{{\\textbf{{{iso.split()[0]}}} \\\\ ({' '.join(iso.split()[1:])})}}" for iso in isotope_labels]
    header_cells += ["\\makecell{\\textbf{Reduced} \\\\ $\\chi^2$}"]

    header_row = " & ".join(header_cells) + " \\\\"

    table = "\\begin{tabular}{|" + "l" + "|c" * n_isotopes + "|c|}\n"
    table += "\\hline\n"
    table += header_row + "\n"
    table += "\\hline\n"
    Depths = Depths[::-1]
    for idx, (fitter, dof) in enumerate(zip(fitters, dofs)):
        values = fitter.values
        errors = fitter.errors
        total = sum(values)
        depth = Depths[idx]

        try:
            redChiSq = np.round(fitter.fval / dof, 3)
        except ZeroDivisionError:
            redChiSq = 'inf'

        cell_contents = [
            f"\\makecell{{{val:.2f} $\\pm$ {err:.2f} \\\\ ({(val/total)*100:.2f}\\%)}}"
            for val, err in zip(values, errors)
        ]
        row = f"{depth} & " + " & ".join(cell_contents) + f" & {redChiSq} \\\\"
        table += row + "\n"
        table += "\\hline\n"

    table += "\\end{tabular}\n"

    with open(output_path, "w") as f:
        f.write(table)

    return table  # return the string for display/debugging if needed




def main():
    """
    Updated main execution function:
    - Parses command-line arguments
    - Loads and groups data
    - Generates histograms with multiple bin widths
    - Applies optional valley deletion and Stalin sort
    - Fits exponential decay models to histogram data
    - Plots results in 2x2 subplots per figure if grouped
    """
    args = parse_arguments()
    geo_channels = prepare_geo_channels()
    imgs_dir = os.path.join(args.output_dir, "")
    ensure_output_dir(imgs_dir)

    print("Loading and preparing data...")
    dfs = load_prepare_and_group_by_coordinate_chunked_grouping(
        args, args.in_dir, args.data_file,
        geo_channels, args.geometry_path,
        args.coord_axis, group_size=args.group_size
    )

    bins_lengths = [5]

    if args.group_size != -1 and len(dfs) > 1:
        # Grouped plotting: 2x2 subplots, 4 per figure
        handle_batched_fitting_and_plotting(args, dfs, bins_lengths)

       

if __name__ == "__main__":
    main()

