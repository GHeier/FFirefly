from . import plot_path
from . import basic
from . import gap_function

import argparse
import sys
import matplotlib
import matplotlib.pyplot as plt
from cycler import cycler

#matplotlib.use('module://matplotlib-backend-kitty')  # Kitty terminal

plt.rcParams["axes.prop_cycle"] = cycler(color=["#9a05fc", "#f00524", "#f0d802"])
plt.rcParams["lines.linewidth"] = 1.0

plt.rcParams['lines.markersize'] = 5.0
plt.rcParams["scatter.marker"] = '.'

# plt.rcParams.update(
#    {
#        "figure.facecolor": "black",
#        "axes.facecolor": "black",
#        "savefig.facecolor": "black",
#        "text.color": "white",
#        "axes.labelcolor": "white",
#        "xtick.color": "white",
#        "ytick.color": "white",
#        "axes.edgecolor": "white",
#        "grid.color": "gray",
#    }
# )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process command-line arguments: flags and files."
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose mode"
    )
    parser.add_argument(
        "--output", type=str, help="Specify output file (use --output=filename)"
    )
    parser.add_argument(
        "-s", "--scatter", action="store_true", help="Perform scatter plot"
    )
    parser.add_argument(
        "-G", "--Gap", action="store_true", help="Plot gap over mesh"
    )
    parser.add_argument(
        "-Gs", "--Gap_Surface", action="store_true", help="Plot gap over surface"
    )

    parser.add_argument("files", nargs="*", help="Input file(s)")

    args, unknown = parser.parse_known_args()
    if unknown:
        parser.error(f"Unknown arguments: {' '.join(unknown)}")

    valid_files = []
    for f in args.files:
        if f.startswith("-"):
            parser.error(f"Invalid file name or unknown flag: {f}")
        valid_files.append(f)

    return {
        "flags": {"verbose": args.verbose, "output": args.output, "scatter": args.scatter, "Gap": args.Gap},
        "files": valid_files,
    }


def get_flags_from_files(result):
    flags = result["flags"]
    for filename in result["files"]:
        lower_name = filename.lower()
        if "band" in lower_name:
            flags["band"] = True
        if "fermi_surface" in lower_name or "fs" in lower_name or "FS" in lower_name:
            flags["fermi_surface"] = True

def sketch(files, plot_type='line', **kwargs):

    # Select plot type
    if plot_type == 'line':
        fig, ax = basic.plot_basic(files, **kwargs)
    elif plot_type == 'scatter':
        fig, ax = basic.plot_scatter(files, **kwargs)
    elif plot_type == 'bar':
        fig, ax = basic.plot_bar(files, **kwargs)
    elif plot_type == 'hist':
        fig, ax = basic.plot_hist(files, **kwargs)
    elif plot_type == "path":
        fig, ax = plot_path.plot_path(files, **kwargs)
    else:
        raise ValueError(f"Unsupported plot type: {plot_type}")

    # Return the figure and axis for further modifications
    return fig, ax


if __name__ == "__main__":
    try:
        result = parse_arguments()
        get_flags_from_files(result)
        #print("Flags:", result['flags'])
        #print("Files:", result['files'])
        plot_type = ''
        if "band" in result["flags"]:
            plot_type = 'band'
        elif result["flags"]["scatter"]:
            plot_type = 'scatter'
        elif result["flags"]["Gap"]:
            pass
        fig, ax = sketch(result["files"], plot_type)
        plt.show()

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
