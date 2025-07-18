from plot import plot_path
from plot import basic
from plot import colorplot
from plot import gap_function

import argparse
import sys
import matplotlib
import matplotlib.pyplot as plt
from cycler import cycler

#matplotlib.use('module://matplotlib-backend-kitty')  # Kitty terminal

# Joey chose colors
plt.rcParams["axes.prop_cycle"] = cycler(color=["#9a05fc", # Purple
                                                "#f00524", # Red
                                                "#f80af1", # Pink
                                                "#0a68f8", # Blue
                                                "#1c841f", # Green
                                                "#865522", # Brown
                                                "#fc9303", # Orange
                                                "black"])
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
        "-l", "--line", action="store_true", help="Perform line plot"
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
        "flags": {"verbose": args.verbose, "output": args.output, "scatter": args.scatter, "line": args.line, "Gap": args.Gap},
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
    if isinstance(files, str):
        files = [files]
    # Select plot type
    if plot_type == 'line':
        fig, ax = basic.plot_basic(files, **kwargs)
    elif plot_type == 'scatter':
        fig, ax = basic.plot_scatter(files, **kwargs)
    elif plot_type == 'bar':
        fig, ax = basic.plot_bar(files, **kwargs)
    elif plot_type == 'hist':
        fig, ax = basic.plot_hist(files, **kwargs)
    elif plot_type == 'colorgrid':
        fig, ax = colorplot.plot_colorgrid(files, **kwargs)
    elif plot_type == 'colorline':
        fig, ax = colorplot.plot_colorline(files, **kwargs)
    elif plot_type == "path":
        fig, ax = plot_path.plot_path(files, **kwargs)
    elif plot_type == "band":
        pass
        #fig, ax = plot_path.plot_band(files, **kwargs)
    else:
        raise ValueError(f"Unsupported plot type: {plot_type}")

    # Return the figure and axis for further modifications
    return fig, ax


if __name__ == "__main__":
    try:
        result = parse_arguments()
        if (not any(result['flags'].values())):
            get_flags_from_files(result)
        print(result)
        #print("Flags:", result['flags'])
        #print("Files:", result['files'])
        plot_type = ''
        if result["flags"]["line"]:
            plot_type = 'line'
        elif "band" in result["flags"]:
            plot_type = 'band'
        elif result["flags"]["scatter"]:
            plot_type = 'scatter'
        elif result["flags"]["Gap"]:
            pass
        print("plot_type: ", plot_type)
        fig, ax = sketch(result["files"], plot_type)
        plt.show()

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
