import mimetypes
import firefly as fly

plt.rcParams["axes.prop_cycle"] = cycler(color=["#9a05fc", # Purple
                                                "#f00524", # Red
                                                "#f80af1", # Pink
                                                "#0a68f8", # Blue
                                                "#1c841f", # Green
                                                "#fc9303", # Orange
                                                "#865522", # Brown
                                                "#00c7a9",  # Teal (adds cool contrast)
                                                "#C9C22A",  # Yellow (bright mid tone)
                                                "#7f7f7f",  # Gray (neutral)
                                                "black"])
plt.rcParams["lines.linewidth"] = 1.0

plt.rcParams['lines.markersize'] = 5.0
plt.rcParams["scatter.marker"] = '.'


def get_file_type(file_path):
    mime_type, encoding = mimetypes.guess_type(file_path)
    return encoding


def main(files, plot_type, ax=None):
    if get_file_type == "h5":
        data_field = fly.Field()
    pass

