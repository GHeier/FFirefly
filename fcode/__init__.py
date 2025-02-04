from .fmodule import *

from .config.load import config


file = config.load_config()
load_c_config(file)

