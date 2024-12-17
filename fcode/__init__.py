from .fmodule import *
from .config.load import config

print("calculation: ", config.calculation)
config.load_config()
print("calculation: ", config.calculation)
