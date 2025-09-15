from .src.config.load import config

globals().update({k: getattr(config, k) for k in dir(config) if not k.startswith("_")})
