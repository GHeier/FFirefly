from .src.objects import diagram

globals().update({k: getattr(diagram, k) for k in dir(diagram) if not k.startswith("_")})
