# Making `import firefly` work

The `firefly` Python package lives at `pypkg/firefly` (declared in `pypkg/pyproject.toml`).
It isn't on `PYTHONPATH` by default — you need to install it (in editable mode, so edits to
the source are picked up immediately) into **every Python interpreter** that needs to
`import firefly`.

This machine has two independent interpreters, each with its own `site-packages`, so the
install has to be repeated for each:

| Interpreter | Path | Why it needs `firefly` |
|---|---|---|
| conda `base` env | `/home/g/miniconda3/bin/python3` (3.13) | normal interactive/script use |
| system Python | `/usr/bin/python3` (3.14) | **PyCall.jl** is built against this exact interpreter (see `~/.julia/packages/PyCall/*/deps/deps.jl` → `pyprogramname = "/usr/bin/python3"`), and the Julia package's `julia_config.jl` does `pyimport("firefly")` |

Installing into one does **not** make the other see it — `pip install -e` only registers
the package in the `site-packages` of whichever interpreter's `pip` ran it.

## Prerequisite: the `src`/`plot` symlinks must resolve

`pypkg/firefly/__init__.py` reaches the C++/Fortran-backed core via relative imports
(`from .src...`, `from . import plot`). These resolve through symlinks
`pypkg/firefly/src` and `pypkg/firefly/plot`, which must point at the repo's top-level
`src/` and `plot/` directories:

```bash
ls -la /home/g/FFirefly/pypkg/firefly/   # should show src -> ../../src, plot -> ../../plot
readlink -f /home/g/FFirefly/pypkg/firefly/src    # -> /home/g/FFirefly/src
readlink -f /home/g/FFirefly/pypkg/firefly/plot   # -> /home/g/FFirefly/plot
```

If they're dangling (e.g. point at a stale absolute path from a different checkout
location), recreate them as relative links so they survive the repo being cloned/moved
anywhere:

```bash
cd /home/g/FFirefly/pypkg/firefly
rm -f src plot
ln -s ../../src src
ln -s ../../plot plot
```

`pypkg/pyproject.toml` also needs a standard package-discovery config (the package lives
at `pypkg/firefly/`, with `pyproject.toml` at `pypkg/`):

```toml
[tool.setuptools.packages.find]
where = ["."]
include = ["firefly*"]
```

## 1. Conda interpreter

With the conda env active (so `pip` resolves to `/home/g/miniconda3/bin/pip`):

```bash
pip install -e /home/g/FFirefly/pypkg
```

Verify:

```bash
python3 -c "import firefly; print('OK', firefly.__file__)"
```

## 2. System interpreter (`/usr/bin/python3`)

The system Python ships with neither `pip` nor `ensurepip`, so install pip via apt first:

```bash
sudo apt update && sudo apt install -y python3-pip
```

Debian/Ubuntu's system Python also refuses plain `pip install` outside apt-managed paths
(PEP 668 "externally-managed-environment", to stop pip and apt from clobbering each
other's files). `--break-system-packages` is the deliberate, narrow override for cases
like this where you specifically want to manage a package outside apt:

```bash
sudo /usr/bin/python3 -m pip install --break-system-packages numpy matplotlib
sudo /usr/bin/python3 -m pip install --break-system-packages -e /home/g/FFirefly/pypkg
```

(`numpy` and `matplotlib` are the third-party packages `firefly`'s eager-import chain
needs — `cpp_imports.py` imports `numpy`, and `plot/fly_plot.py` imports `matplotlib`,
which pulls in `cycler`. Heavier optional deps like `h5py`/`triqs` are only touched by
lazily-loaded attributes such as `firefly.diagram`/`firefly.load_triqs_H` and aren't
required for the bare `import firefly`.)

Verify:

```bash
/usr/bin/python3 -c "import firefly; print('OK', firefly.__file__)"
```

## 3. Julia package (`jlpkg/Firefly`)

The Julia package depends on everything above: `jlpkg/Firefly/src/Firefly.jl` does
`include("./src/config/load/julia_config.jl")`, and that file does

```julia
using PyCall
firefly = pyimport("firefly")
cfg = firefly.config
```

So `using Firefly` only works once **(a)** the `include` paths resolve and **(b)** PyCall's
Python (`/usr/bin/python3` — see the table above) can `import firefly`, i.e. step 2 has
been done.

### Fix the package's own dangling symlink

Just like `pypkg/firefly/{src,plot}`, `jlpkg/Firefly/src/src` is a symlink that must point
at the repo's top-level `src/`:

```bash
readlink -f /home/g/FFirefly/jlpkg/Firefly/src/src   # -> /home/g/FFirefly/src
```

If it's dangling, recreate it as a relative link (3 levels up from
`jlpkg/Firefly/src/src` reaches the repo root):

```bash
cd /home/g/FFirefly/jlpkg/Firefly/src
rm -f src
ln -s ../../../src src
```

### Instantiate / precompile (per `README.md`)

```bash
julia --project=/home/g/FFirefly/jlpkg/Firefly -e 'using Pkg; Pkg.Registry.update(); Pkg.instantiate(); Pkg.precompile()'
```

### Verify

```bash
julia --project=/home/g/FFirefly/jlpkg/Firefly -e '
using Firefly
using PyCall
println(pyimport("firefly").config)   # should print the firefly.config module, not throw
Firefly.printv("Firefly loaded OK")
'
```

If `pyimport("firefly")` raises `ModuleNotFoundError`, PyCall isn't pointed at an
interpreter that has `firefly` installed — check
`~/.julia/packages/PyCall/*/deps/deps.jl` for `pyprogramname`/`python`, and either
install `firefly` into *that* interpreter (as in step 2, adjusted for its path) or
rebuild PyCall against one that already has it:

```julia
ENV["PYTHON"] = "/path/to/desired/python3"
using Pkg
Pkg.build("PyCall")
```
