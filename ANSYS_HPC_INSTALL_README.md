# Ansys Electronics Desktop / HFSS on the Caltech HPC cluster

This documents the installation of Ansys Electronics Desktop 2025 R2 (HFSS)
for the Golwala group on the Caltech Resnick HPC cluster. It assumes that the
installer media has been copied to the shared group directory and that the
user has access to the Astro license server.

## 1. Installer media and installation directory

Keep the installer media separate from the installed software. In this setup:

```text
Installer media: /resnick/groups/golwala/software/ansys_em
Install target:  /resnick/groups/golwala/rshenoy/software/ansys_em
```

The installer media contains `INSTALL`, `252-1.dvd`, and the
`electromagneticssuite/` payload. The install target will contain `v252/`,
`shared_files/`, and `install.log` after installation.

Check that the media is present:

```bash
ls -l /resnick/groups/golwala/software/ansys_em/INSTALL
ls -l /resnick/groups/golwala/software/ansys_em/252-1.dvd
```

Create the destination and run the installer:

```bash
mkdir -p /resnick/groups/golwala/rshenoy/software/ansys_em
cd /resnick/groups/golwala/software/ansys_em
./INSTALL
```

When prompted for the installation directory, enter:

```text
/resnick/groups/golwala/rshenoy/software/ansys_em
```

The installation may take several minutes and approximately 40 GB in this
configuration. Confirm completion with:

```bash
tail -40 /resnick/groups/golwala/rshenoy/software/ansys_em/install.log
```

The log should end with `Installation Complete`.

## 2. Shell environment

Create `/resnick/groups/golwala/rshenoy/ansys_env.sh` with:

```bash
# Astro license server
export ANSYSLMD_LICENSE_FILE=1055@license.astro.caltech.edu
export ANSYSL_LOCK_TYPE=1

# Ansys Electronics Desktop 2025 R2
export ANSYS_ROOT=/resnick/groups/golwala/rshenoy/software/ansys_em/v252
export ANSYSEM_ROOT252="$ANSYS_ROOT/AnsysEM"
export PATH="$ANSYSEM_ROOT252:$PATH"
```

Load it in each shell that will use Ansys:

```bash
source /resnick/groups/golwala/rshenoy/ansys_env.sh
```

Verify the launcher path:

```bash
echo "$ANSYSLMD_LICENSE_FILE"
echo "$ANSYS_ROOT"
which ansysedt
ls -l "$(which ansysedt)"
```

The expected launcher is:

```text
/resnick/groups/golwala/rshenoy/software/ansys_em/v252/AnsysEM/ansysedt
```

It is normally a symlink to `.runtimeexewrapper`; that is expected.

## 3. License diagnostics

The Ansys licensing utility is installed at:

```bash
$ANSYS_ROOT/licensingclient/linx64/ansysli_util
```

Useful checks are:

```bash
ANSYS_LI="$ANSYS_ROOT/licensingclient/linx64/ansysli_util"
"$ANSYS_LI" -version
"$ANSYS_LI" -envvar
"$ANSYS_LI" -printlicpath
"$ANSYS_LI" -hostinfo license.astro.caltech.edu
```

`-printavail` may report `No Feature List Retrieved` with the Astro licensing
backend even when the installation and license path are correct. Do not use
that message alone as proof that the license is broken.

## 4. Login-node display limitation

Running:

```bash
ansysedt -version
```

may produce:

```text
Cannot open X display "(not specified)".
```

This means the login node has no graphical X session. It does not by itself
indicate a failed installation. Use a graphical interactive HPC session for
the AEDT GUI, or use PyAEDT in non-graphical mode for batch HFSS automation.

## 5. Python environment (one-time setup)

The HPC checkout is `/home/rshenoy/BBRSim/Blackbody-Simulations`.
The Python environment lives outside the checkout:

```text
/home/rshenoy/BBRSim/
├── Blackbody-Simulations/    # Git checkout
├── envs/bb/                 # Generated Python virtual environment
└── bb_env.sh                # Combined shell setup
```

This setup uses Python's `venv`, not Conda. The system Python on login3 was
3.9.21, whereas the BBRsim analysis package requires Python >=3.10. We used
the available Python 3.11.6 module. Run the following once on HPC:

```bash
cd /home/rshenoy/BBRSim
module load python/3.11.6-gcc-13.2.0-fh6i4o3
python3 -m venv /home/rshenoy/BBRSim/envs/bb
source /home/rshenoy/BBRSim/envs/bb/bin/activate

python -m pip install --upgrade pip
python -m pip install numpy scipy pandas matplotlib uproot pyaedt seaborn ipykernel
python -m pip check
```

The package selection covers both repositories:

| Packages | Purpose |
|---|---|
| numpy, scipy, pandas, matplotlib | Shared numerical analysis and plotting |
| pyaedt | Control the separately installed AEDT/HFSS application |
| seaborn | Blackbody-Simulations plotting |
| uproot | Read BBRsim ROOT output without installing ROOT's Python bindings |
| ipykernel | Optional notebook kernel support |

Dependencies are installed by pip; these commands do not pin versions. The
BBRSimulation `bbrsim` Python package, the Geant4 modules and the BBRsim build
are covered in section 8.

## 6. Combined setup for each session

Create `/home/rshenoy/BBRSim/bb_env.sh` with this content:

```bash
# bb_env.sh: Python, Ansys, Geant4 and BBRsim on the Caltech HPC.
# Use: source /home/rshenoy/BBRSim/bb_env.sh   (source it; do not run it)

# Drop MATLAB's library directory and empty entries from LD_LIBRARY_PATH:
# they break cmake and can shadow other libraries (ANSYS_HPC_INSTALL_README.md, section 8).
LD_LIBRARY_PATH="$(printf '%s' "${LD_LIBRARY_PATH:-}" | tr ':' '\n' | grep -v -e '^$' -e '/Matlab/' | paste -sd: -)"
if [ -n "$LD_LIBRARY_PATH" ]; then export LD_LIBRARY_PATH; else unset LD_LIBRARY_PATH; fi

# Load the Python installation used to create the environment.
module load python/3.11.6-gcc-13.2.0-fh6i4o3 || return

# Activate Python packages.
source /home/rshenoy/BBRSim/envs/bb/bin/activate || return

# Load Ansys paths and licensing.
source /resnick/groups/golwala/rshenoy/ansys_env.sh || return
export ANSYSEM_ROOT252="$ANSYS_ROOT/AnsysEM"

# Geant4 11.1.2 and gcc 13.2.0, the compiler of both Geant4 and the Python module.
module load gcc/13.2.0-gcc-13.2.0-w55nxkl geant4/11.1.2-gcc-13.2.0-gk5zjkn || return

# The installed BBRsim (BBRSIMDATA, PYTHONPATH, its library), once it is built.
bbrsim_env=/home/rshenoy/BBRSim/BBRSimulation/install/share/BBRsim/bbrsim_env.sh
if [ -r "$bbrsim_env" ]; then source "$bbrsim_env" || return; fi
unset bbrsim_env
```

In each new Bash session, or batch script after the module command is
available, run:

```bash
source /home/rshenoy/BBRSim/bb_env.sh
```

The module supplies the Python installation used by the virtual environment;
activation selects its installed packages. Sourcing `ansys_env.sh` exports
the Ansys paths and licensing settings to Python and its child processes.
No package reinstallation is needed for a new session.

One script serves both the HFSS and the Geant4 work. The cluster's modules do
not change `LD_LIBRARY_PATH` (checked 2026-09-30), so loading Geant4 does not
put its libraries in front of Ansys's. Section 8 covers the Geant4 part and its
verification.

There are two manually maintained scripts (`ansys_env.sh` and `bb_env.sh`).
The third file, `envs/bb/bin/activate`, is generated by `venv`; do not edit it.
Use `source`, rather than executing these setup scripts in a child shell.

The Ansys setup file is in `/resnick/groups/golwala/rshenoy/`, **not**
`/home/rshenoy/`. `ANSYS_ROOT` ends in `/v252`, while `ANSYSEM_ROOT252` must
include `/AnsysEM`, the directory containing `ansysedt`. `PATH` contains
executable directories, not the path to `ansys_env.sh`.

## 7. Verification and current status (2026-09-28)

After sourcing the combined setup, run:

```bash
python -c "import sys, ansys.aedt.core, uproot; print(sys.executable); print('Imports succeeded')"
command -v ansysedt
echo "$ANSYSEM_ROOT252"
```

The following output was confirmed on HPC:

```text
/home/rshenoy/BBRSim/envs/bb/bin/python
Imports succeeded
/resnick/groups/golwala/rshenoy/software/ansys_em/v252/AnsysEM/ansysedt
/resnick/groups/golwala/rshenoy/software/ansys_em/v252/AnsysEM
```

The installer reported completion, first-time AEDT configuration succeeded,
and the Python import/path checks above passed.

Update 2026-09-28: non-graphical AEDT 2025.2 startup on a compute node was
verified with PyAEDT 1.7.0 (`python hpc/aedt_smoke_test.py`, Slurm job
3614903, debug QOS): the desktop started over gRPC in 43 s, created and
saved a project, and released cleanly, so the base AEDT license checkout
works from compute nodes. An HFSS solver license checkout and an actual
simulation have **not yet been verified**; see `hpc/RESULTS.md`.

Before batch execution, Blackbody-Simulations still needs its graphical
`Hfss(..., non_graphical=False)` launches adapted for non-graphical use, plus
an available HFSS project with the expected design. Those script changes
have not been made as part of this environment setup.

## 8. Geant4 and BBRsim (2026-09-30)

BBRsim, the Geant4 package (Mac checkout `/Users/rohanshenoy/BBRsim/BBRSimulation`),
is built against the cluster's own Geant4 module. Nothing Geant4-related is
installed by us, and no administrator rights are needed.

| Need | Module or path |
|---|---|
| Compiler | `gcc/13.2.0-gcc-13.2.0-w55nxkl` (the compiler of Geant4 and of the Python module) |
| Geant4 | `geant4/11.1.2-gcc-13.2.0-gk5zjkn` (multithreaded; loads all datasets via `geant4-data/11.1.0`) |
| CMake | `/usr/bin/cmake` 3.31.8 |
| ROOT | not needed: Geant4 writes ROOT files itself and `uproot` reads them. The only module, `root/6.28.04-gcc-11.3.1`, uses another compiler. |

The module system is Environment Modules: search with `module avail <name>`
(there is no `module spider`).

**`LD_LIBRARY_PATH` trap.** A `~/.bashrc` that appends MATLAB's
`/central/software/Matlab/R2024a/bin/glnxa64` to `LD_LIBRARY_PATH` breaks
`/usr/bin/cmake` (`undefined symbol: uv_fs_get_system_error`, from MATLAB's
older `libuv`) and can shadow Geant4's libraries. Appending to an empty
`LD_LIBRARY_PATH` also leaves an empty entry, which the loader reads as the
current directory. `bb_env.sh` (section 6) removes both.

**Code transfer.** The checkout is `/home/rshenoy/BBRSim/BBRSimulation`, cloned
from a git bundle made on the Mac:

```bash
# Mac
cd /Users/rohanshenoy/BBRsim/BBRSimulation
git bundle create ~/Desktop/bbrsim-main.bundle main bbrsim-V00-01-00
scp ~/Desktop/bbrsim-main.bundle rshenoy@login.hpc.caltech.edu:/home/rshenoy/BBRSim/

# HPC, first time
cd /home/rshenoy/BBRSim && git clone -b main bbrsim-main.bundle BBRSimulation
# HPC, later: copy a new bundle over the old one, then
cd /home/rshenoy/BBRSim/BBRSimulation && git pull --ff-only
```

`bb_env.sh` (section 6) loads the compiler and Geant4 modules, and sources the
installed BBRsim's `bbrsim_env.sh` once `install/` exists.

**Python package**, once, into the `bb` environment:

```bash
source /home/rshenoy/BBRSim/bb_env.sh
python -m pip install -e "/home/rshenoy/BBRSim/BBRSimulation/tools/python[test]"
```

**Build** (in a batch job; `-DWITH_GEANT4_UIVIS=OFF` because there is no display):

```bash
cd /home/rshenoy/BBRSim/BBRSimulation
cmake -S . -B build -DCMAKE_C_COMPILER="$(command -v gcc)" -DCMAKE_CXX_COMPILER="$(command -v g++)" \
      -DCMAKE_BUILD_TYPE=Release -DWITH_GEANT4_UIVIS=OFF -DBUILD_BBRSIM_TESTS=ON \
      -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build -j8 && cmake --install build && ctest --test-dir build -j8
cmake -S examples/testworld -B build/examples/testworld \
      -DCMAKE_C_COMPILER="$(command -v gcc)" -DCMAKE_CXX_COMPILER="$(command -v g++)" \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$PWD/install" -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build/examples/testworld -j8
```

**Verification.** With two compatibility fixes to BBRsim's tests, the build,
install and `ctest` gave 87/87 on Geant4 11.1.2. The library and examples
needed no changes. The fixes:
- `G4iosSetDestination` first appears in Geant4 11.2 (BBRsim commit `a73f5d2`).
- Before 11.4, `G4OpticalPhysics` creates the electron too late in an
  optical-only physics list. This fix was applied on HPC and sent for commit in
  BBRsim.

The first attempt (Slurm job 3691877, before the second fix) passed 75 of 87.

Slurm job 3692888, with the test world built against the install and 8 threads:

| Check | Result |
|---|---|
| `pytest` | 105 passed |
| `planck.mac`, `check_planck_spectrum.py --temp 4` | PASS: ratio obs/theory 1.046; KS p 0.882 |
| `Validation_CrackTransmit.mac`, `check_crack_transmittance.py` | PASS: T 0.50, z = -0.16 sigma; 500 GHz HFSS data on all 40000 crack entries |
| `Validation_CrackTransmit.mac`, `check_crack_ratio.py` | PASS: 2.019 +/- 0.098 vs 1.962 |

Neither run log had any `G4Exception`, `BBR0xx` or `GeomNav` messages.

The merged `bb_env.sh` was then checked in a fresh login shell on login4:
- Geant4 11.1.2, g++ 13.2.0 and cmake 3.31.8 were all available.
- `bbrsim` and PyAEDT both imported.
- `LD_LIBRARY_PATH` held only BBRsim's `install/lib` plus the ROOT entry from
  `~/.bashrc`: no MATLAB, no empty entry and nothing from the modules.
- All of the test-world binary's libraries resolved.

The AEDT smoke test also passed under it (Slurm job 3693291 on hpc-91-13). The
gRPC server was up in 43.9 s, `smoke.aedt` was saved, and the session closed
cleanly.

Not yet run on HPC:
- BBRsim's full `validation/Scripts/run_regression.sh`, which still assumes
  macOS (Apple Clang, `conda run`);
- the remaining fixtures and the light-pipe example.

The fixed-seed reference numbers in BBRsim come from Geant4 11.4 and are not
expected to match 11.1.2.
