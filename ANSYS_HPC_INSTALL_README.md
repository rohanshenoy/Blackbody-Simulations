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

Dependencies are installed by pip; these commands do not pin versions. This
does not yet install the local `bbrsim` Python package from BBRSimulation.
Geant4 module selection, the BBRSimulation build, and its Python package
installation are deferred.

## 6. Combined setup for each session

Create `/home/rshenoy/BBRSim/bb_env.sh` with this content:

```bash
# Load the Python installation used to create the environment.
module load python/3.11.6-gcc-13.2.0-fh6i4o3 || return

# Activate Python packages.
source /home/rshenoy/BBRSim/envs/bb/bin/activate || return

# Load Ansys paths and licensing.
source /resnick/groups/golwala/rshenoy/ansys_env.sh || return
export ANSYSEM_ROOT252="$ANSYS_ROOT/AnsysEM"
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

**When starting work with the Geant4-based BBRSimulation package on HPC,
update `bb_env.sh`** to load the selected, compatible Geant4/compiler modules
and source the installed BBRsim `bbrsim_env.sh` after building and installing
the package. Also install its Python analysis package into `bb`. The current
combined setup covers Python and Ansys only; Geant4 compatibility and the
BBRsim installation paths still need to be established.

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
and the Python import/path checks above passed. Non-graphical AEDT startup,
an actual HFSS license checkout, and a simulation have **not yet been
verified**. An X-display error or successful DNS lookup does not verify a
license checkout. The next step is an AEDT startup test in an interactive
compute allocation.

Before batch execution, Blackbody-Simulations still needs its graphical
`Hfss(..., non_graphical=False)` launches adapted for non-graphical use, plus
an available HFSS project with the expected design. Those script changes
have not been made as part of this environment setup.
