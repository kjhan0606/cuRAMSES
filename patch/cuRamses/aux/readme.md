# cuRAMSES Auxiliary Python Scripts

A collection of auxiliary scripts for the cuRAMSES project, organized into four categories: generators, paper figures, analysis plots, and verification tests.

---

## 1. Generators (Table & Configuration Generation)

### `generate_multi_z_table.py`
Generates a multi-redshift Grackle cooling table used by the Eunha cooling module.
```bash
python3 generate_multi_z_table.py --output grackle_multi_z.bin
```
- Dependencies: numpy, struct, subprocess
- Output: `grackle_multi_z.bin` (98 MB, 10 redshift snapshots)

### `generate_de_table.py`
Generates a dark energy transfer function ratio table R_DE(k,a) using the CAMB Weyl potential method.
```bash
python3 generate_de_table.py --w0=-0.8 --wa=0 --cs2=0.0 -o de_table_wcdm.dat
```
- Dependencies: numpy, camb
- Options: `--w0`, `--wa`, `--cs2`, `-o`

### `generate_neutrino_table.py`
Generates a neutrino transfer function ratio table R_nu(k,a) = T_nu/T_cb using CAMB.
```bash
python3 generate_neutrino_table.py neutrino_table.dat
```
- Dependencies: numpy, camb

### `de_linear_test_setup.py`
Batch-generates CAMB transfer functions, MUSIC IC configs, RAMSES namelists, and SLURM scripts for four DE models (LCDM, wCDM, CPL-thaw, CPL-phantom).
```bash
cd test_ksection/run_de_pk_test/
python3 ../../misc/de_linear_test_setup.py
```
- Dependencies: numpy, camb
- Must be run from the target directory

### `ramses_nml_generator.py`
Interactive and CLI tool for generating RAMSES namelist (.nml) files with preset-based configuration.
```bash
python3 ramses_nml_generator.py                                    # Interactive mode
python3 ramses_nml_generator.py --preset hr5_production -o out.nml # Non-interactive
python3 ramses_nml_generator.py --import old.nml -o new.nml        # Import and re-export
python3 ramses_nml_generator.py --list-presets                     # List available presets
```
- Dependencies: none (Python 3.6+ standard library only)
- 5 presets: cosmo_uniform, cosmo_zoomin, restart, gravity_only, hr5_production
- Features: parameter validation (mutual exclusion, dependencies, range checks), existing .nml import/round-trip

---

## 2. Analysis Plots

### Cooling & Star Formation
| Script | Description | Output |
|--------|-------------|--------|
| `plot_cooling_difference.py` | RAMSES vs Grackle cooling rate comparison on (nH, T) diagram | PNG |
| `plot_grackle_table.py` | Grackle binary cooling table visualization (net rate + mean molecular weight) | PNG |
| `plot_sfrd.py` | Cosmic SFRD vs redshift with Madau & Dickinson (2014) and observations | `sfrd_comparison.png` |
| `plot_sfrd_comparison.py` | Multi-run SFRD comparison (holdback / no-holdback / different g*) | PNG |

### Power Spectrum
| Script | Description | Output |
|--------|-------------|--------|
| `compare_pk_camb.py` | Simulation P(k) vs CAMB linear theory comparison | `pk_ratio_*.png`, `growth_factor.png` |

```bash
python3 compare_pk_camb.py --basedir test_ksection/run_de_pk_test/
```

---


## Dependencies

| Package | Required by |
|---------|-------------|
| numpy | Nearly all scripts |
| camb | `generate_de_table.py`, `generate_neutrino_table.py`, `de_linear_test_setup.py`, `compare_pk_camb.py` |
| scipy | `ksection_domain_viz.py`, `ksection_progressive_viz.py` |
| Pillow (PIL) | `ksection_domain_viz.py`, `ksection_progressive_viz.py` |
| (none) | `ramses_nml_generator.py` (Python 3.6+ standard library only) |

```bash
pip install numpy matplotlib camb scipy Pillow  # Install all
```
