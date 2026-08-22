[1]: https://ramses-organisation.readthedocs.io/en/latest
[5]: https://bitbucket.org/vperret/dice
[6]: https://bitbucket.org/ohahn/music
[7]: https://github.com/osyris-project/osyris
[8]: https://github.com/pynbody/pynbody
[9]: https://yt-project.org

# cuRAMSES

**cuRAMSES** is an in-house, performance-optimised fork of the
[RAMSES](https://bitbucket.org/rteyssie/ramses) cosmological AMR code,
developed at the Korea Institute for Advanced Study (KIAS) for the
Horizon Run 5 and HR6 simulations. The cuRAMSES contributions
(recursive *k*-section domain decomposition, Morton-key hash-table
neighbour lookup, FFTW3 direct Poisson solver, spatial hash-binning
for SN/AGN feedback, variable-rank HDF5 restart, and CUDA-accelerated
multigrid/Godunov kernels) are bundled under [`patch/cuRamses/`](./patch/cuRamses/).

### Upstream lineage

cuRAMSES is derived from
[`bitbucket.org/rteyssie/ramses`](https://bitbucket.org/rteyssie/ramses)
master branch at commit
[`d92cd7d9`](https://bitbucket.org/rteyssie/ramses/commits/d92cd7d9)
(2016-11-08, RAMSES version 3.10), as identified by `git`-blob-hash
intersection of all unmodified base-directory source files. The
sub-grid physics (sink particles, AGN feedback, supernova feedback,
metal yields) trace back to the Horizon-AGN version of RAMSES
maintained by Yohan Dubois and collaborators, which is itself
not publicly distributed.

### Companion paper

J. Kim, *cuRAMSES: Scalable AMR Optimizations for Large-Scale
Cosmological Simulations*, MNRAS (2026, in revision).
Cite this paper when using cuRAMSES.

### Acknowledgements

The cuRAMSES contributions in [`patch/cuRamses/`](./patch/cuRamses/)
rest on contributions from many people:

- **Romain Teyssier** — public release and continued development of RAMSES.
- **Yohan Dubois** — Horizon-AGN version of RAMSES, on which the present
  in-house lineage is based. We thank him for sharing his code.
- **Taysun Kimm** — chemical enrichment model (alpha-element and stellar-wind yields).
- **Christophe Pichon, Sébastien Peirani, and Julien Devriendt** — design
  discussions on the Horizon-AGN code and simulation.
- The upstream RAMSES developer community for continued maintenance and
  improvements (see [`CONTRIBUTORS.md`](./CONTRIBUTORS.md) for the
  full list).

cuRAMSES is supported by the Korea Institute for Advanced Study and the
KIAS Center for Advanced Computation. J.~Kim is supported by KIAS
Individual Grants (KG039603).

### Build

A reference build is provided as [`bin/Makefile`](./bin/Makefile) /
[`patch/cuRamses/Makefile`](./patch/cuRamses/Makefile). The default
configuration uses the Intel `ifx` compiler with HDF5 parallel I/O and
optional CUDA acceleration (`make HDF5=1 USE_CUDA=1`). 128-bit Morton
keys for deep-AMR runs are enabled with `MORTON128=1`.

### Bug-fix history

Notable correctness fixes to the cuRAMSES contributions, kept here
because their symptoms (silent data misrouting, restart deadlocks)
are hard to diagnose from run logs alone.

- **2026-08-22 — CUDA multigrid residual convergence.** With GPU
  restriction/interpolation enabled, the CUDA kernels computed valid
  initial and post-smoothing residual norms, but the Fortran driver then
  overwrote both values by recomputing them from stale host-side
  `f(:,1)`. This made the reported relative error remain exactly one,
  forcing every solve to its iteration limit and potentially triggering a
  strict-solver abort. The device solution arithmetic itself was not
  overwritten. The driver now preserves the GPU norms, limits the host
  residual exchange to the CPU and non-GPU-RI fallback paths, limits host
  norm recomputation to the CPU path, and applies one common MPI reduction
  per norm. In a small
  two-rank CUDA/nDGP smoke run of the lagRamses descendant, the original
  path stalled at `Error=1.000` and aborted; after this correction its GPU
  MG error decreased through approximately `8.40e-2`, `6.15e-3`,
  `4.53e-4`, and `3.35e-5`, and both CPU and GPU variants completed. This
  confirms the convergence fix, not full CPU/GPU numerical parity.

- **2026-08-02 — Grid-load-balanced variable-ncpu AMR restart.** Binary
  and HDF5 restarts that changed `ncpu` could fail with `ERROR: out of
  free grids`: the rebuilt decomposition divided Hilbert key space
  uniformly without considering the actual, clustered grid load. A
  same-ncpu restart was unaffected because it read the already balanced
  `bound_key` values from the checkpoint. The binary reader now builds
  cumulative-load-balanced Hilbert or *k*-section boundaries after reading
  grid positions. The streaming HDF5 reader makes one `xg`-only pre-pass
  and uses the same cumulative-load Hilbert boundary calculation; HDF5
  variable-ncpu or cross-ordering restarts into *k*-section now stop
  explicitly instead of using a tree tied to the old decomposition. These
  changes correspond to lagRamses commits `2b62526` (binary) and `cd3d64b`
  (HDF5). In lagRamses validation, a 16-rank checkpoint containing
  2,409,681 grids split across 8, 32, and 64 domains with imbalances
  1.000003, 1.000020, and 1.000073, respectively, with at most five grids
  difference. `econs`, `epot`, `ekin`, and `eint` were identical to all
  printed digits at 16, 8, 32, and 64 ranks, while `mcons` remained about
  `1.8e-16`. At production scale, a 367,389,712-grid checkpoint (100.7M
  grids at level 10, 122.5M at level 11, and 125.0M at level 12) restarted
  on 32 ranks with 11,471,339--11,490,524 grids per domain, imbalance
  1.000836, and no errors. On the same job and nodes, the unmodified binary
  failed under identical conditions after a receiving rank exceeded 13.8M
  grids.

- **2026-07-30 — leaf-aware work load balancing.** Work mode now uses
  the actual particle occupancy of each AMR leaf instead of charging a
  grid's complete particle list to every child. When SIDM is enabled,
  the cost also includes the number of particle pairs sampled by the
  primitive SIDM scatterer; DMO and other runs pay no SIDM proxy cost.
  Sparse rank-by-level wall times update an exponentially smoothed mesh
  and rank correction, while automatic remapping is rate-limited and
  requires its predicted savings to repay the measured remap cost.
  Memory mode remains independent of subcycling, timing, and SIDM work.
  Hilbert, bisection, and k-section all use the same leaf-cost function.
  Verified with a clean serial `make USE_FFTW=1` build and unit tests
  for leaf occupancy, memory/work costs, SIDM gating, and remap
  economics. (The legacy Makefile lacks complete Fortran module
  dependencies, so a clean parallel `make -j` is not supported.)

- **2026-07-29 (`5f8d0a9`) — HDF5 particle restart free-slot
  initialisation.** The HDF5 reader restored active particles
  contiguously into `levelp(1:npart)` but left
  `levelp(npart+1:npartmax)` uninitialised. The next HDF5 checkpoint
  scanned the complete array for positive levels and could mistake
  dirty free slots for particles. The resulting pack-buffer overrun
  corrupted the heap and commonly appeared as a later segmentation
  fault during deallocation. The reader now clears every free slot.
  The writer also verifies that the active-slot count equals `npart`
  before packing. Production HDF5+FFTW tests completed same-rank
  4-to-4 and variable-rank 4-to-2 restart-and-rewrite sequences. Both
  retained 32768 particles. Particle identities, masses, and positions
  were identical, while velocities and gravity agreed within a
  relative tolerance of `1e-11`.

- **2026-07-22 (`5644c40`) — variable-ncpu restore: `cmp_cpumap`
  argument shapes.** The chunked-*k*-section variable-ncpu restores in
  [`patch/cuRamses/init_hydro.f90`](./patch/cuRamses/init_hydro.f90)
  and [`patch/cuRamses/init_poisson.f90`](./patch/cuRamses/init_poisson.f90)
  declared `xx_father(1:1,1:ndim)` and `c_tmp(1:1)` while `cmp_cpumap`
  takes explicit-shape dummies `x(1:nvector,1:ndim)`, `c(1:nvector)`.
  Under Fortran sequence association the dummy element `x(1,2)` maps to
  element `nvector+1` of the 3-element actual, so the `'hilbert'`
  ordering branch read garbage *y*/*z* coordinates and computed wrong
  destination ranks — silent grid misrouting on any restart with
  `ncpu /= ncpu_file`. The `'ksection'` branch (assumed-shape dummy)
  and `init_amr.f90` (already `1:nvector`) were unaffected. Found while
  root-causing a related defect in the lagRamses descendant, where the
  same call site passed a deallocated read buffer instead of the
  position array and hung 32-rank restarts inside `MPI_Waitall`; the
  corrected destination logic was verified there end to end (16-rank
  checkpoint restarted on 32 ranks, psi/hydro/poisson payloads
  delivered with zero misroutes, mass conservation at machine
  precision).

---

![GitHub logo dark-mode-only](./doc/img/full_project_logo_dark.svg#gh-dark-mode-only)
![GitHub logo light-mode-only](./doc/img/full_project_logo.svg#gh-light-mode-only)

## This is the [ramses](https://github.com/ramses-organisation/ramses/) GitHub repository.

Ramses is an open source code to model astrophysical systems, featuring self-gravitating, magnetised, compressible, radiative fluid flows. It is based  on the Adaptive Mesh Refinement (AMR)  technique on a  fully-threaded graded octree.
[ramses](https://github.com/ramses-organisation/ramses/) is written in  Fortran 90 and is making intensive use of the Message Passing Interface (MPI) library.

### ⬇️ Get the code

Download the code by cloning the git repository using
```
$ git clone https://github.com/ramses-organisation/ramses
```

### ℹ️ Get Support

You can go to the user's guide using [Read The Docs][1].
Please register also to the [mailing list](http://groups.google.com/group/ramses_users).
You can get support on Ramses on the [Github Discussion page](https://github.com/ramses-organisation/ramses/discussions)

### 🛠️ Tools and ressources

To generate idealised initial conditions of galaxies, check out [DICE][5].
To generate cosmological initial conditions, check out [MUSIC][6].

To visualize RAMSES data, we encourage you to use [YT][9], [OSYRIS][7] or [PYNBODY][8].

You'll find a lot of useful ressources, links and news about the Ramses community on https://ramses.cnrs.fr/, a website edited by the French "SNO" Ramses.
