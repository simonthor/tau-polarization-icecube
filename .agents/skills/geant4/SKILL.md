---
name: geant4
description: >-
  Comprehensive guide for writing, building, running, and debugging Geant4
  applications (the C++ Monte Carlo toolkit for simulating particle transport
  through matter, used across HEP, nuclear, medical, and space physics). Use
  this skill whenever the user is working with Geant4 — writing or editing
  detector geometry (G4VUserDetectorConstruction, G4Box/G4Tubs, G4LogicalVolume,
  G4PVPlacement), physics lists (FTFP_BERT, QGSP_BERT, modular physics lists,
  optical/EM physics), primary generators (G4ParticleGun, G4GeneralParticleSource/GPS),
  user action classes, sensitive detectors, scoring/hits, or the G4AnalysisManager;
  setting up main.cc; writing CMakeLists.txt or compiling with geant4-config;
  authoring .mac macro files, visualization (vis.mac, OpenGL/Qt), or /run, /gun,
  /gps, /score, /vis, /control commands; multithreading (G4MTRunManager,
  G4RunManagerFactory, ConstructSDandField); or fixing build/runtime errors like
  missing Geant4 headers, link errors, or G4*DATA dataset errors. ESPECIALLY use
  this skill when Geant4 is installed via CVMFS (/cvmfs/geant4.cern.ch or
  /cvmfs/sft.cern.ch LCG views) and the user needs help sourcing the environment,
  setting include paths, or resolving compiler/ABI mismatches — even if they don't
  say the word "skill".
---

# Geant4

Geant4 is a C++ toolkit for simulating the passage of particles through matter. An application is a C++ program you write and compile against the Geant4 libraries — there is no standalone "geant4" binary. This skill covers the full lifecycle: writing the application, building it (with particular care for CVMFS installs), and steering it at runtime via macro files.

## How to use this skill

The material is split across three reference files so you only load what the task needs. Read the one(s) that match what the user is doing:

- **`references/application-structure.md`** — Writing or editing the C++. The mandatory and optional user classes, a correct v11 `main.cc`, detector construction (materials, solids, volumes, sensitive detectors), physics lists, primary generators (gun & GPS), scoring and the analysis manager, and multithreading rules. Read this for almost any "write/fix Geant4 code" request.
- **`references/building-and-cvmfs.md`** — Compiling and linking. `CMakeLists.txt` templates, `geant4-config`, how `find_package(Geant4)` works, and **the CVMFS workflow in depth**: sourcing the right environment, getting include/library paths right, the manual `g++` compile path, dataset env vars, and the compiler/ABI pitfalls that dominate CVMFS builds. Read this for any build error, CMake question, or anything involving `/cvmfs`.
- **`references/macros-and-visualization.md`** — Runtime control. The UI command tree, batch vs interactive mode, `init_vis.mac`/`vis.mac` templates, `/gun/` vs `/gps/` commands, command-based scoring (`/score/`), visualization drivers, and random-seed reproducibility. Read this for `.mac` files, visualization, or any `/...` command.

When a request spans areas (e.g. "add a calorimeter and write it out to ROOT, then build it on lxplus"), read the relevant files together.

## Orientation: what a Geant4 application is made of

Every application supplies its own `main()` plus three **mandatory** user classes registered on the run manager, and any number of **optional** action classes:

| Role | Base class | Required? |
|---|---|---|
| Geometry & materials | `G4VUserDetectorConstruction` | **Yes** |
| Particles & physics processes | `G4VUserPhysicsList` (usually a reference list like `FTFP_BERT`, or `G4VModularPhysicsList`) | **Yes** |
| Register the action classes | `G4VUserActionInitialization` | **Yes** |
| Generate primary particles | `G4VUserPrimaryGeneratorAction` | **Yes** (registered inside ActionInitialization) |
| Per-run / per-event / per-step hooks | `G4UserRunAction`, `G4UserEventAction`, `G4UserSteppingAction`, `G4UserTrackingAction`, `G4UserStackingAction` | Optional |

The control flow: `main()` creates a run manager (`G4RunManagerFactory::CreateRunManager()` in v11), registers the mandatory classes, optionally starts visualization/UI, then runs events via `/run/beamOn N` (macro) or `runManager->BeamOn(N)` (code). Geometry and physics are frozen once a run starts but may change between runs.

## Key version note (v11 vs v10)

This skill targets **Geant4 v11.x**, the current series. Flag these v11 differences if you see v10-style code:
- Run manager is created via `G4RunManagerFactory::CreateRunManager()` (header `G4RunManagerFactory.hh`), not `new G4RunManager`. The factory auto-selects MT or serial based on the build; override at runtime with the `G4RUN_MANAGER_TYPE` env var (`Serial`/`MT`/`Tasking`/`TBB`).
- `G4AnalysisManager` is a single unified class (alias for `G4GenericAnalysisManager`) via `#include "G4AnalysisManager.hh"` — no more `g4root.hh`/`g4csv.hh`. File format follows the extension or `SetDefaultFileType()`.
- For multithreading, sensitive detectors and fields **must** be built in `DetectorConstruction::ConstructSDandField()` (thread-local), not in `Construct()` (shared). Migrating old single-thread code that built SDs in `Construct()` is a common source of bugs.

## Working principles

- **Prefer the NIST material database** (`G4NistManager::FindOrBuildMaterial("G4_WATER")`) over hand-building materials unless the user needs a custom composition.
- **Geant4 solids take half-lengths**, not full dimensions — `new G4Box("b", dx/2, dy/2, dz/2)`. This trips people up constantly; double-check it when reading or writing geometry.
- **Always enable overlap checking** (`true` as the last arg of `G4PVPlacement`) while developing geometry — overlaps cause silent tracking errors.
- **Use a reference physics list** (`FTFP_BERT` is the default for HEP) unless there's a specific reason to build a custom modular list. None of the reference lists include optical photons — add `G4OpticalPhysics` on top if needed.
- **Don't hardcode dataset paths.** The `G4*DATA` environment variables should come from sourcing `geant4.sh` (or the CVMFS/LCG setup). A "dataset not found" runtime error almost always means the environment wasn't sourced or a version mismatch — see the CVMFS reference.
- When in doubt about an API, the canonical reference is the official examples (`examples/basic/B1`–`B5`, `examples/extended/...` in the Geant4 source tree) and the Book For Application Developers at https://geant4-userdoc.web.cern.ch.
