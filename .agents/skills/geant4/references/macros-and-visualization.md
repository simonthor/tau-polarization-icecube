# Geant4 Macros, UI Commands & Visualization (v11.x)

Runtime control of a Geant4 application: the UI command tree, `.mac` macro files, primary-source commands, command-based scoring, visualization, and reproducible random seeding. Command names and example macros verified against the official Book For Application Developers and the canonical `examples/basic/B1` macros.

## Table of contents
1. [The UI command system](#1-the-ui-command-system)
2. [Batch vs interactive; common commands](#2-batch-vs-interactive-common-commands)
3. [Macro patterns (init_vis.mac, run.mac, loops)](#3-macro-patterns)
4. [/gun/ vs /gps/ primary source commands](#4-gun-vs-gps-primary-source-commands)
5. [Command-based scoring /score/](#5-command-based-scoring-score)
6. [Visualization](#6-visualization)
7. [Random seeding & reproducibility](#7-random-seeding--reproducibility)

---

## 1. The UI command system

Commands reach the kernel three ways:
1. **Interactively** in a (G)UI session.
2. **From a macro file** via `/control/execute <file>`.
3. **From C++**: `G4UImanager::GetUIpointer()->ApplyCommand("/run/beamOn 100");`

Commands are organized in a directory tree by category: `/control/` (flow, aliases), `/run/`, `/event/`, `/tracking/`, `/gun/` (only if the app uses `G4ParticleGun`), `/gps/` (only if it uses `G4GeneralParticleSource`), `/process/`, `/particle/`, `/material/`, `/geometry/`, `/units/`, `/vis/`, `/score/`, `/random/`. In a session, `help` lists everything and `/control/manual /dir/` dumps a directory.

---

## 2. Batch vs interactive; common commands

Mode is chosen by `main()` from `argc` (see `application-structure.md`): a macro argument → batch; no argument → interactive GUI.

**`/control/`**
- `/control/execute <file>` — run a macro.
- `/control/verbose <0|1|2>` — echo applied commands.
- `/control/alias <name> <value>` / use as `{name}`; `/control/listAlias`.
- `/control/loop <macro> <counter> <init> <final> [step]` — repeat a macro, counter usable as `{counter}`.
- `/control/foreach <macro> <counter> <"v1 v2 v3">` — run once per value.
- `/control/getEnv <var>` — import a shell env var as an alias.
- `/control/shell <cmd>`, `/control/echo <str>`, `/control/ifBatch <macro>`, `/control/ifInteractive <macro>`.

**`/run/`**
- `/run/initialize` — initialize geometry + physics.
- `/run/beamOn <N> [macro] [nSelect]` — run N events.
- `/run/verbose <0-4>`, `/run/printProgress <mod>`.
- `/run/numberOfThreads <n>` (MT, PreInit only).
- `/run/setCut <value> [unit]`, `/run/setCutForAGivenParticle <particle> <value> <unit>`.

**Verbosity elsewhere:** `/tracking/verbose <n>`, `/event/verbose <n>`, `/vis/verbose <level>`.

---

## 3. Macro patterns

### `init_vis.mac` (run at startup in interactive mode — verbatim from example B1)
```
# Set some default verbose
/control/verbose 2
/control/saveHistory
/run/verbose 2
#
# Change the default number of threads (in multi-threaded mode)
#/run/numberOfThreads 4
#
/run/initialize
#
/control/execute vis.mac
```

### `run.mac` (typical batch run)
```
/run/verbose 2
/gun/particle e-
/gun/energy 1 GeV
/run/beamOn 100
```

### Energy scan with a loop
```
# scan.mac, invoked as: /control/foreach scan.mac ENERGY "1 2 5 10 20"
/gun/energy {ENERGY} MeV
/run/beamOn 1000
```

---

## 4. /gun/ vs /gps/ primary source commands

`/gun/` and `/gps/` are mutually exclusive — only the one matching your primary generator exists.

### `/gun/` — G4ParticleGun (simple fixed beam)
- `/gun/particle <name>` (e.g. `e-`, `proton`, `gamma`, `geantino`)
- `/gun/energy <E> <unit>` (kinetic energy)
- `/gun/direction <ex> <ey> <ez>` (need not be normalized)
- `/gun/position <X> <Y> <Z> <unit>`
- `/gun/momentum <px> <py> <pz> <unit>`, `/gun/momentumAmp <p> <unit>`
- `/gun/number <N>` (particles per event), `/gun/time <t> <unit>`, `/gun/polarization <Px> <Py> <Pz>`
- `/gun/ion <Z> <A> [Q] [E_keV]`

```
/gun/particle proton
/gun/energy 150 MeV
/gun/position 0 0 -20 cm
/gun/direction 0 0 1
/run/beamOn 1000
```

### `/gps/` — G4GeneralParticleSource (rich distributions)
**Position** `/gps/pos/`: `type <Point|Plane|Beam|Surface|Volume>`, `shape <Circle|Square|Sphere|Cylinder|...>`, `centre <X Y Z unit>`, `halfx/halfy/halfz`, `radius`, `confine <physVolName>`.
**Angular** `/gps/ang/`: `type <iso|cos|beam1d|beam2d|focused|planar>`, `mintheta/maxtheta/minphi/maxphi <val unit>`.
**Energy** `/gps/ene/`: `type <Mono|Lin|Pow|Exp|Gauss|Brem|Bbody|User|Arb>`, `mono <E unit>`, `min/max <E unit>`, `sigma`, `alpha` (power law), `temp`, `gradient`/`intercept` (linear).
**Multiple sources** `/gps/source/`: `add <intensity>`, `set <i>`, `list`.
Shortcuts: `/gps/particle`, `/gps/energy`, `/gps/direction`, `/gps/position`, `/gps/number`.

Point source, isotropic, mono-energetic:
```
/gps/particle gamma
/gps/pos/type Point
/gps/pos/centre 0 0 0 cm
/gps/ang/type iso
/gps/ene/type Mono
/gps/ene/mono 1 MeV
/run/beamOn 100000
```

Isotropic power-law spectrum (e.g. cosmic-ray-like):
```
/gps/particle proton
/gps/ang/type iso
/gps/ene/type Pow
/gps/ene/min 10 MeV
/gps/ene/max 1 GeV
/gps/ene/alpha -2.7
/run/beamOn 100000
```

---

## 5. Command-based scoring /score/

Overlay a scoring mesh without writing a sensitive detector. Workflow: create mesh → set size/bins/position → declare quantities (each optionally followed by a filter that applies to the preceding scorer) → `/score/close` → run → dump/draw.

```
/score/create/boxMesh boxMesh_1
/score/mesh/boxSize 100. 100. 100. cm        # half-sizes
/score/mesh/nBin 30 30 30
#
/score/quantity/energyDeposit eDep
/score/quantity/doseDeposit dose
/score/quantity/nOfStep nOfStepGamma
/score/filter/particle gammaFilter gamma     # applies to nOfStepGamma above
/score/close
#
/run/beamOn 100000
#
/score/dumpQuantityToFile boxMesh_1 eDep eDep.csv
/score/drawProjection boxMesh_1 eDep
```
Other quantities: `cellFlux`, `trackLength`, `nOfSecondary`, `flatSurfaceFlux`. Filters: `/score/filter/kineticEnergy <name> <Emin> <Emax> <unit>`. Cylinder mesh: `/score/create/cylinderMesh` + `/score/mesh/cylinderSize <R> <dZ> <unit>`. Output is CSV by default; enable ROOT ntuple output by instantiating `G4TScoreNtupleWriter` in `main()`.

---

## 6. Visualization

### Drivers and `/vis/open`
Common driver strings: `OGL` (generic OpenGL, auto-picks a variant), `Qt`/`OGLSQt` (interactive Qt GUI — mouse zoom/rotate, image export), `TSG`/`TSG_OFFSCREEN` (ToolsSG, offscreen for batch images), `VtkNative`/`VtkQt`, `HepRepFile` (HepRApp XML), `DAWNFILE` (`.prim` → DAWN), `VRML2FILE` (`.wrl`), `RayTracer` (JPEG), `ASCIITree` (text tree via `/vis/drawTree`). Default selectable via the `G4VIS_DEFAULT_DRIVER` env var.

### `vis.mac` (verbatim-derived from example B1 — safe to ship)
```
# Open a viewer (no arg = system choice; or e.g. /vis/open OGL)
/vis/open
#
# Quieten and defer redraws while building the scene
/vis/viewer/set/autoRefresh false
/vis/verbose errors
#
# Draw the geometry
/vis/drawVolume
#
# View angle and style
/vis/viewer/set/viewpointThetaPhi 120 150
/vis/viewer/set/style surface
/vis/viewer/set/auxiliaryEdge true
/vis/viewer/set/lineSegmentsPerCircle 100
#
# Draw smooth trajectories at end of event, colored by charge, with step points
/vis/scene/add/trajectories smooth
/vis/modeling/trajectories/create/drawByCharge
/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true
/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2
#
# Draw hits at end of event
#/vis/scene/add/hits
#
# Superimpose all events of a run
/vis/scene/endOfEventAction accumulate
#
# Decorations
/vis/scene/add/axes
/vis/scene/add/scale
/vis/scene/add/eventID
/vis/scene/add/logo2D
#
# Re-enable auto refresh
/vis/viewer/set/autoRefresh true
/vis/verbose warnings
# For file-based drivers, force a render:
#/vis/viewer/flush
```

Useful extras:
- Color a volume: `/vis/geometry/set/colour Envelope 0 0 0 1 .3` (name depth r g b alpha).
- Hide a volume: `/vis/geometry/set/visibility World 0 false`.
- Color trajectories by particle: `/vis/modeling/trajectories/create/drawByParticleID` then `/vis/modeling/trajectories/drawByParticleID-0/set e+ yellow`.
- Show only gammas: `/vis/filtering/trajectories/create/particleFilter` then `.../particleFilter-0/add gamma`.
- Note auto-indexed model names: first `drawByCharge` model is `drawByCharge-0`, first particle filter is `particleFilter-0`.
- If too many tracks crash the viewer: `/tracking/storeTrajectory 0`.

---

## 7. Random seeding & reproducibility

UI commands (`/random/`):
- `/random/setSeeds <s1> <s2>` — seed the engine (give at least two integers).
- `/random/setDirectoryName <dir>` and `/random/setSavingFlag true` — save engine status at the start of each run (`currentRun.rndm`) and event (`currentEvent.rndm`).
- `/random/saveEachEventFlag true` — save per-event status to `runXXXevtYYY.rndm`.
- `/random/resetEngineFrom <file>` — restore engine status (re-run a specific event).

C++ equivalents:
```cpp
G4Random::setTheSeed(12345);
G4long seeds[2] = {12345, 67890};
G4Random::setTheSeeds(seeds);
G4Random::saveEngineStatus("myEngine.rndm");
G4Random::restoreEngineStatus("myEngine.rndm");
```

Reproducibility macro:
```
/random/setSeeds 123456789 987654321
/random/setDirectoryName randomStatus
/random/setSavingFlag true
/run/beamOn 100000
# Reproduce one event later:
#/random/resetEngineFrom randomStatus/currentEvent.rndm
#/run/beamOn 1
```

Notes:
- **Command-line seeding is an application convention**, not a built-in UI command — the standard examples parse `argv` and call `G4Random::setTheSeeds`. Don't expect a `/random/` command for it.
- **MT reproducibility:** the master manages per-event seeds, so with the default `seedOnce=0`, results are reproducible regardless of thread count.

---

## Key sources
- Built-in commands: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Control/commands.html
- Visualization: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Visualization/
- General Particle Source: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/GettingStarted/generalParticleSource.html
- Command-based scoring: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Detector/commandScore.html
- Example B1 macros: https://github.com/Geant4/geant4/tree/master/examples/basic/B1
