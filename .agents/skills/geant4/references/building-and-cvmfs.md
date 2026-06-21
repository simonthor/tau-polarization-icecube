# Building Geant4 Applications — CMake, geant4-config, and CVMFS

How to compile and link a Geant4 v11 application, with an in-depth section on CVMFS installs (the dominant setup on CERN/grid/cluster machines). The CVMFS guidance is the heart of this file — read §3–§6 carefully when `/cvmfs` is involved.

## Table of contents
1. [Standard CMakeLists.txt](#1-standard-cmakeliststxt)
2. [How find_package(Geant4) locates the install; geant4-config](#2-how-find_packagegeant4-locates-the-install-geant4-config)
3. [CVMFS: sourcing the right environment](#3-cvmfs-sourcing-the-right-environment)
4. [CVMFS: include/header paths when compiling](#4-cvmfs-includeheader-paths-when-compiling)
5. [CVMFS: pitfalls (compiler/ABI, datasets, CMAKE_PREFIX_PATH)](#5-cvmfs-pitfalls)
6. [Dataset environment variables](#6-dataset-environment-variables)

---

## 1. Standard CMakeLists.txt

A complete, broadly-compatible template (this is the classic example-B1 style; works against any v11 install):

```cmake
cmake_minimum_required(VERSION 3.16...3.27)
project(MyApp)

# Find Geant4, pulling in all available UI and visualization drivers
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Sets Geant4 header search paths + the compile definitions Geant4 was built with
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(myApp myApp.cc ${sources} ${headers})
target_link_libraries(myApp ${Geant4_LIBRARIES})

install(TARGETS myApp DESTINATION bin)
```

Modern target-style alternative (current example B1) — drops `Geant4_USE_FILE`; the imported `Geant4_LIBRARIES` target carries include dirs as usage requirements:
```cmake
find_package(Geant4 REQUIRED ui_all vis_all)
add_executable(myApp myApp.cc ${sources} ${headers})
target_include_directories(myApp PRIVATE include)
target_link_libraries(myApp PRIVATE ${Geant4_LIBRARIES})
```

Key variables set by `find_package(Geant4)`: `Geant4_INCLUDE_DIRS` (resolves to `<prefix>/include/Geant4`), `Geant4_LIBRARIES`, `Geant4_DEFINITIONS`, `Geant4_USE_FILE`. The words after `REQUIRED` (`ui_all`, `vis_all`) are **component requests** that activate UI/vis driver dependencies.

Configure & build:
```bash
mkdir build && cd build
cmake ..                              # if the env is already sourced (CVMFS view etc.)
# or point CMake at the install explicitly:
cmake -DGeant4_DIR=/path/to/install/lib/cmake/Geant4 ..
cmake --build . -j$(nproc)            # or: make -j$(nproc)
```

---

## 2. How find_package(Geant4) locates the install; geant4-config

- `Geant4Config.cmake` is what `find_package` looks for. It lives under `<prefix>/lib/cmake/Geant4/` (on some builds `lib64/...` — varies by platform).
- If CMake can't find it automatically, point it there with **either**:
  - `-DGeant4_DIR=<dir containing Geant4Config.cmake>`, or
  - `-DCMAKE_PREFIX_PATH=<install prefix>`.
- **`geant4.sh`** (in `<prefix>/bin`): sourcing it puts Geant4 on `PATH`/`LD_LIBRARY_PATH`, makes `geant4-config` available, and **exports all the dataset env vars**. This is the recommended runtime setup.
- **`geant4-config`** is the introspection tool. The most useful invocations:
  ```bash
  geant4-config --prefix     # install prefix
  geant4-config --version    # e.g. 11.2.1
  geant4-config --cflags     # -I<include dirs> and -D<definitions>
  geant4-config --libs       # -L<libdir> -lG4run -lG4event ... (the full link line)
  ```
  `Geant4Config.cmake` is then at `$(geant4-config --prefix)/lib/cmake/Geant4`.

---

## 3. CVMFS: sourcing the right environment

There are **two** distinct CVMFS sources. Don't hardcode versions/platform tags — list the directories to discover what's actually present.

### A) The LCG views at `/cvmfs/sft.cern.ch` — recommended (complete, consistent toolchain)
An LCG "view" is a single directory whose `setup.sh` wires up a mutually compatible **gcc + CMake + ROOT + Geant4 + Python** in one shot. This is the safest option because it guarantees a matching compiler (avoiding the #1 CVMFS pitfall):
```bash
# Discover available releases & platforms:
ls /cvmfs/sft.cern.ch/lcg/views/

# Source one (example — pick the release & platform that match your OS/compiler):
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
```
Platform tag format: `x86_64-<os>-<compiler>-<opt|dbg>`, e.g. `x86_64-el9-gcc13-opt` (`el9` = RHEL9/AlmaLinux9). After sourcing a view, `CMAKE_PREFIX_PATH`/`PATH`/`LD_LIBRARY_PATH` are extended, so `find_package(Geant4)` **just works** — you normally don't set `Geant4_DIR` by hand. Use https://lcginfo.cern.ch to see which Geant4/ROOT version ships in which `LCG_<n>`.

On ATLAS/LHC machines the same is often done via `lsetup "views LCG_106 x86_64-el9-gcc13-opt"`.

### B) The dedicated Geant4 repo `/cvmfs/geant4.cern.ch` — Geant4 only
Gives Geant4 + its datasets, but **not** a matching compiler/CMake/ROOT — you must supply a compatible toolchain yourself (see §5).
```bash
# Discover what's available:
ls /cvmfs/geant4.cern.ch/geant4/

# Path pattern: /cvmfs/geant4.cern.ch/geant4/<version>/<platform-tag>/
# platform tag e.g. x86_64-centos9-gcc11-optdeb-MT  (MT = multithreaded build)

# Source the CMake-based env (sets PATH, lib paths, dataset vars, geant4-config):
source /cvmfs/geant4.cern.ch/geant4/11.2/x86_64-el9-gcc11-optdeb-MT/bin/geant4.sh

# Some installs instead expose the legacy GNUmake env:
# source /cvmfs/geant4.cern.ch/geant4/11.1/x86_64-centos9-gcc11-optdeb-MT/share/Geant4/geant4make/geant4make.sh
```

### Always re-source in each fresh shell / batch job
The environment is per-shell. In a cluster job script, source the setup at the top before `cmake`/`make`/running the binary.

---

## 4. CVMFS: include/header paths when compiling

This is the part people most often get wrong. There are two paths — **use CMake whenever you can.**

### CMake path (recommended)
Once the environment is sourced (§3), nothing special is needed — `find_package(Geant4)` resolves `Geant4_INCLUDE_DIRS` to `/cvmfs/.../include/Geant4` and `include(${Geant4_USE_FILE})` (or the imported target) adds the `-I` flags and compile definitions automatically. You do **not** write `-I/cvmfs/...` yourself.

If CMake picks the wrong Geant4 (e.g. a system one shadows the CVMFS one), pin it:
```bash
cmake -DGeant4_DIR=$(geant4-config --prefix)/lib/cmake/Geant4 ..
```

### Manual g++ path (no CMake) — explicit includes into CVMFS
The cleanest way is to let `geant4-config` emit the flags — it already points into CVMFS:
```bash
# After sourcing the CVMFS/LCG env so geant4-config is on PATH:
g++ -std=c++17 $(geant4-config --cflags) myApp.cc -o myApp $(geant4-config --libs)
```
`$(geant4-config --cflags)` expands to include flags like `-I/cvmfs/geant4.cern.ch/geant4/<ver>/<platform>/include/Geant4` plus `-D` definitions; `$(geant4-config --libs)` expands to `-L/cvmfs/.../lib -lG4run -lG4event -lG4tracking ...` (a long, version-dependent list — another reason to prefer it over hand-listing).

If `geant4-config` isn't on `PATH`, point at the include dir explicitly:
```bash
G4PREFIX=/cvmfs/geant4.cern.ch/geant4/11.2/x86_64-el9-gcc11-optdeb-MT
g++ -std=c++17 -I$G4PREFIX/include/Geant4 myApp.cc -o myApp \
    -L$G4PREFIX/lib -Wl,-rpath,$G4PREFIX/lib \
    -lG4run -lG4event -lG4tracking -lG4processes -lG4digits_hits \
    -lG4track -lG4particles -lG4geometry -lG4materials -lG4global -lG4intercoms
```
The headers live in **`<prefix>/include/Geant4`** (a flat directory of all `G4*.hh`), so a single `-I<prefix>/include/Geant4` covers them. Add `-Wl,-rpath,$G4PREFIX/lib` so the binary finds the CVMFS libraries at runtime without needing `LD_LIBRARY_PATH` re-set (though sourcing `geant4.sh` sets it too). Match the C++ standard to the build (v11 needs C++17; check `geant4-config --cflags`).

---

## 5. CVMFS pitfalls

These cause the large majority of CVMFS build/run failures:

1. **Compiler / ABI mismatch — the #1 issue.** You must compile with the *same gcc family/version* the CVMFS Geant4 was built with (the `gccNN` in the platform tag). Mixing your system gcc with a CVMFS gcc11/gcc13 build causes cryptic link errors or runtime libstdc++ ABI breakage. **Using an LCG view (§3A) avoids this** because the view bundles the matching compiler. With the bare `/cvmfs/geant4.cern.ch` repo (§3B), you must source a matching gcc yourself (e.g. from an LCG release/compiler area).

2. **Source the env *before* configuring.** Run the `source ...setup.sh`/`geant4.sh` before `cmake` and before `g++`, so `CMAKE_PREFIX_PATH`, `PATH`, library paths, and dataset vars are present. A common symptom of forgetting: `find_package` fails, or "G4...DATA not found" at runtime.

3. **`CMAKE_PREFIX_PATH` collisions.** If both a system Geant4 and a CVMFS one are visible, CMake may grab the wrong one. Pin with `-DGeant4_DIR=$(geant4-config --prefix)/lib/cmake/Geant4`. Also watch `lib` vs `lib64` differences across builds.

4. **Dataset version mismatch.** A binary built/linked against one Geant4 version expects specific dataset versions (e.g. `G4EMLOW8.5`). Sourcing the *matching* `geant4.sh`/view sets the correct `G4*DATA` vars. Don't mix a binary of one version with data env vars of another.

5. **Stale build dir after re-sourcing.** If you switch LCG views/Geant4 versions, wipe and re-create the CMake `build/` directory — cached paths from the old environment will point at the wrong CVMFS tree.

---

## 6. Dataset environment variables

Geant4 needs external data files (cross sections, EM low-energy data, gamma levels, etc.). Each dataset has its own env var. Sourcing `geant4.sh` / an LCG view exports them all to absolute CVMFS paths — **don't set these by hand** unless overriding intentionally.

| Variable | Dataset (example dir name) |
|---|---|
| `G4NEUTRONHPDATA` | `G4NDL` (high-precision neutron) |
| `G4LEDATA` | `G4EMLOW` (low-energy EM) |
| `G4LEVELGAMMADATA` | `PhotonEvaporation` |
| `G4RADIOACTIVEDATA` | `RadioactiveDecay` |
| `G4PARTICLEXSDATA` | `G4PARTICLEXS` |
| `G4PIIDATA` | `G4PII` |
| `G4REALSURFACEDATA` | `RealSurface` |
| `G4SAIDXSDATA` | `G4SAIDDATA` |
| `G4ABLADATA` | `G4ABLA` |
| `G4INCLDATA` | `G4INCL` |
| `G4ENSDFSTATEDATA` | `G4ENSDFSTATE` |

Resolution rules:
- If `GEANT4_DATA_DIR` is **set**, it overrides the default search (individual `G4*DATA` vars still take precedence over it). No fallback.
- If unset, Geant4 searches a list of prefixes that **includes `/cvmfs/geant4.cern.ch`** as a built-in fallback — which is why a CVMFS-built binary can often find data even without sourcing anything. Still, source the env to be safe and get the matching versions.

**Diagnosing a "dataset not found" error:** check that the env is sourced (`echo $G4LEDATA`), that the path it points to exists (`ls $G4LEDATA`), and that its version matches the Geant4 version you built against.

---

## Key sources
- Installation Guide (build tools, geant4-config, Geant4_DIR, CMAKE_PREFIX_PATH): https://geant4.web.cern.ch/documentation/dev/ig_html/InstallationGuide/buildtools.html
- Postinstall setup (geant4.sh, dataset env vars, search order): https://geant4.web.cern.ch/documentation/dev/ig_html/InstallationGuide/postinstall.html
- Example B1 CMakeLists.txt: https://github.com/Geant4/geant4/blob/master/examples/basic/B1/CMakeLists.txt
- LCG software builds / views on /cvmfs/sft.cern.ch: https://atlas-software.docs.cern.ch/analysis-software/lcg-software/
- LCG release/platform lookup: https://lcginfo.cern.ch/
- Geant4 forum threads on CVMFS (lxplus9 etc.): https://geant4-forum.web.cern.ch/
