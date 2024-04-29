# Tau decay simulations in ice
This directory contains code for simulating tau decays in ice. There are numerous important features with this program:
- The Geant4 version used is 11.2.0, i.e. the absolute newest version of Geant4. This way, it can be tested and compared against other Geant4 versions (e.g. the one used in the DeepCore simulations, which is 10.x.x), to identify potential relevant bug fixes in Geant4 (or lack thereof)
- Tauola has been integrated into Geant4 to handle all of the $\tau$ decays. This makes $\tau$ decays as accurate as possible, even when using Geant4 (which otherwise simulates tau decays quite poorly). See [tauola integration](#tauola-integration) for more details.
- The decay products of the $\tau$ are extracted in the exact same way as in the DeepCore simulations, using a `UserTrackingAction`. This way, by adding or removing the other features, it is possible to investigate the differences that certain features have on the Geant4 output compared to the DeepCore simulations
- The tau 4-momenta are read from a file automatically. This means that GENIE tau outputs can be directly passed into Geant4 as a csv file, which makes it easy to integrate Geant4 with GENIE. Note that this is kind of hard-coded, see [how to run](#how-to-run).

## How to run
The best way to run the software is by using [generate_events.sh](../icecube-tau-decays/generate_events.sh) in [icecube-tau-decays/](../icecube-tau-decays/). `generate_events.sh` is controlled by editing various parameters in the [settings.yml](../icecube-tau-decays/settings.yml) file. More documentation is available in [generate_events.sh](../icecube-tau-decays/generate_events.sh). In short, to only run Geant4, set `start_step: 5` and run
```bash
./generate_events.sh
```
Note that this must be (or at least it is strongly recommended to) executed from the [icecube-tau-decays/](../icecube-tau-decays/) directory. This will change if/when the hard-coding of file names is removed.

## Output
The output is a file in [data](../data/) with the name `geant4_run_e{energy}_tauola.csv` (or some similar variation, see the source code in [TractingAction](./src/TrackingAction.cc), as this is also hard-coded). `{energy}` in the file name is replaced with the energy value specified in [settings.yml](../icecube-tau-decays/settings.yml). 

A log file is also produced, with the name `geant4_run_e{evenrgy}_tauola.log` in [logfiles/](../logfiles/).

## Tauola integration
For Tauola integration into an arbitrary Geant4 application, only a few simple modifications must be done.
1. Copy [TauolaDecayer.cc](./src/TauolaDecayer.cc), [TauolaDecayerPhysics.cc](./src/TauolaDecayerPhysics.cc) to your `src/` directory.
2. Copy [TauolaDecayer.hh](./include/TauolaDecayer.hh), [TauolaDecayerPhysics.hh](./include/TauolaDecayerPhysics.hh) to your `include/` directory.
3. In your `PhysicsList` class (see [PhysicsList.cc](./src/PhysicsList.cc) for an example), add the following two lines. Only having the `RegisterPhysics(new TauolaDecayerPhysics());` at the end of all `RegisterPhysics` calls has been tested, and is known to work.

```c++
// Your existing #include's
#include "TauolaDecayerPhysics.hh"
// your existing code...

PhysicsList::PhysicsList() {
    // your existing code...
    RegisterPhysics(new TauolaDecayerPhysics());
}
// your existing code...
```
4. Add the file [FindTauola++.cmake](./FindTauola++.cmake) to your CMake module path, and add the line `find_package(Tauola++ REQUIRED)` to your CMakeLists.txt file (I have not managed to get this to work).
5. If the above step does not work, copy the content of [FindTauola++.cmake](./FindTauola++.cmake) to your CMakeLists.txt file instead.
6. In CMakeLsits.txt, there should be a line where `target_link_libraries()` is called. In it, add  `${Tauola++_LIBRARIES}`. E.g. `target_link_libraries(exampleB1 ${Geant4_LIBRARIES})` becomes `target_link_libraries(exampleB1 ${Geant4_LIBRARIES}  ${Tauola++_LIBRARIES})`.
7. Compile the Geant4 application using `cmake ... -DTauola++_DIR="path/to/Tauola/"`, where ... is replaced by all command line arguments that is usually passed to cmake when compiling this specific application.

These are the only changes needed. After this, one can compile and run the Geant4 application as usual (i.e. using `make`). All $\tau$ leptons should now be decayed using Tauola.

### Limitations to Tauola integration
Currently, only unpolarized decays with Tauola are supported. It is easy to change the code for the tau leptons to e.g. always have a fully left-handed polarization. This can be done by changing the `polx, poly, polz` parameters to the negative values of the tau momentum. A comment in the code describes how to do this.

Calculating and passing a realistic polarization has not been implemented yet, as the integration on this depends highly on the rest of the software stack and how one can pass a polarization vector to Geant4.

Tauola does not work in with multithreading. It might be possible to circumvent this in some way, but officially there is no support.

<!--TODO add more text about how the rest of the software works-->