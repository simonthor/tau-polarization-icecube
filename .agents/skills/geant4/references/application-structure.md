# Geant4 Application Structure (v11.x)

Everything needed to write or fix the C++ of a Geant4 application. All class names and signatures verified against the official Geant4 v11 *Book For Application Developers*.

## Table of contents
1. [User classes overview](#1-user-classes-overview)
2. [main.cc skeleton](#2-maincc-skeleton)
3. [Detector construction](#3-detector-construction)
4. [Physics lists](#4-physics-lists)
5. [Primary generator](#5-primary-generator)
6. [User action classes](#6-user-action-classes)
7. [Scoring & extracting output](#7-scoring--extracting-output)
8. [Multithreading](#8-multithreading)

---

## 1. User classes overview

A Geant4 application supplies its own `main()` plus user classes derived from Geant4 base classes. Two families:

**Initialization classes** (registered on the run manager via `SetUserInitialization`):
- `G4VUserDetectorConstruction` — geometry, materials, sensitive detectors, fields. Implement `Construct()` (returns the world `G4VPhysicalVolume*`) and, for MT, `ConstructSDandField()`.
- `G4VUserPhysicsList` — particles + physics processes + production cuts. Almost always you use a prebuilt reference list (`FTFP_BERT`) or derive from `G4VModularPhysicsList`.
- `G4VUserActionInitialization` — registers all the action classes. Implement `Build()` (called per worker thread) and `BuildForMaster()` (master only, MT).

**Action classes** (registered inside `ActionInitialization::Build()` via `SetUserAction`):
- `G4VUserPrimaryGeneratorAction` — **mandatory**; `GeneratePrimaries(G4Event*)`.
- `G4UserRunAction` — `BeginOfRunAction`, `EndOfRunAction`, `GenerateRun`.
- `G4UserEventAction` — `BeginOfEventAction`, `EndOfEventAction`.
- `G4UserTrackingAction` — `PreUserTrackingAction`, `PostUserTrackingAction`.
- `G4UserSteppingAction` — `UserSteppingAction(const G4Step*)`.
- `G4UserStackingAction` — `ClassifyNewTrack`, `NewStage`, `PrepareNewEvent`.

The **run manager** (`G4RunManager`, or `G4MTRunManager` for MT) is the only manager you construct explicitly; it creates and owns the other managers and the event loop. `Initialize()` builds geometry and physics; `BeamOn(n)` runs `n` events.

### ActionInitialization pattern
```cpp
void ActionInitialization::Build() const {        // per worker thread (and in serial mode)
  SetUserAction(new PrimaryGeneratorAction);      // mandatory
  SetUserAction(new RunAction);
  SetUserAction(new EventAction);
  SetUserAction(new SteppingAction);
}
void ActionInitialization::BuildForMaster() const { // MT: master thread only
  SetUserAction(new RunAction);                     // typically only RunAction, for merging
}
```

---

## 2. main.cc skeleton

The v11 idiom: create the run manager via `G4RunManagerFactory`, register the three mandatory classes, and branch on `argc` to pick interactive vs batch mode.

```cpp
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "FTFP_BERT.hh"               // a reference physics list

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

int main(int argc, char** argv)
{
  // Interactive if no macro argument was given
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc, argv);

  // v11: factory auto-selects MT or serial based on how the libs were built.
  auto* runManager = G4RunManagerFactory::CreateRunManager();

  // Mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new FTFP_BERT);          // physics list
  runManager->SetUserInitialization(new ActionInitialization);

  // Visualization
  auto* visManager = new G4VisExecutive;
  visManager->Initialize();

  auto* UImanager = G4UImanager::GetUIpointer();
  if (!ui) {                                                 // batch mode
    UImanager->ApplyCommand("/control/execute " + G4String(argv[1]));
  } else {                                                   // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;   // deletes all other managers; do this last
  return 0;
}
```

Run as `./myApp run.mac` (batch) or `./myApp` (interactive GUI).

---

## 3. Detector construction

Implement `Construct()` to build the geometry tree and return the world physical volume.

### Materials — prefer the NIST database
```cpp
#include "G4NistManager.hh"
auto* nist  = G4NistManager::Instance();
auto* air   = nist->FindOrBuildMaterial("G4_AIR");
auto* water = nist->FindOrBuildMaterial("G4_WATER");
auto* pb    = nist->FindOrBuildMaterial("G4_Pb");
auto* vac   = nist->FindOrBuildMaterial("G4_Galactic");   // good for the world
```
Custom materials when needed:
```cpp
auto* elH = nist->FindOrBuildElement("H");
auto* elO = nist->FindOrBuildElement("O");
auto* h2o = new G4Material("MyWater", 1.0*g/cm3, 2 /*ncomponents*/);
h2o->AddElement(elH, 2);     // by number of atoms
h2o->AddElement(elO, 1);
// or by mass fraction: mat->AddElement(el, 0.7);
// single-element: new G4Material("Al", 13., 26.98*g/mole, 2.70*g/cm3);
```

### Solids — REMEMBER: half-lengths
```cpp
auto* box  = new G4Box("Box", dx/2, dy/2, dz/2);                  // half-extents!
auto* tube = new G4Tubs("Tube", rMin, rMax, halfZ, startPhi, deltaPhi);
auto* sph  = new G4Sphere("Sph", rMin, rMax, sPhi, dPhi, sTheta, dTheta);
// Booleans: G4UnionSolid, G4IntersectionSolid, G4SubtractionSolid
```

### Logical & physical volumes
```cpp
auto* boxLV = new G4LogicalVolume(box, water, "BoxLV");

// G4PVPlacement(rotation, translation, logicalVol, name, motherLV, pMany, copyNo, checkOverlaps)
new G4PVPlacement(nullptr, G4ThreeVector(0,0,5*cm), boxLV, "Box",
                  worldLV, false, 0, /*checkOverlaps=*/true);
```

### Full Construct()
```cpp
G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist     = G4NistManager::Instance();
  auto* worldMat = nist->FindOrBuildMaterial("G4_Galactic");

  auto* worldS  = new G4Box("World", 1*m, 1*m, 1*m);          // half-extents
  auto* worldLV = new G4LogicalVolume(worldS, worldMat, "World");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World",
                                    nullptr /*no mother*/, false, 0, true);
  // ... place daughter volumes with worldLV as the mother ...
  return worldPV;                                             // MUST return the world
}
```

### Sensitive detectors
Attach SDs in `ConstructSDandField()` so they are created per worker thread (thread-safe under MT):
```cpp
void DetectorConstruction::ConstructSDandField() {
  auto* sd = new MySensitiveDetector("MySD", "MyHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  SetSensitiveDetector("BoxLV", sd);   // attach by logical-volume name
}
```
A custom SD derives from `G4VSensitiveDetector` and implements `ProcessHits(G4Step*, G4TouchableHistory*)`, filling a `G4THitsCollection<MyHit>`.

---

## 4. Physics lists

### Easiest: use a reference list
```cpp
#include "FTFP_BERT.hh"
runManager->SetUserInitialization(new FTFP_BERT);   // HEP default since v10
```
Common reference lists: `FTFP_BERT` (default), `FTFP_BERT_HP` (with thermal-neutron high-precision), `QGSP_BERT`, `QGSP_BIC`, `QGSP_BIC_HP`, `Shielding`, `QBBC`. **None include optical photons.**

### Select by name via factory
```cpp
#include "G4PhysListFactory.hh"
G4PhysListFactory factory;
auto* pl = factory.GetReferencePhysList("FTFP_BERT_EMZ");
runManager->SetUserInitialization(pl);
```

### EM physics variants (suffixes)
| Suffix | Constructor | Use |
|---|---|---|
| (none) | `G4EmStandardPhysics` | default |
| `_EMV`/`_EMX` | option1/option2 | faster, less precise |
| `_EMY` | option3 | accurate (medical/low-E) |
| `_EMZ` | option4 | most accurate |
| `_LIV` | `G4EmLivermorePhysics` | low-energy Livermore |
| `_PEN` | `G4EmPenelopePhysics` | Penelope |

### Custom modular list (e.g. to add optical photons)
```cpp
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4OpticalPhysics.hh"

class PhysicsList : public G4VModularPhysicsList {
 public:
  PhysicsList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4OpticalPhysics());
  }
};
```
Or add modules to a reference list directly:
```cpp
auto* phys = new FTFP_BERT;
phys->RegisterPhysics(new G4OpticalPhysics());   // optical on top
```
Modular API: `RegisterPhysics`, `ReplacePhysics` (swap, e.g. EM), `RemovePhysics`.

---

## 5. Primary generator

### G4ParticleGun — fixed kinematics
```cpp
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction() {
  fGun = new G4ParticleGun(1 /*particles per event*/);
  auto* p = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fGun->SetParticleDefinition(p);
  fGun->SetParticleEnergy(6.*MeV);
  fGun->SetParticlePosition(G4ThreeVector(0,0,-10*cm));
  fGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
}
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* evt) {
  fGun->GeneratePrimaryVertex(evt);   // can randomize per event here
}
// destructor: delete fGun;
```

### G4GeneralParticleSource (GPS) — macro-driven distributions
Use it identically in code (`fGPS->GeneratePrimaryVertex(evt)`), then configure spatial/angular/energy distributions and multiple sources from macros (see `macros-and-visualization.md`). Preferred when you need energy spectra, isotropic emission, surface/volume sources, etc.

---

## 6. User action classes

These are hooks where you collect data. Common pattern — accumulate energy per event, then per run:

```cpp
// EventAction
void EventAction::BeginOfEventAction(const G4Event*) { fEdep = 0.; }
void EventAction::EndOfEventAction(const G4Event*) {
  fRunAction->AddEdep(fEdep);                 // hand off to RunAction
}
// SteppingAction
void SteppingAction::UserSteppingAction(const G4Step* step) {
  fEventAction->AddEdep(step->GetTotalEnergyDeposit());
}
// RunAction
void RunAction::EndOfRunAction(const G4Run* run) {
  G4int n = run->GetNumberOfEvent();
  if (n == 0) return;
  // report / write accumulated quantities
}
```
Under MT, prefer accumulating into the analysis manager or `G4Accumulable` so results merge correctly across threads.

---

## 7. Scoring & extracting output

Four complementary mechanisms — pick by need:

### (a) Sensitive detectors + hits collections
Full per-step control via a custom `G4VSensitiveDetector` building `G4THitsCollection<MyHit>`. Most flexible, most code.

### (b) G4MultiFunctionalDetector + primitive scorers (low/no custom code)
```cpp
auto* det = new G4MultiFunctionalDetector("crystal");
G4SDManager::GetSDMpointer()->AddNewDetector(det);
det->RegisterPrimitive(new G4PSEnergyDeposit("edep"));
det->RegisterPrimitive(new G4PSTrackLength("trackLength"));
crystalLV->SetSensitiveDetector(det);
// retrieve: G4SDManager::GetSDMpointer()->GetCollectionID("crystal/edep")
```
Scorers include `G4PSEnergyDeposit`, `G4PSDoseDeposit`, `G4PSTrackLength`, `G4PSNofStep`, `G4PSCellFlux`, `G4PSFlatSurfaceCurrent`. Results auto-merge across MT workers. Filters via `scorer->SetFilter(new G4SDParticleFilter(...))`.

### (c) Command-based scoring (no C++) — `/score/...` macros
See `macros-and-visualization.md`. Good for quick dose/flux meshes overlaid on geometry.

### (d) G4AnalysisManager — histograms & ntuples to ROOT/CSV/HDF5/XML
This is usually what physicists want for downstream analysis. Single unified class in v11:
```cpp
#include "G4AnalysisManager.hh"

// RunAction constructor / setup
auto* man = G4AnalysisManager::Instance();
man->SetDefaultFileType("root");        // or csv, hdf5, xml
man->SetNtupleMerging(true);            // MT: merge per-worker ROOT ntuples to master
man->CreateH1("Eabs", "Edep in absorber", 100, 0., 800*MeV);   // returns histogram id
man->CreateNtuple("data", "Edep and TrackL");
man->CreateNtupleDColumn("Eabs");       // also I/F/S column variants
man->FinishNtuple();

// BeginOfRunAction
man->OpenFile("output");                // -> output.root

// EventAction / fill point
man->FillH1(0, energyAbs);
man->FillNtupleDColumn(0, energyAbs);
man->AddNtupleRow();

// EndOfRunAction
man->Write();
man->CloseFile();
```
Note: ntuple merging across threads works for ROOT only (`SetNtupleMerging(true)`); CSV/HDF5/XML produce per-thread files. Histograms merge automatically on `Write()`.

---

## 8. Multithreading

- **Run manager**: `G4MTRunManager` (or just let `G4RunManagerFactory` choose). Set threads with `/run/numberOfThreads N` (PreInit phase) or the `G4FORCENUMBEROFTHREADS` env var. The master cannot use `SetUserAction()` — all actions go through `ActionInitialization`.
- **Shared (built once on master):** geometry (solids, logical/physical volumes) and the material table — built in `Construct()`.
- **Thread-local (built per worker):** sensitive detectors and fields — **must** be in `ConstructSDandField()`. This is the #1 migration issue from old single-thread code.
- **Merging results:** primitive-scorer hits maps and analysis-manager histograms merge automatically; ROOT ntuples merge with `SetNtupleMerging(true)`; for custom run data use `G4Run::Merge()` / `GenerateRun()`, or `G4Accumulable` registered with `G4AccumulableManager`.
- **Reproducibility:** the master manages per-event seeds so results are reproducible regardless of thread count (default `seedOnce=0`).

Switching MT/serial needs no recompile — set `G4RUN_MANAGER_TYPE=Serial` (or `MT`/`Tasking`/`TBB`) at runtime.

---

## Key sources
- Book For Application Developers: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/
- Guide For Physics Lists: https://geant4.web.cern.ch/documentation/dev/plg_html/PhysicsListGuide/physicslistguide.html
- Canonical examples (B1–B5, extended/): https://github.com/Geant4/geant4/tree/master/examples
