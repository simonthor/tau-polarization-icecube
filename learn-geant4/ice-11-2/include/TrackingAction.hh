#ifndef TrackingAction_h
#define TrackingAction_h 1


#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include <mutex>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace B1 {
class EventAction;

class TrackingAction : public G4UserTrackingAction {

  public:  
    TrackingAction(EventAction* eventAction);
    
   ~TrackingAction() {};
   
    void PreUserTrackingAction(const G4Track*);   
    void PostUserTrackingAction(const G4Track*) {};

  private:
    EventAction* fEventAction = nullptr;
    std::mutex fileMutex; // Mutex for file access
    std::string out_filename;
    double energy; // incoming neutrino energy. Read from settings.yaml
    std::string tau_out_filename;
    bool write_tau_info = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}

#endif