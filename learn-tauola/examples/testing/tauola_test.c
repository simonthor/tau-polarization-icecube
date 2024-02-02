/**
 * Example of use of tauola C++ interface. Pythia events are
 * generated with a stable tau. Taus are subsequently decay via
 * tauola.
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */

#include "Tauola/Log.h"
#include "Tauola/Plots.h"
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"

//pythia header files
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

//MC-TESTER header files
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"

using namespace std;
using namespace Pythia8; 
using namespace Tauolapp;

unsigned int NumberOfEvents = 10000;
unsigned int EventsToCheck  = 20;

// elementary test of HepMC typically executed before
// detector simulation based on http://home.fnal.gov/~mrenna/HCPSS/HCPSShepmc.html
// similar test was performed in Fortran
// we perform it before and after Tauola (for the first several events)
void checkMomentumConservationInEvent(HepMC::GenEvent *evt)
{
	//cout<<"List of stable particles: "<<endl;

	double px=0.0,py=0.0,pz=0.0,e=0.0;
	
	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
	      p != evt->particles_end(); ++p )
	{
		if( (*p)->status() == 1 )
		{
			HepMC::FourVector m = (*p)->momentum();
			px+=m.px();
			py+=m.py();
			pz+=m.pz();
			e +=m.e();
			//(*p)->print();
		}
	}
	cout.precision(6);
	cout.setf(ios_base::scientific, ios_base::floatfield);
	cout<<endl<<"Vector Sum: "<<px<<" "<<py<<" "<<pz<<" "<<e<<endl;
}

int main(int argc,char **argv){

  // Program needs at least 4 parameters
  if(argc<5)
  {
    cout<<endl<<"Usage: "<<argv[0]<<" <pythia_conf> <pythia_mode> <no_events> <tauola_mode> <mixing_angle>"<<endl;
    cout<<endl<<"   eg. "<<argv[0]<<" pythia_H.conf 0 10000 4 0.7853"<<endl;
    cout<<endl;
    return -1;
  }

  Log::SummaryAtExit();

  // Initialisation of pythia
  Pythia pythia;
  Event& event = pythia.event;
  
  // Pythia8 HepMC interface depends on Pythia8 version
  HepMC::Pythia8ToHepMC ToHepMC;
  
  // Initial pythia configuration
  pythia.particleData.readString("15:mayDecay = off");

  /*
    Read input parameters from console. List of parameters:
    1. Pythia configuration filename
    2. Are we using pp collisions? (If not - e+ e- collisions)
    3. Number of events
    4. Tauola decay mode (refer to documentation)
    5. Higgs scalar-pseudoscalar mixing angle

    Example where all input parameters are used:

    ./taumain_pythia_example.exe pythia_H.conf 0 100000 4 0.7853
      - use pythia_H.conf
      - initialize using e+ e- collisions
      - generate 100 000 events
      - fix TAUOLA decay to channel 4 (RHO_MODE)
      - Higgs scalar-pseudoscalar mixing angle set to 0.7853
  */

  // 1. Load pythia configuration file (argv[1], from console)
  if(argc>1) pythia.readFile(argv[1]);

  // 2. Initialize pythia to pp or e+e- collisions (argv[2], from console)
  if(atoi(argv[2])==1) // p_bar  p_bar  collisions
  {
    pythia.readString("Beams:idA = -2212");
    pythia.readString("Beams:idB = -2212");
    pythia.readString("Beams:eCM =  14000.0");
  }
  else // e+ e- collisions
  {
    pythia.readString("Beams:idA =  11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("Beams:eCM =  500");
  }
  
  pythia.init();

  // 3. Get number of events (argv[3], from console)
  if(argc>3) NumberOfEvents=atoi(argv[3]);

  // 4. Set Tauola decay mode (argv[4], from console)
  if(argc>4)
  {
    Tauola::setSameParticleDecayMode(atoi(argv[4]));
    Tauola::setOppositeParticleDecayMode(atoi(argv[4]));
  }

  // 5. Set Higgs scalar-pseudoscalar mixing angle (argv[5], from console)
  if(argc>5)
  {
    Tauola::setHiggsScalarPseudoscalarMixingAngle(atof(argv[5]));
    Tauola::setHiggsScalarPseudoscalarPDG(25);
  }

  Tauola::setRadiation(false); // turn off radiation in leptionic decays
  Tauola::initialize();

  MC_Initialize();

  // Begin event loop
  for(unsigned int iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
  {
    cout.setf(ios_base::fixed, ios_base::floatfield);
    if (iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<endl;
    if (!pythia.next()) continue;

    // Convert event record to HepMC
    HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();

    // Conversion needed if HepMC use different momentum units
    // than Pythia. However, requires HepMC 2.04 or higher.
    HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);

    ToHepMC.fill_next_event(event, HepMCEvt);

    if(iEvent<EventsToCheck)
    {
      cout<<"                                          "<<endl;
      cout<<"Momentum conservation chceck BEFORE/AFTER Tauola"<<endl;
      checkMomentumConservationInEvent(HepMCEvt);
    }

    // Run TAUOLA on the event
    TauolaHepMCEvent * t_event = new TauolaHepMCEvent(HepMCEvt);

    // We do not let Pythia decay taus, so we don't have to undecay them
    //t_event->undecayTaus();
    t_event->decayTaus();
    delete t_event;

    if(iEvent<EventsToCheck)
    {
      checkMomentumConservationInEvent(HepMCEvt);
    }

    // Run MC-TESTER on the event
    HepMCEvent temp_event(*HepMCEvt,false);
    MC_Analyze(&temp_event);

    // Clean up HepMC event
    delete HepMCEvt;
  }
  pythia.stat();
  Tauola::summary();
  MC_Finalize();
}

