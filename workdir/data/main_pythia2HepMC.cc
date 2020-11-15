// main30.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Steve Mrenna.

// Example how to create a copy of the event record, where the original one
// is translated to another format, to meet various analysis needs.
// In this specific case the idea is to set up the history information
// of the underlying hard process to be close to the PYTHIA 6 structure.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
//#include "Pythia8/LHEF3.h"
#include "Pythia8/LesHouches.h"

// Make root histograms
//#include "TH1.h"

// ROOT, for saving file.
//#include "TFile.h"

using namespace Pythia8;

//--------------------------------------------------------------------------

void translate(Event&, Event&);

//--------------------------------------------------------------------------

int run(const string & infile, int nevents, const string & cfgfile, const string & outputfile)
{

  //Set output file names
  std::srand(500);
  string outname;
  if (outputfile == ""){
	  size_t lastindex = infile.find_last_of(".");
	  outname = infile.substr(0, lastindex)+".hep";
  }
  else{outname = outputfile;}

  // Generator. Shorthand for the event.
  Pythia pythia("",false); //Set printBanner to false
  Event& event = pythia.event;
  Event& hard = pythia.process;
  pythia.readFile( cfgfile );
  if ( infile.find(".slha") != std::string::npos ){
    cout << "Using SLHA file as input" << endl;
    pythia.readString("SLHA:file = " + infile);
    if ( nevents < 0) {
    	nevents = 100;
    	cout << "Negative number of events for SLHA input. Setting nevents to " << nevents << endl;
    }
  }
  else{
    cout << "Using LHE file as input" << endl;
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = " + infile);
  }


// trying to write hepmc
// Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;

// Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(outname, std::ios::out);

  // Create an LHEF object that can access relevant information in pythia. process is the hard event.

  // LHEF3FromPythia8 myLHEF3(&pythia.process, &pythia.settings, &pythia.info, &pythia.particleData);

  // Open a file on which LHEF events should be stored, and write header.
  // myLHEF3.openLHEF(lhefilename);

  // Initialize.
  pythia.init();

  // Store initialization info in the LHAup object.
  // myLHEF3.setInit();

  // Begin event loop.
  int iAbort = 0;
  int iEvent = 0;
  while (iEvent < nevents or nevents < 0){

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;

      // Generate events. Quit if failure.
      if (!pythia.next()) {
        if (++iAbort < 10) continue;
        cout << " Event generation aborted prematurely, owing to error!\n";
        break;
      }

      ++iEvent;


      // Store event info in the LHAup object.
      // myLHEF3.setEvent();

      // Construct new empty HepMC event and fill it.
      // Units will be as chosen for HepMC build; but can be changed
      // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;
    //
    // if (iEvent < nShow) hard.list();

  // End of event loop.
  }

  // Final statistics.
  pythia.stat();
  // Update the cross section info based on Monte Carlo integration during run.
  // myLHEF3.closeLHEF(true);
  // Done.
  return 0;
}

//--------------------------------------------------------------------------



void help( const char * name )
{
	  cout << "syntax: " << name << " [-h] [-f <input file>] [-o <output file>] [-n <number of events>] [-c <pythia cfg file>]" << endl;
	  cout << "        -f <input file>:  pythia input LHE or SLHA file" << endl;
	  cout << "        -c <pythia config file>:  pythia config file [pythia8_stau.cfg]" << endl;
	  cout << "        -o <output file>:  output filename for naming the output file and histograms [<input file>.out]" << endl;
	  cout << "        -n <number of events>:  Number of events to be generated [100]. If n < 0, it will run over all events in the LHE file" << endl;
  exit( 0 );
};

int main( int argc, const char * argv[] ) {
  int nevents = -1;
  double width = 0.;
  string cfgfile = "pythia8.cfg";
  string outfile = "";
  string infile = "";
  for ( int i=1; i!=argc ; ++i )
  {
    string s = argv[i];
    if ( s== "-h" )
    {
      help ( argv[0] );
    }

    if ( s== "-c" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      cfgfile = argv[i+1];
      i++;
      continue;
    }


    if ( s== "-n" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      nevents = atoi(argv[i+1]);
      i++;
      continue;
    }

    if ( s== "-o" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      outfile = argv[i+1];
      i++;
      continue;
    }



    if ( s== "-f" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      infile = argv[i+1];
      i++;
      continue;
    }



    cout << "Error. Argument " << argv[i] << " unknown." << endl;
    help ( argv[0] );
  };

  int r = run(infile, nevents, cfgfile, outfile);

  return 0;
}
