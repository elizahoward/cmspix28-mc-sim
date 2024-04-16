/*

  Authors: 
  Anthony Badea

  Description:
  Script to generate minbias pythia samples
  Based on pythia main08 example https://pythia.org/latest-manual/examples/main08.html

*/

// c, pythia, root includes
#include <iostream>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "Pythia8Plugins/HepMC3.h"
#include "TH1D.h"

// declare pythia8 namespace
using namespace Pythia8;

// // print a message to the console
// void msg(string m){
//   printf("\r%s",m.c_str());                               
//   std::cout << std::endl;
// }

// // progress bar for the people taken from alex tuna and ann wang
// void pbftp(double time_diff, int nprocessed, int ntotal){
//   if(nprocessed%10 == 0){
//     double rate      = (double)(nprocessed+1)/time_diff;
//     std::cout << "\r > " << nprocessed << " / " << ntotal
// 	      << " | "   << std::fixed << std::setprecision(1) << 100*(double)(nprocessed)/(double)(ntotal) << "%"
// 	      << " | "   << std::fixed << std::setprecision(1) << rate << "Hz"
// 	      << " | "   << std::fixed << std::setprecision(1) << time_diff/60 << "m elapsed"
// 	      << " | "   << std::fixed << std::setprecision(1) << (double)(ntotal-nprocessed)/(rate*60) << "m remaining"
// 	      << std::flush;
//   }
//   // add new line at end of events
//   if (nprocessed+1 == ntotal){
//     msg("");
//   }
// }


int main(int argc, char* argv[]) {

  // check user is inputting correct parameters
  if(argc != 5){
    std::cout << "Usage: ./minbias.exe <outFileName> <maxEvents> <pTHatMin> <pTHatMax>" << std::endl;
    return 1;
  }
  
  // read in user parameters
  std::string outFileName = argv[1];
  int maxEvents = std::stoi(argv[2]);
  double pTHatMin = std::stod(argv[3]);
  double pTHatMax = std::stod(argv[4]);

  // other settings
  double global_pTHatMin = 0.25; // GeV, min pT bin to enable SoftQCD:nonDiffractive
  
  // Generator.
  Pythia pythia;
  
  // Shorthand for some public members of pythia (also static ones).
  Settings& settings = pythia.settings;
  const Info& info = pythia.info;
  // number of bins
  // int nBin = 1;
  // set up pT bins
  // double pTlimitTwo[11] = {0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5};
  // double pTlimitTwo[11] = {0., 0.25};
  
  // // Create root file and tree
  // TFile* f = new TFile((outFileName + ".root").c_str(),"RECREATE");
  // TTree* t = new TTree("t","t");

  // // Create a histogram for mc event weights
  // TH1D* MetaData = new TH1D("MetaData", "", 2, 0, 1);
  
  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC3::Pythia8ToHepMC3 toHepMC;
  // Specify file where HepMC events will be stored.
  HepMC3::WriterAscii ascii_io((outFileName + ".hepmc").c_str());

  // // initialize branches variables
  // int nParticles;
  // double weight;
  // std::vector<int> *id = 0;
  // std::vector<int> *status = 0;
  // std::vector<int> *statusHepMC = 0;
  // std::vector<double> *mass = 0;
  // std::vector<double> *pt = 0;
  // std::vector<double> *eta = 0;
  // std::vector<double> *phi = 0;
  // std::vector<double> *tau = 0;
  // std::vector<double> *xDec = 0;
  // std::vector<double> *yDec = 0;
  // std::vector<double> *zDec = 0;
  // std::vector<double> *tDec = 0;

  // // initialize branches
  // t->Branch("nParticles", &nParticles);
  // t->Branch("weight", &weight);
  // t->Branch("id", &id);
  // t->Branch("status", &status);
  // t->Branch("statusHepMC", &statusHepMC);
  // t->Branch("mass", &mass);
  // t->Branch("pt", &pt);
  // t->Branch("eta", &eta);
  // t->Branch("phi", &phi);
  // t->Branch("tau", &tau);
  // t->Branch("xDec", &xDec);
  // t->Branch("yDec", &yDec);
  // t->Branch("zDec", &zDec);
  // t->Branch("tDec", &tDec);

  // // time keeper for progress bar
  // std::chrono::time_point<std::chrono::system_clock> time_start;
  // std::chrono::duration<double> elapsed_seconds;
  // time_start = std::chrono::system_clock::now();

  if (pTHatMin < global_pTHatMin){
    pythia.readString("HardQCD:all = off");
    pythia.readString("SoftQCD:nonDiffractive = on");
  }
  else{
    pythia.readString("HardQCD:all = on");
    pythia.readString("SoftQCD:nonDiffractive = off");
  }
  settings.parm("PhaseSpace:pTHatMin", pTHatMin);
  settings.parm("PhaseSpace:pTHatMax", pTHatMax);
  
  // Hard processes in one step, but pT-weighted.
  // if (pTHatMin > 0){
  //   pythia.readString("PhaseSpace:bias2Selection = on");
  //   pythia.readString("PhaseSpace:bias2SelectionPow = 4.");
  //   pythia.readString("PhaseSpace:bias2SelectionRef = 20.");
  // }

   // Initialize for LHC at 13.6 TeV
  pythia.readString("Beams:eCM = 13600.");
  pythia.init();

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;
  
  for (int iE = 0; iE < maxEvents; ++iE) {
    
    if(!pythia.next()) continue;
    
    // Construct new empty HepMC event and fill it.
    // Default units are ( HepMC3::Units::GEV, HepMC3::Units::MM)
    // but can be changed in the GenEvent constructor.
    HepMC3::GenEvent hepmcevt;
    toHepMC.fill_next_event( pythia, &hepmcevt );
    
    // Write the HepMC event to file.
    ascii_io.write_event(hepmcevt);
  }
  
  // // Setting
  // // Loop over number of bins, i.e. number of subruns.
  // for (int iBin = 0; iBin < nBin; ++iBin) {

  //   // Normally HardQCD, but in two cases nonDiffractive.
  //   // Need MPI on in nonDiffractive to get first interaction, but not else.
  //   if(iBin == 0){
  //     pythia.readString("HardQCD:all = off");
  //     pythia.readString("SoftQCD:nonDiffractive = on");
  //   }
  //   else{
  //     pythia.readString("HardQCD:all = on");
  //     pythia.readString("SoftQCD:nonDiffractive = off");
  //   }

  //   // Mode 5: hardcoded here. Use settings.parm for non-string input.
  //   // Hard processes in one step, but pT-weighted.
  //   settings.parm("PhaseSpace:pTHatMin", pTlimitTwo[iBin]);
  //   settings.parm("PhaseSpace:pTHatMax", pTlimitTwo[iBin + 1]);
  //   if (iBin == 1) {
  //     pythia.readString("PhaseSpace:bias2Selection = on");
  //     pythia.readString("PhaseSpace:bias2SelectionPow = 4.");
  //     pythia.readString("PhaseSpace:bias2SelectionRef = 20.");
  //   }

  //   // Initialize for LHC at 13.6 TeV
  //   pythia.readString("Beams:eCM = 13600.");
  //   pythia.init();
    
  //   // Shorthand for the event record in pythia.
  //   Event& event = pythia.event;
    
  //   for (int iE = 0; iE < maxEvents; ++iE) {
      
  //     if(!pythia.next()) continue;
      
  //     // Construct new empty HepMC event and fill it.
  //     // Default units are ( HepMC3::Units::GEV, HepMC3::Units::MM)
  //     // but can be changed in the GenEvent constructor.
  //     HepMC3::GenEvent hepmcevt;
  //     toHepMC.fill_next_event( pythia, &hepmcevt );

  //     // Write the HepMC event to file.
  //     ascii_io.write_event(hepmcevt);

  //     // progress bar
  //     // elapsed_seconds = (std::chrono::system_clock::now() - time_start);
  //     // pbftp(elapsed_seconds.count(), iE, maxEvents);

  //     // clear for new event
  //     // id->clear();
  //     // status->clear();
  //     // statusHepMC->clear();
  //     // mass->clear();
  //     // pt->clear();
  //     // eta->clear();
  //     // phi->clear();
  //     // tau->clear();
  //     // xDec->clear();
  //     // yDec->clear();
  //     // zDec->clear();
  //     // tDec->clear();
      
  //     // get event level information
  //     // nParticles = event.size();
  //     // weight = info.weight();

  //     // loop over the particles. available properties listed here https://pythia.org/latest-manual/ParticleProperties.html
  //     // for(int iP=0; iP<nParticles; iP++){

  //     // 	// save particle information
  //     // 	id->push_back(pythia.event[iP].id());
  //     // 	status->push_back(pythia.event[iP].status());
  //     // 	statusHepMC->push_back(pythia.event[iP].statusHepMC());
  //     // 	mass->push_back(pythia.event[iP].m());
  //     // 	pt->push_back(pythia.event[iP].pT());
  //     // 	eta->push_back(pythia.event[iP].eta());
  //     // 	phi->push_back(pythia.event[iP].phi());
  //     // 	tau->push_back(pythia.event[iP].tau());
  //     // 	xDec->push_back(pythia.event[iP].xDec());
  //     // 	yDec->push_back(pythia.event[iP].yDec());
  //     // 	zDec->push_back(pythia.event[iP].zDec());
  //     // 	tDec->push_back(pythia.event[iP].tDec());
  //     // }

  //     // fill tree on each particle loop
  //     // t->Fill();
  //   }
    
  //   // pythia.stat();
  // }

  // // save metadata to root (this is automatic in hepmc file)
  // MetaData->SetBinContent(1, info.sigmaGen());
  // MetaData->GetXaxis()->SetBinLabel(1, "sigmaGen");
  // MetaData->SetBinContent(2, info.weightSum());
  // MetaData->GetXaxis()->SetBinLabel(2, "weightsSum");
  // MetaData->Write();
  
  // // write and cleanup
  // t->Write();  
  // delete t;

  // // close file and cleanup
  // f->Close();
  // delete f;

  return 0;
}
