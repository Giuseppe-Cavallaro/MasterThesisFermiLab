// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a TTree for studying the matching
// between flashes and events
//

#ifndef FlashMatchAnaEdep_H
#define FlashMatchAnaEdep_H 1

// ROOT includes
#include "TH1.h"
#include "TEfficiency.h"
#include "TTree.h"

// C++ includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
//DUNE includes
#include "dune/OpticalDetector/OpFlashSort.h"
#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaUtilsBase.h"
#include "dune/AnaUtils/DUNEAnaHitUtils.h"

#include "/dune/app/users/gcavalla/dunetpc19/srcs/dunetpc/dune/AnaUtils/DUNEAnaEventUtils.h"
#include "/dune/app/users/gcavalla/dunetpc19/srcs/dunetpc/dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "/dune/app/users/gcavalla/dunetpc19/srcs/dunetpc/dune/AnaUtils/DUNEAnaTrackUtils.h"
#include "/dune/app/users/gcavalla/dunetpc19/srcs/dunetpc/dune/AnaUtils/DUNEAnaShowerUtils.h"
#include "/dune/app/users/gcavalla/dunetpc19/srcs/dunetpc/dune/AnaUtils/DUNEAnaHitUtils.h"



namespace opdet {

  class FlashMatchAnaEdep : public art::EDAnalyzer{
  public:

    // Standard constructor and destructor for an ART module.
    FlashMatchAnaEdep(const fhicl::ParameterSet&);
    virtual ~FlashMatchAnaEdep();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The analyzer routine, called once per event.
    void analyze (const art::Event&);

    void endJob();

  private:

    art::ServiceHandle< geo::Geometry > geom;

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fOpFlashModuleLabel;       // Input tag for OpFlash collection
    std::string fOpHitModuleLabel;         // Input tag for OpHit collection
    std::string fSignalLabel;              // Input tag for the signal generator label
    std::string fGeantLabel;               // Input tag for GEANT
    std::string fSimLabel;                 // Input tag for Sim
    std::string fHitsLabel;                //Input tag for Charge hits
    std::string fSimulationProducerLabel;  // Input tag for Sim
    
    //pandora
     std::string fTrackLabel;
     //std::string fShowerLabel;
     std::string fPFParticleLabel;   
     std::string fWireLabel;
            std::string fHitLabel;
            
            std::string fShowerLabel;
            std::string fTrackToHitLabel;
            std::string fShowerToHitLabel;
            std::string fHitToSpacePointLabel;   
            fhicl::ParameterSet fParam;      
    
    
    
    bool fBeam;                            // Simulated events are beam neutrinos

    TTree * fFlashMatchTree;
//    TTree * fLargestFlashTree;
//    TTree * fSelectedFlashTree;

    /* TEfficiency * fRecoEfficiencyVsE;
    TEfficiency * fRecoEfficiencyVsX;
    TEfficiency * fRecoEfficiencyVsXandE;
    TEfficiency * fLargestEfficiencyVsE;
    TEfficiency * fLargestEfficiencyVsX;
    TEfficiency * fLargestEfficiencyVsXandE;
    TEfficiency * fSelectedEfficiencyVsE;
    TEfficiency * fSelectedEfficiencyVsX;
    TEfficiency * fSelectedEfficiencyVsXandE;
    */
    
    // Parameters from the fhicl
    int   fNBinsE;
    float fLowE;
    float fHighE;
    int   fNBinsX;
    float fLowX;
    float fHighX;
    float fDistanceCut;

    Int_t fEventID;
    Int_t fPandoraContained50=0;
   

    Float_t fTrueX;
    Float_t fTrueY;
    Float_t fTrueZ;
    Float_t fTrueT;
    Float_t fDetectedT;
    Float_t fTrueE;
    Int_t   fTruePDG;
    Int_t   fTrueCCNC;
    Float_t fRecoX;
    
//mappa di luce

    Float_t dist1;
    Float_t dist2;
    Float_t sentheta;
    Float_t areaeff;
   // std::vector< Int_t >   fPEperOpDetcorr; 
   //corretti
    
//********Ion and scint
    Float_t fEdep;
    Float_t fTotEdep;
    Float_t fLdepSim; //track length from simulation
    Float_t fLdepGeom; // track lengt from geometrical calculation
    Int_t fTotEion; //number of ionizing electrons
    Int_t fTotGammaScint; //number of scintillation gammas
    Float_t fEiondE; //number of ionizing electrons divided by deposited energy
    Float_t fGammaScintdE; //number of scintillation gammas divided by deposited energy
   
    
//*********Charge from reco hits
    Float_t fTotalCharge; //Integral under the calibrated signal waveform of the hit, in tick x ADC units
    Float_t fTotalChargeCorr;//Total charge corrected for the electron lifetime
    Float_t fTotalChargeFinal;//Total charge corrected for the electron lifetime
    Float_t fADCSum; //The sum of calibrated ADC counts of the hit (0. by default)
    Float_t fPeakAmplitude; //The estimated amplitude of the hit at its peak, in ADC units
    Float_t fFirstHitTime; //First Time recorded among hits, Time of the signal peak,in us
    Float_t fLastHitTime; //Last Time recorded among hits, Time of the signal peak, in us
    Float_t fMeanHitTime; //Last Time recorded among hits, Time of the signal peak, in us
    Float_t fHitDist; //First hit time * 1.6 mm/us --> theor distance travelled by charge (cm)
    
    std::vector< Float_t > fTruePxallpart;
    std::vector< Float_t > fTruePyallpart;
    std::vector< Float_t > fTruePzallpart;
    std::vector< Float_t > fTrueEallpart;
    std::vector< Int_t >   fTrueAllPDG;

    Int_t fNFlashes;
    std::vector< Int_t >   fFlashIDVector;
    std::vector< Float_t > fYCenterVector;
    std::vector< Float_t > fZCenterVector;
    std::vector< Float_t > fYWidthVector;
    std::vector< Float_t > fZWidthVector;
    std::vector< Float_t > fTimeVector;
    std::vector< Float_t > fRecoXVector;
    std::vector< Float_t > fTimeWidthVector;
    std::vector< Float_t > fTimeDiffVector;
    std::vector< Float_t > fTotalPEVector;
    std::vector< Float_t > fPurityVector;
    std::vector< Float_t > fDistanceVector;
    Int_t fNOpDets;
    std::vector<Int_t> fNHitOpDetVector;
    std::vector< Float_t > fOpHitPeakTime; //Time of OpHit in us (?)
    std::vector< Float_t > fOpHitFastToTotal;
    //*******Ion and Scint  
    std::vector< Float_t > fEnergyDepositionVector;
    std::vector< Float_t > fPointX;
    std::vector< Float_t > fPointY;
    std::vector< Float_t > fPointZ;    
    std::vector< Int_t >   fGammaScint; //scint gamma each step
    std::vector< Int_t >   fPEperOpDet; //scint gamma each step
    std::vector< Float_t > fStepLength;
    std::vector< Float_t > fStepLCumVector;
    std::vector< Float_t > fStepEdepCumVector;
    
     std::vector< Float_t >  fTrackPointX;
     std::vector< Float_t >  fTrackPointY;
     std::vector< Float_t >  fTrackPointZ;
     
     std::vector< Float_t > fTrackMiddlePointX_1;
     std::vector< Float_t > fTrackMiddlePointY_1;
     std::vector< Float_t > fTrackMiddlePointZ_1;
     
     std::vector< Float_t > fTrackMiddlePointX;
     std::vector< Float_t > fTrackMiddlePointY;
     std::vector< Float_t > fTrackMiddlePointZ;

    std::vector< Int_t > fHitMultiplicity; //How many hits could this one be shared with. Index of this hit among the Multiplicity() hits in the signal window
    std::vector< Float_t > fHitPeakTime; //Time of the signal peak, converted in us
    //std::vector< Float_t > fHitPeakTimeTicks; //Time of the signal peak, in tick units.
    
    Int_t    fFlashID;
    Float_t  fYCenter;
    Float_t  fZCenter;
    Float_t  fYWidth;
    Float_t  fZWidth;
    Float_t  fTime;
    Float_t  fTimeWidth;
    Float_t  fTimeDiff;
    Float_t  fTotalPE;    
    Float_t  fSumPE;//Sum of PE in all
    Float_t  fPurity;
    Float_t  fDistance;
    Int_t    fNHitOpDets;
    std::vector< Float_t > fPEsPerOpDetVector;

    // For counting waveforms

    std::string fOpDetWaveformLabel;
    float fBaseline;
    float fPE;
    TTree * fCountTree;
    Int_t fnwaveforms1pe;
    Int_t fnwaveforms2pe;
    Int_t fnwaveforms3pe;
    
   //Convert charge to energy
   float Wion = 23.6e-6; //MeV/electron
   float Recombination = 0.64; //from Protodune note
   
   //map
    std::vector< Float_t > fMiddledepx;
    std::vector< Float_t > fMiddledepy;
    std::vector< Float_t > fMiddledepz;  
    
    std::vector< Float_t > fMiddlex;
    std::vector< Float_t > fMiddley;
    std::vector< Float_t > fMiddlez;
    
    
    std::vector< Int_t > fIdDep;
    std::vector< Int_t > fTrackDep;
    std::vector< Float_t > fAreaeff;
   // std::vector< Int_t > fPhotonperDep;  
    std::vector< Float_t > fConteggi;
    std::vector< Float_t > fConteggiTot;
    //std::vector< Int_t >fCountperOpdet;
    std::vector< Int_t > fIdOpdet;
    std::vector< Float_t > fCountperOpdet;
    std::vector< Float_t >  fConteggiTotDep;
    std::vector< Float_t > fConteggiTotDepPE;
    //std::vector< Int_t > fTotalConteggiCorr;
    std::vector< Float_t > fTuttiPE;
    std::vector< Float_t > fCountperOpdetPE;
    std::vector< Float_t > fConteggiTotPE;
    
   std::vector< Float_t > fGammaStep;
   std::vector< Float_t > fPEStep;
  };

}

#endif // FlashMatchAnaEdep_H

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  FlashMatchAnaEdep::FlashMatchAnaEdep(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fOpFlashModuleLabel = pset.get<std::string>("OpFlashModuleLabel");
    fOpHitModuleLabel   = pset.get<std::string>("OpHitModuleLabel");
    fHitsLabel          = pset.get<std::string>("HitsLabel");
    fSignalLabel        = pset.get<std::string>("SignalLabel");
    fGeantLabel         = pset.get<std::string>("GeantLabel");
    fSimLabel           = pset.get<std::string>("SimLabel");
  
    fWireLabel=pset.get<std::string>("WireLabel");
    fHitLabel=pset.get<std::string>("HitLabel");
   
    fShowerLabel=pset.get<std::string>("ShowerLabel");
    fTrackToHitLabel=pset.get<std::string>("TrackToHitLabel");
    fShowerToHitLabel=pset.get<std::string>("ShowerToHitLabel");
    fHitToSpacePointLabel=pset.get<std::string>("HitToSpacePointLabel");
    fParam= pset.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg");
    
    fSimulationProducerLabel = pset.get<std::string>("SimulationLabel");
    
    //roba extra pandora
      fPFParticleLabel = pset.get<std::string>("PFParticleLabel");
      fTrackLabel = pset.get<std::string>("TrackLabel");
      //fShowerLabel = pset.get<std::string>("ShowerLabel");
    
    fBeam               = pset.get<bool>("Beam");
    fNBinsE             = pset.get<int>("NBinsE");
    fLowE               = pset.get<float>("LowE");
    fHighE              = pset.get<float>("HighE");
    fNBinsX             = pset.get<int>("NBinsX");
    fLowX               = pset.get<float>("LowX");
    fHighX              = pset.get<float>("HighX");
    fDistanceCut        = pset.get<float>("DistanceCut");

    fOpDetWaveformLabel = pset.get<std::string>("OpDetWaveformLabel","");
    fBaseline           = pset.get<float>("Baseline", 1500.);
    fPE                 = pset.get<float>("PE", 18.);

    art::ServiceHandle< art::TFileService > tfs;

    fFlashMatchTree = tfs->make<TTree>("FlashMatchTree","FlashMatchTree");
    fFlashMatchTree->Branch("EventID",                     &fEventID,   "EventID/I");
    //fFlashMatchTree->Branch("TrueX",                       &fTrueX,     "TrueX/F");
    //fFlashMatchTree->Branch("TrueY",                       &fTrueY,     "TrueY/F");
    //fFlashMatchTree->Branch("TrueZ",                       &fTrueZ,     "TrueZ/F");
    //fFlashMatchTree->Branch("TrueT",                       &fTrueT,     "TrueT/F");
    //fFlashMatchTree->Branch("DetectedT",                   &fDetectedT, "DetectedT/F");
    //fFlashMatchTree->Branch("TrueE",                       &fTrueE,     "TrueE/F");
    //fFlashMatchTree->Branch("TruePDG",                     &fTruePDG,   "TruePDG/I");
    //fFlashMatchTree->Branch("TrueCCNC",                    &fTrueCCNC,  "TrueCCNC/I");
    //fFlashMatchTree->Branch("NFlashes",                    &fNFlashes,  "NFlashes/I");
    //fFlashMatchTree->Branch("FlashIDVector",               &fFlashIDVector);
    //fFlashMatchTree->Branch("YCenterVector",               &fYCenterVector);
    //fFlashMatchTree->Branch("ZCenterVector",               &fZCenterVector);
    //fFlashMatchTree->Branch("YWidthVector",                &fYWidthVector);
    //fFlashMatchTree->Branch("ZWidthVector",                &fZWidthVector);
    //fFlashMatchTree->Branch("TimeVector",                  &fTimeVector);
    //fFlashMatchTree->Branch("TimeWidthVector",             &fTimeWidthVector);
    //fFlashMatchTree->Branch("TimeDiffVector",              &fTimeDiffVector);
    fFlashMatchTree->Branch("TotalPEVector",               &fTotalPEVector);
    fFlashMatchTree->Branch("SumPE",                       &fSumPE,      "SumPE/F");
    //fFlashMatchTree->Branch("NOpDets",                     &fNOpDets, "NOpDets/I");
    fFlashMatchTree->Branch("NHitOpDetVector",             &fNHitOpDetVector);
    fFlashMatchTree->Branch("OpHitPeakTime",               &fOpHitPeakTime);
    //fFlashMatchTree->Branch("OpHitFastToTotal",            &fOpHitFastToTotal);
    //fFlashMatchTree->Branch("Purity",                      &fPurityVector);
    //fFlashMatchTree->Branch("Distance",                    &fDistanceVector);
    //fFlashMatchTree->Branch("RecoXVector",                 &fRecoXVector);
    fFlashMatchTree->Branch("TruePxallpart",               &fTruePxallpart);
    fFlashMatchTree->Branch("TruePyallpart",               &fTruePyallpart);
    fFlashMatchTree->Branch("TruePzallpart",               &fTruePzallpart);
    fFlashMatchTree->Branch("TrueEallpart",                &fTrueEallpart);
    fFlashMatchTree->Branch("TrueAllPDG",                  &fTrueAllPDG);
    fFlashMatchTree->Branch("PointX",     		     &fPointX);
    fFlashMatchTree->Branch("PointY",     		     &fPointY);
    fFlashMatchTree->Branch("PointZ",     		     &fPointZ);
    fFlashMatchTree->Branch("GammaScint",     	     &fGammaScint);
    fFlashMatchTree->Branch("PEperOpDet",     	     &fPEperOpDet);
    //fFlashMatchTree->Branch("PEperOpDetcorr",     	     &fPEperOpDetcorr);
    fFlashMatchTree->Branch("EnergyDepositionVector",      &fEnergyDepositionVector);
    fFlashMatchTree->Branch("fStepLCumVector",             &fStepLCumVector);
    fFlashMatchTree->Branch("fStepEdepCumVector",          &fStepEdepCumVector);
    fFlashMatchTree->Branch("TotEdep",                     &fTotEdep, "TotEdep/F");
    fFlashMatchTree->Branch("StepLengthEdep",              &fStepLength);
    fFlashMatchTree->Branch("LdepSim",                     &fLdepSim,   "LdepSim/F");
    fFlashMatchTree->Branch("LdepGeom",                    &fLdepGeom,"LdepGeom/F");
    fFlashMatchTree->Branch("TotEion",                     &fTotEion, "TotEion/I");
    fFlashMatchTree->Branch("TotGammaScint",               &fTotGammaScint,"TotGammaScint/I");
    fFlashMatchTree->Branch("TotalCharge",                 &fTotalCharge, "TotalCharge/F");
    fFlashMatchTree->Branch("TotalChargeCorr",            &fTotalChargeCorr, "TotalChargeCorr/F");
    fFlashMatchTree->Branch("TotalChargeFinal",            &fTotalChargeFinal, "TotalChargeFinal/F");
    fFlashMatchTree->Branch("ADCSum",                      &fADCSum,      "ADCSum/F");
    fFlashMatchTree->Branch("PeakAmplitude",               &fPeakAmplitude,"PeakAmplitude/F");
    fFlashMatchTree->Branch("HitMultiplicity",             &fHitMultiplicity);
    fFlashMatchTree->Branch("FirstHitTime",                &fFirstHitTime,  "FirstHitTime/F");
    fFlashMatchTree->Branch("LastHitTime",                 &fLastHitTime,  "LastHitTime/F");
    fFlashMatchTree->Branch("MeanHitTime",                 &fMeanHitTime, "MeanHitTime/F");
    fFlashMatchTree->Branch("HitPeakTime",                 &fHitPeakTime);
    //fFlashMatchTree->Branch("HitPeakTimeTicks",            &fHitPeakTimeTicks);
    fFlashMatchTree->Branch("HitDist",                     &fHitDist,  "HitDist/F");
    
    //mappadiluce
    fFlashMatchTree->Branch("Middledepx",               &fMiddledepx);
    fFlashMatchTree->Branch("Middledepy",               &fMiddledepy);
    fFlashMatchTree->Branch("Middledepz",               &fMiddledepz);
    fFlashMatchTree->Branch("PandoraContained50",                     &fPandoraContained50);
    fFlashMatchTree->Branch("IdDep",                &fIdDep);
    fFlashMatchTree->Branch("IdOpdet",                &fIdOpdet);
    fFlashMatchTree->Branch("TrackDep",                  &fTrackDep);
    fFlashMatchTree->Branch("Areaeff",                  &fAreaeff);
    //fFlashMatchTree->Branch("PhotonperDep",                  &fPhotonperDep);
  
    fFlashMatchTree->Branch("CountperOpdet",                  &fCountperOpdet);
    fFlashMatchTree->Branch("Conteggi",                  &fConteggi);
    fFlashMatchTree->Branch("ConteggiTot",                  &fConteggiTot);
    fFlashMatchTree->Branch("ConteggiTotDep",                  &fConteggiTotDep);
    fFlashMatchTree->Branch("ConteggiTotDepPE",                  &fConteggiTotDepPE);
    fFlashMatchTree->Branch("TuttiPE",                  &fTuttiPE);
    
    fFlashMatchTree->Branch("ConteggiTotPE",                  &fConteggiTotPE);
    fFlashMatchTree->Branch("CountperOpdetPE",                  &fCountperOpdetPE);
    //traccia
    fFlashMatchTree->Branch("GammaStep",                  &fGammaStep);
    fFlashMatchTree->Branch("PEStep",                  &fPEStep);
    fFlashMatchTree->Branch("Middlex",               &fMiddlex);
    fFlashMatchTree->Branch("Middley",               &fMiddley);
    fFlashMatchTree->Branch("Middlez",               &fMiddlez);
    fFlashMatchTree->Branch("TrackPointX",               &fTrackPointX);
    fFlashMatchTree->Branch("TrackPointY",               &fTrackPointY);
    fFlashMatchTree->Branch("TrackPointZ",               &fTrackPointZ);
    
    fFlashMatchTree->Branch("TrackMiddlePointX",               &fTrackMiddlePointX);
    fFlashMatchTree->Branch("TrackMiddlePointY",               &fTrackMiddlePointY);
    fFlashMatchTree->Branch("TrackMiddlePointZ",               &fTrackMiddlePointZ);
    //fFlashMatchTree->Branch("TotalConteggiCorr",                  &fTotalConteggiCorr);

/*    fLargestFlashTree = tfs->make<TTree>("LargestFlashTree","LargestFlashTree");
    fLargestFlashTree->Branch("EventID",                     &fEventID,   "EventID/I");
    fLargestFlashTree->Branch("TrueX",                       &fTrueX,     "TrueX/F");
    fLargestFlashTree->Branch("TrueY",                       &fTrueY,     "TrueY/F");
    fLargestFlashTree->Branch("TrueZ",                       &fTrueZ,     "TrueZ/F");
    fLargestFlashTree->Branch("TrueT",                       &fTrueT,     "TrueT/F");
    fLargestFlashTree->Branch("DetectedT",                   &fDetectedT, "DetectedT/F");
    fLargestFlashTree->Branch("TrueE",                       &fTrueE,     "TrueE/F");
    fLargestFlashTree->Branch("TruePDG",                     &fTruePDG,   "TruePDG/I");
    fLargestFlashTree->Branch("TrueCCNC",                    &fTrueCCNC,  "TrueCCNC/I");
    fLargestFlashTree->Branch("NFlashes",                    &fNFlashes,  "NFlashes/I");
    fLargestFlashTree->Branch("FlashID",                     &fFlashID,   "FlashID/I");
    fLargestFlashTree->Branch("YCenter",                     &fYCenter,   "YCenter/F");
    fLargestFlashTree->Branch("ZCenter",                     &fZCenter,   "ZCenter/F");
    fLargestFlashTree->Branch("YWidth",                      &fYWidth,    "YWidth/F");
    fLargestFlashTree->Branch("ZWidth",                      &fZWidth,    "ZWidth/F");
    fLargestFlashTree->Branch("Time",                        &fTime,      "Time/F");
    fLargestFlashTree->Branch("TimeWidth",                   &fTimeWidth, "TimeWidth/F");
    fLargestFlashTree->Branch("TimeDiff",                    &fTimeDiff,  "TimeDiff/F");
    fLargestFlashTree->Branch("TotalPE",                     &fTotalPE,   "TotalPE/F");
    fLargestFlashTree->Branch("NOpDets",                     &fNOpDets,   "NOpDets/I");
    fLargestFlashTree->Branch("NHitOpDets",                  &fNHitOpDets,"NHitOpDets/I");
    fLargestFlashTree->Branch("PEsPerOpDetVector",           &fPEsPerOpDetVector);
    fLargestFlashTree->Branch("Purity",                      &fPurity,    "Purity/F");
    fLargestFlashTree->Branch("Distance",                    &fDistance,  "Distance/F");
    fLargestFlashTree->Branch("RecoX",                       &fRecoX,     "RecoX/F");
    fLargestFlashTree->Branch("TotEdep",                     &fTotEdep,     "TotEdep/F");
    fLargestFlashTree->Branch("EnergyDepositionVector",      &fEnergyDepositionVector);
    fLargestFlashTree->Branch("TotalCharge",                 &fTotalCharge,  "TotalCharge/F");
    fLargestFlashTree->Branch("ADCSum",                      &fADCSum,       "ADCSum/F");
    fLargestFlashTree->Branch("PeakAmplitude",               &fPeakAmplitude,"PeakAmplitude/F");
    fLargestFlashTree->Branch("HitMultiplicity",             &fHitMultiplicity);
    //fLargestFlashTree->Branch("FirstHitTime",                &fFirstHitTime,  "FirstHitTime/F");
    //fLargestFlashTree->Branch("LastHitTime",                 &fLastHitTime,   "LastHitTime/F");
    fLargestFlashTree->Branch("MeanHitTime",                 &fMeanHitTime,   "MeanHitTime/F");
    fLargestFlashTree->Branch("HitPeakTime",                 &fHitPeakTime);
    fLargestFlashTree->Branch("HitDist",                &fHitDist,  "HitDist/F");
  
    fSelectedFlashTree = tfs->make<TTree>("SelectedFlashTree","SelectedFlashTree");
    fSelectedFlashTree->Branch("EventID",                     &fEventID,   "EventID/I");
    fSelectedFlashTree->Branch("TrueX",                       &fTrueX,     "TrueX/F");
    fSelectedFlashTree->Branch("TrueY",                       &fTrueY,     "TrueY/F");
    fSelectedFlashTree->Branch("TrueZ",                       &fTrueZ,     "TrueZ/F");
    fSelectedFlashTree->Branch("TrueT",                       &fTrueT,     "TrueT/F");
    fSelectedFlashTree->Branch("DetectedT",                   &fDetectedT, "DetectedT/F");
    fSelectedFlashTree->Branch("TrueE",                       &fTrueE,     "TrueE/F");
    fSelectedFlashTree->Branch("TruePDG",                     &fTruePDG,   "TruePDG/I");
    fSelectedFlashTree->Branch("TrueCCNC",                    &fTrueCCNC,  "TrueCCNC/I");
    fSelectedFlashTree->Branch("NFlashes",                    &fNFlashes,  "NFlashes/I");
    fSelectedFlashTree->Branch("FlashID",                     &fFlashID,   "FlashID/I");
    fSelectedFlashTree->Branch("YCenter",                     &fYCenter,   "YCenter/F");
    fSelectedFlashTree->Branch("ZCenter",                     &fZCenter,   "ZCenter/F");
    fSelectedFlashTree->Branch("YWidth",                      &fYWidth,    "YWidth/F");
    fSelectedFlashTree->Branch("ZWidth",                      &fZWidth,    "ZWidth/F");
    fSelectedFlashTree->Branch("Time",                        &fTime,      "Time/F");
    fSelectedFlashTree->Branch("TimeWidth",                   &fTimeWidth, "TimeWidth/F");
    fSelectedFlashTree->Branch("TimeDiff",                    &fTimeDiff,  "TimeDiff/F");
    fSelectedFlashTree->Branch("TotalPE",                     &fTotalPE,   "TotalPE/F");
    fSelectedFlashTree->Branch("NOpDets",                     &fNOpDets,   "NOpDets/I");
    fSelectedFlashTree->Branch("NHitOpDets",                  &fNHitOpDets,"NHitOpDets/I");
    fSelectedFlashTree->Branch("PEsPerOpDetVector",           &fPEsPerOpDetVector);
    fSelectedFlashTree->Branch("Purity",                      &fPurity,    "Purity/F");
    fSelectedFlashTree->Branch("Distance",                    &fDistance,  "Distance/F");
    fSelectedFlashTree->Branch("RecoX",                       &fRecoX,     "RecoX/F");
    fSelectedFlashTree->Branch("TruePxallpart",               &fTruePxallpart);
    fSelectedFlashTree->Branch("TruePyallpart",               &fTruePyallpart);
    fSelectedFlashTree->Branch("TruePzallpart",               &fTruePzallpart);
    fSelectedFlashTree->Branch("TrueEallpart",                &fTrueEallpart);
    fSelectedFlashTree->Branch("TrueAllPDG",                  &fTrueAllPDG);
    fSelectedFlashTree->Branch("TotEdep",                     &fTotEdep,   "TotEdep/F");
    fSelectedFlashTree->Branch("LdepSim",                     &fLdepSim,     "LdepSim/F");
    fSelectedFlashTree->Branch("TotEion",                     &fTotEion,   "TotEion/F");
    fSelectedFlashTree->Branch("TotGammaScint",               &fTotGammaScint,"TotGammaScint/F");
    fSelectedFlashTree->Branch("EnergyDepositionVector",      &fEnergyDepositionVector);
    fSelectedFlashTree->Branch("TotalCharge",                 &fTotalCharge,  "TotalCharge/F");
    fSelectedFlashTree->Branch("ADCSum",                      &fADCSum,       "ADCSum/F");
    fSelectedFlashTree->Branch("PeakAmplitude",               &fPeakAmplitude,"PeakAmplitude/F");
    fSelectedFlashTree->Branch("HitMultiplicity",             &fHitMultiplicity);
    //fSelectedFlashTree->Branch("FirstHitTime",                &fFirstHitTime,  "FirstHitTime/F");
    //fSelectedFlashTree->Branch("LastHitTime",                 &fLastHitTime,   "LastHitTime/F");
    fSelectedFlashTree->Branch("MeanHitTime",                 &fMeanHitTime,   "MeanHitTime/F");
    fSelectedFlashTree->Branch("HitPeakTime",                 &fHitPeakTime);
    fSelectedFlashTree->Branch("HitDist",                &fHitDist,  "HitDist/F");
*/ 
    if (!fOpDetWaveformLabel.empty()) {
      fCountTree = tfs->make<TTree>("CountWaveforms","CountWaveforms");
      fCountTree->Branch("EventID",       &fEventID,      "EventID/I");
      fCountTree->Branch("nwaveforms1pe", &fnwaveforms1pe, "nwaveforms1pe/I");
      fCountTree->Branch("nwaveforms2pe", &fnwaveforms2pe, "nwaveforms2pe/I");
      fCountTree->Branch("nwaveforms3pe", &fnwaveforms3pe, "nwaveforms3pe/I");
    }
    /*
    fRecoEfficiencyVsE         = tfs->make<TEfficiency>("recoEfficiencyVsE",         ";Energy (GeV);Efficiency",  fNBinsE, fLowE, fHighE);
    fRecoEfficiencyVsX         = tfs->make<TEfficiency>("recoEfficiencyVsX",         ";Position (cm);Efficiency", fNBinsX, fLowX, fHighX);
    fRecoEfficiencyVsXandE     = tfs->make<TEfficiency>("recoEfficiencyVsXandE",     ";Position (cm);Energy (GeV);Efficiency", fNBinsX, fLowX, fHighX, fNBinsE, fLowE, fHighE);
    fLargestEfficiencyVsE      = tfs->make<TEfficiency>("largestEfficiencyVsE",      ";Energy (GeV);Efficiency",  fNBinsE, fLowE, fHighE);
    fLargestEfficiencyVsX      = tfs->make<TEfficiency>("largestEfficiencyVsX",      ";Position (cm);Efficiency", fNBinsX, fLowX, fHighX);
    fLargestEfficiencyVsXandE  = tfs->make<TEfficiency>("largestEfficiencyVsXandE",  ";Position (cm);Energy (GeV);Efficiency", fNBinsX, fLowX, fHighX, fNBinsE, fLowE, fHighE);
    fSelectedEfficiencyVsE     = tfs->make<TEfficiency>("selectedEfficiencyVsE",     ";Energy (GeV);Efficiency",  fNBinsE, fLowE, fHighE);
    fSelectedEfficiencyVsX     = tfs->make<TEfficiency>("selectedEfficiencyVsX",     ";Position (cm);Efficiency", fNBinsX, fLowX, fHighX);
    fSelectedEfficiencyVsXandE = tfs->make<TEfficiency>("selectedEfficiencyVsXandE", ";Position (cm);Energy (GeV);Efficiency", fNBinsX, fLowX, fHighX, fNBinsE, fLowE, fHighE);
    */ 
  }

  //-----------------------------------------------------------------------
  // Destructor
  FlashMatchAnaEdep::~FlashMatchAnaEdep()
  {}

  //-----------------------------------------------------------------------
  void FlashMatchAnaEdep::beginJob()
  {}

  //-----------------------------------------------------------------------
  void FlashMatchAnaEdep::analyze(const art::Event& evt)
  {
    // Get the required services
    art::ServiceHandle< cheat::PhotonBackTrackerService > pbt;
    art::ServiceHandle< cheat::ParticleInventoryService > pinv;
    art::ServiceHandle< art::TFileService > tfs;

    //pbt->Rebuild(evt);


    // Record the event ID
    fEventID = evt.id().event();


    ///////////////////////////////////////////////////
    // Count waveforms if a waveform label was given //
    ///////////////////////////////////////////////////

    if (!fOpDetWaveformLabel.empty()) {
      fnwaveforms1pe = 0;
      fnwaveforms2pe = 0;
      fnwaveforms3pe = 0;
      art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
      if (evt.getByLabel(fOpDetWaveformLabel, wfHandle)) {
        fnwaveforms1pe = wfHandle->size();

        for (auto wf: *wfHandle) {
          auto it = max_element(std::begin(wf), std::end(wf));
          double peak = *it - fBaseline;
          if ( peak > (1.5*fPE)) {
            ++fnwaveforms2pe;

            if ( peak > (2.5*fPE) )
              ++fnwaveforms3pe;
          }
        }
        fCountTree->Fill();
      }
    }


    //////////////////////////////////////
    // Access all the Flash Information //
    //////////////////////////////////////

    // Get flashes from event
    art::Handle< std::vector< recob::OpFlash > > FlashHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if (evt.getByLabel(fOpFlashModuleLabel, FlashHandle)) {
      art::fill_ptr_vector(flashlist, FlashHandle);
      std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
    }
    else {
      mf::LogWarning("FlashMatchAnaEdep") << "Cannot load any flashes. Failing";
      return;
    }
    
 
     //////////////////////////////////////
    // Access all the OpHit Information //
    //////////////////////////////////////
       
    art::Handle< std::vector< recob::OpHit > > HitHandle;
    std::vector<art::Ptr<recob::OpHit> > hitlist;
    if (evt.getByLabel(fOpHitModuleLabel, HitHandle)) {
      art::fill_ptr_vector(hitlist, HitHandle);
    }

    // Get total PE in all Ophits
    fSumPE = 0;
    for (auto hit: hitlist) {
        fSumPE += hit->PE();
        fOpHitPeakTime.emplace_back(hit->PeakTime());  
        fOpHitFastToTotal.emplace_back(hit->FastToTotal());                 
        }
        
      fPEperOpDet.clear();
      //fPEperOpDetcorr.clear();
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        fPEperOpDet.emplace_back(0);
       // fPEperOpDetcorr.emplace_back(0);
      }
      
     unsigned int iC = 0;   
     for (auto hit:hitlist) {     
      iC = hit->OpChannel();
      unsigned int iOD = geom->OpDetFromOpChannel(iC);
      fPEperOpDet[iOD] += hit->PE();
     }
     
    // Get assosciations between flashes and hits
    //art::FindManyP< recob::OpHit > Assns(flashlist, evt, fOpFlashModuleLabel);

    /////////////////////////
    // G4 Deposited Energy
    //////////////////////
    fTotEdep = 0;
    
    //IonAndScint
    fTotGammaScint = 0;
    fTotEion = 0;
    fLdepSim=0;
    fLdepGeom=0;
    fGammaScintdE = 0;
    fEiondE = 0;
    int countdep=0;
    int countmedio=0;
    
    std::cout<<"Looking for deposited energy of IonAndScint Algo"<<std::endl;
    int nSimEnergyDeposits = 0;
    double dedxsteps=0;
    double dEdx=0;
    art::Handle< std::vector<sim::SimEnergyDeposit> > energyDepositHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit> > energyDepositlist;
      if(evt.getByLabel(fSimLabel, energyDepositHandle)){
	art::fill_ptr_vector(energyDepositlist, energyDepositHandle);
	nSimEnergyDeposits = energyDepositlist.size();
	//look for start point and end point of total deposition, make a scan in Z (beam dir) and get corresponding x,y
	double edepstartx = 0;
	double edepstarty = 0;
	double edepstartz = energyDepositlist[0]->StartZ();
	double edependx = 0;
	double edependy = 0;
	double edependz =  energyDepositlist[0]->EndZ();
	double thisposition = 0;
	int istart = 0;
	int iend = 0;
	
       //mappatura del muone
       	/*mappatura totale
       	
       	for(int i = 0; i < nSimEnergyDeposits; i++){
       	if(energyDepositlist[i]->PdgCode()==13){
       	fIdDep.emplace_back(i);
       	fTrackDep.emplace_back(energyDepositlist[i]->TrackID());
       	fMiddledepx.emplace_back(energyDepositlist[i]->MidPointX());
       	fMiddledepy.emplace_back(energyDepositlist[i]->MidPointY());
       	fMiddledepz.emplace_back(energyDepositlist[i]->MidPointZ());
       	//fPhotonperDep.emplace_back(energyDepositlist[i]->NumPhotons());
       	countdep++;
       	}
       	
       	}*/
       //mappatura della traccia dividendo in step da 2.4 cm, ossia circa ogni 60 depositi
       for(int i = 0; i < nSimEnergyDeposits; i++){
       	if(energyDepositlist[i]->PdgCode()==13){
       	fIdDep.emplace_back(i);
       	fTrackDep.emplace_back(energyDepositlist[i]->TrackID());
       	fMiddledepx.emplace_back(energyDepositlist[i]->MidPointX());
       	fMiddledepy.emplace_back(energyDepositlist[i]->MidPointY());
       	fMiddledepz.emplace_back(energyDepositlist[i]->MidPointZ());
       	//fPhotonperDep.emplace_back(energyDepositlist[i]->NumPhotons());
       	countdep++;
       	}
       	
       	}
       	
       	for(int i = 29; i < countdep; i+=60){
       	
       	fMiddlex.emplace_back(energyDepositlist[i]->MidPointX());
       	fMiddley.emplace_back(energyDepositlist[i]->MidPointY());
       	fMiddlez.emplace_back(energyDepositlist[i]->MidPointZ());
       	
       	
       	countmedio++;
       	
       	}
       	
	
       std::cout<<"Number of deposits: "<<countdep<<std::endl;
       
       for(int i = 0; i < nSimEnergyDeposits; i++){
          //points of the simulated track
          fPointX.emplace_back(energyDepositlist[i]->StartX());
          fPointY.emplace_back(energyDepositlist[i]->StartY());
          fPointZ.emplace_back(energyDepositlist[i]->StartZ());
          
          fGammaScint.emplace_back(energyDepositlist[i]->NumPhotons());
              
	  fEdep = energyDepositlist[i]->E(); 
	  fEnergyDepositionVector.emplace_back(fEdep);
	  fTotEdep += fEdep;
          fStepEdepCumVector.emplace_back(fTotEdep);
	  fTotGammaScint += energyDepositlist[i]->NumPhotons();
	  fTotEion += energyDepositlist[i]->NumElectrons();
	  //find starting point, scan along z
	  thisposition = energyDepositlist[i]->StartZ();
	  if(thisposition < edepstartz ){
	    edepstartz=thisposition;
	    istart=i;
	  }
	  //find end point, scan along z
	  thisposition = energyDepositlist[i]->EndZ();
	  if(thisposition > edependz ){
	    edependz=thisposition;
	    iend=i;
	  }
	  //save single steps of edep length
	  fStepLength.emplace_back( energyDepositlist[i]->StepLength()); //in cm
	 //sum of single steps of edep length
	  fLdepSim += energyDepositlist[i]->StepLength(); //in cm
	  fStepLCumVector.emplace_back(fLdepSim); //cumulative of energy in each step
	  
	  dedxsteps += fEdep/fLdepSim; //in MeV/cm
	}
	
	//std::cout<<"Number of deposits: "<<nSimEnergyDeposits<<" Single step L(cm): "<<energyDepositlist[2]->StepLength()<<" Total number of emittend photons: "<<fTotGammaScint<<" and ionization electrons: "<<fTotEion<<" ---> Check w.r.t path length: "<<std::endl;

	edepstartx=energyDepositlist[istart]->StartX();
	edepstarty=energyDepositlist[istart]->StartY();
	edependx=energyDepositlist[iend]->EndX();
	edependy=energyDepositlist[iend]->EndY();
	double diffx=TMath::Power((edependx-edepstartx),2);
	double diffy=TMath::Power((edependy-edepstarty),2);
	double diffz=TMath::Power((edependz-edepstartz),2);
	fLdepGeom=TMath::Sqrt(diffx+diffy+diffz);

        std::cout<<"Sum of every deposit step length is "<<fLdepSim<<" cm."<<std::endl;
	std::cout<<"Scan in Z, Start of energy deposition (x,y,z): "<<edepstartx<<" "<<edepstarty<<" "<<edepstartz<<" End (x,y,z): "<<edependx<<" "<<edependy<<" "<<edependz<<" cm --> length: "<<fLdepGeom<<" cm"<<std::endl;
	std::cout<<"Gamma/path length: "<<fTotGammaScint/fLdepSim<<" Ion. electrons/path length: "<<fTotEion/fLdepSim<<std::endl;
	fGammaScintdE=fTotGammaScint/fTotEdep;
	fEiondE=fTotEion/fTotEdep;
	std::cout<<"Gamma/TotEdep: "<<fGammaScintdE<<" Ion. electrons/TotEdep: "<<fEiondE<<std::endl;
	std::cout<<"(Gamma/TotEdep)/dx: "<<fGammaScintdE/fLdepSim<<" (Ion.electrons/TotEdep)/dx: "<<fEiondE/fLdepSim<<std::endl;
	dEdx=fTotEdep/fLdepSim;
	std::cout<<"Sum of every step de/dx: "<<dedxsteps<<" MeV/cm, ...divinding TotEdep by Total L--> "<<std::endl;  
	std::cout<<"Total dE/dX: "<<dEdx<<" (dx is sum of steps L). MeV/cm"<<std::endl;
	
      }else{
	mf::LogWarning("FlashMatchAnaEdep") << "Cannot Find Deposited Energy. Failing";
	return;
      }

    std::cout<<"TotEdep is: "<<fTotEdep<<" (MeV) ";
    
    //mappatura
    
    int inc=0;
    int inx=0;
    int num_1=0;
    int num_2=0;
    float opdetx_1=0.05;
    float opdety_1=591.569;
    float opdetz_1=1357.47;
    float opdetx_2=-0.05;
    float opdety_2=-31.135;
    float opdetz_2=1357.47;
    float totaleperdep_1=0;
    float totaleperdep_2=0;
    float totaleperdep_1PE=0;
    float totaleperdep_2PE=0;
    
      float z[24]={1357.47,1308.67,1245.87,1197.07,1125.08,1076.28,1013.48,964.679,892.689,843.889,781.089,732.289,660.299,611.499,548.699,499.899,427.909,379.109,316.309,267.509,195.519,146.719,83.9188,35.1188};
    
    for (int iOZ = 0; iOZ < countdep; ++iOZ){
    
   inc=0;
   inx=0;
   num_1=0;
   num_2=0;
   opdetx_1=0.05;
   opdety_1=591.569;
   opdetz_1=1357.47;
   opdetx_2=-0.05;
   opdety_2=-31.135;
   opdetz_2=1357.47;
   totaleperdep_1=0;
   totaleperdep_2=0;
   totaleperdep_1PE=0;
   totaleperdep_2PE=0;
    
    for(int iOD = 0; iOD < 240; ++iOD){
    
    if (inc==10){
        num_1++;
        opdetz_1=z[num_1];
        opdety_1=591.569;
        inc=0;
        }
    
    dist1=TMath::Sqrt(TMath::Power((fMiddledepy[iOZ]-opdety_1),2)+TMath::Power((fMiddledepz[iOZ]-opdetz_1),2));
    dist2=TMath::Sqrt(TMath::Power((fMiddledepx[iOZ]-opdetx_1),2)+TMath::Power((fMiddledepy[iOZ]-opdety_1),2)+TMath::Power((fMiddledepz[iOZ]-opdetz_1),2));
    sentheta=TMath::Sqrt(1-TMath::Power((dist1/dist2),2));
    fIdOpdet.emplace_back(iOD);
    areaeff=(435.24*sentheta)/(4*3.14*TMath::Power(dist2,2));
    fAreaeff.emplace_back(areaeff);
    fConteggi.emplace_back(areaeff*fGammaScint[iOZ]);
    fTuttiPE.emplace_back(0.6*areaeff*fGammaScint[iOZ]);
    totaleperdep_1+=areaeff*fGammaScint[iOZ];
    totaleperdep_1PE+=0.6*areaeff*fGammaScint[iOZ];
    
    opdety_1-=62.271;    
    
    inc++;
    }
    
    for(int iOD = 240; iOD < 480; ++iOD){
    
    if (inx==10){
        num_2++;
        opdetz_2=z[num_2];
        inx=0;
        opdety_2=-31.135;
        }
    
    dist1=TMath::Sqrt(TMath::Power((fMiddledepy[iOZ]-opdety_2),2)+TMath::Power((fMiddledepz[iOZ]-opdetz_2),2));
    dist2=TMath::Sqrt(TMath::Power((fMiddledepx[iOZ]-opdetx_2),2)+TMath::Power((fMiddledepy[iOZ]-opdety_2),2)+TMath::Power((fMiddledepz[iOZ]-opdetz_2),2));
    sentheta=TMath::Sqrt(1-TMath::Power((dist1/dist2),2));
    fIdOpdet.emplace_back(iOD);
    areaeff=(435.24*sentheta)/(4*3.14*TMath::Power(dist2,2));
    fAreaeff.emplace_back(areaeff);
    fConteggi.emplace_back(areaeff*fGammaScint[iOZ]);
    fTuttiPE.emplace_back(0.6*areaeff*fGammaScint[iOZ]);
    totaleperdep_2+=areaeff*fGammaScint[iOZ];
    totaleperdep_2PE+=0.6*areaeff*fGammaScint[iOZ];
    
    opdety_2-=62.271;
    
    inx++;
	  }
	  
    fConteggiTotDep.emplace_back(totaleperdep_1+totaleperdep_2);  
    fConteggiTotDepPE.emplace_back(totaleperdep_1PE+totaleperdep_2PE);	  
    }
    
    int totale=0;
    totale=countdep*480;
    float ConteggiTot=0;
    float ConteggiTotPE=0;
    
    for (int iOZ = 0; iOZ < totale; ++iOZ){
    
    ConteggiTot+=fConteggi[iOZ];
    ConteggiTotPE+=fTuttiPE[iOZ];
    
    }
    
    fConteggiTot.emplace_back(ConteggiTot);
    fConteggiTotPE.emplace_back(ConteggiTotPE);
    
    
    float result=0;
    float resultPE=0;
    
    for(int iOD = 0; iOD < 480; ++iOD){
    
    result=0;
    resultPE=0;
    for (int iOZ = 0; iOZ < countdep; ++iOZ){
    
    result+=fConteggi[iOD+480*iOZ];
    resultPE+=fTuttiPE[iOD+480*iOZ];
    
    
    }
    fCountperOpdet.emplace_back(result);
    fCountperOpdetPE.emplace_back(resultPE);
    
    }
    
    float sommagamma=0;
    float sommape=0;
    
    for (int i=0; i<(countmedio-1); ++i){
    sommagamma=0;
    sommape=0;
      for (int z=0; z<60; ++z){
      sommagamma+=fConteggiTotDep[z+60*i];
      sommape+=fConteggiTotDepPE[z+60*i];
      }
      fGammaStep.emplace_back(sommagamma);
      fPEStep.emplace_back(sommape);
    
    }
    
    fPandoraContained50=0;
    
    double lengthvect[10];
    const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fPFParticleLabel); 
    int k=0;
    for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){

  
    if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp,evt,fPFParticleLabel,fTrackLabel)){
      
      art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp,evt,fPFParticleLabel, fTrackLabel); 
      
      if(pfp->PdgCode()==13){
      /*std::cout<<"start x"<< track->Start().X() <<std::endl;
       std::cout<<"start y"<< track->Start().Y() <<std::endl;
        std::cout<<"start z"<< track->Start().Z() <<std::endl;
        std::cout<<"end x"<< track->End().X() <<std::endl;
       std::cout<<"end y"<< track->End().Y() <<std::endl;
        std::cout<<"end z"<< track->End().Z() <<std::endl;*/
        
        
        
      if ((fabs(track->Start().X()))<310&&(fabs(track->Start().Y()))<550&&(fabs(track->Start().Z()))>50&&(fabs(track->Start().Z()))<1340&&(fabs(track->End().X()))<310
      &&(fabs(track->End().Y()))<550&&(fabs(track->End().Z()))>50&&(fabs(track->End().Z())<1340)){
      fPandoraContained50=1;
      lengthvect[k]=track->Length();
      k++;
      }
      else{
      fPandoraContained50=0;
      break;
      }
      }
      
      }
    
    }
    
    //double max=0;
    int zinc=0;
    int punti=0;
    int c=0;
    int ind=0;
    float dist_parz=0;
    float distx=0;
    float disty=0;
    float distz=0;
    int mentat=0;
    
    if (fPandoraContained50==1){
    
    
    zinc=std::distance(lengthvect, std::max_element(lengthvect, lengthvect + sizeof(lengthvect)/sizeof(lengthvect[0])));
    
    
    const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect_1 = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fPFParticleLabel);
    for(const art::Ptr<recob::PFParticle> &pfp_1: pfparticleVect_1){

  
    if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp_1,evt,fPFParticleLabel,fTrackLabel)){
      
      art::Ptr<recob::Track> track_1 = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp_1,evt,fPFParticleLabel, fTrackLabel); 
      
      if(pfp_1->PdgCode()==13){
      
      
      //traccia + lunga che è quella che mi interessa
      if(c==zinc){
      punti=track_1->CountValidPoints();
      for(int i=0; i<punti;++i){
      
      ind++;
      fTrackPointX.emplace_back(track_1->LocationAtPoint(i).X());
      fTrackPointY.emplace_back(track_1->LocationAtPoint(i).Y());
      fTrackPointZ.emplace_back(track_1->LocationAtPoint(i).Z());
      } 
      c++;
     
      
      
      }
      
      }
      
    }
    
 }  }   
 
    if (ind==1){
    for(int i=0; i<punti;++i){
    
    if (i=punti-1){
    break;
    }
    distx= pow((fTrackPointX[i]-fTrackPointX[i+1]),2);
    disty= pow((fTrackPointX[i]-fTrackPointX[i+1]),2);
    distz= pow((fTrackPointX[i]-fTrackPointX[i+1]),2);
    dist_parz= sqrt (distx + distx + distz);
    dist_totale+=dist_parz;
    if (dist_totale>1.2){
      fTrackMiddlePointX_1.emplace_back(track_1->LocationAtPoint(i).X());
      fTrackMiddlePointY_1.emplace_back(track_1->LocationAtPoint(i).Y());
      fTrackMiddlePointZ_1.emplace_back(track_1->LocationAtPoint(i).Z());
      mentat++;
      dist_totale=0;
    }
    
    }
     for for(int i=0; i<mentat;i+=2){
     fTrackMiddlePointX.emplace_back(fTrackMiddlePointX_1[i]);
     fTrackMiddlePointY.emplace_back(fTrackMiddlePointY_1[i]);
     fTrackMiddlePointZ.emplace_back(fTrackMiddlePointZ_1[i]);
    }
    
    /*int totale=0;
    for(int iOD = 0; iOD < 480; ++iOD){
    
    totale+=fCountperOpdet[iOD];}
    fTotalConteggiCorr.emplace_back(totale);*/
    
    //fPEperOpDet[iOD]
    /*int inc=0;
    int inx=0;
    float opdetx_1=0.05;
    float opdety_1=591.569;
    float opdetz_1=1357.47;
    float opdetx_2=-0.05;
    float opdety_2=-31.135;
    float opdetz_2=1357.47;
   for (int iOZ = 0; iOZ < countdep; ++iOZ){
    
   inc=0;
   inx=0; 
   opdetx_1=0.05;
   opdety_1=591.569;
   opdetz_1=1357.47;
   opdetx_2=-0.05;
   opdety_2=-31.135;
   opdetz_2=1357.47;
    
    
    for(int iOD = 0; iOD < 480; ++iOD){
    
    if (iOD<240){
    
    if (inc==10){
        opdetz_1=opdetz_1-48.8;
        inc=0;
        opdety_1=591.569;
        }
    
    inc++;
    dist1=TMath::Sqrt(TMath::Power((fMiddledepy[iOZ]-opdety_1),2)+TMath::Power((fMiddledepz[iOZ]-opdetz_1),2));
    dist2=TMath::Sqrt(TMath::Power((fMiddledepx[iOZ]-opdetx_1),2)+TMath::Power((fMiddledepy[iOZ]-opdety_1),2)+TMath::Power((fMiddledepz[iOZ]-opdetz_1),2));
    sentheta=TMath::Sqrt(1-TMath::Power((dist1/dist2),2));
    fIdOpdet.emplace_back(iOD);
    areaeff=(435.24*sentheta)/(4*3.14*TMath::Power(dist2,2));
    fAreaeff.emplace_back(areaeff);
    fConteggi.emplace_back(areaeff*fGammaScint[iOZ]);
    
    
    opdety_1=opdety_1-(inc*62.27);
    
    }
    
    if (iOD>240){
    
    if (inx==10){
        opdetz_2=opdetz_2-48.8;
        inx=0;
        opdety_2=-31.135;
        }
    
    inx++;
    dist1=TMath::Sqrt(TMath::Power((fMiddledepy[iOZ]-opdety_2),2)+TMath::Power((fMiddledepz[iOZ]-opdetz_2),2));
    dist2=TMath::Sqrt(TMath::Power((fMiddledepx[iOZ]-opdetx_2),2)+TMath::Power((fMiddledepy[iOZ]-opdety_2),2)+TMath::Power((fMiddledepz[iOZ]-opdetz_2),2));
    sentheta=TMath::Sqrt(1-TMath::Power((dist1/dist2),2));
    fIdOpdet.emplace_back(iOD);
    areaeff=(435.24*sentheta)/(4*3.14*TMath::Power(dist2,2));
    fAreaeff.emplace_back(areaeff);
    fConteggi.emplace_back(areaeff*fGammaScint[iOZ]);
    
    opdety_2=opdety_2-(inx*62.27);
    
    }}}*/
    
    //int result=0;
    
    /*for(int iOD = 0; iOD < 480; ++iOD){
    
    result=0;
    for (int iOZ = 0; iOZ < countdep; ++iOZ){
    
    result+=fConteggi[iOD+480*iOZ];
    //fCountperOpdet.emplace_back(result);
    
    }
    fCountperOpdet.emplace_back(result);
    
    }
    
    int totale=0;
    for(int iOD = 0; iOD < 480; ++iOD){
    
    totale+=fCountperOpdet[iOD];}
    fTotalConteggiCorr.emplace_back(totale);*/
    
    /*for(unsigned int iOD = 0; iOD < 240; ++iOD){
        
        
        if (inc==10){
        opdetz_1=opdetz_1-48.8;
        inc=0;
        }
        opdety_1=opdety_1-(inc*62.27);
        
        
        dist1=TMath::Sqrt(TMath::Power((fMiddledepy[i]-opdety_1),2)+TMath::Power((fMiddledepz[i]-opdetz_1),2));
        dist2=TMath::Sqrt(TMath::Power((fMiddledepx[i]-opdetx_1),2)+TMath::Power((fMiddledepy[i]-opdety_1),2)+TMath::Power((fMiddledepz[i]-opdetz_1),2));
        sentheta=TMath::Sqrt(1-TMath::Power((dist1/dist2),2));
        fAreaeff.emplace_back((957.529*sentheta)/(4*3.14*TMath::Power(dist2,2)));
        inc++;
        //fPEperOpDetcorr[iOD]=fPEperOpDet[iOD]*areaeff;
      }
    
    inc=0;
     for(unsigned int iOD = 240; iOD < 480; ++iOD){
        
        if (inc==10){
        opdetz_1=opdetz_2-48.8;
        inc=0;
        }
        opdety_1=opdety_2-(inc*62.27);
        
        dist1=TMath::Sqrt(TMath::Power((fMiddledepy[i]-opdety_2),2)+TMath::Power((fMiddledepz[i]-opdetz_2),2));
        dist2=TMath::Sqrt(TMath::Power((fMiddledepx[i]-opdetx_2),2)+TMath::Power((fMiddledepy[i]-opdety_2),2)+TMath::Power((fMiddledepz[i]-opdetz_2),2));
        sentheta=TMath::Sqrt(1-TMath::Power((dist1/dist2),2));
        fAreaeff.emplace_back((957.529*sentheta)/(4*3.14*TMath::Power(dist2,2)));
        inc++;
        //fPEperOpDetcorr[iOD]=fPEperOpDet[iOD]*areaeff;
      }
    } */
    //////////////////////////////////////
    // Access all the truth information //
    //////////////////////////////////////

    std::set<int> signal_trackids;
    geo::PlaneID planeid;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    try {
      auto MClistHandle = evt.getValidHandle<std::vector<simb::MCTruth> >(fSignalLabel);

      art::Ptr<simb::MCTruth> mctruth(MClistHandle, 0);
      if (mctruth->NParticles() == 0) {
        mf::LogError("FlashMatchAnaEdep") << "No MCTruth Particles";
      }

      // Get all the track ids associated with the signal event.
      art::FindManyP<simb::MCParticle> SignalGeantAssns(MClistHandle,evt,fGeantLabel);
      for ( size_t i = 0; i < SignalGeantAssns.size(); i++) {
        auto parts = SignalGeantAssns.at(i);
        for (auto part = parts.begin(); part != parts.end(); part++) {
          signal_trackids.emplace((*part)->TrackId());
        }
      }

      // Get just the neutrino, entry 0 from the list, and record its properties
      const simb::MCParticle& part(mctruth->GetParticle(0));
      fTrueX     = part.Vx();
      fTrueY     = part.Vy();
      fTrueZ     = part.Vz();
      fTrueT     = part.T()*1000; // ns -> us
      fTrueE     = part.E();
      fTruePDG   = part.PdgCode();
      //CCNC only for the neutrino otherwise set to 1000.
      // 0=CC 1=NC
      int intType=-1000.;
      if (mctruth->Origin() == simb::kBeamNeutrino){
	fTrueCCNC  = mctruth->GetNeutrino().CCNC();
	intType=mctruth->GetNeutrino().InteractionType();
	std::cout<<"Interaction Type is: "<<intType;
	if(intType==1091) std::cout<<" --> CCDIS"<<std::endl;
	if(intType==1001) std::cout<<" --> CCQE"<<std::endl;
	if(intType==1003) std::cout<<" --> ResCCNuProtonPiPlus"<<std::endl;
	if(intType==1004) std::cout<<" --> ResCCNuNeutronPi0"<<std::endl;
	if(intType==1005) std::cout<<" --> ResCCNuNeutronPiPlus"<<std::endl;
	std::cout<<std::endl;
      }else{
	fTrueCCNC=1000;
      }

      if (mctruth->Origin() == simb::kBeamNeutrino && fTruePDG==12 &&	fTrueCCNC==0) std::cout<<" *********** nue CC ********* "<<std::endl;
      std::cout<<"TrueE is: "<<fTrueE<<" (GeV)"<<std::endl;
      
      // Get all the paricle including neutrino, and record its properties
      unsigned int const nParticles = mctruth->NParticles();
      std::cout<<"There are: "<<nParticles<<" secondary particles"<<std::endl;

      for (unsigned int i = 0; i < nParticles; ++i) {
	simb::MCParticle const& particle = mctruth->GetParticle(i);
        fTruePxallpart    .emplace_back(particle.Px());
        fTruePyallpart    .emplace_back(particle.Py());
        fTruePzallpart    .emplace_back(particle.Pz());
        fTrueEallpart     .emplace_back(particle.E());
        fTrueAllPDG       .emplace_back(particle.PdgCode());
      }

      
      // Get the PlaneID which describes the location of the true vertex
      int plane = 0;
      double loc[] = {part.Vx(), part.Vy(), part.Vz()};
      geo::TPCID tpc = geom->FindTPCAtPosition(loc);
      if (! geom->HasTPC(tpc) ) {
        mf::LogInfo("FlashMatchAnaEdep") << "No valid TPC for " << tpc;
        return;
      }
      geo::PlaneID tempid(tpc, plane);
      planeid = tempid;
      std::cout<<" Plane ID from true Vertex: "<<planeid<<" TPC: "<<tpc<<'\n';
      
      // Convert true X to would-be charge arrival time, and convert from ticks to us, add to MC time
      double deltaTicks = detProp.ConvertXToTicks(part.Vx(), planeid);
      double deltaT = clockData.TPCTick2Time(deltaTicks);
      fDetectedT = fTrueT + deltaT;
      std::cout<<"True Vertex X to get the deltaTicks: "<<part.Vx()<<'\n';
      std::cout<<"MC True T: "<< fTrueT<<" ticks and coversion to us: "<<deltaTicks<<" "<<deltaT<< '\n';
      std::cout<<"Detected time (would-be charge arrival time) (us): "<< fDetectedT << '\n';
    }
    catch (art::Exception const& err) 
    {
      // If the error isn't that a product wasn't found, throw it back up.
      if ( err.categoryCode() != art::errors::ProductNotFound ) throw;

      // Otherwise, this event just doesn't have signal. Fill default values
      fTrueX = 0;
      fTrueY = 0;
      fTrueZ = 0;
      fTrueT = 0;
      fTrueE = 0;
      fTruePDG = 0;
      fTrueCCNC = 0;
      fDetectedT = 0;

      mf::LogError("FlashMatchAnaEdep") << "Event doesn't have signal ";
    }

    // Get the maximum possible time difference by getting number of ticks corresponding to
    // one full drift distance, and converting to time.
    double maxT=0;
    if(fBeam) maxT = clockData.TPCTick2Time(detProp.NumberTimeSamples());


    //////////////////////////////////////
    // CHARGE collected                 //
    ////////////////////////////////////// 
//    unsigned short wirePlane = geo::kZ;

    // Total charge collected in the event and time info of the hits
    fTotalCharge = 0.0;
    fTotalChargeCorr = 0.0;
    fTotalChargeFinal = 0.0;
    fADCSum = 0.0;
    fPeakAmplitude = 0.0;
    //fFirstHitTime = -1000;
    fHitDist = -1000;
    fMeanHitTime = -1000;
    //double firstHit=100000;
    double meanHit=0;
    //double lastHit=-100000;
    int collpl=0;
    auto hitListHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitsLabel);

    for(size_t iHit = 0; iHit < hitListHandle->size(); ++iHit){
      //thishit=static_cast<int>(iHit);  
      art::Ptr<recob::Hit> hitPtr(hitListHandle, iHit);
	if(geom->SignalType(hitPtr->Channel()) == geo::kCollection){
     // if(hitPtr->View() == wirePlane){
	collpl++;
	fTotalCharge += hitPtr->Integral();        
	fADCSum += hitPtr->SummedADC();
	fPeakAmplitude += hitPtr->PeakAmplitude();
	fHitMultiplicity.emplace_back(hitPtr->Multiplicity());
	//fHitPeakTimeTicks.emplace_back(hitPtr->PeakTime());
	//convert hit time from ticks to us
	double hitT = clockData.TPCTick2Time(hitPtr->PeakTime());
	fHitPeakTime.emplace_back(hitT);
	meanHit += hitT;
	//	if(thishit==0) firstHit=hitT;
	/*	if(static_cast<int>(iHit) == 0){
	  firstHit=hitT;
	  lastHit=hitT;
	  }*/
	//if(hitT<firstHit) firstHit=hitT;
	//if(hitT>lastHit) lastHit=hitT;
      }
    }
  
    if(collpl==0){
      std::cout<<"********** no hits in the Collection Plane... ***********"<<'\n';
      //  fFirstHitTime = -1000;
      //fLastHitTime = -1000;
      fMeanHitTime = -1000;
    }else{
      // fFirstHitTime = firstHit;
      //fLastHitTime = lastHit;
      fMeanHitTime = meanHit/collpl;
      fHitDist = (fMeanHitTime * 1.6)/10.; //cm
      //std::cout << "First Hit in time is " << firstHit <<" (us), Last Hit in time is " << lastHit <<" (us)" << '\n';
      std::cout << "Nb of hits: " << collpl <<" Mean hit time is " << fMeanHitTime <<" (us)" << " --> should correspond to charge travelling: "<<fHitDist<<" cm"<< '\n';    
    }
        
    // Output the total charge
    std::cout << "Total Charge from Reco Hits: " << fTotalCharge << '\n';    
    
    
   //compute the correction to the collected charge using DUNEAnaHitUtils
   std::vector<art::Ptr<recob::Hit>> HitsInColl = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsLabel),2); //Collection plane has Plane_ID = 2
   fTotalChargeCorr = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, HitsInColl);
    std::cout << "Total Charge from reco hits corrected with electron lifetime: " << fTotalChargeCorr << '\n';    
         
   fTotalChargeFinal = (fTotalChargeCorr*Wion)/Recombination;
  std::cout<<"Total Charge from reco hits corrected with electron lifetime and multiplied by Wion/R: "<<fTotalChargeFinal<<std::endl;    
    
      
    /////////////////////////
    // Analyze the flashes //
    /////////////////////////

    // Set up some flags to fill as we loop
    // through flashes. These will control
    // filling of efficiency plots after the loop.
 
    //   bool AnyReconstructed = false;
    //bool LargestFound     = false;
    // bool LargestRight     = false;
   // bool SelectedFound    = false;
    //  bool SelectedRight    = false;

    // For every OpFlash in the vector
    fNOpDets   = geom->NOpDets();
    fNFlashes  = flashlist.size();
    for(unsigned int i = 0; i < flashlist.size(); ++i)
    {
      // Get OpFlash and associated hits
      recob::OpFlash TheFlash = *flashlist[i];
      art::Ptr<recob::OpFlash> FlashP = flashlist[i];
      std::vector< art::Ptr<recob::OpHit> > hitFromFlash = pbt->OpFlashToOpHits_Ps(FlashP);
      std::vector< art::Ptr<recob::OpHit> > matchedHits = pbt->OpFlashToOpHits_Ps(FlashP);
      //std::vector< art::Ptr<recob::OpHit> > matchedHits = Assns.at(i);

      // Calculate the flash purity
      double purity = pbt->OpHitCollectionPurity(signal_trackids, matchedHits);

      // Calcuate relative detection time
      //if event is not from neutrino beam add true MC time
      double flashrealtime= TheFlash.Time() + fTrueT;
      //       double timeDiff = fDetectedT - TheFlash.Time();
      double timeDiff = fDetectedT - flashrealtime;

      if (!planeid) {
        // planeid isn't valid
        fRecoX = 0;
      }
      else {
        double ticks = clockData.Time2Tick(timeDiff);
        fRecoX = detProp.ConvertTicksToX(ticks, planeid);
      }

      ///********************* USE THIS For Neutrino Events, skip for Single Part muons
      // Check if this is a possible flash (w/in 1 drift window)
      // Otherwise, skip it
      if(fBeam){
        if (timeDiff < -10 || timeDiff > maxT){
	  mf::LogError("FlashMatchAnaEdep") << "Skipping Flash: not w/in 1 drift window (timeDiff < -10 OR timeDiff > maxT)" <<
                       " \n DetectedTime " << fDetectedT << " FlashTime "<< TheFlash.Time() << " --> TimeDiff: " << timeDiff <<
	  "\n MaxT is the maximum possible time difference by getting number of ticks corresponding to one full drift distance, and converting to time: "<< maxT;
	  continue;
	}
      }
      //*************************************************************
      
      // Put flash info into variables
      fFlashID     = i;
      fYCenter     = TheFlash.YCenter();
      fZCenter     = TheFlash.ZCenter();
      fYWidth      = TheFlash.YWidth();
      fZWidth      = TheFlash.ZWidth();
      fTime        = TheFlash.Time();
      fTimeWidth   = TheFlash.TimeWidth();
      fTimeDiff    = timeDiff;
      fTotalPE     = TheFlash.TotalPE();
      fPurity      = purity;
      std::cout<<"Flash: "<<fFlashID<<" PE of this Flash: "<<fTotalPE<<'\n';
      std::cout<<"Flash time: "<<fTime<<'\n';
      std::cout<<"True X: "<<fTrueX<<" RecoX from ticks "<<fRecoX<< '\n';

      // Calculate distance from MC truth vertex in the Y-Z plane
      fDistance = sqrt( pow(fTrueY-fYCenter,2) +  pow(fTrueZ-fZCenter,2) );
      std::cout<<"Flash distance from MC truth vertex in the Y-Z plane: "<<fDistance<< '\n';
      std::cout<<"TrueX+flash distance gives: "<< sqrt( pow(fTrueX,2) +  pow(fDistance,2) )<< '\n';

      // Loop through all the opdets with hits in this flash
      fPEsPerOpDetVector.clear();
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        fPEsPerOpDetVector.emplace_back(0);
      }
      for(unsigned int iC=0; iC < geom->NOpChannels(); ++iC)
      {
        unsigned int iOD = geom->OpDetFromOpChannel(iC);
        fPEsPerOpDetVector[iOD] += TheFlash.PE(iC);
      }

      fNHitOpDets = 0;
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        if (fPEsPerOpDetVector[iOD] > 0) ++fNHitOpDets;
      }
      fNHitOpDetVector.emplace_back(fNHitOpDets);

      // Add flash info to the tree of all possible flashes
      fFlashIDVector    .emplace_back(fFlashID);
      fYCenterVector    .emplace_back(fYCenter);
      fZCenterVector    .emplace_back(fZCenter);
      fYWidthVector     .emplace_back(fYWidth);
      fZWidthVector     .emplace_back(fZWidth);
      fTimeVector       .emplace_back(fTime);
      fTimeWidthVector  .emplace_back(fTimeWidth);
      fTimeDiffVector   .emplace_back(fTimeDiff);
      fTotalPEVector    .emplace_back(fTotalPE);
      fPurityVector     .emplace_back(fPurity);
      fDistanceVector   .emplace_back(fDistance);
      fRecoXVector      .emplace_back(fRecoX);


      // Did we reconstruct any flashes with signal in them?
      // if (fPurity > 0) AnyReconstructed = true;
/*
      // First == Largest, so if this is the first flash it is also the largest.
      // So, fill the LargestFlash tree and the LargestFlash efficiency plots
      if (!LargestFound) {

        // Write out the info into the tree for the largest flash
        fLargestFlashTree->Fill();

        // Record that we found the largest flash
        // and if we got it right
        LargestFound = true;
	//        if (fPurity > 0) LargestRight = true;
      }


      // The first time we get into here we have the largest flash that is
      // within the distance cut. So, fill the SelectedFlash tree and the
      // selected flash efficiency plots
      if (!SelectedFound && fDistance < fDistanceCut) {

        // Write out the info for the selected flash
        fSelectedFlashTree->Fill();

        // Record that we found the selected flash
        // and if we got it right
        SelectedFound = true;
	//        if (fPurity > 0) SelectedRight = true;
      }
*/
    }

    // Fill these TEfficiencies once for every event
    // but use the booleans to decide if it was
    // "selected" or not.
    /*
    fRecoEfficiencyVsE->Fill(AnyReconstructed, fTrueE);
    fRecoEfficiencyVsX->Fill(AnyReconstructed, fTrueX);
    fRecoEfficiencyVsXandE->Fill(AnyReconstructed, fTrueX, fTrueE);

    fLargestEfficiencyVsE->Fill(LargestRight, fTrueE);
    fLargestEfficiencyVsX->Fill(LargestRight, fTrueX);
    fLargestEfficiencyVsXandE->Fill(LargestRight, fTrueX, fTrueE);

    fSelectedEfficiencyVsE->Fill(SelectedRight, fTrueE);
    fSelectedEfficiencyVsX->Fill(SelectedRight, fTrueX);
    fSelectedEfficiencyVsXandE->Fill(SelectedRight, fTrueX, fTrueE);
    */

    ///////////////////////////////////////////////
    // Write out the FlashMatchTree and clean up //
    ///////////////////////////////////////////////

    fFlashMatchTree->Fill();
    fFlashIDVector              .clear();
    fYCenterVector              .clear();
    fZCenterVector              .clear();
    fYWidthVector               .clear();
    fZWidthVector               .clear();
    fTimeVector                 .clear();
    fTimeWidthVector            .clear();
    fTimeDiffVector             .clear();
    fTotalPEVector              .clear();
    fNHitOpDetVector            .clear();
    fPurityVector               .clear();
    fDistanceVector             .clear();
    fRecoXVector                .clear();
    fTruePxallpart              .clear();
    fTruePyallpart              .clear();
    fTruePzallpart              .clear();
    fTrueEallpart               .clear();
    fTrueAllPDG                 .clear();
    fEnergyDepositionVector     .clear();
    fHitMultiplicity            .clear();
    fHitPeakTime                .clear();
    fPointX                     .clear();
    fPointY                     .clear();
    fPointZ                     .clear();
    fGammaScint                 .clear();
    fStepLength		 .clear();
    fStepLCumVector		 .clear();
    fStepEdepCumVector          .clear();
    fOpHitPeakTime  		 .clear();
    fOpHitFastToTotal  	 .clear();    
    fMiddledepx.clear(); 
    fMiddledepy.clear(); 
    fMiddledepz.clear(); 
    fIdDep.clear(); 
    fTrackDep.clear(); 
    fAreaeff.clear(); 
    //fPhotonperDep.clear();
     fIdOpdet.clear();
    //fCountperOpdet.clear();
    //fTotalConteggiCorr.clear();
    fConteggi.clear();
    fConteggiTot.clear();
    fCountperOpdet.clear();
    
    fConteggiTotDep.clear();
    fConteggiTotDepPE.clear();
    fTuttiPE.clear();
    fCountperOpdetPE.clear();
    fConteggiTotPE.clear();
    
    fGammaStep.clear();
    fPEStep.clear();
    
    fMiddlex.clear(); 
    fMiddley.clear(); 
    fMiddlez.clear(); 
    fTrackPointX.clear(); 
    fTrackPointY.clear(); 
    fTrackPointZ.clear(); 
    
    fTrackPointZ.clear();
    fTrackPointY.clear();
    fTrackPointX.clear();
    
    fTrackPointZ_1.clear();
    fTrackPointY_1.clear();
    fTrackPointX_1.clear();
    
  }

  //-----------------------------------------------------------------------

 void FlashMatchAnaEdep::endJob(){

 }

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(FlashMatchAnaEdep)
}
