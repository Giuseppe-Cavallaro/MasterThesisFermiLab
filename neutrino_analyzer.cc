// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a TTree for studying light and charge
//

#ifndef LightChargeAna_H
#define LightChargeAna_H 1

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
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larsim/Utils/TruthMatchUtils.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//DUNE
#include "dune/OpticalDetector/OpFlashSort.h"
#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaHitUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"
 


namespace opdet {

  class LightChargeAna : public art::EDAnalyzer{
  public:

    // Standard constructor and destructor for an ART module.
    LightChargeAna(const fhicl::ParameterSet&);
    virtual ~LightChargeAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The analyzer routine, called once per event.
    void analyze (const art::Event&);

    void endJob();
    //    void findDepositionVertex(std::vector<sim::SimChannel> SCHandle);
    //void checkDepositionVertex(std::vector<sim::SimChannel> SCHandle, double x, double y, double z);
    void readLightMap();
    float getfvis(float x, float y, float z, float q);
    
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
    std::string fSimulationProducerLabel;  // Input tag for Sim
    std::string fHitsLabel;                //Input tag for Charge hits
    std::string fParticleModuleLabel;       //Input tag for Reco Particles
    std::string fTrackLabel;                //Input tag for Reco Tracks
    std::string fShowerLabel;               //Input tag for Reco Showers
    std::string fSpacePointLabel;      //Input Tag to go from hits to point in space
    std::string fHitToSpacePointLabel;
    bool fBeam;                            // Simulated events are beam neutrinos
    bool fSimCh;                      // Simulated events are from official MC 

    TTree * fFlashMatchTree;
  
   //Values needed to convert charge to energy
   float RecombFactor=0.64;
   float Wion=23.6e-6; //MeV/electron
   float Wph=19.5e-6; //MeV/electron

    //Light map
    std::vector<float> Bin_x;  
    std::vector<float> Bin_y; 
    std::vector<float> Bin_z;
    std::vector<float> Bin_w;
   
    // Parameters from the fhicl
    int   fNBinsE;
    float fLowE;
    float fHighE;
    int   fNBinsX;
    float fLowX;
    float fHighX;
    float fDistanceCut;

    Int_t fEventID;

    Int_t fContained50=0;
    Int_t fmuPandoraContained50=0;
    Float_t fTrueX;
    Float_t fTrueY;
    Float_t fTrueZ;
    Float_t fTrueT;
    Float_t fDetectedT;
    Float_t fTrueE;
    Float_t fmuTrueE;
    Int_t   fTruePDG;
    Int_t   fTrueCCNC;
    Float_t fRecoX;
//********Ion and scint
    Float_t fEdep;
    Float_t fTotEdep;
    Float_t fmuTotEdep;
    Int_t fTotGammaScint;
    Int_t fTotEion;
//*******SimChannel if not Ion and Scint
    //std::vector<double> fMCdEdxVec;
    //double fMaxMCdEdx = 0.;
    double dE;
 // distance of energy deposition from conversion point
    //double trackLength;
 // list of event numbers corresponding to unusually low dEdx measurement
    std::vector<unsigned int> fNullList;

 // position of conversion point
    double xi = 0.;
    double yi = 0.;
    double zi = 1000.;
    // for finding conversion point
    double ziFalse = -500.;
    unsigned int nearDeps = 0;
    bool goodDepVertex = false;
    // for conversion distance
    double fConvDist;
    TH1D* fConvDistHist;
    double Vxi;           //
    double Vyi;           // primary vertex coordinates
    double Vzi;           //

    //deposited energy by plane
    double depositU,depositV,depositZ;
    Float_t fDepositedU = 0.;
    Float_t fDepositedV = 0.;
    Float_t fDepositedZ = 0.;
    Float_t fDepVertexX = -1000;
    Float_t fDepVertexY = -1000;
    Float_t fDepVertexZ = -1000;
    
    //*********Charge from reco hits
    Float_t fTotalCharge; //Integral under the calibrated signal waveform of the hit, in tick x ADC units
    Float_t fTotalCorrCharge; //Total Charge corrected for electron lifetime
    Float_t fTotalFinalCharge; //Total Charge corrected for electron lifetime and multiplied by Wion/R
    
    Float_t fmuTotalCorrCharge; //Total Charge corrected for electron lifetime (muon)
    Float_t fmuTotalFinalCharge; //Total Charge corrected for electron lifetime and multiplied by Wion/R (muon)
    
    Float_t fADCSum; //The sum of calibrated ADC counts of the hit (0. by default)
    Float_t fPeakAmplitude; //The estimated amplitude of the hit at its peak, in ADC units
    // Float_t fFirstHitTime; //First Time recorded among hits, Time of the signal peak,in us
    //Float_t fLastHitTime; //Last Time recorded among hits, Time of the signal peak, in us
    Float_t fMeanHitTime; //Last Time recorded among hits, Time of the signal peak, in us
    Float_t fHitDist; //First hit time * 1.6 mm/us --> theor distance travelled by charge (cm)

     //**** Space points
    Float_t fvis;
    Float_t Fvis;
    
    Float_t mufvis;
    Float_t muFvis;
    
    Float_t fLight;
    Float_t fE_QL;
    
    std::vector< Float_t > fTruePxallpart;
    std::vector< Float_t > fTruePyallpart;
    std::vector< Float_t > fTruePzallpart;
    
    //std::vector< Float_t > fspointx;
    //std::vector< Float_t > fspointy;
    //std::vector< Float_t > fspointz;
    //std::vector< Float_t > fspointcharge;
    
    //std::vector< Float_t > fmuspointx;
    //std::vector< Float_t > fmuspointy;
    //std::vector< Float_t > fmuspointz;
    //std::vector< Float_t > fmuspointcharge;
    
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
    std::vector< Int_t >   fFlashFrameVector;
    std::vector< Bool_t >  fInBeamFrameVector;
    std::vector< Int_t >   fOnBeamTimeVector;
    std::vector< Float_t > fTotalPEVector;
    std::vector< Float_t > fPurityVector;
    std::vector< Float_t > fDistanceVector;
    Int_t fNOpDets;
    std::vector<Int_t> fNHitOpDetVector;
    //*******Ion and Scint  
    std::vector< Float_t > fEnergyDepositionVector;

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
    //Int_t    fFlashFrame; // unused
    //Bool_t   fInBeamFrame; // unused
    //Int_t    fOnBeamTime; // unused
    Float_t  fTotalPE;
    //Ophits in the event disregarding the Flash
    Float_t  fSumPE=0;
    Float_t  fOpHitTotalArea=0;
    std::vector< Float_t > fOpHitPEVector;
    std::vector< Float_t > fOpHitPeakTimeVector;
    std::vector< Float_t > fOpHitPeakTimeAbsVector;
    std::vector< Float_t > fOpHitFastToTotalVector;
    
    Float_t  fPurity;
    Float_t  fDistance;
    Int_t    fNHitOpDets;
    std::vector< Float_t > fPEsPerOpDetVector;

    //Reco Showers and Tracks
    Int_t fNShowers;
    std::vector< Float_t > fShwEnVector;
    std::vector< Float_t > fShwdEdxVector;
    //3rd value is the biggest --> collection plane?
    Float_t fShwEn;
    Float_t fShwdEdx;
    Float_t fShwLength;
    Float_t fShwOpenAngle;
    Float_t fShwStartX,fShwStartY,fShwStartZ;
    TVector3 shwstart;

    Float_t fPandoraVtxX,fPandoraVtxY,fPandoraVtxZ;
    
    // For counting waveforms
    std::string fOpDetWaveformLabel;
    float fBaseline;
    float fPE;
    TTree * fCountTree;
    Int_t fnwaveforms1pe;
    Int_t fnwaveforms2pe;
    Int_t fnwaveforms3pe;

  };

}

#endif // LightChargeAna_H

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  LightChargeAna::LightChargeAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fOpFlashModuleLabel = pset.get<std::string>("OpFlashModuleLabel");
    fOpHitModuleLabel   = pset.get<std::string>("OpHitModuleLabel");
    fHitsLabel          = pset.get<std::string>("HitsLabel");
    fParticleModuleLabel= pset.get<std::string>("ParticleModuleLabel");
    fShowerLabel        = pset.get<std::string>("ShowerLabel");
    fTrackLabel         = pset.get<std::string>("TrackLabel");
    fSpacePointLabel    = pset.get<std::string>("SpacePointLabel");
    fSignalLabel        = pset.get<std::string>("SignalLabel");
    fGeantLabel         = pset.get<std::string>("GeantLabel");
    fSimLabel           = pset.get<std::string>("SimLabel");
    fSimulationProducerLabel = pset.get<std::string>("SimulationLabel");
    fBeam               = pset.get<bool>("Beam");
    fSimCh              = pset.get<bool>("SimCh");
    fNBinsE             = pset.get<int>("NBinsE");
    fLowE               = pset.get<float>("LowE");
    fHighE              = pset.get<float>("HighE");
    fNBinsX             = pset.get<int>("NBinsX");
    fLowX               = pset.get<float>("LowX");
    fHighX              = pset.get<float>("HighX");
    fDistanceCut        = pset.get<float>("DistanceCut");
    fHitToSpacePointLabel=pset.get<std::string>("HitToSpacePointLabel");

    fOpDetWaveformLabel = pset.get<std::string>("OpDetWaveformLabel","");
    fBaseline           = pset.get<float>("Baseline", 1500.);
    fPE                 = pset.get<float>("PE", 18.);

    art::ServiceHandle< art::TFileService > tfs;

    fFlashMatchTree = tfs->make<TTree>("FlashMatchTree","FlashMatchTree");
    fFlashMatchTree->Branch("EventID",                     &fEventID,   "EventID/I");
    fFlashMatchTree->Branch("Contained50",                     &fContained50);
    fFlashMatchTree->Branch("muPandoraContained50",                     &fmuPandoraContained50);
    fFlashMatchTree->Branch("TrueX",                       &fTrueX,     "TrueX/F");
    fFlashMatchTree->Branch("TrueY",                       &fTrueY,     "TrueY/F");
    fFlashMatchTree->Branch("TrueZ",                       &fTrueZ,     "TrueZ/F");
    fFlashMatchTree->Branch("TrueT",                       &fTrueT,     "TrueT/F");
    fFlashMatchTree->Branch("DetectedT",                   &fDetectedT, "DetectedT/F");
    fFlashMatchTree->Branch("TrueE",                       &fTrueE,     "TrueE/F");
    fFlashMatchTree->Branch("muTrueE",                       &fmuTrueE,     "muTrueE/F");
    fFlashMatchTree->Branch("TruePDG",                     &fTruePDG,   "TruePDG/I");
    fFlashMatchTree->Branch("TrueCCNC",                    &fTrueCCNC,  "TrueCCNC/I");
    fFlashMatchTree->Branch("NFlashes",                    &fNFlashes,  "NFlashes/I");
    fFlashMatchTree->Branch("FlashIDVector",               &fFlashIDVector);
    fFlashMatchTree->Branch("YCenterVector",               &fYCenterVector);
    fFlashMatchTree->Branch("ZCenterVector",               &fZCenterVector);
    fFlashMatchTree->Branch("YWidthVector",                &fYWidthVector);
    fFlashMatchTree->Branch("ZWidthVector",                &fZWidthVector);
    fFlashMatchTree->Branch("TimeVector",                  &fTimeVector);
    fFlashMatchTree->Branch("TimeWidthVector",             &fTimeWidthVector);
    fFlashMatchTree->Branch("TimeDiffVector",              &fTimeDiffVector);
    fFlashMatchTree->Branch("TotalPEVector",               &fTotalPEVector);
    fFlashMatchTree->Branch("SumPE",                       &fSumPE,   "PE/F");
    fFlashMatchTree->Branch("OpHitTotalArea",         &fOpHitTotalArea,"OpHitTotalArea/F");
    fFlashMatchTree->Branch("OpHitPEVector",               &fOpHitPEVector);
    fFlashMatchTree->Branch("OpHitPeakTimeVector",         &fOpHitPeakTimeVector);
    fFlashMatchTree->Branch("OpHitPeakTimeAbsVector",      &fOpHitPeakTimeAbsVector);
    fFlashMatchTree->Branch("OpHitFastToTotalVector",      &fOpHitFastToTotalVector);
    fFlashMatchTree->Branch("NOpDets",                     &fNOpDets, "NOpDets/I");
    fFlashMatchTree->Branch("NHitOpDetVector",             &fNHitOpDetVector);
    fFlashMatchTree->Branch("Purity",                      &fPurityVector);
    fFlashMatchTree->Branch("Distance",                    &fDistanceVector);
    fFlashMatchTree->Branch("RecoXVector",                 &fRecoXVector);
    fFlashMatchTree->Branch("TruePxallpart",               &fTruePxallpart);
    fFlashMatchTree->Branch("TruePyallpart",               &fTruePyallpart);
    fFlashMatchTree->Branch("TruePzallpart",               &fTruePzallpart);
    
    //spacepoint
    
   // fFlashMatchTree->Branch("spointx",               &fspointx);
   // fFlashMatchTree->Branch("spointy",               &fspointy);
   // fFlashMatchTree->Branch("spointz",               &fspointz);
   // fFlashMatchTree->Branch("spointcharge",               &fspointcharge);
    
   // fFlashMatchTree->Branch("muspointx",               &fmuspointx);
   // fFlashMatchTree->Branch("muspointy",               &fmuspointy);
   // fFlashMatchTree->Branch("muspointz",               &fmuspointz);
   // fFlashMatchTree->Branch("muspointcharge",               &fmuspointcharge);
    
    fFlashMatchTree->Branch("TrueEallpart",                &fTrueEallpart);
    fFlashMatchTree->Branch("TrueAllPDG",                  &fTrueAllPDG);
    // fFlashMatchTree->Branch("Edep",                        &fEdep,       "Edep/F");
    fFlashMatchTree->Branch("TotEdep",                     &fTotEdep,     "TotEdep/F");
    fFlashMatchTree->Branch("muTotEdep",                     &fmuTotEdep,     "muTotEdep/F");
    fFlashMatchTree->Branch("TotEion",                     &fTotEion,     "TotEion/I");
    fFlashMatchTree->Branch("TotGammaScint",               &fTotGammaScint,"TotGammScint/I");
    fFlashMatchTree->Branch("EdepU",                       &fDepositedU,  "EdepU/F");
    fFlashMatchTree->Branch("EdepV",                       &fDepositedV,  "EdepV/F");
    fFlashMatchTree->Branch("EdepZ",                       &fDepositedZ,  "EdepZ/F");
    fFlashMatchTree->Branch("StartEdepX",                  &fDepVertexX,  "StartEdepX/F");
    fFlashMatchTree->Branch("StartEdepY",                  &fDepVertexY,  "StartEdepY/F");
    fFlashMatchTree->Branch("StartEdepZ",                  &fDepVertexZ,  "StartEdepZ/F");
    // fFlashMatchTree->Branch("EnergyDepositionVector",      &fEnergyDepositionVector);
    fFlashMatchTree->Branch("TotalCharge",                 &fTotalCharge,  "TotalCharge/F");
    fFlashMatchTree->Branch("TotalCorrCharge",             &fTotalCorrCharge,"TotalCorrCharge/F");
    fFlashMatchTree->Branch("TotalFinalCharge",            &fTotalFinalCharge,"E_Q/F");
    
    //mu
    fFlashMatchTree->Branch("muTotalFinalCharge",            &fmuTotalFinalCharge,"muE_Q/F");
    fFlashMatchTree->Branch("Light",            &fLight,"Light/F");
    fFlashMatchTree->Branch("E_QL",            &fE_QL,"E_QL/F");
    
    fFlashMatchTree->Branch("ADCSum",                      &fADCSum,       "ADCSum/F");
    fFlashMatchTree->Branch("PeakAmplitude",               &fPeakAmplitude,"PeakAmplitude/F");
    fFlashMatchTree->Branch("HitMultiplicity",             &fHitMultiplicity);
    //fFlashMatchTree->Branch("FirstHitTime",                &fFirstHitTime,  "FirstHitTime/F");
    //fFlashMatchTree->Branch("LastHitTime",                 &fLastHitTime,  "LastHitTime/F");
    fFlashMatchTree->Branch("MeanHitTime",                 &fMeanHitTime,  "MeanHitTime/F");
    fFlashMatchTree->Branch("HitPeakTime",                 &fHitPeakTime);
    //fFlashMatchTree->Branch("HitPeakTimeTicks",            &fHitPeakTimeTicks);
    fFlashMatchTree->Branch("HitDist",                     &fHitDist,  "HitDist/F");
    fFlashMatchTree->Branch("NShowers",                    &fNShowers, "NShowers/I");
    fFlashMatchTree->Branch("ShwLength",                   &fShwLength,"ShwLength/F");
    fFlashMatchTree->Branch("ShwOpenAngle",                &fShwOpenAngle, "ShwOpenAngle/F");
    fFlashMatchTree->Branch("ShwEnVector",                 &fShwEnVector);
    fFlashMatchTree->Branch("ShwdEdxVector",               &fShwdEdxVector);
    fFlashMatchTree->Branch("ShwEn",                       &fShwEn, "ShwEn/F");
    fFlashMatchTree->Branch("ShwdEdx",                     &fShwdEdx, "ShwdEdx/F");
    fFlashMatchTree->Branch("ShwStartX",                   &fShwStartX,"ShwStartX/F");
    fFlashMatchTree->Branch("ShwStartY",                   &fShwStartY,"ShwStartY/F");
    fFlashMatchTree->Branch("ShwStartZ",                   &fShwStartZ,"ShwStartZ/F");
    fFlashMatchTree->Branch("PandoraVtxX",                 &fPandoraVtxX,"PandoraVtxX/F");
    fFlashMatchTree->Branch("PandoraVtxY",                 &fPandoraVtxY,"PandoraVtxY/F");
    fFlashMatchTree->Branch("PandoraVtxZ",                 &fPandoraVtxZ,"PandoraVtxZ/F");
    
    fFlashMatchTree->Branch("Fvis",                        &Fvis,"Fvis/F");
    fFlashMatchTree->Branch("muFvis",                        &muFvis,"muFvis/F");
   
    if (!fOpDetWaveformLabel.empty()) {
      fCountTree = tfs->make<TTree>("CountWaveforms","CountWaveforms");
      fCountTree->Branch("EventID",       &fEventID,      "EventID/I");
      fCountTree->Branch("nwaveforms1pe", &fnwaveforms1pe, "nwaveforms1pe/I");
      fCountTree->Branch("nwaveforms2pe", &fnwaveforms2pe, "nwaveforms2pe/I");
      fCountTree->Branch("nwaveforms3pe", &fnwaveforms3pe, "nwaveforms3pe/I");
    }
 
  }

  //-----------------------------------------------------------------------
  // Destructor
  LightChargeAna::~LightChargeAna()
  {}

 //-----------------------------------------------------------------------
  void LightChargeAna::beginJob()
  {

   std::cout << "******* Reading the Light Map **************";
   readLightMap();
   std::cout << "...done ***************"<<std::endl;

  }

  //-----------------------------------------------------------------------
  void LightChargeAna::analyze(const art::Event& evt)
  {
    // Get the required services
    art::ServiceHandle< cheat::PhotonBackTrackerService > pbt;
    art::ServiceHandle< cheat::BackTrackerService > bt;
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
      mf::LogWarning("LightChargeAna") << "Cannot load any flashes. Failing";
      return;
    }

    // Get assosciations between flashes and hits
    //art::FindManyP< recob::OpHit > Assns(flashlist, evt, fOpFlashModuleLabel);
    art::Handle< std::vector< recob::OpHit > > HitHandle;
    std::vector<art::Ptr<recob::OpHit> > hitlist;
    if (evt.getByLabel(fOpHitModuleLabel, HitHandle)) {
      art::fill_ptr_vector(hitlist, HitHandle);
    }

    // Get total PE in all hits
    fSumPE = 0;
    fOpHitTotalArea=0;
    for (auto hit: hitlist){
      fSumPE += hit->PE();
      fOpHitTotalArea += hit->Area();
      fOpHitPEVector.emplace_back(hit->PE());
      fOpHitPeakTimeVector.emplace_back(hit->PeakTime());
      fOpHitPeakTimeAbsVector.emplace_back(hit->PeakTimeAbs());
      fOpHitFastToTotalVector.emplace_back(hit->FastToTotal());
    }
    /////////////////////////
    // G4 Deposited Energy
    //////////////////////
    fTotEdep = 0;
    fmuTotEdep = 0;
    fTotGammaScint = 0;
    fTotEion = 0;
    double t0=0; //start time of ionization
    double minT=10000000;
    //IonAndScint
    if(!fSimCh){
      std::cout<<"Looking for deposited energy of IonAndScint Algo"<<std::endl;
      int nSimEnergyDeposits = 0;
      art::Handle< std::vector<sim::SimEnergyDeposit> > energyDepositHandle;
      std::vector<art::Ptr<sim::SimEnergyDeposit> > energyDepositlist;
      if(evt.getByLabel(fSimLabel, energyDepositHandle)){
	art::fill_ptr_vector(energyDepositlist, energyDepositHandle);
	nSimEnergyDeposits = energyDepositlist.size();
	for(int i = 0; i < nSimEnergyDeposits; i++){
	  t0=energyDepositlist[i]->StartT()*1000.; //ns-->us
	  if(t0<minT) minT=t0;
	  //	  std::cout<<"start Time of this deposit: "<<t0<<std::endl;
	  fEdep = energyDepositlist[i]->E();
	  fTotEdep += fEdep;
	  if (energyDepositlist[i]->PdgCode()==13){
	  fmuTotEdep += fEdep;
	  }
	  fEnergyDepositionVector.emplace_back(fEdep);
	  fTotGammaScint += energyDepositlist[i]->NumPhotons();
	  fTotEion += energyDepositlist[i]->NumElectrons();
	  
	}
      }else{
	mf::LogWarning("LightChargeAna") << "Cannot Find Deposited Energy. Failing";
	return;
      }
    }
    t0=minT;
    std::cout<<"Time the deposition starts: "<<t0<<std::endl;
    //SimChannel
    if(fSimCh){
      std::cout<<"Looking for total deposited energy in SimChannel"<<std::endl;
      art::Handle<std::vector<sim::SimChannel>> simChanHandle;
      evt.getByLabel(fSimulationProducerLabel, simChanHandle);
      // initialize at the start of every event
      dE = 0.;
      /*fDepVertexX = -1000;
      fDepVertexY = -1000;
      fDepVertexZ = -1000;
      //trackLength = 0.;
       goodDepVertex = false;
      unsigned int nTries = 0;
      ziFalse = -500.;
      */
    // find conversion point
      /*      while(!goodDepVertex && nTries < 300){
	  this->findDepositionVertex( (*simChanHandle) );
	  this->checkDepositionVertex( (*simChanHandle), xi, yi, zi );
	  nTries++;
	  }*/
      // std::cout << "Number of attemps at finding deposition vertex = " << nTries << "." << std::endl;
      // find energy deposited into first 2.5 cm of "track"
      // std::cout<<"Deposition vertex = energy deposited in the first 2.5 cm of the \"track\", found in x,y,z: "<<xi<<" "<<yi<<" "<<zi<<std::endl;
      //using TDCIDEs_t sim::SimChannel::fTDCIDEs = list of energy deposits for each TDC with signal
      //Definition at line 151 of file SimChannel.h.
      //Find total energy deposited in the collection plane 
      /* if(nTries < 300 && goodDepVertex){
	   fDepVertexX=xi;
	   fDepVertexY=yi;
	   fDepVertexZ=zi;*/
	  //	  fStartEdepHist->Fill(xi,yi,zi);
	   for(auto const& channel : (*simChanHandle)){
	     if(geom->SignalType(channel.Channel()) == geo::kCollection){
	       auto const& timeSlices = channel.TDCIDEMap();
	       for(auto const& t : timeSlices){
		 auto const& eDeps = t.second;
		 for(auto const& eDep : eDeps){
			  //			  fShowerEdepHist->Fill(eDep.x,eDep.y,eDep.z);
			  //trackLength = std::sqrt( (eDep.x-xi)*(eDep.x-xi) + (eDep.y-yi)*(eDep.y-yi) + (eDep.z-zi)*(eDep.z-zi) );
			  //std::cout<<"\"TrackLength\" so far: "<<trackLength<<std::endl;
			  //if(trackLength < 2.5 /* && eDep.z >= zi*/)
			  //{
		   dE += eDep.energy;
			      //  fTrackEdepHist->Fill(eDep.x,eDep.y,eDep.z);
			      //}
		 } // every energy deposition
	       } // in every time slice
	     } // for the collection plane
	   } // and for every wire
	   //std::cout<<"TrackLength in eDep steps: "<<trackLength<<std::endl;
	   fTotEdep = dE;
	   //dE /= 2.5;
	   //std::cout<<" --->  ~2.5cm of \"track\" gives dE/dx "<<dE<<std::endl;
	   //fMCdEdxVec.push_back(dE);
	   // if(dE > fMaxMCdEdx) fMaxMCdEdx = dE;
	   //if(dE < 0.5) fNullList.push_back(fEventID);
	   //}
      art::Handle< std::vector<simb::MCParticle> > particleHandle;
      evt.getByLabel(fSimulationProducerLabel, particleHandle);
      for(auto const& particle : (*particleHandle)){
	if(particle.Process() == "primary"){
	  Vxi = particle.Vx();
	  Vyi = particle.Vy();
	  Vzi = particle.Vz();
	}
      }
      fConvDist =  std::sqrt( (Vxi-xi)*(Vxi-xi) + (Vyi-yi)*(Vyi-yi) + (Vzi-zi)*(Vzi-zi) );
      //      fConvDistHist->Fill(fConvDist);
      std::cout<<"Conversion distance: "<<fConvDist<<std::endl;

      // Find the energy deposited on each plane in the TPC
      depositU=0;
      depositV=0;
      depositZ=0;
      // const std::vector<art::Ptr< sim::SimChannel > >& simChannels = backtracker->SimChannels();
      for(auto const& channelIt : (*simChanHandle)){
	//for (auto channelIt = simChannels.begin(); channelIt != simChannels.end(); ++channelIt) {
	int plane = geom->View(channelIt.Channel());
	for (auto const& tdcIt : channelIt.TDCIDEMap()) {
	  for (auto const& ideIt : tdcIt.second) {
	    // std::cout<<"Plane: "<<plane;
	    switch (plane) {
	    case geo::kU:
	      //std::cout<<" --> U "<<" dep en: "<<ideIt.energy<<std::endl;
	      depositU += ideIt.energy;
	      break;
	    case geo::kV:
	      depositV += ideIt.energy;
	      break;
	    case geo::kZ:
	      //std::cout<<" --> Z "<<" dep en: "<<ideIt.energy<<std::endl;
	      depositZ += ideIt.energy;
	      break;
	    }
	  }
	}
      }

      fDepositedU=depositU;
      fDepositedV=depositV;
      fDepositedZ=depositZ;
      std::cout<<"---> dE is (GeV): "<<fTotEdep/1000.<<" U,V,Z: "<<fDepositedU/1000.<<" "<<fDepositedU/1000.<<" "<<fDepositedZ/1000.<<std::endl;
    }//if SimCh
     
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
        mf::LogError("LightChargeAna") << "No MCTruth Particles";
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
      if (mctruth->Origin() == simb::kBeamNeutrino){
	fTrueCCNC  = mctruth->GetNeutrino().CCNC();
      }else{
	fTrueCCNC=1000;
      }

      std::cout<<"TrueE is: "<<fTrueE<<" Interaction CCNC is: "<<fTrueCCNC<<std::endl;
      
      fmuTrueE=0;

      // Get all the paricle including neutrino, and record its properties
      unsigned int const nParticles = mctruth->NParticles();
      for (unsigned int i = 0; i < nParticles; ++i) {
	simb::MCParticle const& particle = mctruth->GetParticle(i);
        fTruePxallpart    .emplace_back(particle.Px());
        fTruePyallpart    .emplace_back(particle.Py());
        fTruePzallpart    .emplace_back(particle.Pz());
        fTrueEallpart     .emplace_back(particle.E());
        fTrueAllPDG       .emplace_back(particle.PdgCode());
        if(particle.PdgCode()==13){
        fmuTrueE+=particle.E();
        
        }
      }

      // Get the PlaneID which describes the location of the true vertex
      int plane = 0;
      double loc[] = {part.Vx(), part.Vy(), part.Vz()};
      geo::TPCID tpc = geom->FindTPCAtPosition(loc);
      if (! geom->HasTPC(tpc) ) {
        mf::LogInfo("LightChargeAna") << "No valid TPC for " << tpc;
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
      mf::LogError("LightChargeAna") << "Event doesn't have signal ";
    }

    // Get the maximum possible time difference by getting number of ticks corresponding to
    // one full drift distance, and converting to time.
    double maxT=0;
    if(fBeam) maxT = clockData.TPCTick2Time(detProp.NumberTimeSamples());


    //////////////////////////////////////
    // CHARGE collected                 //
    ////////////////////////////////////// 
    //unsigned short wirePlane = geo::kZ;
    //std::cout<<"Wire Plane geo kZ is: "<<wirePlane<<std::endl;
    // Total charge collected in the event and time info of the hits
    fTotalCharge = 0.0;
    fTotalCorrCharge = 0.0;
    fTotalFinalCharge = 0.0;
    fmuTotalCorrCharge = 0.0;
    fmuTotalFinalCharge = 0.0;
    
    fADCSum = 0.0;
    fPeakAmplitude = 0.0;
    //fFirstHitTime = -1000;
    fHitDist = -1000;
    fMeanHitTime = -1000;
    //double firstHit=100000;
    double meanHit=0;
    //double lastHit=-100000;
    int collpl=0;
    //std::vector<double> hitXYZ;
    //std::vector<const sim::IDE*> ides;
    
    const double tauLifetime=detProp.ElectronLifetime();

    // here is the copy of associations to hits, based on original hit assns
    auto hitListHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitsLabel);
    std::cout<<"There are a total of: "<<hitListHandle->size()<<" hits in the event (collection+induction)"<<std::endl;
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

	/*	try{
	   ides= bt->HitToSimIDEs_Ps(clockData, hitPtr);
	}
	catch(...){}
	 if(ides.size()>0){
            hitXYZ = bt->SimIDEsToXYZ(ides);
	    std::cout<<"HitToXYZ: "<<hitXYZ[0]<<" "<<hitXYZ[2]<<" "<<hitXYZ[2]<<std::endl;
	    }*/
      }				      
    }

   std::vector<art::Ptr<recob::Hit>> eventHits=dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsLabel),2); //Plane_ID 2 is collection
   fTotalCorrCharge =145.137881*dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, eventHits);
   fTotalFinalCharge = (fTotalCorrCharge*Wion)/RecombFactor;
   std::cout<<"Electron lifetime: "<<tauLifetime<<" Using Recombination Factor: "<<RecombFactor<<" and Wion: "<<Wion<<" MeV/e"<<std::endl;
   std::cout<<"Reco Hits, Total Charge: "<<fTotalCharge<<" Total Charge corrected by electron lifetime: "<<fTotalCorrCharge<<" multiplied by Wion/R: "<<fTotalFinalCharge<<std::endl;
    
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
      std::cout << "Nb of hits in the collection plane: " << collpl <<" Mean hit time is " << fMeanHitTime <<" (us)" << " --> should correspond to charge travelling: "<<fHitDist<<" cm"<< '\n';
    
  }
    // Output the total charge
    std::cout << "Total Charge: " << fTotalCharge << '\n';
    
    ////////////////////////////////////////////////////
    //Space points of collection plane hits
    ////////////////////////////////////////////////
    int collsp=0;
    int indsp=0;
    int notknownsp=0;
    Float_t sppointx, sppointy, sppointz;
    Float_t musppointx, musppointy, musppointz;
    Float_t fvis_i=0;
    Float_t mufvis_i=0;
    Float_t totspcharge=0;
    Float_t mutotspcharge=0;
    //Float_t collspcharge=0;
    Fvis=0;
    muFvis=0;
    fmuTotalFinalCharge=0;
    fmuTotalCorrCharge=0;
    fE_QL=0;
    
    auto spHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fSpacePointLabel);
    art::FindManyP<recob::Hit> hitsFromSp(spHandle, evt, fSpacePointLabel);
    
    std::cout<<"Space points: "<<spHandle->size()<<std::endl;
    for (size_t i = 0; i < spHandle->size(); ++i){
      auto collhits = hitsFromSp.at(i);
      //      for (const auto & h : collhits) // find Collection hit, assume one is enough
      for(auto & h : collhits){
	if (h->SignalType() == geo::kCollection){ //
	  //||h->SignalType() == geo::kInduction){
	  //std::cout<<"This is a spcepoints-collection plane hit"<<std::endl;
	  // if(h->SignalType() == geo::kCollection){
	  collsp++;
	    //collspcharge += h->Integral();
	    // }
	  //if(h->SignalType() == geo::kInduction) indsp++;
	  art::Ptr<recob::SpacePoint> hitSP(spHandle, i);
	  sppointx=hitSP->XYZ()[0];
	  sppointy=hitSP->XYZ()[1];
	  sppointz=hitSP->XYZ()[2];
	  //fspointx    .emplace_back(sppointx);
          //fspointy    .emplace_back(sppointy);
          //fspointz    .emplace_back(sppointz);
	  float spcharge=h->Integral();
	  //fspointcharge.emplace_back(spcharge);
	  fvis_i=getfvis(sppointx, sppointy, sppointz, spcharge);
	  Fvis += fvis_i;
	  //std::cout<<"Position: "<<hitSP->XYZ()[0]<<" "<<hitSP->XYZ()[1]<<" "<<hitSP->XYZ()[2]<<std::endl;
	  totspcharge += spcharge;
	  //TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, h, true));
	  
	  //if (g4ID==13){
	//std::cout<<"Found: beccato!"<<std::endl;
	//}
	 TruthMatchUtils::G4ID g4ID = TruthMatchUtils::TrueParticleID(clockData, h, true);
	 if (TruthMatchUtils::Valid(g4ID)){
	 const simb::MCParticle *particle(pinv->TrackIdToParticle_P(g4ID));
	 if (particle->PdgCode()==13)
	{
	     musppointx=hitSP->XYZ()[0];
	     musppointy=hitSP->XYZ()[1];
	     musppointz=hitSP->XYZ()[2];
	     float muspcharge=h->Integral();
	     mufvis_i=getfvis(musppointx, musppointy, musppointz, muspcharge);
	     muFvis += mufvis_i;
	     mutotspcharge += muspcharge;
	     //fmuTotalCorrCharge += dune_ana::DUNEAnaHitUtils::LifetimeCorrection(clockData, detProp, h);
             //fmuTotalFinalCharge = (fmuTotalCorrCharge*Wion)/RecombFactor;
	  }
	  
	  }
	  }else if(h->SignalType() == geo::kInduction){
	    indsp++;
	  }else{
	  notknownsp++;
	}
	
	
      }
    }
    
    std::cout<<"Found: "<<collsp<<" spacepoints associated to collection hits, for a total charge of: "<<totspcharge<<std::endl;
    std::cout<<"Induction: "<<indsp<<" not known: "<<notknownsp<<" for a total of: "<<collsp+indsp+notknownsp<<std::endl;
    if(totspcharge!=0) {
      Fvis=Fvis/totspcharge;
    }else{
      Fvis=0;
    }
     std::cout<<"Fvis of the event: "<<Fvis<<std::endl; 
     
    if(mutotspcharge!=0) {
      muFvis=muFvis/mutotspcharge;
    }else{
      muFvis=0;
    }
     std::cout<<"muFvis of the event: "<<muFvis<<std::endl;  
     
     fLight=fSumPE/(Fvis*0.03*0.05);
     fE_QL=(fTotalCorrCharge+fLight)*Wph;
     
  //std::vector<art::Ptr<recob::Hit>> eventHits=dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsLabel),2); //Plane_ID 2 is collection
  // fTotalCorrCharge = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, eventHits);
  // fTotalFinalCharge = (fTotalCorrCharge*Wion)/RecombFactor;   
    
    const std::vector<art::Ptr<recob::PFParticle>> pfps = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fParticleModuleLabel);
    for (const art::Ptr<recob::PFParticle> pfp: pfps)
        {
	//const std::vector<art::Ptr<recob::Hit>> hits(dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fSpacePointLabel));
	//std::vector<art::Ptr<recob::Hit>> eventHits_1(dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fSpacePointLabel));
	//for (size_t i = 0; i < eventHits_1.size(); ++i){
	    //auto collhitsus = eventHits_1.at(i);
	    //TruthMatchUtils::G4ID g4IDus = TruthMatchUtils::TrueParticleID(clockData, collhitsus, true);
	    //std::cout<<"found hit"<< g4IDus <<std::endl;}
	std::vector<art::Ptr<recob::Hit>> eventHits_1=dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fSpacePointLabel),2); //Plane_ID 2 is collection
	TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, eventHits_1, true));
	//TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, eventHits_1, true));
	//std::cout<<"found hit"<<std::endl;
	if (TruthMatchUtils::Valid(g4ID))
	{
 //std::cout<<"found G4ID"<<std::endl;
	    const simb::MCParticle *particle(pinv->TrackIdToParticle_P(g4ID));
 //std::cout<<"it's "<<particle->PdgCode()<<std::endl;
 
	      if (particle->PdgCode()==13)
	{
	
		  fmuTotalCorrCharge +=145.137881*dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, eventHits_1);}
		  
		  }
	}
        
        fmuTotalFinalCharge = (fmuTotalCorrCharge*Wion)/RecombFactor;
	
	    //for(auto & h1 : eventHits_1){
	     
	    //std::vector<art::Ptr<recob::SpacePoint>> spHIT = dune_ana::DUNEAnaHitUtils::GetSpacePoints(h1, evt, fHitsLabel, fHitToSpacePointLabel);
	    //for (size_t i = 0; i < spHIT.size(); ++i){
	   // art::Ptr<recob::SpacePoint> split = spHIT.at(i);
	  //  musppointx=split->XYZ()[0];
	  //  musppointy=split->XYZ()[1];
	  //  musppointz=split->XYZ()[2];
	    //fmuspointx    .emplace_back(musppointx);
            //fmuspointy    .emplace_back(musppointy);
            //fmuspointz    .emplace_back(musppointz);
           // float muspcharge=h1->Integral();
	    //fmuspointcharge.emplace_back(muspcharge);
	  //  mufvis_i=getfvis(musppointx, musppointy, musppointz, muspcharge);
	  //  muFvis += mufvis_i;
	  //  mutotspcharge += muspcharge;
            
          //  }	
	
  //  } 
//}
//}
//else 
//{
    //std::cout<<"NOT FOUND G4ID"<<std::endl;
//}
	
	

//}

   fmuTotalFinalCharge = (fmuTotalCorrCharge*Wion)/RecombFactor; 

/* std::vector<art::Ptr<recob::Hit>> eventHits=dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitsLabel),2); //Plane_ID 2 is collection
   fTotalCorrCharge = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, eventHits);
   fTotalFinalCharge = (fTotalCorrCharge*Wion)/RecombFactor;
   std::cout<<"Electron lifetime: "<<tauLifetime<<" Using Recombination Factor: "<<RecombFactor<<" and Wion: "<<Wion<<" MeV/e"<<std::endl;
   std::cout<<"Reco Hits, Total Charge: "<<fTotalCharge<<" Total Charge corrected by electron lifetime: "<<fTotalCorrCharge<<" multiplied by Wion/R: "<<fTotalFinalCharge<<std::endl;*/
    
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
      std::cout << "Nb of hits in the collection plane: " << collpl <<" Mean hit time is " << fMeanHitTime <<" (us)" << " --> should correspond to charge travelling: "<<fHitDist<<" cm"<< '\n';
    
  }
    // Output the total charge
    std::cout << "Total Charge: " << fTotalCharge << '\n';
    
    //std::cout<<"Found: "<<collsp<<" spacepoints associated to collection hits, for a total charge of: "<<totspcharge<<std::endl;
    //std::cout<<"Induction: "<<indsp<<" not known: "<<notknownsp<<" for a total of: "<<collsp+indsp+notknownsp<<std::endl;
    //if(totspcharge!=0) {
     // Fvis=Fvis/totspcharge;
    //}else{
      //Fvis=0;
    //}
     //std::cout<<"Fvis of the event: "<<Fvis<<std::endl; 
      //This works, same nb of spacepoints I get with pandora
    /* auto spListHandle = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointLabel);
    std::cout<<"List of Event Space Points: "<<spListHandle->size()<<std::endl;
    //for(size_t isp = 0; isp < spListHandle->size(); ++isp){
    for(size_t isp = 0; isp < 5; ++isp){
      art::Ptr<recob::SpacePoint> spPtr(spListHandle, isp);
      std::cout<<"x,y,z: "<<spPtr->XYZ()[0]<<" "<<spPtr->XYZ()[0]<<" "<<spPtr->XYZ()[0]<<std::endl;
      }
    */
    
    
    /////////////////////////
    // Analyze the flashes //
    /////////////////////////

    // Set up some flags to fill as we loop
    // through flashes. These will control
    // filling of efficiency plots after the loop.

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
	  mf::LogError("LightChargeAna") << "Skipping Flash: not w/in 1 drift window (timeDiff < -10 OR timeDiff > maxT)" <<
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
      //fPESum       += fTotalPE;
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

    }

    ///////////////////////////////////////
    //Reco tracks and showers            //
    ///////////////////////////////////////

    const std::vector<art::Ptr<recob::PFParticle>> particles = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fParticleModuleLabel);
    std::cout<<"**** check Pandora: "<<std::endl;
    int cPart=0;
    fNShowers=0;
    fShwLength=0;
    fShwOpenAngle=0;
    fShwEn=0;
    fShwdEdx=0;
    fShwStartX=-9999;
    fShwStartY=-9999;
    fShwStartZ=-9999;
    double maxE=0;
    double maxdedx=0;
    double maxL=0;
    double maxA=0;
    int getshwvtx=0;
    unsigned int bp=999;
    int bestplane=999;

    //lar_pandora::PFParticlesToVertices particlesToVertices;
    art::Ptr<recob::SpacePoint> spacePoint;
    int totpfsp=0;
    
    for (const art::Ptr<recob::PFParticle> &particle : particles){
      // int thisShwID = -999;
      cPart++;
      //check Pandora vertex
      /*    if(dune_ana::DUNEAnaPFParticleUtils::GetVertex(particle,evt,fParticleModuleLabel)->isValid()){
	const art::Ptr<recob::Vertex> vtx = dune_ana::DUNEAnaPFParticleUtils::GetVertex(particle,evt,fParticleModuleLabel);
	recob::Track::Point_t vtxpos =vtx->position();
	fPandoraVtxX=vtxpos.X();
	fPandoraVtxY=vtxpos.Y();
	fPandoraVtxZ=vtxpos.Z();
	std::cout<<"Particle: "<<cPart<<" Pandora Vertex X Y Z: "<<fPandoraVtxX<<" "<<fPandoraVtxY<<" "<<fPandoraVtxZ<<std::endl;
      }else{
	std::cout<<"Vertex status is not valid"<<std::endl;
	}*/
      
      //check pandora vertex --> Sempre emtpy!!
      /*      short iv = 0;
      if(particle->IsPrimary()){
	std::cout<<"PFParticle "<<cPart<<" is primary"<<std::endl;
      }else{
	std::cout<<"PFParticle "<<cPart<<" is not primary"<<std::endl;
      }
      lar_pandora::PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
      // if (particlesToVertices.end() != vIter) {
	const lar_pandora::VertexVector &vertexVector = vIter->second;
	if (!vertexVector.empty()) {
	  if (vertexVector.size() == 1) {
	    const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
	    double xyz[3] = {0.0, 0.0, 0.0} ;
	    vertex->XYZ(xyz);
	    std::cout<<"There's one vertex, PFParticle to vertices: "<<iv<<") X:"<<xyz[0]<<" Y: "<<xyz[1]<<" Z: "<<xyz[2]<<std::endl;
	    iv++;
	  }else{
	    std::cout<<"There's more than one vertex: "<< vertexVector.size() << std::endl;
	  }
	}else{
	  std::cout<<"Vector of LAr Pandora vertices is empty"<<std::endl;
	}
	//}
	*/ 
      //Check Spacepoints
      	std::vector<art::Ptr<recob::SpacePoint>> spacePoints= dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(particle,evt,fParticleModuleLabel);
	std::cout<<"Get Spacepoints from Pandora"<<std::endl;
   	if(spacePoints.size()>0){
	  std::cout<<"Found "<<spacePoints.size()<<" space points"<<std::endl;
	  totpfsp=totpfsp+spacePoints.size();
	  //	  for (unsigned int iSpacePoint = 0; iSpacePoint < spacePoints.size(); ++iSpacePoint){
	    for (unsigned int iSpacePoint = 0; iSpacePoint < 5; ++iSpacePoint){
	    spacePoint = spacePoints[iSpacePoint];
	    /*sppointx=spacePoint->XYZ()[0];
	    sppointy=spacePoint->XYZ()[1];
	    sppointz=spacePoint->XYZ()[2];
	    std::cout<<"The space point is in x,y,z: "<<sppointx<<" "<<sppointy<<" "<<sppointz<<std::endl;*/
	  }//loop on spacepoints collection
	}//if spacepoints
	else{
	  std::cout<<"space point vector empty"<<std::endl;
	}
      
      //check reco shower
      if (dune_ana::DUNEAnaPFParticleUtils::IsShower(particle,evt,fParticleModuleLabel,fShowerLabel)){
	const art::Ptr<recob::Shower> shw = dune_ana::DUNEAnaPFParticleUtils::GetShower(particle,evt,fParticleModuleLabel,fShowerLabel);
	fNShowers++;
	std::cout<<"Particle: "<<cPart<<" Shower: "<<fNShowers<<std::endl;
	std::cout<<" Shower length: "<<shw->Length()<<" [cm], and Cone aperture: "<<shw->OpenAngle()<<" [rad]"<<std::endl;
	std::cout<<" Shower Energy at 2 is: "<<shw->Energy().at(2)<<std::endl;
	std::cout<<"Best plane to evaluate energy is: "<<shw->best_plane()<<std::endl;
	std::cout<<"Calculated Energy and dEdx per each plane: "<<std::endl;
	for(unsigned int i=0; i<shw->Energy().size(); i++){
	  std::cout << shw->Energy().at(i)<<" "<<shw->dEdx().at(i)<<std::endl;
	  fShwEnVector.emplace_back(shw->Energy().at(i));
	  fShwdEdxVector.emplace_back(shw->dEdx().at(i));
	}	
	if(shw->Energy().at(2)>maxE){
	  maxE=shw->Energy().at(2);
	  maxdedx=shw->dEdx().at(2);
	  maxL=shw->Length();
	  maxA=shw->OpenAngle();
	  getshwvtx=cPart;
	  bp=shw->best_plane();
	}
	
      }else{
	mf::LogWarning("LightChargeAna") << "Didn't Find any reco shower in this event";
	/*	fShwLength=0;
	fShwOpenAngle=0;
	fShwEn=0;
	fShwdEdx=0;*/
	//return;
      }
     
      //Check reco tracks
      if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(particle,evt,fParticleModuleLabel,fTrackLabel)){
	const art::Ptr<recob::Track> trk = dune_ana::DUNEAnaPFParticleUtils::GetTrack(particle,evt,fParticleModuleLabel,fTrackLabel);
	std::cout<<" RECONSTRUCTED TRACK: Particle: "<<cPart<<" Track ID: "<<trk->ID()<<" Particle ID: "<<trk->ParticleId();
	TVector3 trkvtx = trk->Vertex<TVector3>(); 
	std::cout<<" Track Vertex: "<<trkvtx(0)<<" "<<trkvtx(1)<<" "<<trkvtx(2)<<std::endl;//" Length: "<<trk->Length<<std::endl;
   	
      }else{
	//std::cout<<"Didn't Find any reco track in this event"<<std::endl;
	mf::LogWarning("LightChargeAna") << "Didn't Find any reco track in this event";
      }
      
    }//PFParticles

    //Take best shower and Check reco vertex for saved shower
      if(fNShowers==0){
	std::cout<<"This event doesn't have any reco Shower"<<std::endl;
      }else{
	std::cout<<"Check Pandora Vertex for selected shower (part nb "<<getshwvtx<<")"<<std::endl;
	int c=0;
	for (const art::Ptr<recob::PFParticle> &part : particles){
	  c++;
	  if(c==getshwvtx){
	    const art::Ptr<recob::Shower> bestshw = dune_ana::DUNEAnaPFParticleUtils::GetShower(part,evt,fParticleModuleLabel,fShowerLabel);
	    //std::cout<<" Saving info for plane nb "<<bp<<std::endl;
	    /*	for(unsigned int i=0; i<bestshw->Energy().size(); i++){
		if(i==bp){
		fShwEn=bestshw->Energy().at(i);
		fShwdEdx=bestshw->dEdx().at(i);	
		fShwLength=bestshw->Length();
		fShwOpenAngle=bestshw->OpenAngle();
		shwstart = bestshw->ShowerStart();
		bestplane=bestshw->best_plane();
		}
		}*/
	    std::cout<<" Saving info for plane nb 2 "<<std::endl;
	    fShwEn=bestshw->Energy().at(2);
	    fShwdEdx=bestshw->dEdx().at(2);	
	    fShwLength=bestshw->Length();
	    fShwOpenAngle=bestshw->OpenAngle();
	    shwstart = bestshw->ShowerStart();
	    //bestplane=bestshw->best_plane();
	    bestplane=bp;
	    
	    if(dune_ana::DUNEAnaPFParticleUtils::GetVertex(part,evt,fParticleModuleLabel)->isValid()){
	      const art::Ptr<recob::Vertex> shwvtx = dune_ana::DUNEAnaPFParticleUtils::GetVertex(part,evt,fParticleModuleLabel);
	      recob::Track::Point_t shwvtxpos =shwvtx->position();
	      fPandoraVtxX=shwvtxpos.X();
	      fPandoraVtxY=shwvtxpos.Y();
	      fPandoraVtxZ=shwvtxpos.Z();
	      std::cout<<" Vertex X Y Z: "<<fPandoraVtxX<<" "<<fPandoraVtxY<<" "<<fPandoraVtxZ<<std::endl;
	    }else{
	      std::cout<<"Vertex status is not valid"<<std::endl;
	    }
	  }//c
	}//for part of particles
    
    
	std::cout<<"Pandora has reconstructed: "<<fNShowers<<" showers, saving most energetic shower (part nb "<<getshwvtx<<") info of plane 2,and best plane is (nb "<<bestplane<<")"<<std::endl;
	//Take only shower of best plane
	fShwStartX=shwstart(0);
	fShwStartY=shwstart(1);
	fShwStartZ=shwstart(2);
	std::cout<<"Shower Start: "<<fShwStartX<<" "<<fShwStartY<<" "<<fShwStartZ<<std::endl;
	std::cout<<"Energy, dEdx: "<<fShwEn<<" "<<fShwdEdx<<" Length: "<<fShwLength<<" Opening Angle: "<<fShwOpenAngle<<std::endl;
	std::cout<<"While most energetic plane of most energetic shower has En, dEdx, L and Angle: "<<maxE<<" "<<maxdedx<<" "<<maxL<<" "<<maxA<<std::endl;
	std::cout<<"Pandora Vertex for best shower: "<<fPandoraVtxX<<" "<<fPandoraVtxY<<" "<<fPandoraVtxZ<<std::endl;
      }//there are showers

      std::cout<<"Total Space points from PFParticles: "<<totpfsp<<std::endl;

      
    std::cout<<"*************** This Event has: "<<std::endl;
    std::cout<<"Total PE reco: "<< fSumPE <<" with total area: "<<fOpHitTotalArea<<std::endl;
    std::cout<<"Total Charge: "<< fTotalCharge <<std::endl;
    //std::cout<<"Reco  Shower(s) Energy per plane:";
    //for(unsigned i=0; i<fShwEnVector.size(); i++) std::cout<<" "<<fShwEnVector.at(i);
    //std::cout<<std::endl;
    std::cout<<"*************** End "<<std::endl;
    
    //contenimento mu
    
    fmuPandoraContained50=0;
    
    //double lengthvect[10];
    const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fParticleModuleLabel); 
    int k=0;
    for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){

  
    if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp,evt,fParticleModuleLabel,fTrackLabel)){
      
      art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp,evt,fParticleModuleLabel, fTrackLabel); 
      
      if(pfp->PdgCode()==13){
      /*std::cout<<"start x"<< track->Start().X() <<std::endl;
       std::cout<<"start y"<< track->Start().Y() <<std::endl;
        std::cout<<"start z"<< track->Start().Z() <<std::endl;
        std::cout<<"end x"<< track->End().X() <<std::endl;
       std::cout<<"end y"<< track->End().Y() <<std::endl;
        std::cout<<"end z"<< track->End().Z() <<std::endl;*/
        
        
        
      if ((fabs(track->Start().X()))<310&&(fabs(track->Start().Y()))<550&&(fabs(track->Start().Z()))>50&&(fabs(track->Start().Z()))<1340&&(fabs(track->End().X()))<310
      &&(fabs(track->End().Y()))<550&&(fabs(track->End().Z()))>50&&(fabs(track->End().Z())<1340)){
      fmuPandoraContained50=1;
      //lengthvect[k]=track->Length();
      k++;
      }
      else{
      fmuPandoraContained50=0;
      break;
      }
      }
      
      }
    
    }
    
    //contenimento totale
    
    //fPandoraContained50=0;
    
    //double lengthvect[10];
    //const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect2 = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fParticleModuleLabel); 
    //int k2=0;
    //for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect2){

  
    //if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp,evt,fParticleModuleLabel,fTrackLabel)){
      
      //art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp,evt,fParticleModuleLabel, fTrackLabel); 
      
      
      /*std::cout<<"start x"<< track->Start().X() <<std::endl;
       std::cout<<"start y"<< track->Start().Y() <<std::endl;
        std::cout<<"start z"<< track->Start().Z() <<std::endl;
        std::cout<<"end x"<< track->End().X() <<std::endl;
       std::cout<<"end y"<< track->End().Y() <<std::endl;
        std::cout<<"end z"<< track->End().Z() <<std::endl;*/
        
      fContained50=0;  
      int nSimEnergyDeposits2 = 0;
      art::Handle< std::vector<sim::SimEnergyDeposit> > energyDepositHandle2;
      std::vector<art::Ptr<sim::SimEnergyDeposit> > energyDepositlist2;
      if(evt.getByLabel(fSimLabel, energyDepositHandle2)){
	art::fill_ptr_vector(energyDepositlist2, energyDepositHandle2);
	nSimEnergyDeposits2 = energyDepositlist2.size();
      for(int i = 0; i < nSimEnergyDeposits2; i++){
	  
      if ((fabs(energyDepositlist2[i]->MidPointX()))<310&&(fabs(energyDepositlist2[i]->MidPointY()))<550&&(fabs(energyDepositlist2[i]->MidPointZ()))>50&&(fabs(energyDepositlist2[i]->MidPointZ()))<1340){
      fContained50=1;
    
      }
      
      else{
      fContained50=0;
      break;
      }
      
      
      }
    }
    
    
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
    fOpHitPEVector              .clear();
    fOpHitPeakTimeVector        .clear();
    fOpHitPeakTimeAbsVector     .clear();
    fOpHitFastToTotalVector     .clear();
    fShwEnVector                .clear();
    fShwdEdxVector              .clear();
    //fspointx.clear();
    //fspointy.clear();
    //fspointz.clear();
    //fspointcharge.clear();
    //fmuspointx.clear();
    //fmuspointy.clear();
    //fmuspointz.clear();
    //fmuspointcharge.clear();
    //fHitPeakTimeTicks           .clear();
    //fPESum=0;
    
  }

void LightChargeAna::endJob()
  {
    /*    fMCdEdxHist = tfs->make<TH1D>("fMCdEdxHist",";dE/dx (MeV/cm)",200,0.0,10.);

    for(unsigned int e = 0; e < fMCdEdxVec.size(); e++)
    {
      fMCdEdxHist->Fill(fMCdEdxVec[e]);
    }
    */
    //    std::cout << "Events with < 0.5 dE/dx entry: ";
    // for(unsigned int a = 0; a < fNullList.size(); a++) std::cout << fNullList[a] << ", ";

  }
  
void LightChargeAna::readLightMap(){
    int linenumber=0;
    //std::string tempx, tempy, tempz, tempw;
    //std::string line;
    float  tempx, tempy, tempz, tempw;
    std::ifstream infile;
    infile.open("mappa30000.txt",std::ios::in);
    // infile.open("pippo.txt",std::ifstream::in);
    if(infile.is_open()){
      //      while (std::getline(infile,line)){
      while( infile >> tempx >> tempy >> tempz >> tempw){
	linenumber++;
	//	  std::stringstream ss(line);
	//	  ss >> std::ws >> tempx >> tempy >> tempz >> tempw;
	try{
	  
	  Bin_x.push_back(tempx);
	  Bin_y.push_back(tempy);
	  Bin_z.push_back(tempz);
	  Bin_w.push_back(tempw);
	  /*	  	  Bin_x.push_back(stod(tempx));
	  Bin_y.push_back(stod(tempy));
	  Bin_z.push_back(stod(tempz));
	  Bin_w.push_back(stod(tempw));*/
	}catch (const std::exception& ex){
	  std::cout << "Line # " << linenumber << " is invalid, skipping..." << std::endl;
	}
      }
      infile.close();
    }else{
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		<< "NO FILE FOUND FOR LIGHT MAP\n"
		<< "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

    }
    std::cout<<" filled "<<Bin_w.size()<<" light points from "<<linenumber<<" lines"<<std::endl;
   
  }

  float LightChargeAna::getfvis(float x, float y, float z, float q){
    if(q==0) return 0;

    //PANDORA offset
    //if(x<0) x=x-43.45;
    //if(x>0) x=x+43.45;
    if(x<0) x=x-0;
    if(x>0) x=x+0;
    
    // std::cout<<"Get expected light for position: "<<x<<" "<<y<<" "<<z<<std::endl;
    //BIN SIZE is x=20, y=31.429 , z=38.5714
    float stepx=10;
    float stepy=31.429/2;
    float stepz=38.5714/2;
    bool bin=false;
    float xmin,xmax,ymin,ymax,zmin,zmax;
    float thisf=0;
    float vis=0;
    int mapsize=Bin_w.size();
    //   std::cout<<" Loop over "<<mapsize<<" light points"<<std::endl;
    for (int f = 0; f<mapsize; f++){
    //  for (int f = 1; f<5; f++){  
      //  if (Bin_w.at(f) != 0) {
      xmin=Bin_x.at(f)-stepx;
      xmax=Bin_x.at(f)+stepx;
      ymin=Bin_y.at(f)-stepy;
      ymax=Bin_y.at(f)+stepy;
      zmin=Bin_z.at(f)-stepz;
      zmax=Bin_z.at(f)+stepz;
      if(xmin <= x && x < xmax && ymin <= y && y < ymax && zmin <= z && z < zmax ){
	//std::cout<<"Light map bin: "<<Bin_x.at(f)<<" "<<Bin_y.at(f)<<" "<<Bin_z.at(f)<<" "<<Bin_w.at(f)<<std::endl;
	thisf=Bin_w.at(f);
	bin=true;
      }
    }
  
    if(!bin){
      //      std::cout<<"didn't find a bin in the light map for position: "<<x<<" "<<y<<" "<<z<<std::endl;
      thisf=0;
    }
    //    std::cout<<" ---> f is: "<<thisf<<" and q "<<q;
    vis=q*thisf;
    //std::cout<<" ---> vis: "<<vis<<std::endl;
    return vis;
  }  

 }// namespace opdet
  


namespace opdet {
  DEFINE_ART_MODULE(LightChargeAna)
}
