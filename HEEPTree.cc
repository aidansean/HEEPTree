// -*- C++ -*-
//
// Package:    HEEPTree
// Class:      HEEPTree
// 
/**\class HEEPTree HEEPTree.cc RecoEgamma/HEEPTree/src/HEEPTree.cc

Description: <one line class summary>
Implementation:
<Notes on implementation>
*/
//
// Original Author:  Charaf Otman
//         Created:  Thu Jan 17 14:41:56 CET 2008
//
//Cleaning ladies : Aidan
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "UserCode/HEEPSkims/interface/BranchWrapper.h"
#include "UserCode/HEEPSkims/interface/HEEPTree.h"
#include "UserCode/HEEPSkims/interface/EtSort.h"

#include <iostream>
#include <TMath.h>
#include <vector>

#define PI 3.141592654
#define TWOPI 6.283185308

using namespace std ;
using namespace reco;
using namespace edm ;

enum variableTypes{
  kBool,
  kDouble,
  kFloat,
  kInt,
  kVectorBool,
  kVectorDouble,
  kVectorFloat,
  kVectorInt,
  kVectorVectorBool,
  kVectorVectorDouble,
  kVectorVectorFloat,
  kVectorVectorInt
};

HEEPTree::HEEPTree(const edm::ParameterSet& iConfig){
  current_var_type = -1 ;
  debug = true ;
  
  beamSpotLabel_ = consumes<BeamSpot>(iConfig.getParameter<InputTag>("beamSpot"));
  ScPtMin_       = iConfig.getUntrackedParameter<double>("ScPtMin"      , 10.);
  GsfPtMin_      = iConfig.getUntrackedParameter<double>("GsfPtMin"     ,  5.);
  GsfTrackPtMin_ = iConfig.getUntrackedParameter<double>("GsfTrackPtMin",  5.);
  muPtMin_       = iConfig.getUntrackedParameter<double>("muPtMin"      ,  5.);
  
  EcalHcal1EffAreaEndcaps_ = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaEndcaps", 0.) ;
  EcalHcal1EffAreaBarrel_  = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaBarrel" , 0.) ;
}


HEEPTree::~HEEPTree(){}

bool HEEPTree::branch_exists(std::string name){
  for(unsigned int i=0 ; i<all_vars_.size() ; i++){
    if(all_vars_.at(i)->name()==name) return true ;
  }
  return false ;
}
bool HEEPTree::add_branch(std::string name){ return add_branch(name, current_var_type) ; }

void HEEPTree::set_branch_type(int type){ current_var_type = type ; }
int  HEEPTree::get_branch_type(){ return current_var_type ; }

bool HEEPTree::add_branch(std::string name, int type){
  // First check to see if this branch name has already been used
  bool success = !(branch_exists(name)) ;
  if(success==false){
    return false ;
  }
  list_of_branches.push_back(std::pair<std::string,int>(name,type)) ;
  switch(type){
    case kBool:{
      BranchWrapperB*  bw = new BranchWrapperB(name) ;
      vars_B_  .push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kDouble:{
      BranchWrapperD*  bw = new BranchWrapperD(name) ;
      vars_D_  .push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
    }
    case kFloat:{
      BranchWrapperF*  bw = new BranchWrapperF(name) ;
      vars_F_  .push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kInt:{
      BranchWrapperI*  bw = new BranchWrapperI(name) ;
      vars_I_  .push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorBool:{
      BranchWrapperBV* bw = new BranchWrapperBV(name) ;
      vars_BV_ .push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorDouble:{
      BranchWrapperDV* bw = new BranchWrapperDV(name) ;
      vars_DV_ .push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorFloat:{
      BranchWrapperFV* bw = new BranchWrapperFV(name) ;
      vars_FV_ .push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorInt:{
      BranchWrapperIV* bw = new BranchWrapperIV(name) ;
      vars_IV_ .push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorVectorBool:{
      BranchWrapperBVV* bw = new BranchWrapperBVV(name) ;
      vars_BVV_.push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorVectorDouble:{
      BranchWrapperDVV* bw = new BranchWrapperDVV(name) ;
      vars_DVV_.push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorVectorFloat:{
      BranchWrapperFVV* bw = new BranchWrapperFVV(name) ;
      vars_FVV_.push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorVectorInt:{
      BranchWrapperIVV* bw = new BranchWrapperIVV(name) ;
      vars_IVV_.push_back(bw) ;
      all_vars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    default : return false ; // Bail out if we don't know the type of branch
  }
  
  return true ;
}

// ------------ method called once each job just before starting event loop  ------------
void HEEPTree::beginJob(){
  edm::Service<TFileService> fs;
  myFile = new TFile("outfile.root", "RECREATE") ;
  mytree = new TTree("HEEPTree", "HEEPTree") ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                  Event information                                 //
  ////////////////////////////////////////////////////////////////////////////////////////
  set_branch_type(kInt) ;
  add_branch("eventNumber"    ) ;
  add_branch("runNumber"      ) ;
  add_branch("luminosityBlock") ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                  Primary vertices                                  //
  ////////////////////////////////////////////////////////////////////////////////////////
  add_branch("pv_n", kInt) ;
  set_branch_type(kVectorFloat) ;
  add_branch("pv_x") ;
  add_branch("pv_y") ;
  add_branch("pv_z") ;
  add_branch("pv_isValid", kVectorBool) ;
  add_branch("pv_normChi2") ;
  set_branch_type(kVectorInt) ;
  add_branch("pv_ndof") ;
  add_branch("pv_nTracks") ;
  add_branch("pv_totTrackSize") ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                   Superclusters                                    //
  ////////////////////////////////////////////////////////////////////////////////////////
  set_branch_type(kVectorFloat) ;
  add_branch("sc_energy") ;
  add_branch("sc_eta") ;
  add_branch("sc_etacorr") ;
  add_branch("sc_theta") ;
  add_branch("sc_thetacorr") ;
  add_branch("sc_et") ;
  add_branch("sc_phi") ;
  add_branch("sc_px") ;
  add_branch("sc_py") ;
  add_branch("sc_pz") ;
  add_branch("sc_x") ;
  add_branch("sc_y") ;
  add_branch("sc_z") ;
  add_branch("sc_indexforgsf", kVectorInt) ;
  add_branch("scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", kVectorBool) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                     GSF tracks                                     //
  ////////////////////////////////////////////////////////////////////////////////////////
  add_branch("gsftrack_n", kInt) ;
  set_branch_type(kVectorFloat) ;
  add_branch("gsftrack_eta") ;
  add_branch("gsftrack_phi") ;
  add_branch("gsftrack_p"  ) ;
  add_branch("gsftrack_pt" ) ;
  add_branch("gsftrack_px" ) ;
  add_branch("gsftrack_py" ) ;
  add_branch("gsftrack_pz" ) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                               GSF Electron collection                              //
  ////////////////////////////////////////////////////////////////////////////////////////
  set_branch_type(kVectorVectorInt) ;
  add_branch("gsf_crystal_ietaorix") ;
  add_branch("gsf_crystal_iphioriy") ;
  set_branch_type(kVectorVectorFloat) ;
  add_branch("gsf_crystal_energy") ;
  add_branch("gsf_crystal_eta") ;
  
  set_branch_type(kVectorBool) ;
  add_branch("gsf_isEB") ;
  add_branch("gsf_isEE") ;
  set_branch_type(kVectorFloat) ;
  add_branch("gsf_theta") ;
  add_branch("gsf_deltaEtaATcalo") ;
  add_branch("gsf_deltaPhiATcalo") ;
  add_branch("gsf_sigmaetaeta") ;
  add_branch("gsf_sigmaIetaIeta") ;
  add_branch("gsf_ecalEnergy") ;
  add_branch("gsf_eOVERp") ;
  add_branch("gsf_dxy") ;
  add_branch("gsf_dxy_beamSpot") ;
  add_branch("gsf_dxy_firstPVtx") ;
  add_branch("gsf_dxy_firstPVtxwithBS") ;
  add_branch("gsf_dxyError") ;
  add_branch("gsf_dz") ;
  add_branch("gsf_dz_beamSpot") ;
  add_branch("gsf_dz_firstPVtx") ;
  add_branch("gsf_dz_firstPVtxwithBS") ;
  add_branch("gsf_dzError") ;
  add_branch("gsf_vz") ;
  set_branch_type(kVectorInt) ;
  add_branch("gsf_nHits") ;
  add_branch("gsf_nLostInnerHits") ;
  add_branch("gsf_nLostOuterHits") ;
  add_branch("gsf_convFlags") ;
  set_branch_type(kVectorFloat) ;
  add_branch("gsf_convDist") ;
  add_branch("gsf_convDcot") ;
  add_branch("gsf_convRadius") ;
  add_branch("gsf_fBrem") ;
  add_branch("gsf_e1x5") ;
  add_branch("gsf_e2x5") ;
  add_branch("gsf_e5x5") ;
  add_branch("gsf_e1x3") ;
  add_branch("gsf_p") ;
  add_branch("gsf_e") ;
  add_branch("gsf_deltaeta") ;
  add_branch("gsf_deltaphi") ;
  add_branch("gsf_hovere") ;
  add_branch("gsf_hdepth1overe") ;
  add_branch("gsf_hdepth2overe") ;
  add_branch("gsf_trackiso") ;
  add_branch("gsf_ecaliso") ;
  add_branch("gsf_hcaliso1") ;
  add_branch("gsf_hcaliso2") ;
  add_branch("gsf_isecaldriven"   , kVectorBool) ;
  add_branch("gsf_istrackerdriven", kVectorBool) ;
  add_branch("gsf_hitsinfo", kVectorVectorInt) ;
  
  set_branch_type(kVectorBool) ;
  add_branch("gsfMatch_hltL1sL1SingleEG12") ;
  add_branch("gsfMatch_hltL1sL1Mu3p5EG12") ;
  add_branch("gsfMatch_hltL1sL1SingleEG22") ;
  add_branch("gsfMatch_hltEle33CaloIdLPixelMatchFilter") ;
  add_branch("gsfMatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter") ;
  add_branch("gsfMatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter") ;
  add_branch("gsfMatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter") ;
  add_branch("gsfMatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter") ;
  add_branch("gsfMatch_hltEle27WP80TrackIsoFilter") ;
  add_branch("gsfMatch_hltMu22Photon22CaloIdLHEFilter") ;
  
  set_branch_type(kVectorFloat) ;
  add_branch("gsfsc_e") ;
  add_branch("gsfsc_pt") ;
  add_branch("gsfsc_eta") ;
  add_branch("gsfsc_phi") ;
  add_branch("gsfsc_px") ;
  add_branch("gsfsc_py") ;
  add_branch("gsfsc_pz") ;
  add_branch("gsf_e2x5overe5x5") ;
  add_branch("gsf_e1x5overe5x5") ;
  add_branch("gsf_gsfet") ;
  
  set_branch_type(kVectorBool) ;
  add_branch("gsfpass_ET") ;
  add_branch("gsfpass_PT") ;
  add_branch("gsfpass_DETETA") ;
  add_branch("gsfpass_CRACK") ;
  add_branch("gsfpass_DETAIN") ;
  add_branch("gsfpass_DPHIIN") ;
  add_branch("gsfpass_HADEM") ;
  add_branch("gsfpass_SIGMAIETAIETA") ;
  add_branch("gsfpass_E2X5OVER5X5") ;
  add_branch("gsfpass_ISOLEMHADDEPTH1") ;
  add_branch("gsfpass_ISOLHADDEPTH2") ;
  add_branch("gsfpass_ISOLPTTRKS") ;
  add_branch("gsfpass_ECALDRIVEN") ;
  add_branch("gsfpass_INVALID") ;
  add_branch("gsfpass_NOMISSINGHITS") ;
  add_branch("gsfpass_NOCONVERSION") ;
  add_branch("gsfpass_DXYFIRSTPV") ;
  add_branch("gsfpass_HEEP") ;
  add_branch("gsfpass_ID") ;
  add_branch("gsfpass_ISO") ;

  //charge information
  set_branch_type(kVectorInt) ;
  add_branch("sc_pixcharge") ;
  add_branch("ctf_charge") ;
  add_branch("gsf_charge") ;
  add_branch("gsf_class") ;
  
  set_branch_type(kVectorBool) ;
  add_branch("gsf_ctfscpixconsistent") ;
  add_branch("gsf_scpixconsistent") ;
  add_branch("gsf_ctfconsistent") ;

  //conversion information
  set_branch_type(kVectorFloat) ;
  add_branch("conv_vtxProb") ;
  add_branch("conv_lxy") ;
  set_branch_type(kVectorInt) ;
  add_branch("conv_nHitsMax") ;
  add_branch("conv_eleind") ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                             GED GSF Electron collection                            //
  ////////////////////////////////////////////////////////////////////////////////////////
  add_branch("gsf_n", kInt) ;
  set_branch_type(kVectorFloat) ;
  add_branch("gsf_pt" ) ;
  add_branch("gsf_px" ) ;
  add_branch("gsf_py" ) ;
  add_branch("gsf_pz" ) ;
  add_branch("gsf_E"  ) ;
  add_branch("gsf_eta") ;
  add_branch("gsf_phi") ;
  add_branch("gsf_charge", kVectorInt) ;
  add_branch("gsf_hovere"  ) ;
  add_branch("gsf_preshowerEnergy") ;
  add_branch("gsf_eseffsixix") ;
  add_branch("gsf_eseffsiyiy") ;
  add_branch("gsf_eseffsirir") ;
  set_branch_type(kVectorVectorFloat) ;
  add_branch("gsf_eshitsixix") ;
  add_branch("gsf_eshitsiyiy") ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                  Muon collection                                   //
  ////////////////////////////////////////////////////////////////////////////////////////
  add_branch("muon_n", kInt) ;
  set_branch_type(kVectorFloat) ;
  add_branch("muon_pt") ;
  add_branch("muon_ptError") ;
  add_branch("muon_gTrk_pt") ;
  add_branch("muon_gTrk_ptError") ;
  add_branch("muon_eta") ;
  add_branch("muon_etaError") ;
  add_branch("muon_phi") ;
  add_branch("muon_phiError") ;
  add_branch("muon_theta") ;
  add_branch("muon_thetaError") ;
  add_branch("muon_outerPt") ;
  add_branch("muon_outerEta") ;
  add_branch("muon_outerPhi") ;
  add_branch("muon_outerTheta") ;
  add_branch("muon_px") ;
  add_branch("muon_py") ;
  add_branch("muon_pz") ;
  set_branch_type(kVectorInt) ;
  add_branch("muon_charge") ;
  add_branch("muon_nhitspixel") ;
  add_branch("muon_nhitstrack") ;
  add_branch("muon_nhitsmuons") ;
  add_branch("muon_nhitstotal") ;
  add_branch("muon_nlayerswithhits") ;
  add_branch("muon_nlosthits") ;
  add_branch("muon_nSegmentMatch") ;
  set_branch_type(kVectorBool) ;
  add_branch("muon_isTrackerMuon") ;
  add_branch("muon_isPFMuon") ;
  add_branch("muon_isPFIsolationValid") ;
  add_branch("muon_chi2"    , kVectorFloat) ;
  add_branch("muon_ndof"    , kVectorFloat) ;
  add_branch("muon_normChi2", kVectorFloat) ;
  set_branch_type(kVectorFloat) ;
  add_branch("muon_d0") ;
  add_branch("muon_d0Error") ;
  add_branch("muon_dz_cmsCenter") ;
  add_branch("muon_dz_beamSpot") ;
  add_branch("muon_dz_firstPVtx") ;
  add_branch("muon_dz_firstPVtxwithBS") ;
  add_branch("muon_dzError") ;
  add_branch("muon_dxy_cmsCenter") ;
  add_branch("muon_dxy_beamSpot") ;
  add_branch("muon_dxy_firstPVtx") ;
  add_branch("muon_dxy_firstPVtxwithBS") ;
  add_branch("muon_dxyError") ;
  add_branch("muon_trackIso03") ;
  add_branch("muon_trackIso05") ;
  add_branch("muon_trackIso03_ptInVeto") ;
  add_branch("muon_trackIso05_ptInVeto") ;
  add_branch("muon_emIso03") ;
  add_branch("muon_emIso05") ;
  add_branch("muon_emIso03_ptInVeto") ;
  add_branch("muon_emIso05_ptInVeto") ;
  add_branch("muon_hadIso03") ;
  add_branch("muon_hadIso05") ;
  add_branch("muon_hadIso03_ptInVeto") ;
  add_branch("muon_hadIso05_ptInVeto") ;
  add_branch("muon_innerPosx") ;
  add_branch("muon_innerPosy") ;
  add_branch("muon_innerPosz") ;
  set_branch_type(kVectorBool) ;
  add_branch("muMatch_hltL1sMu16Eta2p1") ;
  add_branch("muMatch_hltL1sL1Mu3p5EG12") ;
  add_branch("muMatch_hltL1Mu3p5EG12L3Filtered22") ;
  add_branch("muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q") ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                               Configure everything                                 //
  ////////////////////////////////////////////////////////////////////////////////////////
  configure_branches() ;
}

void HEEPTree::configure_branches(){
  for(unsigned int i=0 ; i<all_vars_.size() ; i++){
    all_vars_.at(i)->config(mytree) ;
  }
  return ;
}

// ------------ method called to for each event  ------------
void HEEPTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  beginEvent() ;
  // rho variable
  rho = 0 ;
  edm::Handle<double> rho_ ;
  bool isrho = iEvent.getByLabel(edm::InputTag("kt6PFJets:rho"),rho_) ;
  if(isrho) rho = *rho_ ;

  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                  Event information                                 //
  ////////////////////////////////////////////////////////////////////////////////////////
  store("eventNumber"    , (int)(iEvent.id().event()          )) ;
  store("runNumber"      , (int)(iEvent.id().run()            )) ;
  store("luminosityBlock", (int)(iEvent.id().luminosityBlock())) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                Beamspot information                                //
  ////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotLabel_, theBeamSpot);
  
  // Get all the beam spot info
  //float sigmaZ = theBeamSpot->sigmaZ() ;
  //float sigmaZ0Error = theBeamSpot->sigmaZ0Error() ;
  //float sq = sqrt(sigmaZ*sigmaZ+sigmaZ0Error*sigmaZ0Error) ;
  
  //float bsposx = theBeamSpot->position().x();
  //float bsposy = theBeamSpot->position().y();
  //float bsposz = theBeamSpot->position().z();

  math::XYZPoint beamspot(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  math::XYZPoint firstpvertex(0.,0.,0.);
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                Trigger information                                 //
  ////////////////////////////////////////////////////////////////////////////////////////
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT");
  edm::Handle<trigger::TriggerEvent> trigEvent; 
  iEvent.getByLabel(trigEventTag,trigEvent);
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                  Primary vertices                                  //
  ////////////////////////////////////////////////////////////////////////////////////////
  // Retrieve primary vertex collection
  Handle<reco::VertexCollection> primaryVertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVertexColl);
  const reco::VertexCollection* pvcoll = primaryVertexColl.product();

  // We take only the first primary vertex, i.e. the one with the electrons
  if(pvcoll->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcoll->begin();
    firstpvertex.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  // There are two vtx collections (at least) : offlinePrimaryVertices and offlinePrimaryVerticeswithBS
  // Now storing info about the second one
  math::XYZPoint firstpvertexwithBS(0.,0.,0.);
  Handle<reco::VertexCollection> primaryVertexCollwithBS;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primaryVertexCollwithBS);
  const reco::VertexCollection* pvcollwithBS = primaryVertexCollwithBS.product();

  if(pvcollwithBS->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcollwithBS->begin();
    firstpvertexwithBS.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  //Retrieve primary vertex collection
  float pvz = -999 ;
  bool first_pv = true ;
  for(reco::VertexCollection::const_iterator pvIt = pvcoll->begin(); pvIt != pvcoll->end(); ++pvIt){
    store("pv_x"           , pvIt->x()) ;
    store("pv_y"           , pvIt->y()) ;
    store("pv_z"           , pvIt->z()) ;
    store("pv_isValid"     , pvIt->isValid()) ;
    store("pv_ndof"        , (int)pvIt->ndof()) ;
    store("pv_nTracks"     , pvIt->nTracks()) ;
    store("pv_normChi2"    , pvIt->normalizedChi2()) ;
    store("pv_totTrackSize", (int)(pvIt->tracksSize())) ;
    if(first_pv==true){
      // We need the z coordinate of the first vertex later on...
      pvz = pvIt->z() ;
      first_pv = false ;
    }
  }
  store("pv_n", (int)pvcoll->size()) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                   Superclusters                                    //
  ////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<reco::SuperClusterCollection> pHybridSuperClusters;
  edm::Handle<reco::SuperClusterCollection> pIslandSuperClusters;
  
  try{
    iEvent.getByLabel("correctedHybridSuperClusters"               ,"",pHybridSuperClusters);
    iEvent.getByLabel("correctedMulti5x5SuperClustersWithPreshower","",pIslandSuperClusters);
  }
  catch(cms::Exception &ex){}
  
  const reco::SuperClusterCollection *hybridSuperClusters = pHybridSuperClusters.product() ;
  const reco::SuperClusterCollection *islandSuperClusters = pIslandSuperClusters.product() ;
  
  // Merge these two supercluster collections into one (sclusters collection)
  // There's probably a slicker way to do this...
  std::vector<const reco::SuperCluster*>    sclusters ;
  std::vector<reco::SuperClusterRef>     refsclusters ;
  for(reco::SuperClusterCollection::const_iterator hsc = hybridSuperClusters->begin() ; hsc!=hybridSuperClusters->end() ; hsc++ ){ sclusters.push_back(&(*hsc)) ; }
  for(reco::SuperClusterCollection::const_iterator isc = islandSuperClusters->begin() ; isc!=islandSuperClusters->end() ; isc++ ){ sclusters.push_back(&(*isc)) ; }
  for(unsigned int i=0 ; i<hybridSuperClusters->size() ; i++){ refsclusters.push_back(reco::SuperClusterRef(pHybridSuperClusters,i)) ; }
  for(unsigned int i=0 ; i<islandSuperClusters->size() ; i++){ refsclusters.push_back(reco::SuperClusterRef(pIslandSuperClusters,i)) ; }
  
  // Sort all the refSC and SC by transverse energy
  std::sort(refsclusters.begin(),refsclusters.end(),refScEtGreater);
  std::sort(   sclusters.begin(),   sclusters.end(),   scEtGreater);
  
  unsigned int scsize = 0 ;
  for(unsigned int i_sc=0 ; i_sc<sclusters.size() ; i_sc++){
    reco::SuperCluster* sc = (reco::SuperCluster*)sclusters.at(i_sc) ;
    float sc_energy = sc->rawEnergy()+sc->preshowerEnergy() ;
    float sc_et     = sc_energy/cosh(sc->eta()) ;
    if( sc_et <= ScPtMin_) continue ;
    scsize++ ;
      
    store("sc_eta"      , sc->eta()) ;
    store("sc_etacorr"  , etacorr( sc->eta(), pvz, sc->position().z() )) ;
    store("sc_theta"    , 2.*atan(exp(-1.*sc->eta()))) ;
    store("sc_thetacorr", 2.*atan(exp(-1.*etacorr( sc->eta(), pvz, sc->position().z() ) ))) ;
    store("sc_phi"      , sc->phi()) ;
    store("sc_energy"   , sc_energy) ;
    store("sc_et"       , sc_et    );
    store("sc_px"       , sc_et*cos(sc->phi())) ;
    store("sc_py"       , sc_et*sin(sc->phi()));
    store("sc_pz"       , sc_energy*tanh(sc->eta())) ;
    store("sc_x"        , sc->position().x()) ;
    store("sc_y"        , sc->position().y()) ;
    store("sc_z"        , sc->position().z()) ;
    
    // Trigger matching 
    // HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50 2nd leg 
    // The second leg is a sc, not a gsf (T&P trigger)
    std::string filterName = "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter";
    // It is important to specify the right HLT process for the filter, not doing this is a common bug
    bool scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter = false ;
    trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())) ;
    if(filterIndex<trigEvent->sizeFilters()){
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex) ;
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects()) ;
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){
        const trigger::TriggerObject& obj = trigObjColl[*keyIt] ;
        if(deltaR(sc->eta(),sc->phi(),obj.eta(), obj.phi())<0.3){
          scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter = true ;
          break ;
        }
      }
    }
    store("scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter) ;
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                     GSF tracks                                     //
  ////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<GsfTrackCollection> gsfTracksH ;
  iEvent.getByLabel("electronGsfTracks",gsfTracksH) ;
  const GsfTrackCollection *gsftracks = gsfTracksH.product();
  
  edm::Handle<reco::GsfElectronCollection> pGsfElectrons;
  iEvent.getByLabel("gedGsfElectrons","",pGsfElectrons);
  reco::GsfElectronCollection gsfelectrons(pGsfElectrons->begin(),pGsfElectrons->end());
  
  int gsftrack_n = 0 ;
  for(GsfTrackCollection::const_iterator gsftrackiter = gsftracks->begin() ; gsftrackiter!=gsftracks->end() ; ++gsftrackiter){
    if(gsftrackiter->pt()<GsfTrackPtMin_) continue ;
    gsftrack_n++ ;
    store("gsftrack_eta", gsftrackiter->eta()) ;
    store("gsftrack_phi", gsftrackiter->phi()) ;  
    store("gsftrack_p"  , gsftrackiter->p()  ) ;
    store("gsftrack_pt" , gsftrackiter->pt() ) ;
    store("gsftrack_px" , gsftrackiter->px() ) ;
    store("gsftrack_py" , gsftrackiter->py() ) ;
    store("gsftrack_pz" , gsftrackiter->pz() ) ;
  }
  store("gsftrack_n", gsftrack_n) ;
  
  EcalClusterLazyTools lazytool(iEvent,iSetup,InputTag("reducedEcalRecHitsEB"),InputTag("reducedEcalRecHitsEE"),InputTag("reducedEcalRecHitsES"));

  int gsf_size          = 0 ;
  int gsf0_crystal_size = 0 ; 
  int gsf1_crystal_size = 0 ; 
  //int pfele_size = 0 ;

  edm::ESHandle<CaloGeometry> pGeometry ;
  iSetup.get<CaloGeometryRecord>().get(pGeometry) ;
  CaloGeometry* geometry = (CaloGeometry*) pGeometry.product() ;
  
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  
  Handle<EcalRecHitCollection> EBhits;
  Handle<EcalRecHitCollection> EEhits;
  iEvent.getByLabel("reducedEcalRecHitsEB",EBhits);
  iEvent.getByLabel("reducedEcalRecHitsEE",EEhits);
  
  for(reco::GsfElectronCollection::const_iterator gsfiterforptcut = gsfelectrons.begin() ; gsfiterforptcut != gsfelectrons.end() ; ++gsfiterforptcut){
    if(gsfiterforptcut->caloEnergy()*sin(gsfiterforptcut->p4().theta())<GsfPtMin_) continue;
    
    if(fabs((*gsfiterforptcut).superCluster()->eta())<1.479){//First : Barrel
      for(reco::CaloCluster_iterator bcIt = (*gsfiterforptcut).superCluster()->clustersBegin() ; bcIt!=(*gsfiterforptcut).superCluster()->clustersEnd() ; ++bcIt){
        for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin() ; rhIt != (*bcIt)->hitsAndFractions().end() ; ++rhIt){
          // Loop over rec hits in basic cluster
          for(EcalRecHitCollection::const_iterator it = EBhits->begin() ; it!=EBhits->end() ; ++it) { //loop over all rec hits to find the right ones
            if(rhIt->first==(*it).id()){ //found the matching rechit
              if(gsf_size==0) gsf0_crystal_size++;
              if(gsf_size==1) gsf1_crystal_size++; 
            }
          }
        }
      }
    }
    else{//Now looking at endcaps rechits
      for(reco::CaloCluster_iterator bcIt = (*gsfiterforptcut).superCluster()->clustersBegin() ; bcIt!=(*gsfiterforptcut).superCluster()->clustersEnd() ; ++bcIt){
        for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin() ; rhIt!=(*bcIt)->hitsAndFractions().end() ; ++rhIt){
        // Loop over rec hits in basic cluster
          for(EcalRecHitCollection::const_iterator it = EEhits->begin() ; it!=EEhits->end() ; ++it) { //loop over all rec hits to find the right ones
            if(rhIt->first==(*it).id()){ //found the matching rechit
              if(gsf_size==0) gsf0_crystal_size++ ;
              if(gsf_size==1) gsf1_crystal_size++ ;
            }
          }
        }
      }
    }
    gsf_size++;
  }
  
  // Conversions
  int conv_size = 0;
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);
  for(reco::ConversionCollection::const_iterator conv = hConversions->begin(); conv!= hConversions->end(); ++conv){
    reco::Vertex vtx = conv->conversionVertex();
    if(vtx.isValid()){
      for(reco::GsfElectronCollection::const_iterator gsfiterforconv = gsfelectrons.begin(); gsfiterforconv!=gsfelectrons.end(); ++gsfiterforconv){
        if(ConversionTools::matchesConversion(*gsfiterforconv, *conv)){
          conv_size++ ;
          break;
        }
      }
    }
  }
  for(reco::ConversionCollection::const_iterator conv=hConversions->begin() ; conv!=hConversions->end() ; ++conv){
    reco::Vertex vtx = conv->conversionVertex();
    if(vtx.isValid()){
      int iel = -1 ;
      for(reco::GsfElectronCollection::const_iterator gsfiterforconv=gsfelectrons.begin() ; gsfiterforconv!=gsfelectrons.end() ; ++gsfiterforconv){
        iel++;
        if(ConversionTools::matchesConversion(*gsfiterforconv, *conv)){
          store("conv_eleind", iel);
          store("conv_vtxProb", TMath::Prob( vtx.chi2(), vtx.ndof()) ) ;
          math::XYZVector mom(conv->refittedPairMomentum());
          double dbsx = vtx.x() - theBeamSpot->position().x();   
          double dbsy = vtx.y() - theBeamSpot->position().y();
          store("conv_lxy",(mom.x()*dbsx + mom.y()*dbsy)/mom.rho());
          store("conv_nHitsMax", 0);
          int nHitsMax = 0 ;
          for (std::vector<uint8_t>::const_iterator it = conv->nHitsBeforeVtx().begin(); it!=conv->nHitsBeforeVtx().end(); ++it) {
            if ((*it)>nHitsMax) nHitsMax = (int)(*it) ;
          }
          store("conv_nHitsMax", nHitsMax) ;
          break;
        }
      }
    }
  }
  //End of conversion info
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                   GSF electrons                                    //
  ////////////////////////////////////////////////////////////////////////////////////////
  const CaloSubdetectorGeometry* geometryES = geometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower) ;
  CaloSubdetectorTopology* topology_p = 0 ;
  if(geometryES) topology_p = new EcalPreshowerTopology(geometry) ;
   
  int e = 0 ;
  int nHeepEle = 0 ;
  store("gsf_n", (int)gsfelectrons.size()) ;
  for(reco::GsfElectronCollection::const_iterator gsfiter = gsfelectrons.begin() ; gsfiter!=gsfelectrons.end() ; ++gsfiter){
    if(gsfiter->caloEnergy()*sin(gsfiter->p4().theta())<GsfPtMin_) continue;
    int index = -3 ;
    // Try to get the index for the sc assoc to this gsf
    reco::SuperClusterRef gsfrefsc = gsfiter->superCluster();
    for(unsigned int k=0;k<refsclusters.size();k++){
      if(gsfrefsc == refsclusters[k]){
        index = k ;
      }
    }
    store("sc_indexforgsf", index) ;
    
    // Required for preshower variables
    reco::SuperClusterRef cl_ref = gsfiter->superCluster();
    
    //Fill the gsf related variables
    float gsf_deltaeta        = gsfiter->deltaEtaSuperClusterTrackAtVtx() ;
    float gsf_deltaphi        = gsfiter->deltaPhiSuperClusterTrackAtVtx() ;
    float gsf_hovere          = gsfiter->hadronicOverEm() ;
    float gsf_sigmaIetaIeta   = gsfiter->sigmaIetaIeta() ;
    float gsf_e2x5overe5x5    = gsfiter->scE2x5Max()/gsfiter->scE5x5() ;
    float gsf_e1x5overe5x5    = gsfiter->scE1x5()/gsfiter->scE5x5() ;
    float gsf_ecaliso         = gsfiter->dr03EcalRecHitSumEt() ;
    float gsf_hcaliso1        = gsfiter->dr03HcalDepth1TowerSumEt() ;
    float gsf_hcaliso2        = gsfiter->dr03HcalDepth2TowerSumEt() ;
    float gsf_trackiso        = gsfiter->dr03TkSumPt() ;
    bool  gsf_isecaldriven    = gsfiter->ecalDrivenSeed() ;
    bool  gsf_istrackerdriven = gsfiter->trackerDrivenSeed() ;
    float gsf_gsfet           = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
    int gsf_nLostInnerHits    = gsfiter->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ;
    int gsf_nLostOuterHits    = gsfiter->gsfTrack()->trackerExpectedHitsOuter().numberOfLostHits() ;
    int gsf_convFlags         = gsfiter->convFlags() ;
    float gsf_dxy_firstPVtx   = gsfiter->gsfTrack()->dxy(firstpvertex) ;
    store("gsf_e"                  , gsfiter->energy()) ;
    store("gsf_p"                  , gsfiter->p()) ;
    store("gsf_pt"                 , gsfiter->pt()) ;
    store("gsf_class"              , gsfiter->classification()) ;
    store("gsf_e2x5overe5x5"       , gsf_e2x5overe5x5) ;
    store("gsf_e1x5overe5x5"       , gsf_e1x5overe5x5) ;
    store("gsf_eta"                , gsfiter->eta()) ;
    store("gsf_phi"                , gsfiter->phi()) ;
    store("gsf_theta"              , gsfiter->theta());
    store("gsf_px"                 , gsfiter->px()) ;
    store("gsf_py"                 , gsfiter->py()) ;
    store("gsf_pz"                 , gsfiter->pz()) ;
    store("gsf_deltaeta"           , gsf_deltaeta) ;
    store("gsf_deltaphi"           , gsf_deltaphi) ;
    store("gsf_hovere"             , gsf_hovere) ;
    store("gsf_hdepth1overe"       , gsfiter->hcalDepth1OverEcal()) ;
    store("gsf_hdepth2overe"       , gsfiter->hcalDepth2OverEcal()) ;
    store("gsf_trackiso"           , gsf_trackiso) ;
    store("gsf_ecaliso"            , gsf_ecaliso) ;
    store("gsf_hcaliso1"           , gsf_hcaliso1) ;
    store("gsf_hcaliso2"           , gsf_hcaliso2) ;
    store("gsf_charge"             , gsfiter->charge()) ;
    store("gsf_sigmaetaeta"        , gsf_sigmaIetaIeta) ;
    store("gsf_sigmaIetaIeta"      , gsfiter->sigmaIetaIeta()) ;
    store("gsf_isecaldriven"       , gsf_isecaldriven   ) ;
    store("gsf_istrackerdriven"    , gsf_istrackerdriven) ;
    store("gsf_eseffsixix"         , lazytool.eseffsixix(*cl_ref)) ;
    store("gsf_eseffsiyiy"         , lazytool.eseffsiyiy(*cl_ref)) ;
    store("gsf_eseffsirir"         , lazytool.eseffsirir(*cl_ref)) ;
    store("gsf_preshowerEnergy"    , gsfiter->superCluster()->preshowerEnergy()) ;
    store("gsf_gsfet"              , gsf_gsfet) ;
    store("gsf_isEB"               , gsfiter->isEB());
    store("gsf_isEE"               , gsfiter->isEE());
    store("gsf_deltaEtaATcalo"     , gsfiter->deltaEtaSeedClusterTrackAtCalo());
    store("gsf_deltaPhiATcalo"     , gsfiter->deltaPhiSeedClusterTrackAtCalo());
    store("gsf_ecalEnergy"         , gsfiter->ecalEnergy());
    store("gsf_eOVERp"             , gsfiter->eSuperClusterOverP());
    store("gsf_dxy"                , gsfiter->gsfTrack()->dxy());
    store("gsf_dxy_beamSpot"       , gsfiter->gsfTrack()->dxy(beamspot));
    store("gsf_dxy_firstPVtx"      , gsf_dxy_firstPVtx);
    store("gsf_dxy_firstPVtxwithBS", gsfiter->gsfTrack()->dxy(firstpvertexwithBS));
    store("gsf_dxyError"           , gsfiter->gsfTrack()->dxyError());
    store("gsf_dz"                 , gsfiter->gsfTrack()->dz());
    store("gsf_dz_beamSpot"        , gsfiter->gsfTrack()->dz(beamspot));
    store("gsf_dz_firstPVtx"       , gsfiter->gsfTrack()->dz(firstpvertex));
    store("gsf_dz_firstPVtxwithBS" , gsfiter->gsfTrack()->dz(firstpvertexwithBS));
    store("gsf_dzError"            , gsfiter->gsfTrack()->dzError()); 
    store("gsf_vz"                 , gsfiter->gsfTrack()->vz());
    store("gsf_nHits"              , gsfiter->gsfTrack()->numberOfValidHits());
    store("gsf_nLostInnerHits"     , gsf_nLostInnerHits);
    store("gsf_nLostOuterHits"     , gsf_nLostOuterHits);
    store("gsf_convFlags"          , gsf_convFlags);
    store("gsf_convDist"           , gsfiter->convDist());
    store("gsf_convDcot"           , gsfiter->convDcot());
    store("gsf_convRadius"         , gsfiter->convRadius());
    store("gsf_fBrem"              , gsfiter->fbrem());
    store("gsf_e1x5"               , gsfiter->e1x5()) ;
    store("gsf_e2x5"               , gsfiter->e2x5Max()) ;
    store("gsf_e5x5"               , gsfiter->e5x5()) ;
    
    // Get the preshower hits
    double x = gsfiter->superCluster()->x() ;
    double y = gsfiter->superCluster()->y() ;
    double z = gsfiter->superCluster()->z() ;
    std::vector<float> phoESHitsIXIX = lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_p, 0, 1);
    std::vector<float> phoESHitsIYIY = lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_p, 0, 2);
    store("gsf_eshitsixix", phoESHitsIXIX) ;
    store("gsf_eshitsiyiy", phoESHitsIYIY) ;
    
    float gsfsc_eta = gsfiter->superCluster()->eta() ;
    store("gsfsc_e"  , gsfiter->superCluster()->energy()) ;
    store("gsfsc_pt" , (gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy())/cosh(gsfiter->superCluster()->eta())) ;
    store("gsfsc_eta", gsfsc_eta) ;
    store("gsfsc_phi", gsfiter->superCluster()->phi()) ;
    store("gsfsc_px" , gsfiter->pt()*cos(gsfiter->phi())) ;
    store("gsfsc_py" , gsfiter->pt()*sin(gsfiter->phi())) ;
    store("gsfsc_pz" , (gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy())*tanh(gsfiter->superCluster()->eta())) ;
    
    //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DataFormats/TrackReco/interface/HitPattern.h?revision=1.32&view=markup
    reco::HitPattern kfHitPattern = gsfiter->gsfTrack()->hitPattern();
    int nbtrackhits = kfHitPattern.numberOfHits();
    
    std::vector<int> gsf_hitsinfo ;
    for(int hititer=0; hititer<25;hititer++){
      int myhitbin = (hititer<nbtrackhits) ? kfHitPattern.getHitPattern(hititer) : 0 ;
      gsf_hitsinfo.push_back(myhitbin) ;
    }
    store("gsf_hitsinfo", gsf_hitsinfo) ;
    
    // Try to add info about rechit in the SC 
    // Strongly inspired from : http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/DaveC/src/printPhoton.cc
    
    //Crystal variables
    std::vector<float> gsf_crystal_energy  ;
    std::vector<int  > gsf_crystal_ietaorix;
    std::vector<int  > gsf_crystal_iphioriy;
    std::vector<float> gsf_crystal_eta     ;
    
    if(fabs((*gsfiter).superCluster()->eta())<1.479){//First : Barrel
      int iebhit = -1, nclust = 0;
      double amplitot = 0.0;
      double clustot = 0.0;
      
      for(reco::CaloCluster_iterator bcIt = (*gsfiter).superCluster()->clustersBegin() ; bcIt!=(*gsfiter).superCluster()->clustersEnd() ; ++bcIt) {
        // Loop over basic clusters in SC
        // bcIt seems to be a pointer to a pointer
   
        double clusterEnergy = (*bcIt)->energy();
        clustot += clusterEnergy;
        nclust += 1;
        for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin() ; rhIt!=(*bcIt)->hitsAndFractions().end() ; ++rhIt) {
          // Loop over rec hits in basic cluster
          for(EcalRecHitCollection::const_iterator it=EBhits->begin() ; it!=EBhits->end() ; ++it){
            // Loop over all rec hits to find the right ones
            if(rhIt->first==(*it).id()){ // Found the matching rechit
              iebhit +=1; 
              EcalRecHit hit = (*it);
              EBDetId det    = hit.id(); 
              float ampli    = hit.energy();
              amplitot += ampli;
        
              GlobalPoint poseb = geometry->getPosition(hit.detid());
              float eta_eb = poseb.eta();
              int   ieta   = det.ieta() ;      
              int   iphi   = det.iphi() ;
        
              gsf_crystal_energy  .push_back(ampli ) ;
              gsf_crystal_ietaorix.push_back(ieta  ) ;
              gsf_crystal_iphioriy.push_back(iphi  ) ;
              gsf_crystal_eta     .push_back(eta_eb) ;
            }
          }
        }
      }
    }
    //Now looking at endcaps rechits
    else{
      int ieehit = -1, nclustee = 0;
      double amplitotee = 0.0, clustotee = 0.0;
      for(reco::CaloCluster_iterator bcIt = (*gsfiter).superCluster()->clustersBegin() ; bcIt!=(*gsfiter).superCluster()->clustersEnd() ; ++bcIt){
        // Loop over basic clusters in SC
        nclustee +=1;
        double clusterEnergyee = (*bcIt)->energy();
        clustotee += clusterEnergyee;
        for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin() ; rhIt!=(*bcIt)->hitsAndFractions().end() ; ++rhIt){
          // Loop over rec hits in basic cluster
          for(EcalRecHitCollection::const_iterator it = EEhits->begin() ; it!=EEhits->end(); ++it){
            // Loop over all rec hits to find the right ones
            if(rhIt->first==(*it).id() ){ //found the matching rechit
              ieehit += 1;
              EcalRecHit hit = (*it);
              EEDetId det = hit.id(); 
              float ampli = hit.energy();
              amplitotee += ampli;
        
              GlobalPoint posee = geometry->getPosition(hit.detid());
              float eta_ee = posee.eta();
              int   ix     = det.ix();
              int   iy     = det.iy();

              gsf_crystal_energy  .push_back(ampli ) ;
              gsf_crystal_ietaorix.push_back(ix    ) ;
              gsf_crystal_iphioriy.push_back(iy    ) ;
              gsf_crystal_eta     .push_back(eta_ee) ;
            }
          }
        }
      }
    }
    store("gsf_crystal_energy"  , gsf_crystal_energy  ) ;
    store("gsf_crystal_ietaorix", gsf_crystal_ietaorix) ;
    store("gsf_crystal_iphioriy", gsf_crystal_iphioriy) ;
    store("gsf_crystal_eta"     , gsf_crystal_eta     ) ;
    
    const reco::CaloClusterPtr seed = gsfiter->superCluster()->seed();
    store("gsf_e1x3", lazytool.e1x3(*seed)) ;
    
    //////////////////////////////////////////////////////////////////////////////////////
    //                                    HEEP cutflow                                  //
    //////////////////////////////////////////////////////////////////////////////////////
    // HEEP selection v4.0  - 07/05/2012 Laurent 
    bool gsfetbarrel         = gsf_gsfet > 35.0 ;
    bool gsfetendcap         = gsf_gsfet > 35.0 ;
    bool barrelsc            = fabs(gsfsc_eta) < 1.442 ;
    bool endcapsc            = (fabs(gsfsc_eta) > 1.56) && (fabs(gsfsc_eta) < 2.5) ;
    bool deltaetabarrel      = fabs(gsf_deltaeta) < 0.005 ;
    bool deltaetaendcap      = fabs(gsf_deltaeta) < 0.007 ;
    bool deltaphibarrel      = fabs(gsf_deltaphi) < 0.06  ;
    bool deltaphiendcap      = fabs(gsf_deltaphi) < 0.06  ;
    bool hoverebarrel        = gsf_hovere < 0.05 ;
    bool hovereendcap        = gsf_hovere < 0.05 ;
    bool sigmaIetaIetabarrel = true ;
    bool sigmaIetaIetaendcap = gsf_sigmaIetaIeta < 0.03 ;
    bool e2x5overe5x5barrel  = (gsf_e2x5overe5x5 > 0.94) || (gsf_e1x5overe5x5 > 0.83);
    bool e2x5overe5x5endcap  = true ;
    bool ecalisobarrel       = (gsf_ecaliso+gsf_hcaliso1) < (2.+0.03*gsf_gsfet + rho*EcalHcal1EffAreaBarrel_);
    bool ecalisoendcap       = true ;
    if(gsf_gsfet<50.0){ ecalisoendcap = (gsf_ecaliso+gsf_hcaliso1) <  2.5+ rho*EcalHcal1EffAreaEndcaps_ ; }
    else              { ecalisoendcap = (gsf_ecaliso+gsf_hcaliso1) < (2.5+0.03*(gsf_gsfet-50.0)+ rho*EcalHcal1EffAreaEndcaps_ ) ; }
    bool hcaliso2barrel      = true ;
    bool hcaliso2endcap      = true ;
    bool trackisobarrel      = gsf_trackiso<5.0 ;
    bool trackisoendcap      = gsf_trackiso<5.0 ;
    bool noMissingHits       = gsf_nLostInnerHits<=1 ;
    bool noConversion        = gsf_convFlags != 3 ; 
    bool dxyfirstpvbarrel    = fabs(gsf_dxy_firstPVtx) <0.02 ;
    bool dxyfirstpvendcaps   = fabs(gsf_dxy_firstPVtx) <0.05 ;

    //Boolean HEEP cuts
    bool gsfpass_ET              = (gsfetbarrel && barrelsc) || (gsfetendcap && endcapsc) ;
    bool gsfpass_PT              = true ;
    bool gsfpass_DETETA          = true ;
    bool gsfpass_CRACK           = true ;
    bool gsfpass_DETAIN          = (deltaetabarrel && barrelsc) || (deltaetaendcap && endcapsc) ;
    bool gsfpass_DPHIIN          = (deltaphibarrel && barrelsc) || (deltaphiendcap && endcapsc) ;
    bool gsfpass_HADEM           = (hoverebarrel   && barrelsc) || (hovereendcap   && endcapsc) ;
    bool gsfpass_SIGMAIETAIETA   = (sigmaIetaIetabarrel && barrelsc) || (sigmaIetaIetaendcap && endcapsc) ;
    bool gsfpass_E2X5OVER5X5     = (e2x5overe5x5barrel  && barrelsc) || (e2x5overe5x5endcap  && endcapsc) ;
    bool gsfpass_ISOLEMHADDEPTH1 = (ecalisobarrel  && barrelsc) || (ecalisoendcap  && endcapsc) ;
    bool gsfpass_ISOLHADDEPTH2   = (hcaliso2barrel && barrelsc) || (hcaliso2endcap && endcapsc) ;
    bool gsfpass_ISOLPTTRKS      = (trackisobarrel && barrelsc) || (trackisoendcap && endcapsc) ;
    bool gsfpass_ECALDRIVEN      = gsf_isecaldriven ;
    bool gsfpass_INVALID         = true ;
    bool gsfpass_NOMISSINGHITS   = noMissingHits ;
    bool gsfpass_NOCONVERSION    = noConversion ;
    bool gsfpass_DXYFIRSTPV      = (dxyfirstpvbarrel && barrelsc) || (dxyfirstpvendcaps && endcapsc) ;
    bool gsfpass_ID              = (gsfpass_DETAIN && gsfpass_DPHIIN && gsfpass_HADEM && gsfpass_SIGMAIETAIETA && gsfpass_E2X5OVER5X5) ;
    bool gsfpass_ISO             = (gsfpass_ISOLEMHADDEPTH1 && gsfpass_ISOLHADDEPTH2 && gsfpass_ISOLPTTRKS) ;
    bool gsfpass_HEEP = gsfpass_ET && gsfpass_PT && gsfpass_DETETA && gsfpass_CRACK&& gsfpass_DETAIN && gsfpass_DPHIIN && gsfpass_HADEM && gsfpass_SIGMAIETAIETA && gsfpass_E2X5OVER5X5     && gsfpass_ISOLEMHADDEPTH1 && gsfpass_ISOLHADDEPTH2 && gsfpass_ISOLPTTRKS && gsfpass_ECALDRIVEN && gsfpass_INVALID && gsfpass_NOMISSINGHITS && gsfpass_NOCONVERSION && gsfpass_DXYFIRSTPV && gsfpass_ID && gsfpass_ISO ;
    
    store("gsfpass_ET"             , gsfpass_ET              ) ;
    store("gsfpass_PT"             , gsfpass_PT              ) ;
    store("gsfpass_DETETA"         , gsfpass_DETETA          ) ;
    store("gsfpass_CRACK"          , gsfpass_CRACK           ) ;
    store("gsfpass_DETAIN"         , gsfpass_DETAIN          ) ;
    store("gsfpass_DPHIIN"         , gsfpass_DPHIIN          ) ;
    store("gsfpass_HADEM"          , gsfpass_HADEM           ) ;
    store("gsfpass_SIGMAIETAIETA"  , gsfpass_SIGMAIETAIETA   ) ;
    store("gsfpass_E2X5OVER5X5"    , gsfpass_E2X5OVER5X5     ) ;
    store("gsfpass_ISOLEMHADDEPTH1", gsfpass_ISOLEMHADDEPTH1 ) ;
    store("gsfpass_ISOLHADDEPTH2"  , gsfpass_ISOLHADDEPTH2   ) ;
    store("gsfpass_ISOLPTTRKS"     , gsfpass_ISOLPTTRKS      ) ;
    store("gsfpass_ECALDRIVEN"     , gsfpass_ECALDRIVEN      ) ;
    store("gsfpass_INVALID"        , gsfpass_INVALID         ) ;
    store("gsfpass_NOMISSINGHITS"  , gsfpass_NOMISSINGHITS   ) ;
    store("gsfpass_NOCONVERSION"   , gsfpass_NOCONVERSION    ) ;
    store("gsfpass_DXYFIRSTPV"     , gsfpass_DXYFIRSTPV      ) ;
    store("gsfpass_ID"             , gsfpass_ID              ) ;
    store("gsfpass_ISO"            , gsfpass_ISO             ) ;
    store("gsfpass_HEEP"           , gsfpass_HEEP            ) ;
    if(gsfpass_HEEP) ++nHeepEle ;

    //charge info
    store("sc_pixcharge", gsfiter->scPixCharge()) ;
    if(gsfiter->closestCtfTrackRef().isNonnull()) store("ctf_charge", gsfiter->closestCtfTrackRef()->charge());
    store("gsf_charge", gsfiter->gsfTrack()->charge());
    store("gsf_ctfscpixconsistent", gsfiter->isGsfCtfScPixChargeConsistent()) ;
    store("gsf_scpixconsistent"   , gsfiter->isGsfScPixChargeConsistent()   ) ;
    store("gsf_ctfconsistent"     , gsfiter->isGsfCtfChargeConsistent()     ) ;

    //////////////////////////////////////////////////////////////////////////////////////
    //                                 Trigger matching                                 //
    //////////////////////////////////////////////////////////////////////////////////////
    const double barrelEnd       = 1.4791;
    const double regionEtaSizeEB = 0.522 ;
    const double regionEtaSizeEE = 1.0   ;
    const double regionPhiSize   = 1.044 ;
    
    std::string branchPrefix = "gsfMatch_" ;
    std::vector<std::string> filterNames ;
    filterNames.push_back("hltL1sL1SingleEG12") ;
    filterNames.push_back("hltL1sL1Mu3p5EG12" ) ;
    filterNames.push_back("hltL1sL1SingleEG22") ;
    
    // L1 hltL1sL1SingleEG12
    // Careful that L1 triggers only have discrete eta phi. Need to be extremely loose. 
    // See here: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/SHarper/SHNtupliser/src/SHTrigInfo.cc?revision=1.5&view=markup&pathrev=HEAD
    // It is important to specify the right HLT process for the filter, not doing this is a common bug
    for(unsigned int iFilter=0 ; iFilter<filterNames.size() ; iFilter++){
      bool gsfMatch = false ;
      int filterIndex = trigEvent->filterIndex(edm::InputTag(filterNames.at(iFilter),"",trigEventTag.process())); 
      if(filterIndex<trigEvent->sizeFilters()){ 
        const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
        const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
        // Now loop of the trigger objects passing filter
        for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
          const trigger::TriggerObject& obj = trigObjColl[*keyIt];
          // Do what you want with the trigger objects, you have
          // eta,phi,pt,mass,p,px,py,pz,et,energy accessors
          
          float objeta = obj.eta(); 
          float objphi = obj.phi();
          
          double etaBinLow  = 0.0 ;
          double etaBinHigh = 0.0 ;
          
          if(fabs(objeta) < barrelEnd){
            etaBinLow  = objeta - regionEtaSizeEB/2.;
            etaBinHigh = etaBinLow + regionEtaSizeEB;
          }
          else{
            etaBinLow  = objeta - regionEtaSizeEE/2.;
            etaBinHigh = etaBinLow + regionEtaSizeEE;
          }
          
          float deltaPhi = reco::deltaPhi(gsfiter->phi(),objphi);
          
          if(gsfiter->eta() < etaBinHigh && gsfiter->eta() > etaBinLow &&   deltaPhi <regionPhiSize/2. )  {
            gsfMatch = true ;
            break ;
          }
        }
      }
      std::string branchName = branchPrefix + filterNames.at(iFilter) ;
      store(branchName, gsfMatch) ;
    }//end filter size check
    
    filterNames.clear() ;
    filterNames.push_back("hltEle33CaloIdLPixelMatchFilter") ;
    filterNames.push_back("hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter") ;
    filterNames.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter") ;
    filterNames.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter") ;
    filterNames.push_back("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter") ;
    filterNames.push_back("hltEle27WP80TrackIsoFilter") ;
    filterNames.push_back("hltMu22Photon22CaloIdLHEFilter") ;
    for(unsigned iFilter=0 ; iFilter<filterNames.size() ; iFilter++){
      bool gsfMatch = false ;
      // It is important to specify the right HLT process for the filter, not doing this is a common bug
      trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterNames.at(iFilter),"",trigEventTag.process())); 
      if(filterIndex<trigEvent->sizeFilters()){ 
        const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
        const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
        // Now loop over the trigger objects passing filter
        for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); ++keyIt) { 
          const trigger::TriggerObject& obj = trigObjColl[*keyIt];
          if(deltaR(gsfiter->eta(),gsfiter->phi(),obj.eta(), obj.phi())<0.5){
            gsfMatch = true ;
          }
        }
      }//end filter size check
      std::string branchName = branchPrefix + filterNames.at(iFilter) ;
      store(branchName, gsfMatch) ;
    }

    //increment index for gsf
    e++;
  }
  //Conversion info : https://twiki.cern.ch/twiki/bin/viewauth/CMS/ConversionTools
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                  Muon collection                                   //
  ////////////////////////////////////////////////////////////////////////////////////////
  edm::Handle<reco::MuonCollection> muonCollection;
  iEvent.getByLabel("muons",muonCollection);
  const reco::MuonCollection* muons = muonCollection.product();
  store("muon_n", (int)muons->size()) ;
  for(reco::MuonCollection::const_iterator muIt = muons->begin(); muIt != muons->end(); ++muIt){
    if (muIt->isGlobalMuon()) {
      // get TeV optimized track
      reco::Muon::MuonTrackTypePair tevOptMuTrk = muon::tevOptimized(*muIt, 200, 17., 40., 0.25) ;
      if (tevOptMuTrk.first->pt() < muPtMin_) continue;
      store("muon_pt"                 , tevOptMuTrk.first->pt()                                      ) ;
      store("muon_ptError"            , tevOptMuTrk.first->ptError()                                 ) ;
      store("muon_gTrk_pt"            , muIt->globalTrack()->pt()                                    ) ;
      store("muon_gTrk_ptError"       , muIt->globalTrack()->ptError()                               ) ;
      store("muon_eta"                , tevOptMuTrk.first->eta()                                     ) ;
      store("muon_etaError"           , tevOptMuTrk.first->etaError()                                ) ;
      store("muon_phi"                , tevOptMuTrk.first->phi()                                     ) ;
      store("muon_phiError"           , tevOptMuTrk.first->phiError()                                ) ;
      store("muon_theta"              , tevOptMuTrk.first->theta()                                   ) ;
      store("muon_thetaError"         , tevOptMuTrk.first->thetaError()                              ) ;
      store("muon_outerPt"            , muIt->globalTrack()->outerPt()                               ) ;
      store("muon_outerEta"           , muIt->globalTrack()->outerEta()                              ) ;
      store("muon_outerPhi"           , muIt->globalTrack()->outerPhi()                              ) ;
      store("muon_outerTheta"         , muIt->globalTrack()->outerTheta()                            ) ;
      store("muon_px"                 , tevOptMuTrk.first->px()                                      ) ;
      store("muon_py"                 , tevOptMuTrk.first->py()                                      ) ;
      store("muon_pz"                 , tevOptMuTrk.first->pz()                                      ) ;
      store("muon_charge"             , tevOptMuTrk.first->charge()                                  ) ;
      store("muon_nhitspixel"         , muIt->innerTrack()->hitPattern().numberOfValidPixelHits()    ) ;
      store("muon_nhitstrack"         , muIt->globalTrack()->hitPattern().numberOfValidTrackerHits() ) ;                                  
      store("muon_nhitsmuons"         , muIt->globalTrack()->hitPattern().numberOfValidMuonHits()    ) ;                            
      store("muon_nhitstotal"         , muIt->globalTrack()->numberOfValidHits()                     ) ;        
      store("muon_nlayerswithhits"    , muIt->track()->hitPattern().trackerLayersWithMeasurement()   ) ;        
      store("muon_nlosthits"          , muIt->globalTrack()->numberOfLostHits()                      ) ;        
      store("muon_nSegmentMatch"      , muIt->numberOfMatchedStations()                              ) ;        
      store("muon_isTrackerMuon"      , muIt->isTrackerMuon()                                        ) ;        
      store("muon_isPFMuon"           , muIt->isPFMuon()                                             ) ;        
      store("muon_isPFIsolationValid" , muIt->isPFIsolationValid()                                   ) ;        
      store("muon_chi2"               , muIt->globalTrack()->chi2()                                  ) ;
      store("muon_ndof"               , muIt->globalTrack()->ndof()                                  ) ;
      store("muon_normChi2"           , muIt->globalTrack()->normalizedChi2()                        ) ;        
      store("muon_d0"                 , tevOptMuTrk.first->d0()                                      ) ;
      store("muon_d0Error"            , tevOptMuTrk.first->d0Error()                                 ) ;
      store("muon_dz_cmsCenter"       , tevOptMuTrk.first->dz()                                      ) ;
      store("muon_dz_beamSpot"        , tevOptMuTrk.first->dz(beamspot)                              ) ;
      store("muon_dz_firstPVtx"       , tevOptMuTrk.first->dz(firstpvertex)                          ) ;
      store("muon_dz_firstPVtxwithBS" , tevOptMuTrk.first->dz(firstpvertexwithBS)                    ) ;
      store("muon_dzError"            , tevOptMuTrk.first->dzError()                                 ) ;
      store("muon_dxy_cmsCenter"      , tevOptMuTrk.first->dxy()                                     ) ;
      store("muon_dxy_beamSpot"       , tevOptMuTrk.first->dxy(beamspot)                             ) ;
      store("muon_dxy_firstPVtx"      , tevOptMuTrk.first->dxy(firstpvertex)                         ) ;
      store("muon_dxy_firstPVtxwithBS", tevOptMuTrk.first->dxy(firstpvertexwithBS)                   ) ;
      store("muon_dxyError"           , tevOptMuTrk.first->dxyError()                                ) ;
      store("muon_innerPosx"          , muIt->globalTrack()->innerPosition().X()                     ) ;
      store("muon_innerPosy"          , muIt->globalTrack()->innerPosition().Y()                     ) ;
      store("muon_innerPosz"          , muIt->globalTrack()->innerPosition().Z()                     ) ;
      store("muon_trackIso03"         , muIt->isolationR03().sumPt                                   ) ;
      store("muon_trackIso05"         , muIt->isolationR05().sumPt                                   ) ;
      store("muon_trackIso03_ptInVeto", muIt->isolationR03().trackerVetoPt                           ) ;
      store("muon_trackIso05_ptInVeto", muIt->isolationR05().trackerVetoPt                           ) ;
      store("muon_emIso03"            , muIt->isolationR03().emEt                                    ) ;
      store("muon_emIso05"            , muIt->isolationR05().emEt                                    ) ;
      store("muon_emIso03_ptInVeto"   , muIt->isolationR03().emVetoEt                                ) ;
      store("muon_emIso05_ptInVeto"   , muIt->isolationR05().emVetoEt                                ) ;
      store("muon_hadIso03"           , muIt->isolationR03().hadEt                                   ) ;
      store("muon_hadIso05"           , muIt->isolationR05().hadEt                                   ) ;
      store("muon_hadIso03_ptInVeto"  , muIt->isolationR03().hadVetoEt                               ) ;
      store("muon_hadIso05_ptInVeto"  , muIt->isolationR05().hadVetoEt                               ) ;

      std::string branchPrefix = "muMatch_" ;
      std::vector<std::string> filterNames ;
      filterNames.push_back("hltL1sMu16Eta2p1") ;
      filterNames.push_back("hltL1sL1Mu3p5EG12") ;
      filterNames.push_back("hltL1Mu3p5EG12L3Filtered22") ;
      filterNames.push_back("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q") ;
      
      for(unsigned iFilter=0 ; iFilter<filterNames.size() ; iFilter++){
        bool muMatch = false ;
        // It is important to specify the right HLT process for the filter, not doing this is a common bug
        trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterNames.at(iFilter),"",trigEventTag.process())); 
        if(filterIndex<trigEvent->sizeFilters()){ 
          const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
          const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
          // Now loop over the trigger objects passing filter
          for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); ++keyIt) { 
            const trigger::TriggerObject& obj = trigObjColl[*keyIt];
            if(deltaR(tevOptMuTrk.first->eta(), tevOptMuTrk.first->phi(), obj.eta(), obj.phi())<1.){
              muMatch = true ;
            }
          }
        }//end filter size check
        std::string branchName = branchPrefix + filterNames.at(iFilter) ;
        store(branchName, muMatch) ;
      }
    }
  }
  
  mytree->Fill();
  endEvent() ;
}//end of analyze method


void 
HEEPTree::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

void HEEPTree::beginEvent(){
  for(unsigned int i=0 ; i<all_vars_.size() ; i++){
    all_vars_.at(i)->beginEvent() ;
  }
}

void HEEPTree::endEvent(){
  for(unsigned int i=0 ; i<all_vars_.size() ; i++){
    all_vars_.at(i)->endEvent() ;
  }
}


// ------------ method called once each job just after ending the event loop  ------------
void 
HEEPTree::endJob(){
  std::vector<std::string> untouched_branch_names ;
  for(unsigned int i=0 ; i<all_vars_.size() ; i++){
    if(all_vars_.at(i)->is_touched()==false) untouched_branch_names.push_back(all_vars_.at(i)->name()) ;
  }
  if(debug==true){
    if(untouched_branch_names.size()>0){
      std::cout << "The following branches were never touched:" << std::endl ;
      for(unsigned int i=0 ; i<untouched_branch_names.size() ; i++){
        std::cout << "  " << untouched_branch_names.at(i) << std::endl ;
      }
    }
  }
  
  if(myFile){
    myFile->Write() ;
    delete myFile ;
  }
}

bool HEEPTree::store(std::string name, bool value){
  for(unsigned int i=0 ; i<vars_B_.size() ; i++){
    if(vars_B_.at(i)->name()==name){
      vars_B_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_BV_.size() ; i++){
    if(vars_BV_.at(i)->name()==name){
      vars_BV_ .at(i)->push(value) ;
      return true ;
    }
  }
  if(debug) std::cout << "Could not find a (bool) branch named " << name << std::endl ;
  return false ;
}
bool HEEPTree::store(std::string name, double value){
  // Try to fill doubles, then floats
  for(unsigned int i=0 ; i<vars_D_.size() ; i++){
    if(vars_D_.at(i)->name()==name){
      vars_D_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_DV_.size() ; i++){
    if(vars_DV_.at(i)->name()==name){
      vars_DV_ .at(i)->push(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_F_.size() ; i++){
    if(vars_F_.at(i)->name()==name){
      vars_F_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_FV_.size() ; i++){
    if(vars_FV_.at(i)->name()==name){
      vars_FV_ .at(i)->push(value) ;
      return true ;
    }
  }
  if(debug) std::cout << "Could not find a (double) branch named " << name << std::endl ;
  return false ;
}
bool HEEPTree::store(std::string name, float value){
  // Try to fill floats, then doubles
  for(unsigned int i=0 ; i<vars_F_.size() ; i++){
    if(vars_F_.at(i)->name()==name){
      vars_F_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_FV_.size() ; i++){
    if(vars_FV_.at(i)->name()==name){
      vars_FV_ .at(i)->push(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_D_.size() ; i++){
    if(vars_D_.at(i)->name()==name){
      vars_D_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_DV_.size() ; i++){
    if(vars_DV_.at(i)->name()==name){
      vars_DV_ .at(i)->push(value) ;
      return true ;
    }
  }
  if(debug) std::cout << "Could not find a (float) branch named " << name << std::endl ;
  return false ;
}
bool HEEPTree::store(std::string name, unsigned value){ return store(name, (int) value) ; }
bool HEEPTree::store(std::string name, int value){
  for(unsigned int i=0 ; i<vars_I_.size() ; i++){
    if(vars_I_.at(i)->name()==name){
      vars_I_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_IV_.size() ; i++){
    if(vars_IV_.at(i)->name()==name){
      vars_IV_ .at(i)->push(value) ;
      return true ;
    }
  }
  if(debug) std::cout << "Could not find a (int) branch named " << name << std::endl ;
  return false ;
}

bool HEEPTree::store(std::string name, std::vector<bool> values){
  for(unsigned int i=0 ; i<vars_BVV_.size() ; i++){
    if(vars_BVV_.at(i)->name()==name){
      vars_BVV_ .at(i)->push(values) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_BV_.size() ; i++){
    if(vars_BV_.at(i)->name()==name){
      for(unsigned j=0 ; j<values.size() ; j++){
        vars_BV_ .at(i)->push(values.at(j)) ;
      }
      return true ;
    }
  }
  if(debug) std::cout << "Could not find a (vector bool) branch named " << name << std::endl ;
  return false ;
}
bool HEEPTree::store(std::string name, std::vector<int> values){
  for(unsigned int i=0 ; i<vars_IVV_.size() ; i++){
    if(vars_IVV_.at(i)->name()==name){
      vars_IVV_ .at(i)->push(values) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_IV_.size() ; i++){
    if(vars_IV_.at(i)->name()==name){
      for(unsigned j=0 ; j<values.size() ; j++){
        vars_IV_ .at(i)->push(values.at(j)) ;
      }
      return true ;
    }
  }
  if(debug) std::cout << "Could not find a (vector int) branch named " << name << std::endl ;
  return false ;
}
bool HEEPTree::store(std::string name, std::vector<float> values){
  for(unsigned int i=0 ; i<vars_FVV_.size() ; i++){
    if(vars_FVV_.at(i)->name()==name){
      vars_FVV_ .at(i)->push(values) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_FV_.size() ; i++){
    if(vars_FV_.at(i)->name()==name){
      for(unsigned j=0 ; j<values.size() ; j++){
        vars_FV_ .at(i)->push(values.at(j)) ;
      }
      return true ;
    }
  }
  if(debug) std::cout << "Could not find a (vector float) branch named " << name << std::endl ;
  return false ;
}
bool HEEPTree::store(std::string name, std::vector<double> values){
  for(unsigned int i=0 ; i<vars_DVV_.size() ; i++){
    if(vars_DVV_.at(i)->name()==name){
      vars_DVV_ .at(i)->push(values) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_DV_.size() ; i++){
    if(vars_DV_.at(i)->name()==name){
      for(unsigned j=0 ; j<values.size() ; j++){
        vars_DV_ .at(i)->push(values.at(j)) ;
      }
      return true ;
    }
  }
  if(debug) std::cout << "Could not find a (vector double) branch named " << name << std::endl ;
  return false ;
}


DEFINE_FWK_MODULE(HEEPTree);
