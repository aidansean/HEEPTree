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
//Cleaning ladies : Thomas and Laurent
#include <iostream>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "UserCode/HEEPSkims/interface/HEEPTree.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
//ECAL 
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include <TMath.h>
#include <vector>
#define PI 3.141592654
#define TWOPI 6.283185308

using namespace std;
using namespace reco;
using namespace edm;

#include "UserCode/HEEPSkims/interface/HEEPTree.h"
#include "UserCode/HEEPSkims/interface/BranchWrapper.h"


HEEPTree::HEEPTree(const edm::ParameterSet& iConfig)
{
}


HEEPTree::~HEEPTree()
{
}

// ------------ method called once each job just before starting event loop  ------------
void 
HEEPTree::beginJob()
{
  edm::Service<TFileService> fs;
  myFile = new TFile("outfile.root", "RECREATE") ;
  mytree = new TTree("HEEPTree", "HEEPTree") ;
  
  vars_FV_.push_back(new branch_wrapper_FV("el_pt"    )) ;
  vars_FV_.push_back(new branch_wrapper_FV("el_px"    )) ;
  vars_FV_.push_back(new branch_wrapper_FV("el_py"    )) ;
  vars_FV_.push_back(new branch_wrapper_FV("el_px"    )) ;
  vars_FV_.push_back(new branch_wrapper_FV("el_E"     )) ;
  vars_FV_.push_back(new branch_wrapper_FV("el_eta"   )) ;
  vars_FV_.push_back(new branch_wrapper_FV("el_phi"   )) ;
  vars_IV_.push_back(new branch_wrapper_IV("el_charge")) ;
  
  vars_FV_.push_back(new branch_wrapper_FV("el_hovere"  )) ;
  vars_FV_.push_back(new branch_wrapper_FV("el_preshowerEnergy")) ;
  
  vars_I_ .push_back(new branch_wrapper_I ("eventNumber"    )) ;
  vars_I_ .push_back(new branch_wrapper_I ("runNumber"      )) ;
  vars_I_ .push_back(new branch_wrapper_I ("luminosityBlock")) ;
  
  for(unsigned int i=0 ; i<vars_D_ .size() ; i++){ vars_D_ .at(i)->config(mytree) ; }
  for(unsigned int i=0 ; i<vars_F_ .size() ; i++){ vars_F_ .at(i)->config(mytree) ; }
  for(unsigned int i=0 ; i<vars_I_ .size() ; i++){ vars_I_ .at(i)->config(mytree) ; }
  for(unsigned int i=0 ; i<vars_DV_.size() ; i++){ vars_DV_.at(i)->config(mytree) ; }
  for(unsigned int i=0 ; i<vars_FV_.size() ; i++){ vars_FV_.at(i)->config(mytree) ; }
  for(unsigned int i=0 ; i<vars_IV_.size() ; i++){ vars_FV_.at(i)->config(mytree) ; }
}

// ------------ method called to for each event  ------------
void
HEEPTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  beginEvent() ;
  
  //Run and event number
  store("eventNumber"    , (int)(iEvent.id().event()          )) ;
  store("runNumber"      , (int)(iEvent.id().run()            )) ;
  store("luminosityBlock", (int)(iEvent.id().luminosityBlock())) ;
  
  //Final GSF Electron collection
  edm::Handle<reco::GsfElectronCollection> pGsfElectrons;
  iEvent.getByLabel("gsfElectrons","",pGsfElectrons);
  reco::GsfElectronCollection gsfelectrons(pGsfElectrons->begin(),pGsfElectrons->end());
  
  for(reco::GsfElectronCollection::const_iterator gsfiter = gsfelectrons.begin() ; gsfiter != gsfelectrons.end(); ++gsfiter){
    //Fill the gsf related variables
    store("el_pt"      , gsfiter->pt()    ) ;
    store("el_px"      , gsfiter->px()    ) ;
    store("el_py"      , gsfiter->py()    ) ;
    store("el_pz"      , gsfiter->pz()    ) ;
    store("el_E"       , gsfiter->energy()) ;
    store("el_eta"     , gsfiter->eta()   ) ;
    store("el_phi"     , gsfiter->phi()   ) ;
    store("el_charge"  , gsfiter->charge()) ;
    store("el_hovere"  , gsfiter->hadronicOverEm()) ;
    store("el_preshowerEnergy", gsfiter->superCluster()->preshowerEnergy()) ;
  }
  
  mytree->Fill();
  endEvent() ;
 }//end of analyze method


void 
HEEPTree::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

void HEEPTree::beginEvent(){
  for(unsigned int i=0 ; i<vars_D_ .size() ; i++){ vars_D_ .at(i)->event_begin() ; }
  for(unsigned int i=0 ; i<vars_F_ .size() ; i++){ vars_F_ .at(i)->event_begin() ; }
  for(unsigned int i=0 ; i<vars_I_ .size() ; i++){ vars_I_ .at(i)->event_begin() ; }
  for(unsigned int i=0 ; i<vars_DV_.size() ; i++){ vars_DV_.at(i)->event_begin() ; }
  for(unsigned int i=0 ; i<vars_FV_.size() ; i++){ vars_FV_.at(i)->event_begin() ; }
  for(unsigned int i=0 ; i<vars_IV_.size() ; i++){ vars_IV_.at(i)->event_begin() ; }
}

void HEEPTree::endEvent(){
  for(unsigned int i=0 ; i<vars_D_ .size() ; i++){ vars_D_ .at(i)->event_end() ; }
  for(unsigned int i=0 ; i<vars_F_ .size() ; i++){ vars_F_ .at(i)->event_end() ; }
  for(unsigned int i=0 ; i<vars_I_ .size() ; i++){ vars_I_ .at(i)->event_end() ; }
  for(unsigned int i=0 ; i<vars_DV_.size() ; i++){ vars_DV_.at(i)->event_end() ; }
  for(unsigned int i=0 ; i<vars_FV_.size() ; i++){ vars_FV_.at(i)->event_end() ; }
  for(unsigned int i=0 ; i<vars_IV_.size() ; i++){ vars_IV_.at(i)->event_end() ; }
}


// ------------ method called once each job just after ending the event loop  ------------
void 
HEEPTree::endJob() {
  if(myFile){
    myFile->Write() ;
    delete myFile ;
  }
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
  
  return false ;
}
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
  return false ;
}

DEFINE_FWK_MODULE(HEEPTree);
