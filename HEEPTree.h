#ifndef UserCode_HEEPSkims_HEEPTree_h
#define UserCode_HEEPSkims_HEEPTree_h

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
// $Id: HEEPTree.h,v 1.42 2013/03/15 16:23:54 treis Exp $
//
//

// system include files
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "UserCode/HEEPSkims/interface/BranchWrapper.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1F.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

//
// class decleration
//
class HEEPTree : public edm::EDAnalyzer {

private:

public:
  explicit HEEPTree(const edm::ParameterSet& iConfig);
  ~HEEPTree();
  
  bool store(std::string,double);
  bool store(std::string,float );
  bool store(std::string,int   );
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  
  void beginEvent() ;
  void endEvent() ;
  
  std::vector<branch_wrapper_DV*> vars_DV_;
  std::vector<branch_wrapper_FV*> vars_FV_;
  std::vector<branch_wrapper_IV*> vars_IV_;
  std::vector<branch_wrapper_D* > vars_D_;
  std::vector<branch_wrapper_F* > vars_F_;
  std::vector<branch_wrapper_I* > vars_I_;
  // ----------member data ---------------------------

  // config parameters -------------------------------
  TFile* myFile ;
  TTree* mytree ;
};
#endif
//define this as a plug-in

