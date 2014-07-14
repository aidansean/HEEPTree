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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "UserCode/HEEPSkims/interface/BranchWrapper.h"

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
  
  bool store(std::string, bool  );
  bool store(std::string, double);
  bool store(std::string, float );
  bool store(std::string, int   );
  bool store(std::string, unsigned);
  bool store(std::string, std::vector<bool>  );
  bool store(std::string, std::vector<double>);
  bool store(std::string, std::vector<float >);
  bool store(std::string, std::vector<int   >);
  
  bool add_branch(std::string);
  bool add_branch(std::string,int);
  bool branch_exists(std::string);
  
  void set_branch_type(int);
  int  get_branch_type();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  
  void beginEvent() ;
  void endEvent() ;
  
  void configure_branches();
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::BeamSpot> beamSpotLabel_ ;
  double ScPtMin_;
  double GsfPtMin_;
  double GsfTrackPtMin_;
  double muPtMin_;
  
  std::vector<BranchWrapperBase*> all_vars_;
  std::vector<BranchWrapperBVV* > vars_BVV_;
  std::vector<BranchWrapperDVV* > vars_DVV_;
  std::vector<BranchWrapperFVV* > vars_FVV_;
  std::vector<BranchWrapperIVV* > vars_IVV_;
  std::vector<BranchWrapperBV*  > vars_BV_ ;
  std::vector<BranchWrapperDV*  > vars_DV_ ;
  std::vector<BranchWrapperFV*  > vars_FV_ ;
  std::vector<BranchWrapperIV*  > vars_IV_ ;
  std::vector<BranchWrapperB*   > vars_B_  ;
  std::vector<BranchWrapperD*   > vars_D_  ;
  std::vector<BranchWrapperF*   > vars_F_  ;
  std::vector<BranchWrapperI*   > vars_I_  ;
  
  int current_var_type;
  std::vector< std::pair<std::string, int> > list_of_branches ;
  std::vector< std::pair<std::string, int> > missing_branches ;
    
  // Parameters for PU subtraction 
  double EcalHcal1EffAreaBarrel_  ;
  double EcalHcal1EffAreaEndcaps_ ;
  
  bool debug;
  
  float rho ;

  // config parameters -------------------------------
  TFile* myFile ;
  TTree* mytree ;
};
#endif
//define this as a plug-in

