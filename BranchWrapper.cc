#include "UserCode/HEEPSkims/interface/BranchWrapper.h"
#include <iostream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////
//                                  Inherited classes                                 //
////////////////////////////////////////////////////////////////////////////////////////
// Inherited
BranchWrapperBase::BranchWrapperBase(std::string name){
  name_       = name ;
  is_filled_  = false ;
  is_touched_ = false ;
}
BranchWrapperBase::~BranchWrapperBase(){}
void BranchWrapperBase::beginEvent(){
  is_filled_ = false ;
}
void BranchWrapperBase::endEvent(){}

////////////////////////////////////////////////////////////////////////////////////////
//                                    Simple classes                                  //
////////////////////////////////////////////////////////////////////////////////////////
BranchWrapperB::BranchWrapperB(std::string name): BranchWrapperBase(name){
  value_ = false ;
}
int BranchWrapperB::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &value_, Form("%s/O", name().c_str())) ;
  return 0 ;
}
void BranchWrapperB::set(bool value){
  value_ = value ;
  fill() ;
}
void BranchWrapperB::beginEvent(){
  value_ = false ;
  unfill() ;
}
void BranchWrapperB::endEvent(){}

// double
BranchWrapperD::BranchWrapperD(std::string name): BranchWrapperBase(name){
  value_ = -999 ;
}
int BranchWrapperD::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &value_, Form("%s/D", name().c_str())) ;
  return 0 ;
}
void BranchWrapperD::set(double value){
  value_ = value ;
  fill() ;
}
void BranchWrapperD::beginEvent(){
  value_ = -999 ;
  unfill() ;
}
void BranchWrapperD::endEvent(){}


// float
BranchWrapperF::BranchWrapperF(std::string name): BranchWrapperBase(name){
  value_ = -999 ;
}
int BranchWrapperF::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &value_, Form("%s/F", name().c_str())) ;
  return 0 ;
}
void BranchWrapperF::set(float value){
  value_ = value ;
  fill() ;
}
void BranchWrapperF::beginEvent(){
  value_ = -999 ;
  unfill() ;
}
void BranchWrapperF::endEvent(){}

// int
BranchWrapperI::BranchWrapperI(std::string name): BranchWrapperBase(name){
  value_ = -999 ;
}
int BranchWrapperI::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &value_, Form("%s/I", name().c_str())) ;
  return 0 ;
}
void BranchWrapperI::set(int value){
  value_ = value ;
  fill() ;
}
void BranchWrapperI::beginEvent(){
  value_ = -999 ;
  unfill() ;
}
void BranchWrapperI::endEvent(){}

////////////////////////////////////////////////////////////////////////////////////////
//                                    Vector classes                                  //
////////////////////////////////////////////////////////////////////////////////////////
// Vector of bools
BranchWrapperBV::BranchWrapperBV(std::string name): BranchWrapperBase(name){}
BranchWrapperBV::~BranchWrapperBV(){}
int BranchWrapperBV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &values_) ;
  return 0 ;
}
void BranchWrapperBV::push(bool value){
  values_.push_back(value) ;
  fill() ;
}
void BranchWrapperBV::beginEvent(){
  unfill() ;
  values_.clear() ;
}
void BranchWrapperBV::endEvent(){}

// Vector of doubles
BranchWrapperDV::BranchWrapperDV(std::string name): BranchWrapperBase(name){}
BranchWrapperDV::~BranchWrapperDV(){}
int BranchWrapperDV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &values_) ;
  return 0 ;
}
void BranchWrapperDV::push(double value){
  values_.push_back(value) ;
  fill();
}
void BranchWrapperDV::beginEvent(){
  unfill();
  values_.clear() ;
}
void BranchWrapperDV::endEvent(){}

// Vector of floats
BranchWrapperFV::BranchWrapperFV(std::string name): BranchWrapperBase(name){}
BranchWrapperFV::~BranchWrapperFV(){}
int BranchWrapperFV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &values_) ;
  return 0 ;
}
void BranchWrapperFV::push(float value){
  values_.push_back(value) ;
  fill() ;
}
void BranchWrapperFV::beginEvent(){
  unfill();
  values_.clear() ;
}
void BranchWrapperFV::endEvent(){}

// Vector of ints
BranchWrapperIV::BranchWrapperIV(std::string name): BranchWrapperBase(name){}
BranchWrapperIV::~BranchWrapperIV(){}
int BranchWrapperIV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &values_) ;
  return 0 ;
}
void BranchWrapperIV::push(int value){
  values_.push_back(value) ;
  fill() ;
}
void BranchWrapperIV::beginEvent(){
  unfill();
  values_.clear() ;
}
void BranchWrapperIV::endEvent(){}


////////////////////////////////////////////////////////////////////////////////////////
//                                Vector vector classes                               //
////////////////////////////////////////////////////////////////////////////////////////
// Vector of vector of bools
BranchWrapperBVV::BranchWrapperBVV(std::string name): BranchWrapperBase(name){}
BranchWrapperBVV::~BranchWrapperBVV(){}
int BranchWrapperBVV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &values_) ;
  return 0 ;
}
void BranchWrapperBVV::push(std::vector<bool> value){
  values_.push_back(value) ;
  fill();
}
void BranchWrapperBVV::beginEvent(){
  unfill();
  values_.clear() ;
}
void BranchWrapperBVV::endEvent(){}


// Vector of vector of doubles
BranchWrapperDVV::BranchWrapperDVV(std::string name): BranchWrapperBase(name){}
BranchWrapperDVV::~BranchWrapperDVV(){}
int BranchWrapperDVV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &values_) ;
  return 0 ;
}
void BranchWrapperDVV::push(std::vector<double> value){
  values_.push_back(value) ;
  fill();
}
void BranchWrapperDVV::beginEvent(){
  unfill();
  values_.clear() ;
}
void BranchWrapperDVV::endEvent(){}

// Vector of vector of floats
BranchWrapperFVV::BranchWrapperFVV(std::string name): BranchWrapperBase(name){}
BranchWrapperFVV::~BranchWrapperFVV(){}
int BranchWrapperFVV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &values_) ;
  return 0 ;
}
void BranchWrapperFVV::push(std::vector<float> value){
  values_.push_back(value) ;
  fill() ;
}
void BranchWrapperFVV::beginEvent(){
  unfill();
  values_.clear() ;
}
void BranchWrapperFVV::endEvent(){}

// Vector of vector of ints
BranchWrapperIVV::BranchWrapperIVV(std::string name): BranchWrapperBase(name){}
BranchWrapperIVV::~BranchWrapperIVV(){}
int BranchWrapperIVV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &values_) ;
  return 0 ;
}
void BranchWrapperIVV::push(std::vector<int> value){
  values_.push_back(value) ;
  fill() ;
}
void BranchWrapperIVV::beginEvent(){
  unfill() ;
  values_.clear() ;
}
void BranchWrapperIVV::endEvent(){}

////////////////////////////////////////////////////////////////////////////////////////
//                             Templated classes (not used)                           //
////////////////////////////////////////////////////////////////////////////////////////
template <class T>
branchWrapper_simple<T>::branchWrapper_simple(std::string name){
  name_  = name ;
  value_ = -999 ;
}
template <class T>
int branchWrapper_simple<T>::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &value_, Form("%s/F", name().c_str())) ;
  return 0 ;
}
template <class T> void branchWrapper_simple<T>::beginEvent(){ value_ = -999 ; }
template <class T> void branchWrapper_simple<T>::endEvent(){}
template <class T> void branchWrapper_simple<T>::set(T value){ value_ = value ; }
template <class T> T    branchWrapper_simple<T>::get(){ return value_ ; }

template <class T>
branchWrapper_vector<T>::branchWrapper_vector(std::string name){
  name_   = name ;
  values_ = 0 ;
  std::cout << typeid(T).name() << std::endl ;
}
template <class T>
int branchWrapper_vector<T>::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name().c_str())) return 2 ;
  tree->Branch(name().c_str(), &values_) ;
  return 0 ;
}
template <class T> bool branchWrapper_vector<T>::beginEvent(){
  if(!values_) return false ;
  values_->clear() ;
  return true ;
}
template <class T> bool branchWrapper_vector<T>::endEvent(){
  if(!values_) return false ;
  return true ;
}
template <class T> bool branchWrapper_vector<T>::push(T value){
  if(!values_) return false ;
  values_->push_back(value) ;
  return true ;
}
template <class T> T    branchWrapper_vector<T>::get(int index){
  if(!values_) return -999 ;
  return values_->at(index) ;
}
template <class T> int  branchWrapper_vector<T>::nEntries(){
  if(!values_) return -1 ;
  return values_->size() ;
}

