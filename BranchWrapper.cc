#include "UserCode/HEEPSkims/interface/BranchWrapper.h"
#include <iostream>
#include <vector>

// Templated
template <class T>
branchWrapper_simple<T>::branchWrapper_simple(std::string name){
  name_  = name ;
  value_ = -999 ;
  std::cout << typeid(T).name() << std::endl ;
}
template <class T>
int branchWrapper_simple<T>::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name_.c_str())) return 2 ;
  tree->Branch(name_.c_str(), &value_, Form("%s/F", name_.c_str())) ;
  return 0 ;
}
template <class T> void branchWrapper_simple<T>::event_begin(){ value_ = -999 ; }
template <class T> void branchWrapper_simple<T>::event_end(){}
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
  if(tree->GetBranch(name_.c_str())) return 2 ;
  tree->Branch(name_.c_str(), &values_) ;
  return 0 ;
}
template <class T> bool branchWrapper_vector<T>::event_begin(){
  if(!values_) return false ;
  values_->clear() ;
  return true ;
}
template <class T> bool branchWrapper_vector<T>::event_end(){
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

// double
branch_wrapper_D::branch_wrapper_D(std::string name){
  name_  = name   ;
  value_ = -999 ;
}
int branch_wrapper_D::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name_.c_str())) return 2 ;
  tree->Branch(name_.c_str(), &value_, Form("%s/F", name_.c_str())) ;
  return 0 ;
}
void branch_wrapper_D::set(double value){ value_ = value ; }
void branch_wrapper_D::event_begin(){ value_ = -999 ; }
void branch_wrapper_D::event_end(){}


// float
branch_wrapper_F::branch_wrapper_F(std::string name){
  name_  = name   ;
  value_ = -999 ;
}
int  branch_wrapper_F::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name_.c_str())) return 2 ;
  tree->Branch(name_.c_str(), &value_, Form("%s/F", name_.c_str())) ;
  return 0 ;
}
void branch_wrapper_F::set(float value){ value_ = value ; }
void branch_wrapper_F::event_begin(){ value_ = -999 ; }
void branch_wrapper_F::event_end(){}

// int
branch_wrapper_I::branch_wrapper_I(std::string name){
  name_  = name   ;
  value_ = -999 ;
}
int  branch_wrapper_I::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name_.c_str())) return 2 ;
  tree->Branch(name_.c_str(), &value_, Form("%s/I", name_.c_str())) ;
  return 0 ;
}
void branch_wrapper_I::set(int value){ value_ = value ; }
void branch_wrapper_I::event_begin(){ value_ = -999 ; }
void branch_wrapper_I::event_end(){}


// Vector of doubles
branch_wrapper_DV::branch_wrapper_DV(std::string name){
  name_   = name   ;
}
branch_wrapper_DV::~branch_wrapper_DV(){}
int  branch_wrapper_DV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name_.c_str())) return 2 ;
  tree->Branch(name_.c_str(), &values_) ;
  return 0 ;
}
void branch_wrapper_DV::push(double value){ values_.push_back(value) ; }
void branch_wrapper_DV::event_begin(){ values_.clear() ; }
void branch_wrapper_DV::event_end(){}

// Vector of floats
branch_wrapper_FV::branch_wrapper_FV(std::string name){
  name_   = name   ;
}
branch_wrapper_FV::~branch_wrapper_FV(){}
int  branch_wrapper_FV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name_.c_str())) return 2 ;
  tree->Branch(name_.c_str(), &values_) ;
  return 0 ;
}
void branch_wrapper_FV::push(float value){ values_.push_back(value) ; }
void branch_wrapper_FV::event_begin(){ values_.clear() ; }
void branch_wrapper_FV::event_end(){}

// Vector of ints
branch_wrapper_IV::branch_wrapper_IV(std::string name){
  name_   = name   ;
}
branch_wrapper_IV::~branch_wrapper_IV(){}
int  branch_wrapper_IV::config(TTree* tree){
  if(!tree) return 1 ;
  if(tree->GetBranch(name_.c_str())) return 2 ;
  tree->Branch(name_.c_str(), &values_) ;
  return 0 ;
}
void branch_wrapper_IV::push(int value){ values_.push_back(value) ; }
void branch_wrapper_IV::event_begin(){ values_.clear() ; }
void branch_wrapper_IV::event_end(){}
