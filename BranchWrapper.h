#include "TTree.h"
#include <vector>
#include <iostream>

#ifndef BRANCHWRAPPER
#define BRANCHWRAPPER

// Inherited
class BranchWrapperBase{
  public:
    BranchWrapperBase(std::string);
    ~BranchWrapperBase();
    void beginEvent();
    void endEvent();
    std::string name(){ return name_; } ;
    bool  is_filled(){ return is_filled_ ; } ;
    bool is_touched(){ return is_touched_; } ;
    void touch(){ is_touched_ = true ; } ;
    void fill(){ is_filled_ = true ; touch() ; } ;
    void unfill(){ is_filled_ = false ; } ;
    virtual int config(TTree*){return -1; } ;
  private:
    std::string name_ ;
    bool is_filled_;
    bool is_touched_;
};

class BranchWrapperB: public BranchWrapperBase{
  private:
    bool value_;
  public:
    BranchWrapperB(std::string) ;
    ~BranchWrapperB(){} ;
    void set(bool);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};

class BranchWrapperD: public BranchWrapperBase{
  private:
    double value_;
  public:
    BranchWrapperD(std::string) ;
    ~BranchWrapperD(){} ;
    void set(double);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};


class BranchWrapperF: public BranchWrapperBase{
  private:
    float value_;
  public:
    BranchWrapperF(std::string) ;
    ~BranchWrapperF(){} ;
    void set(float);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};

class BranchWrapperI: public BranchWrapperBase{
  private:
    int value_;
  public:
    BranchWrapperI(std::string) ;
    ~BranchWrapperI(){} ;
    void set(int);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};

class BranchWrapperBV: public BranchWrapperBase{
  private:
    std::vector<bool> values_;
  public:
    BranchWrapperBV(std::string) ;
    ~BranchWrapperBV() ;
    void push(bool);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};

class BranchWrapperDV: public BranchWrapperBase{
  private:
    std::vector<double> values_;
  public:
    BranchWrapperDV(std::string) ;
    ~BranchWrapperDV() ;
    void push(double);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};

class BranchWrapperFV: public BranchWrapperBase{
  private:
    std::vector<float> values_;
  public:
    BranchWrapperFV(std::string) ;
    ~BranchWrapperFV() ;
    void push(float);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};

class BranchWrapperIV: public BranchWrapperBase{
  private:
    std::vector<int> values_;
  public:
    BranchWrapperIV(std::string) ;
    ~BranchWrapperIV() ;
    void push(int);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};


class BranchWrapperBVV: public BranchWrapperBase{
  private:
    std::vector<std::vector<bool> > values_;
  public:
    BranchWrapperBVV(std::string) ;
    ~BranchWrapperBVV() ;
    void push(std::vector<bool>);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};

class BranchWrapperDVV: public BranchWrapperBase{
  private:
    std::vector<std::vector<double> > values_;
  public:
    BranchWrapperDVV(std::string) ;
    ~BranchWrapperDVV() ;
    void push(std::vector<double>);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};

class BranchWrapperFVV: public BranchWrapperBase{
  private:
    std::vector<std::vector<float> > values_;
  public:
    BranchWrapperFVV(std::string) ;
    ~BranchWrapperFVV() ;
    void push(std::vector<float>);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};

class BranchWrapperIVV: public BranchWrapperBase{
  private:
    std::vector<std::vector<int> > values_;
  public:
    BranchWrapperIVV(std::string) ;
    ~BranchWrapperIVV() ;
    void push(std::vector<int>);
    int config(TTree*);
    void beginEvent();
    void endEvent();
};


// Templated (not used, yet)
template <class T>
class branchWrapper_simple{
    std::string name_ ;
    T value_ ;
  public:
    std::string name(){ return name_ ; } ;
    branchWrapper_simple(std::string) ;
    ~branchWrapper_simple(){} ;
    void set(T) ;
    T get() ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};
template <class T>
class branchWrapper_vector{
    std::string name_ ;
    std::vector<T> values_ ;
  public:
    std::string name(){ return name_ ; } ;
    branchWrapper_vector(std::string) ;
    ~branchWrapper_vector(){} ;
    bool push(T) ;
    T get(int) ;
    int nEntries() ;
    int config(TTree*) ;
    bool beginEvent() ;
    bool endEvent() ;
};


#endif