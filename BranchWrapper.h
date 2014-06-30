#include "TTree.h"
#include <vector>

#ifndef BRANCHWRAPPER
#define BRANCHWRAPPER

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
    void event_begin() ;
    void event_end() ;
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
    bool event_begin() ;
    bool event_end() ;
};


class branch_wrapper_D{
  private:
    std::string name_;
    float value_;
  public:
    std::string name(){ return name_;};
    branch_wrapper_D(std::string) ;
    ~branch_wrapper_D(){} ;
    void set(double);
    int config(TTree*);
    void event_begin();
    void event_end();
};

class branch_wrapper_F{
  private:
    std::string name_;
    float value_;
  public:
    std::string name(){ return name_;};
    branch_wrapper_F(std::string) ;
    ~branch_wrapper_F(){} ;
    void set(float);
    int config(TTree*);
    void event_begin();
    void event_end();
};

class branch_wrapper_I{
  private:
    std::string name_;
    int value_;
  public:
    std::string name(){ return name_;};
    branch_wrapper_I(std::string) ;
    ~branch_wrapper_I(){} ;
    void set(int);
    int config(TTree*);
    void event_begin();
    void event_end();
};

class branch_wrapper_DV{
  private:
    std::string name_;
    std::vector<double> values_;
  public:
    std::string name(){ return name_;};
    branch_wrapper_DV(std::string) ;
    ~branch_wrapper_DV() ;
    void push(double);
    int config(TTree*);
    void event_begin();
    void event_end();
};

class branch_wrapper_FV{
  private:
    std::string name_;
    std::vector<float> values_;
  public:
    std::string name(){ return name_;};
    branch_wrapper_FV(std::string) ;
    ~branch_wrapper_FV() ;
    void push(float);
    int config(TTree*);
    void event_begin();
    void event_end();
};

class branch_wrapper_IV{
  private:
    std::string name_;
    std::vector<int> values_;
  public:
    std::string name(){ return name_;};
    branch_wrapper_IV(std::string) ;
    ~branch_wrapper_IV() ;
    void push(int);
    int config(TTree*);
    void event_begin();
    void event_end();
};
#endif