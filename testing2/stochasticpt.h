#ifndef STOCHASTICPT_H
#define STOCHASTICPT_H

#include <vector>
#include "global.h"
#include "timer.h"
#include "sampling.h"
#include <bitset>

  class StateQuantaInfo
  {
    struct leftquantainfo
    {
      int first; // The number of left quanta.
      int second;// The number of states in the quanta.
      int third; // After collecting combined quanta from left quanta and dot quanta, the position of this quanta in combined quanta.
      leftquantainfo() :first(0), second(0), third(0){};
      leftquantainfo(int a, int b, int c) :first(a), second(b), third(c){};
    
      friend class boost::serialization::access;
      template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
          ar & first \
            & second\
            & third;
        }
    
    };
    typedef std::pair<int,int> intpair;
    int currentstate;
    std::vector< std::map < intpair, intpair> > quantalink;
    std::vector< std::map < intpair, leftquantainfo> > reversequantalink;
    std::vector< std::vector<Matrix> > RotationMatrices;
    public:
    double norm;
    double mpssampleT=0.0;
    double mpscoeffT=0.0;
    double HdiagonalT=0.0;

    StateQuantaInfo(int state=0): currentstate(state), quantalink(), reversequantalink(), RotationMatrices(){};
    void StoreQuantaInformation(SweepParams &sweepParams, const bool &forward);
    void readRotationandQuanta();
    double getcoeff(const std::vector<int>& ci);
    double getcoeff(const bitstring& ci);
    void preparestate();
    void initstate();
    void get_slater();
    void test_expectation();
    double sampling(std::vector<int>& ci);
    double sampling_approx(std::vector<int>& ci, double& prob);
    double sampling(bitstring& ci);
    double sampling_approx(bitstring& ci, double& prob);
    double local_energy(const std::vector<int>& ci, int integralIndex=0);
    double local_energy(const longbitarray& ci, int integralIndex=0);
    double local_energy(const bitstring& ci, int integralIndex=0);
    //static double local_energy(const std::vector<int>& ci, int integralIndex=0);
  };

template <class T1, class T2>
void insertorsum(std::unordered_map<T1,T2>& hashtable, const T1& key, const T2& value)
{
  if(hashtable.find(key)==hashtable.end())
  {
    hashtable[key] = value;
  }
  else{
    hashtable[key] +=value;
  }
}
#endif
