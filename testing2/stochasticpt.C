#include "MatrixBLAS.h"
#include "spinblock.h"
#include "initblocks.h"
#include "input.h"
#include "timer.h"
#include <ctime>
#include "rotationmat.h"
#include "wavefunction.h"
#include "global.h"
#include "sweep.h"
#include <unordered_map>
#include <unordered_set>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <random>
#include "IntegralMatrix.h"
#include <chrono>
#include "stochasticpt.h"
#include "sampling.h"
#include "heatbath.h"
#include "nonspinmps.h"
using namespace SpinAdapted;
using namespace SpinAdapted::Sweep;
int bitstring::n_orb;
void ReadInput(char* conf);

void getComplementarySites(std::vector<int> &sites, std::vector<int> &complementarySites) ;


void StateQuantaInfo::StoreQuantaInformation(SweepParams &sweepParams, const bool &forward)
{
    bool restart=false, warmUp = false;
    std::vector<int> sites, spinsites;

    int new_site, wave_site;
    new_site = 0;
    sites.push_back(new_site);
    if (dmrginp.spinAdapted())
      spinsites.push_back(new_site);
    else {
      spinsites.push_back(2*new_site);
      spinsites.push_back(2*new_site+1);
    }
    StateInfo stateinfo; 
    std::vector<Matrix> RotationMatrix;

    makeStateInfo(stateinfo, new_site);
    //StateInfo::restore(forward, sites, stateinfo, sweepParams.current_root());
    //LoadRotationMatrix (spinsites, RotationMatrix, currentstate);

    for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
      new_site++;
      sites.push_back(new_site);
      if (dmrginp.spinAdapted())
        spinsites.push_back(new_site);
      else {
        spinsites.push_back(2*new_site);
        spinsites.push_back(2*new_site+1);
      }

      pout <<"Sites: ";
      for (int i : sites)
        pout << i<<" ";
      pout <<endl;

      StateInfo combinedstateinfo, siteState;
      StateInfo leftStateInfo=stateinfo;
      //TensorProduct provides a pointer to left and right StateInfo.
      //Load Stateinfo will overwrite the left StateInfo pointed.

      makeStateInfo(siteState, new_site);
      TensorProduct(leftStateInfo, siteState, combinedstateinfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

      combinedstateinfo.CollectQuanta();
      StateInfo uncollectedstateinfo = *(combinedstateinfo.unCollectedStateInfo);
      pout << uncollectedstateinfo <<endl;
      pout << combinedstateinfo<<endl;

      LoadRotationMatrix (spinsites, RotationMatrix, currentstate);

      SpinAdapted::StateInfo::transform_state(RotationMatrix, combinedstateinfo, stateinfo);
      //StateInfo::restore(forward, sites, stateinfo, sweepParams.current_root());
      pout << stateinfo <<endl;

      //typedef std::unordered_map < std::pair<int,int>, std::pair<int,int> > 
      typedef std::pair<int,int> intpair;
      std::map < intpair, intpair> combinedquanta;
      std::map < intpair, leftquantainfo> reversequanta;
      //Print A^{n_i} dimension.
      for (int i=0; i<combinedstateinfo.quanta.size();i++)
      {
        //Some quanta num discarded after renormalization.
        auto iter = std::find(stateinfo.quanta.begin(), stateinfo.quanta.end(), combinedstateinfo.quanta[i]);
        if (iter !=stateinfo.quanta.end() )
        {
          int index = std::distance( stateinfo.quanta.begin(), iter );
          int quantastate =0;

          for(int k=0; k < combinedstateinfo.oldToNewState[i].size(); k++)
          {
            int oldquanta = combinedstateinfo.oldToNewState[i][k];
            std::pair<int,int> combinepair (uncollectedstateinfo.leftUnMapQuanta[oldquanta],uncollectedstateinfo.rightUnMapQuanta[oldquanta]);


            //TODO
            //Discarded quanta is still in Rotation matrix.
            //std::pair<int,int> pointedpair(index,quantastate);
            std::pair<int,int> pointedpair(index,quantastate);
            combinedquanta[combinepair] = pointedpair;
            std::pair<int,int> d_r_quanta (uncollectedstateinfo.rightUnMapQuanta[oldquanta],index);
            int leftquantastates = leftStateInfo.quantaStates[uncollectedstateinfo.leftUnMapQuanta[oldquanta] ];
            leftquantainfo left_index(uncollectedstateinfo.leftUnMapQuanta[oldquanta], leftquantastates, quantastate);
            reversequanta[d_r_quanta] = left_index;
            quantastate += uncollectedstateinfo.quantaStates[oldquanta];
          }
        }
      }
      //for (auto it =combinedquanta.begin(); it !=combinedquanta.end();it++)
      //  pout <<it->first.first <<" "<<it->first.second<<"-> " <<it->second.first<<" "<<it->second.second<<endl;


      {
        std::string file;
        file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/quantalink-state" % new_site %"-"%currentstate %".tmp" );
        std::ofstream s(file.c_str());
        boost::archive::text_oarchive oa(s);
        oa << combinedquanta;
        s.close();
      }
      {
        std::string file;
        file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/reversequantalink-state" % new_site %"-"%currentstate %".tmp" );
        std::ofstream s(file.c_str());
        boost::archive::text_oarchive oa(s);
        oa << reversequanta;
        s.close();
      }
    

      ++sweepParams.set_block_iter();
    }
    new_site++;

    //TODO
    //The wave function.
    std::vector<int> complementarySites;
    if (dmrginp.spinAdapted())
    {
      complementarySites.assign(1,dmrginp.last_site()-1);
    }
    else
    {
      complementarySites.push_back(dmrginp.last_site()/2-1);

    }
    //wave function and rotation matrix files are named by spin orbitals.
    //StateInfo files are named by spatial orbitals.
    getComplementarySites(spinsites,complementarySites);

    Wavefunction lastsitewave;
    StateInfo waveinfo;
    //TODO
    lastsitewave.LoadWavefunctionInfo (waveinfo, spinsites, currentstate);
    //norm = DotProduct(lastsitewave,lastsitewave);

    pout <<"Last site state info " <<endl;
    pout <<stateinfo<<endl;
    pout <<waveinfo<<endl;
    pout <<*(waveinfo.leftStateInfo)<<endl;
    pout <<*(waveinfo.rightStateInfo)<<endl;
    pout <<lastsitewave<<endl;

    //TODO
    //Assuming quanta num waveinfo is not collected.
    //There is one entangle bond between the last site and other sites in each quanta num.
    typedef std::pair<int,int> intpair;
    std::map < intpair, intpair> combinedquanta;
    std::map < intpair, leftquantainfo> reversequanta;

    StateInfo combinedstateinfo, siteState;
    makeStateInfo(siteState, new_site);
    TensorProduct(stateinfo, siteState, combinedstateinfo, PARTICLE_SPIN_NUMBER_CONSTRAINT);


    for (int i=0; i<combinedstateinfo.quanta.size();i++)
    {
      int quantastate =0;

      std::pair<int,int> combinepair (combinedstateinfo.leftUnMapQuanta[i],combinedstateinfo.rightUnMapQuanta[i]);

      std::pair<int,int> pointedpair(i,quantastate);
      combinedquanta[combinepair] = pointedpair;
      std::pair<int,int> d_r_quanta (combinedstateinfo.rightUnMapQuanta[i],i);
      leftquantainfo leftquanta(combinedstateinfo.leftUnMapQuanta[i],stateinfo.quantaStates[combinedstateinfo.leftUnMapQuanta[i]] , quantastate);
      reversequanta[d_r_quanta] = leftquanta;

    }
    {
    std::string file;
    file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/quantalink-state" % new_site %"-"%currentstate %".tmp" );
    std::ofstream s(file.c_str());
    boost::archive::text_oarchive oa(s);
    oa << combinedquanta;
    s.close();
    }
    {
    std::string file;
    file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/reversequantalink-state" % new_site %"-"%currentstate %".tmp" );
    std::ofstream s(file.c_str());
    boost::archive::text_oarchive oa(s);
    oa << reversequanta;
    s.close();
    }
    //for (auto it =combinedquanta.begin(); it !=combinedquanta.end();it++)
    //  pout <<it->first.first <<" "<<it->first.second<<"-> " <<it->second.first<<" "<<it->second.second<<endl;





    std::vector<Matrix> lastRotationMatrix;
    for (int a=0; a<lastsitewave.nrows(); a++)
      for (int b=0; b<lastsitewave.ncols(); b++)
      {

        if (lastsitewave.allowed(a, b)) {
	        Matrix& lM = lastsitewave.operator_element(a, b);
	        Matrix tM = RotationMatrix[a];
          if (RotationMatrix[a].Ncols()==0) continue;
	        Matrix nM(1,1);
          nM(1,1) = 0.0;
	        MatrixMultiply(tM, 't', lM, 'n', nM, 1.0);
          lastRotationMatrix.push_back(nM);
          pout << nM<<endl;
        }
      }

      sites.push_back(new_site);
      if (dmrginp.spinAdapted())
        spinsites.push_back(new_site);
      else {
        spinsites.push_back(2*new_site);
        spinsites.push_back(2*new_site+1);
      }
      SaveRotationMatrix (spinsites, lastRotationMatrix, currentstate);



}

void StateQuantaInfo::readRotationandQuanta()
{
  quantalink.resize(0);
  reversequantalink.resize(0);
  RotationMatrices.resize(0);
  
  std::vector<Matrix> RotationMatrix;
  std::vector<Matrix> reducedRotationMatrix;//Remove the discarded quanta num.


  int niters;
  int new_site = 0;
  std::vector<int> spinsites;
  if (dmrginp.spinAdapted())
  {
  
    niters = dmrginp.last_site()-1;
    spinsites.push_back(new_site);
  }
  else {
    niters = dmrginp.last_site()/2-1;
    spinsites.push_back(2*new_site);
    spinsites.push_back(2*new_site+1);
  }


  for(int i=0; i<niters; i++)
  {
    new_site++;
    if (dmrginp.spinAdapted())
    {
      spinsites.push_back(new_site);
    }
    else {
      spinsites.push_back(2*new_site);
      spinsites.push_back(2*new_site+1);
    }
    
    reducedRotationMatrix.clear();
    LoadRotationMatrix (spinsites, RotationMatrix, currentstate);
    pout <<"total size: " <<RotationMatrix.size()<<endl;
    for (Matrix m : RotationMatrix)
    {
      if(m.Ncols()) reducedRotationMatrix.push_back(m);
    }
    RotationMatrices.push_back(reducedRotationMatrix);
  
    {
      std::string file;
      file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/quantalink-state" % new_site %"-"%currentstate %".tmp" );
      std::ifstream s(file.c_str());
      boost::archive::text_iarchive oa(s);
      std::map < intpair, intpair> combinedquanta;
      oa >> combinedquanta;
      s.close();
      quantalink.push_back(combinedquanta);
    }
    {
      std::string file;
      file = str(boost::format("%s%s%d%s%d%s") % dmrginp.load_prefix() % "/reversequantalink-state" % new_site %"-"%currentstate %".tmp" );
      std::ifstream s(file.c_str());
      boost::archive::text_iarchive oa(s);
      std::map < intpair, leftquantainfo> reversequanta;
      oa >> reversequanta;
      s.close();
      reversequantalink.push_back(reversequanta);
    }
  }
    norm = 0.0;
    //norm = DotProduct(lastsitewave,lastsitewave);
    for (Matrix m : RotationMatrices.back())
      norm += m(1,1)*m(1,1);

        
}

void StateQuantaInfo::get_slater()
{
  std::vector<int> ci;
  ci.assign(dmrginp.last_site(),0);
  double norm = 0.0;
  
  for(int n=0;n<pow(2,ci.size());n++)
  {
    for(int i=0;i<ci.size();i++)
      ci[i] = (n & ( 1 << i )) >> i;
    int particle_num = 0;
    for(int i=0;i<ci.size();i++)
      particle_num+= ci[i];
    if (particle_num!=dmrginp.total_particle_number()) continue;

    double coeff = getcoeff(ci);
    norm +=coeff*coeff;
    //if (fabs(coeff)>1e-4)
    {
    for(int i=0;i<ci.size();i++)
      pout << ci[i];
    pout <<": ";
    pout << coeff;
    pout <<endl;
    }
  }
  pout <<"END"<<endl;

}

void StateQuantaInfo::test_expectation()
{
  std::vector<int> ci;
  ci.assign(dmrginp.last_site(),0);
  double e0 = 0.0;
  
  for(int n=0;n<pow(2,ci.size());n++)
  {
    for(int i=0;i<ci.size();i++)
      ci[i] = (n & ( 1 << i )) >> i;
    int particle_num = 0;
    for(int i=0;i<ci.size();i++)
      particle_num+= ci[i];
    if (particle_num!=dmrginp.total_particle_number()) continue;

    double coeff = getcoeff(ci);
    e0 +=coeff*coeff*local_energy(ci);
  }
  pout <<"<H_0_E_0>" << e0<<endl;


}

double StateQuantaInfo::getcoeff(const std::vector<int>& ci)
{
    std::clock_t startTime = std::clock();
    int quantanum;
    int statenum=0;

    int niters;
    if (dmrginp.spinAdapted())
    {
    
      niters = dmrginp.last_site()-2;
    }
    else {
      niters = dmrginp.last_site()/2-2;
    }

    if (dmrginp.spinAdapted())
      quantanum = ci[0];
    else {
      quantanum = 2*ci[0]+ci[1];
    }

    //TODO
    //call DGEMV rather than DGEMM
    Matrix lbasis(1,1);
    lbasis(1,1) = 1.0;
    for(int i=0; i< niters+1;i++)
    {
      int localquantanum = dmrginp.spinAdapted()? ci[i+1]: 2*ci[2*i+2]+ci[2*i+3];
      auto iter = quantalink[i].find(std::make_pair(quantanum,localquantanum));
      if(iter==quantalink[i].end()) 
      {
        //pout << "Quanta: "<< quantanum<<" "<<localquantanum<<" was not found"<<endl;
        //pout << "This slater determinant avoids the wave function symmetry"<<endl;
        return 0.0;
      }
      auto pos = iter->second;
      //pout <<"Row position" << pos.first<<" " <<pos.second<<endl;

      Matrix rotationmatrix;
      rotationmatrix.ReSize(lbasis.Ncols(),RotationMatrices[i][pos.first].Ncols());
      //pout <<"total size: " <<RotationMatrices[i].size()<<endl;
      //pout <<"Size: " << rotationmatrix.Nrows() << " " << rotationmatrix.Ncols()<<endl;
      //pout <<"Size: " << RotationMatrices[i][pos.first].Nrows() << " " << RotationMatrices[i][pos.first].Ncols()<<endl;
      for(int k=1;k<=rotationmatrix.Nrows();k++)
      for(int l=1;l<=rotationmatrix.Ncols();l++)
      {
        rotationmatrix(k,l) = RotationMatrices[i][pos.first](pos.second+k,l);
      }

      Matrix newlbasis;
      newlbasis.ReSize(lbasis.Nrows(), rotationmatrix.Ncols());
      newlbasis= 0.0;
      MatrixMultiply(lbasis,'n', rotationmatrix,'n',newlbasis,1.0);
      //MatrixMultiply(lbasis,'n', rotationmatrix,'n',newlbasis,1.0, 0.0);

      //pout <<newlbasis<<endl;
      quantanum = pos.first;
      lbasis = newlbasis;
      //pout <<lbasis<<endl;

    }
    assert(lbasis.Ncols()==1);
    assert(lbasis.Nrows()==1);
    //pout <<lbasis<<endl;
    mpscoeffT += (std::clock() - startTime)/ (double) CLOCKS_PER_SEC;
    return lbasis(1,1);


}

double StateQuantaInfo::getcoeff(const bitstring& ci)
{
    std::clock_t startTime = std::clock();
    int quantanum;
    int statenum=0;

    int niters;
    if (dmrginp.spinAdapted())
    {
    
      niters = dmrginp.last_site()-2;
    }
    else {
      niters = dmrginp.last_site()/2-2;
    }

    if (dmrginp.spinAdapted())
      quantanum = ci[0];
    else {
      quantanum = 2*ci[0]+ci[1];
    }

    //TODO
    //call DGEMV rather than DGEMM
    Matrix lbasis(1,1);
    lbasis(1,1) = 1.0;
    for(int i=0; i< niters+1;i++)
    {
      int localquantanum = dmrginp.spinAdapted()? ci[i+1]: 2*ci[2*i+2]+ci[2*i+3];
      auto iter = quantalink[i].find(std::make_pair(quantanum,localquantanum));
      if(iter==quantalink[i].end()) 
      {
        //pout << "Quanta: "<< quantanum<<" "<<localquantanum<<" was not found"<<endl;
        //pout << "This slater determinant avoids the wave function symmetry"<<endl;
        return 0.0;
      }
      auto pos = iter->second;
      //pout <<"Row position" << pos.first<<" " <<pos.second<<endl;

      Matrix rotationmatrix;
      rotationmatrix.ReSize(lbasis.Ncols(),RotationMatrices[i][pos.first].Ncols());
      //pout <<"total size: " <<RotationMatrices[i].size()<<endl;
      //pout <<"Size: " << rotationmatrix.Nrows() << " " << rotationmatrix.Ncols()<<endl;
      //pout <<"Size: " << RotationMatrices[i][pos.first].Nrows() << " " << RotationMatrices[i][pos.first].Ncols()<<endl;
      for(int k=1;k<=rotationmatrix.Nrows();k++)
      for(int l=1;l<=rotationmatrix.Ncols();l++)
      {
        rotationmatrix(k,l) = RotationMatrices[i][pos.first](pos.second+k,l);
      }

      Matrix newlbasis;
      newlbasis.ReSize(lbasis.Nrows(), rotationmatrix.Ncols());
      newlbasis= 0.0;
      MatrixMultiply(lbasis,'n', rotationmatrix,'n',newlbasis,1.0);
      //MatrixMultiply(lbasis,'n', rotationmatrix,'n',newlbasis,1.0, 0.0);

      //pout <<newlbasis<<endl;
      quantanum = pos.first;
      lbasis = newlbasis;
      //pout <<lbasis<<endl;

    }
    assert(lbasis.Ncols()==1);
    assert(lbasis.Nrows()==1);
    //pout <<lbasis<<endl;
    mpscoeffT += (std::clock() - startTime)/ (double) CLOCKS_PER_SEC;
    return lbasis(1,1);


}

void StateQuantaInfo::initstate()
{
    bool direction=false;
    SweepParams sweepParams;
    sweepParams.current_root() = currentstate;
    Sweep::InitializeStateInfo(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, direction, currentstate);
}

void StateQuantaInfo::preparestate()
{
  if (mpigetrank() == 0) {
   // if (!dmrginp.stochasticpt_restart())
   // {
   //   initstate();
   //   SweepParams sweepParams;
   //   sweepParams.set_sweep_iter() = 0;
   //   sweepParams.current_root() = currentstate;
   //   StoreQuantaInformation(sweepParams, true);
   // }
    readRotationandQuanta();
  }

#ifndef SERIAL
      boost::mpi::communicator world;
      broadcast(world, norm, 0);
      broadcast(world, quantalink, 0);
      broadcast(world, reversequantalink, 0);
      broadcast(world, RotationMatrices, 0);
#endif
}

double StateQuantaInfo::sampling(std::vector<int>& ci)
{
    std::clock_t startTime = std::clock();
    int statenum=0;

    int niters;
    if (dmrginp.spinAdapted())
    {
    
      niters = dmrginp.last_site()-1;
    }
    else {
      niters = dmrginp.last_site()/2-1;
    }
    std::vector<int> reverseci;

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> f(0.0,1.0);



    ColumnVector oldwave;
    int rightquanta;
    for(int i=niters-1; i>=0 ;i--)
    {
      int localquanta =0;
      std::vector<double> coeff2(4,0.0);
      std::vector<ColumnVector> newwave(4);
      for(localquanta=0;localquanta<4;localquanta++)
      {
        int leftquanta;

        if(i==niters-1)
        {
          for(auto quantapair:reversequantalink[i])
          {
            if(quantapair.first.first == localquanta)
            {
              oldwave = RotationMatrices[i][quantapair.first.second];
              auto quanta_l = reversequantalink[i].at(quantapair.first);
              coeff2[localquanta] = oldwave(1)*oldwave(1);
              break;
            }
          }
        }
        else
        {
          auto quanta_d_r = std::make_pair(localquanta, rightquanta);
          auto iter = reversequantalink[i].find(quanta_d_r);
          if (iter==reversequantalink[i].end())
          {
            coeff2[localquanta] = 0.0;
          }
          else{

            auto quanta_l = iter->second;
            Matrix& rotateM1 = RotationMatrices[i][quanta_d_r.second];
            Matrix rotateM2;
            rotateM2.ReSize(quanta_l.second, rotateM1.Ncols());
            for(int k=1;k<=rotateM2.Nrows();k++)
            for(int l=1;l<=rotateM2.Ncols();l++)
              rotateM2(k,l) = rotateM1(k+quanta_l.third,l);
            newwave[localquanta].ReSize(rotateM2.Nrows());
            newwave[localquanta] = 0.0;
            MatrixMultiply(rotateM2,'n', oldwave,'n',newwave[localquanta], 1.0, 0.0);
            coeff2[localquanta] = dotproduct(newwave[localquanta],newwave[localquanta]);
          }
        }

      }
      std::vector<double> sum(4);
      std::partial_sum (coeff2.begin(),coeff2.end() , sum.begin());
      double x = f(generator);
      //cout <<"X: "<<x <<endl;
      //cout <<"coeff: "<<coeff2[0]<< " " <<coeff2[1]<< " "<<coeff2[2]<< " "<<coeff2[3]<< " "<<endl;
      if(x<sum[0]/sum[3])
      {
        localquanta = 0;
      }
      else if(x<sum[1]/sum[3])
      {
        localquanta = 1;
      }
      else if(x<sum[2]/sum[3])
      {
        localquanta = 2;
      }
      else
      {
        localquanta = 3;
      }
      if (i==niters-1)
      {
        for(auto quantapair:reversequantalink[i])
        {
          if(quantapair.first.first == localquanta)
          {
            rightquanta = quantapair.second.first;
            oldwave = RotationMatrices[i][quantapair.first.second];
            break;
          }
        }
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
      }
      else
      {
        auto quanta_d_r = std::make_pair(localquanta, rightquanta);
        auto quanta_l = reversequantalink[i].at(quanta_d_r);
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
        if (i==0)
        {
          reverseci.push_back(quanta_l.first%2);
          reverseci.push_back(quanta_l.first/2);
        }
        rightquanta = quanta_l.first;
        oldwave = newwave[localquanta];
      }

    }
    ci.resize(reverseci.size());
    for(int i=0;i<reverseci.size();i++)
      ci[i] = reverseci[reverseci.size()-1-i];
    assert(oldwave.Nrows()==1);
    mpssampleT += ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    return oldwave(1);



}

double StateQuantaInfo::sampling(bitstring& ci)
{
    std::clock_t startTime = std::clock();
    int statenum=0;

    int niters;
    if (dmrginp.spinAdapted())
    {
    
      niters = dmrginp.last_site()-1;
    }
    else {
      niters = dmrginp.last_site()/2-1;
    }
    std::vector<int> reverseci;

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> f(0.0,1.0);



    ColumnVector oldwave;
    int rightquanta;
    for(int i=niters-1; i>=0 ;i--)
    {
      int localquanta =0;
      std::vector<double> coeff2(4,0.0);
      std::vector<ColumnVector> newwave(4);
      for(localquanta=0;localquanta<4;localquanta++)
      {
        int leftquanta;

        if(i==niters-1)
        {
          for(auto quantapair:reversequantalink[i])
          {
            if(quantapair.first.first == localquanta)
            {
              oldwave = RotationMatrices[i][quantapair.first.second];
              auto quanta_l = reversequantalink[i].at(quantapair.first);
              coeff2[localquanta] = oldwave(1)*oldwave(1);
              break;
            }
          }
        }
        else
        {
          auto quanta_d_r = std::make_pair(localquanta, rightquanta);
          auto iter = reversequantalink[i].find(quanta_d_r);
          if (iter==reversequantalink[i].end())
          {
            coeff2[localquanta] = 0.0;
          }
          else{

            auto quanta_l = iter->second;
            Matrix& rotateM1 = RotationMatrices[i][quanta_d_r.second];
            Matrix rotateM2;
            rotateM2.ReSize(quanta_l.second, rotateM1.Ncols());
            for(int k=1;k<=rotateM2.Nrows();k++)
            for(int l=1;l<=rotateM2.Ncols();l++)
              rotateM2(k,l) = rotateM1(k+quanta_l.third,l);
            newwave[localquanta].ReSize(rotateM2.Nrows());
            newwave[localquanta] = 0.0;
            MatrixMultiply(rotateM2,'n', oldwave,'n',newwave[localquanta], 1.0, 0.0);
            coeff2[localquanta] = dotproduct(newwave[localquanta],newwave[localquanta]);
          }
        }

      }
      std::vector<double> sum(4);
      std::partial_sum (coeff2.begin(),coeff2.end() , sum.begin());
      double x = f(generator);
      //cout <<"X: "<<x <<endl;
      //cout <<"coeff: "<<coeff2[0]<< " " <<coeff2[1]<< " "<<coeff2[2]<< " "<<coeff2[3]<< " "<<endl;
      if(x<sum[0]/sum[3])
      {
        localquanta = 0;
      }
      else if(x<sum[1]/sum[3])
      {
        localquanta = 1;
      }
      else if(x<sum[2]/sum[3])
      {
        localquanta = 2;
      }
      else
      {
        localquanta = 3;
      }
      if (i==niters-1)
      {
        for(auto quantapair:reversequantalink[i])
        {
          if(quantapair.first.first == localquanta)
          {
            rightquanta = quantapair.second.first;
            oldwave = RotationMatrices[i][quantapair.first.second];
            break;
          }
        }
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
      }
      else
      {
        auto quanta_d_r = std::make_pair(localquanta, rightquanta);
        auto quanta_l = reversequantalink[i].at(quanta_d_r);
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
        if (i==0)
        {
          reverseci.push_back(quanta_l.first%2);
          reverseci.push_back(quanta_l.first/2);
        }
        rightquanta = quanta_l.first;
        oldwave = newwave[localquanta];
      }

    }
    ci.reset();
    for(int i=0;i<reverseci.size();i++)
      ci[i] = reverseci[reverseci.size()-1-i];
    assert(oldwave.Nrows()==1);
    mpssampleT += ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    return oldwave(1);
}

double StateQuantaInfo::sampling_approx(std::vector<int>& ci, double& prob)
{
    std::clock_t startTime = std::clock();
    int statenum=0;

    int niters;
    if (dmrginp.spinAdapted())
    {
      niters = dmrginp.last_site()-1;
    }
    else {
      niters = dmrginp.last_site()/2-1;
    }
    std::vector<int> reverseci;

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> f(0.0,1.0);

    std::vector<std::vector<double> > env(niters);
    ////Calculate env
    //
     // for (int i=0;i<niters;i++)
     // {
     //   if(i==0)
     //   {
     //   env[2].resize(4);
     //   for(auto quantapair:reversequantalink[i])
     //   {
     //     int right = quantapair.first.first;
     //     int dot = quantapair.first.second;
     //     int right = quantapair.second.first;
     //     for
     //     env[right] += fabs()
     //     if(quantapair.first.first == localquanta)
     //     {
     //       rightquanta = quantapair.second.first;
     //       oldwave = RotationMatrices[i][quantapair.first.second];
     //       break;
     //     }
     //   }
     //   reverseci.push_back(localquanta%2);
     //     

     //   }

     // }
    
    //



    prob =1.0;
    ColumnVector oldwave;
    int rightquanta;
    for(int i=niters-1; i>=0 ;i--)
    {
      int localquanta =0;
      std::vector<double> coeff2(4,0.0);
      std::vector<ColumnVector> newwave(4);
      for(localquanta=0;localquanta<4;localquanta++)
      {
        int leftquanta;

        if(i==niters-1)
        {
          for(auto quantapair:reversequantalink[i])
          {
            if(quantapair.first.first == localquanta)
            {
              oldwave = RotationMatrices[i][quantapair.first.second];
              auto quanta_l = reversequantalink[i].at(quantapair.first);
              coeff2[localquanta] = fabs(oldwave(1));
              //coeff2[localquanta] = oldwave(1)*oldwave(1);
              break;
            }
          }
        }
        else
        {
          auto quanta_d_r = std::make_pair(localquanta, rightquanta);
          auto iter = reversequantalink[i].find(quanta_d_r);
          if (iter==reversequantalink[i].end())
          {
            coeff2[localquanta] = 0.0;
          }
          else{

            auto quanta_l = iter->second;
            Matrix& rotateM1 = RotationMatrices[i][quanta_d_r.second];
            Matrix rotateM2;
            rotateM2.ReSize(quanta_l.second, rotateM1.Ncols());
            for(int k=1;k<=rotateM2.Nrows();k++)
            for(int l=1;l<=rotateM2.Ncols();l++)
              rotateM2(k,l) = rotateM1(k+quanta_l.third,l);
            newwave[localquanta].ReSize(rotateM2.Nrows());
            newwave[localquanta] = 0.0;
            MatrixMultiply(rotateM2,'n', oldwave,'n',newwave[localquanta], 1.0, 0.0);
            coeff2[localquanta] = 0;
            for (int k=1; k<=rotateM2.Nrows();k++)
            coeff2[localquanta] += fabs(newwave[localquanta](k));
            //coeff2[localquanta] = dotproduct(newwave[localquanta],newwave[localquanta]);
          }
        }

      }
      std::vector<double> sum(4);
      std::partial_sum (coeff2.begin(),coeff2.end() , sum.begin());
      double x = f(generator);
      //cout <<"X: "<<x <<endl;
      //cout <<"coeff: "<<coeff2[0]<< " " <<coeff2[1]<< " "<<coeff2[2]<< " "<<coeff2[3]<< " "<<endl;
      if(x<sum[0]/sum[3])
      {
        localquanta = 0;
        prob *= coeff2[0]/sum[3];
      }
      else if(x<sum[1]/sum[3])
      {
        localquanta = 1;
        prob *= coeff2[1]/sum[3];
      }
      else if(x<sum[2]/sum[3])
      {
        localquanta = 2;
        prob *= coeff2[2]/sum[3];
      }
      else
      {
        localquanta = 3;
        prob *= coeff2[3]/sum[3];
      }
      if (i==niters-1)
      {
        for(auto quantapair:reversequantalink[i])
        {
          if(quantapair.first.first == localquanta)
          {
            rightquanta = quantapair.second.first;
            oldwave = RotationMatrices[i][quantapair.first.second];
            break;
          }
        }
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
      }
      else
      {
        auto quanta_d_r = std::make_pair(localquanta, rightquanta);
        auto quanta_l = reversequantalink[i].at(quanta_d_r);
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
        if (i==0)
        {
          reverseci.push_back(quanta_l.first%2);
          reverseci.push_back(quanta_l.first/2);
        }
        rightquanta = quanta_l.first;
        oldwave = newwave[localquanta];
      }

    }
    ci.resize(reverseci.size());
    for(int i=0;i<reverseci.size();i++)
      ci[i] = reverseci[reverseci.size()-1-i];
    assert(oldwave.Nrows()==1);
    mpssampleT += ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    return oldwave(1);



}

double StateQuantaInfo::sampling_approx(bitstring& ci, double& prob)
{
    std::clock_t startTime = std::clock();
    int statenum=0;

    int niters;
    if (dmrginp.spinAdapted())
    {
      niters = dmrginp.last_site()-1;
    }
    else {
      niters = dmrginp.last_site()/2-1;
    }
    std::vector<int> reverseci;

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> f(0.0,1.0);

    std::vector<std::vector<double> > env(niters);
    ////Calculate env
    //
     // for (int i=0;i<niters;i++)
     // {
     //   if(i==0)
     //   {
     //   env[2].resize(4);
     //   for(auto quantapair:reversequantalink[i])
     //   {
     //     int right = quantapair.first.first;
     //     int dot = quantapair.first.second;
     //     int right = quantapair.second.first;
     //     for
     //     env[right] += fabs()
     //     if(quantapair.first.first == localquanta)
     //     {
     //       rightquanta = quantapair.second.first;
     //       oldwave = RotationMatrices[i][quantapair.first.second];
     //       break;
     //     }
     //   }
     //   reverseci.push_back(localquanta%2);
     //     

     //   }

     // }
    
    //



    prob =1.0;
    ColumnVector oldwave;
    int rightquanta;
    for(int i=niters-1; i>=0 ;i--)
    {
      int localquanta =0;
      std::vector<double> coeff2(4,0.0);
      std::vector<ColumnVector> newwave(4);
      for(localquanta=0;localquanta<4;localquanta++)
      {
        int leftquanta;

        if(i==niters-1)
        {
          for(auto quantapair:reversequantalink[i])
          {
            if(quantapair.first.first == localquanta)
            {
              oldwave = RotationMatrices[i][quantapair.first.second];
              auto quanta_l = reversequantalink[i].at(quantapair.first);
              coeff2[localquanta] = fabs(oldwave(1));
              //coeff2[localquanta] = oldwave(1)*oldwave(1);
              break;
            }
          }
        }
        else
        {
          auto quanta_d_r = std::make_pair(localquanta, rightquanta);
          auto iter = reversequantalink[i].find(quanta_d_r);
          if (iter==reversequantalink[i].end())
          {
            coeff2[localquanta] = 0.0;
          }
          else{

            auto quanta_l = iter->second;
            Matrix& rotateM1 = RotationMatrices[i][quanta_d_r.second];
            Matrix rotateM2;
            rotateM2.ReSize(quanta_l.second, rotateM1.Ncols());
            for(int k=1;k<=rotateM2.Nrows();k++)
            for(int l=1;l<=rotateM2.Ncols();l++)
              rotateM2(k,l) = rotateM1(k+quanta_l.third,l);
            newwave[localquanta].ReSize(rotateM2.Nrows());
            newwave[localquanta] = 0.0;
            MatrixMultiply(rotateM2,'n', oldwave,'n',newwave[localquanta], 1.0, 0.0);
            coeff2[localquanta] = 0;
            for (int k=1; k<=rotateM2.Nrows();k++)
            coeff2[localquanta] += fabs(newwave[localquanta](k));
            //coeff2[localquanta] = dotproduct(newwave[localquanta],newwave[localquanta]);
          }
        }

      }
      std::vector<double> sum(4);
      std::partial_sum (coeff2.begin(),coeff2.end() , sum.begin());
      double x = f(generator);
      //cout <<"X: "<<x <<endl;
      //cout <<"coeff: "<<coeff2[0]<< " " <<coeff2[1]<< " "<<coeff2[2]<< " "<<coeff2[3]<< " "<<endl;
      if(x<sum[0]/sum[3])
      {
        localquanta = 0;
        prob *= coeff2[0]/sum[3];
      }
      else if(x<sum[1]/sum[3])
      {
        localquanta = 1;
        prob *= coeff2[1]/sum[3];
      }
      else if(x<sum[2]/sum[3])
      {
        localquanta = 2;
        prob *= coeff2[2]/sum[3];
      }
      else
      {
        localquanta = 3;
        prob *= coeff2[3]/sum[3];
      }
      if (i==niters-1)
      {
        for(auto quantapair:reversequantalink[i])
        {
          if(quantapair.first.first == localquanta)
          {
            rightquanta = quantapair.second.first;
            oldwave = RotationMatrices[i][quantapair.first.second];
            break;
          }
        }
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
      }
      else
      {
        auto quanta_d_r = std::make_pair(localquanta, rightquanta);
        auto quanta_l = reversequantalink[i].at(quanta_d_r);
        reverseci.push_back(localquanta%2);
        reverseci.push_back(localquanta/2);
        if (i==0)
        {
          reverseci.push_back(quanta_l.first%2);
          reverseci.push_back(quanta_l.first/2);
        }
        rightquanta = quanta_l.first;
        oldwave = newwave[localquanta];
      }

    }
    ci.reset();
    for(int i=0;i<reverseci.size();i++)
      ci[i] = reverseci[reverseci.size()-1-i];
    assert(oldwave.Nrows()==1);
    mpssampleT += ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    return oldwave(1);



}

//double StateQuantaInfo::local_energy(const std::vector<int>& ci, int integralIndex)
//{
//  std::clock_t startTime = std::clock();
//  double energy = 0;
//  energy += coreEnergy[integralIndex];
//  for(int i=0;i<ci.size();i++)
//    energy += ci[i]*v_1[integralIndex](i,i);
//  for(int i=0;i<ci.size();i++)
//  for(int j=i+1;j<ci.size();j++)
//    energy += ci[i]*ci[j]*v_2[integralIndex](i,j,i,j);
//  for(int i=0;i<ci.size();i++)
//  for(int j=i+1;j<ci.size();j++)
//    energy -= ci[i]*ci[j]*v_2[integralIndex](j,i,i,j);
//  //for(int i=0;i<ci.size()/2;i++)
//  //  energy -= ci[2*i]*ci[2*i+1]*v_2[integralIndex](2*i,2*i+1,2*i,2*i+1);
//  HdiagonalT += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
//  return energy;
//
//}

double StateQuantaInfo::local_energy(const std::vector<int>& ci, int integralIndex)
{
  std::clock_t startTime = std::clock();

  double energy = coreEnergy[integralIndex];
  std::vector<int> occupy;
  for(int i=0;i<ci.size();i++)
    if(ci[i])
      occupy.push_back(i);
  for(auto i: occupy)
    energy += v_1[integralIndex](i, i);

  for(int i=0;i<occupy.size();i++)
  for(int j=i+1;j<occupy.size();j++)
  {
    energy += v_2[integralIndex](occupy[i],occupy[j],occupy[i],occupy[j]);
    energy -= v_2[integralIndex](occupy[j],occupy[i],occupy[i],occupy[j]);
  }
  HdiagonalT += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
  return energy;

}

double StateQuantaInfo::local_energy(const longbitarray& ci, int integralIndex)
{
  std::clock_t startTime = std::clock();

  double energy = coreEnergy[integralIndex];
  int n = dmrginp.last_site();

  std::vector<int> occupy;
  for(int i=0;i<n;i++)
    if(ci.getocc(i))
      occupy.push_back(i);
  for(auto i: occupy)
    energy += v_1[integralIndex](i, i);

  for(int i=0;i<occupy.size();i++)
  for(int j=i+1;j<occupy.size();j++)
  {
    energy += v_2[integralIndex](occupy[i],occupy[j],occupy[i],occupy[j]);
    energy -= v_2[integralIndex](occupy[j],occupy[i],occupy[i],occupy[j]);
  }
  HdiagonalT += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
  return energy;

}

double StateQuantaInfo::local_energy(const bitstring& ci, int integralIndex)
{
  std::clock_t startTime = std::clock();

  double energy = coreEnergy[integralIndex];
  int n = dmrginp.last_site();

  std::vector<int> occupy;
  for(int i=0;i<n;i++)
    if(ci[i])
      occupy.push_back(i);
  for(auto i: occupy)
    energy += v_1[integralIndex](i, i);

  for(int i=0;i<occupy.size();i++)
  for(int j=i+1;j<occupy.size();j++)
  {
    energy += v_2[integralIndex](occupy[i],occupy[j],occupy[i],occupy[j]);
    energy -= v_2[integralIndex](occupy[j],occupy[i],occupy[i],occupy[j]);
  }
  HdiagonalT += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
  return energy;

}

void check_sampling_approx(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  std::vector<int> ci;
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();
  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  double PT_onesample = 0.0;
  double PT_onesample_2 = 0.0;
  double PT_twosample = 0.0;
  double PT_twosample_2 = 0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  std::unordered_map<bitstring, double> det_energy;
  int n_spinorb = dmrginp.last_site();

  //double sum_c=0.0;
  //double sum_c_2=0.0;
  //for (int j=0;j<n_H_samples*num_sample;j++)
  //{

  //  std::vector<int> determinant;
  //  double prob;
  //  //double coeff = baseState.sampling(determinant);
  //  double coeff = baseState.sampling_approx(determinant, prob);
  //  double temp = fabs(coeff)/prob;
  //  sum_c += temp/(n_H_samples*num_sample);
  //  sum_c_2 += temp*temp/(n_H_samples*num_sample);
  //  //norm += 1/(coeff*coeff);
  //}
  //{
  //boost::mpi::communicator world;
  //sum_c = all_reduce(world, sum_c, std::plus<double>())/world.size();
  //sum_c_2 = all_reduce(world, sum_c_2, std::plus<double>())/world.size();
  //pout <<"sum of |c|: "<<sum_c<<"+"<< sqrt((sum_c_2*2-sum_c*sum_c)/(num_sample*n_H_samples*world.size()))<<endl;
  //}

  double minele = 1000.0;
  double maxele = 0.0;
  double timer1=0.0;
  double timer2=0.0;
  double timer3=0.0;
  double tol = dmrginp.stochasticpt_tol();
  std::clock_t startTime = std::clock();
  std::clock_t begintime;
  //sd_hashtable.reserve(num_sample);
  for (int i=0;i<num_sample;i++)
  {
    begintime = std::clock();
    double pt2=0.0;
    double pt2_one=0.0;
    double pt2_two=0.0;
    double norm= 0.0;
    double H01_temp=0.0;
    for (int j=0;j<n_H_samples;j++)
    {

      bitstring determinant;
      double prob;
      //double coeff = baseState.sampling(determinant);
      double coeff = baseState.sampling_approx(determinant, prob);
      double local_energy0;
      if(det_energy.find(determinant)==det_energy.end())
      {
        local_energy0= baseState.local_energy(determinant, 0);
        det_energy[determinant] = local_energy0;
      }
      else{
        local_energy0= det_energy[determinant];

      }

      double temp =coeff*coeff/(prob*local_energy0);
      H00 +=temp/(num_sample*n_H_samples);
      H00_2 +=temp*temp/(num_sample*n_H_samples);
      sd_hashtable0[determinant] += 1.0;
      sd_hashtable0_p[determinant] = prob;
      //norm += 1/(coeff*coeff);
    }
    timer1 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = baseState.getcoeff(iter->first);
      double freq = iter->second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1, det_energy, fabs(tol/coeff));
      sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);
      //TODO
      //Remove zeroth order sd from first order wave function.
      //for(auto iter1= sd_hashtable0.begin();iter1!=sd_hashtable0.end();iter1++)
      //  sd_hashtable1.erase(iter1->first);

      double factor = baseState.local_energy(iter->first, 1)/baseState.local_energy(iter->first,0);
      H01_temp +=factor*freq*coeff*coeff/(sd_hashtable0_p[iter->first]*n_H_samples);
      H01_temp -=factor*freq*freq*coeff*coeff/(n_H_samples*(n_H_samples-1)*sd_hashtable0_p[iter->first]*sd_hashtable0_p[iter->first]);
      //H01_temp -=factor*freq/sqrt((n_H_samples*(n_H_samples-1)));

        //pt2 += (freq*(n_H_samples-1)- freq*freq/(coeff*coeff))/(n_H_samples*(n_H_samples-1));
      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        //double e = baseState.local_energy(iter1->first, 0);
        double e = det_energy[iter1->first];
        //TODO
        //if(sd_hashtable0.find(iter1->first)==sd_hashtable0.end())
        //pt2 += (freq*(n_H_samples-1)- freq*freq/(coeff*coeff))*v*v/(n_H_samples*(n_H_samples-1)*e);
        pt2_one +=freq*v*v*coeff*coeff/(sd_hashtable0_p[iter->first]*n_H_samples*e);
        pt2_two -= freq*freq*v*v*coeff*coeff/(sd_hashtable0_p[iter->first]*sd_hashtable0_p[iter->first]*n_H_samples*(n_H_samples-1)*e);
        sd_hashtable1sum[iter1->first] += v*freq*(coeff/sd_hashtable0_p[iter->first]);
        //pout <<"table0";
        //pout <<  std::bitset<64>(iter1->first.data[0])<<" "<<iter1->second<<endl;
        //for(auto newiter = sd_hashtable0.begin();newiter != sd_hashtable0.end(); newiter++)
        //  pout <<  std::bitset<64>(newiter->first.data[0])<<" "<<newiter->second<<endl;

      }
    //  sd_hashtable1.clear();
    //  pout <<"SIZE: " <<sd_hashtable1.size()<<endl;
    //}
      sd_hashtable1.clear();
    }

    timer2 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();
    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();iter++)
    {
      double v = iter->second;
      //double e = baseState.local_energy(iter->first, 0);
      double e = det_energy[iter->first];
      //TODO
      pt2_two += v*v/((n_H_samples-1)*(n_H_samples)*e);
      if(sd_hashtable0.find(iter->first)!=sd_hashtable0.end())
      {
        double coeff =baseState.getcoeff(iter->first);
        //pout << "CI: ";
        //for(int i:iter->first.to_vector(n_spinorb))
        //  pout <<i<<" ";
        //pout<<endl;
        //pout <<"Coeff: " <<coeff<<endl;
        //pout <<"e: " <<e<<endl;
 
        H01_temp  += v*sd_hashtable0[iter->first]*coeff/(sd_hashtable0_p[iter->first]*n_H_samples*(n_H_samples-1)*e);
      }
    }
    pt2 = pt2_one + pt2_two;
    pout <<"PT2: "<< pt2<<endl;
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    PT_onesample +=pt2_one/num_sample;
    PT_onesample_2 +=pt2_one*pt2_one/num_sample;
    PT_twosample +=pt2_two/num_sample;
    PT_twosample_2 +=pt2_two*pt2_two/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable0_p.clear();
    //pout <<"unordered_map elements: "<<sd_hashtable1sum.size()<<endl;
    //pout <<"unordered_map bucket size: "<<sd_hashtable1sum.bucket_count()<<endl;
    //int maxcount = 0;
    //for(int i=0;i<sd_hashtable1sum.bucket_count();i++)
    //{
    //  pout <<sd_hashtable1sum.bucket_size(i)<<" ";
    //  int n = sd_hashtable1sum.bucket_size(i);
    //  maxcount = max(maxcount, n);
    //}
    //pout <<endl;
    //pout <<"unordered_map max bucket: "<<maxcount<<endl;
    sd_hashtable1sum.clear();
    //pout <<"unordered_map bucket size: "<<sd_hashtable1sum.bucket_count()<<endl;

    det_energy.clear();
    timer3 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
  }

  pout <<"Approx sampling |c| time: " <<(std::clock() - startTime)/ (double) CLOCKS_PER_SEC<<endl;


  {
  boost::mpi::communicator world;
  PT_onesample = all_reduce(world, PT_onesample, std::plus<double>())/world.size();
  PT_onesample_2 = all_reduce(world, PT_onesample_2, std::plus<double>())/world.size();
  PT_twosample = all_reduce(world, PT_twosample, std::plus<double>())/world.size();
  PT_twosample_2 = all_reduce(world, PT_twosample_2, std::plus<double>())/world.size();
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"PT one sample" << PT_onesample<<" + "<<sqrt((PT_onesample_2-PT_onesample*PT_onesample)/(num_sample*world.size()))<<endl;
  pout <<"PT two sample" << PT_twosample<<" + "<<sqrt((PT_twosample_2-PT_twosample*PT_twosample)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*n_H_samples*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  if(mpigetrank() ==0)
  {
    cout <<"Timer for head node"<<endl;
    cout <<"MPS sampling: " << baseState.mpssampleT<<endl;
    cout <<"MPS coeff: " << baseState.mpscoeffT<<endl;
    cout <<"Diagonal term: " << baseState.HdiagonalT<<endl;
    cout <<"Single Excitation: " << baseheatbath.singleT<<endl;
    cout <<"Double Excitation: " << baseheatbath.doubleT<<endl;
    cout <<"Excitation: " << baseheatbath.excitation_T<<endl;
    cout <<"Energy: " << baseheatbath.energy_T<<endl;
    cout <<"factor: " << baseheatbath.factor_T<<endl;
    cout <<"Phase1: " << timer1<<endl;
    cout <<"Phase2: " << timer2<<endl;
    cout <<"Phase3: " << timer3<<endl;
  }

  return;


}

void compressedMPS(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();
  StateQuantaInfo perturberState(1000);
  perturberState.preparestate();

  simplemps zeromps, sampledmps;
  zeromps.build(0);
  sampledmps.build(1000);


  //get_slater();
  std::vector<int> ci;
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();
  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  double PT_onesample = 0.0;
  double PT_onesample_2 = 0.0;
  double PT_twosample = 0.0;
  double PT_twosample_2 = 0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  std::unordered_map<bitstring, double> det_energy;
  std::unordered_map<bitstring, double> det_coeff;
  int n_spinorb = dmrginp.last_site();

  double minele = 1000.0;
  double maxele = 0.0;
  double timer1=0.0;
  double timer2=0.0;
  double timer3=0.0;
  double tol = dmrginp.stochasticpt_tol();
  std::clock_t startTime = std::clock();
  std::clock_t begintime;
  //sd_hashtable.reserve(num_sample);
  {
    begintime = std::clock();
    double pt2=0.0;
    double pt2_one=0.0;
    double pt2_two=0.0;
    double norm= 0.0;
    double H01_temp=0.0;
    pout <<"Number of samples: "<< mpigetsize() <<" X "<< num_sample<<endl;
    for (int i=0;i<num_sample;i++)
    {
      bitstring determinant;
      double coeff = zeromps.sampling(determinant);
      //double coeff = baseState.sampling(determinant);
      double local_energy0;
      if(det_energy.find(determinant)==det_energy.end())
      {
        local_energy0= baseState.local_energy(determinant, 0);
        det_energy[determinant] = local_energy0;
      }
      else{
        local_energy0= det_energy[determinant];

      }
      H00 += 1/(local_energy0*num_sample);
      H00_2 += 1/(local_energy0*local_energy0*num_sample);
      //norm += 1/(coeff*coeff);
    }
    timer1 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();
    long num_coeff=0;
    for (int i=0;i<num_sample;i++)
    {

      bitstring determinant;
      //double prob;
      double coeff = sampledmps.sampling(determinant);
      double e = baseState.local_energy(determinant, 0);

      baseheatbath.allexcite(determinant, 1.0,sd_hashtable1, fabs(tol*100));
      //baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1, fabs(tol/coeff));
      sd_hashtable1[determinant] += baseState.local_energy(determinant, 1);
      num_coeff += sd_hashtable1.size();

      double overlap=0.0;
      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        overlap +=  zeromps.getcoeff(iter1->first)*iter1->second;
        //testmps.getcoeff(iter1->first);
      }
      double temp =perturberState.norm*overlap*overlap/(e*coeff*coeff);
      H11 += temp/num_sample;
      H11_2 += temp*temp/num_sample;
      temp = perturberState.norm/e;
      test1 += temp/num_sample;
      test1_2 += temp*temp/num_sample;

      temp = perturberState.norm*overlap*zeromps.getcoeff(determinant)/(e*coeff*coeff);
      H01 += temp/num_sample;
      H01_2 += temp*temp/num_sample;


      sd_hashtable1.clear();

      if(i%(num_sample/10)==0)
      {
        pout <<"Finished " <<i <<" samples"<<endl;
      }

      //sd_hashtable0[determinant] += 1.0;
    }
    timer2 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();

    pout <<"Num of coeff to compute: "<< num_coeff<<endl;
    timer3 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
  }

  pout <<"Approx sampling |c| time: " <<(std::clock() - startTime)/ (double) CLOCKS_PER_SEC<<endl;


  {
  boost::mpi::communicator world;
  PT_onesample = all_reduce(world, PT_onesample, std::plus<double>())/world.size();
  PT_onesample_2 = all_reduce(world, PT_onesample_2, std::plus<double>())/world.size();
  PT_twosample = all_reduce(world, PT_twosample, std::plus<double>())/world.size();
  PT_twosample_2 = all_reduce(world, PT_twosample_2, std::plus<double>())/world.size();
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  //pout <<"PT one sample" << PT_onesample<<" + "<<sqrt((PT_onesample_2-PT_onesample*PT_onesample)/(num_sample*world.size()))<<endl;
  //pout <<"PT two sample" << PT_twosample<<" + "<<sqrt((PT_twosample_2-PT_twosample*PT_twosample)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  if(mpigetrank() ==0)
  {
    cout <<"Timer for head node"<<endl;
    cout <<"MPS sampling: " << baseState.mpssampleT<<endl;
    cout <<"MPS coeff: " << baseState.mpscoeffT<<endl;
    cout <<"MPS coeff: " << zeromps.mpscoeffT<<endl;
    cout <<"Diagonal term: " << baseState.HdiagonalT<<endl;
    cout <<"Single Excitation: " << baseheatbath.singleT<<endl;
    cout <<"Double Excitation: " << baseheatbath.doubleT<<endl;
    cout <<"Excitation: " << baseheatbath.excitation_T<<endl;
    cout <<"Energy: " << baseheatbath.energy_T<<endl;
    cout <<"factor: " << baseheatbath.factor_T<<endl;
    cout <<"Phase1: " << timer1<<endl;
    cout <<"Phase2: " << timer2<<endl;
    cout <<"Phase3: " << timer3<<endl;
  }

  return;


}

void approxsampling_twoH(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();



  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  std::unordered_map<bitstring, double> sd_hashtable2;
  std::unordered_map<bitstring, double> sd_hashtable2sum;
  std::unordered_map<bitstring, double> det_energy;
  int n_spinorb = dmrginp.last_site();
  for (int i=0;i<num_sample;i++)
    for (int j=0;j<n_H_samples;j++)
    {
      bitstring determinant;
      double coeff = baseState.sampling(determinant);
      double local_energy0= baseState.local_energy(determinant, 0);
      test1 +=1.0/(local_energy0*num_sample*n_H_samples);
      test1_2 +=1.0/(local_energy0*local_energy0*num_sample*n_H_samples);
    }
  {
  boost::mpi::communicator world;
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*n_H_samples*world.size()))<<endl;
  }

  double minele = 1000.0;
  double maxele = 0.0;
  double tol = dmrginp.stochasticpt_tol();
  double lowest_weight = 1.0;
  //sd_hashtable.reserve(num_sample);
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> uniformdistribution(0.0,1.0);
  for (int i=0;i<num_sample;i++)
  {
    double pt2_one=0.0;
    double pt2_two=0.0;
    double pt2=0.0;
    double norm= 0.0;
    double H01_temp=0.0;
    double H00_temp=0.0;
    for (int j=0;j<n_H_samples;j++)
    {
      bitstring determinant;
      double prob;
      //double coeff = baseState.sampling(determinant);
      double coeff = baseState.sampling_approx(determinant, prob);
      double local_energy0;
      if(det_energy.find(determinant)==det_energy.end())
      {
        local_energy0= baseState.local_energy(determinant, 0);
        det_energy[determinant] = local_energy0;
      }
      else{
        local_energy0= det_energy[determinant];

      }

      double temp =coeff*coeff/(prob*local_energy0);
      H00_temp +=temp/(n_H_samples);
      sd_hashtable0[determinant] += 1.0;
      sd_hashtable0_p[determinant] = prob;
    }
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = baseState.getcoeff(iter->first);
      double freq = iter->second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1, fabs(tol/coeff));
      sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);

      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        if(fabs(v*freq)<lowest_weight)
        {

          if (uniformdistribution(generator) < abs(v/lowest_weight)) 
          {

            int sign = signbit(v)? -1 : 1;
            sd_hashtable1sum[iter1->first] += lowest_weight*sign*(coeff/sd_hashtable0_p[iter->first]);
          }
        }
        else
            sd_hashtable1sum[iter1->first] += v*freq*(coeff/sd_hashtable0_p[iter->first]);


        //sd_hashtable1sum[iter1->first] += v*freq*(coeff/sd_hashtable0_p[iter->first]);
      }
      sd_hashtable1.clear();
    }
    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();)
    {
      if(fabs(iter->second)<lowest_weight) 
      {
        if (uniformdistribution(generator) < abs(iter->second)/lowest_weight) 
        {
          int sign = signbit(iter->second)? -1 : 1;
          sd_hashtable1sum[iter->first] = lowest_weight*sign;
          iter++;
        }
        else
          iter = sd_hashtable1sum.erase(iter);
      }
      else
        iter++;
    }
    pout <<"Size of sd_hashtable1sum: "<< sd_hashtable1sum.size()<<endl;
    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();iter++)
    {
      double v = iter->second;
      double e = baseState.local_energy(iter->first, 0);
      //double e = det_energy[iter->first];
      baseheatbath.allexcite(iter->first, v/e,sd_hashtable2, tol);
      sd_hashtable2[iter->first] += baseState.local_energy(iter->first, 1)*v/e;
      for(auto iter1=sd_hashtable2.begin();iter1!=sd_hashtable2.end();iter1++)
      {
        double v = iter1->second;
        //if(fabs(v)<lowest_weight)
        //{
        //  if (uniformdistribution(generator) < abs(v/lowest_weight)) 
        //  {
        //    int sign = signbit(v)? -1 : 1;
        //    sd_hashtable2sum[iter1->first] += lowest_weight*sign;
        //  }
        //}
        //else
          sd_hashtable2sum[iter1->first] += iter1->second;
      }
      sd_hashtable2.clear();
    }
    for(auto iter=sd_hashtable2sum.begin();iter!=sd_hashtable2sum.end();)
    {
      if(fabs(iter->second)<lowest_weight) 
      {
        if (uniformdistribution(generator) < abs(iter->second)/lowest_weight) 
        {
          int sign = signbit(iter->second)? -1 : 1;
          sd_hashtable2sum[iter->first] = lowest_weight*sign;
          iter++;
        }
        else
          iter = sd_hashtable2sum.erase(iter);
      }
      else
        iter++;
    }
    pout <<"Size of sd_hashtable2sum: "<< sd_hashtable2sum.size()<<endl;
    for(auto iter=sd_hashtable2sum.begin();iter!=sd_hashtable2sum.end();iter++)
    {
      if(fabs(iter->second)<tol*100) continue;
      double coeff = baseState.getcoeff(iter->first);
      pt2 +=iter->second*coeff/n_H_samples;
    }

    pout <<"PT2 (total): "<<pt2<<endl;
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    H00 += H00_temp/num_sample;
    H00_2 += H00_temp*H00_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable1sum.clear();
    sd_hashtable2sum.clear();
  }

  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*n_H_samples*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  return;


}

void check_sampling_approx_combined(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  std::vector<int> ci;
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();
  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  std::unordered_map<bitstring, double> det_energy;
  int n_spinorb = dmrginp.last_site();

  //double sum_c=0.0;
  //double sum_c_2=0.0;
  //for (int j=0;j<n_H_samples*num_sample;j++)
  //{

  //  std::vector<int> determinant;
  //  double prob;
  //  //double coeff = baseState.sampling(determinant);
  //  double coeff = baseState.sampling_approx(determinant, prob);
  //  double temp = fabs(coeff)/prob;
  //  sum_c += temp/(n_H_samples*num_sample);
  //  sum_c_2 += temp*temp/(n_H_samples*num_sample);
  //  //norm += 1/(coeff*coeff);
  //}
  //{
  //boost::mpi::communicator world;
  //sum_c = all_reduce(world, sum_c, std::plus<double>())/world.size();
  //sum_c_2 = all_reduce(world, sum_c_2, std::plus<double>())/world.size();
  //pout <<"sum of |c|: "<<sum_c<<"+"<< sqrt((sum_c_2*2-sum_c*sum_c)/(num_sample*n_H_samples*world.size()))<<endl;
  //}

  double minele = 1000.0;
  double maxele = 0.0;
  double timer1=0.0;
  double timer2=0.0;
  double timer3=0.0;
  double tol = dmrginp.stochasticpt_tol();
  std::clock_t startTime = std::clock();
  std::clock_t begintime;
  //sd_hashtable.reserve(num_sample);
  for (int i=0;i<num_sample;i++)
    for (int j=0;j<n_H_samples;j++)
    {
      bitstring determinant;
      double coeff = baseState.sampling(determinant);
      double local_energy0= baseState.local_energy(determinant, 0);
      test1 +=1.0/(local_energy0*num_sample*n_H_samples);
      test1_2 +=1.0/(local_energy0*local_energy0*num_sample*n_H_samples);
    }
  {
  boost::mpi::communicator world;
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*n_H_samples*world.size()))<<endl;
  }
  for (int i=0;i<num_sample;i++)
  {
    begintime = std::clock();
    double pt2=0.0;
    double norm= 0.0;
    double H01_temp=0.0;
    double H00_temp=0.0;
    for (int j=0;j<n_H_samples;j++)
    {

      bitstring determinant;
      double prob;
      //double coeff = baseState.sampling(determinant);
      double coeff = baseState.sampling_approx(determinant, prob);
      double local_energy0;
      if(det_energy.find(determinant)==det_energy.end())
      {
        local_energy0= baseState.local_energy(determinant, 0);
        det_energy[determinant] = local_energy0;
      }
      else{
        local_energy0= det_energy[determinant];

      }

      double temp =coeff*coeff/(prob*local_energy0);
      H00_temp +=temp/(n_H_samples);
      sd_hashtable0[determinant] += 1.0;
      sd_hashtable0_p[determinant] = prob;
      //norm += 1/(coeff*coeff);
    }
    timer1 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = baseState.getcoeff(iter->first);
      double freq = iter->second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1, det_energy, fabs(tol/coeff));
      sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);
      //TODO
      //Remove zeroth order sd from first order wave function.
      //for(auto iter1= sd_hashtable0.begin();iter1!=sd_hashtable0.end();iter1++)
      //  sd_hashtable1.erase(iter1->first);

      double factor = baseState.local_energy(iter->first, 1)/baseState.local_energy(iter->first,0);
      //H01_temp +=factor*freq*coeff*coeff/(sd_hashtable0_p[iter->first]*n_H_samples);
      //H01_temp -=factor*freq*freq*coeff*coeff/(n_H_samples*(n_H_samples-1)*sd_hashtable0_p[iter->first]*sd_hashtable0_p[iter->first]);
      //H01_temp -=factor*freq/sqrt((n_H_samples*(n_H_samples-1)));

        //pt2 += (freq*(n_H_samples-1)- freq*freq/(coeff*coeff))/(n_H_samples*(n_H_samples-1));
      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        //double e = baseState.local_energy(iter1->first, 0);
        double e = det_energy[iter1->first];
        //TODO
        //if(sd_hashtable0.find(iter1->first)==sd_hashtable0.end())
        //pt2 += (freq*(n_H_samples-1)- freq*freq/(coeff*coeff))*v*v/(n_H_samples*(n_H_samples-1)*e);
        pt2 +=freq*v*v*coeff*coeff/(sd_hashtable0_p[iter->first]*n_H_samples*e);
        pt2 -= freq*freq*v*v*coeff*coeff/(sd_hashtable0_p[iter->first]*sd_hashtable0_p[iter->first]*n_H_samples*(n_H_samples-1)*e);
        sd_hashtable1sum[iter1->first] += v*freq*(coeff/sd_hashtable0_p[iter->first]);
        //pout <<"table0";
        //pout <<  std::bitset<64>(iter1->first.data[0])<<" "<<iter1->second<<endl;
        //for(auto newiter = sd_hashtable0.begin();newiter != sd_hashtable0.end(); newiter++)
        //  pout <<  std::bitset<64>(newiter->first.data[0])<<" "<<newiter->second<<endl;

      }
    //  sd_hashtable1.clear();
    //  pout <<"SIZE: " <<sd_hashtable1.size()<<endl;
    //}
      sd_hashtable1.clear();
    }

    timer2 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();
    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();iter++)
    {
      double v = iter->second;
      //double e = baseState.local_energy(iter->first, 0);
      double e = det_energy[iter->first];
      //TODO
      pt2 += v*v/((n_H_samples-1)*(n_H_samples)*e);
      //if(sd_hashtable0.find(iter->first)!=sd_hashtable0.end())
      {
        double coeff =baseState.getcoeff(iter->first);
        //pout << "CI: ";
        //for(int i:iter->first.to_vector(n_spinorb))
        //  pout <<i<<" ";
        //pout<<endl;
        //pout <<"Coeff: " <<coeff<<endl;
        //pout <<"e: " <<e<<endl;
 
        H01_temp  += v*coeff/(n_H_samples*e);
      }
    }
    double pt_energy = -pt2 + H01_temp*H01_temp/test1;
    //double pt_energy = -pt2 + H01_temp*H01_temp/(H00_temp);
    test0 += pt_energy/num_sample;
    test0_2 += pt_energy*pt_energy/num_sample;
    pout <<"PT current"<< pt_energy <<endl;
    H00 += H00_temp/num_sample;
    H00_2 += H00_temp*H00_temp/num_sample;
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable0_p.clear();
    //pout <<"unordered_map elements: "<<sd_hashtable1sum.size()<<endl;
    //pout <<"unordered_map bucket size: "<<sd_hashtable1sum.bucket_count()<<endl;
    //int maxcount = 0;
    //for(int i=0;i<sd_hashtable1sum.bucket_count();i++)
    //{
    //  pout <<sd_hashtable1sum.bucket_size(i)<<" ";
    //  int n = sd_hashtable1sum.bucket_size(i);
    //  maxcount = max(maxcount, n);
    //}
    //pout <<endl;
    //pout <<"unordered_map max bucket: "<<maxcount<<endl;
    sd_hashtable1sum.clear();
    //pout <<"unordered_map bucket size: "<<sd_hashtable1sum.bucket_count()<<endl;

    det_energy.clear();
    timer3 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
  }

  pout <<"Approx sampling |c| time: " <<(std::clock() - startTime)/ (double) CLOCKS_PER_SEC<<endl;


  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  //test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  //test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  if(mpigetrank() ==0)
  {
    cout <<"Timer for head node"<<endl;
    cout <<"MPS sampling: " << baseState.mpssampleT<<endl;
    cout <<"MPS coeff: " << baseState.mpscoeffT<<endl;
    cout <<"Diagonal term: " << baseState.HdiagonalT<<endl;
    cout <<"Single Excitation: " << baseheatbath.singleT<<endl;
    cout <<"Double Excitation: " << baseheatbath.doubleT<<endl;
    cout <<"Excitation: " << baseheatbath.excitation_T<<endl;
    cout <<"Energy: " << baseheatbath.energy_T<<endl;
    cout <<"factor: " << baseheatbath.factor_T<<endl;
    cout <<"Phase1: " << timer1<<endl;
    cout <<"Phase2: " << timer2<<endl;
    cout <<"Phase3: " << timer3<<endl;
  }

  return;


}

void exactpt()
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();
  heatbath baseheatbath;
  baseheatbath.precompute();

  double H11=0.0;
  double H00=0.0;
  double H01=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test1=0.0;
  double test2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;


  double pt2 = 0.0;
  double x=0.0;
  std::vector<int> cilist;
  std::vector<double> cicoeff;
  double norm = 0.0;
  for(long n=0;n<pow(2,dmrginp.last_site());n++)
  {
    bitstring ci(n);
    if (ci.count()!=dmrginp.total_particle_number()) continue;

    double coeff = baseState.getcoeff(ci);
    if(fabs(coeff)>1e-13)
    {

      sd_hashtable0[ci] = coeff;
      baseheatbath.allexcite(ci, coeff,sd_hashtable1,0.0);
      sd_hashtable1[ci] += coeff*baseState.local_energy(ci, 1);
      H00 += coeff*coeff/baseState.local_energy(ci, 0);
    }
  }
  double min_e = 1000.0;
  double E0=0;
  for(auto iter=sd_hashtable1.begin();iter!=sd_hashtable1.end();iter++)
  {
      double v = iter->second;
      double e = baseState.local_energy(iter->first, 0);
      min_e = min(e, min_e);
      H11 += v*v/e;
      if(sd_hashtable0.find(iter->first)!=sd_hashtable0.end())
      {
        H01 += sd_hashtable0[iter->first]*v/e;
        E0 +=sd_hashtable0[iter->first]*v;
      }
      //else
      //{
      //  pout <<"default value: "<<sd_hashtable0[iter->first]<<endl;
      //  pout <<sd_hashtable0[iter->first]*v/e<<endl;
      //}
  }
  pout <<"E0: "<<E0<<endl;
  pout <<"Min(e): " <<min_e<<endl;
  pout <<"<0|VQ(H0-E0)^{-1}QV|0>: "<<H11<<endl;
  pout <<"<0|(H0-E0)^{-1}|0>: "<<H00<<endl;
  pout <<"<0|(H0-E0)^{-1}QV|0>: "<<H01<<endl;
  double pt2_energy = -1.0*H11+H01*H01/H00;
  pout <<"PT2 energy: " << pt2_energy<<endl;
  return;
}

void printoutall(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();


  double x=0.0;
  std::vector<int> cilist;
  std::vector<double> cicoeff;
  double norm = 0.0;
  for(long n=0;n<pow(2,dmrginp.last_site());n++)
  {
    bitstring ci(n);
    if (ci.count()!=dmrginp.total_particle_number()) continue;

    double coeff = baseState.getcoeff(ci);
    coeff = fabs(coeff);
    //coeff = coeff*coeff;
    if(fabs(coeff)>1e-13)
    {
      double local_energy0= baseState.local_energy(ci, 0);
      cilist.push_back(n);
      cicoeff.push_back(coeff);
      x += coeff*coeff/local_energy0;
    }
  }
  pout <<"Exact: "<<x<<endl;
  pout <<"Size of sd: "<< cicoeff.size()<<endl;
  double sum_x=std::accumulate(cicoeff.begin(), cicoeff.end(), 0.0);
  pout <<"sum of |c|: " <<sum_x<<endl;
  //std::sort(cicoeff.begin(),cicoeff.end(),std::greater<double>());
  //cout <<"Prob"<<endl;
  //for(auto i:cicoeff)
  //  cout <<i<<endl;


  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  int n_spinorb = dmrginp.last_site();

  double minele = 1000.0;
  double maxele = 0.0;
  double tol = dmrginp.stochasticpt_tol();
  //sd_hashtable.reserve(num_sample);
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::discrete_distribution<int> distribution(cicoeff.begin(),cicoeff.end());
  for (int i=0;i<num_sample;i++)
  {


    double pt2_one=0.0;
    double pt2_two=0.0;
    double norm= 0.0;
    for (int j=0;j<n_H_samples;j++)
    {

      int cin = distribution(generator);
      distribution.reset();
      bitstring ci(cilist[cin]);

      double local_energy0= baseState.local_energy(ci, 0);
      double coeff = baseState.getcoeff(ci);

      //double temp =coeff*coeff/(fabs(coeff)*local_energy0);
      double temp =1.0*fabs(coeff)*sum_x/local_energy0;
      H00 +=temp/(num_sample*n_H_samples);
      H00_2 +=temp*temp/(num_sample*n_H_samples);
      sd_hashtable0[ci] += 1.0;
    }
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    double H01_temp=0.0;
    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = baseState.getcoeff(iter->first);
      double freq = iter->second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1,fabs(tol/coeff));
      sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);
      //TODO
      //Remove zeroth order sd from first order wave function.
      //for(auto iter1= sd_hashtable0.begin();iter1!=sd_hashtable0.end();iter1++)
      //  sd_hashtable1.erase(iter1->first);
      //pout <<"sd_hashtable1.size: " << sd_hashtable1.size()<<endl;

      double factor = baseState.local_energy(iter->first, 1)/baseState.local_energy(iter->first,0);
      H01_temp +=factor*freq*fabs(coeff)*sum_x/n_H_samples;
      H01_temp -=factor*freq*freq*sum_x*sum_x/(n_H_samples*(n_H_samples-1));
      //H01_temp -=factor*freq/sqrt((n_H_samples*(n_H_samples-1)));

        //pt2 += (freq*(n_H_samples-1)- freq*freq/(coeff*coeff))/(n_H_samples*(n_H_samples-1));
      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        double e = baseState.local_energy(iter1->first, 0);
        pt2_one += freq*v*v*fabs(coeff)*sum_x/(n_H_samples*e);
        pt2_two -= freq*freq*sum_x*sum_x*v*v/(n_H_samples*(n_H_samples-1)*e);
        sd_hashtable1sum[iter1->first] += v*freq*sum_x*(coeff/fabs(coeff));
      }
      sd_hashtable1.clear();
    }
    //double temp = pt2;
    pout <<"PT2 (first part): "<<pt2_one<<endl;
    pout <<"Size of sd_hashtable1: "<< sd_hashtable1sum.size()<<endl;

    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();iter++)
    {
      double v = iter->second;
      double e = baseState.local_energy(iter->first, 0);
      //TODO
      //if(sd_hashtable0.find(iter->first)==sd_hashtable0.end())
      //pt2 += v*v/(sqrt((n_H_samples-1)*(n_H_samples))*e);
      //pt2 += v*v/((n_H_samples-1)*(n_H_samples)*e);
      pt2_two += v*v/((n_H_samples-1)*n_H_samples*e);
      if(sd_hashtable0.find(iter->first)!=sd_hashtable0.end())
      {
        double coeff =baseState.getcoeff(iter->first);
        //pout << "CI: ";
        //for(int i:iter->first.to_vector(n_spinorb))
        //  pout <<i<<" ";
        //pout<<endl;
        //pout <<"Coeff: " <<coeff<<endl;
        //pout <<"e: " <<e<<endl;
 
        H01_temp  += v*sd_hashtable0[iter->first]*sum_x*(coeff/fabs(coeff))/(n_H_samples*(n_H_samples-1)*e);
      }
    }
    pout <<"PT2 (second part): "<<pt2_two<<endl;
    double pt2 = pt2_one +pt2_two;
    pout <<"PT2 (total): "<<pt2<<endl;
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable1sum.clear();
  }

  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*n_H_samples*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  return;


}

void printoutall_twobatches(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();


  double x=0.0;
  std::vector<int> cilist;
  std::vector<double> cicoeff;
  double norm = 0.0;
  for(long n=0;n<pow(2,dmrginp.last_site());n++)
  {
    bitstring ci(n);
    if (ci.count()!=dmrginp.total_particle_number()) continue;

    double coeff = baseState.getcoeff(ci);
    coeff = fabs(coeff);
    //coeff = coeff*coeff;
    if(fabs(coeff)>1e-13)
    {
      double local_energy0= baseState.local_energy(ci, 0);
      cilist.push_back(n);
      cicoeff.push_back(coeff);
      x += coeff*coeff/local_energy0;
    }
  }
  pout <<"<0|(H0-E)|0>: "<<x<<endl;
  pout <<"Size of sd: "<< cicoeff.size()<<endl;
  double sum_x=std::accumulate(cicoeff.begin(), cicoeff.end(), 0.0);
  pout <<"sum of |c|: " <<sum_x<<endl;
  //std::sort(cicoeff.begin(),cicoeff.end(),std::greater<double>());
  //cout <<"Prob"<<endl;
  //for(auto i:cicoeff)
  //  cout <<i<<endl;


  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  int n_spinorb = dmrginp.last_site();

  double minele = 1000.0;
  double maxele = 0.0;
  double tol = dmrginp.stochasticpt_tol();
  //sd_hashtable.reserve(num_sample);
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::discrete_distribution<int> distribution(cicoeff.begin(),cicoeff.end());
  for (int i=0;i<num_sample;i++)
  {
    double pt2_one=0.0;
    double pt2_two=0.0;
    double pt2=0.0;
    double norm= 0.0;
    for (int j=0;j<n_H_samples;j++)
    {

      int cin = distribution(generator);
      bitstring ci(cilist[cin]);

      double local_energy0= baseState.local_energy(ci, 0);
      double coeff = baseState.getcoeff(ci);

      //double temp =coeff*coeff/(fabs(coeff)*local_energy0);
      double temp =1.0*fabs(coeff)*sum_x/local_energy0;
      H00 +=temp/(num_sample*n_H_samples);
      H00_2 +=temp*temp/(num_sample*n_H_samples);
      sd_hashtable0[ci] += 1.0;
    }
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    double H01_temp=0.0;
    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = baseState.getcoeff(iter->first);
      double freq = iter->second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1, fabs(tol/coeff));
      sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);

      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        sd_hashtable1sum[iter1->first] += v*freq*sum_x*(coeff/fabs(coeff));
      }
      sd_hashtable1.clear();
    }
    distribution.reset();
    for (int j=0;j<n_H_samples;j++)
    {

      int cin = distribution(generator);
      bitstring ci(cilist[cin]);

      double local_energy0= baseState.local_energy(ci, 0);
      double coeff = baseState.getcoeff(ci);
      baseheatbath.allexcite(ci,1.0, sd_hashtable1,fabs(tol/coeff));
      sd_hashtable1[ci] += baseState.local_energy(ci, 1);
      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        if(sd_hashtable1sum.find(iter1->first)!=sd_hashtable1sum.end())
        {
          double e = baseState.local_energy(iter1->first, 0);
          pt2 += sd_hashtable1sum[iter1->first]*coeff*sum_x*iter1->second/(n_H_samples*n_H_samples*e*fabs(coeff));
        }

      }
      sd_hashtable1.clear();

    }



    pout <<"PT2 (total): "<<pt2<<endl;
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable1sum.clear();
  }

  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*n_H_samples*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  return;


}

/*

void printoutall_twoH(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();


  double x=0.0;
  std::vector<int> cilist;
  std::vector<double> cicoeff;
  double norm = 0.0;
  for(long n=0;n<pow(2,dmrginp.last_site());n++)
  {
    bitstring ci(n);
    if (ci.count()!=dmrginp.total_particle_number()) continue;

    double coeff = baseState.getcoeff(ci);
    coeff = fabs(coeff);
    //coeff = coeff*coeff;
    if(fabs(coeff)>1e-15)
    {
      double local_energy0= baseState.local_energy(ci, 0);
      cilist.push_back(n);
      cicoeff.push_back(coeff);
      x += coeff*coeff/local_energy0;
    }
  }
  pout <<"<0|(H0-E)|0>: "<<x<<endl;
  pout <<"Size of sd: "<< cicoeff.size()<<endl;
  double sum_x=std::accumulate(cicoeff.begin(), cicoeff.end(), 0.0);
  pout <<"sum of |c|: " <<sum_x<<endl;
  //std::sort(cicoeff.begin(),cicoeff.end(),std::greater<double>());
  //cout <<"Prob"<<endl;
  //for(auto i:cicoeff)
  //  cout <<i<<endl;


  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  std::unordered_map<bitstring, double> sd_hashtable2;
  std::unordered_map<bitstring, double> sd_hashtable2sum;
  int n_spinorb = dmrginp.last_site();

  double minele = 1000.0;
  double maxele = 0.0;
  double tol = dmrginp.stochasticpt_tol();
  //sd_hashtable.reserve(num_sample);
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::discrete_distribution<int> distribution(cicoeff.begin(),cicoeff.end());
  for (int i=0;i<num_sample;i++)
  {
    double pt2_one=0.0;
    double pt2_two=0.0;
    double pt2=0.0;
    double norm= 0.0;
    for (int j=0;j<n_H_samples;j++)
    {

      int cin = distribution(generator);
      bitstring ci(cilist[cin]);

      double local_energy0= baseState.local_energy(ci, 0);
      double coeff = baseState.getcoeff(ci);

      //double temp =coeff*coeff/(fabs(coeff)*local_energy0);
      double temp =1.0*fabs(coeff)*sum_x/local_energy0;
      H00 +=temp/(num_sample*n_H_samples);
      H00_2 +=temp*temp/(num_sample*n_H_samples);
      sd_hashtable0[ci] += 1.0;
    }
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    double H01_temp=0.0;
    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = baseState.getcoeff(iter->first);
      double freq = iter->second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1, fabs(tol/coeff));
      sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);

      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        //if(fabs(v)>=tol)
        sd_hashtable1sum[iter1->first] += v*freq*sum_x*(coeff/fabs(coeff))/n_H_samples;
      }
      sd_hashtable1.clear();
    }
    pout <<"Size of sd_hashtable1sum: "<< sd_hashtable1sum.size()<<endl;
    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();iter++)
    {
      double v = iter->second;
      if(fabs(v)<tol*10) continue;
      double e = baseState.local_energy(iter->first, 0);
      baseheatbath.allexcite(iter->first, v/e,sd_hashtable2, 1e-13);
      sd_hashtable2[iter->first] += baseState.local_energy(iter->first, 1)*v/e;
      for(auto iter1=sd_hashtable2.begin();iter1!=sd_hashtable2.end();iter1++)
      {
        sd_hashtable2sum[iter1->first] += iter1->second;
      }
      sd_hashtable2.clear();
    }
    pout <<"Size of sd_hashtable2sum: "<< sd_hashtable2sum.size()<<endl;
    for(auto iter=sd_hashtable2sum.begin();iter!=sd_hashtable2sum.end();iter++)
    {
      if(fabs(iter->second)<tol*10) continue;
      double coeff = baseState.getcoeff(iter->first);
      pt2 +=iter->second*coeff;
    }

    pout <<"PT2 (total): "<<pt2<<endl;
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable1sum.clear();
    sd_hashtable2sum.clear();
  }

  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*n_H_samples*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  return;


}
*/

void printoutall_twoH(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();


  double x=0.0;
  std::vector<int> cilist;
  std::vector<double> cicoeff;
  double norm = 0.0;
  for(long n=0;n<pow(2,dmrginp.last_site());n++)
  {
    bitstring ci(n);
    if (ci.count()!=dmrginp.total_particle_number()) continue;

    double coeff = baseState.getcoeff(ci);
    coeff = fabs(coeff);
    //coeff = coeff*coeff;
    if(fabs(coeff)>1e-15)
    {
      double local_energy0= baseState.local_energy(ci, 0);
      cilist.push_back(n);
      cicoeff.push_back(coeff);
      x += coeff*coeff/local_energy0;
    }
  }
  pout <<"<0|(H0-E)|0>: "<<x<<endl;
  pout <<"Size of sd: "<< cicoeff.size()<<endl;
  double sum_x=std::accumulate(cicoeff.begin(), cicoeff.end(), 0.0);
  pout <<"sum of |c|: " <<sum_x<<endl;
  //std::sort(cicoeff.begin(),cicoeff.end(),std::greater<double>());
  //cout <<"Prob"<<endl;
  //for(auto i:cicoeff)
  //  cout <<i<<endl;


  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  std::unordered_map<bitstring, double> sd_hashtable2;
  std::unordered_map<bitstring, double> sd_hashtable2sum;
  int n_spinorb = dmrginp.last_site();

  double minele = 1000.0;
  double maxele = 0.0;
  double tol = dmrginp.stochasticpt_tol();
  double lowest_weight = 1.0;
  //sd_hashtable.reserve(num_sample);
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::discrete_distribution<int> distribution(cicoeff.begin(),cicoeff.end());
  std::uniform_real_distribution<double> uniformdistribution(0.0,1.0);
  for (int i=0;i<num_sample;i++)
  {
    double pt2_one=0.0;
    double pt2_two=0.0;
    double pt2=0.0;
    double norm= 0.0;
    for (int j=0;j<n_H_samples;j++)
    {

      int cin = distribution(generator);
      bitstring ci(cilist[cin]);

      double local_energy0= baseState.local_energy(ci, 0);
      double coeff = baseState.getcoeff(ci);

      //double temp =coeff*coeff/(fabs(coeff)*local_energy0);
      double temp =1.0*fabs(coeff)*sum_x/local_energy0;
      H00 +=temp/(num_sample*n_H_samples);
      H00_2 +=temp*temp/(num_sample*n_H_samples);
      sd_hashtable0[ci] += 1.0;
    }
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    double H01_temp=0.0;
    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = baseState.getcoeff(iter->first);
      double freq = iter->second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1, fabs(tol/coeff));
      sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);

      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        //if(fabs(v)>=tol)
        sd_hashtable1sum[iter1->first] += v*freq*sum_x*(coeff/fabs(coeff));
      }
      sd_hashtable1.clear();
    }
    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();)
    {
      if(fabs(iter->second)<lowest_weight) 
      {
        if (uniformdistribution(generator) < abs(iter->second)/lowest_weight) 
        {
          int sign = signbit(iter->second)? -1 : 1;
          sd_hashtable1sum[iter->first] = lowest_weight*sign;
          iter++;
        }
        else
          iter = sd_hashtable1sum.erase(iter);
      }
      else
        iter++;
    }
    pout <<"Size of sd_hashtable1sum: "<< sd_hashtable1sum.size()<<endl;
    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();iter++)
    {
      double v = iter->second;
      if(fabs(v)<tol*10) continue;
      double e = baseState.local_energy(iter->first, 0);
      baseheatbath.allexcite(iter->first, v/e,sd_hashtable2, 1e-13);
      sd_hashtable2[iter->first] += baseState.local_energy(iter->first, 1)*v/e;
      for(auto iter1=sd_hashtable2.begin();iter1!=sd_hashtable2.end();iter1++)
      {
        sd_hashtable2sum[iter1->first] += iter1->second;
      }
      sd_hashtable2.clear();
    }
    for(auto iter=sd_hashtable2sum.begin();iter!=sd_hashtable2sum.end();)
    {
      if(fabs(iter->second)<lowest_weight) 
      {
        if (uniformdistribution(generator) < abs(iter->second)/lowest_weight) 
        {
          int sign = signbit(iter->second)? -1 : 1;
          sd_hashtable2sum[iter->first] = lowest_weight*sign;
          iter++;
        }
        else
          iter = sd_hashtable2sum.erase(iter);
      }
      else
        iter++;
    }
    pout <<"Size of sd_hashtable2sum: "<< sd_hashtable2sum.size()<<endl;
    for(auto iter=sd_hashtable2sum.begin();iter!=sd_hashtable2sum.end();iter++)
    {
      if(fabs(iter->second)<tol*10) continue;
      double coeff = baseState.getcoeff(iter->first);
      pt2 +=iter->second*coeff/n_H_samples;
    }

    pout <<"PT2 (total): "<<pt2<<endl;
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable1sum.clear();
    sd_hashtable2sum.clear();
  }

  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*n_H_samples*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  return;


}

void printoutall_combined(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();


  double x=0.0;
  std::vector<int> cilist;
  std::vector<double> cicoeff;
  double norm = 0.0;
  for(long n=0;n<pow(2,dmrginp.last_site());n++)
  {
    bitstring ci(n);
    if (ci.count()!=dmrginp.total_particle_number()) continue;

    double coeff = baseState.getcoeff(ci);
    coeff = fabs(coeff);
    //coeff = coeff*coeff;
    if(fabs(coeff)>1e-13)
    {
      double local_energy0= baseState.local_energy(ci, 0);
      cilist.push_back(n);
      cicoeff.push_back(coeff);
      x += coeff*coeff/local_energy0;
    }
  }
  pout <<"Exact: "<<x<<endl;
  pout <<"Size of sd: "<< cicoeff.size()<<endl;
  double sum_x=std::accumulate(cicoeff.begin(), cicoeff.end(), 0.0);
  pout <<"sum of |c|: " <<sum_x<<endl;
  //std::sort(cicoeff.begin(),cicoeff.end(),std::greater<double>());
  //cout <<"Prob"<<endl;
  //for(auto i:cicoeff)
  //  cout <<i<<endl;


  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  int n_spinorb = dmrginp.last_site();

  double minele = 1000.0;
  double maxele = 0.0;
  double tol = dmrginp.stochasticpt_tol();
  tol = 0.0;
  //sd_hashtable.reserve(num_sample);
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::discrete_distribution<int> distribution(cicoeff.begin(),cicoeff.end());
  for (int i=0;i<num_sample;i++)
  {


    double H01_temp=0.0;
    double H00_temp=0.0;
    double pt2=0.0;
    double norm= 0.0;
    for (int j=0;j<n_H_samples;j++)
    {

      int cin = distribution(generator);
      bitstring ci(cilist[cin]);

      double local_energy0= baseState.local_energy(ci, 0);
      double coeff = baseState.getcoeff(ci);

      //double temp =coeff*coeff/(fabs(coeff)*local_energy0);
      double temp =1.0*fabs(coeff)*sum_x/local_energy0;
      H00_temp +=temp/(n_H_samples);
      sd_hashtable0[ci] += 1.0;
    }
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = baseState.getcoeff(iter->first);
      double freq = iter->second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1,fabs(tol/coeff));
      sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);
      //TODO
      //Remove zeroth order sd from first order wave function.
      //for(auto iter1= sd_hashtable0.begin();iter1!=sd_hashtable0.end();iter1++)
      //  sd_hashtable1.erase(iter1->first);
      //pout <<"sd_hashtable1.size: " << sd_hashtable1.size()<<endl;

      double factor = baseState.local_energy(iter->first, 1)/baseState.local_energy(iter->first,0);
      H01_temp +=factor*freq*fabs(coeff)*sum_x/n_H_samples;
      H01_temp -=factor*freq*freq*sum_x*sum_x/(n_H_samples*(n_H_samples-1));
      //H01_temp -=factor*freq/sqrt((n_H_samples*(n_H_samples-1)));

        //pt2 += (freq*(n_H_samples-1)- freq*freq/(coeff*coeff))/(n_H_samples*(n_H_samples-1));
      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        double e = baseState.local_energy(iter1->first, 0);
        pt2 += freq*v*v*fabs(coeff)*sum_x/(n_H_samples*e);
        pt2 -= freq*freq*sum_x*sum_x*v*v/(n_H_samples*(n_H_samples-1)*e);
        sd_hashtable1sum[iter1->first] += v*freq*sum_x*(coeff/fabs(coeff));
      }
      sd_hashtable1.clear();
    }
    pout <<"Size of sd_hashtable1: "<< sd_hashtable1sum.size()<<endl;

    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();iter++)
    {
      double v = iter->second;
      double e = baseState.local_energy(iter->first, 0);
      //TODO
      //if(sd_hashtable0.find(iter->first)==sd_hashtable0.end())
      //pt2 += v*v/(sqrt((n_H_samples-1)*(n_H_samples))*e);
      //pt2 += v*v/((n_H_samples-1)*(n_H_samples)*e);
      pt2 += v*v/((n_H_samples-1)*n_H_samples*e);
      if(sd_hashtable0.find(iter->first)!=sd_hashtable0.end())
      {
        double coeff =baseState.getcoeff(iter->first);
        //pout << "CI: ";
        //for(int i:iter->first.to_vector(n_spinorb))
        //  pout <<i<<" ";
        //pout<<endl;
        //pout <<"Coeff: " <<coeff<<endl;
        //pout <<"e: " <<e<<endl;
 
        H01_temp  += v*sd_hashtable0[iter->first]*sum_x*(coeff/fabs(coeff))/(n_H_samples*(n_H_samples-1)*e);
      }
    }
    double pt_energy = -pt2 + H01_temp*H01_temp/(H00_temp);
    test0 += pt_energy/num_sample;
    test0_2 += pt_energy*pt_energy/num_sample;
    pout <<"PT current"<< pt_energy <<endl;
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    H00 += H00_temp/num_sample;
    H00_2 += H00_temp*H00_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable1sum.clear();
  }

  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  return;


}

void splitspace(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();


  double x=0.0;
  std::vector<int> cilist;
  std::vector<double> cicoeff;
  double norm = 0.0;
  for(long n=0;n<pow(2,dmrginp.last_site());n++)
  {
    bitstring ci(n);
    if (ci.count()!=dmrginp.total_particle_number()) continue;

    double coeff = baseState.getcoeff(ci);
    coeff = fabs(coeff);
    //coeff = coeff*coeff;
    if(fabs(coeff)>1e-13)
    {
      double local_energy0= baseState.local_energy(ci, 0);
      cilist.push_back(n);
      cicoeff.push_back(coeff);
      x += coeff*coeff/local_energy0;
    }
  }
  pout <<"Exact: "<<x<<endl;
  pout <<"Size of sd: "<< cicoeff.size()<<endl;
  long varsize = cicoeff.size()/10;
  double sum_x=std::accumulate(cicoeff.begin(), cicoeff.begin()+varsize, 0.0);


  std::unordered_set<bitstring> vardets;
  for(auto iter= cilist.begin();iter!=cilist.begin()+varsize;iter++)
    vardets.insert(*iter);
  pout <<"sum of |c|: " <<sum_x<<endl;
  //std::sort(cicoeff.begin(),cicoeff.end(),std::greater<double>());
  //cout <<"Prob"<<endl;
  //for(auto i:cicoeff)
  //  cout <<i<<endl;


  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, double> sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable0_p;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  int n_spinorb = dmrginp.last_site();

  double minele = 1000.0;
  double maxele = 0.0;
  double tol = dmrginp.stochasticpt_tol();
  //sd_hashtable.reserve(num_sample);
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::discrete_distribution<int> distribution(cicoeff.begin(),cicoeff.begin()+varsize);
  for (int i=0;i<num_sample;i++)
  {


    double pt2=0.0;
    double norm= 0.0;
    for (int j=0;j<n_H_samples;j++)
    {

      int cin = distribution(generator);
      bitstring ci(cilist[cin]);

      double local_energy0= baseState.local_energy(ci, 0);
      double coeff = baseState.getcoeff(ci);

      //double temp =coeff*coeff/(fabs(coeff)*local_energy0);
      double temp =1.0*fabs(coeff)*sum_x/local_energy0;
      H00 +=temp/(num_sample*n_H_samples);
      H00_2 +=temp*temp/(num_sample*n_H_samples);
      sd_hashtable0[ci] += 1.0;
    }
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    double H01_temp=0.0;
    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = baseState.getcoeff(iter->first);
      double freq = iter->second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1,fabs(tol/coeff));
      //sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);
      //TODO
      //Remove zeroth order sd from first order wave function.
      //for(auto iter1= sd_hashtable0.begin();iter1!=sd_hashtable0.end();iter1++)
      //  sd_hashtable1.erase(iter1->first);
      for(auto iter1= cilist.begin();iter1!=cilist.begin()+varsize;iter1++)
        sd_hashtable1.erase(*iter1);
      //pout <<"sd_hashtable1.size: " << sd_hashtable1.size()<<endl;

      double factor = baseState.local_energy(iter->first, 1)/baseState.local_energy(iter->first,0);
      H01_temp +=factor*freq*fabs(coeff)*sum_x/n_H_samples;
      H01_temp -=factor*freq*freq*sum_x*sum_x/(n_H_samples*(n_H_samples-1));
      //H01_temp -=factor*freq/sqrt((n_H_samples*(n_H_samples-1)));

        //pt2 += (freq*(n_H_samples-1)- freq*freq/(coeff*coeff))/(n_H_samples*(n_H_samples-1));
      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        double e = baseState.local_energy(iter1->first, 0);
        pt2 += freq*v*v*fabs(coeff)*sum_x/(n_H_samples*e);
        pt2 -= freq*freq*sum_x*sum_x*v*v/(n_H_samples*(n_H_samples-1)*e);
        sd_hashtable1sum[iter1->first] += v*freq*sum_x*(coeff/fabs(coeff));
      }
      sd_hashtable1.clear();
    }
    pout <<"Size of sd_hashtable1: "<< sd_hashtable1sum.size()<<endl;

    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();iter++)
    {
      double v = iter->second;
      double e = baseState.local_energy(iter->first, 0);
      //TODO
      //if(sd_hashtable0.find(iter->first)==sd_hashtable0.end())
      //pt2 += v*v/(sqrt((n_H_samples-1)*(n_H_samples))*e);
      //pt2 += v*v/((n_H_samples-1)*(n_H_samples)*e);
      pt2 += v*v/((n_H_samples-1)*n_H_samples*e);
      if(sd_hashtable0.find(iter->first)!=sd_hashtable0.end())
      {
        double coeff =baseState.getcoeff(iter->first);
        //pout << "CI: ";
        //for(int i:iter->first.to_vector(n_spinorb))
        //  pout <<i<<" ";
        //pout<<endl;
        //pout <<"Coeff: " <<coeff<<endl;
        //pout <<"e: " <<e<<endl;
 
        H01_temp  += v*sd_hashtable0[iter->first]*sum_x*(coeff/fabs(coeff))/(n_H_samples*(n_H_samples-1)*e);
      }
    }
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable1sum.clear();
  }

  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*n_H_samples*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  return;


}

void check_heatbath(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  std::vector<int> ci;
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();
  double expectation0=0.0;
  double expectation1=0.0;
  double expectation2=0.0;
  double expectation0_2=0.0;
  double expectation1_2=0.0;
  double expectation2_2=0.0;
  for (int i=0;i<num_sample;i++)
  {
    double coeff = baseState.sampling(ci);
    double out_coeff;
    double tmp;
    std::vector<int> out_ci;
    //baseheatbath.doubleexcite(ci, 1.0/coeff, out_ci, out_coeff);
    baseheatbath.uniformdoubleexcite(ci, 1.0/coeff, out_ci, out_coeff);
    tmp=baseState.getcoeff(out_ci)*out_coeff;
    expectation2 += tmp;
    expectation2_2 += tmp*tmp;
    //baseheatbath.singleexcite(ci, 1.0/coeff, out_ci, out_coeff);
    baseheatbath.uniformsingleexcite(ci, 1.0/coeff, out_ci, out_coeff);
    //baseheatbath.exactsingleexcite(ci, 1.0/coeff, out_ci, out_coeff);
    tmp=baseState.getcoeff(out_ci)*out_coeff;
    expectation1 +=tmp;
    expectation1_2 +=tmp*tmp;

    tmp=baseState.local_energy(ci, 1);
    expectation0 +=tmp;
    expectation0_2 +=tmp*tmp;
    //expectation +=baseState.getcoeff(ci)*(1.0/coeff)*coreEnergy[0];
  }
  expectation0 /=num_sample;
  expectation1 /=num_sample;
  expectation2 /=num_sample;
  expectation0_2 /=num_sample;
  expectation1_2 /=num_sample;
  expectation2_2 /=num_sample;

  boost::mpi::communicator world;
  expectation0 = all_reduce(world, expectation0, std::plus<double>())/world.size();
  expectation1 = all_reduce(world, expectation1, std::plus<double>())/world.size();
  expectation2 = all_reduce(world, expectation2, std::plus<double>())/world.size();
  expectation0_2 = all_reduce(world, expectation0_2, std::plus<double>())/world.size();
  expectation1_2 = all_reduce(world, expectation1_2, std::plus<double>())/world.size();
  expectation2_2 = all_reduce(world, expectation2_2, std::plus<double>())/world.size();



  pout <<"EXPECT: " <<expectation0 <<" + "<< sqrt((expectation0_2-expectation0*expectation0)/(num_sample*world.size()))<<endl;
  pout <<"EXPECT: " <<expectation1 <<" + "<< sqrt((expectation1_2-expectation1*expectation1)/(num_sample*world.size()))<<endl;
  pout <<"EXPECT: " <<expectation2 <<" + "<< sqrt((expectation2_2-expectation2*expectation2)/(num_sample*world.size()))<<endl;
  pout <<"EXPECT: " <<expectation0+expectation1+expectation2 <<endl;

}

void check_overlap(const unsigned long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  std::vector<int> ci;
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();
  double expectation0=0.0;
  double reverse=0.0;
  double expectation1=0.0;
  double expectation2=0.0;
  for (int i=0;i<num_sample;i++)
  {
    double coeff = baseState.sampling(ci);
    double out_coeff;
    std::vector<int> out_ci;
    baseheatbath.uniformdoubleexcite(ci, 1.0/coeff, out_ci, out_coeff);
    //baseheatbath.doubleexcite(ci, 1.0/coeff, out_ci, out_coeff);
    expectation2 +=baseState.getcoeff(out_ci)*out_coeff/baseState.local_energy(out_ci);
    baseheatbath.uniformsingleexcite(ci, 1.0/coeff, out_ci, out_coeff);
    //baseheatbath.singleexcite(ci, 1.0/coeff, out_ci, out_coeff);
    expectation1 +=baseState.getcoeff(out_ci)*out_coeff/baseState.local_energy(out_ci);

    expectation0 +=baseState.local_energy(ci, 1)/baseState.local_energy(ci,0);
    //reverse +=1.0/baseState.local_energy(ci);
    //expectation +=baseState.getcoeff(ci)*(1.0/coeff)*coreEnergy[0];
  }
  expectation0 /=num_sample;
  expectation1 /=num_sample;
  expectation2 /=num_sample;
  reverse /=num_sample;
  //cout <<"REVERSE: " <<reverse <<endl;
  cout <<"EXPECT: " <<expectation0 <<endl;
  cout <<"EXPECT: " <<expectation1 <<endl;
  cout <<"EXPECT: " <<expectation2 <<endl;
  cout <<"EXPECT: " <<expectation0+expectation1+expectation2 <<endl;

}

void hci(long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  std::vector<int> ci;
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();
  double H11=0.0;
  double H11_2=0.0;
  double H11_new=0.0;
  double H11_2_new=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<bitstring, std::pair<double, unsigned int> > sd_hashtable0;
  std::unordered_map<bitstring, double> sd_hashtable1;
  std::unordered_map<bitstring, double> sd_hashtable1sum;
  int n_spinorb = dmrginp.last_site();

  double minele = 1000.0;
  double maxele = 0.0;
  double tol = dmrginp.stochasticpt_tol();
  //sd_hashtable.reserve(num_sample);
  std::clock_t startTime = std::clock();
  for (int i=0;i<num_sample;i++)
  {

    double pt2=0.0;
    double norm= 0.0;
    for (int j=0;j<n_H_samples;j++)
    {

      bitstring determinant;
      double coeff = baseState.sampling(determinant);
      double local_energy0= baseState.local_energy(determinant, 0);
      H00 +=1.0/(local_energy0*num_sample*n_H_samples);
      H00_2 +=1.0/(local_energy0*local_energy0*num_sample*n_H_samples);
      if(sd_hashtable0.find(determinant)==sd_hashtable0.end())
        sd_hashtable0[determinant] = make_pair(coeff, 1);
      else
        sd_hashtable0[determinant] = make_pair(coeff, sd_hashtable0[determinant].second+1);
    }
    pout <<"Size of sd_hashtable0: "<< sd_hashtable0.size()<<endl;

    double H01_temp=0.0;
    for(auto iter=sd_hashtable0.begin();iter!=sd_hashtable0.end();iter++)
    {
      double coeff = iter->second.first;
      double freq = iter->second.second;
      baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1,fabs(tol/coeff));
      sd_hashtable1[iter->first] += baseState.local_energy(iter->first, 1);

      double factor = baseState.local_energy(iter->first, 1)/baseState.local_energy(iter->first,0);
      H01_temp +=factor*freq/n_H_samples;
      H01_temp -=factor*freq*freq/(n_H_samples*(n_H_samples-1)*coeff*coeff);

      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        double v = iter1->second;
        double e = baseState.local_energy(iter1->first, 0);
        pt2 += freq*v*v/(n_H_samples*e);
        pt2 -= ( freq*freq/(coeff*coeff))*v*v/(n_H_samples*(n_H_samples-1)*e);
        sd_hashtable1sum[iter1->first] += v*freq/coeff;

      }
      sd_hashtable1.clear();
    }
    for(auto iter=sd_hashtable1sum.begin();iter!=sd_hashtable1sum.end();iter++)
    {
      double v = iter->second;
      double e = baseState.local_energy(iter->first, 0);
      pt2 += v*v/((n_H_samples-1)*n_H_samples*e);
      if(sd_hashtable0.find(iter->first)!=sd_hashtable0.end())
      {
        double coeff = sd_hashtable0[iter->first].first;
        double freq = sd_hashtable0[iter->first].second;
 
        H01_temp  += v*freq/(coeff*n_H_samples*(n_H_samples-1)*e);
      }
    }
    H11 +=pt2/num_sample;
    H11_2 +=pt2*pt2/num_sample;
    H01 += H01_temp/num_sample;
    H01_2 += H01_temp*H01_temp/num_sample;
    sd_hashtable0.clear();
    sd_hashtable1sum.clear();
  }

  pout <<"Sampling c^2 time: "<<(std::clock() - startTime)/ (double) CLOCKS_PER_SEC<<endl;

  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*n_H_samples*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  if(mpigetrank() ==0)
  {
    cout <<"Timer for head node"<<endl;
    cout <<"MPS sampling" << baseState.mpssampleT<<endl;
    cout <<"MPS coeff" << baseState.mpscoeffT<<endl;
    cout <<"Diagonal term" << baseState.HdiagonalT<<endl;
    cout <<"Single Excitation" << baseheatbath.singleT<<endl;
    cout <<"Double Excitation" << baseheatbath.doubleT<<endl;
  }

  return;

}

void heatbath_pt(long num_sample)
{
  StateQuantaInfo baseState(0);
  baseState.preparestate();


  //get_slater();
  std::vector<int> ci;
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.precompute();
  double H11=0.0;
  double H11_2=0.0;
  double H00=0.0;
  double H00_2=0.0;
  double H01=0.0;
  double H01_2=0.0;
  double expectation1=0.0;
  double test0=0.0;
  double test0_2=0.0;
  double test1=0.0;
  double test1_2=0.0;
  double test2=0.0;
  double test2_2=0.0;
  int n_H_samples = dmrginp.stochasticpt_Hsamples();
  std::unordered_map<longbitarray, double> sd_hashtable;
  //std::unordered_map<longbitarray, double> sd_hashtable;
  int n_spinorb = dmrginp.last_site();

  double minele = 1000.0;
  double maxele = 0.0;
  //sd_hashtable.reserve(num_sample);
  for (int i=0;i<num_sample;i++)
  {
    std::vector<int> determinant;
    double coeff = baseState.sampling(determinant);
    double local_energy0= baseState.local_energy(determinant, 0);
    H00 +=1.0/(local_energy0*num_sample);
    H00_2 +=1.0/(local_energy0*local_energy0*num_sample);

    double local_energy1= baseState.local_energy(determinant, 1);
    longbitarray cistring(determinant);
    insertorsum(sd_hashtable, cistring, local_energy1/(num_sample*coeff));

    for(int i=0;i<n_H_samples;i++)
    {
      std::vector<int> determinant1;
      double coeff1;
      baseheatbath.singleexcite(determinant, 1.0/(coeff*num_sample*n_H_samples), determinant1, coeff1);
      longbitarray cistring(determinant1);
      insertorsum(sd_hashtable, cistring, coeff1);
    }

    for(int i=0;i<n_H_samples;i++)
    {
      std::vector<int> determinant1;
      double coeff1;
      baseheatbath.doubleexcite(determinant, 1.0/(coeff*num_sample*n_H_samples), determinant1, coeff1);
      longbitarray cistring(determinant1);
      insertorsum(sd_hashtable, cistring, coeff1);
    }

  }

  pout <<"sd_hashtable size: " <<sd_hashtable.size()<<endl;
  for(auto iter = sd_hashtable.begin();iter!=sd_hashtable.end();iter++)
  {
    std::vector<int> determinant1 = iter->first.to_vector(n_spinorb);
    double local_energy0= baseState.local_energy(determinant1, 0);
    pout <<local_energy0<<endl;
    pout <<(iter->second*iter->second)/local_energy0<<endl;

    H11 += (iter->second*iter->second)/local_energy0;
  }
  pout <<"H11 "<<H11<<endl;
  H11_2 = H11 *H11;


  {
  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();
  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  }

  if(mpigetrank() ==0)
  {
    cout <<"Timer for head node"<<endl;
    cout <<"MPS sampling" << baseState.mpssampleT<<endl;
    cout <<"MPS coeff" << baseState.mpscoeffT<<endl;
    cout <<"Diagonal term" << baseState.HdiagonalT<<endl;
    cout <<"Single Excitation" << baseheatbath.singleT<<endl;
    cout <<"Double Excitation" << baseheatbath.doubleT<<endl;
  }

  return;


  //std::vector<double> coeff2;
  //std::vector< std::vector<int> > determinant2;
  for (int i=0;i<num_sample;i++)
  {
    std::vector< std::vector<int> > determinant1(3);
    std::vector<double> coeff1(3);
    double coeff = baseState.sampling(determinant1[0]);
    baseheatbath.doubleexcite(determinant1[0], 1.0/coeff, determinant1[2], coeff1[2]);
    baseheatbath.singleexcite(determinant1[0], 1.0/coeff, determinant1[1], coeff1[1]);
    coeff1[0] = baseState.local_energy(determinant1[0], 1)/coeff;
    double local_energy0 = baseState.local_energy(determinant1[0], 0);
    H00 +=1.0/local_energy0;
    H00_2 +=1.0/(local_energy0*local_energy0);
    
    //test0 += coeff1[0]*coeff;
    //test0_2 += coeff1[0]*coeff*coeff1[0]*coeff;
    //test1 += coeff1[1]*baseState.getcoeff(determinant1[1]);
    //test1_2 += pow(coeff1[1]*baseState.getcoeff(determinant1[1]),2);
    //test2 += coeff1[2]*baseState.getcoeff(determinant1[2]);
    //test2_2 += pow(coeff1[2]*baseState.getcoeff(determinant1[2]),2);
    for (int j=0;j<determinant1.size();j++)
    {
      std::vector<int> newsd;
      double newcoeff;

      double overlap= baseState.getcoeff(determinant1[j]);
      coeff1[j] /= baseState.local_energy(determinant1[j], 0);

      H01 += overlap*coeff1[j];
      H01_2 += (overlap*coeff1[j])*(overlap*coeff1[j]);
      newcoeff = coeff1[j]*overlap*baseState.local_energy(determinant1[j], 1);
      H11 += newcoeff;
      H11_2 += newcoeff*newcoeff;

      baseheatbath.singleexcite(determinant1[j], coeff1[j], newsd, newcoeff);
      newcoeff *= baseState.getcoeff(newsd);
      H11 += newcoeff;
      H11_2 += newcoeff*newcoeff;

      baseheatbath.doubleexcite(determinant1[j], coeff1[j], newsd, newcoeff);
      newcoeff *= baseState.getcoeff(newsd);
      H11 += newcoeff;
      H11_2 += newcoeff*newcoeff;
    }
  }
  
  H11 /=num_sample;
  H11_2 /=num_sample;
  H00 /=num_sample;
  H00_2 /=num_sample;
  H01 /=num_sample;
  H01_2 /=num_sample;
  test1/=num_sample;
  test1_2/=num_sample;
  test0/=num_sample;
  test0_2/=num_sample;
  test2/=num_sample;
  test2_2/=num_sample;

  boost::mpi::communicator world;
  H11 = all_reduce(world, H11, std::plus<double>())/world.size();
  H11_2 = all_reduce(world, H11_2, std::plus<double>())/world.size();
  H00 = all_reduce(world, H00, std::plus<double>())/world.size();
  H00_2 = all_reduce(world, H00_2, std::plus<double>())/world.size();
  H01 = all_reduce(world, H01, std::plus<double>())/world.size();
  H01_2 = all_reduce(world, H01_2, std::plus<double>())/world.size();
  test0 = all_reduce(world, test0, std::plus<double>())/world.size();
  test0_2 = all_reduce(world, test0_2, std::plus<double>())/world.size();
  test1 = all_reduce(world, test1, std::plus<double>())/world.size();
  test1_2 = all_reduce(world, test1_2, std::plus<double>())/world.size();
  test2 = all_reduce(world, test2, std::plus<double>())/world.size();
  test2_2 = all_reduce(world, test2_2, std::plus<double>())/world.size();


  if(mpigetrank() ==0)
  {
    cout <<"Timer for head node"<<endl;
    cout <<"MPS sampling" << baseState.mpssampleT<<endl;
    cout <<"MPS coeff" << baseState.mpscoeffT<<endl;
    cout <<"Diagonal term" << baseState.HdiagonalT<<endl;
    cout <<"Single Excitation" << baseheatbath.singleT<<endl;
    cout <<"Double Excitation" << baseheatbath.doubleT<<endl;
  }

  pout <<"TEST0: " <<test0<<" + "<<sqrt((test0_2-test0*test0)/(num_sample*world.size()))<<endl;
  pout <<"TEST1: " <<test1<<" + "<<sqrt((test1_2-test1*test1)/(num_sample*world.size()))<<endl;
  pout <<"TEST2: " <<test2<<" + "<<sqrt((test2_2-test2*test2)/(num_sample*world.size()))<<endl;
  pout <<"<0|VQ(H_0-E_0)^(-1)QV|0>" << H11<<" + "<<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)|0>" << H00<<" + "<<sqrt((H00_2-H00*H00)/(num_sample*world.size()))<<endl;
  pout <<"<0|(H_0-E_0)^(-1)QV|0>" << H01<<" + "<<sqrt((H01_2-H01*H01)/(num_sample*world.size()))<<endl;
  pout <<"PT energy: " << -H11+H01*H01/H00 << " + " <<sqrt((H11_2-H11*H11)/(num_sample*world.size()))<<endl;
}

int main(int argc, char* argv[])
{
  //test(argc,argv);
//  for(auto i: argv)
//    cout <<string(i)<<endl;
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  //MPI_Comm Calc;
  if (world.rank() == 0) {
    pout << "Runing with " << world.size() << " processors" << endl;
  }
#endif

//  for(auto i: argv)
//    cout <<string(i)<<endl;

  ReadInput(argv[1]);
  //dmrginp.matmultFlops.resize(1, 0.);
  dmrginp.initCumulTimer();
  bitstring::n_orb= dmrginp.last_site();

  dmrginp.initCumulTimer();
  long num_sample = dmrginp.stochasticpt_nsamples();
  //check_heatbath(num_sample);
  //check_sampling_approx(num_sample);
  //exactpt();
  //printoutall_combined(num_sample);
  //printoutall(num_sample);
  //printoutall_twobatches(num_sample);
  //printoutall_twoH(num_sample);
  //splitspace(num_sample);
  //hci(num_sample);
  //check_sampling_approx_combined(num_sample);
  //check_sampling_approx(num_sample);
  compressedMPS(num_sample);
  //approxsampling_twoH(num_sample);

  //check_overlap(num_sample);
  return 0;
  
  StateQuantaInfo baseState(0);
  baseState.preparestate();
  StateQuantaInfo perturberState(1000);
  perturberState.preparestate();
  //baseState.test_expectation();
  //perturberState.test_expectation();

  //get_slater();
  std::vector<int> ci;
  pout << "Num of samples: "<< world.size() <<" X " << num_sample<<endl;
  Timer sampletimer;
  //const unsigned long num_sample = 20000;
  double base_H0=0.0;
  double base_H01=0.0;
  double perturber_H0=0.0;
  double perturber_H01=0.0;
  double overlap0 = 0.0;
  double overlap1 = 0.0;
  double base_H0_2=0.0;
  double base_H01_2=0.0;
  double overlap0_2 = 0.0;
  double overlap1_2 = 0.0;
  double perturber_H0_2=0.0;
  double perturber_H01_2=0.0;
  double cutoverlap = 0.0;
  double testH0 = 0.0;
  double smallest_coeff=1.0;
  unsigned long realcount = 0;
  for (int i=0;i<num_sample;i++)
  {
    double coeff = baseState.sampling(ci);
    double tmp;
    tmp = 1/baseState.local_energy(ci);
    //if (fabs(StateQuantaInfo::local_energy(ci)) < 1e-13) cout <<"ERROR: " <<"Too small E_{ci}"<<endl;
    base_H0+= tmp/num_sample;
    base_H0_2+= tmp*tmp/num_sample;
    tmp = perturberState.getcoeff(ci)/coeff;
    overlap0 += tmp/num_sample;
    overlap0_2 += tmp*tmp/num_sample;
    tmp = perturberState.getcoeff(ci)/(baseState.local_energy(ci)*coeff);
    base_H01 += tmp/num_sample;
    base_H01_2 += tmp*tmp/num_sample;


    testH0 +=baseState.local_energy(ci)/num_sample;
    //testH0 +=StateQuantaInfo::local_energy(ci)/num_sample;

    
    smallest_coeff = min(smallest_coeff, coeff*coeff);
    //if(coeff*coeff > 1.e-6)
    //{
    //  realcount++;
    //  base_H01 += perturberState.getcoeff(ci)/(StateQuantaInfo::local_energy(ci)*num_sample*coeff);
    //  cutoverlap += perturberState.getcoeff(ci)/(coeff*num_sample);
    //
    //}
    //for(auto i:ci)
    //  cout <<i;
    //cout <<":" << coeff<<endl;
  }
  //pout <<"RealCount " <<realcount<<endl;
  //pout <<"testH0" <<testH0 <<endl;
  //pout <<"smallest_coeff" <<smallest_coeff<<endl;
  //pout <<"cutoverlap" <<cutoverlap<<endl;

  //cout <<world.rank()<<" "<<"base_H0: " << base_H0<<endl;

  testH0 = 0.0;
  cutoverlap = 0.0;

  realcount = 0;
  for (int i=0;i<num_sample;i++)
  {
    double tmp;
    double coeff = perturberState.sampling(ci);
    tmp = 1/baseState.local_energy(ci);
    perturber_H0+= tmp/num_sample;
    perturber_H0_2+= tmp*tmp/num_sample;
    //if(coeff*coeff/perturberState.norm > 1.e-6)
    //{
      tmp = baseState.getcoeff(ci)/(baseState.local_energy(ci)*coeff);
      perturber_H01 += tmp/num_sample;
      perturber_H01_2 += tmp*tmp/num_sample;
      realcount++;
      cutoverlap += baseState.getcoeff(ci)/(coeff*num_sample);
    //}
    tmp = baseState.getcoeff(ci)/(coeff);
    overlap1 += tmp/num_sample;
    overlap1_2 += tmp*tmp/num_sample;
    testH0 +=baseState.local_energy(ci)/num_sample;
  }
  //cout <<"RealCount " <<realcount<<endl;
  //cout <<"testH0" <<testH0 <<endl;
  //cout <<"cutoverlap" <<cutoverlap<<endl;


  /*
#ifndef SERIAL

  double global_base_H0;
  double global_base_H01;
  double global_perturber_H0;
  double global_perturber_H01;
  reduce(world, base_H0, global_base_H0, std::plus<double>(),0);
  reduce(world, base_H01, global_base_H01, std::plus<double>(),0);
  reduce(world, perturber_H0, global_perturber_H0, std::plus<double>(),0);
  reduce(world, perturber_H01, global_perturber_H01, std::plus<double>(),0);
  base_H0       =  global_base_H0;
  base_H01      =  global_base_H01;
  perturber_H0  =  global_perturber_H0;
  perturber_H01 =  global_perturber_H01;
#endif
  if(!world.rank())
  {
    base_H0 /= world.size();
    base_H01 /= world.size();
    perturber_H0 /= world.size();
    perturber_H01 /= world.size();
    perturber_H01 *=perturberState.norm;
    perturber_H0 *=perturberState.norm;
    cout <<"base_H0: " << base_H0<<endl;
    cout <<"base_H01: " << base_H01<<endl;
    cout <<"perturber_norm: " << perturberState.norm<<endl;
    cout <<"TEST: " <<test <<endl;
    cout <<"perturber_H0: " << perturber_H0<<endl;
    cout <<"perturber_H01: " << perturber_H01<<endl;

    double pt2_energy = -1.0*perturber_H0+perturber_H01*perturber_H01/(base_H0);
    cout <<"perturbation energy: " << pt2_energy<<endl;
    pt2_energy = -1.0*perturber_H0+base_H01*base_H01/(base_H0);
    cout <<"perturbation energy: " << pt2_energy<<endl;
  }

  */

  perturber_H01 *=perturberState.norm;
  perturber_H0  *=perturberState.norm;
#ifndef SERIAL
  std::vector<double> all_base_H0;
  std::vector<double> all_base_H01;
  std::vector<double> all_perturber_H0;
  std::vector<double> all_perturber_H01;
  std::vector<double> all_overlap0;
  std::vector<double> all_overlap1;
  std::vector<double> all_base_H0_2;
  std::vector<double> all_base_H01_2;
  std::vector<double> all_perturber_H0_2;
  std::vector<double> all_perturber_H01_2;
  std::vector<double> all_overlap0_2;
  std::vector<double> all_overlap1_2;
  gather(world, base_H0, all_base_H0, 0);
  gather(world, base_H01, all_base_H01,0);
  gather(world, perturber_H0, all_perturber_H0, 0);
  gather(world, perturber_H01, all_perturber_H01, 0);
  gather(world, overlap0, all_overlap0, 0);
  gather(world, overlap1, all_overlap1, 0);
  gather(world, base_H0_2, all_base_H0_2, 0);
  gather(world, base_H01_2, all_base_H01_2,0);
  gather(world, perturber_H0_2, all_perturber_H0_2, 0);
  gather(world, perturber_H01_2, all_perturber_H01_2, 0);
  gather(world, overlap0_2, all_overlap0_2, 0);
  gather(world, overlap1_2, all_overlap1_2, 0);
#else
  std::vector<double> all_base_H0(1,base_H0);
  std::vector<double> all_base_H01(1,base_H01);
  std::vector<double> all_perturber_H0(1,perturber_H0);
  std::vector<double> all_perturber_H01(1,perturber_H01);
  std::vector<double> all_overlap0(1,overlap0);
  std::vector<double> all_overlap1(1,overlap1);
#endif


  if(!world.rank())
  {

    pout <<"<0|(H_0-E_0)^{-1}|0>"<<endl;
    for(int i=0;i<world.size();i++)
      pout <<all_base_H0[i]<<"\t\t";
    pout <<endl<<"<1|(H_0-E_0)^{-1}|0>"<<endl;
    for(int i=0;i<world.size();i++)
      pout <<all_base_H01[i]<<"\t\t";
    pout <<endl<<"<0|(H_0-E_0)^{-1}|1>"<<endl;
    for(int i=0;i<world.size();i++)
      pout <<all_perturber_H01[i]<<"\t\t";
    pout <<endl<<"<1|(H_0-E_0)^{-1}|1>"<<endl;
    for(int i=0;i<world.size();i++)
      pout <<all_perturber_H0[i]<<"\t\t";



    std::vector<double> PTenergy(world.size());
    pout <<endl<<"PT2 energy in each process"<<endl;
    for(int i=0;i<world.size();i++)
    {
      PTenergy[i] = -1.0*all_perturber_H0[i]+all_perturber_H01[i]*all_perturber_H01[i]/(all_base_H0[0]);
      pout << PTenergy[i] <<"\t\t";
    }
    pout <<endl<<"Average: " << std::accumulate(PTenergy.begin(),PTenergy.end(),0.0)/world.size()<<endl;
    pout <<"PT2 energy in each process"<<endl;
    for(int i=0;i<world.size();i++)
    {
      PTenergy[i] = -1.0*all_perturber_H0[i]+all_base_H01[i]*all_base_H01[i]/(all_base_H0[0]);
      pout << PTenergy[i] <<"\t\t";
    }
    pout <<endl<<"Average: " << std::accumulate(PTenergy.begin(),PTenergy.end(),0.0)/world.size()<<endl;
    pout << "=============================================="<<endl;
    pout << "Average each component to compute PT2 energy" <<endl;

    base_H0   = std::accumulate(all_base_H0.begin(),all_base_H0.end(),0.0)/world.size();
    base_H01  = std::accumulate(all_base_H01.begin(),all_base_H01.end(),0.0)/world.size();
    perturber_H01   = std::accumulate(all_perturber_H01.begin(),all_perturber_H01.end(),0.0)/world.size();
    perturber_H0    = std::accumulate(all_perturber_H0.begin(),all_perturber_H0.end(),0.0)/world.size();
    overlap0    = std::accumulate(all_overlap0.begin(),all_overlap0.end(),0.0)/world.size();
    overlap1    = std::accumulate(all_overlap1.begin(),all_overlap1.end(),0.0)/world.size();
    base_H0_2   = std::accumulate(all_base_H0_2.begin(),all_base_H0_2.end(),0.0)/world.size();
    base_H01_2  = std::accumulate(all_base_H01_2.begin(),all_base_H01_2.end(),0.0)/world.size();
    perturber_H01_2   = std::accumulate(all_perturber_H01_2.begin(),all_perturber_H01_2.end(),0.0)/world.size();
    perturber_H0_2    = std::accumulate(all_perturber_H0_2.begin(),all_perturber_H0_2.end(),0.0)/world.size();
    overlap0_2    = std::accumulate(all_overlap0_2.begin(),all_overlap0_2.end(),0.0)/world.size();
    overlap1_2    = std::accumulate(all_overlap1_2.begin(),all_overlap1_2.end(),0.0)/world.size();
    cout <<"<0|1>: " << overlap0<<" "<<overlap0_2<<" "<<sqrt((overlap0_2-overlap0*overlap0)/(num_sample*world.size()))<<endl;
    cout <<"<1|0>: " << overlap1<<" "<<overlap1_2<<" "<<sqrt((overlap1_2-overlap1*overlap1)/(num_sample*world.size()))<<endl;
    cout <<"<0|(H_0-E_0)^(-1)|0>" << base_H0<<" "<<base_H0_2<<" "<<sqrt((base_H0_2-base_H0*base_H0)/(num_sample*world.size()))<<endl;
    cout <<"<1|(H_0-E_0)^(-1)|1>" << perturber_H0<<" "<<perturber_H0_2<<" "<<sqrt((perturber_H0_2-perturber_H0*perturber_H0)/(num_sample*world.size()))<<endl;
    cout <<"<1|(H_0-E_0)^(-1)|0>" << base_H01<<" "<<base_H01_2<<" "<<sqrt((base_H01_2-base_H01*base_H01)/(num_sample*world.size()))<<endl;
    cout <<"<0|(H_0-E_0)^(-1)|1>" << perturber_H01<<" "<<perturber_H01_2<<" "<<sqrt((perturber_H01_2-perturber_H01*perturber_H01)/(num_sample*world.size()))<<endl;
    double pt2_energy = -1.0*perturber_H0+perturber_H01*perturber_H01/(base_H0);
    cout <<"PT2 energy: " << pt2_energy<<endl;
    pt2_energy = -1.0*perturber_H0+base_H01*base_H01/(base_H0);
    cout <<"PT2 energy: " << pt2_energy<<endl;
  }
    double tcpu=sampletimer.totalcputime();
    double twall=sampletimer.totalwalltime();
    pout << setprecision(3) <<"\n\n\t\t\t SAMPLE CPU  Time (seconds): " << tcpu << endl;
    pout << setprecision(3) <<"\t\t\t SAMPLE Wall Time (seconds): " << twall << endl;



//  if(world.rank())
//  {
//  MPI_Reduce(&base_H0      , &test, 1, MPI_DOUBLE, MPI_SUM,0, Calc);
//  MPI_Reduce(&base_H01     , &base_H01, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(&perturber_H0 , &perturber_H0, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(&perturber_H01, &perturber_H01, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  }
//  else
//  {
//  MPI_Reduce(&base_H0      , &test, 1, MPI_DOUBLE, MPI_SUM,0, Calc);
//  //MPI_Allreduce(MPI_IN_PLACE, &base_H0, 1, MPI_DOUBLE, MPI_SUM, Calc);
//  //MPI_Reduce(MPI_IN_PLACE, &base_H0, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(MPI_IN_PLACE, &base_H01, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(MPI_IN_PLACE, &perturber_H0, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  MPI_Reduce(MPI_IN_PLACE, &perturber_H01, 1, MPI_DOUBLE, MPI_SUM, 0, Calc);
//  cout <<test<<"TEST"<<endl;
//
//  }
//  cout <<world.rank()<<" "<<"base_H0: " << base_H0<<endl;
  return 0;
}

