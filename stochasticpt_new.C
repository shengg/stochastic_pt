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
#include "stochasticpt_new.h"
#include "sampling.h"
#include "heatbath.h"
#include "simplemps.h"
using namespace SpinAdapted;
using namespace SpinAdapted::Sweep;
int bitstring::n_orb;
void ReadInput(char* conf);

void getComplementarySites(std::vector<int> &sites, std::vector<int> &complementarySites) ;



void compressedMPS(const unsigned long num_sample)
{
  //std::shared_ptr<NonAbelianmps> zeromps, sampledmps;
  std::shared_ptr<simplemps> zeromps, sampledmps;
  if(dmrginp.spinAdapted())
  {
    std::shared_ptr<NonAbelianmps> mps1 = std::make_shared<NonAbelianmps>();
    std::shared_ptr<NonAbelianmps> mps2 = std::make_shared<NonAbelianmps>();
    //Use their own shared_ptr for serialization in boost broadcast;
    mps1->build(0);
    mps2->build(1000);
    zeromps = mps1;
    sampledmps = mps2;
  }
  else{
    std::shared_ptr<Abelianmps> mps1 = std::make_shared<Abelianmps>();
    std::shared_ptr<Abelianmps> mps2 = std::make_shared<Abelianmps>();
    //Use their own shared_ptr for serialization in boost broadcast;
    mps1->build(0);
    mps2->build(1000);
    zeromps = mps1;
    sampledmps = mps2;
  }
  pout <<"QV|0> norm"<<sampledmps->get_norm()<<endl;
  /*
  {
  heatbath baseheatbath;
  baseheatbath.precompute();
    //sampledmps.build(1000);
    double norm = 0.0;
    double largest = 0.0;
    long largest_n = 0;
    double h00 = 0.0;
    double h11 = 0.0;
    double h10 = 0.0;
    for(long n=0;n<pow(2,bitstring::n_orb);n++)
    {
      bitstring ci(n);
    int particle_num_a = 0;
    int particle_num_b = 0;
    for(int i=0;i<ci.size()/2;i++)
    {
      particle_num_a+= ci[2*i];
      particle_num_b+= ci[2*i+1];
    }

    //if (particle_num_a*2!=dmrginp.total_particle_number()) continue;
    //if (particle_num_a-particle_num_b!=dmrginp.molecule_quantum().get_s().getirrep() ) continue;

    double coeff = zeromps->getcoeff(ci);
    if(abs(coeff)>4e-3 && ci.Sz()==dmrginp.molecule_quantum().get_s().getirrep() && ci[0]==0 && ci[1]==1)
      cout <<ci<<": "<<coeff<<endl;
    if(abs(coeff)>4e-3 && ci.Sz()==dmrginp.molecule_quantum().get_s().getirrep() && ci[0]==1 && ci[1]==0)
    //if(abs(coeff)>4e-2 && ci.Sz()==dmrginp.molecule_quantum().get_s().getirrep())
      cout <<ci<<": "<<coeff<<endl;
    //double coeff = zeromps->getcoeff(ci);
    norm +=coeff*coeff;
    if(coeff*coeff>largest)
    {
      largest = coeff*coeff;
      largest_n = n;
    }
    h00 += coeff*coeff/local_energy(ci, 0);

      //std::unordered_map<bitstring, double> sd_hashtable1;
      //baseheatbath.allexcite(ci, 1.0,sd_hashtable1, 1e-13);
      //sd_hashtable1[ci] += local_energy(ci, 1);
      //double overlap=0.0;
      //for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      //{
      //  overlap += zeromps->getcoeff(iter1->first)*iter1->second;
      //}
      //h11 += overlap*overlap/local_energy(ci, 0);
      //h10 += overlap*zeromps->getcoeff(ci)/local_energy(ci, 0);
      //sd_hashtable1.clear();

    }
    cout <<"Norm: " <<norm<<endl;
    cout <<"Largest: " <<largest_n<<": "<<bitstring(largest_n)<<" :"<<largest<<endl;
    cout <<"h00"<<h00<<endl;
    cout <<"h11"<<h11<<endl;
    cout <<"h10"<<h10<<endl;
    cout <<"Begin sample"<<endl;
    for(int i=0;i<10;i++)
    {

      bitstring determinant;
      double coeff = sampledmps->sampling(determinant);
      cout <<"ul"<<determinant<<endl;
      cout <<"coeff*coeff"<< coeff*coeff<<endl;
      coeff = sampledmps->getcoeff(determinant);
      cout <<"coeff*coeff"<< coeff*coeff<<endl;
      coeff = sampledmps->getcoeff(determinant);
      cout <<"coeff*coeff"<< coeff*coeff/sampledmps->get_norm()<<endl;

    }
    abort();

  }
*/

  //get_slater();
  std::vector<int> ci;
  //const unsigned long num_sample = 20000;
  heatbath baseheatbath;
  baseheatbath.set_tol() = dmrginp.stochasticpt_tol();
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
      double coeff = zeromps->sampling(determinant);
      //if(fabs(coeff)<1e-15) continue;
      //cout <<"coeff"<<coeff<<endl;
      //cout <<"coeff^2"<<coeff*coeff<<endl;
      //double coeff = baseState.sampling(determinant);
      double local_energy0;
      if(det_energy.find(determinant)==det_energy.end())
      {
        local_energy0= local_energy(determinant, 0);
        det_energy[determinant] = local_energy0;
      }
      else{
        local_energy0= det_energy[determinant];

      }
      //baseheatbath.allexcite(determinant, 1.0,sd_hashtable1, fabs(tol*100));
      //sd_hashtable1[determinant] += local_energy(determinant, 1);
      //
      //double overlap=0.0;
      //for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      //{
      //  overlap +=  zeromps->getcoeff(iter1->first)*iter1->second;
      //}
      //sd_hashtable1.clear();
      //double temp = local_energy0;
      ////double temp = overlap/coeff;
      //H00 += temp/num_sample;
      //H00_2 += temp*temp/num_sample;
      H00 += 1/(local_energy0*num_sample);
      H00_2 += 1/(local_energy0*local_energy0*num_sample);
      //norm += 1/(coeff*coeff);
    }
    timer1 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();
    long num_coeff=0;
    //cout <<"Step 2"<<endl;
    for (int i=0;i<num_sample;i++)
    {

      bitstring determinant;
      //double prob;
      double coeff = sampledmps->sampling(determinant);
      coeff = sampledmps->getcoeff(determinant);
      //cout <<"small coeff" <<coeff<<endl;
      //if(fabs(coeff)<1e-15) continue;

      //cout <<"coeff"<<coeff<<endl;
      //cout <<"coeff^2"<<coeff*coeff/sampledmps->get_norm()<<endl;
      double e = local_energy(determinant, 0);

      tol = 1e-15;
      baseheatbath.allexcite(determinant, 1.0,sd_hashtable1, fabs(tol*100));
      //baseheatbath.allexcite(iter->first, 1.0,sd_hashtable1, fabs(tol/coeff));
      sd_hashtable1[determinant] += local_energy(determinant, 1);
      num_coeff += sd_hashtable1.size();

      double overlap=0.0;
      for(auto iter1=sd_hashtable1.begin();iter1!=sd_hashtable1.end();iter1++)
      {
        //assert(iter1->first.Sz()==determinant.Sz());
        overlap +=  zeromps->getcoeff(iter1->first)*iter1->second;
        //testmps.getcoeff(iter1->first);
      }
      double temp =sampledmps->get_norm()*(overlap*overlap/(coeff*coeff))/e;
      //cout <<"Norm: "<<sampledmps->get_norm()<<endl;
      //cout <<"Sz: "<<determinant.Sz()<<endl;
      //pout <<"Overlap" <<overlap<<endl;
      //pout <<"Coeff"<<coeff<<endl;
      //pout <<"Coeff"<<sampledmps->getcoeff(determinant)<<endl;
      //if(coeff/overlap >2.0 || coeff/overlap<0.5)
      //  pout <<"ERROR: "<<determinant<<endl;
      //cout <<"Overlap" <<overlap*overlap*sampledmps->get_norm()/(coeff*coeff)<<endl;
      H11 += temp/num_sample;
      H11_2 += temp*temp/num_sample;
      temp = sampledmps->get_norm()/e;
      test1 += temp/num_sample;
      test1_2 += temp*temp/num_sample;

      //temp = sampledmps->get_norm()*overlap/coeff;
      //temp = sampledmps->get_norm()*zeromps->getcoeff(determinant)/coeff;
      temp = sampledmps->get_norm()*overlap*zeromps->getcoeff(determinant)/(e*coeff*coeff);
      H01 += temp/num_sample;
      H01_2 += temp*temp/num_sample;


      sd_hashtable1.clear();

      if(i%(max(num_sample/10,(const unsigned long)1))==0)
      {
        cout <<"Process " <<mpigetrank()<<": Finished " <<i <<" samples"<<endl;
      }

      //sd_hashtable0[determinant] += 1.0;
    }
    timer2 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
    begintime = std::clock();

    pout <<"Num of coeff to compute: "<< num_coeff<<endl;
    timer3 += (std::clock() - begintime)/ (double) CLOCKS_PER_SEC;
  }

  cout <<"Processor "<<mpigetrank()<<": Approx sampling |c| time: " <<(std::clock() - startTime)/ (double) CLOCKS_PER_SEC<<endl;


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
    zeromps->print_timer();
    sampledmps->print_timer();
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
  bitstring::n_orb= dmrginp.spinAdapted()?dmrginp.last_site()*2:dmrginp.last_site();

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
}

