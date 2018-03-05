#include <vector>
#include "global.h"
#include "input.h"
#include "heatbath.h"
#include <ctime>
#include <random>
#include "IntegralMatrix.h"
#include <chrono>

using namespace SpinAdapted;
using namespace std;

void heatbath::precompute(){

  int num_spatial_orbs = dmrginp.spinAdapted()?dmrginp.last_site():dmrginp.last_site()/2;

  //Direct.ReSize(num_spatial_orbs, num_spatial_orbs);
  //Exchange.ReSize(num_spatial_orbs, num_spatial_orbs);
  Direct =Matrix(num_spatial_orbs, num_spatial_orbs);
  Exchange =Matrix(num_spatial_orbs, num_spatial_orbs);
  for(int p=0;p<num_spatial_orbs;p++)
    for(int q=0;q<num_spatial_orbs;q++)
    {
      Direct(p+1,q+1) = v_2[0](2*p,2*q,2*p,2*q);
      Exchange(p+1,q+1) = v_2[0](2*p,2*q,2*q,2*p);
    }

  p1_single.resize(num_spatial_orbs*2, 0.0);
  p1.resize(num_spatial_orbs*2, 0.0);
  p2.resize(num_spatial_orbs*2);
  p3.resize(num_spatial_orbs*2);
  p4.resize(num_spatial_orbs*2);


  for(int p=0;p<num_spatial_orbs*2;p++)
  {
    for(int q=0;q<num_spatial_orbs*2;q++)
      if (q!=p)
        p1_single[p] += fabs(v_1[integralIndex](p,q));
  }

  for(int p=0;p<num_spatial_orbs*2;p++)
  {
    for(int q=0;q<num_spatial_orbs*2;q++)
      if (q!=p)
      {
        for(int r=0;r<num_spatial_orbs*2;r++)
        {
          if (r!=p && r!=q)
          {
            for(int s=0;s<num_spatial_orbs*2;s++)
              if (s!=r && s!=p && s!=q)
              {
                p1[p] += fabs(v_2[integralIndex](r,s,p,q));
              }

          }
        }

      }
  }

  for(int p=0;p<num_spatial_orbs*2;p++)
  {
    p2[p].resize(num_spatial_orbs*2, 0.0);
    for(int q=0;q<num_spatial_orbs*2;q++)
      if(q!=p)
      {
        for(int r=0;r<num_spatial_orbs*2;r++)
          if (r!=p && r!=q)
          {
            for(int s=0;s<num_spatial_orbs*2;s++)
              if (s!=r && s!=p && s!=q)
                p2[p][q] += fabs(v_2[integralIndex](r,s,p,q));
          }

      }

  }

  for(int p=0;p<num_spatial_orbs*2;p++)
  {
    p3[p].resize(num_spatial_orbs*2);
    for(int q=0;q<num_spatial_orbs*2;q++)
      if(q!=p)
      {
        p3[p][q].resize(num_spatial_orbs*2, 0.0);
        //double sum = 0;
        for(int r=0;r<num_spatial_orbs*2;r++)
          if (r!=p && r!=q)
          {
            for(int s=0;s<num_spatial_orbs*2;s++)
            {
              if (s!=r && s!=p && s!=q)
                p3[p][q][r] += fabs(v_2[integralIndex](r,s,p,q));
            }
            //sum += p3[p][q][r];
          }
        //for(int r=0;r<num_spatial_orbs*2;r++)
        //{
        //    p3[p][q][r] /= sum;
        //}
      }
  }

  for(int p=0;p<num_spatial_orbs*2;p++)
  {
    p4[p].resize(num_spatial_orbs*2);
    for(int q=0;q<num_spatial_orbs*2;q++)
      {
        p4[p][q].resize(num_spatial_orbs*2, 0.0);
        for(int r=0;r<num_spatial_orbs*2;r++)
        {
          for(int s=0;s<num_spatial_orbs*2;s++)
          {
            if (s!=r && s!=p && s!=q)
              p4[p][q][r] += fabs(v_2[integralIndex](r,s,p,q));
          }
        }
      }
  }

  doubleH.resize(num_spatial_orbs*2);
  for(int p=0;p<num_spatial_orbs*2;p++)
  {
    doubleH[p].resize(num_spatial_orbs*2);
    for(int q=0;q<p;q++)
    {
      //doubleH.resize(num_spatial_orbs*2*num_spatial_orbs*2);
      for(int r=0;r<num_spatial_orbs*2;r++)
      for(int s=0;s<r;s++)
      {
        if(s!=p && s!=q && r!=p && r!=q )
        //if(s!=r &&s!=p && s!=q && r!=p && r!=q )
        {
          double h = v_2[integralIndex](r,s,p,q) -v_2[integralIndex](r,s,q,p)+v_2[integralIndex](s,r,q,p) -v_2[integralIndex](s,r,p,q) ;
        //doubleH[p][q].push_back(std::make_tuple(v_2[integralIndex](r,s,p,q), r, s));
          if(fabs(h)>tol)
            doubleH[p][q].push_back(std::make_tuple(h, r, s));

        }

      }
      std::sort(doubleH[p][q].begin(),doubleH[p][q].end(),[](std::tuple<double, int, int> a, std::tuple<double, int, int> b){return fabs(std::get<0>(a)) > abs(std::get<0>(b));});
    }
  }
}

    void heatbath::doubleexcite(const std::vector<int>& sd_in, double coeff, std::vector< std::vector<int> >& sd_out, std::vector<double>& out_coeff){
      return;
    }

    void heatbath::doubleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff){
      std::clock_t startTime = std::clock();
      double sample_p = 1.0;
      sd_out = sd_in;
      std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      std::vector<double> Sp;
      Sp.resize(occupy.size());
      for(int i=0;i<occupy.size();i++)
        Sp[i] = p1[occupy[i]];
      std::discrete_distribution<int> occupy_random(Sp.begin(),Sp.end());
      int p_i = occupy_random(generator);
      int p = occupy[p_i];
      //cout <<"Step 1: " << std::accumulate(Sp.begin(),Sp.end(), 0.0)<<endl;
      if(fabs(Sp[p_i] < NUMERICAL_ZERO))
      {
        out_coeff = 0.0;
        return;
      }
      else
        sample_p *=Sp[p_i]/std::accumulate(Sp.begin(),Sp.end(), 0.0);

      std::vector<double> Sq;
      Sq.resize(occupy.size());
      for(int i=0;i<Sq.size();i++)
        Sq[i] = p2[p][occupy[i]];
      Sq[p_i] = 0.0;
      std::discrete_distribution<int> occupy_random2(Sq.begin(),Sq.end());
      int q_i = occupy_random2(generator);
      int q = occupy[q_i];
      //cout <<"Step 2: " << std::accumulate(Sq.begin(),Sq.end(), 0.0)<<endl;
      if(fabs(Sq[q_i] < NUMERICAL_ZERO))
      {
        out_coeff = 0.0;
        return;
      }
        sample_p *=Sq[q_i]/std::accumulate(Sq.begin(),Sq.end(), 0.0);

      std::vector<double> Sr(unoccupy.size());
      for(int i=0;i<Sr.size();i++)
        Sr[i] = p3[p][q][unoccupy[i]];
      std::discrete_distribution<int> unoccupy_random(Sr.begin(),Sr.end());
      int r_i = unoccupy_random(generator);
      int r = unoccupy[r_i];
      //cout <<"Step 3: " << std::accumulate(Sr.begin(),Sr.end(), 0.0)<<endl;

      if(fabs(Sr[r_i] < NUMERICAL_ZERO))
      {
        out_coeff = 0.0;
        return;
      }
        sample_p *=Sr[r_i]/std::accumulate(Sr.begin(),Sr.end(), 0.0);

      std::vector<double> Ss(unoccupy.size(), 0.0);
      for(int i=0;i<Ss.size();i++)
        Ss[i] = fabs(v_2[integralIndex](r, unoccupy[i], p, q));
      Ss[r_i] = 0.0;
      std::discrete_distribution<int> unoccupy_random2(Ss.begin(),Ss.end());
      int s_i = unoccupy_random2(generator);
      int s = unoccupy[s_i];
      //cout <<"Step 4: " << Ss[s_i] <<" "<<std::accumulate(Ss.begin(),Ss.end(), 0.0)<<endl;
      if(fabs(Ss[s_i] < NUMERICAL_ZERO))
      {
        out_coeff = 0.0;
        return;
      }
        sample_p *=Ss[s_i]/std::accumulate(Ss.begin(),Ss.end(), 0.0);
      //double h = v_2[integralIndex](r, s, p, q) -v_2[integralIndex](r, s, q, p);
      double h = 0.5*v_2[integralIndex](r, s, p, q);
      out_coeff = coeff*h/sample_p;
      int nelec = 0;
      nelec = std::accumulate(sd_out.begin()+min(p,q)+1, sd_out.begin()+max(p,q), 0);
      sd_out[p] = 0;
      sd_out[q] = 0;
      nelec = std::accumulate(sd_out.begin()+min(r,s)+1, sd_out.begin()+max(r,s), nelec);
      sd_out[r] = 1;
      sd_out[s] = 1;
      if (nelec%2)
        out_coeff *=-1;
      //cout <<"Sample_p double: " << sample_p<<endl;
      //cout <<"Sample_p double: " << out_coeff<<endl;
      //cout <<"Sample_p double: " << h<<endl;
      //cout << r<<' '<<s <<"<-"<<p<<' '<<q<<endl;
      if(r>s) out_coeff *=-1;
      if(p>q) out_coeff *=-1;
      doubleT+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
    }

    void heatbath::uniformdoubleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff){
      std::clock_t startTime = std::clock();
      sd_out = sd_in;
      std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }

      std::uniform_int_distribution<int> occupy_random1(0,occupy.size()-1);
      std::uniform_int_distribution<int> occupy_random2(0,occupy.size()-2);
      std::uniform_int_distribution<int> unoccupy_random1(0,unoccupy.size()-1);
      std::uniform_int_distribution<int> unoccupy_random2(0,unoccupy.size()-2);
      int p_i = occupy_random1(generator);
      int p = occupy[p_i];
      occupy.erase(occupy.begin()+p_i);
      int q_i = occupy_random2(generator);
      int q = occupy[q_i];
      int r_i = unoccupy_random1(generator);
      int r = unoccupy[r_i];
      unoccupy.erase(unoccupy.begin()+r_i);
      int s_i = unoccupy_random2(generator);
      int s = unoccupy[s_i];

      double h = 0.5*v_2[integralIndex](r, s, p, q);
      out_coeff = coeff*h*occupy.size()*(occupy.size()+1)*unoccupy.size()*(unoccupy.size()+1);
      int nelec = 0;
      nelec = std::accumulate(sd_out.begin()+min(p,q)+1, sd_out.begin()+max(p,q), 0);
      sd_out[p] = 0;
      sd_out[q] = 0;
      nelec = std::accumulate(sd_out.begin()+min(r,s)+1, sd_out.begin()+max(r,s), nelec);
      sd_out[r] = 1;
      sd_out[s] = 1;
      if (nelec%2)
        out_coeff *=-1;
      //cout << r<<' '<<s <<"<-"<<p<<' '<<q<<endl;
      if(r>s) out_coeff *=-1;
      if(p>q) out_coeff *=-1;
      //cout <<"Sample_p double: " << out_coeff<<endl;
      doubleT+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
    }

/*
    void heatbath::doubleandsingleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff){
      sd_out = sd_in;
      double sample_p = 1.0;
      int integralIndex = 0;
      std::default_random_engine generator;
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==0) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      std::vector<double> Sp;
      Sp.resize(occupy.size());
      for(int i=0;i<occupy.size();i++)
        Sp[i] = p1[occupy[i]];
      std::discrete_distribution<int> occupy_random(Sp.begin(),Sp.end());
      int p_i = occupy_random(generator);
      int p = occupy[p_i];
      std::vector<double> Sq(occupy.size());
      sample_p *= Sp[p_i]/std::accumulate(Sp.begin(),Sp.end());

      for(int i=0;i<Sq.size();i++)
        Sq[i] = p2[p][occupy[i]];
      Sq[p_i] = 0.0;
      std::discrete_distribution<int> occupy_random2(Sq.begin(),Sq.end());
      int q_i = occupy_random2(generator);
      int q = occupy[q_i];
      sample_p *= Sq[q_i]/std::accumulate(Sq.begin(),Sq.end());

      sd_out[p] = 0;
      sd_out[q] = 0;
      if (p>q) out_coeff *=-1;

      std::vector<double> Sr(sd_in.size());
      for(int i=0;i<Sr.size();i++)
        Sr[i] = p3[p][q][unoccupy[i]];
      std::discrete_distribution<int> unoccupy_random(Sr.begin(),Sr.end());
      int r_i = unoccupy_random(generator);
      int r = unoccupy[r_i];
      sample_p *= Sr[r_i]/std::accumulate(Sr.begin(),Sr.end());

      std::vector<double> Ss(sd_in.size(), 0.0);
      for(int i=0;i<Ss.size();i++)
      {
        if (sd_out[i]==0)
          Ss[i] = v_2[integralIndex](r, i, p, q);
      }
      Ss[r] = 0.0;
      std::discrete_distribution<int> unoccupy_random2(Ss.begin(),Ss.end());
      int s = unoccupy_random2(generator);
      sample_p *= Ss[s_i]/std::accumulate(Ss.begin(),Ss.end(), 0.0);
      out_coeff = coeff*(v_2[integralIndex](r, s, p, q)/sample_p);
      if (r>s) out_coeff *=-1;
      sd_out[r] = 1;
      sd_out[s] = 1;
    }
    */
    void heatbath::singleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff){
      std::clock_t startTime = std::clock();
      double sample_p = 1.0;
      sd_out = sd_in;
      std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      std::vector<double> Sp;
      Sp.resize(occupy.size());
      for(int i=0;i<occupy.size();i++)
        Sp[i] = p1_single[occupy[i]];
      std::discrete_distribution<int> occupy_random(Sp.begin(),Sp.end());
      int p_i = occupy_random(generator);
      int p = occupy[p_i];
      if (fabs(Sp[p_i])<NUMERICAL_ZERO)
      {
        out_coeff =0.0;
        return;
      }
      sample_p *=Sp[p_i]/std::accumulate(Sp.begin(),Sp.end(), 0.0);

      std::vector<double> Sq(sd_out.size(), 0.0);
      for(int i=0;i<Sq.size();i++)
        if (sd_out[i]==0)
          Sq[i] = fabs(v_1[integralIndex](i,p));

      sd_out[p] = 0;

      std::discrete_distribution<int> unoccupy_random(Sq.begin(),Sq.end());
      int q = unoccupy_random(generator);
      sd_out[q] = 1;
      if (fabs(Sq[q])< NUMERICAL_ZERO)
      {
        out_coeff =0.0;
        return;
      }
      sample_p *= Sq[q]/std::accumulate(Sq.begin(),Sq.end(), 0.0);

      double h = v_1[integralIndex](q,p);
      //cout <<"h" << h<<endl;
      //pout <<"h" <<h<<endl;
      for(int i=0; i< occupy.size();i++)
      {
        if (i!=p_i)
        {
        h += v_2[integralIndex](q, occupy[i], p, occupy[i]);
        h -= v_2[integralIndex](q, occupy[i], occupy[i], p);
        }
      }
      //pout <<"h" <<h<<endl;

      out_coeff = coeff*h/sample_p;
      //cout <<"pq" <<p << " " <<q<<endl;
      if (p!=q)
      {
        int nelec = std::accumulate(sd_out.begin()+min(p,q)+1, sd_out.begin()+max(p,q), 0);
        if (nelec%2)
          out_coeff *=-1;
      }
      //cout <<"h" << h<<endl;
      //cout <<"Sample_p single: " << sample_p<<endl;
      //cout <<"out_coeff" <<out_coeff<<endl;
      singleT+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
    }
    void heatbath::uniformsingleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff){
      std::clock_t startTime = std::clock();
      sd_out = sd_in;
      std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      std::uniform_int_distribution<int> occupy_random(0,occupy.size()-1);
      std::uniform_int_distribution<int> unoccupy_random(0,unoccupy.size()-1);
      int p_i = occupy_random(generator);
      int q_i = unoccupy_random(generator);
      int p = occupy[p_i];
      int q = unoccupy[q_i];
      sd_out[p] = 0;
      sd_out[q] = 1;

      double h = v_1[integralIndex](q,p);
      for(int i=0; i< occupy.size();i++)
      {
        if (i!=p_i)
        {
        h += v_2[integralIndex](q, occupy[i], p, occupy[i]);
        h -= v_2[integralIndex](occupy[i], q, p, occupy[i]);
        }
      }
      //cout <<"h" << h<<endl;

      out_coeff = coeff*h*occupy.size()*unoccupy.size();
      //cout <<"pq" <<p << " " <<q<<endl;
      if (p!=q)
      {
        int nelec = std::accumulate(sd_out.begin()+min(p,q)+1, sd_out.begin()+max(p,q), 0);
        if (nelec%2)
          out_coeff *=-1;
      }
      //cout <<"h" << h<<endl;
      //if (fabs(h)>NUMERICAL_ZERO)
      //{
      ////cout <<"in_coeff" <<coeff<<endl;
      ////cout <<"out_coeff" <<out_coeff<<endl;
      //  for (auto i: sd_in)
      //    cout <<i;
      //  cout <<endl;
      //  for (auto i: sd_out)
      //    cout <<i;
      //  cout <<endl;

      //}
      singleT+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
    }
    void heatbath::exactsingleexcite(const std::vector<int>& sd_in, double coeff,  std::vector<int>& sd_out, double& out_coeff){
      std::clock_t startTime = std::clock();
      sd_out = sd_in;
      std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      std::vector<double> P(occupy.size()*unoccupy.size());
      for (int i=0;i <occupy.size();i++)
      for (int j=0;j <unoccupy.size();j++)
        P[i*unoccupy.size()+j] =fabs(v_1[integralIndex](occupy[i],unoccupy[j]));
      std::discrete_distribution<int> randomij(P.begin(),P.end());
      int ij = randomij(generator);
      int p_i = ij/unoccupy.size();
      int q_i = ij%unoccupy.size();
      int p = occupy[p_i];
      int q = unoccupy[q_i];
      sd_out[p] = 0;
      sd_out[q] = 1;

      double h = v_1[integralIndex](q,p);
      for(int i=0; i< occupy.size();i++)
      {
        if (i!=p_i)
        {
        h += v_2[integralIndex](q, occupy[i], p, occupy[i]);
        h -= v_2[integralIndex](occupy[i], q, p, occupy[i]);
        }
      }
      //cout <<"h" << h<<endl;

      //cout <<"h " <<h<<endl;
      //cout <<"P " <<P[ij]<<endl;
      out_coeff = coeff*(h/P[ij])*std::accumulate(P.begin(),P.end(),0.0);
      //out_coeff = coeff*(h/fabs(h))*std::accumulate(P.begin(),P.end(),0.0);
      //cout <<"pq" <<p << " " <<q<<endl;
      if (p!=q)
      {
        int nelec = std::accumulate(sd_out.begin()+min(p,q)+1, sd_out.begin()+max(p,q), 0);
        if (nelec%2)
          out_coeff *=-1;
      }
      //cout <<"h" << h<<endl;
      //if (fabs(h)>NUMERICAL_ZERO)
      //{
      ////cout <<"in_coeff" <<coeff<<endl;
      ////cout <<"out_coeff" <<out_coeff<<endl;
      //  for (auto i: sd_in)
      //    cout <<i;
      //  cout <<endl;
      //  for (auto i: sd_out)
      //    cout <<i;
      //  cout <<endl;

      //}
      singleT+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
    }

    void heatbath::singledoubleexcite(const std::vector<int>& sd_in, double coeff,  std::vector< std::vector<int> >& sd_out, std::vector<double>& out_coeff){
      double sample_p = 1.0;
      std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      std::vector<double> Sp;
      Sp.resize(occupy.size());
      for(int i=0;i<occupy.size();i++)
        Sp[i] = p1[occupy[i]];
      std::discrete_distribution<int> occupy_random(Sp.begin(),Sp.end());
      int p_i = occupy_random(generator);
      int p = occupy[p_i];
      //cout <<"Step 1: " << std::accumulate(Sp.begin(),Sp.end(), 0.0)<<endl;
      if(fabs(Sp[p_i] < NUMERICAL_ZERO))
      {
        return;
      }
      else
        sample_p *=Sp[p_i]/std::accumulate(Sp.begin(),Sp.end(), 0.0);

      std::vector<double> Sq;
      Sq.resize(occupy.size());
      for(int i=0;i<Sq.size();i++)
        Sq[i] = p2[p][occupy[i]];
      Sq[p_i] = 0.0;
      std::discrete_distribution<int> occupy_random2(Sq.begin(),Sq.end());
      int q_i = occupy_random2(generator);
      int q = occupy[q_i];
      //cout <<"Step 2: " << std::accumulate(Sq.begin(),Sq.end(), 0.0)<<endl;
      if(fabs(Sq[q_i] < NUMERICAL_ZERO))
      {
        return;
      }
        sample_p *=Sq[q_i]/std::accumulate(Sq.begin(),Sq.end(), 0.0);

      std::vector<double> Sr(unoccupy.size());
      for(int i=0;i<Sr.size();i++)
        Sr[i] = p3[p][q][unoccupy[i]];
      std::discrete_distribution<int> unoccupy_random(Sr.begin(),Sr.end());
      int r_i = unoccupy_random(generator);
      int r = unoccupy[r_i];
      //cout <<"Step 3: " << std::accumulate(Sr.begin(),Sr.end(), 0.0)<<endl;

      if(fabs(Sr[r_i] < NUMERICAL_ZERO))
      {
        return;
      }
      sample_p *=Sr[r_i]/std::accumulate(Sr.begin(),Sr.end(), 0.0);

      double H2sum=0.0;
      for (int i=0;i<unoccupy.size();i++)
      {
        if (i!=r_i)
          H2sum += 0.5*fabs(v_2[integralIndex](r,unoccupy[i], p, q));
      }
      double H1sum = v_1[integralIndex](r,p);
      for (int i=0;i<occupy.size();i++)
      {
        if(occupy[i]!=p)
        {
          H1sum += v_2[integralIndex](r,occupy[i], p, occupy[i]);
          H1sum -= v_2[integralIndex](r,occupy[i], occupy[i], p);
        }
      }

      if (true){
      //if (fabs(H1sum)>H2sum){
        //cout <<"SINGLE and DOUBLE"<<endl;
        sd_out.push_back(sd_in);
        out_coeff.push_back(coeff);
        out_coeff.back() /= sample_p;
        out_coeff.back() *=Sq[q_i]/std::accumulate(Sq.begin(),Sq.end(), 0.0);
        out_coeff.back() *=H1sum*changedeterm(sd_out.back(),r,p);
        //TODO
        //out_coeff.back() = 0.0;

        std::vector<double> Ss(unoccupy.size(), 0.0);
        for(int i=0;i<Ss.size();i++)
          Ss[i] = fabs(v_2[integralIndex](r, unoccupy[i], p, q));
        Ss[r_i] = 0.0;
        std::discrete_distribution<int> unoccupy_random2(Ss.begin(),Ss.end());
        int s_i = unoccupy_random2(generator);
        int s = unoccupy[s_i];
        if(fabs(Ss[s_i] < NUMERICAL_ZERO))
        {
          cout <<"Invalid double excitation"<<endl;
          return;
        }
        sample_p *=Ss[s_i]/std::accumulate(Ss.begin(),Ss.end(), 0.0);
        sd_out.push_back(sd_in);
        out_coeff.push_back(coeff);
        out_coeff.back() /=sample_p;
        out_coeff.back() *= 0.5*v_2[integralIndex](r, s, p, q);
        out_coeff.back() *= changedeterm(sd_out.back(), r,s,q,p);
      }
      else{
        std::uniform_real_distribution<double> uniformrandom(0.0,1.0);
        //if (uniformrandom(generator)<(fabs(H1sum))/(abs(H1sum)+H2sum))
        if (false)
        {
          cout <<"SINGLE"<<endl;
          sd_out.push_back(sd_in);
          out_coeff.push_back(coeff);
          out_coeff.back() /= sample_p;
          out_coeff.back() *=Sq[q_i]/std::accumulate(Sq.begin(),Sq.end(), 0.0);
          out_coeff.back() *=H1sum*changedeterm(sd_out.back(),r,p);
          out_coeff.back()*= (fabs(H1sum)+H2sum)/abs(H1sum);
          out_coeff.back() = 0.0;
        }
        else{
          cout <<"DOUBLE"<<endl;
          //sample_p *= 1-(fabs(H1sum)/(abs(H1sum)+H2sum));
          std::vector<double> Ss(unoccupy.size(), 0.0);
          for(int i=0;i<Ss.size();i++)
            Ss[i] = fabs(v_2[integralIndex](r, unoccupy[i], p, q));
          Ss[r_i] = 0.0;
          std::discrete_distribution<int> unoccupy_random2(Ss.begin(),Ss.end());
          int s_i = unoccupy_random2(generator);
          int s = unoccupy[s_i];
          if(fabs(Ss[s_i] < NUMERICAL_ZERO))
          {
            return;
          }
          sample_p *=Ss[s_i]/std::accumulate(Ss.begin(),Ss.end(), 0.0);
          sd_out.push_back(sd_in);
          out_coeff.push_back(coeff);
          out_coeff.back() /=sample_p;
          out_coeff.back() *= 0.5*v_2[integralIndex](r, s, p, q);
          out_coeff.back() *= changedeterm(sd_out.back(), r,s,q,p);
          }
        return;
      }

    }

    void heatbath::allexcite(const std::vector<int>& sd_in, double coeff,  std::unordered_map<longbitarray, double>& sd_table, double tol){
      std::clock_t startTime = std::clock();
      double sample_p = 1.0;
      std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      longbitarray sdbits(sd_in);
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      for(int p: occupy)
        for(int q: occupy)
          if(p!=q)
          {
            for(int i=0;i<doubleH[p][q].size();i++)
            {
              if(fabs(std::get<0>(doubleH[p][q][i])) < abs(tol/coeff))
                break;
              if(sdbits.getocc(std::get<1>(doubleH[p][q][i]))) continue;
              if(sdbits.getocc(std::get<2>(doubleH[p][q][i]))) continue;
              longbitarray newsdbits = sdbits;
              newsdbits.unset(p);
              newsdbits.unset(q);
              newsdbits.set(std::get<1>(doubleH[p][q][i]));
              newsdbits.set(std::get<2>(doubleH[p][q][i]));
              int factor = permutefactor(sd_in, std::get<1>(doubleH[p][q][i]), std::get<2>(doubleH[p][q][i]), q, p);
              sd_table[newsdbits] += 0.5*coeff*std::get<0>(doubleH[p][q][i])*factor;
            }
          }
      for(int p: occupy)
      {
        for(int r: unoccupy)
        {

          double h = v_1[integralIndex](r,p);
          for(int q: occupy)
          {
            h += q==p? 0.0: v_2[integralIndex](r,q, p, q);
            h -= p==q? 0.0: v_2[integralIndex](r,q, q, p);
          }
          if(fabs(h)> abs(tol/coeff))
          {
            longbitarray newsdbits = sdbits;
            newsdbits.unset(p);
            newsdbits.set(r);
            int factor = permutefactor(sd_in, r, p);
            sd_table[newsdbits] += coeff*h*factor;
          }
        }
      }
      excitation_T+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;

    }

    void heatbath::allexcite(const bitstring& sd_in, double coeff,  std::unordered_map<bitstring, double>& sd_table, double tol){
      std::clock_t startTime = std::clock();
      double sample_p = 1.0;
      //std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.size();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      for(int p: occupy)
        for(int q: occupy)
          if(p!=q)
          {
            for(int i=0;i<doubleH[p][q].size();i++)
            {
              if(fabs(std::get<0>(doubleH[p][q][i])) < abs(tol/coeff))
                break;
              if(sd_in[std::get<1>(doubleH[p][q][i])]) continue;
              if(sd_in[std::get<2>(doubleH[p][q][i])]) continue;
              bitstring newsdbits = sd_in;
              int factor = newsdbits.excitation(std::get<1>(doubleH[p][q][i]), std::get<2>(doubleH[p][q][i]), q, p);
              sd_table[newsdbits] += 0.5*coeff*std::get<0>(doubleH[p][q][i])*factor;
            }
          }
      for(int p: occupy)
      {
        for(int r: unoccupy)
        {

          double h = v_1[integralIndex](r,p);
          for(int q: occupy)
          {
            h += q==p? 0.0: v_2[integralIndex](r,q, p, q);
            h -= p==q? 0.0: v_2[integralIndex](r,q, q, p);
          }
          if(fabs(h)> abs(tol/coeff))
          {
            bitstring newsdbits = sd_in;
            int factor = newsdbits.excitation(r, p);
            sd_table[newsdbits] += coeff*h*factor;
          }
        }
      }

      excitation_T+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
    }

    void heatbath::allexcite(const bitstring& sd_in, double coeff,  std::unordered_map<bitstring, double>& sd_table, std::unordered_map<bitstring, double>& det_energy, double tol){
      std::clock_t startTime = std::clock();
      double energy = det_energy[sd_in];
      assert(fabs(energy)>1.0e-8);
      //double sample_p = 1.0;
      //std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
      std::vector<int> occupy;
      std::vector<int> unoccupy;
      for(int i=0;i<sd_in.length();i++)
      {
        if(sd_in[i]==1) occupy.push_back(i);
        else unoccupy.push_back(i);
      }
      for(int p: occupy)
        for(int q: occupy)
          if(p>q)
          {
            for(int i=0;i<doubleH[p][q].size();i++)
            {
              //if(fabs(std::get<0>(doubleH[p][q][i])) < abs(tol/coeff))
              if(fabs(std::get<0>(doubleH[p][q][i])) < abs(tol))
                break;
              int r = std::get<1>(doubleH[p][q][i]);
              int s = std::get<2>(doubleH[p][q][i]);
              if(sd_in[r] || sd_in[s]) continue;
              bitstring newsdbits = sd_in;
              std::clock_t startTime1 = std::clock();
              int factor = newsdbits.excitation(r, s, q, p);
              factor_T+= (std::clock() -startTime1)/(double) CLOCKS_PER_SEC;
              sd_table[newsdbits] += 0.5*coeff*std::get<0>(doubleH[p][q][i])*factor;
              if(det_energy.find(newsdbits)==det_energy.end())
              {

                if(r%2!=p%2)
                  swap(r,s);
                if(r%2!=p%2) {cout <<"pqrs"<<p<<" "<<q<<" "<<r<<" "<<s<<endl;abort();}
                assert(r%2==p%2);
                //det_energy[newsdbits] = local_energy(newsdbits);
                det_energy[newsdbits] = EnergyAfterExcitation(occupy, p, q, r, s, energy);
              }
              /*
              else{
                if(r%2!=p%2)
                  swap(r,s);
                assert(r%2==p%2);
                assert(fabs(det_energy[newsdbits] - EnergyAfterExcitation(occupy, p, q, r, s, energy))< 1.0e-8);
                if(fabs(det_energy[newsdbits] - EnergyAfterExcitation(occupy, p, q, r, s, energy))> 1.0e-8) {cout <<EnergyAfterExcitation(occupy, p, q, r, s, energy)<<det_energy[newsdbits]<<endl;abort();}

              }
              */
            }
          }

      excitation_T+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
      doubleT += (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
      startTime = std::clock();
      for(int p: occupy)
      {
        for(int r: unoccupy)
        {

          double h = v_1[integralIndex](r,p);
          for(int q: occupy)
          {
            h += q==p? 0.0: v_2[integralIndex](r,q, p, q);
            h -= p==q? 0.0: v_2[integralIndex](r,q, q, p);
          }
          if(fabs(h)> abs(tol/coeff))
          {
            bitstring newsdbits = sd_in;
            std::clock_t startTime1 = std::clock();
            int factor = newsdbits.excitation(r, p);
            factor_T+= (std::clock() -startTime1)/(double) CLOCKS_PER_SEC;
            sd_table[newsdbits] += coeff*h*factor;
            if(det_energy.find(newsdbits)==det_energy.end())
            {
              //det_energy[newsdbits] = local_energy(newsdbits);

              det_energy[newsdbits] = EnergyAfterExcitation(occupy, p, r, energy);
            }
            /*
            else
            {
                assert(fabs(det_energy[newsdbits] - EnergyAfterExcitation(occupy, p, r, energy))< 1.0e-8);
                if(fabs(det_energy[newsdbits] - EnergyAfterExcitation(occupy, p, r, energy))> 1.0e-8) {cout <<EnergyAfterExcitation(occupy, p, r, energy)<<det_energy[newsdbits]<<endl;abort();}
            }
            */
          }
        }
      }

      excitation_T+= (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
      singleT += (std::clock() -startTime)/(double) CLOCKS_PER_SEC;
    }

double heatbath::local_energy(const bitstring& ci, int integralIndex)
{
  std::clock_t startTime = std::clock();

  double energy = coreEnergy[integralIndex];
  int n = ci.size();

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
  energy_T += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
  return energy;

}


double heatbath::EnergyAfterExcitation(vector<int>& occupy, int i, int j, int a, int b, double Energyd) {
  /*!
     Calculates the new energy of a determinant after double excitation.

     .. note:: Assumes that i!=j and a!=b

     :Arguments:

       vector<int>& closed:
           Occupied orbitals in a vector.
       int i:
           Orbital index for destruction operator.
       int j:
           Orbital index for destruction operator.
       int a:
           Orbital index for creation operator.
       int b:
           Orbital index for creation operator.
       double Energyd:
           Old determinant energy.

     :Returns:

        double E:
            Energy after excitation.
   */

  std::clock_t startTime = std::clock();
  assert(i%2==a%2);
  assert(j%2==b%2);
  double E = Energyd;
  E += v_1[integralIndex](a, a)+v_1[integralIndex](b, b) -v_1[integralIndex](i, i)-v_1[integralIndex](j, j)   ;

/*
  for (int I: occupy) {
    if (I == i|| I ==j) continue;
    E +=  v_2[integralIndex](I,a,I,a)-v_2[integralIndex](I,i,I,i)+v_2[integralIndex](I,b,I,b)-v_2[integralIndex](I,j,I,j);
    E +=  -v_2[integralIndex](I,a, a, I)+v_2[integralIndex](I,i,i, I)-v_2[integralIndex](I,b, b, I)+v_2[integralIndex](I,j,j, I);
}
  E += v_2[integralIndex](a, b, a, b) - v_2[integralIndex](a, b, b, a) -v_2[integralIndex](i, j, i, j) +v_2[integralIndex](i, j, j, i);
*/
  for(int I: occupy)
  {
    if (I==i || I==j) continue;
    E += Direct(a/2+1, I/2+1) - Direct(i/2+1, I/2+1);
    if (i%2 == I%2)
      E += -Exchange(a/2+1, I/2+1) + Exchange(i/2+1, I/2+1);
  }
  for(int I: occupy)
  {
    if (I==i || I==j) continue;
    E += Direct(b/2+1, I/2+1) - Direct(j/2+1, I/2+1);
    if (j%2 == I%2)
      E += -Exchange(b/2+1, I/2+1) + Exchange(j/2+1, I/2+1);
  }
  E += Direct(a/2+1, b/2+1) - Direct(i/2+1, j/2+1);
  if (i%2==j%2)
  E += -Exchange(a/2+1, b/2+1) + Exchange(i/2+1, j/2+1);

  energy_T += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
  return E;
}

double heatbath::EnergyAfterExcitation(vector<int>& occupy, int i, int a, double Energyd) {
  /*!
     Calculates the new energy of a determinant after single excitation.

     .. note:: Assumes that the spin of i and a orbitals is the same

     :Arguments:

       vector<int>& closed:
           Occupied orbitals in a vector.
       int i:
           Orbital index for destruction operator.
       int a:
           Orbital index for creation operator.
       double Energyd:
           Old determinant energy.

     :Returns:

        double E:
            Energy after excitation.
   */

  std::clock_t startTime = std::clock();
  double E = Energyd;
  E += v_1[integralIndex](a, a)-v_1[integralIndex](i, i);

  for (int I :occupy) {
    if (I == i) continue;
    //E +=  v_2[integralIndex](I,a,I,a)-v_2[integralIndex](I,i,I,i);
    E +=  Direct(I/2+1, a/2+1) - Direct(I/2+1, i/2+1);
    if ( (I%2) == (i%2) )
    //E +=  -v_2[integralIndex](I,a, a, I)+v_2[integralIndex](I,i,i, I);
    E +=  -Exchange(I/2+1, a/2+1) + Exchange(I/2+1, i/2+1);
  }
  energy_T += (std::clock() -startTime)/ (double) CLOCKS_PER_SEC;
  return E;
}


