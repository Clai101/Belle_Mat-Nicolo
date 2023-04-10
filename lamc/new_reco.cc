#include "my_belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  using namespace std;
  void User_reco::hist_def( void )
  { extern BelleTupleManager* BASF_Histogram;    
    t1 = BASF_Histogram->ntuple ("lam_p_k_pi",
				 "ml mach p chu chl chlt chlt chdt chdsc en ecm ntr");
  };
  
  
  int fill_tup(Particle lamc, /*vector<Particle> all,*/ double elec, double posi, double ecm, double r2, BelleTuple *t)
  {    
    //int chb = dynamic_cast<UserInfo&>(B.userInfo()).channel();
    
    return 1;
  };
  
  
  void User_reco::event ( BelleEvent* evptr, int* status ) {
    
    *status=0;
    
    static int nevent=0;
    static int nwritt=0;
    if(++nevent<2 || !(nevent%1000)) cout << "Event number " << nevent
					  << " selected" << nwritt << endl;
    
    Belle_runhead_Manager& rhdmgr = Belle_runhead_Manager::get_manager();
    Belle_runhead_Manager::const_iterator rhd = rhdmgr.begin();
    Evtcls_hadron_info_Manager&  ehimgr =
      Evtcls_hadron_info_Manager::get_manager();
    
    Evtcls_hadronic_flag_Manager&  ehadfl =
      Evtcls_hadronic_flag_Manager::get_manager();
    
    HepPoint3D ip_position = IpProfile::position();
    const HepSymMatrix& runIp_err = IpProfile::position_err();
    Mdst_vee2_Manager &vee2_mgr = Mdst_vee2_Manager::get_manager();
    Evtcls_hadron_info_Manager::iterator iti = ehimgr.begin();
    Evtcls_hadronic_flag_Manager::iterator eti = ehadfl.begin();
    
    double r2=0;
    int ntrk=0;
    double evis=0;
    double Pz=0;
    double hjmass=0;

    
    if (iti!=ehimgr.end()){
      r2 = (*iti).R2();
      //  ntrk = (*iti).Ntrk();
      //      evis = (*iti).Evis();
      //      Pz   = (*iti).Pz();
      //      hjmass = (*iti).HeavyJetMass();
    }
        
    double ecm = BeamEnergy::Ecm();
    double elec = BeamEnergy::E_HER();
    double posi = BeamEnergy::E_LER();
    
    /*************** Make particle lists ********************************/

    //Base particles

    std::vector<Particle> p, ap, k_p, k_m, pi_p, pi_m, pi0, gamma, all;

    //fill vectors
    makeProton(p, ap, 1);
    makeKPi(k_p, k_m, pi_p, pi_m, 1);
    makePi0(pi0);
    withEminCutPi0(pi0, 0.05);
    makeGamma(gamma);

    //Cuts
    withDrDzCut(p, 1., 2.);
    withDrDzCut(k_m, 1., 2.);
    withDrDzCut(pi_p, 1., 2.);
    withDrDzCut(ap, 1., 2.);
    withDrDzCut(k_p, 1., 2.);
    withDrDzCut(pi_m, 1., 2.);

    withKaonIdCut(k_p, k_m, 0.6);
    withProtonIdCut(p, ap, 0.6);
    
    deepCopy(pi_p, all);
    deepCopy(pi_m, all);

    setGenHepInfoF(p);
    setGenHepInfoF(k_m);
    setGenHepInfoF(pi_p);
    
    //Undetected particles
    std::vector<Particle> lamc_p, lamc_m;
    std::vector<Particle> lam, alam;
    std::vector<Particle> ups, rho, rho_2m, rho_2p, rho4;
    std::vector<Particle> D0, aD0, D_p, D_m;
    std::vector<Particle> k_s;
    
    combination(rho, m_ptypeRHO0, pi_p, pi_m);

    combination(rho_2p, m_ptypeRHO0, pi_p, pi_p);
    combination(rho_2m, m_ptypeRHO0, pi_m, pi_m);
    combination(rho_ppm, m_ptypeRHO0, pi_m, pi_m, pi_p);
    combination(rho_mmp, m_ptypeRHO0, pi_p, pi_p, pi_m);
    combination(rho4, m_ptypeRHO0, rho_2p, rho_2m);

    makeKs(k_s);
    makeLambda(lam,alam);

  for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l) {
    HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    if (abs(l->mass()-0.4977)>0.03 || V.perp()<0.1 ||
    cos(V.angle(P))<0.99 || l->mdstVee2().z_dist()>1. ) {
      k_s.erase(l); --l; continue;
    }
  }

   for(std::vector<Particle>::iterator l = lam.begin(); l!=lam.end(); ++l) {
      HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
      Vector3 P(l->px(),l->py(),0);
      V=V-ip_position;
      V.setZ(0.);
      double p_id;
      if (l->child(0).pType().mass()>l->child(1).pType().mass()) p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(0).mdstCharged()));
      else p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(1).mdstCharged()));                               
      if (abs(l->mass()-1.1157)>0.01 || l->mdstVee2().z_dist()>1. || p_id<0.6 ) {
        lam.erase(l); --l;
      }
    }


   for(std::vector<Particle>::iterator l = alam.begin(); l!=lam.end(); ++l) {
      HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
      Vector3 P(l->px(),l->py(),0);
      V=V-ip_position;
      V.setZ(0.);
      double p_id;
      if (l->child(0).pType().mass()>l->child(1).pType().mass()) p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(0).mdstCharged()));
      else p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(1).mdstCharged()));                               
      if (abs(l->mass()-1.1157)>0.01 || l->mdstVee2().z_dist()>1. || p_id<0.6 ) {
        lam.erase(l); --l;
      }
    }

    /*D*/

    combination(D0, m_ptypeD0, k_m, pi_p, 0.05);
    combination(aD0, m_ptypeD0B, k_p, pi_m, 0.05);
    setUserInfo(D0,  1);
    setUserInfo(aD0,  1);

    combination(D0, m_ptypeD0, k_m, rho_2p, pi_m, 0.05);
    combination(aD0, m_ptypeD0B, k_p, rho_2m, pi_p, 0.05);
    setUserInfo(D0,  2);
    setUserInfo(aD0,  2);

    combination(D0, m_ptypeD0, k_m, pi0, pi_p, 0.05);
    combination(aD0, m_ptypeD0B, k_p, pi0, pi_m, 0.05);
    setUserInfo(D0, 3);
    setUserInfo(aD0, 3);

    combination(D0, m_ptypeD0, k_m, k_p, 0.05);
    combination(aD0, m_ptypeD0B, k_p, k_m, 0.05);
    setUserInfo(D0, 4);
    setUserInfo(aD0, 4);

    combination(D0, m_ptypeD0, k_s, pi_p, pi_m, 0.05);
    combination(aD0, m_ptypeD0B, k_s, pi_p, pi_m, 0.05);
    setUserInfo(D0, 5);
    setUserInfo(aD0, 5);

    combination(D0, m_ptypeD0, k_s, pi0, 0.05);
    combination(aD0, m_ptypeD0B, k_s, pi0, 0.05);
    setUserInfo(D0, 6);
    setUserInfo(aD0, 6);

    combination(D_p, m_ptypeDP, k_m, pi_p, pi_p, 0.05);
    combination(D_m, m_ptypeDM, k_p, pi_m, pi_m, 0.05);
    setUserInfo(D_p,  11);
    setUserInfo(D_m,  11);

    combination(D_p, m_ptypeD0, k_s, pi_p, 0.05);
    combination(D_m, m_ptypeD0B, k_s, pi_m, 0.05);
    setUserInfo(D_m, 12);
    setUserInfo(D_p, 12);

    /*
    Lambda mesons
    */

    combination(lamc_p, m_ptypeLAMC, p, k_m, pi_p, 0.1);
    combination(lamc_m, m_ptypeLAMC, ap, k_p, pi_m, 0.1);
    setUserInfo(lamc_p,  1);
    setUserInfo(lamc_m,  1);

    combination(lamc_p, m_ptypeLAMC, lam, pi_p, 0.1);
    combination(lamc_m, m_ptypeLAMC, alam, pi_m, 0.1);
    setUserInfo(lamc_p,  2);
    setUserInfo(lamc_m,  2);

    combination(lamc_p, m_ptypeLAMC, p, k_s, 0.1);
    combination(lamc_m, m_ptypeLAMC, ap, k_s, 0.1);
    setUserInfo(lamc_p,  3);
    setUserInfo(lamc_m,  3);

    combination(lamc_p, m_ptypeLAMC, p, k_s, rho, 0.1);
    combination(lamc_m, m_ptypeLAMC, ap, k_s, rho, 0.1);
    setUserInfo(lamc_p,  4);
    setUserInfo(lamc_m,  4);

    /*
    Epsilon 4 pi
    */

    combination(ups, m_ptypeUPS4, lamc_p, lamc_m, 2.0);
    setUserInfo(ups, 1);

    combination(ups, m_ptypeUPS4, lamc_p, lamc_m, rho, 2.0);
    setUserInfo(ups, 2);

    combination(ups, m_ptypeUPS4, lamc_p, lamc_m, rho4, 2.0);
    setUserInfo(ups, 3);chl chlt chlt chdt chdst chdsc en ntr

    combination(ups, m_ptypeUPS4, lamc_m, D0, p, 2.0);
    combination(ups, m_ptypeUPS4, lamc_p, aD0, ap, 2.0);
    setUserInfo(ups, 4);
  
    combination(ups, m_ptypeUPS4, lamc_m, D0, p, rho, 2.0);
    combination(ups, m_ptypeUPS4, lamc_p, aD0, ap, rho, 2.0);
    setUserInfo(ups, 5);
    
    combination(ups, m_ptypeUPS4, lamc_m, D_p, p, pi_m, 2.0);
    combination(ups, m_ptypeUPS4, lamc_p, D_m, ap, pi_p, 2.0);
    setUserInfo(ups,  6);
    
    combination(ups, m_ptypeUPS4, lamc_m, D_p, p, pi_mmp, 2.0);
    combination(ups, m_ptypeUPS4, lamc_p, D_m, ap, pi_ppm, 2.0);
    setUserInfo(ups,  7);

  
    
    for(int j=0; j<ups.size(); ++j) {
        Particle u=ups[j];	
        int ntr=0;
        for(int jj=0; jj<all.size(); ++jj) 
	      if (!checkSame(all[jj],u)) ntr++;
        Particle lamc = u.child(0);
        Particle ach = u.child(1);
        double en = pStar(u, elec, posi).e();
        double p = pStar(u, elec, posi).vect().mag();
        double mass_lamc = lamc.mass();
        int chu = dynamic_cast<UserInfo&>(u.userInfo()).channel();
        int chl = dynamic_cast<UserInfo&>(lamc.userInfo()).channel();
        int chach = dynamic_cast<UserInfo&>(ach.userInfo()).channel();
        
        int chdsc = -1;
        int chdt = -1;
        int chlt = -1;
        
	double mass_ach = ach.mass();
	
        switch(chu){
          case 1:
            chlt = chach;
            break;
          case 2:
            chlt = chach;
            break;
          case 3:
            chlt = chach;
            break;
          case 4:
            chdt = chach;
            break;
          case 5:
            chdt = chach;
            break;
          case 6:
            chdt = chach;
            break;
          case 7:
            chdt = chach;
            break;
        }

        t1->column("ml", mass_lamc);     
        t1->column("mach", mass_ach);
        t1->column("mdst", mass_dst);

        t1->column("p", p);

        t1->column("chu", chu);     
        t1->column("chl", chl);
        t1->column("chlt", chlt);
        t1->column("chdt", chdt);
        t1->column("chdst", chdst);
        t1->column("chdsc", chdsc);

        t1->column("en", en);
        t1->column("ecm", ecm);
        t1->column("ntr", ntr);

        t1->dumpData();


	  *status = 1;
    }
    
    if (*status==1) nwritt++;
    
  }
#if defined(BELLE_NAMESPACE)
}
#endif
