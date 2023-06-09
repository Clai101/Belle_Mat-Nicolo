#include "my_belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  using namespace std;
  void User_reco::hist_def( void )
  { extern BelleTupleManager* BASF_Histogram;    
    t1 = BASF_Histogram->ntuple ("lam_p_k_pi",
				 "en p ml mach mdst chu chl chdt chdsc chlt chdst ntr");
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

    std::vector<Particle> p, ap, k_p, k_m, pi_p, pi_m, pi0, gamma, all;

    makeProton(p, ap, 1);
    makeKPi(k_p, k_m, pi_p, pi_m, 1);
    makePi0(pi0);
    withEminCutPi0(pi0, 0.05);
    makeGamma(gamma);
    withDrDzCut(p, 1., 2.);
    withDrDzCut(k_m, 1., 2.);
    withDrDzCut(pi_p, 1., 2.);
    withDrDzCut(ap, 1., 2.);
    withDrDzCut(k_p, 1., 2.);
    withDrDzCut(pi_m, 1., 2.);

    deepCopy(pi_p, all);
    deepCopy(pi_m, all);

    withKaonIdCut(k_p, k_m, 0.6);
    withProtonIdCut(p, ap, 0.6);
    //    withPionIdCut(pi_p, pi_m, 0.1);
    
    setGenHepInfoF(p);
    setGenHepInfoF(k_m);
    setGenHepInfoF(pi_p);
    
    /*makeGamma(gam);
    withEminCut(gam, 0.05);
    deepCopy(pi_p, all);
    deepCopy(pi_m, all);
    deepCopy(gam, all);*/
    
    std::vector<Particle> lamc_p, lamc_m;
    std::vector<Particle> lam, alam;
    std::vector<Particle> ups, rho, rho_2m, rho_2p, rho4, rho_m, rho_p;
    std::vector<Particle> D0, aD0, D_p, D_m;
    std::vector<Particle> D_star_p, D_star_m, D_star0, aD_star0;
    
    makeLambda(lam, alam);
    /* 

    Mesons

    */

    combination(rho, m_ptypeRHO0, pi_p, pi_m);

    combination(rho_2p, m_ptypeRHO0, pi_p, pi_p);
    combination(rho_2m, m_ptypeRHO0, pi_m, pi_m);
    combination(rho4, m_ptypeRHO0, rho_2p, rho_2m);

    //combination(rho4, m_ptypeRHO0, pi_p, pi_p, pi_m, pi_m);
    

    /*for(std::vector<Particle>::iterator l = rho4.begin(); l!=rho4.end(); ++l) {
      for(std::vector<Particle>::iterator j = l; j!=rho4.end(); ++j) {
        if (rho4[j].child(1).child(1) == rho4[l].child(0).child(1) && rho4[j].child(1).child(0) == rho4[l].child(0).child(0) &&
	    rho4[j].child(1).child(1) == rho4[l].child(0).child(1) && rho4[j].child(1).child(0) == rho4[l].child(0).child(0)) { 
          rho4.erase(j); --j; continue;}
        }
      }
    */

    combination(rho_p, m_ptypeRHO0, rho_2p, pi_m);
    
    /*for(std::vector<Particle>::iterator l = rho_p.begin(); l!=rho_p.end(); ++l) {
      for(std::vector<Particle>::iterator j = l; j!=rho_p.end(); ++j) {
        if (rho_p[j].child(0).child(0) == rho_p[l].child(1)) { 
          rho_p.erase(j); --j; continue;}
        }
      }*/

    combination(rho_m, m_ptypeRHO0, rho_2m, pi_p);   
    /*for(std::vector<Particle>::iterator l = rho_m.begin(); l!=rho_m.end(); ++l) {
      for(std::vector<Particle>::iterator j = l; j!=rho_m.end(); ++j) {
        if (rho_m[j].child(0).child(1) == rho_m[l].child(1)) { 
          rho_m.erase(j); --j; continue;}
        }
      }
    */

    /*D*/

    combination(D0, m_ptypeD0, k_m, pi_p, 0.05);
    combination(aD0, m_ptypeD0B, k_p, pi_m, 0.05);
    setUserInfo(D0,  1);
    setUserInfo(aD0,  1);

    combination(D0, m_ptypeD0, k_m, rho_p, 0.05);
    combination(aD0, m_ptypeD0B, k_p, rho_m, 0.05);
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

    /* D0 -> K+ K-*/

    combination(D_p, m_ptypeDP, k_m, pi_p, pi_p, 0.05);
    combination(D_m, m_ptypeDM, k_p, pi_m, pi_m, 0.05);
    setUserInfo(D_p,  11);
    setUserInfo(D_m,  11);

    /*D_star*/

    combination(D_star_p, m_ptypeDstarP, D0, pi_p, 0.05);
    combination(D_star_m, m_ptypeDstarM, D0, pi_m, 0.05);
    combination(D_star0, m_ptypeDstar0, D0, pi0, 0.05);
    combination(aD_star0, m_ptypeDstarB, aD0, pi0, 0.05);
    setUserInfo(D_star_p, 11);
    setUserInfo(D_star_m, 11);
    setUserInfo(D_star0, 1);
    setUserInfo(aD_star0, 1);

    combination(D_star_p, m_ptypeDstarP, D_p, pi0, 0.01);
    combination(D_star_m, m_ptypeDstarM, D_m, pi0, 0.01);
    combination(D_star0, m_ptypeDstar0, D0, gamma, 0.01);
    combination(aD_star0, m_ptypeDstarB, aD0, gamma, 0.01);
    setUserInfo(D_star_p, 12);
    setUserInfo(D_star_m, 12);
    setUserInfo(D_star0, 2);
    setUserInfo(aD_star0, 2);

    /*Lambda_b_0

    combination(Lm_b_0, m_ptypeDP, p, pi_m, 0.1);
    combination(aLm_b_0, m_ptypeDM, ap, pi_p, 0.1);
    setUserInfo(Delta_0, 1);
    setUserInfo(aDelta_p, 1);


    
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

    /*
    combination(lamc_p, m_ptypeLAMC, lam, pi_p, pi0, 0.1);
    combination(lamc_m, m_ptypeLAMC, alam, pi_m, pi0, 0.1);
    setUserInfo(lamc_p,  3);
    setUserInfo(lamc_m,  3);
    */

    /*combination(lam, m_ptypeLAM, p, pi_m, 0.1);
    combination(alam, m_ptypeLAM, ap, pi_p, 0.1);
    setUserInfo(lam,  1);
    setUserInfo(alam,  1);
    */
    /*
    Epsilon 4 pi
    */

    combination(ups, m_ptypeUPS4, lamc_p, lamc_m, 2.0);
    setUserInfo(ups, 1);

    combination(ups, m_ptypeUPS4, lamc_p, lamc_m, rho, 2.0);
    setUserInfo(ups, 2);

    combination(ups, m_ptypeUPS4, lamc_p, lamc_m, rho4, 2.0);
    setUserInfo(ups, 3);

    combination(ups, m_ptypeUPS4, lamc_m, D_star0, p , 2.0);
    combination(ups, m_ptypeUPS4, lamc_p, aD_star0, ap , 2.0);
    setUserInfo(ups, 4);

    combination(ups, m_ptypeLAMC, lamc_m, D_star_p, p, pi_m, 2.0);
    combination(ups, m_ptypeLAMC, lamc_p, D_star_m, ap, pi_p, 2.0);
    setUserInfo(ups, 5);

    combination(ups, m_ptypeUPS4, lamc_m, D0, p, 2.0);
    combination(ups, m_ptypeUPS4, lamc_p, aD0, ap, 2.0);
    setUserInfo(ups, 6);

    combination(ups, m_ptypeUPS4, lamc_m, D_p, p, pi_m, 2.0);
    combination(ups, m_ptypeUPS4, lamc_p, D_m, ap, pi_p, 2.0);
    setUserInfo(ups,  7);
    

    /*
    ach - anti c hadron
    chu - channel of ups decay
    chl - channel of lamc decay
    chach - channel of ach decay
    chdsc - channel of D* (d star) child (D meson)
    chdt - channel of D meson deacay (t for tag)
    chdst - channel of D* meson deacay (t for tag)
    chlt - channel of lamc, which is ach, decay
    */
    for(int j=0; j<ups.size(); ++j) {
        Particle u=ups[j];	
        int ntr=0;
        for(int jj=0; jj<all.size(); ++jj) 
	      if (!checkSame(all[jj],u)) ntr++;
        Particle lamc = u.child(0);
        Particle ach = u.child(1);
        Particle dsc = u.child(1).child(0);
        double en = pStar(u, elec, posi).e();
        double p = pStar(u, elec, posi).vect().mag();
        double mass_lamc = lamc.mass();
	      // double mass = ach.mass();
        int chu = dynamic_cast<UserInfo&>(u.userInfo()).channel();
        int chl = dynamic_cast<UserInfo&>(lamc.userInfo()).channel();
        int chach = dynamic_cast<UserInfo&>(ach.userInfo()).channel();
        
        int chdsc = -1;
        int chdt = -1;
        int chdst = -1;
        int chlt = -1;
        
	      double mass_dst = 0;
	      double mass_ach = 0;

      	if (ch == 4 or ch == 5) {
             chdsc = dynamic_cast<UserInfo&>(dsc.userInfo()).channel();
      	     mass_dst = ach.mass();
      	     mass_ach = dsc.mass();
      	}
      	else{
      	     mass_ach = ach.mass();
      	}

	      //double ecm = (en + p)*(en + p) - p*p;

        // cout << "in lamc cicle"  << mass << endl;

        /* double m12, m23, m13;
        m12 = (lamc.child(0).p() + lamc.child(1).p()).m();
        m13 = (lamc.child(0).p() + lamc.child(2).p()).m();
        m23 = (lamc.child(1).p() + lamc.child(2).p()).m();
        */

        switch(ch){
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
            chdst = chach;
            break;
          case 5:
            chdst = chach;
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
        t1->column("ntr", ntr);

        t1->dumpData();


	  *status = 1;
    }
    
    if (*status==1) nwritt++;
    
  }
#if defined(BELLE_NAMESPACE)
}
#endif
