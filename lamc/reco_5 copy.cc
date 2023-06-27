#include "my_belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  using namespace std;
  void User_reco::hist_def( void )
  { extern BelleTupleManager* BASF_Histogram;    
    t1 = BASF_Histogram->ntuple ("sigma",
         "m_lamc m_sigm mass_ach momentum channel_ups channel_tagging channel_lam_tag channel_d_tag energy ecm ntr mrec2_r mrec2_p");
    t2 = BASF_Histogram->ntuple ("lam_lept",
      "m_lamc mass_ach momentum channel_ups channel_tagging channel_lam_tag channel_d_tag energy ecm ntr mrec2_r mrec2_p");
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

    std::vector<Particle> p, ap, k_p, k_m, pi_p, pi_m, pi0, gamma, all, e_p, e_m, mu_m, mu_p;

    //fill vectors
    makeProton(p, ap, 1);
    makeKPi(k_p, k_m, pi_p, pi_m, 1);
    makePi0(pi0);
    withEminCutPi0(pi0, 0.05);
    makeGamma(gamma);
    makeLepton(e_p, e_m, mu_p, mu_m, 1);

    //Cuts
    withDrDzCut(p, 1., 2.);
    withDrDzCut(k_m, 1., 2.);
    withDrDzCut(pi_p, 1., 2.);
    withDrDzCut(ap, 1., 2.);
    withDrDzCut(k_p, 1., 2.);
    withDrDzCut(pi_m, 1., 2.);
    withLeptonIdCut(e_p, e_m, mu_p, mu_m, 0.01, 0.1);
    withPSCut(mu_p, 1.);
    withPSCut(mu_m, 1.);
    withPSCut(e_p, 1.);
    withPSCut(e_m, 1.);

    withKaonIdCut(k_p, k_m, 0.6);
    withProtonIdCut(p, ap, 0.6);
    
    deepCopy(pi_p, all);
    deepCopy(pi_m, all);

    setGenHepInfoF(p);
    setGenHepInfoF(k_m);
    setGenHepInfoF(pi_p);
    
    //Undetected particles
    std::vector<Particle> lamc_p, lamc_m, lamct_p, lamct_m;
    std::vector<Particle> sigc_pp, sigc_mm, sigc0, asigc0;
    std::vector<Particle> lam, alam;
    std::vector<Particle> ups, rho, rho_2m, rho_2p, rho4, rho_ppm, rho_mmp;
    std::vector<Particle> D0, aD0, D_p, D_m;
    std::vector<Particle> k_s;
    
    combination(rho, m_ptypeRHO0, pi_p, pi_m);

    combination(rho_2p, m_ptypeRHO0, pi_p, pi_p);
    combination(rho_2m, m_ptypeRHO0, pi_m, pi_m);
    combination(rho_mmp, m_ptypeRHO0, pi_m, pi_m, pi_p);
    combination(rho_ppm, m_ptypeRHO0, pi_p, pi_p, pi_m);
    combination(rho4, m_ptypeRHO0, rho_2p, rho_2m);

    makeKs(k_s);
    makeLambda(lam,alam);

    for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l) {
      HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
      Vector3 P(l->px(),l->py(),0);
      V=V-ip_position;
      V.setZ(0.);
      if (abs(l->mass()-0.4977)>0.03 || V.perp()<0.1 ||
      cos(V.angle(P))<0.99 || l->mdstVee2().z_dist()>1.) {
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


    for(std::vector<Particle>::iterator l = alam.begin(); l!=alam.end(); ++l) {
      HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
      Vector3 P(l->px(),l->py(),0);
      V=V-ip_position;
      V.setZ(0.);
      double p_id;
      if (l->child(0).pType().mass()>l->child(1).pType().mass()) p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(0).mdstCharged()));
      else p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(1).mdstCharged()));                               
      if (abs(l->mass()-1.1157)>0.01 || l->mdstVee2().z_dist()>1. || p_id<0.6 ) {
        alam.erase(l); --l;
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

    combination(lamct_p, m_ptypeLAMC, p, k_m, pi_p, 0.1);
    combination(lamct_m, m_ptypeLAMC, ap, k_p, pi_m, 0.1);
    setUserInfo(lamct_p,  1);
    setUserInfo(lamct_m,  1);

    combination(lamct_p, m_ptypeLAMC, lam, pi_p, 0.1);
    combination(lamct_m, m_ptypeLAMC, alam, pi_m, 0.1);
    setUserInfo(lamct_p,  2);
    setUserInfo(lamct_m,  2);

    combination(lamct_p, m_ptypeLAMC, p, k_s, 0.1);
    combination(lamct_m, m_ptypeLAMC, ap, k_s, 0.1);
    setUserInfo(lamct_p,  3);
    setUserInfo(lamct_m,  3);

    combination(lamct_p, m_ptypeLAMC, p, k_s, rho, 0.1);
    combination(lamct_m, m_ptypeLAMC, ap, k_s, rho, 0.1);
    setUserInfo(lamct_p,  4);
    setUserInfo(lamct_m,  4);


    /*
    semileptonics mode
    */

    combination(lamc_p, m_ptypeLAMC, lam, e_p);
    combination(lamc_m, m_ptypeLAMC, alam, e_m);
    setUserInfo(lamc_p,  10);
    setUserInfo(lamc_m,  10);

    combination(lamc_p, m_ptypeLAMC, lam, mu_p);
    combination(lamc_m, m_ptypeLAMC, alam, mu_m);
    setUserInfo(lamc_p,  11);
    setUserInfo(lamc_m,  11);

    /*
    combination(lamc_p, m_ptypeLAMC, lam, pi_p, 0.1);
    combination(lamc_m, m_ptypeLAMC, alam, pi_m, 0.1);
    setUserInfo(lamc_p,  12);
    setUserInfo(lamc_m,  12);
    */

    /*
    Sigma_c
    */
   
    combination(sigc_pp, m_ptypeSIGC0, lamc_p, pi_p, 0.5);
    combination(sigc_mm, m_ptypeSIGC0, lamc_m, pi_m, 0.5);
    setUserInfo(sigc_pp,  11);
    setUserInfo(sigc_mm,  11);

    combination(sigc0, m_ptypeSIGC0, lamc_p, pi_m, 0.5);
    combination(asigc0, m_ptypeSIGC0, lamc_m, pi_p, 0.5);
    setUserInfo(sigc0,  12);
    setUserInfo(asigc0,  12);


    /*
    elec posi
    */

    combination(ups, m_ptypeUPS4, lamc_p, lamct_m);
    combination(ups, m_ptypeUPS4, lamc_m, lamct_p);
    setUserInfo(ups, 1);

    combination(ups, m_ptypeUPS4, lamc_p, lamct_m, rho);
    combination(ups, m_ptypeUPS4, lamc_m, lamct_p, rho);
    setUserInfo(ups, 2);

    combination(ups, m_ptypeUPS4, lamc_p, lamct_m, rho4);
    combination(ups, m_ptypeUPS4, lamc_m, lamct_p, rho4);
    setUserInfo(ups, 3);

    combination(ups, m_ptypeUPS4, lamc_m, D0, p);
    combination(ups, m_ptypeUPS4, lamc_p, aD0, ap);
    setUserInfo(ups, 4);
  
    combination(ups, m_ptypeUPS4, lamc_m, D0, p, rho);
    combination(ups, m_ptypeUPS4, lamc_p, aD0, ap, rho);
    setUserInfo(ups, 5);
    
    combination(ups, m_ptypeUPS4, lamc_m, D_p, p, pi_m);
    combination(ups, m_ptypeUPS4, lamc_p, D_m, ap, pi_p);
    setUserInfo(ups,  6);
    
    combination(ups, m_ptypeUPS4, lamc_m, D_p, p, rho_mmp);
    combination(ups, m_ptypeUPS4, lamc_p, D_m, ap, rho_ppm);
    setUserInfo(ups,  7);
    
    /*
    elec posi through sigmac
    */

    combination(ups, m_ptypeUPS4, sigc_pp, lamc_m, pi_m);
    combination(ups, m_ptypeUPS4, sigc_mm, lamc_p, pi_p);
    setUserInfo(ups, 11);

    combination(ups, m_ptypeUPS4, sigc0, lamc_m, pi_p);
    combination(ups, m_ptypeUPS4, asigc0, lamc_p, pi_m);
    setUserInfo(ups, 12);

    combination(ups, m_ptypeUPS4, sigc_pp, lamc_m, pi_m, rho);
    combination(ups, m_ptypeUPS4, sigc_mm, lamc_p, pi_p, rho);
    setUserInfo(ups, 13);

    combination(ups, m_ptypeUPS4, sigc0, lamc_m, pi_p, rho);
    combination(ups, m_ptypeUPS4, sigc0, lamc_p, pi_m, rho);
    setUserInfo(ups, 14);
  
    combination(ups, m_ptypeUPS4, sigc_mm, D0, p, pi_p);
    combination(ups, m_ptypeUPS4, sigc_pp, aD0, ap, pi_m);
    setUserInfo(ups, 15);
    
    combination(ups, m_ptypeUPS4, sigc_mm, D_p, p);
    combination(ups, m_ptypeUPS4, sigc_pp, D_m, ap);
    setUserInfo(ups,  16);
    
    combination(ups, m_ptypeUPS4, sigc_mm, D_p, p, rho);
    combination(ups, m_ptypeUPS4, sigc_pp, D_m, ap, rho);
    setUserInfo(ups,  17);
  
    
    for(int j=0; j<ups.size(); ++j) 
    {
      Particle u=ups[j];  
      short ntr=0;
      for(int jj=0; jj<all.size(); ++jj) 
      if (!checkSame(all[jj],u)) ntr++;
      Particle charm_tagging = u.child(0);
      Particle ach = u.child(1);
      double energy = pStar(u, elec, posi).e();
      double momentum = pStar(u, elec, posi).vect().mag();
      double m_lamc = 0;
      double m_sigm = 0;


      short channel = dynamic_cast<UserInfo&>(u.userInfo()).channel();
      short chach = dynamic_cast<UserInfo&>(ach.userInfo()).channel();
      
      short channel_d_tag = -1;
      short channel_lam_tag = -1;
      short channel_tagging = -1;
      
      double mass_ach = ach.mass();

      /*
      Почему-то этот кусок работает неправильно
      */

      VectorL beam = VectorL(elec + posi, 0, 0, elec - posi); 
      VectorL prec = u.p();
      VectorL mis = beam - prec;
      double mrec2 = mis.m2();

      VectorL pl = charm_tagging.p();
      VectorL mis2 = beam - (prec - pl);
      double mrec2_2 = mis2.m2();

      //

      if (channel > 10){
        Particle lam = charm_tagging.child(0);
        channel_tagging = dynamic_cast<UserInfo&>(lam.userInfo()).channel();
        m_lamc = lam.mass();
        m_sigm = charm_tagging.mass();
      }
      else {
        short channel_tagging = dynamic_cast<UserInfo&>(charm_tagging.userInfo()).channel();
        m_lamc = charm_tagging.mass();
      }

      if (channel <= 3 or (channel >= 11 and channel < 15)){
        channel_lam_tag = chach;
      }
      else {
        channel_d_tag = chach;
      }


      if (channel < 10){
        t2->column("m_lamc", m_lamc);     
        t2->column("mass_ach", mass_ach);
        t2->column("momentum", momentum);
        t2->column("channel_ups", channel);     
        t2->column("channel_tagging", channel_tagging);
        t2->column("channel_lam_tag", channel_lam_tag);
        t2->column("channel_d_tag", channel_d_tag);
        t2->column("energy", energy);
        t2->column("ecm", ecm);
        t2->column("ntr", ntr);
        t2->column("mrec2_r", mrec2);
        t2->column("mrec2_p", mrec2_2);
        t2->dumpData();
      }
      else{
        t1->column("m_sigm", m_sigm);
        t1->column("m_lamc", m_lamc);     
        t1->column("mass_ach", mass_ach);
        t1->column("momentum", momentum);
        t1->column("channel_ups", channel);     
        t1->column("channel_tagging", channel_tagging);
        t1->column("channel_lam_tag", channel_lam_tag);
        t1->column("channel_d_tag", channel_d_tag);
        t1->column("energy", energy);
        t1->column("ecm", ecm);
        t1->column("ntr", ntr);
        t1->column("mrec2_r", mrec2);
        t1->column("mrec2_p", mrec2_2);
        t1->dumpData();
      }

  *status = 1;
  }
  
  if (*status==1) nwritt++;
  
  }
#if defined(BELLE_NAMESPACE)
}
#endif
