#ifndef trackVars2_h
#define trackVars2_h



class trackVars2 {
      public:
    double pt, eta, phi, dz, dxy,  dist, dsig;
    double mass, t3Dip, t3Dsip, t2Dip, t2Dsip;
//     double t3DipSigned, t3DsipSigned, t2DipSigned, t2DsipSigned;
    double PCA_sx, PCA_sy, PCA_sz, PCA_sxerr, PCA_syerr, PCA_szerr;
    double PCA_tx, PCA_ty, PCA_tz, PCA_txerr, PCA_tyerr, PCA_tzerr;
    double dotprodTrack, dotprodSeed;
    double dotprodTrackSeed2D, dotprodTrackSeed3D;
    double dotprodTrackSeedVectors2D, dotprodTrackSeedVectors3D;
    double seedPCA_pv, trackPCA_pv;
    double seedMass;
    
    
//     double PCA_J;
    double PCA_JetAxis_distance, PCAPair_Jet_dotprod, PCAAxis_JetAxis_DEta, PCAAxis_JetAxis_DPhi; 
    
    int index, seed_index;
    
    void set_values (double, double, double, double,  double, double, double,
    double, double, double, double,  double, double,
    double, double, double, double,  double, double,
    double, double
    );
    
    void set_signedIPs(double, double, double, double);
    
    void set_vars (double, double, double, double,
    double, double, double, double, double
    );
    
    void set_index ( int );
    void set_SeedIndex(int);
    void set_distances ( double, double );
    
     void set_JetAxisVars(double, double, double, double);
     void setSeedMass(double);

};

inline void trackVars2::setSeedMass(double sm){
    seedMass=sm;
};

inline void trackVars2::set_values (double pt2, double eta2, double phi2, double dz2, double dxy2, double distaaa, double dsig2,
double PCA_sx2, double PCA_sy2, double PCA_sz2, double PCA_sxerr2, double PCA_syerr2, double PCA_szerr2, 
double PCA_tx2, double PCA_ty2, double PCA_tz2, double PCA_txerr2, double PCA_tyerr2, double PCA_tzerr2,
double dotprodTrack2, double dotprodSeed2) {

    pt=pt2;
    eta=eta2;
    phi=phi2; 
    dz=dz2; 
    dxy=dxy2; 
    dist=distaaa;
    dsig=dsig2;
    PCA_sx=PCA_sx2;
    PCA_sy=PCA_sy2;
    PCA_sz=PCA_sz2; 
    PCA_sxerr=PCA_sxerr2;
    PCA_syerr=PCA_syerr2;
    PCA_szerr=PCA_szerr2;
    PCA_tx=PCA_tx2;
    PCA_ty=PCA_ty2;
    PCA_tz=PCA_tz2; 
    PCA_txerr=PCA_txerr2;
    PCA_tyerr=PCA_tyerr2;
    PCA_tzerr=PCA_tzerr2;
    dotprodTrack=dotprodTrack2;
    dotprodSeed=dotprodSeed2;
    
//    std::cout << "filling   "<< pt << " " << eta << " " << phi << " " << dz << " " << dxy << " " << dist << " " << dsig << " " << std::endl;
//    std::cout << "filling   "<< PCA_sx << " " << PCA_sy << " " << PCA_sz << " " << PCA_tx << " " << PCA_ty << " " << PCA_tz << " "  << std::endl;
//    std::cout << "filling   "<< PCA_sxerr << " " << PCA_syerr << " " << PCA_szerr << " " << PCA_txerr << " " << PCA_tyerr << " " << PCA_tzerr << " "  << std::endl;
//        std::cout << "filling " << dotprodTrack << "  " << dotprodSeed << std::endl;
}

inline void trackVars2::set_vars ( double m, double t2dip, double t2dsip, double t3dip, double t3dsip, double t2dTS, double t3dTS, double t2dTSV, double t3dTSV){
mass=m;
t3Dip=t3dip;
t3Dsip=t3dsip;
t2Dip=t2dip;
t2Dsip=t2dsip;
dotprodTrackSeed2D=t2dTS;
dotprodTrackSeed3D=t3dTS;
dotprodTrackSeedVectors2D=t2dTSV;
dotprodTrackSeedVectors3D=t3dTSV;

//  std::cout << "filling  myTrack "<< std::endl;


}

// inline void trackVars2::set_signedIPs ( double a, double b, double c, double d){
// t3DipSigned=c;
// t3DsipSigned=d;
// t2DipSigned=a;
// t2DsipSigned=b;
// }


inline void trackVars2::set_index ( int a){
index=a;
}

inline void trackVars2::set_SeedIndex ( int a){
seed_index=a;
}


inline void trackVars2::set_distances ( double a, double b){
seedPCA_pv=a;
trackPCA_pv=b;
}


inline void trackVars2::set_JetAxisVars(double jadist, double dotprod, double d_eta, double d_phi){
    
//     std::cout<<one<<std::endl;
//     PCA_J=one;
    
    PCA_JetAxis_distance=jadist;
    PCAPair_Jet_dotprod=dotprod;
    PCAAxis_JetAxis_DEta=d_eta;
    PCAAxis_JetAxis_DPhi=d_phi; 
}

struct sortfunction2
{
    inline bool operator() (const trackVars2& struct1, const trackVars2& struct2)
    {
        return (struct1.dist < struct2.dist);
    }
};


class trackGenMatch2 {
      public:
    double chi_square;
    int numberOfDaughters;
    int MomFlav;
    int BChain;
    int GenIndex;
    int Status;
    
    void set_chi ( double );
    void set_numberOfDaughters ( int );
    void set_MomFlav ( int );
    void set_BChain ( int );
    void set_GenIndex ( int );
    void set_Status ( int );
};


inline void trackGenMatch2::set_chi ( double a){
chi_square=a;
}

inline void trackGenMatch2::set_numberOfDaughters ( int a){
numberOfDaughters=a;
}

inline void trackGenMatch2::set_MomFlav ( int a){
MomFlav=a;
}

inline void trackGenMatch2::set_BChain ( int a){
BChain=a;
}

inline void trackGenMatch2::set_GenIndex ( int a){
GenIndex=a;
}

inline void trackGenMatch2::set_Status ( int a){
Status=a;
}




struct sortgen2
{
    inline bool operator() (const trackGenMatch2& struct1, const trackGenMatch2& struct2)
    {
        return (struct1.chi_square < struct2.chi_square);
    }
};

#endif