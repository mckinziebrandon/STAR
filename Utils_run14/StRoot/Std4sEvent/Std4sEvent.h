#ifndef __STD4SEVENT_H__
#define __STD4SEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "StPhysicalHelixD.hh"
#include "StarClassLibrary/StThreeVectorF.hh"

// A. Schmah 05/31/2011
// A. Schmah 06/18/2013 added track structure



//------------------------------------------------------------------------------------
class Std4sTrack : public TObject
{
private:
    // Track properties
    StPhysicalHelixD THelix_Track; // helix of this particle
    Float_t          dca_to_prim;  // distance of closest approach to mother particle decay vertex (or primary vertex)
    Float_t          m2;           // squared mass
    Float_t          Pi_nSigma;    // pion nsgima dE/dx
    Float_t          K_nSigma;     // kaon nsigma dE/dx
    Float_t          P_nSigma;     // proton nsigma dE/dx
    Float_t          qp;           // charge times momentum
    Int_t            nhits;        // number of TPC points
    Int_t            trackId;      // trackId in the event
    Int_t            track_SE_ME;  // same event, mixed event flag

public:
    Std4sTrack() :
        THelix_Track(),dca_to_prim(-1),m2(-1),Pi_nSigma(-1),
        K_nSigma(-1),P_nSigma(-1),qp(-1),nhits(-1),trackId(-1),track_SE_ME(-1)
    {
    }
        ~Std4sTrack() {}

        // setters
        void set_THelix_Track(StPhysicalHelixD thel)        { THelix_Track = thel;        }
        void set_dca_to_prim(Float_t f)                     { dca_to_prim = f;            }
        void set_m2(Float_t f)                              { m2 = f;                     }
        void set_Pi_nSigma(Float_t f)                       { Pi_nSigma = f;              }
        void set_K_nSigma(Float_t f)                        { K_nSigma = f;               }
        void set_P_nSigma(Float_t f)                        { P_nSigma = f;               }
        void set_qp(Float_t f)                              { qp = f;                     }
        void set_nhits(Int_t i)                             { nhits = i;                  }
        void set_trackId(Int_t i)                           { trackId = i;                }
        void set_track_SE_ME(Int_t i)                       { track_SE_ME = i;            }


        // getters
        StPhysicalHelixD get_THelix_Track() const           { return THelix_Track;        }
        Float_t get_dca_to_prim()           const           { return dca_to_prim;         }
        Float_t get_m2 ()                   const           { return m2;                  }
        Float_t get_Pi_nSigma()             const           { return Pi_nSigma;           }
        Float_t get_K_nSigma()              const           { return K_nSigma;            }
        Float_t get_P_nSigma()              const           { return P_nSigma;            }
        Float_t get_qp()                    const           { return qp;                  }
        Int_t   get_nhits()                 const           { return nhits;               }
        Int_t   get_trackId()               const           { return trackId;             }
        Int_t   get_track_SE_ME()           const           { return track_SE_ME;         }


        ClassDef(Std4sTrack,1)  // A simple track of a particle
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class Std4sLambda : public TObject
{
private:
    // Track properties
    TLorentzVector TLV_Lambda; // Lorentz vector properties of this particle
    StThreeVectorF TV3_Lambda; // decay vertex of this particle
    StThreeVectorF TV3_Lambda_plane; // cross product of the two daughter particles at decay vertex

    Float_t        dca_to_prim; // distance of closest approach to mother particle decay vertex (or primary vertex)
    Float_t        vertex_dist_to_prim; // distance of decay vertex to mother particle decay vertex (or primary vertex)
    Float_t        dca_daughters; // distance of of closest approach between daughters

    UShort_t       fIndexPi;
    UShort_t       fIndexP;

public:
    Std4sLambda() :
        TLV_Lambda(-1),TV3_Lambda(),TV3_Lambda_plane(),dca_to_prim(-1),vertex_dist_to_prim(-1),dca_daughters(-1),
        fIndexPi(-1),fIndexP(-1)
    {
    }
        ~Std4sLambda() {}

        // setters
        void set_dca_to_prim(Float_t f)                     { dca_to_prim = f;            }
        void set_vertex_dist_to_prim(Float_t f)             { vertex_dist_to_prim = f;    }
        void set_dca_daughters(Float_t f)                   { dca_daughters = f;          }
        void set_TLV_Lambda(TLorentzVector tlv)             { TLV_Lambda = tlv;           }
        void set_TV3_Lambda(StThreeVectorF tv3)             { TV3_Lambda = tv3;           }
        void set_TV3_Lambda_plane(StThreeVectorF tv3)       { TV3_Lambda_plane = tv3;     }
        void set_IndexPi(Int_t i)                           { fIndexPi = i;               }
        void set_IndexP(Int_t i)                            { fIndexP = i;                }

        // getters
        Float_t get_dca_to_prim()             const         { return dca_to_prim;         }
        Float_t get_vertex_dist_to_prim()     const         { return vertex_dist_to_prim; }
        Float_t get_dca_daughters()           const         { return dca_daughters;       }
        TLorentzVector get_TLV_Lambda()       const         { return TLV_Lambda;          }
        StThreeVectorF get_TV3_Lambda()       const         { return TV3_Lambda;          }
        StThreeVectorF get_TV3_Lambda_plane() const         { return TV3_Lambda_plane;    }
        Int_t   get_IndexPi()                 const         { return fIndexPi;            }
        Int_t   get_IndexP()                  const         { return fIndexP;             }

        ClassDef(Std4sLambda,1)  // A simple track of a particle
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class Std4sXi : public TObject
{
private:
    // Xi properties
    StThreeVectorF TV3_Xi; // decay vertex of this particle
    StPhysicalHelixD THelix_Xi; // helix of this particle

    Float_t        dca_to_prim; // distance of closest approach to mother particle decay vertex (or primary vertex)
    Float_t        vertex_dist_to_prim; // distance of decay vertex to mother particle decay vertex (or primary vertex)
    Float_t        dca_daughters; // distance of of closest approach between daughters

    UShort_t       fIndexLambdaXi;
    UShort_t       fIndexXiPi;

public:
    Std4sXi() :
        TV3_Xi(),THelix_Xi(),dca_to_prim(-1),vertex_dist_to_prim(-1),dca_daughters(-1),fIndexLambdaXi(-1),fIndexXiPi()
    {
    }
        ~Std4sXi() {}

        // setters
        void set_dca_to_prim(Float_t f)                     { dca_to_prim = f;            }
        void set_vertex_dist_to_prim(Float_t f)             { vertex_dist_to_prim = f;    }
        void set_dca_daughters(Float_t f)                   { dca_daughters = f;          }
        void set_IndexLambdaXi(Int_t i)                     { fIndexLambdaXi = i;         }
        void set_IndexXiPi(Int_t i)                         { fIndexXiPi = i;             }
        void set_TV3_Xi(StThreeVectorF tv3)                 { TV3_Xi = tv3;               }
        void set_THelix_Xi(StPhysicalHelixD thel)           { THelix_Xi = thel;           }

        // getters
        Float_t get_dca_to_prim()         const             { return dca_to_prim;         }
        Float_t get_vertex_dist_to_prim() const             { return vertex_dist_to_prim; }
        Float_t get_dca_daughters()       const             { return dca_daughters;       }
        Int_t   get_IndexLambdaXi()       const             { return fIndexLambdaXi;      }
        Int_t   get_IndexXiPi()           const             { return fIndexXiPi;          }
        StThreeVectorF get_TV3_Xi()       const             { return TV3_Xi;              }
        StPhysicalHelixD get_THelix_Xi()  const             { return THelix_Xi;           }


        ClassDef(Std4sXi,1)  // A simple track of a particle
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class Std4sMother : public TObject
{
private:
    // Track properties
    TLorentzVector TLV_d4s; // Lorentz vector properties of this particle
    TLorentzVector TLV_Xi;  // Lorentz vector properties of this particle
    StThreeVectorF TV3_d4s; // decay vertex of this particle
    StPhysicalHelixD THelix_d4s; // helix of this particle

    Float_t        dca_to_prim; // distance of closest approach to mother particle decay vertex (or primary vertex)
    Float_t        vertex_dist_to_prim; // distance of decay vertex to mother particle decay vertex (or primary vertex)
    Float_t        dca_daughters; // distance of of closest approach between daughters
    Float_t        scalar_product; // scalar product of momentum vector and vector between decay vertex and primary vertex
    Float_t        angle_Lambdas; // angle between Lambda decay planes

    UShort_t       fIndexLambda;
    UShort_t       fIndexXi;

public:
    Std4sMother() :
        TLV_d4s(-1),TLV_Xi(-1),TV3_d4s(),THelix_d4s(),dca_to_prim(-1),vertex_dist_to_prim(-1),dca_daughters(-1),scalar_product(-1),
        angle_Lambdas(-1),fIndexLambda(-1),fIndexXi(-1)
    {
    }
        ~Std4sMother() {}

        // setters
        void set_dca_to_prim(Float_t f)                     { dca_to_prim = f;            }
        void set_vertex_dist_to_prim(Float_t f)             { vertex_dist_to_prim = f;    }
        void set_dca_daughters(Float_t f)                   { dca_daughters = f;          }
        void set_scalar_product(Float_t f)                  { scalar_product = f;         }
        void set_angle_Lambdas(Float_t f)                   { angle_Lambdas = f;          }
        void set_IndexLambda(Int_t i)                       { fIndexLambda = i;           }
        void set_IndexXi(Int_t i)                           { fIndexXi = i;               }
        void set_TLV_d4s(TLorentzVector tlv)                { TLV_d4s = tlv;              }
        void set_TLV_Xi(TLorentzVector tlv)                 { TLV_Xi = tlv;               }
        void set_TV3_d4s(StThreeVectorF tv3)                { TV3_d4s = tv3;              }
        void set_THelix_d4s(StPhysicalHelixD thel)          { THelix_d4s = thel;          }

        // getters
        Float_t get_dca_to_prim()         const             { return dca_to_prim;         }
        Float_t get_vertex_dist_to_prim() const             { return vertex_dist_to_prim; }
        Float_t get_dca_daughters()       const             { return dca_daughters;       }
        Float_t get_scalar_product()      const             { return scalar_product;      }
        Float_t get_angle_Lambdas()       const             { return angle_Lambdas;       }
        Int_t   get_IndexLambda()         const             { return fIndexLambda;        }
        Int_t   get_IndexXi()             const             { return fIndexXi;            }
        TLorentzVector get_TLV_d4s()      const             { return TLV_d4s;             }
        TLorentzVector get_TLV_Xi()       const             { return TLV_Xi;              }
        StThreeVectorF get_TV3_d4s()      const             { return TV3_d4s;             }
        StPhysicalHelixD get_THelix_d4s() const             { return THelix_d4s;          }


        ClassDef(Std4sMother,1)  // A simple track of a particle
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class Std4sEvent : public TObject
{
private:
    Float_t x;
    Float_t y;
    Float_t z;
    Int_t   id;
    Float_t mult;
    Int_t   n_prim;
    Int_t   n_non_prim;
    Int_t   n_tof_prim;
    Int_t   SE_ME_flag;

    Float_t ZDCx;
    Float_t BBCx;
    Float_t vzVpd;
    Float_t Psi2;
    Float_t Bfield;

    UShort_t      fNumd4s;
    UShort_t      fNumXi;
    UShort_t      fNumLambda;
    UShort_t      fNumTrack;

    TClonesArray* fd4s;
    TClonesArray* fXi;
    TClonesArray* fLambda;
    TClonesArray* fTrack;

public:
    Std4sEvent() :
        x(-1),y(-1),z(-1),id(-1),mult(0),n_prim(0),n_non_prim(0),
        n_tof_prim(0),SE_ME_flag(-1),ZDCx(-1),BBCx(-1),vzVpd(-1),Psi2(-1),Bfield(-1),fNumd4s(0),fNumXi(0),fNumLambda(0),fNumTrack(0)
    {
        fd4s         = new TClonesArray( "Std4sMother", 10 );
        fXi          = new TClonesArray( "Std4sXi", 10 );
        fLambda      = new TClonesArray( "Std4sLambda", 10 );
        fTrack       = new TClonesArray( "Std4sTrack", 10 );
    }
        ~Std4sEvent()
        {
            delete fd4s;
            fd4s = NULL;

            delete fXi;
            fXi = NULL;

            delete fLambda;
            fLambda = NULL;

            delete fTrack;
            fTrack = NULL;
        }

        void       setx(Float_t r)                    { x = r;                         }
        Float_t    getx() const                       { return x;                      }

        void       sety(Float_t r)                    { y = r;                         }
        Float_t    gety() const                       { return y;                      }

        void       setz(Float_t r)                    { z = r;                         }
        Float_t    getz() const                       { return z;                      }

        void       setid(Int_t  r)                    { id = r;                        }
        Int_t      getid() const                      { return id;                     }

        void       setmult(Float_t r)                 { mult = r;                      }
        Float_t    getmult() const                    { return mult;                   }

        void       setn_prim(Int_t r)                 { n_prim = r;                    }
        Int_t      getn_prim() const                  { return n_prim;                 }

        void       setn_non_prim(Int_t r)             { n_non_prim = r;                }
        Int_t      getn_non_prim() const              { return n_non_prim;             }

        void       setn_tof_prim(Int_t r)             { n_tof_prim = r;                }
        Int_t      getn_tof_prim() const              { return n_tof_prim;             }

        void       setSE_ME_flag(Int_t r)             { SE_ME_flag = r;                }
        Int_t      getSE_ME_flag() const              { return SE_ME_flag;             }

        void       setZDCx(Float_t r)                 { ZDCx = r;                      }
        Float_t    getZDCx() const                    { return ZDCx;                   }

        void       setBBCx(Float_t r)                 { BBCx = r;                      }
        Float_t    getBBCx() const                    { return BBCx;                   }

        void       setvzVpd(Float_t r)                { vzVpd = r;                     }
        Float_t    getvzVpd() const                   { return vzVpd;                  }

        void       setPsi2(Float_t r)                 { Psi2 = r;                      }
        Float_t    getPsi2() const                    { return Psi2;                   }

        void       setBfield(Float_t r)               { Bfield = r;                    }
        Float_t    getBfield() const                  { return Bfield;                 }


        //--------------------------------------
        // d4s
        Std4sMother* created4s()
        {
            if (fNumd4s == fd4s->GetSize())
                fd4s->Expand( fNumd4s + 10 );
            if (fNumd4s >= 200)
            {
                Fatal( "Std4sEvent::created4s()", "ERROR: Too many tracks (>200)!" );
                exit( 2 );
            }

            new((*fd4s)[fNumd4s++]) Std4sMother;
            return (Std4sMother*)((*fd4s)[fNumd4s - 1]);
        }
        void cleard4sList()
        {
            fNumd4s   = 0;
            fd4s      ->Clear();
        }
        UShort_t getNumd4s() const
        {
            return fNumd4s;
        }
        Std4sMother* getd4s(UShort_t i) const
        {
            return i < fNumd4s ? (Std4sMother*)((*fd4s)[i]) : NULL;
        }
        //--------------------------------------


        //--------------------------------------
        // Xi
        Std4sXi* createXi()
        {
            if (fNumXi == fXi->GetSize())
                fXi->Expand( fNumXi + 5 );
            if (fNumXi >= 200)
            {
                Fatal( "Std4sMother::createXi()", "ERROR: Too many Xis (>200)!" );
                exit( 2 );
            }

            new((*fXi)[fNumXi++]) Std4sXi;
            return (Std4sXi*)((*fXi)[fNumXi - 1]);
        }
        void clearXiList()
        {
            fNumXi   = 0;
            fXi      ->Clear();
        }
        UShort_t getNumXi() const
        {
            return fNumXi;
        }
        Std4sXi* getXi(UShort_t i) const
        {
            return i < fNumXi ? (Std4sXi*)((*fXi)[i]) : NULL;
        }
        //--------------------------------------


        //--------------------------------------
        // Lambda
        Std4sLambda* createLambda()
        {
            if (fNumLambda == fLambda->GetSize())
                fLambda->Expand( fNumLambda + 5 );
            if (fNumLambda >= 200)
            {
                Fatal( "Std4sMother::createLambda()", "ERROR: Too many Lambdas (>200)!" );
                exit( 2 );
            }

            new((*fLambda)[fNumLambda++]) Std4sLambda;
            return (Std4sLambda*)((*fLambda)[fNumLambda - 1]);
        }
        void clearLambdaList()
        {
            fNumLambda   = 0;
            fLambda      ->Clear();
        }
        UShort_t getNumLambda() const
        {
            return fNumLambda;
        }
        Std4sLambda* getLambda(UShort_t i) const
        {
            return i < fNumLambda ? (Std4sLambda*)((*fLambda)[i]) : NULL;
        }
        //--------------------------------------


        //--------------------------------------
        // Track
        Std4sTrack* createTrack()
        {
            if (fNumTrack == fTrack->GetSize())
                fTrack->Expand( fNumTrack + 5 );
            if (fNumTrack >= 2000)
            {
                Fatal( "Std4sMother::createTrack()", "ERROR: Too many Tracks (>2000)!" );
                exit( 2 );
            }

            new((*fTrack)[fNumTrack++]) Std4sTrack;
            return (Std4sTrack*)((*fTrack)[fNumTrack - 1]);
        }
        void clearTrackList()
        {
            fNumTrack   = 0;
            fTrack      ->Clear();
        }
        UShort_t getNumTrack() const
        {
            return fNumTrack;
        }
        Std4sTrack* getTrack(UShort_t i) const
        {
            return i < fNumTrack ? (Std4sTrack*)((*fTrack)[i]) : NULL;
        }
        //--------------------------------------


        //--------------------------------------
        Bool_t isBaryon(UShort_t iLambda) const
        {
            if(iLambda < fNumLambda)
            {
                Std4sLambda* Lambda = (Std4sLambda*)((*fLambda)[iLambda]);
                Int_t iP = Lambda->get_IndexP();
                if(iP < fNumTrack)
                {
                    Std4sTrack* Track = (Std4sTrack*)((*fTrack)[iP]);
                    Float_t qp = Track ->get_qp();
                    if(qp > 0.0)
                    {
                        return 1;
                    }
                    else return 0;
                }
                else cout << "WARNING: number of Tracks exceeded in isBaryon!" << endl;
            }
            else cout << "WARNING: number of Lambdas exceeded in isBaryon!" << endl;
        }
        //--------------------------------------


        //--------------------------------------
        // calculate decay plane angle of Lambda candidate
        StThreeVectorF calc_decay_plane_vec(UShort_t iLambda) const
        {
            if(iLambda < fNumLambda)
            {
                Std4sLambda* Lambda = (Std4sLambda*)((*fLambda)[iLambda]);
                Int_t iP  = Lambda->get_IndexP();
                Int_t iPi = Lambda->get_IndexPi();
                if(iP < fNumTrack && iPi < fNumTrack)
                {
                    Std4sTrack* Track_P  = (Std4sTrack*)((*fTrack)[iP]);
                    Std4sTrack* Track_Pi = (Std4sTrack*)((*fTrack)[iPi]);

                    StPhysicalHelixD helix_P  = Track_P  ->get_THelix_Track();
                    StPhysicalHelixD helix_Pi = Track_Pi ->get_THelix_Track();
                    pair <Double_t, Double_t> pair_lengths;
                    pair_lengths = helix_P.pathLengths(helix_Pi);
                    StThreeVectorF dir_P  = helix_P.cat(pair_lengths.first);
                    StThreeVectorF dir_Pi = helix_Pi.cat(pair_lengths.second);
                    StThreeVectorF Stv3_perpPPi = dir_P.cross(dir_Pi);
                    return Stv3_perpPPi;
                }
                else cout << "WARNING: number of Tracks exceeded in isBaryon!" << endl;
            }
            else cout << "WARNING: number of Lambdas exceeded in isBaryon!" << endl;
        }
        //--------------------------------------


        ClassDef(Std4sEvent,1)  // A simple event compiled of tracks
};
//------------------------------------------------------------------------------------



#endif // __STD4SEVENT_H__
