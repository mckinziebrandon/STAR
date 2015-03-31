#ifndef __STD0EVENT_H__
#define __STD0EVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"

// A. Schmah 12.03.2015

class StD0Track : public TObject
{
private:
    // Track properties
    Float_t m2A; // mass2 of particle A
    Float_t m2B;
    Float_t nsA; // nsigma dE/dx of particle A
    Float_t nsB;
    Float_t dcaA; // distance of closest approach of particle A
    Float_t dcaB;
    Float_t iQxA; // Q-vector x-component of particle A
    Float_t iQyA;
    Float_t iQxB;
    Float_t iQyB;
    Float_t etaA;  // pseudo rapidity of particle A
    Float_t etaB;

    Float_t InvAB;  // invariant mass of particles A and B
    Float_t pt; // transverse momentum of particles AB(C)
    Float_t rap; // rapidity of particles AB(C)
    Float_t phi; // azimuth angle of particles AB(C)
    Float_t theta; // polar angle of particles AB(C)

    Float_t qpA; // momentum times charge
    Float_t qpB;

    Float_t VerdistX;
    Float_t VerdistY;
    Float_t dcaAB;

public:
    StD0Track() :
        m2A(-1),m2B(-1),nsA(-1),nsB(-1),dcaA(-1),dcaB(-1),iQxA(-1),iQyA(-1),iQxB(-1),iQyB(-1),
        etaA(-1),etaB(-1),InvAB(-1),pt(-1),
        rap(-1),phi(-1),theta(-1),qpA(-1),qpB(-1),
        VerdistX(-1),VerdistY(-1),dcaAB(-1)
    {
    }
        ~StD0Track() {}

        // setters
        void setm2A(Float_t f)                     { m2A = f;         }
        void setm2B(Float_t f)                     { m2B = f;         }
        void setnsA(Float_t f)                     { nsA = f;         }
        void setnsB(Float_t f)                     { nsB = f;         }
        void setdcaA(Float_t f)                    { dcaA = f;        }
        void setdcaB(Float_t f)                    { dcaB = f;        }
        void setiQxA(Float_t f)                    { iQxA = f;        }
        void setiQyA(Float_t f)                    { iQyA = f;        }
        void setiQxB(Float_t f)                    { iQxB = f;        }
        void setiQyB(Float_t f)                    { iQyB = f;        }
        void setetaA(Float_t f)                    { etaA = f;        }
        void setetaB(Float_t f)                    { etaB = f;        }
        void setInvAB(Float_t f)                   { InvAB = f;       }
        void setpt(Float_t f)                      { pt = f;          }
        void setrap(Float_t f)                     { rap = f;         }
        void setphi(Float_t f)                     { phi = f;         }
        void settheta(Float_t f)                   { theta = f;       }
        void setqpA(Float_t f)                     { qpA = f;         }
        void setqpB(Float_t f)                     { qpB = f;         }
        void setVerdistX(Float_t f)                { VerdistX = f;    }
        void setVerdistY(Float_t f)                { VerdistY = f;    }
        void setdcaAB(Float_t f)                   { dcaAB = f;       }


        // getters
        Float_t getm2A() const                     { return m2A;         }
        Float_t getm2B() const                     { return m2B;         }
        Float_t getnsA() const                     { return nsA;         }
        Float_t getnsB() const                     { return nsB;         }
        Float_t getdcaA() const                    { return dcaA;        }
        Float_t getdcaB() const                    { return dcaB;        }
        Float_t getiQxA() const                    { return iQxA;        }
        Float_t getiQyA() const                    { return iQyA;        }
        Float_t getiQxB() const                    { return iQxB;        }
        Float_t getiQyB() const                    { return iQyB;        }
        Float_t getetaA() const                    { return etaA;        }
        Float_t getetaB() const                    { return etaB;        }
        Float_t getInvAB() const                   { return InvAB;       }
        Float_t getpt() const                      { return pt;          }
        Float_t getrap() const                     { return rap;         }
        Float_t getphi() const                     { return phi;         }
        Float_t gettheta() const                   { return theta;       }
        Float_t getqpA() const                     { return qpA;         }
        Float_t getqpB() const                     { return qpB;         }
        Float_t getVerdistX() const                { return VerdistX;    }
        Float_t getVerdistY() const                { return VerdistY;    }
        Float_t getdcaAB() const                   { return dcaAB;       }


        ClassDef(StD0Track,1)  // A simple track of a particle
};

class StD0Event : public TObject
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
    Float_t EP_Qx_eta_pos_ptw;
    Float_t EP_Qy_eta_pos_ptw;
    Float_t EP_Qx_eta_neg_ptw;
    Float_t EP_Qy_eta_neg_ptw;
    Float_t EP_Qx_ptw;
    Float_t EP_Qy_ptw;
    Int_t   Qtracks_eta_pos;
    Int_t   Qtracks_eta_neg;
    Int_t   Qtracks_full;

    Float_t ZDCx;
    Float_t BBCx;
    Float_t vzVpd;

    UShort_t      fNumTracks;

    TClonesArray* fTracks;      //->

public:
    StD0Event() :
        x(-1),y(-1),z(-1),id(-1),mult(0),n_prim(0),n_non_prim(0),
        n_tof_prim(0),EP_Qx_eta_pos_ptw(-1),EP_Qy_eta_pos_ptw(-1),EP_Qx_eta_neg_ptw(-1),EP_Qy_eta_neg_ptw(-1),
        EP_Qx_ptw(-1),EP_Qy_ptw(-1),Qtracks_eta_pos(0),Qtracks_eta_neg(0),Qtracks_full(0),ZDCx(-1),BBCx(-1),vzVpd(-1),fNumTracks(0)
    {
        fTracks      = new TClonesArray( "StD0Track", 10 );
    }
        ~StD0Event()
        {
            delete fTracks;
            fTracks = NULL;
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

        void       setEP_Qx_eta_pos_ptw(Float_t r)    { EP_Qx_eta_pos_ptw = r;         }
        Float_t    getEP_Qx_eta_pos_ptw() const       { return EP_Qx_eta_pos_ptw;      }

        void       setEP_Qy_eta_pos_ptw(Float_t r)    { EP_Qy_eta_pos_ptw = r;         }
        Float_t    getEP_Qy_eta_pos_ptw() const       { return EP_Qy_eta_pos_ptw;      }

        void       setEP_Qx_eta_neg_ptw(Float_t r)    { EP_Qx_eta_neg_ptw = r;         }
        Float_t    getEP_Qx_eta_neg_ptw() const       { return EP_Qx_eta_neg_ptw;      }

        void       setEP_Qy_eta_neg_ptw(Float_t r)    { EP_Qy_eta_neg_ptw = r;         }
        Float_t    getEP_Qy_eta_neg_ptw() const       { return EP_Qy_eta_neg_ptw;      }

        void       setEP_Qx_ptw(Float_t r)            { EP_Qx_ptw = r;                 }
        Float_t    getEP_Qx_ptw() const               { return EP_Qx_ptw;              }

        void       setEP_Qy_ptw(Float_t r)            { EP_Qy_ptw = r;                 }
        Float_t    getEP_Qy_ptw() const               { return EP_Qy_ptw;              }

        void       setQtracks_eta_pos(Int_t r)        { Qtracks_eta_pos = r;           }
        Int_t      getQtracks_eta_pos() const         { return Qtracks_eta_pos;        }

        void       setQtracks_eta_neg(Int_t r)        { Qtracks_eta_neg = r;           }
        Int_t      getQtracks_eta_neg() const         { return Qtracks_eta_neg;        }

        void       setQtracks_full(Int_t r)           { Qtracks_full = r;              }
        Int_t      getQtracks_full() const            { return Qtracks_full;           }

        void       setZDCx(Float_t r)                 { ZDCx = r;                      }
        Float_t    getZDCx() const                    { return ZDCx;                   }

        void       setBBCx(Float_t r)                 { BBCx = r;                      }
        Float_t    getBBCx() const                    { return BBCx;                   }

        void       setvzVpd(Float_t r)                { vzVpd = r;                     }
        Float_t    getvzVpd() const                   { return vzVpd;                  }

        StD0Track* createTrack()
        {
            if (fNumTracks == fTracks->GetSize())
                fTracks->Expand( fNumTracks + 10 );
            if (fNumTracks >= 10000)
            {
                Fatal( "StD0Event::createTrack()", "ERROR: Too many tracks (>1000)!" );
                exit( 2 );
            }

            new((*fTracks)[fNumTracks++]) StD0Track;
            return (StD0Track*)((*fTracks)[fNumTracks - 1]);
        }
        void clearTrackList()
        {
            fNumTracks   = 0;
            fTracks      ->Clear();
        }
        UShort_t getNumTracks() const
        {
            return fNumTracks;
        }
        StD0Track* getTrack(UShort_t i) const
        {
            return i < fNumTracks ? (StD0Track*)((*fTracks)[i]) : NULL;
        }

        ClassDef(StD0Event,1)  // A simple event compiled of tracks
};


#endif // __STD0EVENT_H__
