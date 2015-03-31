#ifndef __STALEXEVENT_H__
#define __STALEXEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"

// A. Schmah 12.10.2011

class StAlexTrack : public TObject
{
private:
    // Track properties
    Float_t phi;
    Float_t theta;
    Float_t p_t;
    Float_t eta;
    Float_t dca;
    Float_t pq;
    Float_t m2;
    Float_t nsK;
    Float_t nsPi;
    Float_t nsP;
    Float_t nsEl;
    Float_t dEdx;
    Float_t iQx;
    Float_t iQy;
    Float_t phi_EP;
    Float_t phi_EP_raw;
    Float_t phi_EP_eta_gap;

public:
    StAlexTrack() :
        phi(-1),theta(-1),p_t(-1),eta(-1),dca(-1),pq(-1),m2(-1),nsK(-1),nsPi(-1),nsP(-1),nsEl(-1),dEdx(-1),iQx(-1),iQy(-1),
        phi_EP(-1),phi_EP_raw(-1),phi_EP_eta_gap(-1)
    {
    }
        ~StAlexTrack() {}

        void     setphi(Float_t d)           { phi = d;              }
        Float_t  getphi() const              { return phi;           }

        void     settheta(Float_t d)         { theta = d;            }
        Float_t  gettheta() const            { return theta;         }

        void     setp_t(Float_t d)           { p_t = d;              }
        Float_t  getp_t() const              { return p_t;           }

        void     seteta(Float_t d)           { eta = d;              }
        Float_t  geteta() const              { return eta;           }

        void     setdca(Float_t d)           { dca = d;              }
        Float_t  getdca() const              { return dca;           }

        void     setpq(Float_t d)            { pq = d;               }
        Float_t  getpq() const               { return pq;            }

        void     setm2(Float_t d)            { m2 = d;               }
        Float_t  getm2() const               { return m2;            }

        void     setnsK(Float_t d)           { nsK = d;              }
        Float_t  getnsK() const              { return nsK;           }

        void     setnsPi(Float_t d)          { nsPi = d;             }
        Float_t  getnsPi() const             { return nsPi;          }

        void     setnsP(Float_t d)           { nsP = d;              }
        Float_t  getnsP() const              { return nsP;           }

        void     setnsEl(Float_t d)          { nsEl = d;             }
        Float_t  getnsEl() const             { return nsEl;          }

        void     setdEdx(Float_t d)          { dEdx = d;             }
        Float_t  getdEdx() const             { return dEdx;          }

        void     setiQx(Float_t d)           { iQx = d;              }
        Float_t  getiQx() const              { return iQx;           }

        void     setiQy(Float_t d)           { iQy = d;              }
        Float_t  getiQy() const              { return iQy;           }

        void     setphi_EP(Float_t r)        { phi_EP = r;           }
        Float_t  getphi_EP() const           { return phi_EP;        }

        void     setphi_EP_raw(Float_t r)    { phi_EP_raw = r;       }
        Float_t  getphi_EP_raw() const       { return phi_EP_raw;    }

        void     setphi_EP_eta_gap(Float_t r)  { phi_EP_eta_gap = r;    }
        Float_t  getphi_EP_eta_gap() const     { return phi_EP_eta_gap; }

        ClassDef(StAlexTrack,1)  // A simple track of a particle
};

class StAlexEvent : public TObject
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
    StAlexEvent() :
        x(-1),y(-1),z(-1),id(-1),mult(0),n_prim(0),n_non_prim(0),
        n_tof_prim(0),EP_Qx_eta_pos_ptw(-1),EP_Qy_eta_pos_ptw(-1),EP_Qx_eta_neg_ptw(-1),EP_Qy_eta_neg_ptw(-1),
        EP_Qx_ptw(-1),EP_Qy_ptw(-1),Qtracks_eta_pos(0),Qtracks_eta_neg(0),Qtracks_full(0),ZDCx(-1),BBCx(-1),vzVpd(-1),fNumTracks(0)
    {
        fTracks      = new TClonesArray( "StAlexTrack", 10 );
    }
        ~StAlexEvent()
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

        StAlexTrack* createTrack()
        {
            if (fNumTracks == fTracks->GetSize())
                fTracks->Expand( fNumTracks + 10 );
            if (fNumTracks >= 10000)
            {
                Fatal( "StAlexEvent::createTrack()", "ERROR: Too many tracks (>1000)!" );
                exit( 2 );
            }

            new((*fTracks)[fNumTracks++]) StAlexTrack;
            return (StAlexTrack*)((*fTracks)[fNumTracks - 1]);
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
        StAlexTrack* getTrack(UShort_t i) const
        {
            return i < fNumTracks ? (StAlexTrack*)((*fTracks)[i]) : NULL;
        }

        ClassDef(StAlexEvent,1)  // A simple event compiled of tracks
};


#endif // __STALEXEVENT_H__
