#ifndef __STLAMBDAEVENT_H__
#define __STLAMBDAEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "StarClassLibrary/StThreeVectorF.hh"

// A. Schmah 05/31/2011



//------------------------------------------------------------------------------------
class StLambdaLambda : public TObject
{
private:
    // Track properties
    TLorentzVector TLV_Lambda; // Lorentz vector properties of this particle
    StThreeVectorF TV3_Lambda; // decay vertex of this particle

    StThreeVectorF TV3_proton_origin; //
    StThreeVectorF TV3_proton_gMom; //

    StThreeVectorF TV3_pion_origin; //
    StThreeVectorF TV3_pion_gMom; // 

    Float_t        dca_to_prim; // distance of closest approach to mother particle decay vertex (or primary vertex)
    Float_t        vertex_dist_to_prim; // distance of decay vertex to mother particle decay vertex (or primary vertex)
    Float_t        dca_daughters; // distance of of closest approach between daughters

    Float_t        Pi_dca;
    Float_t        Pi_m2;
    Float_t        Pi_nSigma;
    Float_t        Pi_qp;
    Int_t          Pi_Id;

    Float_t        P_dca;
    Float_t        P_m2;
    Float_t        P_nSigma;
    Float_t        P_qp;
    Int_t          P_Id;

public:
    StLambdaLambda() :
        TLV_Lambda(-1),TV3_Lambda(),TV3_proton_origin(),TV3_proton_gMom(),TV3_pion_origin(),TV3_pion_gMom(),dca_to_prim(-1),vertex_dist_to_prim(-1),dca_daughters(-1),Pi_dca(-1),Pi_m2(-1),Pi_nSigma(-1),Pi_qp(-1),Pi_Id(-1),
        P_dca(-1),P_m2(-1),P_nSigma(-1),P_qp(-1),P_Id(-1)
    {
    }
        ~StLambdaLambda() {}

        // setters
        void set_Pi_dca(Float_t f)                          { Pi_dca = f;                 }
        void set_Pi_m2(Float_t f)                           { Pi_m2 = f;                  }
        void set_Pi_nSigma(Float_t f)                       { Pi_nSigma = f;              }
        void set_Pi_qp(Float_t f)                           { Pi_qp = f;                  }
        void set_Pi_Id(Int_t i)                             { Pi_Id = i;                  }
        void set_P_dca(Float_t f)                           { P_dca = f;                  }
        void set_P_m2(Float_t f)                            { P_m2 = f;                   }
        void set_P_nSigma(Float_t f)                        { P_nSigma = f;               }
        void set_P_qp(Float_t f)                            { P_qp = f;                   }
        void set_P_Id(Int_t i)                              { P_Id = i;                   }
        void set_dca_to_prim(Float_t f)                     { dca_to_prim = f;            }
        void set_vertex_dist_to_prim(Float_t f)             { vertex_dist_to_prim = f;    }
        void set_dca_daughters(Float_t f)                   { dca_daughters = f;          }
        void set_TLV_Lambda(TLorentzVector tlv)             { TLV_Lambda = tlv;           }
        void set_TV3_Lambda(StThreeVectorF tv3)             { TV3_Lambda = tv3;           }
        void set_TV3_proton_origin(StThreeVectorF tv3)      { TV3_proton_origin = tv3;    }
        void set_TV3_proton_gMom(StThreeVectorF tv3)        { TV3_proton_gMom   = tv3;    }
        void set_TV3_pion_origin(StThreeVectorF tv3)        { TV3_pion_origin   = tv3;    }
        void set_TV3_pion_gMom(StThreeVectorF tv3)          { TV3_pion_gMom     = tv3;    }

        // getters
        Float_t get_Pi_dca()              const             { return Pi_dca;              }
        Float_t get_Pi_m2 ()              const             { return Pi_m2;               }
        Float_t get_Pi_nSigma()           const             { return Pi_nSigma;           }
        Float_t get_Pi_qp()               const             { return Pi_qp;               }
        Int_t   get_Pi_Id()               const             { return Pi_Id;               }
        Float_t get_P_dca()               const             { return P_dca;               }
        Float_t get_P_m2 ()               const             { return P_m2;                }
        Float_t get_P_nSigma()            const             { return P_nSigma;            }
        Float_t get_P_qp()                const             { return P_qp;                }
        Int_t   get_P_Id()                const             { return P_Id;                }
        Float_t get_dca_to_prim()         const             { return dca_to_prim;         }
        Float_t get_vertex_dist_to_prim() const             { return vertex_dist_to_prim; }
        Float_t get_dca_daughters()       const             { return dca_daughters;       }
        TLorentzVector get_TLV_Lambda()   const             { return TLV_Lambda;          }
        StThreeVectorF get_TV3_Lambda()   const             { return TV3_Lambda;          }
        StThreeVectorF get_TV3_proton_origin() const        { return TV3_proton_origin;   }
        StThreeVectorF get_TV3_proton_gMom()   const        { return TV3_proton_gMom;     }
        StThreeVectorF get_TV3_pion_origin()   const        { return TV3_pion_origin;     }
        StThreeVectorF get_TV3_pion_gMom()     const        { return TV3_pion_gMom;       }

        ClassDef(StLambdaLambda,1)  // A simple track of a particle
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class StLambdaEvent : public TObject
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

    UShort_t      fNumLambda;

    TClonesArray* fLambda;

public:
    StLambdaEvent() :
        x(-1),y(-1),z(-1),id(-1),mult(0),n_prim(0),n_non_prim(0),
        n_tof_prim(0),SE_ME_flag(-1),ZDCx(-1),BBCx(-1),vzVpd(-1),fNumLambda(0)
    {
        fLambda      = new TClonesArray( "StLambdaLambda", 10 );
    }
        ~StLambdaEvent()
        {
            delete fLambda;
            fLambda = NULL;
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


        //--------------------------------------
        // Lambda
        StLambdaLambda* createLambda()
        {
            if (fNumLambda == fLambda->GetSize())
                fLambda->Expand( fNumLambda + 5 );
            if (fNumLambda >= 200)
            {
                Fatal( "StLambdaMother::createLambda()", "ERROR: Too many Lambdas (>200)!" );
                exit( 2 );
            }

            new((*fLambda)[fNumLambda++]) StLambdaLambda;
            return (StLambdaLambda*)((*fLambda)[fNumLambda - 1]);
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
        StLambdaLambda* getLambda(UShort_t i) const
        {
            return i < fNumLambda ? (StLambdaLambda*)((*fLambda)[i]) : NULL;
        }
        //--------------------------------------


        ClassDef(StLambdaEvent,1)  // A simple event compiled of tracks
};
//------------------------------------------------------------------------------------



#endif // __STLAMBDAEVENT_H__
