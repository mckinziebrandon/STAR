#ifndef __STPARTICLEEVENT_H__
#define __STPARTICLEEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TVector2.h"
//#include "TLorentzVector.h"

// Modified by A. Schmah 07.04.2009: Sim info added
// Modified by A. Schmah 17.04.2009: Forward wall info added
// Modified by A. Schmah 18.06.2009: New structure with Sim track and WallHit
// Modified by A. Schmah 05.04.2010: Adapted to STAR structure
// Modified by A. Schmah 23.04.2010: Added TriggerId
// Modified by P. Huck   13.08.2010: Added Primary Tracks and nHitsdEdx

//class StParticleTrack : public TLorentzVector
class StParticleTrack : public TObject
{
private:

    // Track properties
    Short_t fCharge;            // particle's charge (+1,-1)
    Float_t fBeta;              // particle's beta
    Float_t fMomentum;          // particle's global momentum [MeV]
    Float_t fMass2;             // particle's mass^2 [MeV^2]
    Float_t fTPCdEdx;           // Mdc dE/dx for inner and outer segment

    // TPC dE/dx PID
    Float_t fnSigmaEl;
    Float_t fnSigmaPion;
    Float_t fnSigmaKaon;
    Float_t fnSigmaP;

    Float_t fDCA;  // distance of closest approach to primary vertex
    Short_t fnHitsPoss;
    Short_t fnHitsFit;

    // STAR helix information
    Float_t  fdipAngle;
    Float_t  fcurvature;
    Float_t  fphase;
    Int_t    fh;
    Float_t  foriginX;
    Float_t  foriginY;
    Float_t  foriginZ;
    Short_t fnHitsdEdx;
    Float_t  fdipAnglePrim;
    Float_t  fcurvaturePrim;
    Float_t  fphasePrim;
    Int_t    fhPrim;
    Float_t  foriginXPrim;
    Float_t  foriginYPrim;
    Float_t  foriginZPrim;
    Float_t fMomentumPrim;      // particle's primary momentum [MeV]
    Float_t  fTofHitX;
    Float_t  fTofHitY;
    Float_t  fTofHitZ;
    Float_t  fTof;
    Float_t  fPathLength;
    Float_t  fChi2Prob;
    Float_t  fyLocal;
    TVector2 fQVecFullTPC;
    TVector2 fQVecEtaPos;
    TVector2 fQVecEtaNeg;
    Bool_t fFlagFullTPC;
    Bool_t fFlagEtaPos;
    Bool_t fFlagEtaNeg;

public:
    StParticleTrack() :
        fCharge(0),fBeta(-1),fMomentum(-1),fMass2(-1),
        fTPCdEdx(-1),fnSigmaEl(-1),fnSigmaPion(-1),fnSigmaKaon(-1),fnSigmaP(-1), fDCA(-1),
        fnHitsPoss(0),fnHitsFit(0),fdipAngle(-1),fcurvature(-1),fphase(-1),fh(-1),foriginX(-1),
        foriginY(-1),foriginZ(-1),fnHitsdEdx(0),fdipAnglePrim(-1),fcurvaturePrim(-1),fphasePrim(-1),fhPrim(-1),foriginXPrim(-1),
	foriginYPrim(-1),foriginZPrim(-1),fMomentumPrim(-1),fTofHitX(-1),fTofHitY(-1),fTofHitZ(-1),fTof(-1),fPathLength(-1),fChi2Prob(-1),fyLocal(-1),
	fQVecFullTPC(-999,-999), fQVecEtaPos(-999,-999), fQVecEtaNeg(-999,-999),
	fFlagFullTPC(0), fFlagEtaPos(0), fFlagEtaNeg(0)
    {
    }
        ~StParticleTrack() {}

        void     setTPCdEdx(Float_t d)           { fTPCdEdx = d;              }
        Float_t  getTPCdEdx() const              { return fTPCdEdx;           }

        void     setCharge(Short_t c)            { fCharge = c;               }
        Short_t  getCharge() const               { return fCharge;            }

        void     setBeta(Float_t b)              { fBeta = b;                 }
        Float_t  getBeta() const                 { return fBeta;              }

        void     setMomentum(Float_t m)          { fMomentum = m;             }
        Float_t  getMomentum() const             { return fMomentum;          }

        void     setMass2(Float_t m)             { fMass2 = m;                }
        Float_t  getMass2() const                { return fMass2;             }

        void     setnSigmaEl(Float_t a)          { fnSigmaEl = a;             }
        Float_t  getnSigmaEl() const             { return fnSigmaEl;          }

        void     setnSigmaPion(Float_t a)        { fnSigmaPion = a;           }
        Float_t  getnSigmaPion() const           { return fnSigmaPion;        }

        void     setnSigmaKaon(Float_t a)        { fnSigmaKaon = a;           }
        Float_t  getnSigmaKaon() const           { return fnSigmaKaon;        }

        void     setnSigmaP(Float_t a)           { fnSigmaP = a;              }
        Float_t  getnSigmaP() const              { return fnSigmaP;           }

        void     setDCA(Float_t a)               { fDCA = a;                  }
        Float_t  getDCA() const                  { return fDCA;               }

        void     setnHitsPoss(Short_t a)         { fnHitsPoss = a;            }
        Short_t  getnHitsPoss() const            { return fnHitsPoss;         }

        void     setnHitsFit(Short_t a)          { fnHitsFit = a;             }
        Short_t  getnHitsFit() const             { return fnHitsFit;          }

        void     setnHitsdEdx(Short_t a)          { fnHitsdEdx = a;           }
        Short_t  getnHitsdEdx() const             { return fnHitsdEdx;        }

        void     setdipAngle(Float_t a)          { fdipAngle = a;             }
        Float_t  getdipAngle() const             { return fdipAngle;          }
 
        void     setcurvature(Float_t a)         { fcurvature = a;            }
        Float_t  getcurvature() const            { return fcurvature;         }

        void     setphase(Float_t a)             { fphase = a;                }
        Float_t  getphase() const                { return fphase;             }

        void     seth(Int_t a)                   { fh = a;                    }
        Int_t    geth() const                    { return fh;                 }

        void     setoriginX(Float_t a)           { foriginX = a;              }
        Float_t  getoriginX() const              { return foriginX;           }

        void     setoriginY(Float_t a)           { foriginY = a;              }
        Float_t  getoriginY() const              { return foriginY;           }

        void     setoriginZ(Float_t a)           { foriginZ = a;              }
        Float_t  getoriginZ() const              { return foriginZ;           }

        void     setdipAnglePrim(Float_t a)      { fdipAnglePrim = a;         }
        Float_t  getdipAnglePrim() const         { return fdipAnglePrim;      }
 
        void     setcurvaturePrim(Float_t a)     { fcurvaturePrim = a;        }
        Float_t  getcurvaturePrim() const        { return fcurvaturePrim;     }

        void     setphasePrim(Float_t a)         { fphasePrim = a;            }
        Float_t  getphasePrim() const            { return fphasePrim;         }

        void     sethPrim(Int_t a)               { fhPrim = a;                }
        Int_t    gethPrim() const                { return fhPrim;             }

        void     setoriginXPrim(Float_t a)       { foriginXPrim = a;          }
        Float_t  getoriginXPrim() const          { return foriginXPrim;       }

        void     setoriginYPrim(Float_t a)       { foriginYPrim = a;          }
        Float_t  getoriginYPrim() const          { return foriginYPrim;       }

        void     setoriginZPrim(Float_t a)       { foriginZPrim = a;          }
        Float_t  getoriginZPrim() const          { return foriginZPrim;       }

        void     setMomentumPrim(Float_t m)      { fMomentumPrim = m;         }
        Float_t  getMomentumPrim() const         { return fMomentumPrim;      }

        void     setTofHitX(Float_t a)           { fTofHitX = a;              }
        Float_t  getTofHitX() const              { return fTofHitX;           }

        void     setTofHitY(Float_t a)           { fTofHitY = a;              }
        Float_t  getTofHitY() const              { return fTofHitY;           }

        void     setTofHitZ(Float_t a)           { fTofHitZ = a;              }
        Float_t  getTofHitZ() const              { return fTofHitZ;           }

        void     setTof(Float_t a)           { fTof = a;              }
        Float_t  getTof() const              { return fTof;           }

        void     setPathLength(Float_t a)           { fPathLength = a;              }
        Float_t  getPathLength() const              { return fPathLength;           }

        void     setChi2Prob(Float_t a)           { fChi2Prob = a;              }
        Float_t  getChi2Prob() const              { return fChi2Prob;           }

        void     setYLocal(Float_t a)           { fyLocal = a;              }
        Float_t  getYLocal() const              { return fyLocal;           }

        void     setQVecFullTPC(TVector2* vec)   { fQVecFullTPC.Set(vec->X(),vec->Y());        }
        const TVector2* getQVecFullTPC() const         { return &fQVecFullTPC;     }

        void     setQVecEtaPos(TVector2* vec)   { fQVecEtaPos.Set(vec->X(),vec->Y());        }
        const TVector2* getQVecEtaPos() const         { return &fQVecEtaPos;     }

        void     setQVecEtaNeg(TVector2* vec)   { fQVecEtaNeg.Set(vec->X(),vec->Y());        }
        const TVector2* getQVecEtaNeg() const         { return &fQVecEtaNeg;     }

        void    setFlagFullTPC(Bool_t a)       { fFlagFullTPC = a;              }
        Bool_t  getFlagFullTPC() const          { return fFlagFullTPC;           }

        void    setFlagEtaPos(Bool_t a)       { fFlagEtaPos = a;              }
        Bool_t  getFlagEtaPos() const          { return fFlagEtaPos;           }

        void    setFlagEtaNeg(Bool_t a)       { fFlagEtaNeg = a;              }
        Bool_t  getFlagEtaNeg() const          { return fFlagEtaNeg;           }
        ClassDef(StParticleTrack,2)  // A simple track of a particle
};

class StParticleEvent : public TObject
{
private:
    Int_t  fRunId;
    Int_t  fEventId;
    Int_t  fEventNumber;
    Int_t  fRunNumber;

    Float_t fVertexX;
    Float_t fVertexY;
    Float_t fVertexZ;

    Float_t fVertex_frac; // fraction of track combinations used for vertex calculation

    Int_t fNVertices;
    UShort_t      fNumTracks;
    Int_t fNTriggerIds;
    Int_t fTriggerId;
    Int_t fTriggerId2;
    Int_t fTriggerId3;
    Int_t fTriggerId4;
    Int_t fTriggerId5;

    Short_t frefMult;
    Float_t fBBC_delta_t;
    Float_t fVPD_delta_t;
    Float_t fEP_Q_vectors_rc[6];
    Float_t fQtracks_used_Arr[6];


    TClonesArray* fTracks;      //->

public:
    StParticleEvent() :
        fRunId(-1), fEventId(-1), fEventNumber(-1), fRunNumber(-1), fVertexX(-1), fVertexY(-1), fVertexZ(-1),fVertex_frac(0),fNVertices(0),fNumTracks(0),fNTriggerIds(0)
        ,fTriggerId(0),fTriggerId2(0),fTriggerId3(0),fTriggerId4(0),fTriggerId5(0),frefMult(0)
        ,fBBC_delta_t(0),fVPD_delta_t(0)
    {
        fTracks      = new TClonesArray( "StParticleTrack", 10 );
	for ( Int_t i = 0; i < 6; ++i ) { 
	  fEP_Q_vectors_rc[i] = -999.; 
	  fQtracks_used_Arr[i] = -999.; 
	}
    }
        ~StParticleEvent()
        {
            delete fTracks;
            fTracks = NULL;
        }

        void     setRunId(Int_t r)           { fRunId = r;              }
        Int_t    getRunId() const            { return fRunId;           }

        void     setEventId(Int_t r)         { fEventId = r;            }
        Int_t    getEventId() const          { return fEventId;         }

        void     setEventNumber(Int_t r)     { fEventNumber = r;        }
        Int_t    getEventNumber() const      { return fEventNumber;     }

        void     setRunNumber(Int_t r)       { fRunNumber = r;          }
        Int_t    getRunNumber() const        { return fRunNumber;       }

        void     setVertexX(Float_t x)       { fVertexX = x;            }
        Float_t  getVertexX() const          { return fVertexX;         }

        void     setVertexY(Float_t y)       { fVertexY = y;            }
        Float_t  getVertexY() const          { return fVertexY;         }

        void     setVertexZ(Float_t z)       { fVertexZ = z;            }
        Float_t  getVertexZ() const          { return fVertexZ;         }

        void     setVertex_frac(Float_t f)   { fVertex_frac = f;        }
        Float_t  getVertex_frac() const      { return fVertex_frac;     }

        void     setNVertices(Int_t n)       { fNVertices = n;          }
        Int_t    getNVertices() const        { return fNVertices;       }

        void     setNTriggerIds(Int_t r)     { fNTriggerIds = r;        }
        Int_t    getNTriggerIds() const      { return fNTriggerIds;     }

        void     setTriggerId(Int_t r)       { fTriggerId = r;          }
        Int_t    getTriggerId() const        { return fTriggerId;       }

        void     setTriggerId2(Int_t r)      { fTriggerId2 = r;         }
        Int_t    getTriggerId2() const       { return fTriggerId2;      }

        void     setTriggerId3(Int_t r)      { fTriggerId3 = r;         }
        Int_t    getTriggerId3() const       { return fTriggerId3;      }

        void     setTriggerId4(Int_t r)      { fTriggerId4 = r;         }
        Int_t    getTriggerId4() const       { return fTriggerId4;      }

        void     setTriggerId5(Int_t r)      { fTriggerId5 = r;         }
        Int_t    getTriggerId5() const       { return fTriggerId5;      }

        void     setrefMult(Short_t r)       { frefMult = r;            }
        Short_t  getrefMult() const          { return frefMult;         }

        void     setBBC_delta_t(Float_t r)   { fBBC_delta_t = r;        }
        Float_t  getBBC_delta_t() const      { return fBBC_delta_t;     }

        void     setVPD_delta_t(Float_t r)   { fVPD_delta_t = r;        }
        Float_t  getVPD_delta_t() const      { return fVPD_delta_t;     }

        void     set_EP_Q_vectors_rc(Float_t arr[])   { for ( Int_t i = 0; i < 6; ++i ) fEP_Q_vectors_rc[i] = arr[i]; }
        const Float_t*  get_EP_Q_vectors_rc() const      { return &(fEP_Q_vectors_rc[0]);     }

        void     set_Qtracks_used_Arr(Float_t arr[])   { for ( Int_t i = 0; i < 6; ++i ) fQtracks_used_Arr[i] = arr[i]; }
        const Float_t*  get_Qtracks_used_Arr() const      { return &(fQtracks_used_Arr[0]);     }

        StParticleTrack* createTrack()
        {
            if (fNumTracks == fTracks->GetSize())
                fTracks->Expand( fNumTracks + 10 );
            if (fNumTracks >= 10000)
            {
                Fatal( "StParticleEvent::createTrack()", "ERROR: Too many tracks (>1000)!" );
                exit( 2 );
            }

            new((*fTracks)[fNumTracks++]) StParticleTrack;
            return (StParticleTrack*)((*fTracks)[fNumTracks - 1]);
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
        StParticleTrack* getTrack(UShort_t i) const
        {
            return i < fNumTracks ? (StParticleTrack*)((*fTracks)[i]) : NULL;
        }

        ClassDef(StParticleEvent,2)  // A simple event compiled of tracks
};


#endif // __STPARTICLEEVENT_H__
