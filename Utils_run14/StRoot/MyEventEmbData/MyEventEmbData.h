#ifndef MyEventEmbData_h
#define MyEventEmbData_h

#include "TObject.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "Riostream.h"
#include <vector>

using namespace std;

class MyTrackEmbData : public TObject {

  private:
    Float_t        Mass2     ;
    Float_t        nSigmaEl  ;
    Float_t        Momentum  ;
    Float_t        Beta      ;
    Float_t        dca       ;
    Float_t        pt        ;
    Float_t        TPCdEdx   ;
    Float_t        eta       ;
    Float_t        y         ;
    Float_t        phi       ;
    Float_t        nHitsPoss  ;
    Float_t        nHitsFit   ;
    Float_t        nHitsdEdx  ;
    Float_t        yLocal  ;

  public:
    MyTrackEmbData() :
      Mass2(-999),nSigmaEl(-999),Momentum(-999),Beta(-999),dca(-999),
      pt(-999),TPCdEdx(-999),eta(-999),y(-999),phi(-999),nHitsPoss(-999),
      nHitsFit(-999),nHitsdEdx(-999),yLocal(-999)
  {
  } 
    ~MyTrackEmbData() { }

    void       SetM2(Float_t xx)            { Mass2     = xx;     }
    void       SetnSigmaEl(Float_t xx)      { nSigmaEl  = xx;     }
    void       SetMomentum(Float_t xx)      { Momentum  = xx;     }
    void       SetBeta(Float_t xx)          { Beta      = xx;     }
    void       Setdca(Float_t xx)           { dca       = xx;     }
    void       Setpt(Float_t xx)            { pt        = xx;     }
    void       SetTPCdEdx(Float_t xx)       { TPCdEdx   = xx;     }
    void       Seteta(Float_t xx)           { eta       = xx;     }
    void       Sety(Float_t xx)             { y         = xx;     }
    void       Setphi(Float_t xx)           { phi       = xx;     }
    void       SetnHitsPoss(Float_t xx)     { nHitsPoss = xx;     }
    void       SetnHitsFit(Float_t xx)      { nHitsFit  = xx;     }
    void       SetnHitsdEdx(Float_t xx)     { nHitsdEdx = xx;     }
    void       SetyLocal(Float_t xx)        { yLocal    = xx;     }

    Float_t    GetM2() const { return Mass2; }
    Float_t    GetnSigmaEl () const { return nSigmaEl ; }
    Float_t    GetMomentum () const { return Momentum ; }
    Float_t    GetBeta     () const { return Beta     ; }
    Float_t    Getdca      () const { return dca      ; }
    Float_t    Getpt       () const { return pt       ; }
    Float_t    GetTPCdEdx  () const { return TPCdEdx  ; }
    Float_t    Geteta      () const { return eta      ; }
    Float_t    Gety        () const { return y        ; }
    Float_t    Getphi      () const { return phi      ; }
    Float_t    GetnHitsPoss() const { return nHitsPoss; }
    Float_t    GetnHitsFit () const { return nHitsFit ; }
    Float_t    GetnHitsdEdx() const { return nHitsdEdx; }
    Float_t    GetyLocal() const { return yLocal; }
    Short_t    GetCharge() const { return (Momentum<0) ? -1 : 1; }

    ClassDef(MyTrackEmbData,2)  
};

class MyEventEmbData : public TObject {

  private:
    Float_t         refMult  ;
    Float_t         RunId    ;
    Float_t        NrElectrons ;          
    Float_t        NrPions ;          
    Float_t        NrKaons ;          
    Float_t        NrProtons ;          
    Float_t        Prim_X   ;
    Float_t        Prim_Y   ;
    Float_t        Prim_Z   ;
    TClonesArray  *ElectronTracks;            //->
    TClonesArray  *PionTracks;            //->
    TClonesArray  *KaonTracks;            //->
    TClonesArray  *ProtonTracks;            //->

  public:
    MyEventEmbData() : refMult(-999),RunId(-999),NrElectrons(0),NrPions(0),NrKaons(0),NrProtons(0),Prim_X(-999),Prim_Y(-999),Prim_Z(-999) { 
      ElectronTracks = new TClonesArray("MyTrackEmbData", 10); 
      PionTracks = new TClonesArray("MyTrackEmbData", 10); 
      KaonTracks = new TClonesArray("MyTrackEmbData", 10); 
      ProtonTracks = new TClonesArray("MyTrackEmbData", 10); 
    }
    ~MyEventEmbData() {
      delete ElectronTracks; ElectronTracks = NULL;
      delete PionTracks; PionTracks = NULL;
      delete KaonTracks; KaonTracks = NULL;
      delete ProtonTracks; ProtonTracks = NULL;
    }

    void   SetHeader(Float_t rm, Float_t ri) { refMult = rm; RunId = ri; }
    void       SetPrim_X(Float_t px) { Prim_X = px; }
    void       SetPrim_Y(Float_t py) { Prim_Y = py; }
    void       SetPrim_Z(Float_t pz) { Prim_Z = pz; }

    MyTrackEmbData*   AddElectronTrack() { 
      if (NrElectrons == ElectronTracks->GetSize()) ElectronTracks->Expand( NrElectrons + 10 );
      new((*ElectronTracks)[NrElectrons++]) MyTrackEmbData;
      return (MyTrackEmbData*)((*ElectronTracks)[NrElectrons - 1]);
    } 
    MyTrackEmbData*   AddPionTrack() { 
      if (NrPions == PionTracks->GetSize()) PionTracks->Expand( NrPions + 10 );
      new((*PionTracks)[NrPions++]) MyTrackEmbData;
      return (MyTrackEmbData*)((*PionTracks)[NrPions - 1]);
    } 
    MyTrackEmbData*   AddKaonTrack() { 
      if (NrKaons == KaonTracks->GetSize()) KaonTracks->Expand( NrKaons + 10 );
      new((*KaonTracks)[NrKaons++]) MyTrackEmbData;
      return (MyTrackEmbData*)((*KaonTracks)[NrKaons - 1]);
    } 
    MyTrackEmbData*   AddProtonTrack() { 
      if (NrProtons == ProtonTracks->GetSize()) ProtonTracks->Expand( NrProtons + 10 );
      new((*ProtonTracks)[NrProtons++]) MyTrackEmbData;
      return (MyTrackEmbData*)((*ProtonTracks)[NrProtons - 1]);
    } 

    void clearTrackLists() {
      NrElectrons = 0; ElectronTracks->Clear();
      NrPions = 0; PionTracks->Clear();
      NrKaons = 0; KaonTracks->Clear();
      NrProtons = 0; ProtonTracks->Clear();
    }

    Float_t  GetRefMult  () const { return refMult  ; }
    Float_t  GetRunId    () const { return RunId    ; }
    Int_t         GetNrElectrons() const { return NrElectrons; }
    Int_t         GetNrPions() const { return NrPions; }
    Int_t         GetNrKaons() const { return NrKaons; }
    Int_t         GetNrProtons() const { return NrProtons; }
    Float_t       GetPrim_X() const { return Prim_X; }
    Float_t       GetPrim_Y() const { return Prim_Y; }
    Float_t       GetPrim_Z() const { return Prim_Z; }
    TClonesArray *GetElectronTracks() const {return ElectronTracks;}
    TClonesArray *GetPionTracks() const {return PionTracks;}
    TClonesArray *GetKaonTracks() const {return KaonTracks;}
    TClonesArray *GetProtonTracks() const {return ProtonTracks;}

    ClassDef(MyEventEmbData,1) 
};


#endif
