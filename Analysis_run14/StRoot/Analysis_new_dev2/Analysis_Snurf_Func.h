
#ifndef ANALYSIS_SNURF_FUNC_H
#define ANALYSIS_SNURF_FUNC_H
#include "StPhysicalHelixD.hh"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVector.hh"
#include "TLorentzVector.h"
//#include "StBTofUtil/StV0TofCorrection.h"
#include <algorithm>

//using namespace std;
//using namespace fastjet;

//StV0TofCorrection* v0tofcorr = new StV0TofCorrection();

// Analysis function are listed below
// Some functios like fPID or fDCA_Helix_Estimate are called several times in the particle analysis functions
// All particle analysis functions have a similar structure and can be called in parallel



void fPID(StPicoAlexTrack* track)
{
    // proton = 14, anti-Proton = 15, pi+ = 8, pi- = 9, K+ = 11, K- = 12, e+ = 2, e- = 3  (GEANT PID definitions)
    // PIDs_tracks[Ana][i] Array (global definition) i: Protons == 0, Kaons == 1, Pions == 2, Electrons == 3
    // If System == -1 there is no META hit! MetaQual is also negative
    for(Int_t i = 0; i < N_Ana; i++)
    {
        for(Int_t j = 0; j < N_max_PIDs_per_track; j++)
        {
            PIDs_track[i][j] = 0;
        }
    }
    Float_t nHitsFit   = track->nHitsFit();
    Float_t nHitsPoss  = track->nHitsMax();
    Float_t nHitsdEdx  = track->nHitsDedx();
    Float_t yLocal     = track->btofYLocal();
    Float_t dca        = track->dca();
    //Float_t nHitsFit = static_cast<Float_t>(track->nHitsFit()); // sign stores the charge
    //Float_t Polarity   = (nHitsFit >= 0) ? 1.0 : -1.0; // charge of the particle
    Float_t Polarity = static_cast<Float_t>(track->charge());
    StThreeVectorF p(track->pMom());  // primary momentum
    Float_t Momentum = p.mag(); // momentum
    Float_t Beta     = static_cast<Float_t>(track->btofBeta()); // beta

    Float_t Mass2 = -100.0;
    // calculate mass2
    if(track->btofMatchFlag() > 0 && track->btof() != 0 && Beta != 0)
    {
        Mass2 = Momentum*Momentum*(1.0/(Beta*Beta) - 1.0);
    }

    Float_t nSigmaEl   = (track->nSigmaElectron())*nsigma_scaling_fac;
    Float_t nSigmaPion = (track->nSigmaPion())*nsigma_scaling_fac;
    Float_t nSigmaKaon = (track->nSigmaKaon())*nsigma_scaling_fac;
    Float_t nSigmaP    = (track->nSigmaProton())*nsigma_scaling_fac;

    // PID for analysis number [1]
    PIDs_track[1][0] = 14; // proton

    // PID for analysis number [0]
    if(
       nSigmaP > -4.0
       && nSigmaP < 4.0
      )
    {
        if(
           ((track->btofMatchFlag() > 0 && track->btof() != 0 && Beta != 0)
            && Mass2 > 0.7*0.7   // 0.7
            && Mass2 < 1.2*1.2)  // 1.2
           || Beta < 0.0
          )
        {
            if(Polarity > 0)
            {
                PIDs_track[0][0] = 14; // proton
            }
            if(Polarity < 0)
            {
                PIDs_track[0][0] = 15; // anti-proton
            }
        }
    }

    if(
       nSigmaPion > -4.0
       && nSigmaPion < 4.0
      )
    {
        if(
           ((track->btofMatchFlag() > 0 && track->btof() != 0 && Beta != 0)
            && Mass2 > -0.2*0.2   // -0.2
            && Mass2 < 0.3*0.3)   // 0.3
           || Beta < 0.0
          )
        {
            if(Polarity > 0)
            {
                PIDs_track[0][2] = 8; // pi+
            }
            if(Polarity < 0)
            {
                PIDs_track[0][2] = 9; // pi-
            }
        }
    }

    if(
       nSigmaKaon > -2.0
       && nSigmaKaon < 2.0
      )
    {
        if(
           ((track->btofMatchFlag() > 0 && track->btof() != 0 && Beta != 0)
            && Mass2 > 0.4*0.4
            && Mass2 < 0.6*0.6)
           || Beta <= 0.0
          )
        {

            if(Polarity > 0)
            {
                PIDs_track[0][1] = 11; // K+
            }
            if(Polarity < 0)
            {
                PIDs_track[0][1] = 12; // K-
            }
        }
    }


    if(
       nSigmaEl > -2.0
       && nSigmaEl < 2.0
      )
    {
        if(
           ((track->btofMatchFlag() > 0 && track->btof() != 0 && Beta != 0)
            && Mass2 > -0.005
            && Mass2 < 0.006)
           //|| Beta < 0.0
          )
        {
            if(Polarity > 0)
            {
                PIDs_track[0][3] = 2; // e+
            }
            if(Polarity < 0)
            {
                PIDs_track[0][3] = 3; // e-
            }
        }
    }

    // PID for analysis number [2]
    if(
       nSigmaP > -3.0
       && nSigmaP < 3.0
      )
    {
        if(Polarity > 0)
        {
            PIDs_track[2][0] = 14; // proton
        }
        if(Polarity < 0)
        {
            PIDs_track[2][0] = 15; // anti-proton
        }
    }

    if(
       nSigmaPion > -3.0
       && nSigmaPion < 3.0
      )
    {
        if(Polarity > 0)
        {
            PIDs_track[2][2] = 8; // pi+
        }
        if(Polarity < 0)
        {
            PIDs_track[2][2] = 9; // pi-
        }
    }

    //cout << "nSigmaPion = " << nSigmaPion << ", nSigmaP = " << nSigmaP << ", PIDs_track[2][2] = " << PIDs_track[2][2] << ", Polarity = " << Polarity << endl;

    if(
       nSigmaKaon > -3.0
       && nSigmaKaon < 3.0
      )
    {
        if(Polarity > 0)
        {
            PIDs_track[2][1] = 11; // K+
        }
        if(Polarity < 0)
        {
            PIDs_track[2][1] = 12; // K-
        }
    }


    if(
       nSigmaEl > -3.0
       && nSigmaEl < 3.0
      )
    {
        if(Polarity > 0)
        {
            PIDs_track[2][3] = 2; // e+
        }
        if(Polarity < 0)
        {
            PIDs_track[2][3] = 3; // e-
        }
    }


    // PID for analysis number [3]
    // For Patrick's di-lepton analysis
    if( nHitsFit > 14 && nHitsPoss > 0 && (nHitsFit/nHitsPoss) > 0.52 && nHitsdEdx > 14 ) // eta cut on lepton tree!
    {
        if( Momentum > 0.0 && Momentum < 10.0 && Beta > 0.0 ) // min pt cut on lepton tree!
        {
            if( yLocal > -1.8 && yLocal < 1.8 )
            {
                Float_t BetaInv = 1.0/Beta-1.0;
                if ( BetaInv > -0.045 && BetaInv < 0.06 &&
                    ( Momentum < 0.23
                     || ( Momentum>=0.23 && Momentum<0.29 && BetaInv<-0.5*Momentum+0.175 )
                     || ( Momentum>=0.29 && BetaInv<0.03 ) )
                   )
                {// TOF Cut for Hadron reduction

                    if ( ( Momentum<0.4 && nSigmaEl>-2.6 && nSigmaEl<4.0 )
                        || ( Momentum>=0.4 && nSigmaEl<3.0 && nSigmaEl>0.0 )
                        || ( Momentum>=0.4 && Momentum<1.0 && nSigmaEl<=0.0 && nSigmaEl>Momentum-3.0 )
                        || ( Momentum>=1.0 && nSigmaEl<=0.0 && nSigmaEl>-2.0 ) )
                    {
                        if(Polarity > 0)
                        {
                            PIDs_track[3][3] = 2; // e+
                        }
                        if(Polarity < 0)
                        {
                            PIDs_track[3][3] = 3; // e-
                        }
                    }
                }
            }
        }
    }


    // For Patrick's single track tree analysis
    if( fAnalysisNum == 124 )
    {
        // electrons
        if ( nSigmaEl > -2 && nSigmaEl < 2 && nHitsFit > 14 && nHitsPoss > 0 ) {//&& (nHitsFit/nHitsPoss) > 0.52 )
            if(Polarity > 0) {
                PIDs_track[0][3] = 2; // global e+
            }
            if(Polarity < 0 && dca < 2 && Momentum > 0 ) {
                hnSigElDist[0]->Fill(nSigmaEl);
                PIDs_track[0][3] = 3; // primary e-
            }
        }
        // pi K p
        if ( nHitsFit > 14 && nHitsPoss > 0 && dca < 2 && Momentum > 0 ) {
            Float_t nsig = 0.2;
            if( nSigmaPion > -nsig && nSigmaPion < nsig ) {
                if(Polarity > 0) PIDs_track[0][2] = 8; // pi+
                if(Polarity < 0) PIDs_track[0][2] = 9; // pi-
            }
            if( nSigmaKaon > -nsig && nSigmaKaon < nsig ) {
                if(Polarity > 0) PIDs_track[0][1] = 11; // K+
                if(Polarity < 0) PIDs_track[0][1] = 12; // K-
            }
            if( nSigmaP > -nsig && nSigmaP < nsig ) {
                if(Polarity > 0) PIDs_track[0][0] = 14; // proton
                if(Polarity < 0) PIDs_track[0][0] = 15; // anti-proton
            }
        }
        return;
    }
}


void fillAlexEvent(StPicoDst* PicoDst_in, StPicoAlexEvent* PicoAlexEvent_out)
{
    StPicoEvent* event_in = PicoDst_in->event();


    // setters for event information
    PicoAlexEvent_out->setrunId(event_in->runId());
    PicoAlexEvent_out->seteventId(event_in->eventId());
    PicoAlexEvent_out->setfillId(event_in->fillId());
    PicoAlexEvent_out->setbField(event_in->bField());
    PicoAlexEvent_out->setprimaryVertex(event_in->primaryVertex());
    PicoAlexEvent_out->settriggerWord(event_in->triggerWord());
    PicoAlexEvent_out->setrefMultPos(event_in->refMultPos());
    PicoAlexEvent_out->setrefMultNeg(event_in->refMultNeg());
    PicoAlexEvent_out->setrefMultFtpcEast(event_in->refMultFtpcEast());
    PicoAlexEvent_out->setrefMultFtpcWest(event_in->refMultFtpcWest());
    PicoAlexEvent_out->setnVpdHitsEast(event_in->nVpdHitsEast());
    PicoAlexEvent_out->setnVpdHitsWest(event_in->nVpdHitsWest());
    PicoAlexEvent_out->setnT0(event_in->nT0());
    PicoAlexEvent_out->setvzVpd(event_in->vzVpd());
    PicoAlexEvent_out->setZDCx(event_in->ZDCx());
    PicoAlexEvent_out->setBBCx(event_in->BBCx());
    //PicoAlexEvent_out->setVpd();
    PicoAlexEvent_out->setZdcSumAdcEast(event_in->ZdcSumAdcEast());
    PicoAlexEvent_out->setZdcSumAdcWest(event_in->ZdcSumAdcWest());
    //PicoAlexEvent_out->setZdcSmdEastHorizontal();
    //PicoAlexEvent_out->setZdcSmdEastVertical();
    //PicoAlexEvent_out->setZdcSmdWestHorizontal();
    //PicoAlexEvent_out->setZdcSmdWestVertical();
    PicoAlexEvent_out->setbackgroundRate(event_in->backgroundRate());
    PicoAlexEvent_out->setbbcBlueBackgroundRate(event_in->bbcBlueBackgroundRate());
    PicoAlexEvent_out->setbbcYellowBackgroundRate(event_in->bbcYellowBackgroundRate());
    PicoAlexEvent_out->setbbcEastRate(event_in->bbcEastRate());
    PicoAlexEvent_out->setbbcWestRate(event_in->bbcWestRate());
    PicoAlexEvent_out->setzdcEastRate(event_in->zdcEastRate());
    PicoAlexEvent_out->setzdcWestRate(event_in->zdcWestRate());
    PicoAlexEvent_out->setspaceCharge(event_in->spaceCharge());
    PicoAlexEvent_out->setbtofTrayMultiplicity(event_in->btofTrayMultiplicity());
    PicoAlexEvent_out->setnumberOfGlobalTracks(event_in->numberOfGlobalTracks());
    //cout << "N global tracks = " << event_in->numberOfGlobalTracks() << endl;
    PicoAlexEvent_out->setranking(event_in->ranking());
    PicoAlexEvent_out->setnBEMCMatch(event_in->nBEMCMatch());
    for(Int_t ibbc = 0; ibbc < 24; ibbc++)
    {
        PicoAlexEvent_out->setbbcAdcEast(ibbc,event_in->bbcAdcEast(ibbc));
        PicoAlexEvent_out->setbbcAdcWest(ibbc,event_in->bbcAdcWest(ibbc));
    }


    // other user's functions
    PicoAlexEvent_out->setyear(event_in->year());
    PicoAlexEvent_out->setday(event_in->day());
    PicoAlexEvent_out->setenergy(event_in->energy());
#if 1
    PicoAlexEvent_out->setisMinBias(event_in->isMinBias());
    PicoAlexEvent_out->setisMBSlow(event_in->isMBSlow());
    PicoAlexEvent_out->setisCentral(event_in->isCentral());
    PicoAlexEvent_out->setisHT(event_in->isHT());
    PicoAlexEvent_out->setisHT11(event_in->isHT11());
    PicoAlexEvent_out->setisHT15(event_in->isHT15());
#endif
    StPicoTrack*     track_in;
    StPicoAlexTrack* track_out;
    Int_t erun_nevents_in  = PicoDst_in->numberOfTracks(); // number of tracks in this event
    //cout << "erun_nevents_in = " << erun_nevents_in << endl;


    for(UShort_t i = 0; i < erun_nevents_in; ++i) // loop over all tracks of the actual event
    {
        track_in = PicoDst_in->track( i ); // take the track
        track_out = PicoAlexEvent_out->createTrack();

        Int_t HFTmap = track_in->nHitsMapHFT();

        StDcaGeometry dcaG;
        dcaG.set(track_in->params(),track_in->errMatrix());
        StPhysicalHelixD helix_in = dcaG.helix();

        Int_t btofindex = track_in->bTofPidTraitsIndex();
        Int_t emcindex  = track_in->emcPidTraitsIndex();
        StPicoBTofPidTraits* btofpidtraits;
        StPicoEmcPidTraits*  emcpidtraits;
        if(btofindex >= 0)
        {
            btofpidtraits = PicoDst_in->btofPidTraits(btofindex);
        }
        if(emcindex >= 0)
        {
            emcpidtraits = PicoDst_in->emcPidTraits(emcindex);
        }


        // setters for track information
        track_out->setid(track_in->id());
        track_out->setchi2(track_in->chi2());
        track_out->setchi2prob(0.0);
        track_out->setgMom(helix_in.momentum(MAGFIELDFACTOR*event_in->bField()));
        track_out->setpMom(track_in->pMom());
        track_out->setorigin(helix_in.origin().x(), helix_in.origin().y(), helix_in.origin().z());
        track_out->setflowFlag(0);
        track_out->setQi(0,0);
        track_out->setdca(helix_in.geometricSignedDistance(event_in->primaryVertex()));
        track_out->setnHitsFit(track_in->nHitsFit()*track_in->charge());
        track_out->setnHitsMax(track_in->nHitsMax());
        track_out->setnHitsDedx(track_in->nHitsDedx());
        track_out->setdEdx(track_in->dEdx());
        track_out->setnSigmaPion(track_in->nSigmaPion());
        track_out->setnSigmaKaon(track_in->nSigmaKaon());
        track_out->setnSigmaProton(track_in->nSigmaProton());
        track_out->setnSigmaElectron(track_in->nSigmaElectron());
        track_out->setnHitsMapHFT(HFTmap);
        if(btofindex >= 0)
        {
            track_out->setbtofCellId(btofpidtraits->btofCellId());
            track_out->setbtofMatchFlag(btofpidtraits->btofMatchFlag());
            track_out->setbtof(btofpidtraits->btof());
            track_out->setbtofBeta(btofpidtraits->btofBeta());
            track_out->setbtofYLocal(btofpidtraits->btofYLocal());
            track_out->setbtofZLocal(btofpidtraits->btofZLocal());
            track_out->setbtofHisPos(btofpidtraits->btofHitPos().x(), btofpidtraits->btofHitPos().y(), btofpidtraits->btofHitPos().z());
        }
        if(emcindex >= 0)
        {
            track_out->setbemcId(emcpidtraits->bemcId());
            track_out->setadc0(emcpidtraits->adc0());
            track_out->sete0(emcpidtraits->e0());
            track_out->sete(emcpidtraits->e());
            track_out->setzDist(emcpidtraits->zDist());
            track_out->setphiDist(emcpidtraits->phiDist());
            track_out->setnEta(emcpidtraits->nEta());
            track_out->setnPhi(emcpidtraits->nPhi());
            track_out->setbtowId(emcpidtraits->btowId());
            track_out->setbtowId2(emcpidtraits->btowId2());
            track_out->setbtowId3(emcpidtraits->btowId3());
            track_out->sete1(emcpidtraits->e1());
            track_out->sete2(emcpidtraits->e2());
            track_out->sete3(emcpidtraits->e3());
            track_out->setetaTowDist(emcpidtraits->etaTowDist());
            track_out->setphiTowDist(emcpidtraits->phiTowDist());
        }
    }
}



Double_t calcDeterminant(StThreeVectorF& v1,StThreeVectorF& v2,StThreeVectorF& v3)
{
  // calculating the Determinant of a 3 x 3 Matrix 
  // with the column vectors [v1, v2, v3]
  // using the RULE of SARRUS
  //
  // | v1(0)   v2(0)   v3(0) |      | v1(0) v2(0) v3(0)|v1(0) v2(0)  .
  // |                       |      |  \\     \\     X   |  /     /    . 
  // |                       |      |   \\     \\   / \\  | /     /     . 
  // |                       |      |    \\     \\ /   \\ |/     /      . 
  // |                       |      |     \\     X     \\/     /       . 
  // |                       |      |      \\   / \\    /\\    /        .  
  // |                       |      |       \\ /   \\  / |\\  /         . 
  // | v1(1)   v2(1)   v3(1) |   =  | v1(1) v2(1) v3(1)|v1(1) v2(1)  .
  // |                       |      |       / \\    /\\  | /\\          . 
  // |                       |      |      /   \\  /  \\ |/  \\         . 
  // |                       |      |     /     \\/    \\/    \\        . 
  // |                       |      |    /      /\\    /\\     \\       . 
  // |                       |      |   /      /  \\  / |\\     \\      .  
  // |                       |      |  /      /    \\/  | \\     \\     . 
  // | v1(2)   v2(2)   v3(2) |      | v1(2) v2(2) v3(2)| v1(2) v2(2) .  
  //                                 /      /     /  \\     \\     \\   .
  //                                                                
  //                                -      -     -    +     +     +  .

  return ( v1(0) * v2(1) * v3(2) 
	   + v2(0) * v3(1) * v1(2) 
	   + v3(0) * v1(1) * v2(2) 
	   - v3(0) * v2(1) * v1(2) 
	   - v1(0) * v3(1) * v2(2) 
	   - v2(0) * v1(1) * v3(2)); 
}


StThreeVectorF calculatePointOfClosestApproach(StThreeVectorF &base1, StThreeVectorF &dir1,
								    StThreeVectorF &base2, StThreeVectorF &dir2)
{
  //  calculating point of closest approach
  //        
  //        from the equations of the straight lines of g and h 
  //        g: x1 = base1 + l * dir1 
  //        h: x2 = base2 + m * dir2 
  //        
  //        you can construct the following planes:
  //        
  //        E1: e1 = base1  +  a * dir1  +  b * (dir1 x dir2)
  //        E2: e2 = base2  +  s * dir2  +  t * (dir1 x dir2)
  //        
  //        now the intersection point of E1 with g2 = {P1} 
  //        and the intersection point of E2 with g1 = {P2}
  //        
  //        form the base points of the perpendicular to both straight lines.
  //        
  //        The point of closest approach is the middle point between P1 and P2: 
  //        
  //        vertex = (p2 - p1)/2
  // 
  //        E1 ^ g2:
  //
  //           e1 = x2
  //    -->    base1  +  a * dir1  +  b * (dir1 x dir2) = base2 + m * dir2 
  //    -->    base1 - base2 = m * dir2  -  a * dir1  -  b * (dir1 x dir2)       
  //                                          (m)
  //    -->    [ dir2, -dir1, -(dir1 x dir2)] (a) = base1 - base2        
  //                                          (b)
  //           
  //           using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D12 = det [dir2, -dir1, -(dir1 x dir2)] 
  //               = det [dir2,  dir1,  (dir1 x dir2)]
  //           
  //           Dm  = det [base1 - base2, -dir1, -(dir1 x dir2)]
  //               = det [base1 - base2,  dir1,  (dir1 x dir2)]
  //  
  //            m  = Dm/D12
  //           
  //           P1: p1 = x2(m)
  //                  = base2 + Dm/D12 * dir2
  //
  //        E2 ^ g1:
  //
  //           e2 = x1
  //    -->    base2  +  s * dir2  +  t * (dir1 x dir2) = base1 + l * dir1 
  //    -->    base2 - base1 = l * dir1  -  s * dir2  -  t * (dir1 x dir2)       
  //                                          (l)
  //    -->    [ dir1, -dir2, -(dir1 x dir2)] (s) = base2 - base1        
  //                                          (t)
  //           
  //           again using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D21 =  det [dir1, -dir2, -(dir1 x dir2)] 
  //               =  det [dir1,  dir2,  (dir1 x dir2)]
  //               = -det [dir2,  dir1,  (dir1 x dir2)]
  //               = -D12
  //           
  //           Dl  =  det [base2 - base1, -dir2, -(dir1 x dir2)]
  //               =  det [base2 - base1,  dir1,  (dir1 x dir2)]
  //               = -det [base1 - base2,  dir1,  (dir1 x dir2)]
  //
  //            l  =   Dl/D21
  //               = - Dl/D12
  //           
  //           P2: p2 = x1(m)
  //                  = base1 - Dl/D12 * dir1
  //           
  //           
  //           vertex = p1 + 1/2 * (p2 - p1)
  //                  = 1/2 * (p2 + p1)
  //                  = 1/2 *( (base1 + base2) +  1/D12 * ( Dm * dir2 - Dl * dir1) )
  //                      

  StThreeVectorF cross = dir1.cross(dir2); // cross product: dir1 x dir2

  // straight lines are either skew or have a cross point
	      
  StThreeVectorF diff = base1;
  diff-=base2; // Difference of two base vectors base1 - base2
		
  Double_t D;
  D =  calcDeterminant(dir2, dir1 ,cross);

  if (!(fabs(D) > 0.))
    {
      ::Warning(":calculatePointOfClosestApproach","Dirs and cross-product are lin. dependent: returning default Vertex (-20000,-20000,-20000)");
      return StThreeVectorF(-20000.,-20000.,-20000.);
    }

  Double_t Dm =  calcDeterminant(diff , dir1, cross);
  Double_t Dl = -calcDeterminant(diff , dir2, cross);

  StThreeVectorF vertex;
  StThreeVectorF dm;
  StThreeVectorF dl;

  dm = dir2;
  dm *= Dm;

  dl = dir1;
  dl *= Dl;

  vertex = dm - dl;

  vertex *= ((1.)/D);

  vertex+=base1;
  vertex+=base2;
  vertex*=0.5;

  return StThreeVectorF(vertex);
}



StThreeVectorF calculateCrossPoint(StThreeVectorF &base1, StThreeVectorF &dir1,
							StThreeVectorF &base2, StThreeVectorF &dir2)
{ 
  // calculating cross point 
  // taking all three equations into account solving the overdetermined set of lin. equations
  // of 
  // base1 + l * dir2 =  base1 + m * dir2 
  //
  // set of lin. equations:
  //  
  //   base1(0) + l * dir1(0) = base2(0) + m * dir2(0) 
  //   base1(1) + l * dir1(1) = base2(1) + m * dir2(1)
  //   base1(2) + l * dir1(2) = base2(2) + m * dir2(2) this line is ignored
  //
  //   written in matrix form
  //
  //        l
  //   M * |   | = base2 - base1
  //       \\ m /
  //
  //   M is a 3x2 matrix
  //     
  // to solve multiply the equation by the transposed Matrix of M from the left: M 
  //     
  //  T      /  l \\                                                               .
  // M * M * |    | = M  * (base2 - base1)
  //         \\ -m /
  // MIND THE '-' of m
  //
  //     / dir1(0) dir2(0) \\                                                      .
  //     |                 |    T   / dir1(0) dir1(1) dir1(2) \\                   .
  // M = | dir1(1) dir2(1) |,  M  = |                         |
  //     |                 |        \\ dir2(0) dir2(1) dir2(2) /                   .
  //     \\ dir1(2) dir2(2) /                                    
  //
  //  T      / (dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2))   (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))  \\ .
  // M * M = |                                                                                                                |
  //         \\ (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))   (dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2))  /                        
  //
  //  T       / d1d1 d1d2 \\                           .
  // M  * M = |           |
  //          \\ d1d2 d2d2 /
  //
  // diff = base2 - base1
  //
  //  T           /  (dir1(0)*diff(0) + dir1(1)*diff(1) + dir1(2)*diff(2)) \\         .
  // M  * diff =  |                                                        |
  //              \\  (dir2(0)*diff(0) + dir2(1)*diff(1) + dir2(2)*diff(2)) /
  //
  //  T           /  d1diff  \\                                          .
  // M  * diff =  |          |
  //              \\  d2diff  /
  // 
  // now the new Matrix set is to be solved by CRAMER'S Rule:
  // 
  // / d1d1 d1d2 \\   /  l \\   /  d1diff \\                   .
  // |           | * |    | = |          |
  // \\ d1d2 d2d2 /   \\ -m /   \\  d2diff /
  //
  //     | d1d1 d1d2 |
  // D = |           | = d1d1*d2d2 - d1d2*d1d2;
  //     | d1d2 d2d2 |
  // 
  //     | d1diff d1d2 |
  // Dl= |              | = d1diff*d2d2 - d1d2*d2diff;
  //     | d2diff d2d2 |              
  //
  // l = Dl/D = l_cross
  // 
  // vertex = base1 + l_cross * dir1
  //

  Double_t d1d1 = dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2);
  Double_t d2d2 = dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2);
  Double_t d1d2 = dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2);
  
  Double_t D = d1d1*d2d2 - (d1d2*d1d2);
  
  if (!(fabs(D) > 0.))
    {
      ::Warning("calculateCrossPoint","Error while calculating cross point ... eqns are lin. dependent:returning default Vertex (-20000,-20000,-20000)");
      return StThreeVectorF(-20000.,-20000.,-20000.);
    }

  Double_t d1diff = dir1(0)*(base2(0)-base1(0))+dir1(1)*(base2(1)-base1(1))+dir1(2)*(base2(2)-base1(2));
  Double_t d2diff = dir2(0)*(base2(0)-base1(0))+dir2(1)*(base2(1)-base1(1))+dir2(2)*(base2(2)-base1(2));

  Double_t Dlambda = d1diff*d2d2-d1d2*d2diff;
  
  Double_t lambda = Dlambda/D;
  
  StThreeVectorF vertex;
  vertex += dir1;
  vertex *= lambda;
  vertex += base1;

  //cout << "Cross point calculated" << endl;
  return StThreeVectorF(vertex);

 // return StThreeVectorF(-20000.,-20000.,-20000.);
}


StThreeVectorF calcVertexAnalytical(StThreeVectorF &base1, StThreeVectorF &dir1,
							 StThreeVectorF &base2, StThreeVectorF &dir2)
{
  // Calculates the Vertex of two straight lines define by the vectors base and dir
  //
  //      g: x1 = base1 + l * dir1 
  //      h: x2 = base2 + m * dir2 , where l,m are real numbers 
  //
  // 1. are g and h
  //       parallel / identical, i.e. are dir1 and dir2 linear dependent?
  //       
  //                                        /-                               
  //                                        |
  //                                        |   = 0    linear dependent, no unique solution, returning dummy  
  //      => cross product : dir1 x dir2 = -|  
  //                                        |  != 0    linear independent
  //                                        |
  //                                        \\-         
  //
  // 2. are g and h 
  //       skew or do they have a crossing point, i.e are dir1, dir2 and (base1 - base2) linear dependent ?
  //
  //                                                    /-                               
  //                                                    |
  //                                                    |   = 0    linear dependent
  //                                                    |          g and h are intersecting
  //                                                    |          calculating vertex as point of intersection
  //                                                    |
  //    => determinant: det[ dir1, dir2, base1-base2]= -|
  //                                                    |  != 0    linear independent
  //                                                    |          g and h are skew
  //                                                    |          calulating vertex as point of closest approach
  //                                                    |
  //                                                    \\-         
  //  
  // 3.
  //    (a) calculating intersection point
  //    (b) calculating point of closest approach



  // 1. exists a unique solution ?

  if ((dir1.cross(dir2)).mag()> 0.) // dir1 and dir2 linear independent
    {
      // straight lines are either skew or have a cross point

      StThreeVectorF diff = base1;
      diff-=base2; // Difference of two base vectors base1 - base2
      
      // 2. skew or intersecting ?
	
      if (fabs(calcDeterminant(dir2, dir1 ,diff))>0.) 
	{
	  // 3. (b) skew 
	  return StThreeVectorF(calculatePointOfClosestApproach(base1, dir1, base2, dir2));
	}
      else
	{
	  // 3. (a) intersection 
	  return StThreeVectorF(calculateCrossPoint(base1 ,dir1, base2 ,dir2));
	}
    }
  else
    {
      // dir1 and dir2 linear dependent -> g1 and g2 identical or parallel
      return StThreeVectorF(-10000000.,-10000000.,-10000000.);
    }
  return StThreeVectorF(-10000000.,-10000000.,-10000000.);
}

Double_t calculateMinimumDistanceStraightToPoint(StThreeVectorF &base, StThreeVectorF &dir,
									 StThreeVectorF &point)
{
  // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.mag()>0))
    {
      return -1000000.;
    }
  
  StThreeVectorF diff = base-point;

  StThreeVectorF cross = dir.cross(diff);
  
  return cross.mag()/dir.mag();
}

Double_t calculateMinimumDistance(StThreeVectorF &base1, StThreeVectorF &dir1,
							  StThreeVectorF &base2, StThreeVectorF &dir2)
{
  // calculates the minimum distance of two tracks given as parametric straights x = base + n * dir

  StThreeVectorF cross = dir1.cross(dir2);

  StThreeVectorF ab = base1 - base2;

  if ( !( fabs(cross.mag())>0.)) // dir1 || dir2
    {
      return calculateMinimumDistanceStraightToPoint(base1, dir1, base2);
    }
 
  return fabs(ab.dot(cross)/cross.mag());
}




Int_t fCross_points_Circles(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,
                           Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c)
{
    // (x1,y1) -> center of circle 1, r1 -> radius of circle 1
    // (x2,y2) -> center of circle 2, r2 -> radius of circle 2
    // (x1_c,y1_c) -> crossing point 1
    // (x2_c,y2_c) -> crossing point 2
    // Solution see Wikipedia

    Double_t s = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));  // distance between the center of the circles

    x1_c = 0;
    y1_c = 0;
    x2_c = 0;
    y2_c = 0;

    if(x1 != x2 && y1 != y2 && s < (r1+r2))
    {
        Double_t m  = (x1-x2)/(y2-y1);
        Double_t n  = (-r2*r2 + r1*r1 + y2*y2 - y1*y1 + x2*x2 - x1*x1)/(2.0*(y2-y1));
        Double_t p  = (2.0*(-x1 + m*(n-y1)))/(1.0 + m*m);
        Double_t q  = (x1*x1 + (n-y1)*(n-y1) -r1*r1)/(1.0 + m*m);
        Double_t p2 = (p/2.0)*(p/2.0);

        if(p2 >= q)
        {
            x1_c = (-p/2.0) + sqrt(p2 - q);
            x2_c = (-p/2.0) - sqrt(p2 - q);
            y1_c = m*x1_c + n;
            y2_c = m*x2_c + n;
            return 1;
        }
        else return 0;
    }
    else return 0;

}


Int_t fDCA_Helix_Estimate(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB)
{

    // Calculates the 2D crossing point, calculates the corresponding 3D point and returns pathA and pathB

    Double_t x1 = helixA.xcenter();
    Double_t y1 = helixA.ycenter();
    Double_t x2 = helixB.xcenter();
    Double_t y2 = helixB.ycenter();
    Double_t c1 = helixA.curvature();
    Double_t c2 = helixB.curvature();
    Double_t r1 = 0.0;
    Double_t r2 = 0.0;
    if(c1 != 0 && c2 != 0)
    {
        r1 = 1.0/c1;
        r2 = 1.0/c2;
    }

    Double_t x1_c = 0.0;
    Double_t y1_c = 0.0;
    Double_t x2_c = 0.0;
    Double_t y2_c = 0.0;

    Int_t bCross_points = fCross_points_Circles(x1,y1,r1,x2,y2,r2,x1_c,y1_c,x2_c,y2_c);

    //cout << "bCross_points = " << bCross_points << ", xyr(1) = {" << x1 << ", " << y1 << ", " << r1
    //    << "}, xyr(2) = {"  << x2 << ", " << y2 << ", " << r2 << "}, p1 = {" << x1_c << ", " << y1_c << "}, p2 = {" << x2_c << ", " << y2_c << "}" << endl;

    if(bCross_points == 0) return 0;

    StThreeVectorF pointA,pointB,pointA1,pointB1,pointA2,pointB2;

    Double_t path_lengthA_c1,path_lengthA_c2,path_lengthB_c1,path_lengthB_c2;

    // first crossing point for helix A
    pair< double, double > path_lengthA = helixA.pathLength(sqrt(x1_c*x1_c+y1_c*y1_c));
    Double_t path_lengthA1 = path_lengthA.first;
    Double_t path_lengthA2 = path_lengthA.second;
    pointA1 = helixA.at(path_lengthA1);
    pointA2 = helixA.at(path_lengthA2);
    if( ((x1_c-pointA1.x())*(x1_c-pointA1.x()) + (y1_c-pointA1.y())*(y1_c-pointA1.y())) <
       ((x1_c-pointA2.x())*(x1_c-pointA2.x()) + (y1_c-pointA2.y())*(y1_c-pointA2.y())))
    {
        path_lengthA_c1 = path_lengthA1;
    }
    else
    {
        path_lengthA_c1 = path_lengthA2;
    }

    // second crossing point for helix A
    path_lengthA = helixA.pathLength(sqrt(x2_c*x2_c+y2_c*y2_c));
    path_lengthA1 = path_lengthA.first;
    path_lengthA2 = path_lengthA.second;
    pointA1 = helixA.at(path_lengthA1);
    pointA2 = helixA.at(path_lengthA2);
    if( ((x2_c-pointA1.x())*(x2_c-pointA1.x()) + (y2_c-pointA1.y())*(y2_c-pointA1.y())) <
       ((x2_c-pointA2.x())*(x2_c-pointA2.x()) + (y2_c-pointA2.y())*(y2_c-pointA2.y())))
    {
        path_lengthA_c2 = path_lengthA1;
    }
    else
    {
        path_lengthA_c2 = path_lengthA2;
    }

    // first crossing point for helix B
    pair< double, double > path_lengthB = helixB.pathLength(sqrt(x1_c*x1_c+y1_c*y1_c));
    Double_t path_lengthB1 = path_lengthB.first;
    Double_t path_lengthB2 = path_lengthB.second;
    pointB1 = helixB.at(path_lengthB1);
    pointB2 = helixB.at(path_lengthB2);
    if( ((x1_c-pointB1.x())*(x1_c-pointB1.x()) + (y1_c-pointB1.y())*(y1_c-pointB1.y())) <
       ((x1_c-pointB2.x())*(x1_c-pointB2.x()) + (y1_c-pointB2.y())*(y1_c-pointB2.y())))
    {
        path_lengthB_c1 = path_lengthB1;
    }
    else
    {
        path_lengthB_c1 = path_lengthB2;
    }

    // second crossing point for helix B
    path_lengthB = helixB.pathLength(sqrt(x2_c*x2_c+y2_c*y2_c));
    path_lengthB1 = path_lengthB.first;
    path_lengthB2 = path_lengthB.second;
    pointB1 = helixB.at(path_lengthB1);
    pointB2 = helixB.at(path_lengthB2);
    if( ((x2_c-pointB1.x())*(x2_c-pointB1.x()) + (y2_c-pointB1.y())*(y2_c-pointB1.y())) <
       ((x2_c-pointB2.x())*(x2_c-pointB2.x()) + (y2_c-pointB2.y())*(y2_c-pointB2.y())))
    {
        path_lengthB_c2 = path_lengthB1;
    }
    else
    {
        path_lengthB_c2 = path_lengthB2;
    }

    pointA1 = helixA.at(path_lengthA_c1);
    pointA2 = helixA.at(path_lengthA_c2);

    pointB1 = helixB.at(path_lengthB_c1);
    pointB2 = helixB.at(path_lengthB_c2);

    if((pointA1-pointB1).mag() < (pointA2-pointB2).mag())
    {
        pathA = path_lengthA_c1;
        pathB = path_lengthB_c1;
        dcaAB = (pointA1-pointB1).mag();
    }
    else
    {
        pathA = path_lengthA_c2;
        pathB = path_lengthB_c2;
        dcaAB = (pointA2-pointB2).mag();
    }

    //cout << "3D-point (A1) = (" << pointA1.x() << ", " << pointA1.y() << ", " << pointA1.z() << "), calc = (" << x1_c << ", " << y1_c << ")" << endl;
    //cout << "3D-point (A2) = (" << pointA2.x() << ", " << pointA2.y() << ", " << pointA2.z() << "), calc = (" << x2_c << ", " << y2_c << ")" << endl;
    //cout << "3D-point (B1) = (" << pointB1.x() << ", " << pointB1.y() << ", " << pointB1.z() << "), calc = (" << x1_c << ", " << y1_c << ")" << endl;
    //cout << "3D-point (B2) = (" << pointB2.x() << ", " << pointB2.y() << ", " << pointB2.z() << "), calc = (" << x2_c << ", " << y2_c << ")" << endl;

    /*
    StThreeVectorF pointA,pointB,pointA1,pointB1,pointA2,pointB2;
    pointA = helixA.at(0); // 3D-vector of helixA at path 0
    pointB = helixB.at(0); // 3D-vector of helixA at path 0
    Double_t Delta_phiA = TMath::ATan2(pointA.y()-y1,pointA.x()-x1); // polar angle for path length 0 -> starting point for helix A
    Double_t Delta_phiB = TMath::ATan2(pointB.y()-y2,pointB.x()-x2); // polar angle for path length 0 -> starting point for helix B

    //Double_t phi_valA = 0.0;
    //Double_t x_new = x1 + r1*TMath::Cos(Delta_phi+phi_valA);
    //Double_t y_new = y1 + r1*TMath::Sin(Delta_phi+phi_valA);
    //cout << "3D-point = (" << pointA.x() << ", " << pointA.y() << "), calc = (" << x_new << ", " << y_new << ")" << endl;

    Double_t path_lengthA1 = (TMath::ATan2((y1_c-y1),(x1_c-x1)) - Delta_phiA)*r1;  // crossing point 1 for helix A
    Double_t path_lengthA2 = (TMath::ATan2((y2_c-y1),(x2_c-x1)) - Delta_phiA)*r1;  // crossing point 2 for helix A

    Double_t path_lengthB1 = (TMath::ATan2((y1_c-y2),(x1_c-x2)) - Delta_phiB)*r2;  // crossing point 1 for helix B
    Double_t path_lengthB2 = (TMath::ATan2((y2_c-y2),(x2_c-x2)) - Delta_phiB)*r2;  // crossing point 2 for helix B

    pointA1 = helixA.at(path_lengthA1);
    pointA2 = helixA.at(path_lengthA2);

    pointB1 = helixB.at(path_lengthB1);
    pointB2 = helixB.at(path_lengthB2);

    Double_t path_lengthTest = (TMath::ATan2((pointA.y()-y1),(pointA.x()-x1)) - Delta_phiA)*r1;  // crossing point 1 for helix A
    cout << "path_lengthTest = " << path_lengthTest << endl;


    cout << "3D-point (A1) = (" << pointA1.x() << ", " << pointA1.y() << ", " << pointA1.z() << "), calc = (" << x1_c << ", " << y1_c << ")" << endl;
    cout << "3D-point (A2) = (" << pointA2.x() << ", " << pointA2.y() << ", " << pointA2.z() << "), calc = (" << x2_c << ", " << y2_c << ")" << endl;
    cout << "3D-point (B1) = (" << pointB1.x() << ", " << pointB1.y() << ", " << pointB1.z() << "), calc = (" << x1_c << ", " << y1_c << ")" << endl;
    cout << "3D-point (B2) = (" << pointB2.x() << ", " << pointB2.y() << ", " << pointB2.z() << "), calc = (" << x2_c << ", " << y2_c << ")" << endl;
    */

    return 1;

}










void fHelixAtoPointdca(StThreeVectorF space_vec, StPhysicalHelixD helixA, Float_t &pathA, Float_t &dcaAB)
{
    // V1.1
    Float_t pA[2] = {0.0,-0.5}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA;
    for(Int_t r = 0; r < 2; r++)
    {
        testA     = helixA.at(pA[r]); // 3D-vector of helixA at path pA[r]
        distarray[r] = (testA-space_vec).mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.0005 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
        //    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pA[1]-pA[0])*scale; // the next length interval
            pA[0]     = pA[1] + scale_length; // the new path
            testA     = helixA.at(pA[0]); // 3D-vector of helixA at path pA[0]
            distarray[0] = (testA-space_vec).mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pA[0]-pA[1])*scale;
            pA[1]     = pA[0] + scale_length;
            testA     =  helixA.at(pA[1]); // 3D-vector of helixA at path pA[0]
            distarray[1] = (testA-space_vec).mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathA = pA[1];
        dcaAB = distarray[1];
    }
    //cout << "pathA = " << pathA << ", dcaAB = " << dcaAB << endl;
}

// V1.0
void fHelixABdca(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB)
{
    //cout << "Standard fHelixABdca called..." << endl;
    Float_t pA[2] = {0.0,0.0};
    Float_t pB[2] = {0.0,-0.1}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA, testB;
    for(Int_t r = 0; r < 2; r++)
    {
        testB     = helixB.at(pB[r]);  // 3D-vector of helixB point at path pB[r]

        Float_t pathA_dca = -999.0;
        Float_t dcaAB_dca = -999.0;
        fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
        testA = helixA.at(pathA_dca);
        //testA     = helixA.at(helixA.pathLength(testB)); // 3D-vector of helixA point at dca to testB

        distarray[r] = (testA-testB).mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 10.0;
    while(fabs(scale_length) > 0.0005 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pB[0] = " << pB[0]
        //    << ", pB[1] = " << pB[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pB[1]-pB[0])*scale; // the next length interval
            pB[0]     = pB[1] + scale_length; // the new path
            testB     = helixB.at(pB[0]); // new vector testB

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[0] = pathA_dca;
            //pA[0]     = helixA.pathLength(testB); // pathA at dca to helixB

            testA     = helixA.at(pA[0]); // new vector testA
            distarray[0] = (testA-testB).mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pB[0]-pB[1])*scale;
            pB[1]     = pB[0] + scale_length;
            testB     = helixB.at(pB[1]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[1] = pathA_dca;
            //pA[1]     = helixA.pathLength(testB); // pathA at dca to helixB

            testA     = helixA.at(pA[1]); // pathA at dca to helixB
            distarray[1] = (testA-testB).mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathB = pB[0];
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathB = pB[1];
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}

// V1.1
void fHelixAtoPointdca_start_params(StThreeVectorF space_vec, StPhysicalHelixD helixA, Float_t &pathA, Float_t &dcaAB, Float_t path_in_A)
{
    Float_t pA[2] = {(Float_t)(path_in_A+0.1),(Float_t)(path_in_A-0.1)}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA;
    for(Int_t r = 0; r < 2; r++)
    {
        testA     = helixA.at(pA[r]); // 3D-vector of helixA at path pA[r]
        distarray[r] = (testA-space_vec).mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.0005 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
        //    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pA[1]-pA[0])*scale; // the next length interval
            pA[0]     = pA[1] + scale_length; // the new path
            testA     = helixA.at(pA[0]); // 3D-vector of helixA at path pA[0]
            distarray[0] = (testA-space_vec).mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pA[0]-pA[1])*scale;
            pA[1]     = pA[0] + scale_length;
            testA     =  helixA.at(pA[1]); // 3D-vector of helixA at path pA[0]
            distarray[1] = (testA-space_vec).mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}

// V1.1
void fHelixABdca_start_params(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB, Float_t path_in_A, Float_t path_in_B)
{
    Float_t pA[2] = {(Float_t)(path_in_A+0.1),(Float_t)(path_in_A-0.1)};
    Float_t pB[2] = {(Float_t)(path_in_B+0.1),(Float_t)(path_in_B-0.1)}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA, testB;
    for(Int_t r = 0; r < 2; r++)
    {
        testB     = helixB.at(pB[r]);  // 3D-vector of helixB point at path pB[r]

        Float_t pathA_dca = -999.0;
        Float_t dcaAB_dca = -999.0;
        fHelixAtoPointdca_start_params(testB,helixA,pathA_dca,dcaAB_dca,path_in_A); // new helix to point dca calculation
        testA = helixA.at(pathA_dca);
        //testA     = helixA.at(helixA.pathLength(testB)); // 3D-vector of helixA point at dca to testB

        distarray[r] = (testA-testB).mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 10.0;
    while(fabs(scale_length) > 0.0005 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pB[0] = " << pB[0]
        //    << ", pB[1] = " << pB[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pB[1]-pB[0])*scale; // the next length interval
            pB[0]     = pB[1] + scale_length; // the new path
            testB     = helixB.at(pB[0]); // new vector testB

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca_start_params(testB,helixA,pathA_dca,dcaAB_dca,path_in_A); // new helix to point dca calculation
            pA[0] = pathA_dca;
            //pA[0]     = helixA.pathLength(testB); // pathA at dca to helixB

            testA     = helixA.at(pA[0]); // new vector testA
            distarray[0] = (testA-testB).mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pB[0]-pB[1])*scale;
            pB[1]     = pB[0] + scale_length;
            testB     = helixB.at(pB[1]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca_start_params(testB,helixA,pathA_dca,dcaAB_dca,path_in_A); // new helix to point dca calculation
            pA[1] = pathA_dca;
            //pA[1]     = helixA.pathLength(testB); // pathA at dca to helixB

            testA     = helixA.at(pA[1]); // pathA at dca to helixB
            distarray[1] = (testA-testB).mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathB = pB[0];
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathB = pB[1];
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}



void fHelixAtoLinedca(StThreeVectorF dirB, StThreeVectorF spaceB, StPhysicalHelixD helixA, Float_t &pathA, Float_t &pathB, Float_t &dcaAB)
{
    Float_t pA[2] = {0.0,0.0};
    Float_t pB[2] = {0.0,-5.0}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA, testB;
    for(Int_t r = 0; r < 2; r++)
    {
        testB     = spaceB+pB[r]*dirB;  // 3D-vector of helixB point at path pB[r]
        testA     = helixA.at(helixA.pathLength(testB)); // 3D-vector of helixA point at dca to testB
        distarray[r] = (testA-testB).mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.05 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pB[0] = " << pB[0]
        //    << ", pB[1] = " << pB[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pB[1]-pB[0])*scale; // the next length interval
            pB[0]     = pB[1] + scale_length; // the new path
            testB     = spaceB+pB[0]*dirB;  // 3D-vector of helixB point at path pB[r]
            pA[0]     = helixA.pathLength(testB); // pathA at dca to helixB
            testA     = helixA.at(pA[0]); // new vector testA
            distarray[0] = (testA-testB).mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pB[0]-pB[1])*scale;
            pB[1]     = pB[0] + scale_length;
            testB     = spaceB+pB[1]*dirB;  // 3D-vector of helixB point at path pB[r]
            pA[1]     = helixA.pathLength(testB); // pathA at dca to helixB
            testA     = helixA.at(pA[1]); // pathA at dca to helixB
            distarray[1] = (testA-testB).mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathB = pB[0];
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathB = pB[1];
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}

Float_t CalcTruncMean(Float_t data_array[][2], Float_t n_rms, Float_t n_percent, Int_t entries)
{
    // V1.2 Weights added
    // V1.1 Bug fixed, wrong mean value was used in second loop
    // First calculate the mean and rms of the original distribution
    Float_t Mean         = 0.0;
    Float_t Weight       = 0.0;
    Float_t RMS          = 0.0;
    Float_t RMS_weight   = 0.0;
    Int_t   n_above_zero = 0;
    for(Int_t i = 0; i < entries; i++)
    {
        if(data_array[i][0] > -900.0)
        {
            n_above_zero++;
            Weight += data_array[i][1];
            Mean   += data_array[i][0]*data_array[i][1];
        }
    }
    if(n_above_zero > 0 && Weight > 0.0)
    {
        Mean /= Weight;  // first mean value of original distribution
        for(Int_t i = 0; i < entries; i++)
        {
            if(data_array[i][0] > -900.0)
            {
                RMS_weight += pow(data_array[i][1],2);
                RMS        += pow((data_array[i][0]-Mean)*data_array[i][1],2);
            }
        }
        if(RMS_weight > 0.0) RMS = sqrt(RMS/RMS_weight); // first rms value of original distribution
    }
    else  // no entry -> return -1
    {
        return -999.0;
    }

    // Loop until n_percent is fullfilled
    Int_t n_above_zero_itt = n_above_zero;
    while(
          (Float_t)n_above_zero_itt/(Float_t)n_above_zero > n_percent
         )
    {
        Float_t Mean_itt       = 0.0;
        Float_t Weight_itt     = 0.0;
        Float_t RMS_itt        = 0.0;
        Float_t RMS_itt_weight = 0.0;
        n_above_zero_itt       = 0;

        for(Int_t i = 0; i < entries; i++)
        {
            if(data_array[i][0] > -900.0 && fabs(data_array[i][0] - Mean) < n_rms*RMS)
            {
                n_above_zero_itt++;
                Weight_itt += data_array[i][1];
                Mean_itt   += data_array[i][0]*data_array[i][1];
            }
        }
        if(n_above_zero_itt > 0 && Weight_itt > 0.0)
        {
            Mean_itt /= Weight_itt;  // itterated mean value
            for(Int_t i = 0; i < entries; i++)
            {
                if(data_array[i][0] > -900.0 && fabs(data_array[i][0] - Mean) < n_rms*RMS) // use the same selection as before (Mean,RMS)
                {
                    RMS_itt_weight += pow(data_array[i][1],2);
                    RMS_itt        += pow((data_array[i][0]-Mean_itt)*data_array[i][1],2); // calculate the new RMS_itt with new Mean_itt
                }
            }
            if(RMS_itt_weight > 0) RMS_itt = sqrt(RMS_itt/RMS_itt_weight); // itterated rms value
        }
        else  // no entry -> return last Mean value
        {
            return Mean;
        }
        if(Mean == Mean_itt) // no difference to value before
        {
            return Mean;
        }
        Mean = Mean_itt;
        RMS  = RMS_itt;
    }
    return Mean;
}



void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y1_val);
    Zero_line -> SetY2(y2_val);
    Zero_line -> SetLineWidth(LineWidth);
    Zero_line -> SetLineStyle(LineStyle);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}



Double_t calc_event_plane_weight(Double_t phiA, Double_t p_t, Double_t eta, Int_t RunId, Double_t EventVertexZ, Int_t Polarity, Double_t &phi_w)
{
    Double_t total_weight = 1.0;
    Double_t MaxEntry = 0.0;

    Double_t sign = 0.0;
    if(eta > 0.0) sign = 1.0;
    if(eta < 0.0) sign = 1.0; // -1.0 for 1st order event plane

    Double_t weight_phi = 1.0;

    //******************* Phi correction ****************************
    // determine the file time
    TString name;
    Long64_t val;
    Double_t day = 0.0;
    char NoP[50];
    sprintf(NoP,"%u",RunId);
    name = NoP[2];
    sscanf(name.Data(),"%Li",&val);
    day = 100.0 * val;
    name = NoP[3];
    sscanf(name.Data(),"%Li",&val);
    day = day + 10.0 * val;
    name = NoP[4];
    sscanf(name.Data(),"%Li",&val);
    day = day + 1.0 * val;

    Int_t day_bin     = (Int_t)day;
    Int_t z_bin       = (Int_t)((EventVertexZ-phi_corr_z_start)/phi_corr_delta_z);
    Int_t eta_bin     = (Int_t)((eta-phi_corr_eta_start)/phi_corr_delta_eta);
    Int_t pt_bin      = (Int_t)((p_t-phi_corr_pt_start)/phi_corr_delta_pt);

    if(pt_bin >= nPhi_corr_pt_bins) pt_bin = nPhi_corr_pt_bins-1;

    Int_t pol_bin        = 0;
    Int_t Phi_corr_bin   = 0;
    Int_t weight_counter = 0;
    if(Polarity > 0) pol_bin = 0;
    if(Polarity < 0) pol_bin = 1;
    if(
       day_bin    >= 0
       && eta_bin >= 0
       && z_bin   >= 0
       && pt_bin  >= 0
       && day_bin < nPhi_corr_days
       && eta_bin < nPhi_corr_eta_bins
       && pt_bin  < nPhi_corr_pt_bins
       && z_bin   < nPhi_corr_z_bins
       && hPhi_days_in[day_bin] == 1
      )
    {
        MaxEntry           = hPhi_corr_in_max_entry[day_bin][z_bin][eta_bin][pol_bin][pt_bin];
        //MaxEntry = 99; ???
        //cout << "MaxEntry = " << MaxEntry << endl;
        if(MaxEntry < 100) // phi correction file has too few entries, take correction from different day
        {
            Double_t MaxEntry_day_before   = 0.0;
            Double_t MaxEntry_day_after    = 0.0;
            Double_t weight_phi_day_before = 0.0;
            Double_t weight_phi_day_after  = 0.0;
            if((day_bin - 1) >= 0)
            {
                if(hPhi_days_in[day_bin - 1] == 1)
                {
                    MaxEntry_day_before = hPhi_corr_in_max_entry[day_bin-1][z_bin][eta_bin][pol_bin][pt_bin];
                    if(MaxEntry_day_before > 100)
                    {
                        Phi_corr_bin = hPhi_corr_in[day_bin-1][z_bin][eta_bin][pol_bin][pt_bin] ->FindBin(phiA);
                        weight_phi_day_before   = hPhi_corr_in[day_bin-1][z_bin][eta_bin][pol_bin][pt_bin] ->GetBinContent(Phi_corr_bin);
                        weight_counter++;
                    }
                }

            }
            if((day_bin + 1) < nPhi_corr_days)
            {
                if(hPhi_days_in[day_bin + 1] == 1)
                {
                    MaxEntry_day_after = hPhi_corr_in_max_entry[day_bin+1][z_bin][eta_bin][pol_bin][pt_bin];
                    if(MaxEntry_day_after > 100)
                    {
                        Phi_corr_bin = hPhi_corr_in[day_bin+1][z_bin][eta_bin][pol_bin][pt_bin] ->FindBin(phiA);
                        weight_phi_day_after   = hPhi_corr_in[day_bin+1][z_bin][eta_bin][pol_bin][pt_bin] ->GetBinContent(Phi_corr_bin);
                        weight_counter++;
                    }
                }
            }

            if(weight_counter > 0)
            {
                weight_phi = (weight_phi_day_before + weight_phi_day_after)/((Double_t)weight_counter);
                //cout << "weight_phi_day_before = " << weight_phi_day_before << ", weight_phi_day_after = " << weight_phi_day_after << endl;
            }
            else weight_phi = 1.0;
            //cout << "weight_phi = " << weight_phi << endl;

        }
        else // phi correction file has enough entries, everything ok
        {
            Phi_corr_bin = hPhi_corr_in[day_bin][z_bin][eta_bin][pol_bin][pt_bin] ->FindBin(phiA);
            weight_phi   = hPhi_corr_in[day_bin][z_bin][eta_bin][pol_bin][pt_bin] ->GetBinContent(Phi_corr_bin);
            //cout << "weight_phi = " << weight_phi << endl;
        }
    }
    //************** End phi correction ****************************

    if(weight_phi > 0.0) weight_phi = 1.0/weight_phi;
    else weight_phi = 1.0;

    if(weight_phi > 5.0) weight_phi = 5.0;  // don't correct for too much for small entries

    phi_w = weight_phi;

    Double_t p_t_weight = 1.0;
    if(p_t < 2.0)  p_t_weight = p_t;
    if(p_t >= 2.0) p_t_weight = 2.0;

    total_weight = p_t_weight*weight_phi*sign;

    //cout << "p_t_weight = " << p_t_weight << ", weight_phi = " << weight_phi << ", sign = "
    //   << sign << ", total_weight = " << total_weight << endl;

    return total_weight;
}


//************************** Functions for time-of-flight correction **************************

Float_t SquareRoot(Float_t square) {
  return (square>0) ? TMath::Sqrt(square) : -TMath::Sqrt(-square);
}


Float_t calcTof(TLorentzVector trackAB, Float_t PathAB) {
  // calculates new Tof of trackAB for AB -> A + B with Pathlength PathAB
  Float_t InvMassAB  = trackAB.M();
  Float_t MomentumAB = trackAB.P();
  Float_t BetaAB = TMath::Sqrt(1./(1+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));
  Float_t TofAB = PathAB / (BetaAB*clight);
  return TofAB;
}

Float_t calcPath(StPhysicalHelixD helixA, StThreeVectorF vectorAB, StThreeVectorF vectorTof) {
  // calculates PathLength of charged trackA (helixA) from vectorAB to vectorTof
  Float_t path_tof, dca_tof, path_dec, dca_dec;
  fHelixAtoPointdca(vectorTof,helixA,path_tof,dca_tof);
  fHelixAtoPointdca(vectorAB,helixA,path_dec,dca_dec);
  Float_t PathCorr = path_tof-path_dec;
  return PathCorr;
}

Float_t correctBeta4SingleDecay(StPicoAlexTrack trackA, TLorentzVector trackAB, StPhysicalHelixD helixA, StThreeVectorF vectorAB, Float_t PathAB) {
  // corrects beta in case of AB -> A + B
  Float_t TofAB = calcTof(trackAB, PathAB);
  Float_t TofCorr = trackA.btof() - TofAB;
  StThreeVectorF vectorTof = trackA.btofHisPos();
  Float_t PathCorr = calcPath(helixA,vectorAB,vectorTof);
  Float_t BetaCorr = PathCorr / (TofCorr*clight);

  //Double_t BetaA       = trackA.getBeta();  // Velocity after time-of-flight reconstruction
  //Double_t Mass2ACorr = trackA.getMomentum() * trackA.getMomentum() * (1.0/(BetaCorr*BetaCorr)-1.0);
  //if(Mass2ACorr > 100)
  //{
  //    cout << "Mass2ACorr = " << Mass2ACorr << ", PathCorr = " << PathCorr << ", t1 = " << trackA.btof() << ", TofAB = " << TofAB << ", TofCorr = " << TofCorr
  //        << ", BetaA = " << BetaA << endl;
  //}

  if ( debug_flag ) cout << "T) Tof = " << trackA.btof() << "\tTofCorr = " << TofCorr
			 //<< "\tPath = " << trackA.getPathLength() << "\tPathCorr = " << PathCorr
			 << "\tBeta = " << trackA.btofBeta() << "\tBetaCorr = " << BetaCorr << endl;
  return  BetaCorr;
}

Float_t correctBeta4DoubleDecay(StPicoAlexTrack trackA, TLorentzVector trackAB, TLorentzVector trackABC, StPhysicalHelixD helixA, StThreeVectorF vectorAB, Float_t PathAB, Float_t PathABC)
{
  // corrects beta in case of ABC -> AB + C -> A + B + C
  Float_t TofAB     = calcTof(trackAB,PathAB);    // Lambda
  Float_t TofABC    = calcTof(trackABC,PathABC);  // Omega
  Float_t TofCorr   = trackA.btof() - TofAB - TofABC;
  StThreeVectorF vectorTof = trackA.btofHisPos();
  Float_t PathCorr  = calcPath(helixA,vectorAB,vectorTof);
  Float_t BetaCorr  = PathCorr / (TofCorr*clight);
  if ( debug_flag ) cout << "T) Tof = " << trackA.btof() << "\tTofCorr = " << TofCorr
			 //<< "\tPath = " << trackA.getPathLength() << "\tPathCorr = " << PathCorr
			 << "\tBeta = " << trackA.btofBeta() << "\tBetaCorr = " << BetaCorr << endl;
  return  BetaCorr;
}
//******************************************************************************************


Double_t calc_phi_event_plane_2nd(Double_t mQx, Double_t mQy)
{
    // calculates the angle out of a 2D vector (mQx,mQy)
    // angle is in the range from 0..pi
    if(mQx == 0.0 && mQy == 0.0) return -400.0;
    Double_t Psi = 0.0;

    /*
    if(mQx != 0.0)
    {
        Psi = TMath::ATan(mQy/mQx);
        if(mQx > 0.0 && mQy > 0.0) Psi = Psi;
        if(mQx > 0.0 && mQy < 0.0) Psi = 2.0*TMath::Pi()+Psi;
        if(mQx < 0.0 && mQy > 0.0) Psi = 1.0*TMath::Pi()+Psi;
        if(mQx < 0.0 && mQy < 0.0) Psi = 1.0*TMath::Pi()+Psi;
        Psi = Psi/2.0;
    }
    else Psi = TMath::Pi()/4.0;
    */
    Psi = TMath::ATan2(mQy,mQx);
    Psi = Psi/2.0;

    return Psi;
}



//----------------------------------------------------------------------------------------------------------------------------------------------
Int_t calc_Decay_Properties(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t massA, Float_t massB, Float_t pA, Float_t pB,
                   StThreeVectorF vector_prim, Float_t cut_VerdistX, Float_t cut_dcaAB,
                            StThreeVectorF &vectorAB, TLorentzVector &ltrackA, TLorentzVector &ltrackB,
                            TLorentzVector &trackAB, Float_t &VerdistX, Float_t &dcaAB)
{
    Float_t dcaAB_f;
    Float_t pathA_test = -1.0;
    Float_t pathB_test = -1.0;
    Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);

    //cout << "pathA_test = " << pathA_test << ", pathB_test = " << pathB_test << endl;

    StThreeVectorF baseA          = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
    StThreeVectorF baseB          = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
    vectorAB  = baseA + baseB;
    vectorAB *= 0.5; // decay vertex

    StThreeVectorF dirA,dirB;
    dirA  = helixA.at(pathA_test-2.0) - helixA.at(pathA_test+2.0);
    dirB  = helixB.at(pathB_test-2.0) - helixB.at(pathB_test+2.0);

    StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

    dcaAB = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
    StThreeVectorF vectorABtoPrim = vectorAB_lin - vector_prim; // vector primary vertex to decay vertex
    VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay verter


    // calculate the scalar product with the approximated secondary vertex position
    dirA     = helixA.cat(pathA_test); // direction vector at dca for helixA
    dirB     = helixB.cat(pathB_test); // direction vector at dca for helixB
    dirA    *= pA/dirA.mag(); // new momentum vector at decay vertex
    dirB    *= pB/dirB.mag(); // new momentum vector at decay vertex

    ltrackA.SetXYZM(dirA.x(),dirA.y(),dirA.z(),massA);
    ltrackB.SetXYZM(dirB.x(),dirB.y(),dirB.z(),massB);
    trackAB = ltrackA + ltrackB; // mother particle
    StThreeVectorF dirY_lin;
    dirY_lin.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
    dirY_lin *= 1.0/dirY_lin.mag();
    Float_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim/vectorABtoPrim.mag());

    //cout << "i = " << i << ", j = " << j << ", VerdistX_lin = " << VerdistX_lin << ", dcaAB_lin = " << dcaAB_lin << ", scalarProduct_lin = " << scalarProduct_lin << endl;

    Float_t pathA_f, pathB_f;
    if( VerdistX > (cut_VerdistX*0.7) && dcaAB < (cut_dcaAB*1.3) && scalarProduct_lin > 0.0 )
    {

        if(fDCA_Helix_out == 1)
        {
            fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB,pathA_test,pathB_test); // calculate dca between two helices
        }
        else
        {
            fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB); // calculate dca between two helices
        }

        baseA    = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
        baseB    = helixB.at(pathB_f);  // space vector of helixB at dca to helixA

        vectorAB       = baseA + baseB;
        vectorAB       *= 0.5; // decay vertex

        vectorABtoPrim = vectorAB - vector_prim; // vector primary vertex to decay vertex
        VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

        if(VerdistX > cut_VerdistX && dcaAB < cut_dcaAB)
        {
            dirA     = helixA.cat(pathA_f); // direction vector at dca for helixA
            dirB     = helixB.cat(pathB_f); // direction vector at dca for helixB

            dirA *= pA/dirA.mag(); // new momentum vector at decay vertex
            dirB *= pB/dirB.mag(); // new momentum vector at decay vertex

            ltrackA.SetXYZM(dirA.x(),dirA.y(),dirA.z(),massA);
            ltrackB.SetXYZM(dirB.x(),dirB.y(),dirB.z(),massB);

            trackAB = ltrackA + ltrackB; // mother particle

            dirY_lin.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
            dirY_lin *= 1.0/dirY_lin.mag();
            scalarProduct_lin = dirY_lin.dot(vectorABtoPrim/vectorABtoPrim.mag());

            if(scalarProduct_lin > 0.0)
            {
                return 1;
            }
            else
            {
                return 0;
            }
            //cout << "massA = " << massA << ", massB = " << massB << ", pA = " << pA << " pB = " << pB << ", M = " << trackAB.M() << endl;
        }
        else return 0;
    }
    else return 0;
}
//----------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------
// This is the old event plane analysis function
Int_t EventPlane_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                          Int_t ParticleA, Int_t Ana_Num,Int_t run_events,Int_t event_number, Float_t &EP_Qx, Float_t &EP_Qy,
                          Float_t &EP_Qx_eta_pos, Float_t &EP_Qy_eta_pos, Float_t &EP_Qx_eta_neg, Float_t &EP_Qy_eta_neg)
{
    // Event vertex information
    event_A_ana            = picoDst_A;
    StThreeVectorF pVertex = event_A_ana   ->primaryVertex();
    Float_t EventVertexX   = pVertex.x();
    Float_t EventVertexY   = pVertex.y();
    Float_t EventVertexZ   = pVertex.z();
    Int_t   refMult        = event_A_ana->refMult();
    Int_t   RunId          = event_A_ana->runId();
    //Int_t TriggerWord      = event_A_ana->triggerWord();
    Float_t ZDCx           = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx           = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd          = event_A_ana->vzVpd();

    EP_Qx                = 0.0; // x value of event plane vector
    EP_Qy                = 0.0; // y value of event plane vector
    EP_Qx_eta_pos        = 0.0;   // sub event vectors for event plane resolution with pt and phi weight
    EP_Qy_eta_pos        = 0.0;
    EP_Qx_eta_neg        = 0.0;
    EP_Qy_eta_neg        = 0.0;
    EP_Qx_eta_pos_ptw    = 0.0;   // sub event vectors for event plane only with pt weight
    EP_Qy_eta_pos_ptw    = 0.0;
    EP_Qx_eta_neg_ptw    = 0.0;
    EP_Qy_eta_neg_ptw    = 0.0;
    EP_Qx1_eta_pos_ptw    = 0.0;   // sub event vectors for event plane only with pt weight - first order harmonic
    EP_Qy1_eta_pos_ptw    = 0.0;
    EP_Qx1_eta_neg_ptw    = 0.0;
    EP_Qy1_eta_neg_ptw    = 0.0;
    EP_Qx_ptw            = 0.0;   // event plane vector only with pt weight
    EP_Qy_ptw            = 0.0;
    Qtracks_used_eta_pos = 0;
    Qtracks_used_eta_neg = 0;
    Qtracks_used         = 0;
    EP_Qx_subA_ptw       = 0.0;   // sub event vectors for event plane resolution
    EP_Qy_subA_ptw       = 0.0;
    EP_Qx_subB_ptw       = 0.0;
    EP_Qy_subB_ptw       = 0.0;

    StThreeVectorF vector_prim, vector_prim_new;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    if(
       PID_counter_Array[Ana_Num][ParticleA] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && fabs(EventVertexZ) < vertex_z_cut
       //&& event_A_ana->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " particles *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA;
        StThreeVectorF vectorA, vectoratsA, vectorAB, vectorprimAB, vectornewA, vectornewB, vector_sum, dirA;
        StThreeVectorF testA, testAB, vectorABtoPrim, baseY, dirY;

        Double_t EP_Qx_phi_weight  = 0.0; // x value of event plane vector
        Double_t EP_Qy_phi_weight  = 0.0; // y value of event plane vector
        Double_t EP_Qx_no_weight   = 0.0; // x value of event plane vector
        Double_t EP_Qy_no_weight   = 0.0; // y value of event plane vector
        Double_t EP_Qx_subA        = 0.0;   // sub event vectors for event plane resolution
        Double_t EP_Qy_subA        = 0.0;
        Double_t EP_Qx_subB        = 0.0;
        Double_t EP_Qy_subB        = 0.0;
        Float_t phi_mean           = 0.0;
        Float_t mean_px            = 0.0;
        Float_t mean_py            = 0.0;
        Float_t mean_pz            = 0.0;
        Float_t mean_pxw           = 0.0;
        Float_t mean_pyw           = 0.0;
        Float_t mean_pzw           = 0.0;
        Float_t mean_py_pos        = 0.0;
        Float_t mean_py_neg        = 0.0;
        Float_t mean_px_pos        = 0.0;
        Float_t mean_px_neg        = 0.0;
        Int_t N_py_pos             = 0;
        Int_t N_py_neg             = 0;
        Int_t N_px_pos             = 0;
        Int_t N_px_neg             = 0;

        // For random sub event calculation
        TRandom ran_number;
        seed_number += run_events; // changes from event to event
        ran_number.SetSeed((UInt_t)seed_number);
        const Int_t N_max_EP_track = 3000;
        Double_t Qx_vector_array[N_max_EP_track];
        Double_t Qy_vector_array[N_max_EP_track];
        Double_t Qx_vector_array_ptw[N_max_EP_track];
        Double_t Qy_vector_array_ptw[N_max_EP_track];
        Int_t N_EP_tracks = 0;

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            //cout << "i = " << i << endl;
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];

            //StPicoAlexTrack trackA = *event->track( trackA_num );
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );

            StThreeVectorF vectorA = trackA.origin();

            // Requires: momentum, origin, signed Magnetic Field
            //           and Charge of particle (+/- 1)
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
            //helixA.setParameters(trackA.getcurvature(),trackA.getdipAngle(),trackA.getphase(),vectorA,trackA.geth());

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t BetaA       = trackA.btofBeta();

            if(
               nHitsFitA     > nHitsFitA_EP_cut  // 14
               //&& nHitsPossA > nHitsPossA_EP_cut
               && MomentumA  > MomentumA_EP_low_cut
               && MomentumA  < MomentumA_EP_high_cut
              )
            {
                if(
                   ((Float_t)nHitsFitA)/((Float_t)nHitsPossA)    > nHits_ratio_EP_cut
                  )
                {

                    Float_t pathA = -999.0;
                    Float_t dcaAB = -999.0;

                    //pathLength (const StThreeVector< double > &p, bool scanPeriods=true) const
                    //Double_t pathA_st = helixA.pathLength(vector_prim);
                    //Double_t dcaAB_st = helixA.distance(vector_prim);


                    fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);

                    //cout << "pathA_st = " << pathA_st << ", pathA = " << pathA
                    //<< ", dcaAB_st = " << dcaAB_st << ", dcaAB = " << dcaAB << endl;

                    //pathA = pathA_st;
                    //dcaAB = dcaAB_st;

                    vectoratsA  = helixA.cat(pathA);
                    Float_t eta = vectoratsA.pseudoRapidity();

                    if(
                       dcaA          < dcaAB_EP_cut  // 2.0
                       && fabs(eta)  < eta_EP_cut
                      )
                    {
                        //cout << "i = " << i << ", dcaAB = " << dcaAB << endl;
                        dirA = helixA.cat(pathA); // direction vector
                        Double_t phiA = dirA.phi(); // phiA has values from -pi..pi

                        //cout << "Length dirA = " << dirA.mag() << endl;

                        Float_t p_x = MomentumA*dirA.x();
                        Float_t p_y = MomentumA*dirA.y();
                        Float_t p_z = MomentumA*dirA.z();
                        Float_t p_t = sqrt(p_x*p_x + p_y*p_y);

                        //cout << "track = " << i << ", pathA = " << pathA << ", dcaAB = " << dcaAB << ", p_t = " << p_t << endl;

                        mean_px = mean_px + p_x;
                        mean_py = mean_py + p_y;
                        mean_pz = mean_pz + p_z;

                        // Old, only for 2nd order event plane
                        Float_t iQx   = TMath::Cos(2.0*phiA);
                        Float_t iQy   = TMath::Sin(2.0*phiA);

                        // Old, only for 2nd order event plane
                        Float_t iQx1   = TMath::Cos(1.0*phiA);
                        Float_t iQy1   = TMath::Sin(1.0*phiA);

                        Double_t phi_w;
                        Double_t total_weight = calc_event_plane_weight(phiA,p_t,eta,RunId,EventVertexZ,(Int_t)PolarityA,phi_w);

                        Double_t p_t_weight = 1.0;
                        if(p_t < 2.0)  p_t_weight = p_t;
                        if(p_t >= 2.0) p_t_weight = 2.0;
                        Double_t phi_weight = 1.0;
                        if(p_t_weight > 0) phi_weight = total_weight/p_t_weight;

                        //total_weight = phi_weight;

                        //cout << "Qtracks_used = " << Qtracks_used << ", total_weight = " << total_weight << ", phi_weight = " << phi_weight << ", p_t_weight = "
                        //    << p_t_weight << ", p_x = " << p_x << ", mean_px = " << mean_px << endl;

                        /*
                        if(iQy != 0.0)
                        {
                            Float_t Psi = TMath::ATan(iQy/iQx);
                            if(Psi < 0) Psi = TMath::Pi()+Psi;
                            cout << "phiA = " << TMath::RadToDeg()*phiA << ", 2*phiA = " << TMath::RadToDeg()*2.0*phiA << ", Psi = " << TMath::RadToDeg()*Psi
                                << ", Psi/2 = " << TMath::RadToDeg()*Psi/2.0 << endl;
                        }
                        */

                        mean_pxw = mean_pxw + phi_weight*p_x;
                        mean_pyw = mean_pyw + phi_weight*p_y;
                        mean_pzw = mean_pzw + phi_weight*p_z;

                        if(p_x > 0.0)
                        {
                            N_px_pos++;
                            mean_px_pos += phi_weight*p_x;
                        }
                        if(p_x < 0.0)
                        {
                            N_px_neg++;
                            mean_px_neg += phi_weight*p_x;
                        }
                        if(p_y > 0.0)
                        {
                            N_py_pos++;
                            mean_py_pos += phi_weight*p_y;
                        }
                        if(p_y < 0.0)
                        {
                            N_py_neg++;
                            mean_py_neg += phi_weight*p_y;
                        }

                        //phi_mean = phi_mean + sign*weight_phi*phiA;
                        phi_mean = phi_mean + phi_weight*phiA;

                        EP_Qx_no_weight  += iQx;
                        EP_Qy_no_weight  += iQy;

                        EP_Qx_phi_weight += iQx*phi_weight;
                        EP_Qy_phi_weight += iQy*phi_weight;

                        Float_t iQx_add = total_weight*iQx;
                        Float_t iQy_add = total_weight*iQy;

                        EP_Qx    += iQx_add;
                        EP_Qy    += iQy_add;

                        EP_Qx_ptw += iQx*p_t_weight;
                        EP_Qy_ptw += iQy*p_t_weight;


                        // For eta sub-event method
                        if(
                           fabs(eta) > eta_gap
                          )
                        {
                            if(eta >= 0.0)
                            {
                                EP_Qx_eta_pos     += iQx_add;
                                EP_Qy_eta_pos     += iQy_add;
                                EP_Qx_eta_pos_ptw += iQx*p_t_weight;
                                EP_Qy_eta_pos_ptw += iQy*p_t_weight;
                                EP_Qx1_eta_pos_ptw += iQx1*p_t_weight;
                                EP_Qy1_eta_pos_ptw += iQy1*p_t_weight;
                                Qtracks_used_eta_pos++;
                            }
                            if(eta < 0.0)
                            {
                                EP_Qx_eta_neg     += iQx_add;
                                EP_Qy_eta_neg     += iQy_add;
                                EP_Qx_eta_neg_ptw += iQx*p_t_weight;
                                EP_Qy_eta_neg_ptw += iQy*p_t_weight;
                                EP_Qx1_eta_neg_ptw += iQx1*p_t_weight;
                                EP_Qy1_eta_neg_ptw += iQy1*p_t_weight;
                                Qtracks_used_eta_neg++;
                            }
                        }


                        if(N_EP_tracks < N_max_EP_track)
                        {
                            Qx_vector_array[N_EP_tracks]     = iQx_add;
                            Qy_vector_array[N_EP_tracks]     = iQy_add;
                            Qx_vector_array_ptw[N_EP_tracks] = iQx*p_t_weight;
                            Qy_vector_array_ptw[N_EP_tracks] = iQy*p_t_weight;
                            N_EP_tracks++;
                        }

                        //cout << "i = " << i << ", iQx_add = " << iQx_add << ", iQy_add = " << iQy_add <<
                        //    ", EP_Qx = " << EP_Qx << ", EP_Qy = " << EP_Qy << endl;

                        //hphi  ->Fill(phiA,weight_phi);
                        //hsign ->Fill(sign);

                        Qtracks_used++;

                        Phi_Corr_NTDataArray[0]     =(Float_t)phiA;
                        Phi_Corr_NTDataArray[1]     =(Float_t)p_t;
                        Phi_Corr_NTDataArray[2]     =(Float_t)refMult;
                        Phi_Corr_NTDataArray[3]     =(Float_t)eta;
                        Phi_Corr_NTDataArray[4]     =(Float_t)EventVertexX;
                        Phi_Corr_NTDataArray[5]     =(Float_t)EventVertexY;
                        Phi_Corr_NTDataArray[6]     =(Float_t)EventVertexZ;
                        Phi_Corr_NTDataArray[7]     =(Float_t)RunId;
                        Phi_Corr_NTDataArray[8]     =(Float_t)BetaA;
                        Phi_Corr_NTDataArray[9]     =(Float_t)PolarityA;
                        Phi_Corr_NTDataArray[10]    =(Float_t)phi_weight;

                        //Phi_Corr_NT->Fill(Phi_Corr_NTDataArray);

                    }
                }
            }
        }


        //cout << "N_EP_tracks = " << N_EP_tracks << ", seed_number = " << seed_number << endl;
        //Int_t Event_ran_number = ran_number.Integer(run_events);  // will generate a random integer from 0..(run_events-1)
        // random generation of sub event A and sub event B
        // Event numbers -> A, odd events -> B
        // A random generated number selects the track, the Qx, Qy data array is changed by
        // copying the last value into the position of the used one
        for(Int_t r = 0; r < N_EP_tracks; r++)
        {
            Int_t track_ran_number = ran_number.Integer(N_EP_tracks-r);  // will generate a random integer from 0..(val-1)
            Float_t iQx_add     = Qx_vector_array[track_ran_number];
            Float_t iQy_add     = Qy_vector_array[track_ran_number];
            Float_t iQx_add_ptw = Qx_vector_array_ptw[track_ran_number];
            Float_t iQy_add_ptw = Qy_vector_array_ptw[track_ran_number];
            Qx_vector_array[track_ran_number] = Qx_vector_array[N_EP_tracks-r-1];
            Qy_vector_array[track_ran_number] = Qy_vector_array[N_EP_tracks-r-1];
            Qx_vector_array_ptw[track_ran_number] = Qx_vector_array_ptw[N_EP_tracks-r-1];
            Qy_vector_array_ptw[track_ran_number] = Qy_vector_array_ptw[N_EP_tracks-r-1];
            if((r % 2) == 0) // even event -> used for sub event analysis
            {
                EP_Qx_subA      += iQx_add;
                EP_Qy_subA      += iQy_add;
                EP_Qx_subA_ptw  += iQx_add_ptw;
                EP_Qy_subA_ptw  += iQy_add_ptw;
            }
            else
            {
                EP_Qx_subB      += iQx_add;
                EP_Qy_subB      += iQy_add;
                EP_Qx_subB_ptw  += iQx_add_ptw;
                EP_Qy_subB_ptw  += iQy_add_ptw;
            }
            //cout << "r = " << r << ", EP_Qx_subA = " << EP_Qx_subA << ", EP_Qx_subB = " << EP_Qx_subB
            //    << ", track_ran_number " << track_ran_number << ", iQx_add = "  << iQx_add << endl;
            //for(Int_t h = 0; h < N_EP_tracks-r; h++)
            //{
            //    cout << "h = " << h << ", Qx = " << Qx_vector_array[h] << endl;
            //}
        }

        //cout << "EP_Qx_ptw = " << EP_Qx_ptw << ", EP_Qx_subA_ptw+EP_Qx_subB_ptw = " << EP_Qx_subA_ptw+EP_Qx_subB_ptw << ", EP_Qx_subA_ptw = " << EP_Qx_subA_ptw << ", EP_Qx_subB_ptw = " << EP_Qx_subB_ptw
        //    << ", Qtracks_used = " << Qtracks_used << ", N_EP_tracks = "  << N_EP_tracks << endl;


        if(Qtracks_used > 0.0)
        {
            phi_mean = phi_mean/((Float_t)Qtracks_used);
            mean_px  = mean_px/((Float_t)Qtracks_used);
            mean_py  = mean_py/((Float_t)Qtracks_used);
            mean_pz  = mean_pz/((Float_t)Qtracks_used);
            mean_pxw = mean_pxw/((Float_t)Qtracks_used);
            mean_pyw = mean_pyw/((Float_t)Qtracks_used);
            mean_pzw = mean_pzw/((Float_t)Qtracks_used);
        }
        if(N_px_pos > 0)
        {
            mean_px_pos /=((Float_t)N_px_pos);
        }
        if(N_px_neg > 0)
        {
            mean_px_neg /=((Float_t)N_px_neg);
        }
        if(N_py_pos > 0)
        {
            mean_py_pos /=((Float_t)N_py_pos);
        }
        if(N_py_neg > 0)
        {
            mean_py_neg /=((Float_t)N_py_neg);
        }

        Double_t Psi            = calc_phi_event_plane_2nd(EP_Qx,EP_Qy);
        Double_t Psi_eta_pos    = calc_phi_event_plane_2nd(EP_Qx_eta_pos,EP_Qy_eta_pos);
        Double_t Psi_eta_neg    = calc_phi_event_plane_2nd(EP_Qx_eta_neg,EP_Qy_eta_neg);
        Double_t Psi_phi_weight = calc_phi_event_plane_2nd(EP_Qx_phi_weight,EP_Qy_phi_weight);
        Double_t Psi_no_weight  = calc_phi_event_plane_2nd(EP_Qx_no_weight,EP_Qy_no_weight);

        // cout << "Fill event plane ntuple" << endl;

        EventPlane_NTDataArray[0]     =(Float_t)EP_Qx;
        EventPlane_NTDataArray[1]     =(Float_t)EP_Qy;
        EventPlane_NTDataArray[2]     =(Float_t)EventVertexX;
        EventPlane_NTDataArray[3]     =(Float_t)EventVertexY;
        EventPlane_NTDataArray[4]     =(Float_t)EventVertexZ;
        EventPlane_NTDataArray[5]     =(Float_t)refMult;
        EventPlane_NTDataArray[6]     =(Float_t)event_number;
        EventPlane_NTDataArray[7]     =(Float_t)Qtracks_used;
        EventPlane_NTDataArray[8]     =(Float_t)phi_mean;
        EventPlane_NTDataArray[9]     =(Float_t)EP_Qx_subA;
        EventPlane_NTDataArray[10]    =(Float_t)EP_Qy_subA;
        EventPlane_NTDataArray[11]    =(Float_t)EP_Qx_subB;
        EventPlane_NTDataArray[12]    =(Float_t)EP_Qy_subB;
        EventPlane_NTDataArray[13]    =(Float_t)mean_px;
        EventPlane_NTDataArray[14]    =(Float_t)mean_py;
        EventPlane_NTDataArray[15]    =(Float_t)mean_pz;
        EventPlane_NTDataArray[16]    =(Float_t)mean_pxw;
        EventPlane_NTDataArray[17]    =(Float_t)mean_pyw;
        EventPlane_NTDataArray[18]    =(Float_t)mean_pzw;
        EventPlane_NTDataArray[19]    =(Float_t)Psi;
        EventPlane_NTDataArray[20]    =(Float_t)Psi_phi_weight;
        EventPlane_NTDataArray[21]    =(Float_t)Psi_no_weight;
        EventPlane_NTDataArray[22]    =(Float_t)EP_Qx_eta_pos;
        EventPlane_NTDataArray[23]    =(Float_t)EP_Qy_eta_pos;
        EventPlane_NTDataArray[24]    =(Float_t)EP_Qx_eta_neg;
        EventPlane_NTDataArray[25]    =(Float_t)EP_Qy_eta_neg;
        EventPlane_NTDataArray[26]    =(Float_t)Psi_eta_pos;
        EventPlane_NTDataArray[27]    =(Float_t)Psi_eta_neg;
        EventPlane_NTDataArray[28]    =(Float_t)Qtracks_used_eta_pos;
        EventPlane_NTDataArray[29]    =(Float_t)Qtracks_used_eta_neg;
        EventPlane_NTDataArray[30]    =(Float_t)mean_py_pos;
        EventPlane_NTDataArray[31]    =(Float_t)mean_py_neg;
        EventPlane_NTDataArray[32]    =(Float_t)mean_px_pos;
        EventPlane_NTDataArray[33]    =(Float_t)mean_px_neg;
        EventPlane_NTDataArray[34]    =(Float_t)EP_Qx_eta_pos_ptw;
        EventPlane_NTDataArray[35]    =(Float_t)EP_Qy_eta_pos_ptw;
        EventPlane_NTDataArray[36]    =(Float_t)EP_Qx_eta_neg_ptw;
        EventPlane_NTDataArray[37]    =(Float_t)EP_Qy_eta_neg_ptw;
        EventPlane_NTDataArray[38]    =(Float_t)EP_Qx_ptw;
        EventPlane_NTDataArray[39]    =(Float_t)EP_Qy_ptw;
        EventPlane_NTDataArray[40]    =(Float_t)n_tofmatch_prim;
        EventPlane_NTDataArray[41]    =(Float_t)n_non_primaries;
        EventPlane_NTDataArray[42]    =(Float_t)RunId;
        EventPlane_NTDataArray[43]    =(Float_t)ZDCx;
        EventPlane_NTDataArray[44]    =(Float_t)BBCx;
        EventPlane_NTDataArray[45]    =(Float_t)vzVpd;
        EventPlane_NTDataArray[46]    =(Float_t)EP_Qx_subA_ptw;
        EventPlane_NTDataArray[47]    =(Float_t)EP_Qy_subA_ptw;
        EventPlane_NTDataArray[48]    =(Float_t)EP_Qx_subB_ptw;
        EventPlane_NTDataArray[49]    =(Float_t)EP_Qy_subB_ptw;


        if(fAnalysisNum == 11)
        {
            EventPlane_NT->Fill(EventPlane_NTDataArray);
        }

        //cout << "EP_Qx_eta_pos_ptw = " << EP_Qx_eta_pos_ptw << ", EP_Qy_eta_pos_ptw = " <<  EP_Qy_eta_pos_ptw
        //    << ", Qtracks_used_eta_pos = " << Qtracks_used_eta_pos << ", refMult = " << refMult << endl;

        //cout << "mean_px = " << mean_px << endl;
        //cout << "Event plane (orig) = " << Psi << endl;

        return 1;
    }
    else return 0;
}
//----------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------
void BBC_analysis(StPicoAlexEvent* picoDst_A)
{
    event_A_ana            = picoDst_A;
    StThreeVectorF pVertex = event_A_ana   ->primaryVertex();
    BBC_x   = pVertex.x();
    BBC_y   = pVertex.y();
    BBC_z   = pVertex.z();
    BBC_refMult    = (Int_t)event_A_ana->refMult();
    BBC_runId              = event_A_ana->runId();
    //Float_t ZDCx           = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    //Float_t BBCx           = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    //Float_t vzVpd          = event_A_ana->vzVpd();
    //Int_t   eventId        = event_A_ana->eventId();

    if(
       (BBC_x*BBC_x + BBC_y*BBC_y) < Event_radius_cut*Event_radius_cut
       && fabs(BBC_z) < vertex_z_cut
       && event_A_ana->isMinBias()
      )
    {
        StStarBbcUtilities* bbc = new StStarBbcUtilities();
        for(UInt_t ipmt = 0; ipmt < 48; ipmt++)
        {
            const UInt_t pmtId = bbc->PmtId(ipmt);
            const UInt_t ew    = bbc->EastWest(ipmt);
            const UShort_t adc = (ew==0) ? event_A_ana->bbcAdcEast(pmtId) : event_A_ana->bbcAdcWest(pmtId);
            //const UInt_t innerOuter = bbc->InnerOuter(pmtId);

            BBC_ADC[ipmt] = adc;

            //cout << "ipmt = " << ipmt << ", pmtId = " << pmtId << ", ew = " << ew << ", innerOuter = "<< innerOuter << ", adc = " << adc << endl;
        }

        //cout << "BBC_refMult = " << BBC_refMult << ", BBC_x = " << BBC_x << ", BBC_y = " << BBC_y << ", BBC_z = " << BBC_z << endl;
        BBC_tree->Fill();
    }
}
//----------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------
// Event plane function to calculate Q-vectors for different harmonics + full TPC EP + eta-sub (gap) EP with different eta gap values + BBC EP
Int_t EventPlane_analysis_V2(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                          Int_t ParticleA, Int_t Ana_Num,Int_t run_events,Int_t event_number)
{
    //cout << "New event" << endl;
    // Event vertex information
    event_A_ana            = picoDst_A;
    StThreeVectorF pVertex = event_A_ana   ->primaryVertex();
    Float_t EventVertexX   = pVertex.x();
    Float_t EventVertexY   = pVertex.y();
    Float_t EventVertexZ   = pVertex.z();
    Int_t   refMult        = event_A_ana->refMult();
    Int_t   RunId          = event_A_ana->runId();
    Float_t ZDCx           = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx           = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd          = event_A_ana->vzVpd();
    Int_t   eventId        = event_A_ana->eventId();

    Int_t Qtracks_used_full = 0;

    StThreeVectorF vector_prim, vector_prim_new;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    Int_t   N_FTPC = 0;
    Float_t N_BBC  = 0;

    //cout << "tracks = " << PID_counter_Array[Ana_Num][ParticleA] << ", radius = " << (EventVertexX*EventVertexX + EventVertexY*EventVertexY) << ", z = " << EventVertexZ << ", mb = " << event_A_ana->isMinBias() <<  endl;

    if(
       PID_counter_Array[Ana_Num][ParticleA] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && fabs(EventVertexZ) < vertex_z_cut
       && event_A_ana->isMinBias()
      )
    {

        //----------------------------------------------------------------------------------------------------
        // Calculate the BBC event plane for the different harmonics
        // Based on Yadav's and Hiroshi's codes
        Float_t fQxBbcEastRaw[n_harmonics];
        Float_t fQyBbcEastRaw[n_harmonics];
        Float_t fQxBbcWestRaw[n_harmonics];
        Float_t fQyBbcWestRaw[n_harmonics];


        //StStarBbcUtilities* bbc = StStarBbcUtilities::Instance();
        //cout << "BBC Util" << endl;
        StStarBbcUtilities* bbc = new StStarBbcUtilities();


        //------------------------------------------------------------------------------
        // BBC event plane
        //   - ADC QA
        //   - Event plane reconstruction for east and west (2), and inner/outer (2)
        //     separately
        //------------------------------------------------------------------------------
        Float_t qw2e_BBC = 0.0;
        Float_t qw2w_BBC = 0.0;
        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics
        {
            fQxBbcEastRaw[n_harm] = 0.0; //  mQxBbcEast[har] = 0.0 ;
            fQyBbcEastRaw[n_harm] = 0.0; //  mQyBbcEast[har] = 0.0 ;
            fQxBbcWestRaw[n_harm] = 0.0; //  mQxBbcWest[har] = 0.0 ;
            fQyBbcWestRaw[n_harm] = 0.0; //  mQyBbcWest[har] = 0.0 ;
        }

        // Require at least 2 PMT hits for each arm
        UInt_t pmtSum[] = {0,0};

        Int_t BBC_hits_EW[2] = {0,0};

        for(UInt_t ipmt = 0; ipmt < 48; ipmt++)
        {
            const UInt_t pmtId = bbc->PmtId(ipmt) ;
            const UInt_t ew    = bbc->EastWest(ipmt) ;
            const UShort_t adc = (ew==0) ? event_A_ana->bbcAdcEast(pmtId) : event_A_ana->bbcAdcWest(pmtId) ;

            const UInt_t innerOuter = bbc->InnerOuter(pmtId) ;

            // Inner only
            if( innerOuter == 1 ) continue;

            // Remove 0 adc
            if( adc == 0 ) continue;

            // Remove pile-up
            if( adc > 4000 ) continue;

            if ( adc < 25 && fBeamTimeNum != 3 ) continue; // for all BES energies except 62.4 GeV

            if ( adc < 15 && fBeamTimeNum == 3 ) continue; // for 62.4 GeV

            if(fAnalysisNum == 111)
            {
                if(ipmt < 24)
                {
                    if(pmtId == 12)
                    {
                        N_BBC = adc;
                        //cout << "N_BBC = " << N_BBC << endl;
                        //cout << "phi = " << bbc->GetPhi(ew,pmtId)<< endl;
                    }
                    hBBC_ADC_tiles_Eeast_West[ew][ipmt] ->Fill(adc);
                    pBBC_ADC_tiles_vs_RefMult_Eeast_West[ew][ipmt] ->Fill(refMult,adc);
                    hBBC_ADC_tiles_vs_RefMult_Eeast_West[ew][ipmt] ->Fill(refMult,adc);
                }
                else
                {
                    hBBC_ADC_tiles_Eeast_West[ew][ipmt-24] ->Fill(adc);
                    pBBC_ADC_tiles_vs_RefMult_Eeast_West[ew][ipmt-24] ->Fill(refMult,adc);
                    hBBC_ADC_tiles_vs_RefMult_Eeast_West[ew][ipmt-24] ->Fill(refMult,adc);
                }
            }

            BBC_hits_EW[ew]++;


            // Remove background and pile-up
            //  - Lower cut off has been determined by looking at the ADC distributions (11.5 GeV)
            //  ADC cut off should be made during the analysis
            //  by using the BBC ADC sum
            //    if( adc <= 40 ) continue ; // 11.5 GeV
            //    if( adc <= 10 ) continue ; // 200 GeV
            //    if( adc <= 40 || adc > 4000 ) continue ; // 11.5 GeV
            //    if( adc <= 10 || adc > 4000 ) continue ; // 200 GeV (not final)
            //    if( adc <= bbc->GetAdcCut(mEnergy) ) continue ;

            pmtSum[ew]++;

            const Float_t phi       = bbc->GetPhi(ew, pmtId);
            const Double_t weight   = static_cast<Double_t>(adc);

            //cout << "phi = " << phi << ", weight = " << weight << endl;

            // Calculate q-vectors
            for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics
            {
                const Float_t harmonics = n_harm+1.0;
                if( ew == 0 ) // East
                {
                    fQxBbcEastRaw[n_harm] += weight * TMath::Cos(harmonics*phi);
                    fQyBbcEastRaw[n_harm] += weight * TMath::Sin(harmonics*phi);
                    if( n_harm == 0 ) qw2e_BBC += weight * weight;
                }
                else // West
                {
                    fQxBbcWestRaw[n_harm] += weight * TMath::Cos(harmonics*phi);
                    fQyBbcWestRaw[n_harm] += weight * TMath::Sin(harmonics*phi);
                    if( n_harm == 0 ) qw2w_BBC += weight * weight;
                }
            }
        }

        hrefMult_BBC_hits_East_West[0] ->Fill(refMult,BBC_hits_EW[0]);
        hrefMult_BBC_hits_East_West[1] ->Fill(refMult,BBC_hits_EW[1]);
        //cout << "QxE = " << fQxBbcEastRaw[1] << ", QyE = " << fQyBbcEastRaw[1]
        //    << ", QxW = " << fQxBbcWestRaw[1] << ", QyW = " << fQyBbcWestRaw[1] << endl;

#if 0
        // Qw for BBC
        Float_t mQwBbcEast = 0.0;
        Float_t mQwBbcWest = 0.0;

        Float_t mQxBbcEastRaw[n_harmonics];
        Float_t mQyBbcEastRaw[n_harmonics];
        Float_t mQxBbcWestRaw[n_harmonics];
        Float_t mQyBbcWestRaw[n_harmonics];

        if( qw2e_BBC > 0.0 ) mQwBbcEast = TMath::Sqrt(qw2e_BBC);
        if( qw2w_BBC > 0.0 ) mQwBbcWest = TMath::Sqrt(qw2w_BBC);

        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics
        {
            // Check # of pmt hits
            for(Int_t ew=0; ew<2; ew++)
            {
                if( pmtSum[ew] == 0 )
                {
                    if(ew==0)
                    {
                        mQxBbcEastRaw[n_harm] = -9999.0;
                        mQyBbcEastRaw[n_harm] = -9999.0;
                    }
                    else
                    {
                        mQxBbcWestRaw[n_harm] = -9999.0;
                        mQyBbcWestRaw[n_harm] = -9999.0;
                    }
                }
            }
        }
#endif

        // Reduced q-vectors
        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics
        {
            if( qw2e_BBC != 0.0 )
            {
                fQxBbcEastRaw[n_harm] /= TMath::Sqrt(qw2e_BBC);
                fQyBbcEastRaw[n_harm] /= TMath::Sqrt(qw2e_BBC);
            }
            else{
                fQxBbcEastRaw[n_harm] = -9999.0;
                fQyBbcEastRaw[n_harm] = -9999.0;
            }

            if( qw2w_BBC != 0.0 )
            {
                fQxBbcWestRaw[n_harm] /= TMath::Sqrt(qw2w_BBC);
                fQyBbcWestRaw[n_harm] /= TMath::Sqrt(qw2w_BBC);
            }
            else
            {
                fQxBbcWestRaw[n_harm] = -9999.0;
                fQyBbcWestRaw[n_harm] = -9999.0;
            }

            // Check # of pmt hits
            for(Int_t ew=0; ew<2; ew++)
            {
                if( pmtSum[ew] == 0 )
                {
                    if(ew==0)
                    {
                        fQxBbcEastRaw[n_harm] = -9999.0;
                        fQyBbcEastRaw[n_harm] = -9999.0;
                    }
                    else
                    {
                        fQxBbcWestRaw[n_harm] = -9999.0;
                        fQyBbcWestRaw[n_harm] = -9999.0;
                    }
                }
            }
        }
        //----------------------------------------------------------------------------------------------------



        //----------------------------------------------------------------------------------------------------
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " particles *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA;
        StThreeVectorF vectorA, vectoratsA, vectorAB, vectorprimAB, vectornewA, vectornewB, vector_sum, dirA;
        StThreeVectorF testA, testAB, vectorABtoPrim, baseY, dirY;

        // For random sub event calculation
        TRandom ran_number;
        seed_number += run_events; // changes from event to event
        ran_number.SetSeed((UInt_t)seed_number);
        const Int_t N_max_EP_track = 3000;
        Int_t N_EP_tracks = 0;
        //----------------------------------------------------------------------------------------------------



        //----------------------------------------------------------------------------------------------------
        // Define variables for harmonics and eta gap values
        Double_t iQxy[2][n_harmonics]; // [x,y][harmonic]

        // For eta (gap) sub event plane method
        Double_t EP_Qxy_eta_pos_neg_array[n_harmonics][n_eta_gap_values][2][2][3]; // [n harmonic][eta gap value][positive eta, negative eta][Qx, Qy][all charges,+,-]
        Double_t EP_Qxy_eta_pos_neg_ptw_array[n_harmonics][n_eta_gap_values][2][2][3]; // [n harmonic][eta gap value][positive eta, negative eta][Qx, Qy][all charges,+,-]

        // For full TPC event plane method
        Double_t EP_Qxy_array[n_harmonics][2]; // [n harmonic][Qx, Qy]
        Double_t EP_Qxy_ptw_array[n_harmonics][2]; // [n harmonic][Qx, Qy]

        // For full TPC event plane resolution -> randum sub events
        Double_t EP_Qxy_array_sub[n_harmonics][2][N_max_EP_track]; // [n harmonic][Qx, Qy][track number]
        Double_t EP_Qxy_ptw_array_sub[n_harmonics][2][N_max_EP_track]; // [n harmonic][Qx, Qy][track number]

        Int_t Qtracks_used_eta_pos_neg_array[n_harmonics][n_eta_gap_values][2][3]; // [n harmonic][eta gap value][positive eta, negative eta][all charges,+,-]
        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics
        {
            for(Int_t n_eta_gap = 0; n_eta_gap < n_eta_gap_values; n_eta_gap++)  // loop over different eta gap values
            {
                for(Int_t n_pos_neg = 0; n_pos_neg < 2; n_pos_neg++) // loop over positive and negative eta hemisphere
                {
                    for(Int_t n_qxy = 0; n_qxy < 2; n_qxy++) // loop over Qx and Qy
                    {
                        for(Int_t n_charge = 0; n_charge < 3; n_charge++) // loop over charges: all, +, -
                        {
                            EP_Qxy_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg][n_qxy][n_charge]      = 0.0;
                            EP_Qxy_eta_pos_neg_ptw_array[n_harm][n_eta_gap][n_pos_neg][n_qxy][n_charge]  = 0.0;
                        }
                    }
                    for(Int_t n_charge = 0; n_charge < 3; n_charge++) // loop over charges: all, +, -
                    {
                        Qtracks_used_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg][n_charge] = 0;
                    }
                }
            }
        }

        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++)
        {
            for(Int_t n_qxy = 0; n_qxy < 2; n_qxy++) // loop over Qx and Qy
            {
                EP_Qxy_array[n_harm][n_qxy]     = 0.0;
                EP_Qxy_ptw_array[n_harm][n_qxy] = 0.0;
            }
        }
        //----------------------------------------------------------------------------------------------------



        //----------------------------------------------------------------------------------------------------
        // Loop over all particles
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            //cout << "i = " << i << endl;
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];

            //StPicoAlexTrack trackA = *event->track( trackA_num );
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            StThreeVectorF vectorA = trackA.origin();

            // Requires: momentum, origin, signed Magnetic Field
            //           and Charge of particle (+/- 1)
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
            //helixA.setParameters(trackA.getcurvature(),trackA.getdipAngle(),trackA.getphase(),vectorA,trackA.geth());

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            //Float_t BetaA       = trackA.btofBeta();

            if(
                nHitsPossA > 0
               //nHitsFitA     > nHitsFitA_EP_cut  // 14
               //&& nHitsPossA > nHitsPossA_EP_cut
               //&& MomentumA  > MomentumA_EP_low_cut
               //&& MomentumA  < MomentumA_EP_high_cut
              )
            {
                if(
                   ((Float_t)nHitsFitA)/((Float_t)nHitsPossA)    > nHits_ratio_EP_cut
                  )
                {

                    Float_t pathA = -999.0;
                    Float_t dcaAB = -999.0;

                    fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);

                    vectoratsA  = helixA.cat(pathA);
                    Float_t eta = vectoratsA.pseudoRapidity();

                    heta_EP ->Fill(eta);
                    if(eta > 3.44 && eta < 3.95)
                    {
                        dirA = helixA.cat(pathA); // direction vector
                        Double_t phiA = dirA.phi(); // phiA has values from -pi..pi
                        //cout << "phiA = " << phiA << endl;
                        h_FTPC_phi ->Fill(phiA);
                        if(phiA < (-1.5708+0.26) && phiA > (-1.5708-0.26))
                        {
                            N_FTPC++;
                            //cout << "N_FTPC = " << N_FTPC << endl;
                        }
                    }

                    if(
                       dcaA          < dcaAB_EP_cut  // 2.0
                       && fabs(eta)  < eta_EP_cut    // 1.0
                      )
                    {
                        //cout << "i = " << i << ", dcaAB = " << dcaAB << endl;
                        dirA = helixA.cat(pathA); // direction vector
                        Double_t phiA = dirA.phi(); // phiA has values from -pi..pi

                        Float_t p_x = MomentumA*dirA.x();
                        Float_t p_y = MomentumA*dirA.y();
                        //Float_t p_z = MomentumA*dirA.z();
                        Float_t p_t = sqrt(p_x*p_x + p_y*p_y);

                        //------------------------------------------------------------------------
                        // Calculate the Q-vectors for the different harmonics
                        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++)
                        {
                            iQxy[0][n_harm] = TMath::Cos( ((Double_t)(n_harm+1))*phiA );
                            iQxy[1][n_harm] = TMath::Sin( ((Double_t)(n_harm+1))*phiA );
                        }
                        //------------------------------------------------------------------------



                        //------------------------------------------------------------------------
                        // Calculate event plane vector weights
                        Double_t phi_w;
                        Double_t total_weight = calc_event_plane_weight(phiA,p_t,eta,RunId,EventVertexZ,(Int_t)PolarityA,phi_w);

                        Double_t p_t_weight = 1.0;
                        if(p_t < 2.0)  p_t_weight = p_t;
                        if(p_t >= 2.0) p_t_weight = 2.0;
                        Double_t phi_weight = 1.0;
                        if(p_t_weight > 0) phi_weight = total_weight/p_t_weight;
                        //------------------------------------------------------------------------



                        //------------------------------------------------------------------------
                        // Calculate the Q-vectors for the different harmonics for full TPC event plane method
                        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++)
                        {
                            for(Int_t n_qxy = 0; n_qxy < 2; n_qxy++) // loop over Qx and Qy
                            {
                                EP_Qxy_array[n_harm][n_qxy]     += iQxy[n_qxy][n_harm]*total_weight;
                                EP_Qxy_ptw_array[n_harm][n_qxy] += iQxy[n_qxy][n_harm]*p_t_weight;
                            }
                        }
                        //------------------------------------------------------------------------



                        //------------------------------------------------------------------------
                        // For eta sub-event method, new, with different eta-gap values
                        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics
                        {
                            for(Int_t n_eta_gap = 0; n_eta_gap < n_eta_gap_values; n_eta_gap++)  // loop over different eta gap values
                            {
                                if(
                                   fabs(eta) > eta_gap_array[n_eta_gap]
                                  )
                                {
                                    Int_t n_pos_neg = 0; // 0 == positive eta, 1 == negative eta
                                    if(eta >= 0.0) {n_pos_neg = 0;}
                                    else {n_pos_neg = 1;}

                                    for(Int_t n_qxy = 0; n_qxy < 2; n_qxy++) // loop over Qx and Qy
                                    {
                                        EP_Qxy_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg][n_qxy][0]          += iQxy[n_qxy][n_harm]*total_weight;
                                        EP_Qxy_eta_pos_neg_ptw_array[n_harm][n_eta_gap][n_pos_neg][n_qxy][0]      += iQxy[n_qxy][n_harm]*p_t_weight;
                                        if(PolarityA > 0.0) // use only positive charged particles for the event plane
                                        {
                                            EP_Qxy_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg][n_qxy][1]      += iQxy[n_qxy][n_harm]*total_weight;
                                            EP_Qxy_eta_pos_neg_ptw_array[n_harm][n_eta_gap][n_pos_neg][n_qxy][1]  += iQxy[n_qxy][n_harm]*p_t_weight;
                                        }
                                        if(PolarityA < 0.0) // use only negative charged particles for the event plane
                                        {
                                            EP_Qxy_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg][n_qxy][2]      += iQxy[n_qxy][n_harm]*total_weight;
                                            EP_Qxy_eta_pos_neg_ptw_array[n_harm][n_eta_gap][n_pos_neg][n_qxy][2]  += iQxy[n_qxy][n_harm]*p_t_weight;
                                        }
                                    }
                                    Qtracks_used_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg][0]++;
                                    if(PolarityA > 0.0) Qtracks_used_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg][1]++;
                                    if(PolarityA < 0.0) Qtracks_used_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg][2]++;
                                }
                            }
                        }
                        //------------------------------------------------------------------------



                        //------------------------------------------------------------------------
                        // Create array for the full TPC event plane resolution calculation
                        // For random sub events
                        if(N_EP_tracks < N_max_EP_track)
                        {
                            for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics
                            {
                                for(Int_t n_qxy = 0; n_qxy < 2; n_qxy++) // loop over Qx and Qy
                                {
                                    EP_Qxy_array_sub[n_harm][n_qxy][N_EP_tracks]     = iQxy[n_qxy][n_harm]*total_weight;
                                    EP_Qxy_ptw_array_sub[n_harm][n_qxy][N_EP_tracks] = iQxy[n_qxy][n_harm]*p_t_weight;
                                }
                            }
                            N_EP_tracks++;
                        }
                        //------------------------------------------------------------------------


                        Qtracks_used_full++;
                    }
                }
            }
        }
        // End of particle loop
        //----------------------------------------------------------------------------------------------------


        hFTPC_BBC_corr->Fill(N_FTPC,N_BBC);


        //------------------------------------------------------------------------
        // cout << "N_EP_tracks = " << N_EP_tracks << ", seed_number = " << seed_number << endl;
        // Int_t Event_ran_number = ran_number.Integer(run_events);  // will generate a random integer from 0..(run_events-1)
        // random generation of sub event A and sub event B
        // Event numbers -> A, odd events -> B
        // A random generated number selects the track, the Qx, Qy data array is changed by
        // copying the last value into the position of the used one
        Double_t EP_Qxy_subAB[n_harmonics][2][2];
        Double_t EP_Qxy_ptw_subAB[n_harmonics][2][2]; // [n harmonic][Qx, Qy][sub A, subB]
        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics
        {
            for(Int_t n_qxy = 0; n_qxy < 2; n_qxy++) // loop over Qx and Qy
            {
                for(Int_t n_pos_neg = 0; n_pos_neg < 2; n_pos_neg++) // loop over positive and negative eta hemisphere
                {
                    EP_Qxy_subAB[n_harm][n_qxy][n_pos_neg]     = 0.0;
                    EP_Qxy_ptw_subAB[n_harm][n_qxy][n_pos_neg] = 0.0;
                }
            }
        }

        for(Int_t r = 0; r < N_EP_tracks; r++)
        {
            Int_t track_ran_number = ran_number.Integer(N_EP_tracks-r);  // will generate a random integer from 0..(val-1)

            for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics
            {
                for(Int_t n_qxy = 0; n_qxy < 2; n_qxy++) // loop over Qx and Qy
                {
                    Double_t iQxy_add_array     = EP_Qxy_array_sub[n_harm][n_qxy][track_ran_number];
                    Double_t iQxy_ptw_add_array = EP_Qxy_ptw_array_sub[n_harm][n_qxy][track_ran_number];
                    EP_Qxy_array_sub[n_harm][n_qxy][track_ran_number]      = EP_Qxy_array_sub[n_harm][n_qxy][N_EP_tracks-r-1];
                    EP_Qxy_ptw_array_sub[n_harm][n_qxy][track_ran_number]  = EP_Qxy_ptw_array_sub[n_harm][n_qxy][N_EP_tracks-r-1];

                    Int_t n_pos_neg = 0; // 0 == positive sub, 1 == negative sub
                    if((r % 2) == 0) // even event -> used for sub event analysis
                    {n_pos_neg = 0;}
                    else {n_pos_neg = 1;}
                    EP_Qxy_subAB[n_harm][n_qxy][n_pos_neg]     += iQxy_add_array;
                    EP_Qxy_ptw_subAB[n_harm][n_qxy][n_pos_neg] += iQxy_ptw_add_array;
                }
            }
        }
        //------------------------------------------------------------------------



        //----------------------------------------------------------------------------------------------------
        // Write the output to different files (harmonics)

        // EP_Qxy_array[n_harm][n_qxy] // 2
        // EP_Qxy_ptw_array[n_harm][n_qxy] // 2
        // EP_Qxy_subAB[n_harm][n_qxy][n_pos_neg] // 4
        // EP_Qxy_ptw_subAB[n_harm][n_qxy][n_pos_neg] // 4
        // EP_Qxy_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg][n_qxy] // 24
        // EP_Qxy_eta_pos_neg_ptw_array[n_harm][n_eta_gap][n_pos_neg][n_qxy] // 24
        // Qtracks_used_eta_pos_neg_array[n_harm][n_eta_gap][n_pos_neg] 12

        for(Int_t n_harm = 0; n_harm < n_harmonics; n_harm++) // loop over different harmonics (5)
        {
            EventPlane_harm_NTDataArray[0]     =(Float_t)EventVertexX;
            EventPlane_harm_NTDataArray[1]     =(Float_t)EventVertexY;
            EventPlane_harm_NTDataArray[2]     =(Float_t)EventVertexZ;
            EventPlane_harm_NTDataArray[3]     =(Float_t)refMult;
            EventPlane_harm_NTDataArray[4]     =(Float_t)event_number;
            EventPlane_harm_NTDataArray[5]     =(Float_t)eventId;
            EventPlane_harm_NTDataArray[6]     =(Float_t)n_tofmatch_prim;
            EventPlane_harm_NTDataArray[7]     =(Float_t)n_non_primaries;
            EventPlane_harm_NTDataArray[8]     =(Float_t)RunId;
            EventPlane_harm_NTDataArray[9]     =(Float_t)ZDCx;
            EventPlane_harm_NTDataArray[10]    =(Float_t)BBCx;
            EventPlane_harm_NTDataArray[11]    =(Float_t)vzVpd;
            EventPlane_harm_NTDataArray[12]    =(Float_t)EP_Qxy_array[n_harm][0];        // Qx full TPC event plane, pT+phi-weight
            EventPlane_harm_NTDataArray[13]    =(Float_t)EP_Qxy_array[n_harm][1];        // Qy full TPC event plane, pT+phi-weight
            EventPlane_harm_NTDataArray[14]    =(Float_t)EP_Qxy_ptw_array[n_harm][0];    // Qx full TPC event plane, only pT-weight
            EventPlane_harm_NTDataArray[15]    =(Float_t)EP_Qxy_ptw_array[n_harm][1];    // Qy full TPC event plane, only pT-weight
            EventPlane_harm_NTDataArray[16]    =(Float_t)EP_Qxy_subAB[n_harm][0][0];     // Qx full TPC event plane, pT+phi-weight, sub event A
            EventPlane_harm_NTDataArray[17]    =(Float_t)EP_Qxy_subAB[n_harm][1][0];     // Qy full TPC event plane, pT+phi-weight, sub event A
            EventPlane_harm_NTDataArray[18]    =(Float_t)EP_Qxy_subAB[n_harm][0][1];     // Qx full TPC event plane, pT+phi-weight, sub event B
            EventPlane_harm_NTDataArray[19]    =(Float_t)EP_Qxy_subAB[n_harm][1][1];     // Qy full TPC event plane, pT+phi-weight, sub event B
            EventPlane_harm_NTDataArray[20]    =(Float_t)EP_Qxy_ptw_subAB[n_harm][0][0]; // Qx full TPC event plane, only pT-weight, sub event A
            EventPlane_harm_NTDataArray[21]    =(Float_t)EP_Qxy_ptw_subAB[n_harm][1][0]; // Qy full TPC event plane, only pT-weight, sub event A
            EventPlane_harm_NTDataArray[22]    =(Float_t)EP_Qxy_ptw_subAB[n_harm][0][1]; // Qx full TPC event plane, only pT-weight, sub event B
            EventPlane_harm_NTDataArray[23]    =(Float_t)EP_Qxy_ptw_subAB[n_harm][1][1]; // Qy full TPC event plane, only pT-weight, sub event B
            EventPlane_harm_NTDataArray[24]    =(Float_t)Qtracks_used_full; // Number of used tracks for full TPC event plane
            EventPlane_harm_NTDataArray[25]    =(Float_t)fQxBbcEastRaw[n_harm]; // Qx BBC East
            EventPlane_harm_NTDataArray[26]    =(Float_t)fQyBbcEastRaw[n_harm]; // Qy BBC East
            EventPlane_harm_NTDataArray[27]    =(Float_t)fQxBbcWestRaw[n_harm]; // Qx BBC West
            EventPlane_harm_NTDataArray[28]    =(Float_t)fQyBbcWestRaw[n_harm]; // Qy BBC West
            EventPlane_harm_NTDataArray[29]    =(Float_t)pmtSum[0]; // BBC hits East
            EventPlane_harm_NTDataArray[30]    =(Float_t)pmtSum[1]; // BBC hits West

            for(Int_t n_eta_gap = 0; n_eta_gap < n_eta_gap_values; n_eta_gap++)  // loop over different eta gap values (6)
            {
                for(Int_t n_charge = 0; n_charge < 3; n_charge++) // loop over charges: all, +, -
                {
                    EventPlane_harm_NTDataArray[31+n_eta_gap*30+n_charge*10]    =(Float_t)EP_Qxy_eta_pos_neg_array[n_harm][n_eta_gap][0][0][n_charge]; // Qx eta gap event plane, pT+phi-weight, pos eta
                    EventPlane_harm_NTDataArray[32+n_eta_gap*30+n_charge*10]    =(Float_t)EP_Qxy_eta_pos_neg_array[n_harm][n_eta_gap][0][1][n_charge]; // Qy eta gap event plane, pT+phi-weight, pos eta
                    EventPlane_harm_NTDataArray[33+n_eta_gap*30+n_charge*10]    =(Float_t)EP_Qxy_eta_pos_neg_array[n_harm][n_eta_gap][1][0][n_charge]; // Qx eta gap event plane, pT+phi-weight, neg eta
                    EventPlane_harm_NTDataArray[34+n_eta_gap*30+n_charge*10]    =(Float_t)EP_Qxy_eta_pos_neg_array[n_harm][n_eta_gap][1][1][n_charge]; // Qy eta gap event plane, pT+phi-weight, neg eta
                    EventPlane_harm_NTDataArray[35+n_eta_gap*30+n_charge*10]    =(Float_t)EP_Qxy_eta_pos_neg_ptw_array[n_harm][n_eta_gap][0][0][n_charge]; // Qx eta gap event plane, only pT-weight, pos eta
                    EventPlane_harm_NTDataArray[36+n_eta_gap*30+n_charge*10]    =(Float_t)EP_Qxy_eta_pos_neg_ptw_array[n_harm][n_eta_gap][0][1][n_charge]; // Qy eta gap event plane, only pT-weight, pos eta
                    EventPlane_harm_NTDataArray[37+n_eta_gap*30+n_charge*10]    =(Float_t)EP_Qxy_eta_pos_neg_ptw_array[n_harm][n_eta_gap][1][0][n_charge]; // Qx eta gap event plane, only pT-weight, neg eta
                    EventPlane_harm_NTDataArray[38+n_eta_gap*30+n_charge*10]    =(Float_t)EP_Qxy_eta_pos_neg_ptw_array[n_harm][n_eta_gap][1][1][n_charge]; // Qy eta gap event plane, only pT-weight, neg eta
                    EventPlane_harm_NTDataArray[39+n_eta_gap*30+n_charge*10]    =(Float_t)Qtracks_used_eta_pos_neg_array[n_harm][n_eta_gap][0][n_charge]; // Number of used tracks for positive eta gap
                    EventPlane_harm_NTDataArray[40+n_eta_gap*30+n_charge*10]    =(Float_t)Qtracks_used_eta_pos_neg_array[n_harm][n_eta_gap][1][n_charge]; // Number of used tracks for negative eta gap
                }
            }

            if(fAnalysisNum == 111)
            {
                EventPlane_array_NT[n_harm]->Fill(EventPlane_harm_NTDataArray);
            }
        }
        //----------------------------------------------------------------------------------------------------



        return 1;
    }
    else return 0;
}
//----------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------
void Get_day_file_id(Int_t runId, Int_t &day, Int_t &file_of_day, Int_t &file_id)
{
    // determine the file time
    TString name;
    Long64_t val;
    Float_t f_day  = 0.0;
    Float_t f_file_of_day = 0.0;
    char NoP[50];
    sprintf(NoP,"%u",(Int_t)runId);
    name = NoP[2];
    sscanf(name.Data(),"%Li",&val);
    f_day = (Float_t)(100.0 * val);
    name = NoP[3];
    sscanf(name.Data(),"%Li",&val);
    f_day += (Float_t)(10.0 * val);
    name = NoP[4];
    sscanf(name.Data(),"%Li",&val);
    f_day += (Float_t)(1.0 * val);

    name = NoP[5];
    sscanf(name.Data(),"%Li",&val);
    f_file_of_day = (Float_t)(100.0 * val);
    name = NoP[6];
    sscanf(name.Data(),"%Li",&val);
    f_file_of_day += (Float_t)(10.0 * val);
    name = NoP[7];
    sscanf(name.Data(),"%Li",&val);
    f_file_of_day += (Float_t)(1.0 * val);

    day         = (Int_t)f_day;
    file_of_day = (Int_t)f_file_of_day;
    file_id     = day*1000+file_of_day;
}
//----------------------------------------------------------------------------------------------------------------------------------------------



void Get_ReCentering_Correction(Int_t RunId, Int_t eBeamTimeNum, Int_t file_bin_rc, Float_t EventVertexZ, Float_t refMult, Float_t &rc_Qx_eta_pos, Float_t &rc_Qy_eta_pos, Float_t &rc_Qx_eta_neg,
                                Float_t &rc_Qy_eta_neg, Float_t &rc_Qx_full, Float_t &rc_Qy_full)
{
    // eBeamTimeNum:
    // 0 = 7.7
    // 1 = 11.5
    // 2 = 39
    // 3 = 62.4
    // 4 = 19.6
    // 5 = 27
    // 6 = 200



    //----------------------------------------------------------------------------------------------------
    const Double_t z_acceptance[N_Beamtime]   = {70.0,50.0,40.0,40.0,70.0,70.0,40.0};
    const Int_t n_z_vertex_bins = 10;
    Int_t z_bin = -1;
    Float_t start_z = -z_acceptance[eBeamTimeNum];
    Float_t delta_z = 2.0*z_acceptance[eBeamTimeNum]/((Float_t)n_z_vertex_bins);
    z_bin = (Int_t)((EventVertexZ-start_z)/delta_z); // z-vertex bin for re-centering correction
    //----------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------
    // Determine re-centering correction
    // h_rc_QxQy_etapm_z_vs_index_[q_id][z_bin][p_id]:  [Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full][z-bin][parameter]
    Float_t corr_params_rc[6][3]; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full, [p_id]: par0, par1, par2

    for(Int_t q_id = 0; q_id < 6; q_id++)  // q_id: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
    {
        for(Int_t p_id = 0; p_id < 3; p_id++) // p_id: par0, par1, par2
        {
            corr_params_rc[q_id][p_id] = h_rc_QxQy_etapm_z_vs_index_[q_id][z_bin][p_id]->GetBinContent(h_rc_QxQy_etapm_z_vs_index_[q_id][z_bin][p_id]->FindBin(file_bin_rc));
        }
    }


    // Correction parameters (always normalized to number of tracks/event)
    Float_t corr_Q_vectors_rc[6]; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
    for(Int_t q_id = 0; q_id < 6; q_id++)  // q_id: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
    {
        corr_Q_vectors_rc[q_id] = corr_params_rc[q_id][0] + corr_params_rc[q_id][1]*refMult + corr_params_rc[q_id][2]*refMult*refMult;
    }

    rc_Qx_eta_pos = corr_Q_vectors_rc[0];
    rc_Qy_eta_pos = corr_Q_vectors_rc[1];
    rc_Qx_eta_neg = corr_Q_vectors_rc[2];
    rc_Qy_eta_neg = corr_Q_vectors_rc[3];
    rc_Qx_full    = corr_Q_vectors_rc[4];
    rc_Qy_full    = corr_Q_vectors_rc[5];
    //----------------------------------------------------------------------------------------------------
}



Float_t Get_shift_Psi(Double_t f_Psi, Int_t EP_type, TH1F* f_hEP_shift_params[][n_shift_par], Int_t f_file_bin_rc)
{
    Float_t shift_Psi = 0.0;

    // EP_type
    // 0 = Full TPC event plane sub A, phi-weight method
    // 1 = Full TPC event plane sub B, phi-weight method
    // 2 = Full TPC event plane all, phi-weight method
    // 3 = Full TPC event plane sub A, re-centering method
    // 4 = Full TPC event plane sub B, re-centering method
    // 5 = Full TPC event plane all, re-centering method
    // 6 = Eta sub event plane eta plus, phi-weight method
    // 7 = Eta sub event plane eta neg, phi-weight method
    // 8 = Eta sub event plane eta plus, re-centering method
    // 9 = Eta sub event plane eta neg, re-centering method

    Float_t c2 = f_hEP_shift_params[EP_type][1]->GetBinContent(f_hEP_shift_params[EP_type][1]->FindBin(f_file_bin_rc));
    Float_t c4 = f_hEP_shift_params[EP_type][2]->GetBinContent(f_hEP_shift_params[EP_type][2]->FindBin(f_file_bin_rc));
    Float_t s2 = f_hEP_shift_params[EP_type][3]->GetBinContent(f_hEP_shift_params[EP_type][3]->FindBin(f_file_bin_rc));
    Float_t s4 = f_hEP_shift_params[EP_type][4]->GetBinContent(f_hEP_shift_params[EP_type][4]->FindBin(f_file_bin_rc));

    shift_Psi = -s2*TMath::Cos(2.0*f_Psi) + c2*TMath::Sin(2.0*f_Psi) + 0.5*(-s4*TMath::Cos(4.0*f_Psi) + c4*TMath::Sin(4.0*f_Psi));

    return shift_Psi;
}



Int_t EventPlane_Patrick_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                                  Int_t ParticleA, Int_t kBeamTimeNum, Int_t Ana_Num,Int_t run_events
                                 )
{
    // Event vertex information
    event_A_ana            = picoDst_A;
    StThreeVectorF pVertex = event_A_ana   ->primaryVertex();
    Float_t EventVertexX   = pVertex.x();
    Float_t EventVertexY   = pVertex.y();
    Float_t EventVertexZ   = pVertex.z();
    Int_t   refMult        = event_A_ana->refMult();
    Int_t   RunId          = event_A_ana->runId();
    //Int_t TriggerWord      = event_A_ana->triggerWord();
    //Float_t ZDCx           = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    //Float_t BBCx           = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    //Float_t vzVpd          = event_A_ana->vzVpd();



    //----------------------------------------------------------------------------------------------------
    Int_t day, file_of_day, file_id;
    Get_day_file_id((Int_t)RunId,day,file_of_day,file_id);
    Int_t file_bin_rc = (Int_t)(h_runId_index_rc_inv->GetBinContent(h_runId_index_rc_inv->FindBin(file_id))); // index bin for re-centering correction
    //----------------------------------------------------------------------------------------------------



    EP_Qx                = 0.0; // x value of event plane vector
    EP_Qy                = 0.0; // y value of event plane vector
    EP_Qx_eta_pos        = 0.0;   // sub event vectors for event plane resolution with pt and phi weight
    EP_Qy_eta_pos        = 0.0;
    EP_Qx_eta_neg        = 0.0;
    EP_Qy_eta_neg        = 0.0;
    EP_Qx_eta_pos_ptw    = 0.0;   // sub event vectors for event plane only with pt weight
    EP_Qy_eta_pos_ptw    = 0.0;
    EP_Qx_eta_neg_ptw    = 0.0;
    EP_Qy_eta_neg_ptw    = 0.0;
    EP_Qx_ptw            = 0.0;   // event plane vector only with pt weight
    EP_Qy_ptw            = 0.0;
    Qtracks_used_eta_pos = 0;
    Qtracks_used_eta_neg = 0;
    Qtracks_used         = 0;
    EP_Qx_subA_ptw       = 0.0;   // sub event vectors for event plane resolution
    EP_Qy_subA_ptw       = 0.0;
    EP_Qx_subB_ptw       = 0.0;
    EP_Qy_subB_ptw       = 0.0;

    for ( Int_t i = 0; i < 6; ++i )
    {
        EP_Q_vectors_rc_glob[i]  = 0;
        Qtracks_used_Arr_glob[i] = 0;
    }
    
    // reset all iQadd arrays to NULL pointers before EventPlane_Patrick_analysis()
    for(UShort_t i = 0; i < N_max_tracks; ++i) { // loop over possible maximum number of tracks
      if ( Arr_iQadd_FullTPC[i] ) Arr_iQadd_FullTPC[i]->Delete();
      if ( Arr_iQadd_EtaPos[i] ) Arr_iQadd_EtaPos[i]->Delete();
      if ( Arr_iQadd_EtaNeg[i] ) Arr_iQadd_EtaNeg[i]->Delete();
    }
    memset(&Arr_iQadd_FullTPC[0],NULL,sizeof(TVector2*)*N_max_tracks);
    memset(&Arr_iQadd_EtaPos[0],NULL,sizeof(TVector2*)*N_max_tracks);
    memset(&Arr_iQadd_EtaNeg[0],NULL,sizeof(TVector2*)*N_max_tracks);


    StThreeVectorF vector_prim, vector_prim_new;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    if(
       PID_counter_Array[Ana_Num][ParticleA] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && fabs(EventVertexZ) < vertex_z_cut
       && event_A_ana->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " particles *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA;
        StThreeVectorF vectorA, vectoratsA, vectorAB, vectorprimAB, vectornewA, vectornewB, vector_sum, dirA;
        StThreeVectorF testA, testAB, vectorABtoPrim, baseY, dirY;

        Double_t EP_Qx_phi_weight  = 0.0; // x value of event plane vector
        Double_t EP_Qy_phi_weight  = 0.0; // y value of event plane vector
        Double_t EP_Qx_no_weight   = 0.0; // x value of event plane vector
        Double_t EP_Qy_no_weight   = 0.0; // y value of event plane vector
        Double_t EP_Qx_subA        = 0.0;   // sub event vectors for event plane resolution
        Double_t EP_Qy_subA        = 0.0;
        Double_t EP_Qx_subB        = 0.0;
        Double_t EP_Qy_subB        = 0.0;
        Float_t phi_mean           = 0.0;
        Float_t mean_px            = 0.0;
        Float_t mean_py            = 0.0;
        Float_t mean_pz            = 0.0;
        Float_t mean_pxw           = 0.0;
        Float_t mean_pyw           = 0.0;
        Float_t mean_pzw           = 0.0;
        Float_t mean_py_pos        = 0.0;
        Float_t mean_py_neg        = 0.0;
        Float_t mean_px_pos        = 0.0;
        Float_t mean_px_neg        = 0.0;
        Int_t N_py_pos             = 0;
        Int_t N_py_neg             = 0;
        Int_t N_px_pos             = 0;
        Int_t N_px_neg             = 0;

        // For random sub event calculation
        TRandom ran_number;
        seed_number += run_events; // changes from event to event
        ran_number.SetSeed((UInt_t)seed_number);
        const Int_t N_max_EP_track = 3000;
        Double_t Qx_vector_array[N_max_EP_track];
        Double_t Qy_vector_array[N_max_EP_track];
        Double_t Qx_vector_array_ptw[N_max_EP_track];
        Double_t Qy_vector_array_ptw[N_max_EP_track];
        Int_t N_EP_tracks = 0;

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];

            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            StThreeVectorF vectorA = trackA.origin();

            // Requires: momentum, origin, signed Magnetic Field
            //           and Charge of particle (+/- 1)
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            //Float_t BetaA       = trackA.btofBeta();

            if(
               nHitsFitA     > nHitsFitA_EP_cut  // 14
               && nHitsPossA > nHitsPossA_EP_cut
               && MomentumA  > MomentumA_EP_low_cut
               && MomentumA  < MomentumA_EP_high_cut
              )
            {
                if(
                   ((Float_t)nHitsFitA)/((Float_t)nHitsPossA)    > nHits_ratio_EP_cut
                  )
                {

                    Float_t pathA = -999.0;
                    Float_t dcaAB = -999.0;

                    fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);

                    vectoratsA = helixA.cat(pathA);

                    Float_t eta = vectoratsA.pseudoRapidity();

                    if(
                       dcaA          < dcaAB_EP_cut  // 2.0
                       && fabs(eta)  < eta_EP_cut
                      )
                    {
                        dirA = helixA.cat(pathA); // direction vector
                        Double_t phiA = dirA.phi(); // phiA has values from -pi..pi


                        Float_t p_x = MomentumA*dirA.x();
                        Float_t p_y = MomentumA*dirA.y();
                        Float_t p_z = MomentumA*dirA.z();
                        Float_t p_t = sqrt(p_x*p_x + p_y*p_y);


                        mean_px = mean_px + p_x;
                        mean_py = mean_py + p_y;
                        mean_pz = mean_pz + p_z;

                        Float_t iQx = TMath::Cos(2.0*phiA);
                        Float_t iQy = TMath::Sin(2.0*phiA);

                        Double_t phi_w;
                        Double_t total_weight = calc_event_plane_weight(phiA,p_t,eta,RunId,EventVertexZ,(Int_t)PolarityA,phi_w);

                        Double_t p_t_weight = 1.0;
                        if(p_t < 2.0)  p_t_weight = p_t;
                        if(p_t >= 2.0) p_t_weight = 2.0;
                        Double_t phi_weight = 1.0;
                        if(p_t_weight > 0) phi_weight = total_weight/p_t_weight;


                        mean_pxw = mean_pxw + phi_weight*p_x;
                        mean_pyw = mean_pyw + phi_weight*p_y;
                        mean_pzw = mean_pzw + phi_weight*p_z;

                        if(p_x > 0.0)
                        {
                            N_px_pos++;
                            mean_px_pos += phi_weight*p_x;
                        }
                        if(p_x < 0.0)
                        {
                            N_px_neg++;
                            mean_px_neg += phi_weight*p_x;
                        }
                        if(p_y > 0.0)
                        {
                            N_py_pos++;
                            mean_py_pos += phi_weight*p_y;
                        }
                        if(p_y < 0.0)
                        {
                            N_py_neg++;
                            mean_py_neg += phi_weight*p_y;
                        }

                        //phi_mean = phi_mean + sign*weight_phi*phiA;
                        phi_mean = phi_mean + phi_weight*phiA;

                        EP_Qx_no_weight  += iQx;
                        EP_Qy_no_weight  += iQy;

                        EP_Qx_phi_weight += iQx*phi_weight;
                        EP_Qy_phi_weight += iQy*phi_weight;

                        Float_t iQx_add = total_weight*iQx;
                        Float_t iQy_add = total_weight*iQy;

                        EP_Qx    += iQx_add;
                        EP_Qy    += iQy_add;

                        EP_Qx_ptw += iQx*p_t_weight;
                        EP_Qy_ptw += iQy*p_t_weight;

                        // save array of iQx/y_add for later auto-correlation removal
			Arr_iQadd_FullTPC[trackA_num] = new TVector2(iQx*p_t_weight,iQy*p_t_weight);

                        // For eta sub-event method
                        if(
                           fabs(eta) > eta_gap
                          )
                        {
                            if(eta >= 0.0)
                            {
                                EP_Qx_eta_pos     += iQx_add;
                                EP_Qy_eta_pos     += iQy_add;
                                EP_Qx_eta_pos_ptw += iQx*p_t_weight;
                                EP_Qy_eta_pos_ptw += iQy*p_t_weight;
                                Arr_iQadd_EtaPos[trackA_num] = new TVector2(iQx*p_t_weight,iQy*p_t_weight);
                                Qtracks_used_eta_pos++;
                            }
                            if(eta < 0.0)
                            {
                                EP_Qx_eta_neg     += iQx_add;
                                EP_Qy_eta_neg     += iQy_add;
                                EP_Qx_eta_neg_ptw += iQx*p_t_weight;
                                EP_Qy_eta_neg_ptw += iQy*p_t_weight;
                                Arr_iQadd_EtaNeg[trackA_num] = new TVector2(iQx*p_t_weight,iQy*p_t_weight);
                                Qtracks_used_eta_neg++;
                            }
                        }


                        if(N_EP_tracks < N_max_EP_track)
                        {
                            Qx_vector_array[N_EP_tracks]     = iQx_add;
                            Qy_vector_array[N_EP_tracks]     = iQy_add;
                            Qx_vector_array_ptw[N_EP_tracks] = iQx*p_t_weight;
                            Qy_vector_array_ptw[N_EP_tracks] = iQy*p_t_weight;
                            N_EP_tracks++;
                        }

                        Qtracks_used++;
                    }
                }
            }
        }



        //----------------------------------------------------------------------------------------------------
        // Apply re-centering correction
        Float_t rc_Qx_eta_pos, rc_Qy_eta_pos, rc_Qx_eta_neg, rc_Qy_eta_neg, rc_Qx_full, rc_Qy_full;
        Get_ReCentering_Correction(RunId,kBeamTimeNum,file_bin_rc,EventVertexZ,(Float_t)refMult,rc_Qx_eta_pos,rc_Qy_eta_pos,rc_Qx_eta_neg,
                                   rc_Qy_eta_neg,rc_Qx_full,rc_Qy_full);

        Float_t corr_Q_vectors_rc[6]     = {rc_Qx_eta_pos,rc_Qy_eta_pos,rc_Qx_eta_neg,rc_Qy_eta_neg,rc_Qx_full,rc_Qy_full}; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
        Float_t EP_Q_vectors[6]          = {EP_Qx_eta_pos_ptw,EP_Qy_eta_pos_ptw,EP_Qx_eta_neg_ptw,EP_Qy_eta_neg_ptw,EP_Qx_ptw,EP_Qy_ptw};
        Float_t N_tracks_EP_Q_vectors[6] = {(Float_t)Qtracks_used_eta_pos,(Float_t)Qtracks_used_eta_pos,(Float_t)Qtracks_used_eta_neg,(Float_t)Qtracks_used_eta_neg,(Float_t)Qtracks_used,(Float_t)Qtracks_used};
        // Re-centering corrected Q-vectors
        Float_t EP_Q_vectors_rc[6]; // [q_id]: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full

        for(Int_t q_id = 0; q_id < 6; q_id++)  // q_id: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
        {
            if(N_tracks_EP_Q_vectors[q_id] > 0.0)
            {
                EP_Q_vectors_rc[q_id] = (EP_Q_vectors[q_id]/N_tracks_EP_Q_vectors[q_id]) - corr_Q_vectors_rc[q_id];
            }
            else
            {
                EP_Q_vectors_rc[q_id] = EP_Q_vectors[q_id];
            }
        }

        Float_t phi_event_plane_rc             = -100.0;
        Float_t phi_event_plane_eta_pos_gap_rc = -100.0;
        Float_t phi_event_plane_eta_neg_gap_rc = -100.0;

        if(N_tracks_EP_Q_vectors[0] > 0.0 && N_tracks_EP_Q_vectors[2] > 0.0)
        {
            phi_event_plane_eta_pos_gap_rc = calc_phi_event_plane_2nd(EP_Q_vectors_rc[0],EP_Q_vectors_rc[1]); //
            phi_event_plane_eta_neg_gap_rc = calc_phi_event_plane_2nd(EP_Q_vectors_rc[2],EP_Q_vectors_rc[3]); //
        }
        if(N_tracks_EP_Q_vectors[4] > 0.0 && N_tracks_EP_Q_vectors[5] > 0.0)
        {
            phi_event_plane_rc             = calc_phi_event_plane_2nd(EP_Q_vectors_rc[4],EP_Q_vectors_rc[5]); // subtract the Q-vector of the track to avoid auto correlations
        }
        //----------------------------------------------------------------------------------------------------



        // random generation of sub event A and sub event B
        // Event numbers -> A, odd events -> B
        // A random generated number selects the track, the Qx, Qy data array is changed by
        // copying the last value into the position of the used one
        for(Int_t r = 0; r < N_EP_tracks; r++)
        {
            Int_t track_ran_number = ran_number.Integer(N_EP_tracks-r);  // will generate a random integer from 0..(val-1)
            Float_t iQx_add     = Qx_vector_array[track_ran_number];
            Float_t iQy_add     = Qy_vector_array[track_ran_number];
            Float_t iQx_add_ptw = Qx_vector_array_ptw[track_ran_number];
            Float_t iQy_add_ptw = Qy_vector_array_ptw[track_ran_number];
            Qx_vector_array[track_ran_number] = Qx_vector_array[N_EP_tracks-r-1];
            Qy_vector_array[track_ran_number] = Qy_vector_array[N_EP_tracks-r-1];
            Qx_vector_array_ptw[track_ran_number] = Qx_vector_array_ptw[N_EP_tracks-r-1];
            Qy_vector_array_ptw[track_ran_number] = Qy_vector_array_ptw[N_EP_tracks-r-1];
            if((r % 2) == 0) // even event -> used for sub event analysis
            {
                EP_Qx_subA      += iQx_add;
                EP_Qy_subA      += iQy_add;
                EP_Qx_subA_ptw  += iQx_add_ptw;
                EP_Qy_subA_ptw  += iQy_add_ptw;
            }
            else
            {
                EP_Qx_subB      += iQx_add;
                EP_Qy_subB      += iQy_add;
                EP_Qx_subB_ptw  += iQx_add_ptw;
                EP_Qy_subB_ptw  += iQy_add_ptw;
            }
        }

        if(Qtracks_used > 0.0)
        {
            phi_mean = phi_mean/((Float_t)Qtracks_used);
            mean_px  = mean_px/((Float_t)Qtracks_used);
            mean_py  = mean_py/((Float_t)Qtracks_used);
            mean_pz  = mean_pz/((Float_t)Qtracks_used);
            mean_pxw = mean_pxw/((Float_t)Qtracks_used);
            mean_pyw = mean_pyw/((Float_t)Qtracks_used);
            mean_pzw = mean_pzw/((Float_t)Qtracks_used);
        }
        if(N_px_pos > 0)
        {
            mean_px_pos /=((Float_t)N_px_pos);
        }
        if(N_px_neg > 0)
        {
            mean_px_neg /=((Float_t)N_px_neg);
        }
        if(N_py_pos > 0)
        {
            mean_py_pos /=((Float_t)N_py_pos);
        }
        if(N_py_neg > 0)
        {
            mean_py_neg /=((Float_t)N_py_neg);
        }

        // Event plane angles with phi-weight correction
        Double_t Psi            = calc_phi_event_plane_2nd(EP_Qx,EP_Qy);
        Double_t Psi_eta_pos    = calc_phi_event_plane_2nd(EP_Qx_eta_pos,EP_Qy_eta_pos);
        Double_t Psi_eta_neg    = calc_phi_event_plane_2nd(EP_Qx_eta_neg,EP_Qy_eta_neg);
        //Double_t Psi_phi_weight = calc_phi_event_plane_2nd(EP_Qx_phi_weight,EP_Qy_phi_weight);
        //Double_t Psi_no_weight  = calc_phi_event_plane_2nd(EP_Qx_no_weight,EP_Qy_no_weight);



        //----------------------------------------------------------------------------------------------------
        // Determine shift correction
        Float_t shift_phi_event_plane                = Get_shift_Psi(Psi,2,hEP_shift_params,file_bin_rc); // Full TPC event plane all, phi-weight method
        Float_t shift_phi_event_plane_eta_gap_pos    = Get_shift_Psi(Psi_eta_pos,6,hEP_shift_params,file_bin_rc); // Eta sub event plane eta plus/neg, phi-weight method
        Float_t shift_phi_event_plane_eta_gap_neg    = Get_shift_Psi(Psi_eta_neg,7,hEP_shift_params,file_bin_rc); // Eta sub event plane eta plus/neg, phi-weight method
        Float_t shift_phi_event_plane_rc             = Get_shift_Psi(phi_event_plane_rc,5,hEP_shift_params,file_bin_rc); // Full TPC event plane all, re-centering method
        Float_t shift_phi_event_plane_eta_gap_rc_pos = Get_shift_Psi(phi_event_plane_eta_pos_gap_rc,6+2,hEP_shift_params,file_bin_rc); // Eta sub event plane eta plus/neg, re-centering method
        Float_t shift_phi_event_plane_eta_gap_rc_neg = Get_shift_Psi(phi_event_plane_eta_neg_gap_rc,7+2,hEP_shift_params,file_bin_rc); // Eta sub event plane eta plus/neg, re-centering method

        // Apply shift correction
        Psi                              += shift_phi_event_plane;
        Psi_eta_pos                      += shift_phi_event_plane_eta_gap_pos;
        Psi_eta_neg                      += shift_phi_event_plane_eta_gap_neg;
        phi_event_plane_rc               += shift_phi_event_plane_rc; // full TPC event plane angle with re-centering and shift correction
        phi_event_plane_eta_pos_gap_rc   += shift_phi_event_plane_eta_gap_rc_pos;
        phi_event_plane_eta_neg_gap_rc   += shift_phi_event_plane_eta_gap_rc_neg;
        //----------------------------------------------------------------------------------------------------

        for(Int_t q_id = 0; q_id < 6; q_id++)  // q_id: Qx-eta_pos,Qy-eta_pos,Qx-eta_neg,Qy-eta_neg,Qx_full,Qy_full
        {
            EP_Q_vectors_rc_glob[q_id]  = EP_Q_vectors_rc[q_id];
            Qtracks_used_Arr_glob[q_id] = N_tracks_EP_Q_vectors[q_id];
        }

        return 1;
    }
    else return 0;
}



void calc_event_plane_angles(Int_t nParticles,Float_t rap_c, StPicoAlexTrack trackA_c,StPicoAlexTrack trackB_c,StPicoAlexTrack trackC_c,StThreeVectorF vectornewA_c,StThreeVectorF vectornewB_c,StThreeVectorF vectornewC_c,Int_t ME_Flag_c, Int_t SE_ME_Flag_c,
                             Int_t RunIdA_c,Float_t EventVertexXA_c,Float_t EventVertexYA_c,Float_t EventVertexZA_c,Int_t RunIdB_c,Float_t EventVertexXB_c,Float_t EventVertexYB_c,Float_t EventVertexZB_c,
                             Float_t &phi_event_plane_c,Float_t &phi_event_plane_eta_gap_c, Float_t &delta_phi_ME_AB_weight_c, Float_t &delta_phi_ME_AB_weight_eta_gap_c)
{
    // Input: Only global track information!
    // nParticles is either 2 or 3, if 2 than put any information for third track
    // rap_c is the rapidity of the mother particle
    // This function calculates the event planes for the full TPC method -> phi_event_plane_c,
    // the corresponding difference between the angles of event A and B (mixed event) -> phi_event_plane_eta_gap_c.
    // The same is done for the eta gap method.

    phi_event_plane_c                = -400.0;
    phi_event_plane_eta_gap_c        = -400.0;
    delta_phi_ME_AB_weight_c         = 0.0;
    delta_phi_ME_AB_weight_eta_gap_c = 0.0;

    if(nParticles == 2 || nParticles == 3)
    {

        Float_t p_xA_c   = vectornewA_c.x();
        Float_t p_yA_c   = vectornewA_c.y();
        Float_t p_tA_c   = sqrt(p_xA_c*p_xA_c + p_yA_c*p_yA_c);
        Float_t etaA_c   = vectornewA_c.pseudoRapidity();

        Float_t p_xB_c   = vectornewB_c.x();
        Float_t p_yB_c   = vectornewB_c.y();
        Float_t p_tB_c   = sqrt(p_xB_c*p_xB_c + p_yB_c*p_yB_c);
        Float_t etaB_c   = vectornewB_c.pseudoRapidity();

        Float_t p_xC_c   = vectornewC_c.x();
        Float_t p_yC_c   = vectornewC_c.y();
        Float_t p_tC_c   = sqrt(p_xC_c*p_xC_c + p_yC_c*p_yC_c);
        Float_t etaC_c   = vectornewC_c.pseudoRapidity();

        Float_t phiA_c  = vectornewA_c.phi();
        Float_t phiB_c  = vectornewB_c.phi();
        Float_t phiC_c  = vectornewC_c.phi();

        Float_t iQxA   = TMath::Cos(2.0*phiA_c);
        Float_t iQyA   = TMath::Sin(2.0*phiA_c);
        Float_t iQxB   = TMath::Cos(2.0*phiB_c);
        Float_t iQyB   = TMath::Sin(2.0*phiB_c);
        Float_t iQxC   = TMath::Cos(2.0*phiC_c);
        Float_t iQyC   = TMath::Sin(2.0*phiC_c);

        Float_t PolarityA_c   = trackA_c.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
        Float_t PolarityB_c   = trackB_c.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
        Float_t PolarityC_c   = trackC_c.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle

        Double_t phi_w;
        Float_t total_weight_A  = calc_event_plane_weight(phiA_c,p_tA_c,etaA_c,RunIdA_c,EventVertexZA_c,(Int_t)PolarityA_c,phi_w);
        Float_t total_weight_eta_gap_A = total_weight_A;
        Float_t total_weight_B  = calc_event_plane_weight(phiB_c,p_tB_c,etaB_c,RunIdB_c,EventVertexZB_c,(Int_t)PolarityB_c,phi_w);
        Float_t total_weight_eta_gap_B = total_weight_B;
        Float_t total_weight_C  = calc_event_plane_weight(phiC_c,p_tC_c,etaC_c,RunIdA_c,EventVertexZA_c,(Int_t)PolarityC_c,phi_w);
        Float_t total_weight_eta_gap_C = total_weight_C;

        Float_t MomentumA_c   = trackA_c.gMom().mag();
        Float_t dcaA_c        = trackA_c.dca();   // distance of closest approach to primary vertex
        Float_t nHitsPossA_c  = trackA_c.nHitsMax();
        Float_t nHitsFitA_c   = trackA_c.nHitsFit();

        Float_t MomentumB_c   = trackB_c.gMom().mag();
        Float_t dcaB_c        = trackB_c.dca();   // distance of closest approach to primary vertex
        Float_t nHitsPossB_c  = trackB_c.nHitsMax();
        Float_t nHitsFitB_c   = trackB_c.nHitsFit();

        Float_t MomentumC_c   = trackC_c.gMom().mag();
        Float_t dcaC_c        = trackC_c.dca();   // distance of closest approach to primary vertex
        Float_t nHitsPossC_c  = trackC_c.nHitsMax();
        Float_t nHitsFitC_c   = trackC_c.nHitsFit();


        // check wether the tracks were used for the event plane calculation, else set weight to 0
        if(!(
             dcaA_c                        < dcaAB_EP_cut
             && fabs(etaA_c)               < eta_EP_cut
             && MomentumA_c                > MomentumA_EP_low_cut
             && MomentumA_c                < MomentumA_EP_high_cut
             && (nHitsFitA_c/nHitsPossA_c) > nHits_ratio_EP_cut
             && nHitsFitA_c                > nHitsFitA_EP_cut
             && nHitsPossA_c               > nHitsPossA_EP_cut
             && fabs(EventVertexZA_c)      < vertex_z_cut
             && (EventVertexXA_c*EventVertexXA_c + EventVertexYA_c*EventVertexYA_c) < Event_radius_cut*Event_radius_cut
            )
          )
        {
            total_weight_A         = 0.0;
            total_weight_eta_gap_A = 0.0;  // track was not used in eta gap EP calculation
        }
        if(
           !(fabs(etaA_c) > eta_gap)
          )
        {
            total_weight_eta_gap_A = 0.0;  // track was not used in eta gap EP calculation
        }

        if(!(
             dcaB_c                        < dcaAB_EP_cut
             && fabs(etaB_c)               < eta_EP_cut
             && MomentumB_c                > MomentumA_EP_low_cut
             && MomentumB_c                < MomentumA_EP_high_cut
             && (nHitsFitB_c/nHitsPossB_c) > nHits_ratio_EP_cut
             && nHitsFitB_c                > nHitsFitA_EP_cut
             && nHitsPossB_c               > nHitsPossA_EP_cut
             && fabs(EventVertexZB_c)      < vertex_z_cut
             && (EventVertexXB_c*EventVertexXB_c + EventVertexYB_c*EventVertexYB_c) < Event_radius_cut*Event_radius_cut
            )
          )
        {
            total_weight_B         = 0.0;
            total_weight_eta_gap_B = 0.0;  // track was not used in eta gap EP calculation
        }
        if(
           !(fabs(etaB_c) > eta_gap)
          )
        {
            total_weight_eta_gap_B = 0.0;  // track was not used in eta gap EP calculation
        }

        if(!(
             dcaC_c                        < dcaAB_EP_cut
             && fabs(etaC_c)               < eta_EP_cut
             && MomentumC_c                > MomentumA_EP_low_cut
             && MomentumC_c                < MomentumA_EP_high_cut
             && (nHitsFitC_c/nHitsPossC_c) > nHits_ratio_EP_cut
             && nHitsFitC_c                > nHitsFitA_EP_cut
             && nHitsPossC_c               > nHitsPossA_EP_cut
             && fabs(EventVertexZA_c)      < vertex_z_cut  // particle A and C are from one decay, e.g. Lambda -> p + pi-
             && (EventVertexXA_c*EventVertexXA_c + EventVertexYA_c*EventVertexYA_c) < Event_radius_cut*Event_radius_cut
            )
          )
        {
            total_weight_C         = 0.0;
            total_weight_eta_gap_C = 0.0;  // track was not used in eta gap EP calculation
        }
        if(
           !(fabs(etaC_c) > eta_gap)
          )
        {
            total_weight_eta_gap_C = 0.0;  // track was not used in eta gap EP calculation
        }

        if(nParticles == 2) // third track does not exist
        {
            total_weight_C         = 0.0;
            total_weight_eta_gap_C = 0.0;
        }


        // Subtract the Q-vector from the event plane vector to avoid auto correlation
        Float_t EP_Qx_new = EP_Qx - total_weight_A*iQxA - total_weight_B*iQxB - total_weight_C*iQxC;
        Float_t EP_Qy_new = EP_Qy - total_weight_A*iQyA - total_weight_B*iQyB - total_weight_C*iQyC;


        // The same for eta gap method
        Float_t EP_Qx_new_eta_gap   = 0;
        Float_t EP_Qy_new_eta_gap   = 0;

        Float_t EP_Qx_new_eta_gap_A = 0; // for mixed event
        Float_t EP_Qy_new_eta_gap_A = 0;
        Float_t EP_Qx_new_eta_gap_B = 0;
        Float_t EP_Qy_new_eta_gap_B = 0;

        if(rap_c < 0.0) // Use positive eta EP for netagive eta particle track
        {
            EP_Qx_new_eta_gap   = EP_Qx_eta_pos;
            EP_Qy_new_eta_gap   = EP_Qy_eta_pos;
            EP_Qx_new_eta_gap_A = EP_Qx_eta_pos; // for mixed event
            EP_Qy_new_eta_gap_A = EP_Qy_eta_pos;
            EP_Qx_new_eta_gap_B = EP_Qx_eta_pos_B;
            EP_Qy_new_eta_gap_B = EP_Qy_eta_pos_B;
        }
        if(rap_c > 0.0) // Use negative eta EP for positive eta particle track
        {
            EP_Qx_new_eta_gap   = EP_Qx_eta_neg;
            EP_Qy_new_eta_gap   = EP_Qy_eta_neg;
            EP_Qx_new_eta_gap_A = EP_Qx_eta_neg; // for mixed event
            EP_Qy_new_eta_gap_A = EP_Qy_eta_neg;
            EP_Qx_new_eta_gap_B = EP_Qx_eta_neg_B;
            EP_Qy_new_eta_gap_B = EP_Qy_eta_neg_B;
        }

        if((rap_c < 0.0 && etaA_c > 0.0) || (rap_c > 0.0 && etaA_c < 0.0))
        {
            EP_Qx_new_eta_gap = EP_Qx_new_eta_gap - total_weight_eta_gap_A*iQxA;
            EP_Qy_new_eta_gap = EP_Qy_new_eta_gap - total_weight_eta_gap_A*iQyA;
        }
        else total_weight_eta_gap_A = 0.0;

        if((rap_c < 0.0 && etaB_c > 0.0) || (rap_c > 0.0 && etaB_c < 0.0))
        {
            EP_Qx_new_eta_gap = EP_Qx_new_eta_gap - total_weight_eta_gap_B*iQxB;
            EP_Qy_new_eta_gap = EP_Qy_new_eta_gap - total_weight_eta_gap_B*iQyB;
        }
        else total_weight_eta_gap_B = 0.0;

        if((rap_c < 0.0 && etaC_c > 0.0) || (rap_c > 0.0 && etaC_c < 0.0))
        {
            EP_Qx_new_eta_gap = EP_Qx_new_eta_gap - total_weight_eta_gap_C*iQxC;
            EP_Qy_new_eta_gap = EP_Qy_new_eta_gap - total_weight_eta_gap_C*iQyC;
        }
        else total_weight_eta_gap_C = 0.0;

        phi_event_plane_c         = -400.0;
        phi_event_plane_eta_gap_c = -400.0;
        phi_event_plane_c         = calc_phi_event_plane_2nd(EP_Qx_new,EP_Qy_new);
        phi_event_plane_eta_gap_c = calc_phi_event_plane_2nd(EP_Qx_new_eta_gap,EP_Qy_new_eta_gap);
        //******************************************************************************************************


        //***************************** Event plane for mixed event analysis ***********************************
        delta_phi_ME_AB_weight_c         = 0.0;
        delta_phi_ME_AB_weight_eta_gap_c = 0.0;
        if(ME_Flag_c == 1 && SE_ME_Flag_c == 1)
        {
            // Full EP
            Float_t phi_event_planeA = -400.0;
            phi_event_planeA = calc_phi_event_plane_2nd(EP_Qx - total_weight_A*iQxA - total_weight_C*iQxC,EP_Qy - total_weight_A*iQyA - total_weight_C*iQyC);

            Float_t phi_event_planeB = -400.0;
            phi_event_planeB = calc_phi_event_plane_2nd(EP_Qx_B - total_weight_B*iQxB,EP_Qy_B - total_weight_B*iQyB);

            delta_phi_ME_AB_weight_c = phi_event_planeA - phi_event_planeB; // difference between the two event planes

            // mean value of the two event planes (full EP)
            if(fabs(delta_phi_ME_AB_weight_c) < TMath::Pi()/2.0)
            {
                phi_event_plane_c = (phi_event_planeA + phi_event_planeB)/2.0;
            }
            else
            {
                phi_event_plane_c = (phi_event_planeA + (phi_event_planeB+TMath::Pi()))/2.0;
            }
            if(fabs(phi_event_plane_c) > TMath::Pi()/2.0 && fabs(phi_event_plane_c) < 1.5*TMath::Pi())
            {
                phi_event_plane_c = phi_event_plane_c-TMath::Pi();
            }
            // difference between the two angles can not be larger than pi/2
            if(fabs(delta_phi_ME_AB_weight_c) > TMath::Pi()/2.0)
            {
                delta_phi_ME_AB_weight_c = TMath::Pi() - fabs(delta_phi_ME_AB_weight_c);
            }



            // eta gap EP
            Float_t phi_event_plane_eta_gap_A = -400.0;
            phi_event_plane_eta_gap_A = calc_phi_event_plane_2nd(EP_Qx_new_eta_gap_A - total_weight_eta_gap_A*iQxA - total_weight_eta_gap_C*iQxC,EP_Qy_new_eta_gap_A - total_weight_eta_gap_A*iQyA - total_weight_eta_gap_C*iQyC);
            Float_t phi_event_plane_eta_gap_B = -400.0;
            phi_event_plane_eta_gap_B = calc_phi_event_plane_2nd(EP_Qx_new_eta_gap_B - total_weight_eta_gap_B*iQxB,EP_Qy_new_eta_gap_B - total_weight_eta_gap_B*iQyB);

            delta_phi_ME_AB_weight_eta_gap_c = phi_event_plane_eta_gap_A - phi_event_plane_eta_gap_B; // difference between the two event planes

            // mean value of the two event planes (eta gap)
            if(fabs(delta_phi_ME_AB_weight_eta_gap_c) < TMath::Pi()/2.0)
            {
                phi_event_plane_eta_gap_c = (phi_event_plane_eta_gap_A + phi_event_plane_eta_gap_B)/2.0;
            }
            else
            {
                phi_event_plane_eta_gap_c = (phi_event_plane_eta_gap_A + (phi_event_plane_eta_gap_B+TMath::Pi()))/2.0;
            }
            if(fabs(phi_event_plane_eta_gap_c) > TMath::Pi()/2.0 && fabs(phi_event_plane_eta_gap_c) < 1.5*TMath::Pi())
            {
                phi_event_plane_eta_gap_c = phi_event_plane_eta_gap_c-TMath::Pi();
            }
            // difference between the two angles can not be larger than pi/2
            if(fabs(delta_phi_ME_AB_weight_eta_gap_c) > TMath::Pi()/2.0)
            {
                delta_phi_ME_AB_weight_eta_gap_c = TMath::Pi() - fabs(delta_phi_ME_AB_weight_eta_gap_c);
            }


            //cout << "phi_event_planeA = " << phi_event_planeA*TMath::RadToDeg() << ", phi_event_planeB = " << phi_event_planeB*TMath::RadToDeg()
            //    << ", phi_event_plane = " << phi_event_plane*TMath::RadToDeg() << ", delta_phi_ME_AB_weight = " << delta_phi_ME_AB_weight*TMath::RadToDeg() << endl;
        }
    }
    else cout << "WARNING! nParticles was not set correctly in calc_event_plane_angles function" << endl;

}


Int_t MyLeptonTreeAna(Int_t PID_counter_Array[][N_max_PIDs],Int_t PID_Array[][N_max_PIDs][N_max_tracks],
                      StPicoAlexEvent* picoDst_A, Int_t Ana_Num)
{
    Float_t EventVertexXA,EventVertexYA,EventVertexZA;
    Int_t refMultA;

    event_A_ana      = picoDst_A;
    EventVertexXA    = event_A_ana->primaryVertex().x();
    EventVertexYA    = event_A_ana->primaryVertex().y();
    EventVertexZA    = event_A_ana->primaryVertex().z();
    refMultA         = event_A_ana->refMult();
    Int_t   RunIdA   = event_A_ana->runId();
    //Float_t ZDCx     = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    //Float_t BBCx     = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    //Float_t vzVpd    = event_A_ana->vzVpd();

    Float_t radius_cut           = Event_radius_cut*Event_radius_cut; // 2.0 cm radius cut for good events

    if(
       event_A_ana   ->isMinBias()
       && fabs(EventVertexZA) < vertex_z_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
      )
    {
        my_event->clearTrackList();
        my_event->setRunId  (RunIdA);
        my_event->setEventId(event_A_ana->eventId());
        my_event->setVertexX(EventVertexXA);
        my_event->setVertexY(EventVertexYA);
        my_event->setVertexZ(EventVertexZA);
        my_event->setrefMult(refMultA);
        my_event->setNTriggerIds( 0 );
        my_event->setTriggerId  ( 0 );
        my_event->setTriggerId2 ( 0 );
        my_event->setTriggerId3 ( 0 );
        my_event->setTriggerId4 ( 0 );
        my_event->setTriggerId5 ( 0 );
        my_event->setVertex_frac( 0 );
        my_event->setNVertices  ( 0 );
        my_event->setBBC_delta_t( 0 );
        my_event->setVPD_delta_t( 0 );
        my_event->set_EP_Q_vectors_rc(EP_Q_vectors_rc_glob);
        my_event->set_Qtracks_used_Arr(Qtracks_used_Arr_glob);

        const Int_t NrPart        = 2;
        Int_t ParticleNrs[NrPart] = {2,3};

        for(Int_t nr = 0; nr < NrPart; ++nr)
        {
            if(PID_counter_Array[Ana_Num][ParticleNrs[nr]] > 0)
            {
                StParticleTrack* mytrackA;

                for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleNrs[nr]]; i++)
                { // particle candidates
                    Int_t trackA_num = PID_Array[Ana_Num][ParticleNrs[nr]][i];
                    StPicoAlexTrack trackA = *picoDst_A->track( trackA_num);

                    mytrackA = my_event->createTrack();

                    Double_t Mass2      = -100.0;
                    if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && trackA.btofBeta() != 0)
                    {
                        Mass2 = trackA.pMom().mag()*trackA.pMom().mag()*(1.0/(trackA.btofBeta()*trackA.btofBeta()) - 1.0);
                    }

                    StPhysicalHelixD helixA_prim, helixA_glob;
                    helixA_prim   = StPhysicalHelixD(trackA.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());
                    helixA_glob   = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());


                    mytrackA->setTPCdEdx       (trackA.dEdx());
                    mytrackA->setCharge        (trackA.charge       ());
                    mytrackA->setBeta          (trackA.btofBeta          ());
                    mytrackA->setMomentum      (trackA.gMom().mag       ());
                    mytrackA->setMass2         (Mass2);
                    mytrackA->setnSigmaEl      (trackA.nSigmaElectron       ()*nsigma_scaling_fac);
                    mytrackA->setnSigmaPion    (trackA.nSigmaPion     ()*nsigma_scaling_fac);
                    mytrackA->setnSigmaKaon    (trackA.nSigmaKaon     ()*nsigma_scaling_fac);
                    mytrackA->setnSigmaP       (trackA.nSigmaProton        ()*nsigma_scaling_fac);
                    mytrackA->setDCA           (trackA.dca            ());
                    mytrackA->setnHitsPoss     (trackA.nHitsMax      ());
                    mytrackA->setnHitsFit      (trackA.nHitsFit       ());
                    mytrackA->setnHitsdEdx     (trackA.nHitsDedx      ());
                    mytrackA->setdipAngle      (helixA_glob.dipAngle       ());
                    mytrackA->setcurvature     (helixA_glob.curvature      ());
                    mytrackA->setphase         (helixA_glob.phase          ());
                    mytrackA->seth             (helixA_glob.h              ());
                    mytrackA->setoriginX       (helixA_glob.origin().x        ());
                    mytrackA->setoriginY       (helixA_glob.origin().y        ());
                    mytrackA->setoriginZ       (helixA_glob.origin().z        ());
                    mytrackA->setdipAnglePrim  (helixA_prim.dipAngle   ());
                    mytrackA->setcurvaturePrim (helixA_prim.curvature  ());
                    mytrackA->setphasePrim     (helixA_prim.phase      ());
                    mytrackA->sethPrim         (helixA_prim.h          ());
                    mytrackA->setoriginXPrim   (helixA_prim.origin().x    ());
                    mytrackA->setoriginYPrim   (helixA_prim.origin().y    ());
                    mytrackA->setoriginZPrim   (helixA_prim.origin().z    ());
                    mytrackA->setMomentumPrim  (trackA.pMom().mag   ());
                    mytrackA->setYLocal        (trackA.btofYLocal         ());

                    if ( Arr_iQadd_FullTPC[trackA_num] ) {
                        mytrackA->setFlagFullTPC(1);
                        mytrackA->setQVecFullTPC(Arr_iQadd_FullTPC[trackA_num]);
                    }
                    if ( Arr_iQadd_EtaPos[trackA_num] ) {
                        mytrackA->setFlagEtaPos(1);
                        mytrackA->setQVecEtaPos(Arr_iQadd_EtaPos[trackA_num]);
                    }
                    if ( Arr_iQadd_EtaNeg[trackA_num] ) {
                        mytrackA->setFlagEtaNeg(1);
                        mytrackA->setQVecEtaNeg(Arr_iQadd_EtaNeg[trackA_num]);
                    }

                }

            }
        }

        my_tree->Fill();
        //my_tree->Show(good_event);
        my_event->Clear();

        return 1;
    }
    else
    {
        return 0;
    }
}


Int_t SingleTrackTree(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A, MyEventEmbData* my_event_single, Int_t ParticleA, Int_t ParticleB, Int_t Ana_Num)
{
    // Patricks function to generate single lepton trees
    if ( PID_counter_Array[Ana_Num][ParticleA] > 0 || PID_counter_Array[Ana_Num][ParticleB] > 0 )
    {
        event_A_ana      = picoDst_A;

        StThreeVectorF vectorprim;
        vectorprim.set(event_A_ana->primaryVertex().x(),event_A_ana->primaryVertex().y(),event_A_ana->primaryVertex().z());

        StPhysicalHelixD helixA_Prim, helixB_Prim;
        TLorentzVector ltrackA_Prim, ltrackB_Prim;
        StThreeVectorF vectorA_Prim, primdirA, vectorB_Prim, primdirB;
        StPhysicalHelixD helixA_Glob, helixB_Glob;
        TLorentzVector ltrackA_Glob, ltrackB_Glob, ltrackAB;
        Float_t pathA_f, pathB_f,pathA_Est,pathB_Est;
        Float_t dcaAB_f, dcaPhPrim, pL;
        StThreeVectorF vectorA_Glob, vectorB_Glob, vectoratsA, vectoratsB, vectorAB, phvec;
        StThreeVectorF primdirA_Glob, primdirB_Glob, decvec;

        Bool_t conv_flag;
        MyTrackEmbData* mytrackA;

        Int_t nInvMass;
        if ( ParticleA == 3 && ParticleB == 2 ) nInvMass = 0; // unlike s
        if ( ParticleA == 3 && ParticleB == 3 ) nInvMass = 1; // neg ls
        if ( ParticleA == 2 && ParticleB == 2 ) nInvMass = 2; // pos ls

        // Loop over all electrons
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++)
        { // e- candidates

            conv_flag = kFALSE;

            // set trackA
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num);

            // primary helix & lorentzvector trackA
            helixA_Prim   = StPhysicalHelixD(trackA.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

            Float_t MomentumA   = trackA.pMom().mag();
            primdirA = helixA_Prim.cat(helixA_Prim.pathLength(vectorprim)); //
            primdirA = MomentumA*primdirA/primdirA.mag();
            ltrackA_Prim.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),0.00051099892);

            // e- pt cut
            if ( ltrackA_Prim.Pt() < 0.2 || ltrackA_Prim.Pt() > 10 ) continue;

            // determine vertex
            if ( PID_counter_Array[Ana_Num][ParticleA] > 0 && PID_counter_Array[Ana_Num][ParticleB] > 0 )
            {
                for(Int_t j = 0; j < PID_counter_Array[Ana_Num][ParticleB]; j++)
                { // e+ candidates
                    // avoid double counting ( double vertex determination! )
                    if ( ParticleA != ParticleB || ( ParticleA == ParticleB && j > i ) )
                    {
                        Int_t trackB_num = PID_Array[Ana_Num][ParticleB][j];
                        if( trackA_num != trackB_num )
                        {
                            StPicoAlexTrack trackB = *picoDst_A->track( trackB_num);

                            // global helices
                            helixA_Glob   = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());
                            Float_t MomentumA_Glob   = trackA.gMom().mag();

                            helixB_Glob   = StPhysicalHelixD(trackB.gMom(),trackB.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackB.charge());
                            Float_t MomentumB_Glob   = trackB.gMom().mag();

                            // calculate decay vertex
                            Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA_Glob,helixB_Glob,pathA_Est,pathB_Est,dcaAB_f);
                            Float_t dcaAB_est = dcaAB_f;
                            if ( dcaAB_est >= 2 ) continue;
                            if(fDCA_Helix_out == 1) fHelixABdca_start_params(helixA_Glob,helixB_Glob,pathA_f,pathB_f,dcaAB_f,pathA_Est,pathB_Est);
                            else fHelixABdca(helixA_Glob,helixB_Glob,pathA_f,pathB_f,dcaAB_f);// calculate dca between two helices
                            if ( dcaAB_f >= 1 ) continue;
                            hDcaDiff[nInvMass]->Fill(dcaAB_est);
                            vectoratsA     = helixA_Glob.at(pathA_f);  // space vector of helixA at dca to helixB
                            vectoratsB     = helixB_Glob.at(pathB_f);  // space vector of helixB at dca to helixA
                            vectorAB       = vectoratsA+vectoratsB;
                            vectorAB       = vectorAB/2.0; // decay vertex
                            // global lorentzvectors at decay vertex
                            primdirA_Glob = helixA_Glob.cat(helixA_Glob.pathLength(vectorAB)); //
                            primdirA_Glob = MomentumA_Glob*primdirA_Glob/primdirA_Glob.mag();
                            ltrackA_Glob.SetXYZM(primdirA_Glob.x(),primdirA_Glob.y(),primdirA_Glob.z(),0.00051099892);
                            primdirB_Glob = helixB_Glob.cat(helixB_Glob.pathLength(vectorAB)); //
                            primdirB_Glob = MomentumB_Glob*primdirB_Glob/primdirB_Glob.mag();
                            ltrackB_Glob.SetXYZM(primdirB_Glob.x(),primdirB_Glob.y(),primdirB_Glob.z(),0.00051099892);
                            // momenta
                            ltrackAB = ltrackA_Glob+ltrackB_Glob;
                            phvec.set(ltrackAB.Px(),ltrackAB.Py(),ltrackAB.Pz());
                            decvec = vectorAB-vectorprim;
                            pL = phvec.dot(decvec);
                            if ( pL <= 0 ) continue;
                            dcaPhPrim = calculateMinimumDistanceStraightToPoint(vectorAB,phvec,vectorprim);
                            if ( dcaPhPrim < 3 ) {
                                if ( ltrackAB.M() < 0.25 ) hInvMassAB[nInvMass]->Fill(ltrackAB.M());
                                if ( ltrackAB.M() < 0.005 ) {
                                    conv_flag = kTRUE;
                                    if ( ParticleA != ParticleB ) {
                                        // fill decay vertex histos
                                        hDecVertex_YX->Fill(vectorAB.x(),vectorAB.y());
                                        hDecVertex_YZ->Fill(vectorAB.z(),vectorAB.y());
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if ( conv_flag ) hnSigElDist[nInvMass+1]->Fill(trackA.nSigmaElectron()*nsigma_scaling_fac);

            if ( ParticleA != ParticleB && conv_flag )
            {
                Double_t Mass2      = -100.0;
                if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && trackA.btofBeta() != 0)
                {
                    Mass2 = trackA.pMom().mag()*trackA.pMom().mag()*(1.0/(trackA.btofBeta()*trackA.btofBeta()) - 1.0);
                }

                mytrackA = my_event_single->AddElectronTrack();
                // variables
                mytrackA->SetM2       ((Float_t)Mass2);  // Squared mass
                mytrackA->SetnSigmaEl ((Float_t)trackA.nSigmaElectron       ()*nsigma_scaling_fac);
                mytrackA->SetMomentum ((Float_t)trackA.charge()*(Float_t)MomentumA);
                mytrackA->SetBeta     ((Float_t)trackA.btofBeta          ());  // Velocity after time-of-flight reconstruction
                mytrackA->Setdca      ((Float_t)trackA.dca            ());
                mytrackA->Setpt       ((Float_t)ltrackA_Prim.Pt());
                mytrackA->SetTPCdEdx  ((Float_t)trackA.dEdx());
                if(ltrackA_Prim.Pt() > 0)
                {
                    mytrackA->Seteta      ((Float_t)ltrackA_Prim.PseudoRapidity());
                }
                else
                {
                    mytrackA->Seteta      (1e20);
                }
                mytrackA->Sety        ((Float_t)ltrackA_Prim.Rapidity());
                mytrackA->Setphi      ((Float_t)ltrackA_Prim.Phi());
                mytrackA->SetnHitsPoss((Float_t)trackA.nHitsMax      ());
                mytrackA->SetnHitsFit ((Float_t)trackA.nHitsFit       ());
                mytrackA->SetnHitsdEdx((Float_t)trackA.nHitsDedx      ());
                mytrackA->SetyLocal   ((Float_t)trackA.btofYLocal         ());
            }
        }

    }
    return 1;
}

Int_t PiKP_SingleTrackTree(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A, MyEventEmbData* my_event_single, Int_t Ana_Num)
{

    // Patricks function to generate single lepton trees
    const Int_t NrPart = 6;
    Int_t ParticleNrs[NrPart] = {8,9,11,12,14,15};
    Float_t masses[NrPart] = {0.13957018,0.13957018,0.493677,0.493677,0.938272013,0.938272013};

    event_A_ana      = picoDst_A;

    StThreeVectorF vectorprim;
    vectorprim.set(event_A_ana->primaryVertex().x(),event_A_ana->primaryVertex().y(),event_A_ana->primaryVertex().z());

    StPhysicalHelixD helixA_Prim;
    TLorentzVector ltrackA_Prim;
    StThreeVectorF vectorA_Prim, primdirA;

    for ( Int_t nr = 0; nr < NrPart; ++nr ) {
        if ( PID_counter_Array[Ana_Num][ParticleNrs[nr]] > 0 )
        {
            //cout << " # particles [" << ParticleNrs[nr] << "] = " <<  PID_counter_Array[Ana_Num][ParticleNrs[nr]] << endl;
            MyTrackEmbData* mytrackA;
            for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleNrs[nr]]; i++)
            { // particle candidates
                Int_t trackA_num = PID_Array[Ana_Num][ParticleNrs[nr]][i];

                StPicoAlexTrack trackA = *picoDst_A->track( trackA_num);

                // primary helix & lorentzvector trackA
                helixA_Prim   = StPhysicalHelixD(trackA.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());
                Float_t MomentumA   = trackA.pMom().mag();

                primdirA = helixA_Prim.cat(helixA_Prim.pathLength(vectorprim));
                primdirA = MomentumA*primdirA/primdirA.mag();
                ltrackA_Prim.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),masses[nr]);
                // pt cut
                if ( ltrackA_Prim.Pt() < 0.2 || ltrackA_Prim.Pt() > 10 ) continue;

                if ( ParticleNrs[nr] == 8  || ParticleNrs[nr] == 9  ) mytrackA = my_event_single->AddPionTrack();
                if ( ParticleNrs[nr] == 11 || ParticleNrs[nr] == 12 ) mytrackA = my_event_single->AddKaonTrack();
                if ( ParticleNrs[nr] == 14 || ParticleNrs[nr] == 15 ) mytrackA = my_event_single->AddProtonTrack();
                // variables
                Double_t Mass2      = -100.0;
                if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && trackA.btofBeta() != 0)
                {
                    Mass2 = trackA.pMom().mag()*trackA.pMom().mag()*(1.0/(trackA.btofBeta()*trackA.btofBeta()) - 1.0);
                }

                mytrackA->SetM2       ((Float_t)Mass2);  // Squared mass
                mytrackA->SetnSigmaEl ((Float_t)trackA.nSigmaElectron       ()*nsigma_scaling_fac);
                mytrackA->SetMomentum ((Float_t)trackA.charge()*(Float_t)MomentumA);
                mytrackA->SetBeta     ((Float_t)trackA.btofBeta          ());  // Velocity after time-of-flight reconstruction
                mytrackA->Setdca      ((Float_t)trackA.dca            ());
                mytrackA->Setpt       ((Float_t)ltrackA_Prim.Pt());
                mytrackA->SetTPCdEdx  ((Float_t)trackA.dEdx());
                if(ltrackA_Prim.Pt() > 0)
                {
                    mytrackA->Seteta      ((Float_t)ltrackA_Prim.PseudoRapidity());
                }
                else
                {
                    mytrackA->Seteta      (1e20);
                }
                mytrackA->Sety        ((Float_t)ltrackA_Prim.Rapidity());
                mytrackA->Setphi      ((Float_t)ltrackA_Prim.Phi());
                mytrackA->SetnHitsPoss((Float_t)trackA.nHitsMax      ());
                mytrackA->SetnHitsFit ((Float_t)trackA.nHitsFit       ());
                mytrackA->SetnHitsdEdx((Float_t)trackA.nHitsDedx      ());
                mytrackA->SetyLocal   ((Float_t)trackA.btofYLocal         ());
            }

        }
    }
#if 0
    cout << "==> " << my_event_single->GetNrPions() << endl;
    cout << "==> " << my_event_single->GetNrKaons() << endl;
    cout << "==> " << my_event_single->GetNrProtons() << endl;
#endif

    return 1;
}



Int_t PhiKPKM_analysis(Int_t PID_counter_Array[][N_max_PIDs],Int_t PID_counter_Array_B[][N_max_PIDs],
                       Int_t PID_Array[][N_max_PIDs][N_max_tracks],Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                       StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                       Int_t ParticleA, Int_t ParticleB, Int_t Ana_Num,
                       Int_t run_events,Int_t run_events_B,Int_t event_number, Int_t SE_ME_Flag)
{
    // Au + Au -> phi + N + N
    //             |
    //             -> K- + K+

    // Event vertex information
    StThreeVectorF vector_prim,vector_primB,vectordiff;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB;

    event_A_ana   = picoDst_A;
    event_B_ana   = picoDst_B;

    EventVertexXA    = event_A_ana   ->primaryVertex().x();
    EventVertexYA    = event_A_ana   ->primaryVertex().y();
    EventVertexZA    = event_A_ana   ->primaryVertex().z();
    EventVertexXB    = EventVertexXA;
    EventVertexYB    = EventVertexYA;
    EventVertexZB    = EventVertexZA;
    refMultA         = event_A_ana   ->refMult();
    refMultB         = refMultA;
    Int_t   RunIdA   = event_A_ana   ->runId();
    Int_t   RunIdB   = event_B_ana   ->runId();
    Float_t ZDCx     = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx     = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd    = event_A_ana->vzVpd();
    vector_prim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    //cout << "Entered phi analysis, vertex = {" << EventVertexXA << ", " << EventVertexYA << ", " << EventVertexZA << "}" << endl;

    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        //cout << "Start mixing, check if events can be mixed..." << endl;
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }

    Float_t radius_cut           = Event_radius_cut*Event_radius_cut; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
    }


    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]   > 0   // K+
       && PID_counter_Array_B[Ana_Num][ParticleB] > 0   // K-
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       //&& event_A_ana   ->isMinBias()
       //&& event_B_ana   ->isMinBias()
      )
    {
        alexPhiMeson_event.clearTrackList();
        alexPhiMeson_event.setx(EventVertexXA);
        alexPhiMeson_event.sety(EventVertexYA);
        alexPhiMeson_event.setz(EventVertexZA);
        alexPhiMeson_event.setid(RunIdA);
        alexPhiMeson_event.setmult(refMultA);
        alexPhiMeson_event.setn_prim(n_primaries);
        alexPhiMeson_event.setn_non_prim(n_non_primaries);
        alexPhiMeson_event.setn_tof_prim(n_tofmatch_prim);
        alexPhiMeson_event.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        alexPhiMeson_event.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        alexPhiMeson_event.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        alexPhiMeson_event.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        alexPhiMeson_event.setEP_Qx_ptw(EP_Qx_ptw);
        alexPhiMeson_event.setEP_Qy_ptw(EP_Qy_ptw);
        alexPhiMeson_event.setQtracks_eta_pos(Qtracks_used_eta_pos);
        alexPhiMeson_event.setQtracks_eta_neg(Qtracks_used_eta_neg);
        alexPhiMeson_event.setQtracks_full(Qtracks_used);
        alexPhiMeson_event.setZDCx(ZDCx);
        alexPhiMeson_event.setBBCx(BBCx);
        alexPhiMeson_event.setvzVpd(vzVpd);

        Int_t flag_tof_KM, flag_tof_KP;

        dummy_counter++;
        //cout << "************** START of event *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB, helixA_glob, helixB_glob;
        TLorentzVector ltrackA, ltrackB;
        StThreeVectorF primdirA, primdirB, vectorA, vectorB, vectorAB;
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // K- candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );

            Float_t MomentumA   = trackA.pMom().mag();
            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            //Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
            Float_t nSigmaKaonA = trackA.nSigmaKaon();

            // calculate mass2
            Float_t Mass2A      = -100.0;
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
                flag_tof_KM = 1;
            }
            else
            {
                flag_tof_KM = 0;
            }

            if(
               dcaA           < 1.5
               && nHitsFitA   > 14
               && nHitsPossA  > 0
               && (nHitsFitA/nHitsPossA) > 0.52
               && MomentumA   > 0.1
               && MomentumA   < 10.0
               && nSigmaKaonA < 2.5/nsigma_scaling_fac
               && nSigmaKaonA > -2.5/nsigma_scaling_fac
               && ((MomentumA <= 0.65 && flag_tof_KM == 0) ||
                   (flag_tof_KM == 1 && ((MomentumA < 1.5 && Mass2A > 0.16 && Mass2A < 0.36) || (MomentumA >= 1.5 && Mass2A > 0.125 && Mass2A < 0.36)) )
                  )
               //&& (
               //    flag_tof_KM == 0 ||
               //    (flag_tof_KM == 1 && ((MomentumA < 1.5 && Mass2A > 0.16 && Mass2A < 0.36) || (MomentumA >= 1.5 && Mass2A > 0.125 && Mass2A < 0.36)) )
               //   )
              )
            {
                helixA   = StPhysicalHelixD(trackA.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());    // gMom or pMom?
                primdirA = helixA.cat(helixA.pathLength(vector_prim)); //
                primdirA = MomentumA*primdirA/primdirA.mag();
                ltrackA.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),0.493677);

                Float_t ptA         = ltrackA.Pt();
                Float_t pzA         = ltrackA.Pz();
                Float_t etaA        = primdirA.pseudoRapidity();

                if(
                   fabs(etaA)     < 1.0
                   && ptA         > 0.1
                   && ptA         < 10.0
                  )
                {
                    for(Int_t j = 0; j < PID_counter_Array_B[Ana_Num][ParticleB]; j++)  // K+ candidates
                    {
                        //cout << "i = " << i << ", j = " << j << endl;
                        Int_t trackB_num = PID_Array_B[Ana_Num][ParticleB][j];

                        if(
                           (trackA_num != trackB_num)
                           || SE_ME_Flag == 1
                          )
                        {
                            StPicoAlexTrack trackB = *picoDst_B->track( trackB_num );

                            Float_t MomentumB     = trackB.pMom().mag();
                            Float_t dcaB          = trackB.dca();
                            Float_t nHitsPossB    = trackB.nHitsMax();
                            Float_t nHitsFitB     = trackB.nHitsFit();
                            Float_t BetaB         = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                            Float_t PolarityB     = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                            //Float_t TPCdEdxB      = trackB.dEdx(); // Combined inner and outer MDC dE/dx
                            Float_t nSigmaKaonB   = trackB.nSigmaKaon();

                            // calculate mass2
                            Float_t Mass2B        = -100.0;
                            if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                            {
                                Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                                flag_tof_KP = 1;
                            }
                            else
                            {
                                flag_tof_KP = 0;
                            }

                            if(
                               dcaB           < 1.5
                               && nHitsFitB   > 14
                               && nHitsPossB  > 0
                               && (nHitsFitB/nHitsPossB) > 0.52
                               && MomentumB   > 0.1
                               && MomentumB   < 10.0
                               && nSigmaKaonB < 2.5/nsigma_scaling_fac
                               && nSigmaKaonB > -2.5/nsigma_scaling_fac
                               && ((MomentumB <= 0.65 && flag_tof_KP == 0) ||
                                   (flag_tof_KP == 1 && ((MomentumB < 1.5 && Mass2B > 0.16 && Mass2B < 0.36) || (MomentumB >= 1.5 && Mass2B > 0.125 && Mass2B < 0.36)) )
                                  )
                               //&& (
                               //    flag_tof_KP == 0 ||
                               //    (flag_tof_KP == 1 && ((MomentumB < 1.5 && Mass2B > 0.16 && Mass2B < 0.36) || (MomentumB >= 1.5 && Mass2B > 0.125 && Mass2B < 0.36)) )
                               //   )
                              )
                            {
                                helixB = StPhysicalHelixD(trackB.pMom(),event_B_ana->primaryVertex() + vectordiff,event_B_ana->bField()*MAGFIELDFACTOR,trackB.charge());
                                primdirB = helixB.cat(helixB.pathLength(vector_prim)); //
                                primdirB = MomentumB*primdirB/primdirB.mag();
                                ltrackB.SetXYZM(primdirB.x(),primdirB.y(),primdirB.z(),0.493677);

                                Float_t ptB           = ltrackB.Pt();
                                Float_t pzB           = ltrackB.Pz();
                                Float_t etaB          = primdirB.pseudoRapidity();

                                if(
                                   fabs(etaB)     < 1.0
                                   && ptB         > 0.1
                                   && ptB         < 10.0
                                  )
                                {
                                    // Invariant mass calculations
                                    TLorentzVector trackAB      = ltrackA+ltrackB;
                                    Double_t InvMassAB          = trackAB.M();
                                    vectorAB.set(trackAB.X(),trackAB.Y(),trackAB.Z());
                                    Float_t phiAB   = vectorAB.phi();
                                    Float_t thetaAB = vectorAB.theta();

                                    Float_t DeltaDipAngle = TMath::ACos((ptA*ptB+pzA*pzB)/(MomentumA*MomentumB));  // To remove electron conversion background, see Phys.Rev.C79:064903,2009
                                    Float_t pt            = trackAB.Pt();  // Transverse momentum of mother particle
                                    Float_t rap           = trackAB.Rapidity(); // Rapidity of mother particle

                                    //cout << "Event number = " << event_number << ", InvMassAB = " << InvMassAB << ", dcaA = " << dcaA << ", dcaB = " << dcaB
                                    //    << ", BetaA = " << BetaA << ", BetaB = " << BetaB << ", FitA = " << nHitsFitA << ", FitB = " << nHitsFitB
                                    //    << ", ptA = " << ptA << ", ptB = " << ptB << endl;

                                    if(
                                       InvMassAB      < 1.15
                                      )
                                    {
                                        Float_t pathA_prim,dcaA_prim,pathB_prim,dcaB_prim;

                                        helixA_glob = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
                                        helixB_glob = StPhysicalHelixD(trackB.gMom(),trackB.origin() + vectordiff,event_B_ana->bField()*MAGFIELDFACTOR, trackB.charge());
                                        fHelixAtoPointdca(vector_prim,helixA_glob,pathA_prim,dcaA_prim);
                                        fHelixAtoPointdca(vector_prim,helixB_glob,pathB_prim,dcaB_prim);

                                        StThreeVectorF vectornewA_prim,vectornewB_prim;
                                        vectornewA_prim     = helixA_glob.cat(pathA_prim);
                                        vectornewB_prim     = helixB_glob.cat(pathB_prim);

                                        vectornewA_prim     = trackA.gMom().mag()*vectornewA_prim/vectornewA_prim.mag(); // momentum vector at primary vertex
                                        vectornewB_prim     = trackB.gMom().mag()*vectornewB_prim/vectornewB_prim.mag(); // momentum vector at primary vertex

                                        Float_t phi_event_plane                 = -400.0;
                                        Float_t phi_event_plane_eta_gap         = -400.0;
                                        Float_t delta_phi_ME_AB_weight          = 0.0;
                                        Float_t delta_phi_ME_AB_weight_eta_gap  = 0.0;

                                        Int_t delta_phi_ME = 0;

                                        // Calculate the event plane anlges
                                        calc_event_plane_angles(2,rap,trackA,trackB,trackA,vectornewA_prim,vectornewB_prim,vectornewA_prim,ME_Flag,SE_ME_Flag,RunIdA,EventVertexXA,EventVertexYA,
                                                                EventVertexZA,RunIdB,EventVertexXB,EventVertexYB,EventVertexZB,
                                                                phi_event_plane,phi_event_plane_eta_gap,delta_phi_ME_AB_weight,delta_phi_ME_AB_weight_eta_gap);


                                        // check whether the event planes between event A and B are close to each other
                                        if(
                                           fabs(delta_phi_ME_AB_weight) < TMath::DegToRad()*30.0
                                           || (SE_ME_Flag == 0)
                                          )
                                        {
                                            delta_phi_ME = 1;
                                        }

                                        if(delta_phi_ME == 1)
                                        {
                                            Float_t p_xA_c   = vectornewA_prim.x();
                                            Float_t p_yA_c   = vectornewA_prim.y();
                                            Float_t p_tA_c   = sqrt(p_xA_c*p_xA_c + p_yA_c*p_yA_c);
                                            Float_t etaA_c   = vectornewA_prim.pseudoRapidity();
                                            Float_t phiA_c   = vectornewA_prim.phi();
                                            Double_t p_t_weightA = 1.0;
                                            if(p_tA_c < 2.0)  p_t_weightA = p_tA_c;
                                            if(p_tA_c >= 2.0) p_t_weightA = 2.0;
                                            Float_t iQxA     = p_t_weightA*TMath::Cos(2.0*phiA_c);
                                            Float_t iQyA     = p_t_weightA*TMath::Sin(2.0*phiA_c);

                                            Float_t p_xB_c   = vectornewB_prim.x();
                                            Float_t p_yB_c   = vectornewB_prim.y();
                                            Float_t p_tB_c   = sqrt(p_xB_c*p_xB_c + p_yB_c*p_yB_c);
                                            Float_t etaB_c   = vectornewB_prim.pseudoRapidity();
                                            Float_t phiB_c   = vectornewB_prim.phi();
                                            Double_t p_t_weightB = 1.0;
                                            if(p_tB_c < 2.0)  p_t_weightB = p_tB_c;
                                            if(p_tB_c >= 2.0) p_t_weightB = 2.0;
                                            Float_t iQxB     = p_t_weightB*TMath::Cos(2.0*phiB_c);
                                            Float_t iQyB     = p_t_weightB*TMath::Sin(2.0*phiB_c);


                                            // Filling Ntuple
                                            alexPhiMeson_track = alexPhiMeson_event.createTrack();
                                            alexPhiMeson_track->setm2A(Mass2A);
                                            alexPhiMeson_track->setm2B(Mass2B);
                                            alexPhiMeson_track->setnsA(nSigmaKaonA);
                                            alexPhiMeson_track->setnsB(nSigmaKaonB);
                                            alexPhiMeson_track->setdcaA(dcaA);
                                            alexPhiMeson_track->setdcaB(dcaB);
                                            alexPhiMeson_track->setiQxA(iQxA);
                                            alexPhiMeson_track->setiQyA(iQyA);
                                            alexPhiMeson_track->setiQxB(iQxB);
                                            alexPhiMeson_track->setiQyB(iQyB);
                                            alexPhiMeson_track->setetaA(etaA_c);
                                            alexPhiMeson_track->setetaB(etaB_c);
                                            alexPhiMeson_track->setInvAB(InvMassAB);
                                            alexPhiMeson_track->setpt(pt);
                                            alexPhiMeson_track->setrap(rap);
                                            alexPhiMeson_track->setphi(phiAB);
                                            alexPhiMeson_track->settheta(thetaAB);
                                            alexPhiMeson_track->setPsi_ep(phi_event_plane);
                                            alexPhiMeson_track->setPsi_ep_eta(phi_event_plane_eta_gap);
                                            alexPhiMeson_track->setPsi_diff_ME(delta_phi_ME_AB_weight);
                                            alexPhiMeson_track->setDelta_dip_angle(DeltaDipAngle);
                                            alexPhiMeson_track->setqpA(MomentumA*PolarityA);
                                            alexPhiMeson_track->setqpB(MomentumB*PolarityB);

                                            //cout << "InvMassAB = " << InvMassAB << endl;
                                            dummy_counter_A++;
                                        }
                                        //*********************************************************************************************
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(fAnalysisNum == 1) Tree_PhiMeson_v2  ->Fill();
        return 1;
    }
    else return 0;
}



Int_t LambdaCplus_analysis(Int_t PID_counter_Array[][N_max_PIDs],Int_t PID_counter_Array_B[][N_max_PIDs],
                       Int_t PID_Array[][N_max_PIDs][N_max_tracks],Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                       StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                       Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t Ana_Num,
                       Int_t run_events,Int_t run_events_B,Int_t event_number, Int_t SE_ME_Flag)
{
    // Au + Au -> LambdaC+ + N + N
    //             |
    //             -> p + pi+ + K-

    // Event vertex information
    StThreeVectorF vector_prim,vector_primB,vectordiff;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB;

    event_A_ana   = picoDst_A;
    event_B_ana   = picoDst_B;

    EventVertexXA    = event_A_ana   ->primaryVertex().x();
    EventVertexYA    = event_A_ana   ->primaryVertex().y();
    EventVertexZA    = event_A_ana   ->primaryVertex().z();
    EventVertexXB    = EventVertexXA;
    EventVertexYB    = EventVertexYA;
    EventVertexZB    = EventVertexZA;
    refMultA         = event_A_ana   ->refMult();
    refMultB         = refMultA;
    Int_t   RunIdA   = event_A_ana   ->runId();
    Int_t   RunIdB   = event_B_ana   ->runId();
    Float_t ZDCx     = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx     = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd    = event_A_ana->vzVpd();
    vector_prim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    //cout << "Entered phi analysis, vertex = {" << EventVertexXA << ", " << EventVertexYA << ", " << EventVertexZA << "}" << endl;

    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        //cout << "Start mixing, check if events can be mixed..." << endl;
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }

    Float_t radius_cut           = Event_radius_cut*Event_radius_cut; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
    }


    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]   > 0   // p
       && PID_counter_Array[Ana_Num][ParticleB]   > 0   // pi+
       && PID_counter_Array_B[Ana_Num][ParticleC] > 0   // K-
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana   ->isMinBias()
       && event_B_ana   ->isMinBias()
      )
    {
        Int_t flag_tof_KM, flag_tof_p, flag_tof_piM;

        dummy_counter++;
        //cout << "************** START of event *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB, helixC;
        TLorentzVector ltrackA, ltrackB, ltrackC;
        StThreeVectorF primdirA, primdirB, primdirC, vectorA, vectorB, vectorC, vectorABC;
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // p candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );

            Float_t MomentumA    = trackA.pMom().mag();
            Float_t dcaA         = trackA.dca();
            Float_t nHitsPossA   = trackA.nHitsMax();
            Float_t nHitsFitA    = trackA.nHitsFit();
            Float_t BetaA        = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t PolarityA    = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t nSigmaProton = trackA.nSigmaProton();

            // calculate mass2
            Float_t Mass2A      = -100.0;
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
                flag_tof_p = 1;
            }
            else
            {
                flag_tof_p = 0;
            }

            if(
               dcaA           < 1.5
               && nHitsFitA   > 14
               && nHitsPossA  > 0
               && (nHitsFitA/nHitsPossA) > 0.52
               && MomentumA   > 0.1
               && MomentumA   < 10.0
               && nSigmaProton < 2.5/nsigma_scaling_fac
               && nSigmaProton > -2.5/nsigma_scaling_fac
               && ((MomentumA <= 1.0 && flag_tof_p == 0) ||
                   (flag_tof_p == 1 && Mass2A > 0.45 && Mass2A < 1.5 )
                  )
              )
            {
                helixA   = StPhysicalHelixD(trackA.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());    // gMom or pMom?
                primdirA = helixA.cat(helixA.pathLength(vector_prim)); //
                primdirA = MomentumA*primdirA/primdirA.mag();
                ltrackA.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),mass_array[ParticleA]);

                Float_t ptA         = ltrackA.Pt();
                Float_t pzA         = ltrackA.Pz();
                Float_t etaA        = primdirA.pseudoRapidity();

                if(
                   fabs(etaA)     < 1.0
                   && ptA         > 0.1
                   && ptA         < 10.0
                  )
                {
                    for(Int_t j = 0; j < PID_counter_Array[Ana_Num][ParticleB]; j++)  // pi+ candidates
                    {
                        //cout << "i = " << i << ", j = " << j << endl;
                        Int_t trackB_num = PID_Array[Ana_Num][ParticleB][j];

                        if(
                           (trackA_num != trackB_num)
                           //|| SE_ME_Flag == 1
                          )
                        {
                            StPicoAlexTrack trackB = *picoDst_A->track( trackB_num );

                            Float_t MomentumB     = trackB.pMom().mag();
                            Float_t dcaB          = trackB.dca();
                            Float_t nHitsPossB    = trackB.nHitsMax();
                            Float_t nHitsFitB     = trackB.nHitsFit();
                            Float_t BetaB         = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                            Float_t PolarityB     = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                            Float_t nSigmaPion    = trackB.nSigmaPion();

                            // calculate mass2
                            Float_t Mass2B        = -100.0;
                            if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                            {
                                Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                                flag_tof_piM = 1;
                            }
                            else
                            {
                                flag_tof_piM = 0;
                            }

                            if(
                               dcaB           < 1.5
                               && nHitsFitB   > 14
                               && nHitsPossB  > 0
                               && (nHitsFitB/nHitsPossB) > 0.52
                               && MomentumB   > 0.1
                               && MomentumB   < 10.0
                               && nSigmaPion < 2.5/nsigma_scaling_fac
                               && nSigmaPion > -2.5/nsigma_scaling_fac
                               && ((MomentumB <= 0.65 && flag_tof_piM == 0) ||
                                   (flag_tof_piM == 1 && Mass2B > -0.1 && Mass2B < 0.1)
                                  )
                              )
                            {
                                helixB = StPhysicalHelixD(trackB.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackB.charge());
                                primdirB = helixB.cat(helixB.pathLength(vector_prim)); //
                                primdirB = MomentumB*primdirB/primdirB.mag();
                                ltrackB.SetXYZM(primdirB.x(),primdirB.y(),primdirB.z(),mass_array[ParticleB]);

                                Float_t ptB           = ltrackB.Pt();
                                Float_t pzB           = ltrackB.Pz();
                                Float_t etaB          = primdirB.pseudoRapidity();

                                if(
                                   fabs(etaB)     < 1.0
                                   && ptB         > 0.1
                                   && ptB         < 10.0
                                  )
                                {
                                    for(Int_t k = 0; k < PID_counter_Array_B[Ana_Num][ParticleC]; k++)  // K- candidates
                                    {
                                        Int_t trackC_num = PID_Array_B[Ana_Num][ParticleC][k];

                                        if(
                                           ((trackC_num != trackA_num) && (trackC_num != trackB_num))
                                           || SE_ME_Flag == 1
                                          )
                                        {
                                            StPicoAlexTrack trackC = *picoDst_B->track( trackC_num );

                                            Float_t MomentumC     = trackC.pMom().mag();
                                            Float_t dcaC          = trackC.dca();
                                            Float_t nHitsPossC    = trackC.nHitsMax();
                                            Float_t nHitsFitC     = trackC.nHitsFit();
                                            Float_t BetaC         = trackC.btofBeta();  // Velocity after time-of-flight reconstruction
                                            Float_t PolarityC     = trackC.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                            Float_t nSigmaKaon    = trackC.nSigmaKaon();

                                            // calculate mass2
                                            Float_t Mass2C        = -100.0;
                                            if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                                            {
                                                Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                                                flag_tof_KM = 1;
                                            }
                                            else
                                            {
                                                flag_tof_KM = 0;
                                            }

                                            if(
                                               dcaC           < 1.5
                                               && nHitsFitC   > 14
                                               && nHitsPossC  > 0
                                               && (nHitsFitC/nHitsPossC) > 0.52
                                               && MomentumC   > 0.1
                                               && MomentumC   < 10.0
                                               && nSigmaKaon < 2.5/nsigma_scaling_fac
                                               && nSigmaKaon > -2.5/nsigma_scaling_fac
                                               && ((MomentumC <= 0.65 && flag_tof_KM == 0) ||
                                                   (flag_tof_KM == 1 && Mass2C > 0.125 && Mass2C < 0.36)
                                                  )
                                              )
                                            {
                                                helixC = StPhysicalHelixD(trackC.pMom(),event_B_ana->primaryVertex() + vectordiff,event_B_ana->bField()*MAGFIELDFACTOR,trackC.charge());
                                                primdirC = helixC.cat(helixC.pathLength(vector_prim)); //
                                                primdirC = MomentumC*primdirC/primdirC.mag();
                                                ltrackC.SetXYZM(primdirC.x(),primdirC.y(),primdirC.z(),mass_array[ParticleC]);

                                                Float_t ptC           = ltrackC.Pt();
                                                Float_t pzC           = ltrackC.Pz();
                                                Float_t etaC          = primdirC.pseudoRapidity();

                                                if(
                                                   fabs(etaC)     < 1.0
                                                   && ptC         > 0.1
                                                   && ptC         < 10.0
                                                  )
                                                {
                                                    // Invariant mass calculations
                                                    TLorentzVector trackABC      = ltrackA+ltrackB+ltrackC;
                                                    Double_t InvMassABC          = trackABC.M();
                                                    vectorABC.set(trackABC.X(),trackABC.Y(),trackABC.Z());
                                                    Float_t phiABC   = vectorABC.phi();
                                                    Float_t thetaABC = vectorABC.theta();

                                                    Float_t pt            = trackABC.Pt();  // Transverse momentum of mother particle
                                                    Float_t rap           = trackABC.Rapidity(); // Rapidity of mother particle

                                                    if(
                                                       InvMassABC      < 5.0
                                                      )
                                                    {
                                                        h_2D_InvMass_LambdaCplus[0] ->Fill(pt,InvMassABC);

                                                        if(flag_tof_KM == 1 && flag_tof_piM == 1 && flag_tof_p == 1)
                                                        {
                                                            h_2D_InvMass_LambdaCplus[1] ->Fill(pt,InvMassABC);
                                                            if(dcaA < 1.0 && dcaB < 1.0  && dcaC < 1.0)
                                                            {
                                                                h_2D_InvMass_LambdaCplus[2] ->Fill(pt,InvMassABC);

                                                                if(erefMult_bin < 9 && erefMult_bin >= 0)
                                                                {
                                                                    h_2D_InvMass_LambdaCplus[3+erefMult_bin] ->Fill(pt,InvMassABC);
                                                                }
                                                            }
                                                        }
                                                        dummy_counter_A++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return 1;
    }
    else return 0;
}



Int_t JPsi_analysis(Int_t PID_counter_Array[][N_max_PIDs],Int_t PID_counter_Array_B[][N_max_PIDs],
                       Int_t PID_Array[][N_max_PIDs][N_max_tracks],Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                       StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                       Int_t ParticleA, Int_t ParticleB, Int_t Ana_Num,
                       Int_t run_events,Int_t run_events_B,Int_t event_number, Int_t SE_ME_Flag)
{
    // Au + Au -> J/Psi + N + N
    //             |
    //             -> e- + e+

    // Event vertex information
    StThreeVectorF vector_prim,vector_primB,vectordiff;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB;

    event_A_ana   = picoDst_A;
    event_B_ana   = picoDst_B;

    EventVertexXA    = event_A_ana   ->primaryVertex().x();
    EventVertexYA    = event_A_ana   ->primaryVertex().y();
    EventVertexZA    = event_A_ana   ->primaryVertex().z();
    EventVertexXB    = EventVertexXA;
    EventVertexYB    = EventVertexYA;
    EventVertexZB    = EventVertexZA;
    refMultA         = event_A_ana   ->refMult();
    refMultB         = refMultA;
    Int_t   RunIdA   = event_A_ana   ->runId();
    Int_t   RunIdB   = event_B_ana   ->runId();
    Float_t ZDCx     = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx     = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd    = event_A_ana->vzVpd();
    vector_prim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    //cout << "Entered phi analysis, vertex = {" << EventVertexXA << ", " << EventVertexYA << ", " << EventVertexZA << "}" << endl;

    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        //cout << "Start mixing, check if events can be mixed..." << endl;
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }

    Float_t radius_cut           = Event_radius_cut*Event_radius_cut; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
    }


    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]   > 0   // e-
       && PID_counter_Array_B[Ana_Num][ParticleB] > 0   // e+
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana   ->isMinBias()
       && event_B_ana   ->isMinBias()
      )
    {
        alexPhiMeson_event.clearTrackList();
        alexPhiMeson_event.setx(EventVertexXA);
        alexPhiMeson_event.sety(EventVertexYA);
        alexPhiMeson_event.setz(EventVertexZA);
        alexPhiMeson_event.setid(RunIdA);
        alexPhiMeson_event.setmult(refMultA);
        alexPhiMeson_event.setn_prim(n_primaries);
        alexPhiMeson_event.setn_non_prim(n_non_primaries);
        alexPhiMeson_event.setn_tof_prim(n_tofmatch_prim);
        alexPhiMeson_event.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        alexPhiMeson_event.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        alexPhiMeson_event.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        alexPhiMeson_event.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        alexPhiMeson_event.setEP_Qx_ptw(EP_Qx_ptw);
        alexPhiMeson_event.setEP_Qy_ptw(EP_Qy_ptw);
        alexPhiMeson_event.setQtracks_eta_pos(Qtracks_used_eta_pos);
        alexPhiMeson_event.setQtracks_eta_neg(Qtracks_used_eta_neg);
        alexPhiMeson_event.setQtracks_full(Qtracks_used);
        alexPhiMeson_event.setZDCx(ZDCx);
        alexPhiMeson_event.setBBCx(BBCx);
        alexPhiMeson_event.setvzVpd(vzVpd);

        Int_t flag_tof_KM, flag_tof_KP;

        dummy_counter++;
        //cout << "************** START of event *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB, helixA_glob, helixB_glob;
        TLorentzVector ltrackA, ltrackB;
        StThreeVectorF primdirA, primdirB, vectorA, vectorB, vectorAB;
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // K- candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );

            Float_t MomentumA   = trackA.pMom().mag();
            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            //Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
            Float_t nSigmaKaonA = trackA.nSigmaElectron();

            // calculate mass2
            Float_t Mass2A      = -100.0;
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
                flag_tof_KM = 1;
            }
            else
            {
                flag_tof_KM = 0;
            }

            if(
               dcaA           < 1.5
               && nHitsFitA   > 14
               && nHitsPossA  > 0
               && (nHitsFitA/nHitsPossA) > 0.52
               && MomentumA   > 0.1
               && MomentumA   < 10.0
               && nSigmaKaonA < 2.5/nsigma_scaling_fac
               && nSigmaKaonA > -2.5/nsigma_scaling_fac
               && Mass2A > -0.1 && Mass2A < 0.01
               //&& ((MomentumA <= 0.65 && flag_tof_KM == 0) ||
               //    (flag_tof_KM == 1 && ((MomentumA < 1.5 && Mass2A > -0.1 && Mass2A < 0.1) || (MomentumA >= 1.5 && Mass2A > -0.1 && Mass2A < 0.1)) )
               //   )
               //&& (
               //    flag_tof_KM == 0 ||
               //    (flag_tof_KM == 1 && ((MomentumA < 1.5 && Mass2A > 0.16 && Mass2A < 0.36) || (MomentumA >= 1.5 && Mass2A > 0.125 && Mass2A < 0.36)) )
               //   )
              )
            {
                helixA   = StPhysicalHelixD(trackA.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());    // gMom or pMom?
                primdirA = helixA.cat(helixA.pathLength(vector_prim)); //
                primdirA = MomentumA*primdirA/primdirA.mag();
                ltrackA.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),0.00051099892);

                Float_t ptA         = ltrackA.Pt();
                Float_t pzA         = ltrackA.Pz();
                Float_t etaA        = primdirA.pseudoRapidity();

                if(
                   fabs(etaA)     < 1.0
                   && ptA         > 0.1
                   && ptA         < 10.0
                  )
                {
                    for(Int_t j = 0; j < PID_counter_Array_B[Ana_Num][ParticleB]; j++)  // K+ candidates
                    {
                        //cout << "i = " << i << ", j = " << j << endl;
                        Int_t trackB_num = PID_Array_B[Ana_Num][ParticleB][j];

                        if(
                           (trackA_num != trackB_num)
                           || SE_ME_Flag == 1
                          )
                        {
                            StPicoAlexTrack trackB = *picoDst_B->track( trackB_num );

                            Float_t MomentumB     = trackB.pMom().mag();
                            Float_t dcaB          = trackB.dca();
                            Float_t nHitsPossB    = trackB.nHitsMax();
                            Float_t nHitsFitB     = trackB.nHitsFit();
                            Float_t BetaB         = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                            Float_t PolarityB     = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                            //Float_t TPCdEdxB      = trackB.dEdx(); // Combined inner and outer MDC dE/dx
                            Float_t nSigmaKaonB   = trackB.nSigmaElectron();

                            // calculate mass2
                            Float_t Mass2B        = -100.0;
                            if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                            {
                                Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                                flag_tof_KP = 1;
                            }
                            else
                            {
                                flag_tof_KP = 0;
                            }

                            if(
                               dcaB           < 1.5
                               && nHitsFitB   > 14
                               && nHitsPossB  > 0
                               && (nHitsFitB/nHitsPossB) > 0.52
                               && MomentumB   > 0.1
                               && MomentumB   < 10.0
                               && nSigmaKaonB < 2.5/nsigma_scaling_fac
                               && nSigmaKaonB > -2.5/nsigma_scaling_fac
                               && Mass2B > -0.1 && Mass2B < 0.01
                               //&& ((MomentumB <= 0.65 && flag_tof_KP == 0) ||
                               //    (flag_tof_KP == 1 && ((MomentumB < 1.5 && Mass2B > -0.1 && Mass2B < 0.1) || (MomentumB >= 1.5 && Mass2B > -0.1 && Mass2B < 0.1)) )
                               //   )
                               //&& (
                               //    flag_tof_KP == 0 ||
                               //    (flag_tof_KP == 1 && ((MomentumB < 1.5 && Mass2B > 0.16 && Mass2B < 0.36) || (MomentumB >= 1.5 && Mass2B > 0.125 && Mass2B < 0.36)) )
                               //   )
                              )
                            {
                                helixB = StPhysicalHelixD(trackB.pMom(),event_B_ana->primaryVertex() + vectordiff,event_B_ana->bField()*MAGFIELDFACTOR,trackB.charge());
                                primdirB = helixB.cat(helixB.pathLength(vector_prim)); //
                                primdirB = MomentumB*primdirB/primdirB.mag();
                                ltrackB.SetXYZM(primdirB.x(),primdirB.y(),primdirB.z(),0.00051099892);

                                Float_t ptB           = ltrackB.Pt();
                                Float_t pzB           = ltrackB.Pz();
                                Float_t etaB          = primdirB.pseudoRapidity();

                                if(
                                   fabs(etaB)     < 1.0
                                   && ptB         > 0.1
                                   && ptB         < 10.0
                                  )
                                {
                                    // Invariant mass calculations
                                    TLorentzVector trackAB      = ltrackA+ltrackB;
                                    Double_t InvMassAB          = trackAB.M();
                                    vectorAB.set(trackAB.X(),trackAB.Y(),trackAB.Z());
                                    Float_t phiAB   = vectorAB.phi();
                                    Float_t thetaAB = vectorAB.theta();

                                    Float_t DeltaDipAngle = TMath::ACos((ptA*ptB+pzA*pzB)/(MomentumA*MomentumB));  // To remove electron conversion background, see Phys.Rev.C79:064903,2009
                                    Float_t pt            = trackAB.Pt();  // Transverse momentum of mother particle
                                    Float_t rap           = trackAB.Rapidity(); // Rapidity of mother particle

                                    //cout << "Event number = " << event_number << ", InvMassAB = " << InvMassAB << ", dcaA = " << dcaA << ", dcaB = " << dcaB
                                    //    << ", BetaA = " << BetaA << ", BetaB = " << BetaB << ", FitA = " << nHitsFitA << ", FitB = " << nHitsFitB
                                    //    << ", ptA = " << ptA << ", ptB = " << ptB << endl;

                                    if(
                                       InvMassAB      < 10.0
                                      )
                                    {
                                        Float_t pathA_prim,dcaA_prim,pathB_prim,dcaB_prim;

                                        helixA_glob = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
                                        helixB_glob = StPhysicalHelixD(trackB.gMom(),trackB.origin() + vectordiff,event_B_ana->bField()*MAGFIELDFACTOR, trackB.charge());
                                        fHelixAtoPointdca(vector_prim,helixA_glob,pathA_prim,dcaA_prim);
                                        fHelixAtoPointdca(vector_prim,helixB_glob,pathB_prim,dcaB_prim);

                                        StThreeVectorF vectornewA_prim,vectornewB_prim;
                                        vectornewA_prim     = helixA_glob.cat(pathA_prim);
                                        vectornewB_prim     = helixB_glob.cat(pathB_prim);

                                        vectornewA_prim     = trackA.gMom().mag()*vectornewA_prim/vectornewA_prim.mag(); // momentum vector at primary vertex
                                        vectornewB_prim     = trackB.gMom().mag()*vectornewB_prim/vectornewB_prim.mag(); // momentum vector at primary vertex

                                        Float_t phi_event_plane                 = -400.0;
                                        Float_t phi_event_plane_eta_gap         = -400.0;
                                        Float_t delta_phi_ME_AB_weight          = 0.0;
                                        Float_t delta_phi_ME_AB_weight_eta_gap  = 0.0;

                                        Int_t delta_phi_ME = 0;

                                        // Calculate the event plane anlges
                                        calc_event_plane_angles(2,rap,trackA,trackB,trackA,vectornewA_prim,vectornewB_prim,vectornewA_prim,ME_Flag,SE_ME_Flag,RunIdA,EventVertexXA,EventVertexYA,
                                                                EventVertexZA,RunIdB,EventVertexXB,EventVertexYB,EventVertexZB,
                                                                phi_event_plane,phi_event_plane_eta_gap,delta_phi_ME_AB_weight,delta_phi_ME_AB_weight_eta_gap);


                                        // check whether the event planes between event A and B are close to each other
                                        if(
                                           fabs(delta_phi_ME_AB_weight) < TMath::DegToRad()*30.0
                                           || (SE_ME_Flag == 0)
                                          )
                                        {
                                            delta_phi_ME = 1;
                                        }

                                        if(delta_phi_ME == 1)
                                        {
                                            Float_t p_xA_c   = vectornewA_prim.x();
                                            Float_t p_yA_c   = vectornewA_prim.y();
                                            Float_t p_tA_c   = sqrt(p_xA_c*p_xA_c + p_yA_c*p_yA_c);
                                            Float_t etaA_c   = vectornewA_prim.pseudoRapidity();
                                            Float_t phiA_c   = vectornewA_prim.phi();
                                            Double_t p_t_weightA = 1.0;
                                            if(p_tA_c < 2.0)  p_t_weightA = p_tA_c;
                                            if(p_tA_c >= 2.0) p_t_weightA = 2.0;
                                            Float_t iQxA     = p_t_weightA*TMath::Cos(2.0*phiA_c);
                                            Float_t iQyA     = p_t_weightA*TMath::Sin(2.0*phiA_c);

                                            Float_t p_xB_c   = vectornewB_prim.x();
                                            Float_t p_yB_c   = vectornewB_prim.y();
                                            Float_t p_tB_c   = sqrt(p_xB_c*p_xB_c + p_yB_c*p_yB_c);
                                            Float_t etaB_c   = vectornewB_prim.pseudoRapidity();
                                            Float_t phiB_c   = vectornewB_prim.phi();
                                            Double_t p_t_weightB = 1.0;
                                            if(p_tB_c < 2.0)  p_t_weightB = p_tB_c;
                                            if(p_tB_c >= 2.0) p_t_weightB = 2.0;
                                            Float_t iQxB     = p_t_weightB*TMath::Cos(2.0*phiB_c);
                                            Float_t iQyB     = p_t_weightB*TMath::Sin(2.0*phiB_c);


                                            // Filling Ntuple
                                            alexPhiMeson_track = alexPhiMeson_event.createTrack();
                                            alexPhiMeson_track->setm2A(Mass2A);
                                            alexPhiMeson_track->setm2B(Mass2B);
                                            alexPhiMeson_track->setnsA(nSigmaKaonA);
                                            alexPhiMeson_track->setnsB(nSigmaKaonB);
                                            alexPhiMeson_track->setdcaA(dcaA);
                                            alexPhiMeson_track->setdcaB(dcaB);
                                            alexPhiMeson_track->setiQxA(iQxA);
                                            alexPhiMeson_track->setiQyA(iQyA);
                                            alexPhiMeson_track->setiQxB(iQxB);
                                            alexPhiMeson_track->setiQyB(iQyB);
                                            alexPhiMeson_track->setetaA(etaA_c);
                                            alexPhiMeson_track->setetaB(etaB_c);
                                            alexPhiMeson_track->setInvAB(InvMassAB);
                                            alexPhiMeson_track->setpt(pt);
                                            alexPhiMeson_track->setrap(rap);
                                            alexPhiMeson_track->setphi(phiAB);
                                            alexPhiMeson_track->settheta(thetaAB);
                                            alexPhiMeson_track->setPsi_ep(phi_event_plane);
                                            alexPhiMeson_track->setPsi_ep_eta(phi_event_plane_eta_gap);
                                            alexPhiMeson_track->setPsi_diff_ME(delta_phi_ME_AB_weight);
                                            alexPhiMeson_track->setDelta_dip_angle(DeltaDipAngle);
                                            alexPhiMeson_track->setqpA(MomentumA*PolarityA);
                                            alexPhiMeson_track->setqpB(MomentumB*PolarityB);

                                            //cout << "InvMassAB = " << InvMassAB << endl;
                                            dummy_counter_A++;
                                        }
                                        //*********************************************************************************************
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(fAnalysisNum == 25) Tree_PhiMeson_v2  ->Fill();
        return 1;
    }
    else return 0;
}




Int_t rho_analysis(Int_t PID_counter_Array[][N_max_PIDs],Int_t PID_counter_Array_B[][N_max_PIDs],
                       Int_t PID_Array[][N_max_PIDs][N_max_tracks],Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                       StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                       Int_t ParticleA, Int_t ParticleB, Int_t Ana_Num,
                       Int_t run_events,Int_t run_events_B,Int_t event_number, Int_t SE_ME_Flag)
{
    // Au + Au ->  rho + N + N
    //              |
    //              -> pi- + pi+

    // Event vertex information
    StThreeVectorF vector_prim,vector_primB,vectordiff;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB;

    event_A_ana   = picoDst_A;
    event_B_ana   = picoDst_B;

    EventVertexXA    = event_A_ana   ->primaryVertex().x();
    EventVertexYA    = event_A_ana   ->primaryVertex().y();
    EventVertexZA    = event_A_ana   ->primaryVertex().z();
    EventVertexXB    = EventVertexXA;
    EventVertexYB    = EventVertexYA;
    EventVertexZB    = EventVertexZA;
    refMultA         = event_A_ana   ->refMult();
    refMultB         = refMultA;
    Int_t   RunIdA   = event_A_ana   ->runId();
    Int_t   RunIdB   = event_B_ana   ->runId();
    //Float_t ZDCx     = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    //Float_t BBCx     = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    //Float_t vzVpd    = event_A_ana->vzVpd();
    vector_prim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    //cout << "Entered phi analysis, vertex = {" << EventVertexXA << ", " << EventVertexYA << ", " << EventVertexZA << "}" << endl;

    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        //cout << "Start mixing, check if events can be mixed..." << endl;
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }

    Float_t radius_cut           = Event_radius_cut*Event_radius_cut; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
    }


    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]   > 0   // pi+
       && PID_counter_Array_B[Ana_Num][ParticleB] > 0   // pi-
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana   ->isMinBias()
       && event_B_ana   ->isMinBias()
      )
    {
        Int_t flag_tof_KM, flag_tof_KP;

        dummy_counter++;
        //cout << "************** START of event *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB, helixA_glob, helixB_glob;
        TLorentzVector ltrackA, ltrackB;
        StThreeVectorF primdirA, primdirB, vectorA, vectorB, vectorAB;
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // K- candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );

            Float_t MomentumA   = trackA.pMom().mag();
            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            //Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
            Float_t nSigmaPionA = trackA.nSigmaPion();

            // calculate mass2
            Float_t Mass2A      = -100.0;
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
                flag_tof_KM = 1;
            }
            else
            {
                flag_tof_KM = 0;
            }

            if(
               dcaA           < 1.0
               && nHitsFitA   > 14
               && nHitsPossA  > 0
               && (nHitsFitA/nHitsPossA) > 0.52
               && MomentumA   > 0.1
               && MomentumA   < 10.0
               && nSigmaPionA < 2.5/nsigma_scaling_fac
               && nSigmaPionA > -2.5/nsigma_scaling_fac
               && ((MomentumA <= 0.65 && flag_tof_KM == 0) ||
                   (flag_tof_KM == 1 && ((MomentumA < 1.5 && Mass2A > -0.05 && Mass2A < 0.05) || (MomentumA >= 1.5 && Mass2A > -0.05 && Mass2A < 0.05)) )
                  )
               //&& (
               //    flag_tof_KM == 0 ||
               //    (flag_tof_KM == 1 && ((MomentumA < 1.5 && Mass2A > 0.16 && Mass2A < 0.36) || (MomentumA >= 1.5 && Mass2A > 0.125 && Mass2A < 0.36)) )
               //   )
              )
            {
                helixA   = StPhysicalHelixD(trackA.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());    // gMom or pMom?
                primdirA = helixA.cat(helixA.pathLength(vector_prim)); //
                primdirA = MomentumA*primdirA/primdirA.mag();
                ltrackA.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),0.13957018);

                Float_t ptA         = ltrackA.Pt();
                //Float_t pzA         = ltrackA.Pz();
                Float_t etaA        = primdirA.pseudoRapidity();

                if(
                   fabs(etaA)     < 1.0
                   && ptA         > 0.1
                   && ptA         < 10.0
                  )
                {
                    for(Int_t j = 0; j < PID_counter_Array_B[Ana_Num][ParticleB]; j++)  // K+ candidates
                    {
                        //cout << "i = " << i << ", j = " << j << endl;
                        Int_t trackB_num = PID_Array_B[Ana_Num][ParticleB][j];

                        if(
                           (trackA_num != trackB_num)
                           || SE_ME_Flag == 1
                          )
                        {
                            StPicoAlexTrack trackB = *picoDst_B->track( trackB_num );

                            Float_t MomentumB     = trackB.pMom().mag();
                            Float_t dcaB          = trackB.dca();
                            Float_t nHitsPossB    = trackB.nHitsMax();
                            Float_t nHitsFitB     = trackB.nHitsFit();
                            Float_t BetaB         = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                            //Float_t PolarityB     = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                            //Float_t TPCdEdxB      = trackB.dEdx(); // Combined inner and outer MDC dE/dx
                            Float_t nSigmaPionB   = trackB.nSigmaPion();

                            // calculate mass2
                            Float_t Mass2B        = -100.0;
                            if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                            {
                                Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                                flag_tof_KP = 1;
                            }
                            else
                            {
                                flag_tof_KP = 0;
                            }

                            if(
                               dcaB           < 1.0
                               && nHitsFitB   > 14
                               && nHitsPossB  > 0
                               && (nHitsFitB/nHitsPossB) > 0.52
                               && MomentumB   > 0.1
                               && MomentumB   < 10.0
                               && nSigmaPionB < 2.5/nsigma_scaling_fac
                               && nSigmaPionB > -2.5/nsigma_scaling_fac
                               && ((MomentumB <= 0.65 && flag_tof_KP == 0) ||
                                   (flag_tof_KP == 1 && ((MomentumB < 1.5 && Mass2B > -0.05 && Mass2B < 0.05) || (MomentumB >= 1.5 && Mass2B > -0.05 && Mass2B < 0.05)) )
                                  )
                               //&& (
                               //    flag_tof_KP == 0 ||
                               //    (flag_tof_KP == 1 && ((MomentumB < 1.5 && Mass2B > 0.16 && Mass2B < 0.36) || (MomentumB >= 1.5 && Mass2B > 0.125 && Mass2B < 0.36)) )
                               //   )
                              )
                            {
                                helixB = StPhysicalHelixD(trackB.pMom(),event_B_ana->primaryVertex() + vectordiff,event_B_ana->bField()*MAGFIELDFACTOR,trackB.charge());
                                primdirB = helixB.cat(helixB.pathLength(vector_prim)); //
                                primdirB = MomentumB*primdirB/primdirB.mag();
                                ltrackB.SetXYZM(primdirB.x(),primdirB.y(),primdirB.z(),0.13957018);

                                Float_t ptB           = ltrackB.Pt();
                                //Float_t pzB           = ltrackB.Pz();
                                Float_t etaB          = primdirB.pseudoRapidity();

                                if(
                                   fabs(etaB)     < 1.0
                                   && ptB         > 0.1
                                   && ptB         < 10.0
                                  )
                                {
                                    // Invariant mass calculations
                                    TLorentzVector trackAB      = ltrackA+ltrackB;
                                    Double_t InvMassAB          = trackAB.M();
                                    vectorAB.set(trackAB.X(),trackAB.Y(),trackAB.Z());
                                    //Float_t phiAB   = vectorAB.phi();
                                    //Float_t thetaAB = vectorAB.theta();

                                    //Float_t DeltaDipAngle = TMath::ACos((ptA*ptB+pzA*pzB)/(MomentumA*MomentumB));  // To remove electron conversion background, see Phys.Rev.C79:064903,2009
                                    Float_t pt            = trackAB.Pt();  // Transverse momentum of mother particle
                                    Float_t rap           = trackAB.Rapidity(); // Rapidity of mother particle

                                    //cout << "Event number = " << event_number << ", InvMassAB = " << InvMassAB << ", dcaA = " << dcaA << ", dcaB = " << dcaB
                                    //    << ", BetaA = " << BetaA << ", BetaB = " << BetaB << ", FitA = " << nHitsFitA << ", FitB = " << nHitsFitB
                                    //    << ", ptA = " << ptA << ", ptB = " << ptB << endl;

                                    if(
                                       InvMassAB      < 200.0
                                      )
                                    {
                                        Float_t pathA_prim,dcaA_prim,pathB_prim,dcaB_prim;

                                        helixA_glob = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
                                        helixB_glob = StPhysicalHelixD(trackB.gMom(),trackB.origin() + vectordiff,event_B_ana->bField()*MAGFIELDFACTOR, trackB.charge());
                                        fHelixAtoPointdca(vector_prim,helixA_glob,pathA_prim,dcaA_prim);
                                        fHelixAtoPointdca(vector_prim,helixB_glob,pathB_prim,dcaB_prim);

                                        StThreeVectorF vectornewA_prim,vectornewB_prim;
                                        vectornewA_prim     = helixA_glob.cat(pathA_prim);
                                        vectornewB_prim     = helixB_glob.cat(pathB_prim);

                                        vectornewA_prim     = trackA.gMom().mag()*vectornewA_prim/vectornewA_prim.mag(); // momentum vector at primary vertex
                                        vectornewB_prim     = trackB.gMom().mag()*vectornewB_prim/vectornewB_prim.mag(); // momentum vector at primary vertex

                                        Float_t phi_event_plane                 = -400.0;
                                        Float_t phi_event_plane_eta_gap         = -400.0;
                                        Float_t delta_phi_ME_AB_weight          = 0.0;
                                        Float_t delta_phi_ME_AB_weight_eta_gap  = 0.0;

                                        Int_t delta_phi_ME = 0;

                                        // Calculate the event plane anlges
                                        calc_event_plane_angles(2,rap,trackA,trackB,trackA,vectornewA_prim,vectornewB_prim,vectornewA_prim,ME_Flag,SE_ME_Flag,RunIdA,EventVertexXA,EventVertexYA,
                                                                EventVertexZA,RunIdB,EventVertexXB,EventVertexYB,EventVertexZB,
                                                                phi_event_plane,phi_event_plane_eta_gap,delta_phi_ME_AB_weight,delta_phi_ME_AB_weight_eta_gap);


                                        // check whether the event planes between event A and B are close to each other
                                        if(
                                           fabs(delta_phi_ME_AB_weight) < TMath::DegToRad()*30.0
                                           || (SE_ME_Flag == 0)
                                          )
                                        {
                                            delta_phi_ME = 1;
                                        }

                                        if(delta_phi_ME == 1)
                                        {
                                            Float_t p_xA_c   = vectornewA_prim.x();
                                            Float_t p_yA_c   = vectornewA_prim.y();
                                            Float_t p_tA_c   = sqrt(p_xA_c*p_xA_c + p_yA_c*p_yA_c);
                                            //Float_t etaA_c   = vectornewA_prim.pseudoRapidity();
                                            //Float_t phiA_c   = vectornewA_prim.phi();
                                            Double_t p_t_weightA = 1.0;
                                            if(p_tA_c < 2.0)  p_t_weightA = p_tA_c;
                                            if(p_tA_c >= 2.0) p_t_weightA = 2.0;
                                            //Float_t iQxA     = p_t_weightA*TMath::Cos(2.0*phiA_c);
                                            //Float_t iQyA     = p_t_weightA*TMath::Sin(2.0*phiA_c);

                                            Float_t p_xB_c   = vectornewB_prim.x();
                                            Float_t p_yB_c   = vectornewB_prim.y();
                                            Float_t p_tB_c   = sqrt(p_xB_c*p_xB_c + p_yB_c*p_yB_c);
                                            //Float_t etaB_c   = vectornewB_prim.pseudoRapidity();
                                            //Float_t phiB_c   = vectornewB_prim.phi();
                                            Double_t p_t_weightB = 1.0;
                                            if(p_tB_c < 2.0)  p_t_weightB = p_tB_c;
                                            if(p_tB_c >= 2.0) p_t_weightB = 2.0;
                                            //Float_t iQxB     = p_t_weightB*TMath::Cos(2.0*phiB_c);
                                            //Float_t iQyB     = p_t_weightB*TMath::Sin(2.0*phiB_c);

                                            Int_t pt_bin = 0;
                                            if(pt > 0.0 && pt < 0.2) pt_bin = 1;
                                            if(pt > 0.2 && pt < 0.4) pt_bin = 2;
                                            if(pt > 0.4 && pt < 0.8) pt_bin = 3;
                                            if(pt > 0.8 && pt < 1.2) pt_bin = 4;
                                            if(pt > 1.2 && pt < 1.7) pt_bin = 5;
                                            if(pt > 1.7 && pt < 2.3) pt_bin = 6;
                                            if(pt > 2.3 && pt < 3.0) pt_bin = 7;
                                            if(pt > 3.0 && pt < 4.0) pt_bin = 8;
                                            if(pt > 4.0 && pt < 5.0) pt_bin = 9;

                                            Int_t pt_bin_array[2] = {0,pt_bin};
                                            // Filling histograms
                                            for(Int_t pt_bin_loop = 0; pt_bin_loop < 2; pt_bin_loop++)
                                            {
                                                Int_t pt_bin_sel = pt_bin_array[pt_bin_loop];
                                                h_inv_mass_pipi[0][pt_bin_sel]->Fill(InvMassAB);
                                                if(dcaA < 0.5 && dcaB < 0.5) h_inv_mass_pipi[1][pt_bin_sel]->Fill(InvMassAB);
                                                if(flag_tof_KP == 1 && flag_tof_KM == 1) h_inv_mass_pipi[2][pt_bin_sel]->Fill(InvMassAB);
                                                if(dcaA < 0.5 && dcaB < 0.5 && flag_tof_KP == 1 && flag_tof_KM == 1) h_inv_mass_pipi[3][pt_bin_sel]->Fill(InvMassAB);
                                                if(flag_tof_KP == 1 && flag_tof_KM == 1 && fabs(nSigmaPionA) < 2.0 && fabs(nSigmaPionB) < 2.0) h_inv_mass_pipi[4][pt_bin_sel]->Fill(InvMassAB);
                                                if(dcaA < 0.5 && dcaB < 0.5 && flag_tof_KP == 1 && flag_tof_KM == 1 && fabs(nSigmaPionA) < 2.0 && fabs(nSigmaPionB) < 2.0) h_inv_mass_pipi[5][pt_bin_sel]->Fill(InvMassAB);
                                            }

                                            dummy_counter_A++;
                                        }
                                        //*********************************************************************************************
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(fAnalysisNum == 1) Tree_PhiMeson_v2  ->Fill();
        return 1;
    }
    else return 0;
}



Int_t Lambda_ppim_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs], Int_t PID_counter_Array_C[][N_max_PIDs],
                           Int_t PID_Array[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks], Int_t PID_Array_C[][N_max_PIDs][N_max_tracks],
                           StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B, StPicoAlexEvent* picoDst_C,
                           Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t Ana_Num,
                           Int_t SE_ME_Flag)
{
    // Au + Au -> Lambda + N + N
    //              |
    //              -> p + pi-

    // Event vertex information
    //cout << "Lambda analysis started" << endl;
    StThreeVectorF vector_prim,vector_primB,vectordiff;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist,ref_mult_frac;
    Int_t refMultA,refMultB,RunIdA,RunIdB;

    event_A_ana   = picoDst_A;
    event_B_ana   = picoDst_B;
    event_C_ana   = picoDst_C;

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();
    EventVertexZA     = event_A_ana->primaryVertex().z();
    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    refMultB          = refMultA;
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_A_ana->vzVpd();
    vector_prim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);
    Double_t pion_low_cut        = -0.5; //
    Double_t pion_high_cut       = 0.1;  // 0.05
    Double_t proton_low_cut      = 0.45;  // 0.76
    Double_t proton_high_cut     = 1.5;  // 1.12

    //cout << "refMultA = " << refMultA <<  endl;

    if(
       SE_ME_Flag == 1  // mixed event analysis for Lambdas
       && (fAnalysisNum == 2 || fAnalysisNum == 14)
      )
    {
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        RunIdB         = event_B_ana->runId();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex

        ref_mult_frac = fabs(100*2.0*((Float_t)(refMultA - refMultB))/((Float_t)(refMultA + refMultB)));
    }
    if(
       SE_ME_Flag == 1  // mixed event analysis for Sigma(1385)
       && (fAnalysisNum == 21)
      )
    {
        EventVertexXB  = event_C_ana->primaryVertex().x();
        EventVertexYB  = event_C_ana->primaryVertex().y();
        EventVertexZB  = event_C_ana->primaryVertex().z();
        refMultB       = event_C_ana->refMult();
        RunIdB         = event_C_ana->runId();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex

        ref_mult_frac = fabs(100*2.0*((Float_t)(refMultA - refMultB))/((Float_t)(refMultA + refMultB)));
    }

    //cout << "XA = " << EventVertexXA << ", XB = " << EventVertexXB << endl;

    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm
    //Float_t ME_ref_mult_frac_cut = 100.0; // 100% difference in reference multiplicity is allowed

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
       //&& ref_mult_frac < ME_ref_mult_frac_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
        //cout << "Triggers checked, Vertex checked..." << endl;
    }

    //cout << "Trigger word (A) = " << event_A_ana->triggerWord() << ", MB (A) = " << event_A_ana->isMinBias() << ", MB (B) = " << event_B_ana->isMinBias() <<
    //    ", vertex = {" << EventVertexXA << ", " << EventVertexYA << ", "<< EventVertexZA << "}" << endl;

    //

    //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
    //    << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;

    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]   > 0   // p
       && PID_counter_Array_B[Ana_Num][ParticleB] > 0   // pi-
       //&& (
       //    fAnalysisNum == 102
       //    && PID_counter_Array[Ana_Num][2] > 0   // e+
       //    && PID_counter_Array[Ana_Num][3] > 0   // e-
       //   )
       //&& event->getNVertices() >= 1
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       //&& event_A_ana   ->isMinBias()
       //&& event_B_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
        //    << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB, helixC, helixD;
        StThreeVectorF vectorA, vectorB, vectorC, vectoratsA, vectoratsB, vectorAB, vector_primAB, vectornewA, vectornewB, vectornewA_lin, vectornewB_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY, dirY_lin;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackA_lin, ltrackB_lin, ltrackB2;
        TLorentzVector ltrackA_pip;

        if(fAnalysisNum == 102)
        {
            //cout << "Fill event" << endl;
            alexV0_event_A.clearTrackList();
            alexV0_event_A.setx(EventVertexXA);
            alexV0_event_A.sety(EventVertexYA);
            alexV0_event_A.setz(EventVertexZA);
            alexV0_event_A.setid(RunIdA);
            alexV0_event_A.setmult(refMultA);
            alexV0_event_A.setn_prim(n_primaries);
            alexV0_event_A.setn_non_prim(n_non_primaries);
            alexV0_event_A.setn_tof_prim(n_tofmatch_prim);
            alexV0_event_A.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
            alexV0_event_A.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
            alexV0_event_A.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
            alexV0_event_A.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
            alexV0_event_A.setEP_Qx_ptw(EP_Qx_ptw);
            alexV0_event_A.setEP_Qy_ptw(EP_Qy_ptw);
            alexV0_event_A.setQtracks_eta_pos(Qtracks_used_eta_pos);
            alexV0_event_A.setQtracks_eta_neg(Qtracks_used_eta_neg);
            alexV0_event_A.setQtracks_full(Qtracks_used);
            alexV0_event_A.setZDCx(ZDCx);
            alexV0_event_A.setBBCx(BBCx);
            alexV0_event_A.setvzVpd(vzVpd);
        }

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // p candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction

            //cout << "nHitsFitA = " << nHitsFitA << ", chargeA = " << trackA.charge() << endl;

            Float_t Mass2A        = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }

            Float_t MassA       = SquareRoot(Mass2A);
            //Float_t TPCdEdxA    = trackA.dEdx(); // TPC dE/dx
            Float_t nSigmaPA    = trackA.nSigmaProton();
            //Float_t TofA        = trackA.btof();
            if(nHitsPossA <= 0) nHitsPossA = 10000.0;

            Double_t dcaA_cut     = 0.4;
            Double_t dcaB_cut     = 0.6;
            Double_t VerdistX_cut = 2.0;
            Double_t VerdistY_cut = 1.5;
            Int_t    flag_proton  = 0;
            Int_t    flag_pion    = 0;
            Int_t    flag_tof_proton  = 0;  // proton
            Int_t    flag_tof_pion    = 0;  // pion
            Int_t    flag_tof_ep      = 0;  // e+
            Int_t    flag_tof_em      = 0;  // e-

            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                flag_tof_proton = 1;
            }
            else
            {
                flag_tof_proton = 0;
            }

            if(flag_tof_proton == 1 && Mass2A < 1.5 && Mass2A > 0.45)
            {
                dcaA_cut     = 0.1;
                flag_proton  = 1;
            }
            else
            {
                flag_proton = 0;
            }

            //cout << "i = " << i << ", dcaA = " << dcaA << ", fitA = " << nHitsFitA << ", possA = " << nHitsPossA << ", pA =" << MomentumA << endl;

            if(
               dcaA         > dcaA_cut // 0.5
               && nHitsFitA > 14 // 14
               && MomentumA > 0.1
               && MomentumA < 10.0
               && (nHitsFitA/nHitsPossA) > 0.52
              )
            {

                for(Int_t j = 0; j < PID_counter_Array_B[Ana_Num][ParticleB]; j++)  // pi- candidates
                {
                    Int_t trackB_num = PID_Array_B[Ana_Num][ParticleB][j];
                    StPicoAlexTrack trackB    = *picoDst_B->track( trackB_num );
                    vectorB = trackB.origin();
                    if(fAnalysisNum == 2 || fAnalysisNum == 14) vectorB = vectorB + vectordiff; // vectordiff == {0,0,0} for same event, don't do it for Sigma(1385) analysis
                    helixB = StPhysicalHelixD(trackB.gMom(),vectorB,event_B_ana->bField()*MAGFIELDFACTOR,trackB.charge());
                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    //Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction

                    //cout << "nHitsFitB = " << nHitsFitA << ", chargeB = " << trackB.charge() << endl;

                    Float_t Mass2B        = -100.0;
                    // calculate mass2
                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                    }


                    Float_t MassB       = SquareRoot(Mass2B);
                    //Float_t TPCdEdxB    = trackB.dEdx(); // TPC dE/dx
                    Float_t nSigmaPionB = trackB.nSigmaPion();
                    //Float_t TofB        = trackB.btof();
                    if(nHitsPossB <= 0) nHitsPossB = 10000.0;

                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        flag_tof_pion = 1;
                    }
                    else
                    {
                        flag_tof_pion = 0;
                    }

                    if(flag_tof_pion == 1 && Mass2B < 0.05 && Mass2B > -0.05) // time-of-flight is available and particle B is a pion
                    {
                        flag_pion = 1;
                    }
                    else
                    {
                        flag_pion = 0;
                    }
                    if(flag_proton == 0 && flag_pion == 0)
                    {
                        dcaA_cut     = 0.5;
                        dcaB_cut     = 1.5;
                        VerdistX_cut = 3.5;
                        VerdistY_cut = 0.8;
                    }
                    if(flag_proton == 0 && flag_pion == 1)
                    {
                        dcaA_cut     = 0.4;
                        dcaB_cut     = 1.2;
                        VerdistX_cut = 3.0;
                        VerdistY_cut = 0.9;
                    }
                    if(flag_proton == 1)
                    {
                        dcaA_cut     = 0.1;
                        dcaB_cut     = 0.6;
                        VerdistX_cut = 2.0;
                        VerdistY_cut = 1.5;
                    }

                    //cout << "j = " << j << ", dcaB = " << dcaB << ", fitB = " << nHitsFitB << ", possB = " << nHitsPossB << ", pB =" << MomentumB << endl;

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       //&& dcaB > dcaA
                       && dcaA      > dcaA_cut // 1.0
                       && dcaB      > dcaB_cut // 1.0
                       && nHitsFitB > 14 // 14
                       && MomentumB > 0.1
                       && MomentumB < 10.0
                       && (nHitsFitB/nHitsPossB) > 0.52
                      )
                    {
                        Float_t pathA_f, pathB_f, dcaAB_f;
                        Float_t pathA_test = -1.0;
                        Float_t pathB_test = -1.0;
                        Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);

                        //cout << "pathA_test = " << pathA_test << ", pathB_test = " << pathB_test << endl;

                        vectoratsA     = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
                        vectoratsB     = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
                        vectorAB       = vectoratsA+vectoratsB;
                        vectorAB       = vectorAB/2.0; // decay vertex

                        Double_t vex_est_x = vectorAB.x();
                        Double_t vex_est_y = vectorAB.y();
                        Double_t vex_est_z = vectorAB.z();

                        StThreeVectorF baseA,dirA,baseB,dirB;
                        baseA = helixA.at(pathA_test);
                        baseB = helixB.at(pathB_test);
                        dirA  = helixA.at(pathA_test-2.0) - helixA.at(pathA_test+2.0);
                        dirB  = helixB.at(pathB_test-2.0) - helixB.at(pathB_test+2.0);

                        StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

                        Double_t vex_lin_x = vectorAB_lin.x();
                        Double_t vex_lin_y = vectorAB_lin.y();
                        Double_t vex_lin_z = vectorAB_lin.z();

                        Double_t dcaAB_lin = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
                        StThreeVectorF vectorABtoPrim_lin = vectorAB_lin - vector_prim; // vector primary vertex to decay vertex
                        Float_t VerdistX_lin = vectorABtoPrim_lin.mag(); // distance between primary vertex and decay verte


                        // calculate the scalar product with the approximated secondary vertex position
                        vectornewA_lin     = helixA.cat(pathA_test); // direction vector at dca for helixA
                        vectornewB_lin     = helixB.cat(pathB_test); // direction vector at dca for helixB
                        vectornewA_lin = MomentumA*vectornewA_lin/vectornewA_lin.mag(); // new momentum vector at decay vertex
                        vectornewB_lin = MomentumB*vectornewB_lin/vectornewB_lin.mag(); // new momentum vector at decay vertex
                        ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),0.93827203);
                        ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),0.13957018);
                        TLorentzVector trackAB_lin      = ltrackA_lin+ltrackB_lin; // mother particle
                        dirY_lin.set(trackAB_lin.Px(),trackAB_lin.Py(),trackAB_lin.Pz());
                        dirY_lin = dirY_lin/dirY_lin.mag();
                        Double_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim_lin/vectorABtoPrim_lin.mag());

                        //cout << "i = " << i << ", j = " << j << ", VerdistX_lin = " << VerdistX_lin << ", dcaAB_lin = " << dcaAB_lin << ", scalarProduct_lin = " << scalarProduct_lin << endl;

                        //if( (VerdistX_lin > 3.0 && dcaAB_lin < 1.5 && fDCA_Helix_out == 1) ||  fDCA_Helix_out == 0 )
                        if( VerdistX_lin > VerdistX_cut && dcaAB_lin < 1.5 && scalarProduct_lin > 0.0 )
                        {

                            if(fDCA_Helix_out == 1)
                            {
                                fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
                            }
                            else
                            {
                                fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
                            }

                            vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                            vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                            vectorAB       = vectoratsA+vectoratsB;
                            vectorAB       = vectorAB/2.0; // decay vertex

                            Double_t vex_x = vectorAB.x();
                            Double_t vex_y = vectorAB.y();
                            Double_t vex_z = vectorAB.z();

                            //cout << "Vertex = {" << vex_x << ", " << vex_y << ", " << vex_z << "}" <<
                            //    ", Vertex_est = {" << vex_est_x << ", " << vex_est_y << ", " << vex_est_z << "}"<<
                            //    ", Vertex_lin = {" << vex_lin_x << ", " << vex_lin_y << ", " << vex_lin_z << "}"
                            //    << ", fDCA_Helix_out = " << fDCA_Helix_out << endl;

                            Vertex_lin_NTDataArray[0]    =(Float_t)vex_x;
                            Vertex_lin_NTDataArray[1]    =(Float_t)vex_y;
                            Vertex_lin_NTDataArray[2]    =(Float_t)vex_z;
                            Vertex_lin_NTDataArray[3]    =(Float_t)vex_est_x;
                            Vertex_lin_NTDataArray[4]    =(Float_t)vex_est_y;
                            Vertex_lin_NTDataArray[5]    =(Float_t)vex_est_z;
                            Vertex_lin_NTDataArray[6]    =(Float_t)vex_lin_x;
                            Vertex_lin_NTDataArray[7]    =(Float_t)vex_lin_y;
                            Vertex_lin_NTDataArray[8]    =(Float_t)vex_lin_z;
                            Vertex_lin_NTDataArray[9]    =(Float_t)fDCA_Helix_out;
                            Vertex_lin_NTDataArray[10]   =(Float_t)MomentumA;
                            Vertex_lin_NTDataArray[11]   =(Float_t)MomentumB;
                            Vertex_lin_NTDataArray[12]   =(Float_t)dcaAB_f;
                            //Vertex_lin_NT->Fill(Vertex_lin_NTDataArray);

                            vectorABtoPrim = vectorAB - vector_prim; // vector primary vertex to decay vertex
                            Float_t VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

                            vectornewA     = helixA.cat(pathA_f); // direction vector at dca for helixA
                            vectornewB     = helixB.cat(pathB_f); // direction vector at dca for helixB

                            vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
                            vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex

                            ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.93827203);
                            ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),0.13957018);
                            ltrackA_pip.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.13957018);

                            //Float_t thetaA = vectoratsA.theta();
                            //Float_t thetaB = vectoratsB.theta();

                            // Missing mass and invariant mass calculations
                            TLorentzVector trackAB      = ltrackA+ltrackB; // mother particle
                            TLorentzVector trackAB_K0S  = ltrackA_pip+ltrackB; // mother particle
                            Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                            Double_t InvMassAB_K0S      = trackAB_K0S.M(); // invariant mass of mother particle
                            Float_t MomentumAB          = trackAB.P(); // momentum of mother particle
                            Float_t BetaAB = TMath::Sqrt(1.0/(1.0+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            dirY = dirY/dirY.mag();

                            Double_t scalarProduct = dirY.dot(vectorABtoPrim/vectorABtoPrim.mag());
                            //cout << "scalarProduct = " << scalarProduct << endl;

                            baseY = vectorAB;
                            Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vector_prim);

                            Float_t pt          = trackAB.Pt();  // Transverse momentum of mother particle
                            Float_t rap         = trackAB.Rapidity(); // Rapidity of mother particle

                            Float_t phiAB   = dirY.phi();
                            Float_t thetaAB = dirY.theta();

                            // beta correction Lambda ========================================================================================
                            Float_t BetaACorr  = BetaA;
                            Float_t Mass2ACorr = Mass2A;
                            Float_t BetaBCorr  = BetaB;
                            Float_t Mass2BCorr = Mass2B;
                            Float_t MassACorr  = MassA;
                            Float_t MassBCorr  = MassB;
                            Float_t PathAB = vectorABtoPrim.mag(); // Lambda uncharged !!


                            if((trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0) || (trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0))
                            {
                                if( debug_flag )
                                {
                                    cout << "------------------------------------------------- LAMBDA --------------------------------------------------------------------------" << endl;
                                    cout << "PathAB = " << PathAB << "\tInvMassAB = " << InvMassAB << "\tMomentumAB = " << MomentumAB
                                        << "\tBetaAB = " << BetaAB << "\tTofAB = " << PathAB / (BetaAB*29.9792458)<< endl;
                                }
                            }
                            if(flag_tof_proton == 1)
                            {
                                BetaACorr = correctBeta4SingleDecay(trackA,trackAB,helixA,vectorAB,PathAB);
                                Mass2ACorr = MomentumA * MomentumA * (1.0/(BetaACorr*BetaACorr)-1.0);
                                MassACorr = SquareRoot(Mass2ACorr);
                                if( debug_flag ) cout << "A) MassA = " << MassA << "\tMassACorr = " <<  MassACorr << endl;
                            }
                            else Mass2ACorr = -100;

                            if(flag_tof_pion == 1)
                            {
                                BetaBCorr = correctBeta4SingleDecay(trackB,trackAB,helixB,vectorAB,PathAB);
                                Mass2BCorr = MomentumB * MomentumB * (1.0/(BetaBCorr*BetaBCorr)-1.0);
                                MassBCorr = SquareRoot(Mass2BCorr);
                                if( debug_flag ) cout << "B) MassB = " << MassB << "\tMassBCorr = " <<  MassBCorr << endl;
                            }
                            else Mass2BCorr = -100;


                            if(
                               InvMassAB        > 1.06
                               && InvMassAB     < 1.19
                               && VerdistX      > VerdistX_cut  // 3.5
                               && dcaAB_f       < 1.5  // 1.5
                               && VerdistY      < VerdistY_cut  // 1.5
                               && scalarProduct > 0.0
                               && ((flag_tof_proton == 0) || (Mass2ACorr > 100.0) || (Mass2ACorr  > proton_low_cut && Mass2ACorr < proton_high_cut))
                               && ((flag_tof_pion   == 0) || (Mass2BCorr > 100.0) || (Mass2BCorr  > pion_low_cut && Mass2BCorr   < pion_high_cut))

                              )
                            {
                                //cout << "Mass checked, topology checked --> event accepted" << endl;
                                //****************** Event plane calculations ***********************************************************
                                // Calculate the same vectors as for the event plane -> no V0!
                                Float_t pathA_prim,dcaA_prim,pathB_prim,dcaB_prim;
                                fHelixAtoPointdca(vector_prim,helixA,pathA_prim,dcaA_prim);
                                fHelixAtoPointdca(vector_prim,helixB,pathB_prim,dcaB_prim);

                                StThreeVectorF vectornewA_prim,vectornewB_prim;
                                vectornewA_prim     = helixA.cat(pathA_prim);
                                vectornewB_prim     = helixB.cat(pathB_prim);

                                vectornewA_prim     = MomentumA*vectornewA_prim/vectornewA_prim.mag(); // momentum vector at primary vertex
                                vectornewB_prim     = MomentumB*vectornewB_prim/vectornewB_prim.mag(); // momentum vector at primary vertex

                                Int_t delta_phi_ME = 0;
                                Float_t phi_event_plane                 = -400.0;
                                Float_t phi_event_plane_eta_gap         = -400.0;
                                Float_t delta_phi_ME_AB_weight          = 0.0;
                                Float_t delta_phi_ME_AB_weight_eta_gap  = 0.0;

                                if(fAnalysisNum == 2 || fAnalysisNum == 14)
                                {
                                    // Calculate the event plane anlges
                                    calc_event_plane_angles(2,rap,trackA,trackB,trackA,vectornewA_prim,vectornewB_prim,vectornewA_prim,ME_Flag,SE_ME_Flag,RunIdA,EventVertexXA,EventVertexYA,
                                                            EventVertexZA,RunIdB,EventVertexXB,EventVertexYB,EventVertexZB,
                                                            phi_event_plane,phi_event_plane_eta_gap,delta_phi_ME_AB_weight,delta_phi_ME_AB_weight_eta_gap);
                                }
                                else delta_phi_ME = 1;

                                // check whether the event planes between event A and B are close to each other
                                if(
                                   fabs(delta_phi_ME_AB_weight) < TMath::DegToRad()*30.0
                                   || (SE_ME_Flag == 0)
                                  )
                                {
                                    delta_phi_ME = 1;
                                }
                                //******************************************************************************************************


                                if(delta_phi_ME == 1)
                                {
                                    Float_t p_xA_c   = vectornewA_prim.x();
                                    Float_t p_yA_c   = vectornewA_prim.y();
                                    Float_t p_tA_c   = sqrt(p_xA_c*p_xA_c + p_yA_c*p_yA_c);
                                    Float_t etaA_c   = vectornewA_prim.pseudoRapidity();
                                    Float_t phiA_c   = vectornewA_prim.phi();
                                    Double_t p_t_weightA = 1.0;
                                    if(p_tA_c < 2.0)  p_t_weightA = p_tA_c;
                                    if(p_tA_c >= 2.0) p_t_weightA = 2.0;
                                    Float_t iQxA     = p_t_weightA*TMath::Cos(2.0*phiA_c);
                                    Float_t iQyA     = p_t_weightA*TMath::Sin(2.0*phiA_c);

                                    Float_t p_xB_c   = vectornewB_prim.x();
                                    Float_t p_yB_c   = vectornewB_prim.y();
                                    Float_t p_tB_c   = sqrt(p_xB_c*p_xB_c + p_yB_c*p_yB_c);
                                    Float_t etaB_c   = vectornewB_prim.pseudoRapidity();
                                    Float_t phiB_c   = vectornewB_prim.phi();
                                    Double_t p_t_weightB = 1.0;
                                    if(p_tB_c < 2.0)  p_t_weightB = p_tB_c;
                                    if(p_tB_c >= 2.0) p_t_weightB = 2.0;
                                    Float_t iQxB     = p_t_weightB*TMath::Cos(2.0*phiB_c);
                                    Float_t iQyB     = p_t_weightB*TMath::Sin(2.0*phiB_c);

                                    // Filling Ntuple
                                    Lambda_X_NTDataArray[0]      =(Float_t)InvMassAB;
                                    Lambda_X_NTDataArray[1]      =(Float_t)InvMassAB_K0S;

                                    Lambda_X_NTDataArray[2]      =(Float_t)Mass2ACorr;
                                    Lambda_X_NTDataArray[3]      =(Float_t)Mass2BCorr;
                                    Lambda_X_NTDataArray[4]      =(Float_t)nSigmaPA;
                                    Lambda_X_NTDataArray[5]      =(Float_t)nSigmaPionB;
                                    Lambda_X_NTDataArray[6]      =(Float_t)MomentumA;
                                    Lambda_X_NTDataArray[7]      =(Float_t)MomentumB;
                                    Lambda_X_NTDataArray[8]      =(Float_t)dcaA;
                                    Lambda_X_NTDataArray[9]      =(Float_t)dcaB;
                                    Lambda_X_NTDataArray[10]     =(Float_t)iQxA;
                                    Lambda_X_NTDataArray[11]     =(Float_t)iQyA;
                                    Lambda_X_NTDataArray[12]     =(Float_t)iQxB;
                                    Lambda_X_NTDataArray[13]     =(Float_t)iQyB;
                                    Lambda_X_NTDataArray[14]     =(Float_t)etaA_c;
                                    Lambda_X_NTDataArray[15]     =(Float_t)etaB_c;

                                    Lambda_X_NTDataArray[16]     =(Float_t)refMultA; //1
                                    Lambda_X_NTDataArray[17]     =(Float_t)dcaAB_f;
                                    Lambda_X_NTDataArray[18]     =(Float_t)VerdistX;
                                    Lambda_X_NTDataArray[19]     =(Float_t)VerdistY;
                                    Lambda_X_NTDataArray[20]     =(Float_t)pt;
                                    Lambda_X_NTDataArray[21]     =(Float_t)rap;
                                    Lambda_X_NTDataArray[22]     =(Float_t)phiAB;
                                    Lambda_X_NTDataArray[23]     =(Float_t)thetaAB;
                                    Lambda_X_NTDataArray[24]     =(Float_t)phi_event_plane;
                                    Lambda_X_NTDataArray[25]     =(Float_t)phi_event_plane_eta_gap;
                                    Lambda_X_NTDataArray[26]     =(Float_t)delta_phi_ME_AB_weight;
                                    Lambda_X_NTDataArray[27]     =(Float_t)RunIdA;  //1
                                    Lambda_X_NTDataArray[28]     =(Float_t)n_primaries; //1
                                    Lambda_X_NTDataArray[29]     =(Float_t)n_non_primaries; //1
                                    Lambda_X_NTDataArray[30]     =(Float_t)n_tofmatch_prim; //1
                                    Lambda_X_NTDataArray[31]     =(Float_t)scalarProduct;
                                    Lambda_X_NTDataArray[32]     =(Float_t)EventVertexXA; //1
                                    Lambda_X_NTDataArray[33]     =(Float_t)EventVertexYA; //1
                                    Lambda_X_NTDataArray[34]     =(Float_t)EventVertexZA; //1
                                    Lambda_X_NTDataArray[35]     =(Float_t)EP_Qx_eta_pos_ptw;  //1
                                    Lambda_X_NTDataArray[36]     =(Float_t)EP_Qy_eta_pos_ptw;  //1
                                    Lambda_X_NTDataArray[37]     =(Float_t)EP_Qx_eta_neg_ptw;  //1
                                    Lambda_X_NTDataArray[38]     =(Float_t)EP_Qy_eta_neg_ptw;  //1
                                    Lambda_X_NTDataArray[39]     =(Float_t)EP_Qx_ptw; //1
                                    Lambda_X_NTDataArray[40]     =(Float_t)EP_Qy_ptw; //1
                                    Lambda_X_NTDataArray[41]     =(Float_t)Qtracks_used_eta_pos; //1
                                    Lambda_X_NTDataArray[42]     =(Float_t)Qtracks_used_eta_neg; //1
                                    Lambda_X_NTDataArray[43]     =(Float_t)Qtracks_used; //1
                                    Lambda_X_NTDataArray[44]     =(Float_t)ZDCx; //1
                                    Lambda_X_NTDataArray[45]     =(Float_t)vzVpd; //1

                                    if(fAnalysisNum == 2 || fAnalysisNum == 14)
                                    {
                                        Lambda_X_NT->Fill(Lambda_X_NTDataArray);
                                    }

                                    if(ME_Flag == 1)
                                    {
                                        Lambda_ME_counter++;
                                    }



                                    if(fAnalysisNum == 102)
                                    {
                                        StThreeVectorF primdirC, primdirD, vectorC, vectorD, vectorCD;
                                        TLorentzVector ltrackC, ltrackD;
                                        StPhysicalHelixD helixC_glob, helixD_glob;

                                        for(Int_t iep = 0; iep < PID_counter_Array[Ana_Num][2]; iep++) // e+ candidates
                                        {
                                            Int_t trackC_num = PID_Array[Ana_Num][2][iep];
                                            StPicoAlexTrack trackC = *picoDst_A->track( trackC_num );

                                            Float_t MomentumC   = trackC.pMom().mag();
                                            Float_t dcaC        = trackC.dca();
                                            Float_t nHitsPossC  = trackC.nHitsMax();
                                            Float_t nHitsFitC   = trackC.nHitsFit();
                                            Float_t BetaC       = trackC.btofBeta();  // Velocity after time-of-flight reconstruction
                                            //Float_t PolarityC   = trackC.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                            //Float_t TPCdEdxC    = trackC.dEdx(); // Combined inner and outer MDC dE/dx
                                            Float_t nSigmaElecC = trackC.nSigmaElectron();

                                            // calculate mass2
                                            Float_t Mass2C      = -100.0;
                                            if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                                            {
                                                Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                                                flag_tof_ep = 1;
                                            }
                                            else
                                            {
                                                flag_tof_ep = 0;
                                            }

                                            if(
                                               dcaC           < 1.5
                                               && nHitsFitC   > 14
                                               && nHitsPossC  > 0
                                               && (nHitsFitC/nHitsPossC) > 0.52
                                               && MomentumC   > 0.02
                                               && MomentumC   < 10.0
                                               && nSigmaElecC < 2.5/nsigma_scaling_fac
                                               && nSigmaElecC > -2.5/nsigma_scaling_fac
                                               && Mass2C > -0.1 && Mass2C < 0.01
                                               && (trackC_num != trackA_num)
                                               && (trackC_num != trackB_num)
                                              )
                                            {
                                                helixC   = StPhysicalHelixD(trackC.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackC.charge());    // gMom or pMom?
                                                primdirC = helixC.cat(helixC.pathLength(vector_prim)); //
                                                primdirC = MomentumC*primdirC/primdirC.mag();
                                                ltrackC.SetXYZM(primdirC.x(),primdirC.y(),primdirC.z(),0.00051099892);

                                                Float_t ptC         = ltrackC.Pt();
                                                Float_t pzC         = ltrackC.Pz();
                                                Float_t etaC        = primdirC.pseudoRapidity();

                                                if(
                                                   fabs(etaC)     < 1.0
                                                   && ptC         > 0.1
                                                   && ptC         < 10.0
                                                  )
                                                {
                                                    for(Int_t iem = 0; iem < PID_counter_Array[Ana_Num][3]; iem++)  // e- candidates
                                                    {
                                                        Int_t trackD_num = PID_Array[Ana_Num][3][iem];

                                                        if(
                                                           (trackD_num != trackA_num)
                                                           && (trackD_num != trackB_num)
                                                           && (trackD_num != trackC_num)
                                                          )
                                                        {
                                                            StPicoAlexTrack trackD = *picoDst_A->track( trackD_num );

                                                            Float_t MomentumD     = trackD.pMom().mag();
                                                            Float_t dcaD          = trackD.dca();
                                                            Float_t nHitsPossD    = trackD.nHitsMax();
                                                            Float_t nHitsFitD     = trackD.nHitsFit();
                                                            Float_t BetaD         = trackD.btofBeta();  // Velocity after time-of-flight reconstruction
                                                            //Float_t PolarityD     = trackD.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                                            //Float_t TPCdEdxD      = trackD.dEdx(); // Combined inner and outer MDC dE/dx
                                                            Float_t nSigmaElecD   = trackD.nSigmaElectron();

                                                            // calculate mass2
                                                            Float_t Mass2D        = -100.0;
                                                            if(trackD.btofMatchFlag() > 0 && trackD.btof() != 0 && BetaD != 0)
                                                            {
                                                                Mass2D = MomentumD*MomentumD*(1.0/(BetaD*BetaD) - 1.0);
                                                                flag_tof_em = 1;
                                                            }
                                                            else
                                                            {
                                                                flag_tof_em = 0;
                                                            }

                                                            if(
                                                               dcaD           < 1.5
                                                               && nHitsFitD   > 14
                                                               && nHitsPossD  > 0
                                                               && (nHitsFitD/nHitsPossD) > 0.52
                                                               && MomentumD   > 0.02
                                                               && MomentumD   < 10.0
                                                               && nSigmaElecD < 2.5/nsigma_scaling_fac
                                                               && nSigmaElecD > -2.5/nsigma_scaling_fac
                                                               && Mass2D > -0.1 && Mass2D < 0.01
                                                              )
                                                            {
                                                                helixD = StPhysicalHelixD(trackD.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackD.charge());
                                                                primdirD = helixD.cat(helixD.pathLength(vector_prim)); //
                                                                primdirD = MomentumD*primdirD/primdirD.mag();
                                                                ltrackD.SetXYZM(primdirD.x(),primdirD.y(),primdirD.z(),0.00051099892);

                                                                Float_t ptD           = ltrackD.Pt();
                                                                Float_t pzD           = ltrackD.Pz();
                                                                Float_t etaD          = primdirD.pseudoRapidity();

                                                                if(
                                                                   fabs(etaD)     < 1.0
                                                                   && ptD         > 0.1
                                                                   && ptD         < 10.0
                                                                  )
                                                                {
                                                                    // Invariant mass calculations
                                                                    TLorentzVector trackCD      = ltrackC+ltrackD;
                                                                    Double_t InvMassCD          = trackCD.M();
                                                                    TLorentzVector trackABCD    = ltrackA+ltrackB+ltrackC+ltrackD;
                                                                    Double_t InvMassABCD        = trackABCD.M();
                                                                    vectorCD.set(trackCD.X(),trackCD.Y(),trackCD.Z());
                                                                    //Float_t phiCD   = vectorCD.phi();
                                                                    //Float_t thetaCD = vectorCD.theta();

                                                                    //Float_t DeltaDipAngleCD = TMath::ACos((ptC*ptD+pzC*pzD)/(MomentumC*MomentumD));  // To remove electron conversion background, see Phys.Rev.C79:064903,2009
                                                                    Float_t pt              = trackCD.Pt();  // Transverse momentum of mother particle
                                                                    Float_t rap             = trackCD.Rapidity(); // Rapidity of mother particle

                                                                    if(
                                                                       InvMassCD      < 10.0
                                                                      )
                                                                    {
                                                                        Float_t pathC_prim,dcaC_prim,pathD_prim,dcaD_prim;

                                                                        helixC_glob = StPhysicalHelixD(trackC.gMom(),trackC.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackC.charge());
                                                                        helixD_glob = StPhysicalHelixD(trackD.gMom(),trackD.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackD.charge());
                                                                        fHelixAtoPointdca(vector_prim,helixC_glob,pathC_prim,dcaC_prim);
                                                                        fHelixAtoPointdca(vector_prim,helixD_glob,pathD_prim,dcaD_prim);

                                                                        StThreeVectorF vectornewC_prim,vectornewD_prim;
                                                                        vectornewC_prim     = helixC_glob.cat(pathC_prim);
                                                                        vectornewD_prim     = helixD_glob.cat(pathD_prim);

                                                                        vectornewC_prim     = trackC.gMom().mag()*vectornewC_prim/vectornewC_prim.mag(); // momentum vector at primary vertex
                                                                        vectornewD_prim     = trackD.gMom().mag()*vectornewD_prim/vectornewD_prim.mag(); // momentum vector at primary vertex


                                                                        Float_t p_xC_c   = vectornewC_prim.x();
                                                                        Float_t p_yC_c   = vectornewC_prim.y();
                                                                        Float_t p_tC_c   = sqrt(p_xC_c*p_xC_c + p_yC_c*p_yC_c);
                                                                        Float_t etaC_c   = vectornewC_prim.pseudoRapidity();
                                                                        Float_t phiC_c   = vectornewC_prim.phi();
                                                                        Double_t p_t_weightC = 1.0;
                                                                        if(p_tC_c < 2.0)  p_t_weightC = p_tC_c;
                                                                        if(p_tC_c >= 2.0) p_t_weightC = 2.0;
                                                                        Float_t iQxC     = p_t_weightC*TMath::Cos(2.0*phiC_c);
                                                                        Float_t iQyC     = p_t_weightC*TMath::Sin(2.0*phiC_c);

                                                                        Float_t p_xD_c   = vectornewD_prim.x();
                                                                        Float_t p_yD_c   = vectornewD_prim.y();
                                                                        Float_t p_tD_c   = sqrt(p_xD_c*p_xD_c + p_yD_c*p_yD_c);
                                                                        Float_t etaD_c   = vectornewD_prim.pseudoRapidity();
                                                                        Float_t phiD_c   = vectornewD_prim.phi();
                                                                        Double_t p_t_weightD = 1.0;
                                                                        if(p_tD_c < 2.0)  p_t_weightD = p_tD_c;
                                                                        if(p_tD_c >= 2.0) p_t_weightD = 2.0;
                                                                        Float_t iQxD     = p_t_weightD*TMath::Cos(2.0*phiD_c);
                                                                        Float_t iQyD     = p_t_weightD*TMath::Sin(2.0*phiD_c);

                                                                        //cout << "Fill track, InvMassABCD = " << InvMassABCD << endl;
                                                                        alexV0_track_A = alexV0_event_A.createTrack();
                                                                        alexV0_track_A->setm2A(Mass2ACorr);
                                                                        alexV0_track_A->setm2B(Mass2BCorr);
                                                                        alexV0_track_A->setm2C(Mass2C);
                                                                        alexV0_track_A->setm2D(Mass2D);
                                                                        alexV0_track_A->setnsA(nSigmaPA);
                                                                        alexV0_track_A->setnsB(nSigmaPionB);
                                                                        alexV0_track_A->setnsC(nSigmaElecC);
                                                                        alexV0_track_A->setnsD(nSigmaElecD);
                                                                        alexV0_track_A->setdcaA(dcaA);
                                                                        alexV0_track_A->setdcaB(dcaB);
                                                                        alexV0_track_A->setdcaC(dcaC);
                                                                        alexV0_track_A->setdcaD(dcaD);
                                                                        alexV0_track_A->setiQxA(iQxA);
                                                                        alexV0_track_A->setiQyA(iQyA);
                                                                        alexV0_track_A->setiQxB(iQxB);
                                                                        alexV0_track_A->setiQyB(iQyB);
                                                                        alexV0_track_A->setiQxC(iQxC);
                                                                        alexV0_track_A->setiQyC(iQyC);
                                                                        alexV0_track_A->setiQxD(iQxD);
                                                                        alexV0_track_A->setiQyD(iQyD);
                                                                        alexV0_track_A->setetaA(etaA_c);
                                                                        alexV0_track_A->setetaB(etaB_c);
                                                                        alexV0_track_A->setetaC(etaC_c);
                                                                        alexV0_track_A->setetaD(etaD_c);
                                                                        alexV0_track_A->setInvAB(InvMassAB);
                                                                        alexV0_track_A->setInvCD(InvMassCD);
                                                                        alexV0_track_A->setInvABC(-1.0);
                                                                        alexV0_track_A->setInvABCD(InvMassABCD);
                                                                        alexV0_track_A->setInvAB_miss(InvMassAB_K0S);
                                                                        alexV0_track_A->setInvABC_miss(-1.0);
                                                                        alexV0_track_A->setdcaAB(dcaAB_f);
                                                                        alexV0_track_A->setdcaBC(-1.0);
                                                                        alexV0_track_A->setdcaCD(-1.0);
                                                                        alexV0_track_A->setdcaABC(-1.0);
                                                                        alexV0_track_A->setVerdistX(VerdistX);
                                                                        alexV0_track_A->setVerdistY(VerdistY);
                                                                        alexV0_track_A->setVerdistX2(-1.0);
                                                                        alexV0_track_A->setVerdistY2(-1.0);
                                                                        alexV0_track_A->setpt(pt);
                                                                        alexV0_track_A->setrap(rap);
                                                                        alexV0_track_A->setphi(phiAB);
                                                                        alexV0_track_A->settheta(thetaAB);
                                                                        alexV0_track_A->setPsi_ep(phi_event_plane);
                                                                        alexV0_track_A->setPsi_ep_eta(phi_event_plane_eta_gap);
                                                                        alexV0_track_A->setPsi_diff_ME(delta_phi_ME_AB_weight);
                                                                        alexV0_track_A->setscal_prod(scalarProduct);
                                                                        alexV0_track_A->setscal_prod2(-1.0);
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }








                                    if(fAnalysisNum == 21)
                                    {
                                        //********************************************************************************************************************
                                        // Sigma(1385) analysis
                                        if(
                                           fAnalysisNum == 21
                                           && InvMassAB > 1.1157-3.0*0.00129 // 1.11223
                                           && InvMassAB < 1.1157+3.0*0.00129 // 1.11881
                                          )
                                        {
                                            for(Int_t c_sw = 0; c_sw < 2; c_sw++) // switch between pi+ and pi-
                                            {
                                                Int_t PID_C = 8;
                                                if(c_sw == 0) PID_C = 8; // pi+
                                                if(c_sw == 1) PID_C = 9; // pi-

                                                Double_t inv_mass_sign = 1.0;
                                                if(PID_C == 8)  // pi+
                                                {
                                                    inv_mass_sign = 1.0;
                                                }
                                                if(PID_C == 9)  // pi-
                                                {
                                                    inv_mass_sign = -1.0;
                                                }

                                                for(Int_t ic = 0; ic < PID_counter_Array_C[Ana_Num][PID_C]; ic++) // pi+ candidates
                                                {
                                                    Int_t trackC_num = PID_Array_C[Ana_Num][PID_C][ic];
                                                    if((trackC_num != trackA_num) && (trackC_num != trackB_num)) // Prevent that a track is used twice
                                                    {
                                                        StPicoAlexTrack trackC = *picoDst_C->track( trackC_num );
                                                        Float_t MomentumC   = trackC.gMom().mag();
                                                        Float_t dcaC        = trackC.dca();   // distance of closest approach to primary vertex
                                                        Float_t nHitsPossC  = trackC.nHitsMax();
                                                        Float_t nHitsFitC   = trackC.nHitsFit();
                                                        //Float_t PolarityC   = trackC.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                                        Float_t BetaC       = trackC.btofBeta();  // Velocity after time-of-flight reconstruction

                                                        Float_t Mass2C        = -100.0;
                                                        // calculate mass2
                                                        if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                                                        {
                                                            Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                                                        }

                                                        //Float_t MassC       = SquareRoot(Mass2C);
                                                        //Float_t TPCdEdxC    = trackC.dEdx(); // TPC dE/dx
                                                        Float_t nSigmaPionC = trackC.nSigmaPion();
                                                        //Float_t TofC        = trackC.btof();

                                                        if(
                                                           dcaC < 1.0
                                                           && nHitsFitC > 14
                                                           && MomentumC > 0.1
                                                           && MomentumC < 10.0
                                                           && (nHitsFitC/nHitsPossC) > 0.52
                                                           && nSigmaPionC > -2.0
                                                           && nSigmaPionC < 2.0
                                                           && ((BetaC > 0.0 && Mass2C < 0.05 && Mass2C > -0.05) || BetaC <= 0.0)
                                                          )
                                                        {
                                                            helixC = StPhysicalHelixD(trackC.gMom(),trackC.origin()+vectordiff,event_C_ana->bField()*MAGFIELDFACTOR,trackC.charge());
                                                            StThreeVectorF primdirC = helixC.cat(helixC.pathLength(vector_prim)); //
                                                            primdirC = MomentumC*primdirC/primdirC.mag();
                                                            ltrackC.SetXYZM(primdirC.x(),primdirC.y(),primdirC.z(),0.13957018);

                                                            // Missing mass and invariant mass calculations
                                                            TLorentzVector trackABC     = trackAB+ltrackC;
                                                            Double_t InvMassABC         = trackABC.M();

                                                            //Float_t pt2          = trackABC.Pt();  // Transverse momentum of mother particle
                                                            Float_t rap2         = trackABC.Rapidity(); // Rapidity of mother particle

                                                            Float_t pathC_prim,dcaC_prim;
                                                            fHelixAtoPointdca(vector_prim,helixC,pathC_prim,dcaC_prim);

                                                            StThreeVectorF vectornewC_prim;
                                                            vectornewC_prim     = helixC.cat(pathC_prim);
                                                            vectornewC_prim     = MomentumC*vectornewC_prim/vectornewC_prim.mag(); // momentum vector at primary vertex

                                                            delta_phi_ME = 0;
                                                            phi_event_plane                 = -400.0;
                                                            phi_event_plane_eta_gap         = -400.0;
                                                            delta_phi_ME_AB_weight          = 0.0;
                                                            delta_phi_ME_AB_weight_eta_gap  = 0.0;

                                                            calc_event_plane_angles(3,rap2,trackA,trackB,trackC,vectornewA_prim,vectornewB_prim,vectornewC_prim,ME_Flag,SE_ME_Flag,RunIdA,EventVertexXA,EventVertexYA,
                                                                                    EventVertexZA,RunIdB,EventVertexXB,EventVertexYB,EventVertexZB,
                                                                                    phi_event_plane,phi_event_plane_eta_gap,delta_phi_ME_AB_weight,delta_phi_ME_AB_weight_eta_gap);

                                                            // check whether the event planes between event A and B are close to each other
                                                            if(
                                                               fabs(delta_phi_ME_AB_weight) < TMath::DegToRad()*30.0
                                                               || (SE_ME_Flag == 0)
                                                              )
                                                            {
                                                                delta_phi_ME = 1;
                                                            }
                                                            //******************************************************************************************************

                                                            if(delta_phi_ME == 1)
                                                            {
                                                                h_inv_mass_Sigma1385_pt[0] ->Fill(inv_mass_sign*InvMassABC);
                                                                if(
                                                                   dcaA        > 0.3
                                                                   && dcaB     > 1.2
                                                                   && VerdistX > 4.0
                                                                   && VerdistY < 0.9
                                                                  )
                                                                {
                                                                    h_inv_mass_Sigma1385_pt[1] ->Fill(inv_mass_sign*InvMassABC);

                                                                    /*
                                                                    if(pt2 >= 0.4 && pt2 < 5.0)
                                                                    {
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][0]) h_inv_mass_Sigma1385_pt[16]->Fill(inv_mass_sign*InvMassABC); // 0-80%
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][0] && refMultA <= Ref_mult_table[fBeamTimeNum][3]) h_inv_mass_Sigma1385_pt[3]->Fill(inv_mass_sign*InvMassABC); // 50-80
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][3] && refMultA <= Ref_mult_table[fBeamTimeNum][6]) h_inv_mass_Sigma1385_pt[4]->Fill(inv_mass_sign*InvMassABC); // 20-50
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][6] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[5]->Fill(inv_mass_sign*InvMassABC); // 0-20
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][5] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[6]->Fill(inv_mass_sign*InvMassABC); // 0-30
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][7] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[7]->Fill(inv_mass_sign*InvMassABC); // 0-10
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][8] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[8]->Fill(inv_mass_sign*InvMassABC); // 0-5
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][6] && refMultA <= Ref_mult_table[fBeamTimeNum][7]) h_inv_mass_Sigma1385_pt[9]->Fill(inv_mass_sign*InvMassABC); // 10-20
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][5] && refMultA <= Ref_mult_table[fBeamTimeNum][6]) h_inv_mass_Sigma1385_pt[10]->Fill(inv_mass_sign*InvMassABC); // 20-30
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][4] && refMultA <= Ref_mult_table[fBeamTimeNum][5]) h_inv_mass_Sigma1385_pt[11]->Fill(inv_mass_sign*InvMassABC); // 30-40
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][3] && refMultA <= Ref_mult_table[fBeamTimeNum][4]) h_inv_mass_Sigma1385_pt[12]->Fill(inv_mass_sign*InvMassABC); // 40-50
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][2] && refMultA <= Ref_mult_table[fBeamTimeNum][3]) h_inv_mass_Sigma1385_pt[13]->Fill(inv_mass_sign*InvMassABC); // 50-60
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][1] && refMultA <= Ref_mult_table[fBeamTimeNum][2]) h_inv_mass_Sigma1385_pt[14]->Fill(inv_mass_sign*InvMassABC); // 60-70
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][0] && refMultA <= Ref_mult_table[fBeamTimeNum][1]) h_inv_mass_Sigma1385_pt[15]->Fill(inv_mass_sign*InvMassABC); // 70-80
                                                                    }
                                                                    if(pt2 >= 1.5 && pt2 < 5.0)
                                                                    {
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][0]) h_inv_mass_Sigma1385_pt[17]->Fill(inv_mass_sign*InvMassABC); // 0-80%
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][0] && refMultA <= Ref_mult_table[fBeamTimeNum][3]) h_inv_mass_Sigma1385_pt[18]->Fill(inv_mass_sign*InvMassABC); // 50-80
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][3] && refMultA <= Ref_mult_table[fBeamTimeNum][6]) h_inv_mass_Sigma1385_pt[19]->Fill(inv_mass_sign*InvMassABC); // 20-50
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][6] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[20]->Fill(inv_mass_sign*InvMassABC); // 0-20
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][5] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[21]->Fill(inv_mass_sign*InvMassABC); // 0-30
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][7] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[22]->Fill(inv_mass_sign*InvMassABC); // 0-10
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][8] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[23]->Fill(inv_mass_sign*InvMassABC); // 0-5
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][6] && refMultA <= Ref_mult_table[fBeamTimeNum][7]) h_inv_mass_Sigma1385_pt[24]->Fill(inv_mass_sign*InvMassABC); // 10-20
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][5] && refMultA <= Ref_mult_table[fBeamTimeNum][6]) h_inv_mass_Sigma1385_pt[25]->Fill(inv_mass_sign*InvMassABC); // 20-30
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][4] && refMultA <= Ref_mult_table[fBeamTimeNum][5]) h_inv_mass_Sigma1385_pt[26]->Fill(inv_mass_sign*InvMassABC); // 30-40
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][3] && refMultA <= Ref_mult_table[fBeamTimeNum][4]) h_inv_mass_Sigma1385_pt[27]->Fill(inv_mass_sign*InvMassABC); // 40-50
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][2] && refMultA <= Ref_mult_table[fBeamTimeNum][3]) h_inv_mass_Sigma1385_pt[28]->Fill(inv_mass_sign*InvMassABC); // 50-60
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][1] && refMultA <= Ref_mult_table[fBeamTimeNum][2]) h_inv_mass_Sigma1385_pt[29]->Fill(inv_mass_sign*InvMassABC); // 60-70
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][0] && refMultA <= Ref_mult_table[fBeamTimeNum][1]) h_inv_mass_Sigma1385_pt[30]->Fill(inv_mass_sign*InvMassABC); // 70-80
                                                                    }
                                                                    if(pt2 >= 0.4 && pt2 < 1.5)
                                                                    {
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][0]) h_inv_mass_Sigma1385_pt[31]->Fill(inv_mass_sign*InvMassABC); // 0-80%
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][0] && refMultA <= Ref_mult_table[fBeamTimeNum][3]) h_inv_mass_Sigma1385_pt[32]->Fill(inv_mass_sign*InvMassABC); // 50-80
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][3] && refMultA <= Ref_mult_table[fBeamTimeNum][6]) h_inv_mass_Sigma1385_pt[33]->Fill(inv_mass_sign*InvMassABC); // 20-50
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][6] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[34]->Fill(inv_mass_sign*InvMassABC); // 0-20
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][5] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[35]->Fill(inv_mass_sign*InvMassABC); // 0-30
                                                                        if(refMultA > Ref_mult_table[fhBeamTimeNum][7] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[36]->Fill(inv_mass_sign*InvMassABC); // 0-10
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][8] && refMultA <= Ref_mult_table[fBeamTimeNum][9]) h_inv_mass_Sigma1385_pt[37]->Fill(inv_mass_sign*InvMassABC); // 0-5
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][6] && refMultA <= Ref_mult_table[fBeamTimeNum][7]) h_inv_mass_Sigma1385_pt[38]->Fill(inv_mass_sign*InvMassABC); // 10-20
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][5] && refMultA <= Ref_mult_table[fBeamTimeNum][6]) h_inv_mass_Sigma1385_pt[39]->Fill(inv_mass_sign*InvMassABC); // 20-30
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][4] && refMultA <= Ref_mult_table[fBeamTimeNum][5]) h_inv_mass_Sigma1385_pt[40]->Fill(inv_mass_sign*InvMassABC); // 30-40
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][3] && refMultA <= Ref_mult_table[fBeamTimeNum][4]) h_inv_mass_Sigma1385_pt[41]->Fill(inv_mass_sign*InvMassABC); // 40-50
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][2] && refMultA <= Ref_mult_table[fBeamTimeNum][3]) h_inv_mass_Sigma1385_pt[42]->Fill(inv_mass_sign*InvMassABC); // 50-60
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][1] && refMultA <= Ref_mult_table[fBeamTimeNum][2]) h_inv_mass_Sigma1385_pt[43]->Fill(inv_mass_sign*InvMassABC); // 60-70
                                                                        if(refMultA > Ref_mult_table[fBeamTimeNum][0] && refMultA <= Ref_mult_table[fBeamTimeNum][1]) h_inv_mass_Sigma1385_pt[44]->Fill(inv_mass_sign*InvMassABC); // 70-80
                                                                    }
                                                                    */
                                                                }

                                                                if(
                                                                   dcaA        > 0.4
                                                                   && dcaB     > 1.5
                                                                   && VerdistX > 4.5
                                                                   && VerdistY < 0.7
                                                                  )
                                                                {
                                                                    h_inv_mass_Sigma1385_pt[2] ->Fill(inv_mass_sign*InvMassABC);
                                                                }

                                                            }

                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    //********************************************************************************************************************



                                }
                            }
                        }
                    }
                }
            }
        }
        if(fAnalysisNum == 102)
        {
            //cout << "Fill tree" << endl;
            Tree_Sigma0_v2  ->Fill();
        }

        //if(fAnalysisNum == 2 || fAnalysisNum == 14) Tree_V0_v2  ->Fill();
        return 1;
    }
    else return 0;
}



Int_t K0S_pippim_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs], Int_t PID_counter_Array_C[][N_max_PIDs],
                           Int_t PID_Array[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks], Int_t PID_Array_C[][N_max_PIDs][N_max_tracks],
                           StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B, StPicoAlexEvent* picoDst_C,
                           Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t Ana_Num,
                           Int_t SE_ME_Flag)
{
    //  mass = 497.648 MeV/c2
    //  ctau = 2.6842 cm
    //
    //  Au + Au -> K0S + N + N
    //              |
    //              -> pi+(A) + pi-(B)
    // Event vertex information
    //cout << "Lambda analysis started" << endl;
    StThreeVectorF vector_prim,vector_primB,vectordiff;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist,ref_mult_frac;
    Int_t refMultA,refMultB,RunIdA,RunIdB;

    event_A_ana   = picoDst_A;
    event_B_ana   = picoDst_B;
    event_C_ana   = picoDst_C;

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();
    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    refMultB          = refMultA;
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    //Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_A_ana->vzVpd();
    vector_prim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);
    Double_t pion_low_cut        = -0.5; //
    Double_t pion_high_cut       = 0.1;  // 0.05
    Double_t proton_low_cut      = -0.5;  // 0.76
    Double_t proton_high_cut     = 0.1;  // 1.12

    //cout << "refMultA = " << refMultA <<  endl;

    if(
       SE_ME_Flag == 1  // mixed event analysis for Lambdas
      )
    {
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        RunIdB         = event_B_ana->runId();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex

        ref_mult_frac = fabs(100*2.0*((Float_t)(refMultA - refMultB))/((Float_t)(refMultA + refMultB)));
    }

    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm
    //Float_t ME_ref_mult_frac_cut = 100.0; // 100% difference in reference multiplicity is allowed

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
       //&& ref_mult_frac < ME_ref_mult_frac_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
        //cout << "Triggers checked, Vertex checked..." << endl;
    }

    //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
    //    << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;

    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]   > 0   // pi+
       && PID_counter_Array_B[Ana_Num][ParticleB] > 0   // pi-
       //&& event->getNVertices() >= 1
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       //&& event_A_ana   ->isMinBias()
       //&& event_B_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
        //    << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB, helixC;
        StThreeVectorF vectorA, vectorB, vectorC, vectoratsA, vectoratsB, vectorAB, vector_primAB, vectornewA, vectornewB, vectornewA_lin, vectornewB_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY, dirY_lin;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackA_lin, ltrackB_lin, ltrackB2;
        TLorentzVector ltrackA_pip;

        /*
        alexV0_event.clearTrackList();
        alexV0_event.setx(EventVertexXA);
        alexV0_event.sety(EventVertexYA);
        alexV0_event.setz(EventVertexZA);
        alexV0_event.setid(RunIdA);
        alexV0_event.setmult(refMultA);
        alexV0_event.setn_prim(n_primaries);
        alexV0_event.setn_non_prim(n_non_primaries);
        alexV0_event.setn_tof_prim(n_tofmatch_prim);
        alexV0_event.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        alexV0_event.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        alexV0_event.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        alexV0_event.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        alexV0_event.setEP_Qx_ptw(EP_Qx_ptw);
        alexV0_event.setEP_Qy_ptw(EP_Qy_ptw);
        alexV0_event.setQtracks_eta_pos(Qtracks_used_eta_pos);
        alexV0_event.setQtracks_eta_neg(Qtracks_used_eta_neg);
        alexV0_event.setQtracks_full(Qtracks_used);
        alexV0_event.setZDCx(ZDCx);
        alexV0_event.setBBCx(BBCx);
        alexV0_event.setvzVpd(vzVpd);
        */

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // p candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t Mass2A        = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }

            Float_t MassA       = SquareRoot(Mass2A);
            //Float_t TPCdEdxA    = trackA.dEdx(); // TPC dE/dx
            Float_t nSigmaPA    = trackA.nSigmaPion();
            //Float_t TofA        = trackA.btof();
            if(nHitsPossA <= 0) nHitsPossA = 10000.0;

            Double_t dcaA_cut     = 0.5;
            Double_t dcaB_cut     = 0.5;
            Double_t VerdistX_cut = 2.5;
            Double_t VerdistY_cut = 1.5;
            Int_t    flag_proton  = 0;
            Int_t    flag_pion    = 0;
            Int_t    flag_tof_proton  = 0;
            Int_t    flag_tof_pion    = 0;

            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                flag_tof_proton = 1; // pi+
            }
            else
            {
                flag_tof_proton = 0;
            }

            if(flag_tof_proton == 1 && Mass2A < 0.1 && Mass2A > -0.5)
            {
                dcaA_cut     = 0.3;
                flag_proton  = 1;
            }
            else
            {
                flag_proton = 0;
            }

            //cout << "i = " << i << ", dcaA = " << dcaA << ", fitA = " << nHitsFitA << ", possA = " << nHitsPossA << ", pA =" << MomentumA << endl;

            if(
               dcaA         > dcaA_cut // 0.5
               && nHitsFitA > 14 // 14
               && MomentumA > 0.1
               && MomentumA < 10.0
               && (nHitsFitA/nHitsPossA) > 0.52
              )
            {

                for(Int_t j = 0; j < PID_counter_Array_B[Ana_Num][ParticleB]; j++)  // pi- candidates
                {
                    Int_t trackB_num = PID_Array_B[Ana_Num][ParticleB][j];
                    StPicoAlexTrack trackB    = *picoDst_B->track( trackB_num );
                    vectorB = trackB.origin();
                    vectorB = vectorB + vectordiff; // vectordiff == {0,0,0} for same event, don't do it for Sigma(1385) analysis
                    helixB = StPhysicalHelixD(trackB.gMom(),vectorB,event_B_ana->bField()*MAGFIELDFACTOR,trackB.charge());
                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    //Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction

                    Float_t Mass2B        = -100.0;
                    // calculate mass2
                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                    }


                    Float_t MassB       = SquareRoot(Mass2B);
                    //Float_t TPCdEdxB    = trackB.dEdx(); // TPC dE/dx
                    Float_t nSigmaPionB = trackB.nSigmaPion();
                    //Float_t TofB        = trackB.btof();
                    if(nHitsPossB <= 0) nHitsPossB = 10000.0;

                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        flag_tof_pion = 1;
                    }
                    else
                    {
                        flag_tof_pion = 0;
                    }

                    if(flag_tof_pion == 1 && Mass2B < 0.1 && Mass2B > -0.5) // time-of-flight is available and particle B is a pion
                    {
                        flag_pion = 1;
                    }
                    else
                    {
                        flag_pion = 0;
                    }
                    if(flag_proton == 0 && flag_pion == 0)
                    {
                        dcaA_cut     = 0.5;
                        dcaB_cut     = 0.5;
                        VerdistX_cut = 3.0;
                        VerdistY_cut = 0.8;
                    }
                    if((flag_proton == 0 && flag_pion == 1) || (flag_proton == 1 && flag_pion == 0))
                    {
                        dcaA_cut     = 0.4;
                        dcaB_cut     = 0.4;
                        VerdistX_cut = 3.0;
                        VerdistY_cut = 0.9;
                    }
                    if(flag_proton == 1 && flag_pion == 1)
                    {
                        dcaA_cut     = 0.3;
                        dcaB_cut     = 0.3;
                        VerdistX_cut = 2.0;
                        VerdistY_cut = 1.5;
                    }

                    //cout << "j = " << j << ", dcaB = " << dcaB << ", fitB = " << nHitsFitB << ", possB = " << nHitsPossB << ", pB =" << MomentumB << endl;

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       //&& dcaB > dcaA
                       && dcaA      > dcaA_cut // 1.0
                       && dcaB      > dcaB_cut // 1.0
                       && nHitsFitB > 14 // 14
                       && MomentumB > 0.1
                       && MomentumB < 10.0
                       && (nHitsFitB/nHitsPossB) > 0.52
                      )
                    {
                        Float_t pathA_f, pathB_f, dcaAB_f;
                        Float_t pathA_test = -1.0;
                        Float_t pathB_test = -1.0;
                        Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);

                        //cout << "pathA_test = " << pathA_test << ", pathB_test = " << pathB_test << endl;

                        vectoratsA     = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
                        vectoratsB     = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
                        vectorAB       = vectoratsA+vectoratsB;
                        vectorAB       = vectorAB/2.0; // decay vertex

                        //Double_t vex_est_x = vectorAB.x();
                        //Double_t vex_est_y = vectorAB.y();
                        //Double_t vex_est_z = vectorAB.z();

                        StThreeVectorF baseA,dirA,baseB,dirB;
                        baseA = helixA.at(pathA_test);
                        baseB = helixB.at(pathB_test);
                        dirA  = helixA.at(pathA_test-2.0) - helixA.at(pathA_test+2.0);
                        dirB  = helixB.at(pathB_test-2.0) - helixB.at(pathB_test+2.0);

                        StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

                        //Double_t vex_lin_x = vectorAB_lin.x();
                        //Double_t vex_lin_y = vectorAB_lin.y();
                        //Double_t vex_lin_z = vectorAB_lin.z();

                        Double_t dcaAB_lin = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
                        StThreeVectorF vectorABtoPrim_lin = vectorAB_lin - vector_prim; // vector primary vertex to decay vertex
                        Float_t VerdistX_lin = vectorABtoPrim_lin.mag(); // distance between primary vertex and decay verte


                        // calculate the scalar product with the approximated secondary vertex position
                        vectornewA_lin     = helixA.cat(pathA_test); // direction vector at dca for helixA
                        vectornewB_lin     = helixB.cat(pathB_test); // direction vector at dca for helixB
                        vectornewA_lin = MomentumA*vectornewA_lin/vectornewA_lin.mag(); // new momentum vector at decay vertex
                        vectornewB_lin = MomentumB*vectornewB_lin/vectornewB_lin.mag(); // new momentum vector at decay vertex
                        ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),0.13957018);
                        ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),0.13957018);
                        TLorentzVector trackAB_lin      = ltrackA_lin+ltrackB_lin; // mother particle
                        dirY_lin.set(trackAB_lin.Px(),trackAB_lin.Py(),trackAB_lin.Pz());
                        dirY_lin = dirY_lin/dirY_lin.mag();
                        Double_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim_lin/vectorABtoPrim_lin.mag());

                        //cout << "i = " << i << ", j = " << j << ", VerdistX_lin = " << VerdistX_lin << ", dcaAB_lin = " << dcaAB_lin << ", scalarProduct_lin = " << scalarProduct_lin << endl;

                        //if( (VerdistX_lin > 3.0 && dcaAB_lin < 1.5 && fDCA_Helix_out == 1) ||  fDCA_Helix_out == 0 )
                        if( VerdistX_lin > VerdistX_cut && dcaAB_lin < 1.5 && scalarProduct_lin > 0.0 )
                        {

                            if(fDCA_Helix_out == 1)
                            {
                                fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
                            }
                            else
                            {
                                fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
                            }

                            vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                            vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                            vectorAB       = vectoratsA+vectoratsB;
                            vectorAB       = vectorAB/2.0; // decay vertex

                            //Double_t vex_x = vectorAB.x();
                            //Double_t vex_y = vectorAB.y();
                            //Double_t vex_z = vectorAB.z();
                            //cout << "Vertex = {" << vex_x << ", " << vex_y << ", " << vex_z << "}" <<
                            //    ", Vertex_est = {" << vex_est_x << ", " << vex_est_y << ", " << vex_est_z << "}"<<
                            //    ", Vertex_lin = {" << vex_lin_x << ", " << vex_lin_y << ", " << vex_lin_z << "}"
                            //    << ", fDCA_Helix_out = " << fDCA_Helix_out << endl;

                            vectorABtoPrim = vectorAB - vector_prim; // vector primary vertex to decay vertex
                            Float_t VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

                            vectornewA     = helixA.cat(pathA_f); // direction vector at dca for helixA
                            vectornewB     = helixB.cat(pathB_f); // direction vector at dca for helixB

                            vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
                            vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex

                            ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.13957018);
                            ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),0.13957018);
                            ltrackA_pip.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.93827203);

                            //Float_t thetaA = vectoratsA.theta();
                            //Float_t thetaB = vectoratsB.theta();

                            // Missing mass and invariant mass calculations
                            TLorentzVector trackAB      = ltrackA+ltrackB; // mother particle
                            TLorentzVector trackAB_K0S  = ltrackA_pip+ltrackB; // mother particle
                            Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                            Double_t InvMassAB_K0S      = trackAB_K0S.M(); // invariant mass of mother particle
                            Float_t MomentumAB          = trackAB.P(); // momentum of mother particle
                            Float_t BetaAB = TMath::Sqrt(1.0/(1.0+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            dirY = dirY/dirY.mag();

                            Double_t scalarProduct = dirY.dot(vectorABtoPrim/vectorABtoPrim.mag());
                            //cout << "scalarProduct = " << scalarProduct << endl;

                            baseY = vectorAB;
                            Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vector_prim);

                            Float_t pt          = trackAB.Pt();  // Transverse momentum of mother particle
                            Float_t rap         = trackAB.Rapidity(); // Rapidity of mother particle

                            Float_t phiAB   = dirY.phi();
                            Float_t thetaAB = dirY.theta();

                            // beta correction Lambda ========================================================================================
                            Float_t BetaACorr  = BetaA;
                            Float_t Mass2ACorr = Mass2A;
                            Float_t BetaBCorr  = BetaB;
                            Float_t Mass2BCorr = Mass2B;
                            Float_t MassACorr  = MassA;
                            Float_t MassBCorr  = MassB;
                            Float_t PathAB = vectorABtoPrim.mag(); // Lambda uncharged !!


                            if((trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0) || (trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0))
                            {
                                if( debug_flag )
                                {
                                    cout << "------------------------------------------------- LAMBDA --------------------------------------------------------------------------" << endl;
                                    cout << "PathAB = " << PathAB << "\tInvMassAB = " << InvMassAB << "\tMomentumAB = " << MomentumAB
                                        << "\tBetaAB = " << BetaAB << "\tTofAB = " << PathAB / (BetaAB*29.9792458)<< endl;
                                }
                            }
                            if(flag_tof_proton == 1)
                            {
                                BetaACorr = correctBeta4SingleDecay(trackA,trackAB,helixA,vectorAB,PathAB);
                                Mass2ACorr = MomentumA * MomentumA * (1.0/(BetaACorr*BetaACorr)-1.0);
                                MassACorr = SquareRoot(Mass2ACorr);
                                if( debug_flag ) cout << "A) MassA = " << MassA << "\tMassACorr = " <<  MassACorr << endl;
                            }
                            else Mass2ACorr = -100;

                            if(flag_tof_pion == 1)
                            {
                                BetaBCorr = correctBeta4SingleDecay(trackB,trackAB,helixB,vectorAB,PathAB);
                                Mass2BCorr = MomentumB * MomentumB * (1.0/(BetaBCorr*BetaBCorr)-1.0);
                                MassBCorr = SquareRoot(Mass2BCorr);
                                if( debug_flag ) cout << "B) MassB = " << MassB << "\tMassBCorr = " <<  MassBCorr << endl;
                            }
                            else Mass2BCorr = -100;


                            if(
                               InvMassAB        > 0.4
                               && InvMassAB     < 0.6
                               && VerdistX      > VerdistX_cut  // 3.5
                               && dcaAB_f       < 1.5  // 1.5
                               && VerdistY      < VerdistY_cut  // 1.5
                               && scalarProduct > 0.0
                               && ((flag_tof_proton == 0) || (Mass2ACorr > 100.0) || (Mass2ACorr  > proton_low_cut && Mass2ACorr < proton_high_cut))
                               && ((flag_tof_pion   == 0) || (Mass2BCorr > 100.0) || (Mass2BCorr  > pion_low_cut && Mass2BCorr   < pion_high_cut))
                              )
                            {
                                //cout << "Mass checked, topology checked --> event accepted" << endl;
                                //****************** Event plane calculations ***********************************************************
                                // Calculate the same vectors as for the event plane -> no V0!
                                Float_t pathA_prim,dcaA_prim,pathB_prim,dcaB_prim;
                                fHelixAtoPointdca(vector_prim,helixA,pathA_prim,dcaA_prim);
                                fHelixAtoPointdca(vector_prim,helixB,pathB_prim,dcaB_prim);

                                StThreeVectorF vectornewA_prim,vectornewB_prim;
                                vectornewA_prim     = helixA.cat(pathA_prim);
                                vectornewB_prim     = helixB.cat(pathB_prim);

                                vectornewA_prim     = MomentumA*vectornewA_prim/vectornewA_prim.mag(); // momentum vector at primary vertex
                                vectornewB_prim     = MomentumB*vectornewB_prim/vectornewB_prim.mag(); // momentum vector at primary vertex

                                Int_t delta_phi_ME = 0;
                                Float_t phi_event_plane                 = -400.0;
                                Float_t phi_event_plane_eta_gap         = -400.0;
                                Float_t delta_phi_ME_AB_weight          = 0.0;
                                Float_t delta_phi_ME_AB_weight_eta_gap  = 0.0;

                                if(fAnalysisNum == 3)
                                {
                                    // Calculate the event plane anlges
                                    calc_event_plane_angles(2,rap,trackA,trackB,trackA,vectornewA_prim,vectornewB_prim,vectornewA_prim,ME_Flag,SE_ME_Flag,RunIdA,EventVertexXA,EventVertexYA,
                                                            EventVertexZA,RunIdB,EventVertexXB,EventVertexYB,EventVertexZB,
                                                            phi_event_plane,phi_event_plane_eta_gap,delta_phi_ME_AB_weight,delta_phi_ME_AB_weight_eta_gap);
                                }
                                else delta_phi_ME = 1;

                                // check whether the event planes between event A and B are close to each other
                                if(
                                   fabs(delta_phi_ME_AB_weight) < TMath::DegToRad()*30.0
                                   || (SE_ME_Flag == 0)
                                  )
                                {
                                    delta_phi_ME = 1;
                                }
                                //******************************************************************************************************


                                if(delta_phi_ME == 1)
                                {
                                    Float_t p_xA_c   = vectornewA_prim.x();
                                    Float_t p_yA_c   = vectornewA_prim.y();
                                    Float_t p_tA_c   = sqrt(p_xA_c*p_xA_c + p_yA_c*p_yA_c);
                                    Float_t etaA_c   = vectornewA_prim.pseudoRapidity();
                                    Float_t phiA_c   = vectornewA_prim.phi();
                                    Double_t p_t_weightA = 1.0;
                                    if(p_tA_c < 2.0)  p_t_weightA = p_tA_c;
                                    if(p_tA_c >= 2.0) p_t_weightA = 2.0;
                                    Float_t iQxA     = p_t_weightA*TMath::Cos(2.0*phiA_c);
                                    Float_t iQyA     = p_t_weightA*TMath::Sin(2.0*phiA_c);

                                    Float_t p_xB_c   = vectornewB_prim.x();
                                    Float_t p_yB_c   = vectornewB_prim.y();
                                    Float_t p_tB_c   = sqrt(p_xB_c*p_xB_c + p_yB_c*p_yB_c);
                                    Float_t etaB_c   = vectornewB_prim.pseudoRapidity();
                                    Float_t phiB_c   = vectornewB_prim.phi();
                                    Double_t p_t_weightB = 1.0;
                                    if(p_tB_c < 2.0)  p_t_weightB = p_tB_c;
                                    if(p_tB_c >= 2.0) p_t_weightB = 2.0;
                                    Float_t iQxB     = p_t_weightB*TMath::Cos(2.0*phiB_c);
                                    Float_t iQyB     = p_t_weightB*TMath::Sin(2.0*phiB_c);

                                    /*
                                    alexV0_track = alexV0_event.createTrack();
                                    alexV0_track->setm2A(Mass2ACorr);
                                    alexV0_track->setm2B(Mass2BCorr);
                                    alexV0_track->setm2C(-1.0);
                                    alexV0_track->setnsA(nSigmaPA);
                                    alexV0_track->setnsB(nSigmaPionB);
                                    alexV0_track->setnsC(-1.0);
                                    alexV0_track->setdcaA(dcaA);
                                    alexV0_track->setdcaB(dcaB);
                                    alexV0_track->setdcaC(-1.0);
                                    alexV0_track->setiQxA(iQxA);
                                    alexV0_track->setiQyA(iQyA);
                                    alexV0_track->setiQxB(iQxB);
                                    alexV0_track->setiQyB(iQyB);
                                    alexV0_track->setiQxC(-1.0);
                                    alexV0_track->setiQyC(-1.0);
                                    alexV0_track->setetaA(etaA_c);
                                    alexV0_track->setetaB(etaB_c);
                                    alexV0_track->setetaC(-1.0);
                                    alexV0_track->setInvAB(InvMassAB);
                                    alexV0_track->setInvABC(-1.0);
                                    alexV0_track->setInvAB_miss(InvMassAB_K0S);
                                    alexV0_track->setInvABC_miss(-1.0);
                                    alexV0_track->setdcaAB(dcaAB_f);
                                    alexV0_track->setdcaBC(-1.0);
                                    alexV0_track->setdcaABC(-1.0);
                                    alexV0_track->setVerdistX(VerdistX);
                                    alexV0_track->setVerdistY(VerdistY);
                                    alexV0_track->setVerdistX2(-1.0);
                                    alexV0_track->setVerdistY2(-1.0);
                                    alexV0_track->setpt(pt);
                                    alexV0_track->setrap(rap);
                                    alexV0_track->setphi(phiAB);
                                    alexV0_track->settheta(thetaAB);
                                    alexV0_track->setPsi_ep(phi_event_plane);
                                    alexV0_track->setPsi_ep_eta(phi_event_plane_eta_gap);
                                    alexV0_track->setPsi_diff_ME(delta_phi_ME_AB_weight);
                                    alexV0_track->setscal_prod(scalarProduct);
                                    alexV0_track->setscal_prod2(-1.0);
                                    */


                                    // Filling Ntuple
                                    Lambda_X_NTDataArray[0]      =(Float_t)InvMassAB;
                                    Lambda_X_NTDataArray[1]      =(Float_t)InvMassAB_K0S;

                                    Lambda_X_NTDataArray[2]      =(Float_t)Mass2ACorr;
                                    Lambda_X_NTDataArray[3]      =(Float_t)Mass2BCorr;
                                    Lambda_X_NTDataArray[4]      =(Float_t)nSigmaPA;
                                    Lambda_X_NTDataArray[5]      =(Float_t)nSigmaPionB;
                                    Lambda_X_NTDataArray[6]      =(Float_t)MomentumA;
                                    Lambda_X_NTDataArray[7]      =(Float_t)MomentumB;
                                    Lambda_X_NTDataArray[8]      =(Float_t)dcaA;
                                    Lambda_X_NTDataArray[9]      =(Float_t)dcaB;
                                    Lambda_X_NTDataArray[10]     =(Float_t)iQxA;
                                    Lambda_X_NTDataArray[11]     =(Float_t)iQyA;
                                    Lambda_X_NTDataArray[12]     =(Float_t)iQxB;
                                    Lambda_X_NTDataArray[13]     =(Float_t)iQyB;
                                    Lambda_X_NTDataArray[14]     =(Float_t)etaA_c;
                                    Lambda_X_NTDataArray[15]     =(Float_t)etaB_c;

                                    Lambda_X_NTDataArray[16]     =(Float_t)refMultA; //1
                                    Lambda_X_NTDataArray[17]     =(Float_t)dcaAB_f;
                                    Lambda_X_NTDataArray[18]     =(Float_t)VerdistX;
                                    Lambda_X_NTDataArray[19]     =(Float_t)VerdistY;
                                    Lambda_X_NTDataArray[20]     =(Float_t)pt;
                                    Lambda_X_NTDataArray[21]     =(Float_t)rap;
                                    Lambda_X_NTDataArray[22]     =(Float_t)phiAB;
                                    Lambda_X_NTDataArray[23]     =(Float_t)thetaAB;
                                    Lambda_X_NTDataArray[24]     =(Float_t)phi_event_plane;
                                    Lambda_X_NTDataArray[25]     =(Float_t)phi_event_plane_eta_gap;
                                    Lambda_X_NTDataArray[26]     =(Float_t)delta_phi_ME_AB_weight;
                                    Lambda_X_NTDataArray[27]     =(Float_t)RunIdA;  //1
                                    Lambda_X_NTDataArray[28]     =(Float_t)n_primaries; //1
                                    Lambda_X_NTDataArray[29]     =(Float_t)n_non_primaries; //1
                                    Lambda_X_NTDataArray[30]     =(Float_t)n_tofmatch_prim; //1
                                    Lambda_X_NTDataArray[31]     =(Float_t)scalarProduct;
                                    Lambda_X_NTDataArray[32]     =(Float_t)EventVertexXA; //1
                                    Lambda_X_NTDataArray[33]     =(Float_t)EventVertexYA; //1
                                    Lambda_X_NTDataArray[34]     =(Float_t)EventVertexZA; //1
                                    Lambda_X_NTDataArray[35]     =(Float_t)EP_Qx_eta_pos_ptw;  //1
                                    Lambda_X_NTDataArray[36]     =(Float_t)EP_Qy_eta_pos_ptw;  //1
                                    Lambda_X_NTDataArray[37]     =(Float_t)EP_Qx_eta_neg_ptw;  //1
                                    Lambda_X_NTDataArray[38]     =(Float_t)EP_Qy_eta_neg_ptw;  //1
                                    Lambda_X_NTDataArray[39]     =(Float_t)EP_Qx_ptw; //1
                                    Lambda_X_NTDataArray[40]     =(Float_t)EP_Qy_ptw; //1
                                    Lambda_X_NTDataArray[41]     =(Float_t)Qtracks_used_eta_pos; //1
                                    Lambda_X_NTDataArray[42]     =(Float_t)Qtracks_used_eta_neg; //1
                                    Lambda_X_NTDataArray[43]     =(Float_t)Qtracks_used; //1
                                    Lambda_X_NTDataArray[44]     =(Float_t)ZDCx; //1
                                    Lambda_X_NTDataArray[45]     =(Float_t)vzVpd; //1

                                    if(fAnalysisNum == 3)
                                    {
                                        Lambda_X_NT->Fill(Lambda_X_NTDataArray);
                                    }

                                    if(ME_Flag == 1)
                                    {
                                        Lambda_ME_counter++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        //if(fAnalysisNum == 2 || fAnalysisNum == 14) Tree_V0_v2  ->Fill();
        return 1;
    }
    else return 0;
}




Int_t D0_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs], Int_t PID_counter_Array_C[][N_max_PIDs],
                           Int_t PID_Array[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks], Int_t PID_Array_C[][N_max_PIDs][N_max_tracks],
                           StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B, StPicoAlexEvent* picoDst_C,
                           Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t Ana_Num,
                           Int_t SE_ME_Flag)
{
    //  mass = 1864.5 MeV/c2
    //  ctau = 122.9 mum
    //
    //  Au + Au -> D0 + N + N
    //              |
    //              -> K-(A) + pi+(B)
    // Event vertex information
    StThreeVectorF vector_prim,vector_primB,vectordiff;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist,ref_mult_frac;
    Int_t refMultA,refMultB,RunIdA,RunIdB;

    event_A_ana   = picoDst_A;
    event_B_ana   = picoDst_B;
    event_C_ana   = picoDst_C;

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();
    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    refMultB          = refMultA;
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_A_ana->vzVpd();
    vector_prim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    //cout << "refMultA = " << refMultA <<  endl;

    //cout << "D0 analysis" << endl;

    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        RunIdB         = event_B_ana->runId();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex

        ref_mult_frac = fabs(100*2.0*((Float_t)(refMultA - refMultB))/((Float_t)(refMultA + refMultB)));
    }

    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;
    Float_t ME_vertex_dist_cut   = 1.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
       //&& ref_mult_frac < ME_ref_mult_frac_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
        //cout << "Triggers checked, Vertex checked..." << endl;
    }

    //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
    //    << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;

    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]   > 0   // K-
       && PID_counter_Array_B[Ana_Num][ParticleB] > 0   // pi+
       //&& event->getNVertices() >= 1
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       //&& event_A_ana   ->isMinBias()
       //&& event_B_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
        //    << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB, helixC;
        StThreeVectorF vectorA, vectorB, vectorC, vectoratsA, vectoratsB, vectorAB, vector_primAB, vectornewA, vectornewB, vectornewA_lin, vectornewB_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY, dirY_lin;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackA_lin, ltrackB_lin, ltrackB2;
        TLorentzVector ltrackA_pip;


        D0_event.clearTrackList();
        D0_event.setx(EventVertexXA);
        D0_event.sety(EventVertexYA);
        D0_event.setz(EventVertexZA);
        D0_event.setid(RunIdA);
        D0_event.setmult(refMultA);
        D0_event.setn_prim(n_primaries);
        D0_event.setn_non_prim(n_non_primaries);
        D0_event.setn_tof_prim(n_tofmatch_prim);
        D0_event.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        D0_event.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        D0_event.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        D0_event.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        D0_event.setEP_Qx_ptw(EP_Qx_ptw);
        D0_event.setEP_Qy_ptw(EP_Qy_ptw);
        D0_event.setQtracks_eta_pos(Qtracks_used_eta_pos);
        D0_event.setQtracks_eta_neg(Qtracks_used_eta_neg);
        D0_event.setQtracks_full(Qtracks_used);
        D0_event.setZDCx(ZDCx);
        D0_event.setBBCx(BBCx);
        D0_event.setvzVpd(vzVpd);


        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // K- candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );

            Int_t HFTmapA = trackA.nHitsMapHFT();
            Int_t flag_HFTA = 0;
            if((HFTmapA>>0 & 0x1) && (HFTmapA>>1 & 0x3) && (HFTmapA>>3 & 0x3)) flag_HFTA = 1;
            if(!flag_HFTA) continue;

            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t Mass2A      = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }

            Float_t nSigmaKaonA    = trackA.nSigmaKaon();
            if(fabs(nSigmaKaonA) > 2.5) continue;
            if(Mass2A > -90.0) // time-of-flight info exists
            {
                if(!(Mass2A > 0.14 && Mass2A < 0.4)) continue; // not in kaon mass range
            }

            Double_t dcaA_cut         = 0.005;
            Double_t dcaB_cut         = 0.005;
            Double_t VerdistX_cut     = 0.008;
            Double_t VerdistY_cut     = 0.006;
            Int_t    flag_kaon        = 0;
            Int_t    flag_pion        = 0;
            Int_t    flag_tof_kaon    = 0;
            Int_t    flag_tof_pion    = 0;

            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                flag_tof_kaon = 1;
            }
            else
            {
                flag_tof_kaon = 0;
            }

            //cout << "i = " << i << ", dcaA = " << dcaA << ", fitA = " << nHitsFitA << ", possA = " << nHitsPossA << ", pA =" << MomentumA << endl;

            if(
               fabs(dcaA)    > dcaA_cut
               && fabs(dcaA) < 0.04
               && nHitsFitA  > 14
               && MomentumA  > 0.1
               && MomentumA  < 10.0
               //&& nHitsPossA > 0
               //&& (nHitsFitA/nHitsPossA) > 0.52
              )
            {
                //cout << "Track A accepted" << endl;
                for(Int_t j = 0; j < PID_counter_Array_B[Ana_Num][ParticleB]; j++)  // pi+ candidates
                {
                    Int_t trackB_num = PID_Array_B[Ana_Num][ParticleB][j];
                    StPicoAlexTrack trackB    = *picoDst_B->track( trackB_num );

                    Int_t HFTmapB = trackB.nHitsMapHFT();
                    Int_t flag_HFTB = 0;
                    if((HFTmapB>>0 & 0x1) && (HFTmapB>>1 & 0x3) && (HFTmapB>>3 & 0x3)) flag_HFTB = 1;
                    if(!flag_HFTB) continue;

                    vectorB = trackB.origin();
                    vectorB = vectorB + vectordiff; // vectordiff == {0,0,0} for same event
                    helixB = StPhysicalHelixD(trackB.gMom(),vectorB,event_B_ana->bField()*MAGFIELDFACTOR,trackB.charge());
                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                    Float_t Mass2B      = -100.0;
                    // calculate mass2
                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                    }

                    Float_t nSigmaPionB = trackB.nSigmaPion();
                    if(fabs(nSigmaPionB) > 2.5) continue;
                    if(Mass2B > -90.0) // time-of-fligh info exists
                    {
                        if(!(Mass2B > -0.4 && Mass2B < 0.14)) continue; // not in pion mass range
                    }

                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        flag_tof_pion = 1;
                    }
                    else
                    {
                        flag_tof_pion = 0;
                    }

                    //cout << "j = " << j << ", dcaB = " << dcaB << ", fitB = " << nHitsFitB << ", possB = " << nHitsPossB << ", pB =" << MomentumB << endl;

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       //&& dcaB > dcaA
                       && fabs(dcaA) > dcaA_cut
                       && fabs(dcaB) > dcaB_cut
                       && fabs(dcaA) < 0.04
                       && fabs(dcaB) < 0.04
                       && nHitsFitB  > 14
                       && MomentumB  > 0.1
                       && MomentumB  < 10.0
                       //&& (nHitsFitB/nHitsPossB) > 0.52
                      )
                    {
                        //cout << "Track B accepted" << endl;
                        Float_t pathA_f, pathB_f, dcaAB_f;
                        Float_t pathA_test = -1.0;
                        Float_t pathB_test = -1.0;
                        Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);

                        //cout << "fDCA_Helix_out: " << fDCA_Helix_out << ", pathA_test: " << pathA_test << ", pathB_test: " << pathB_test << ", dcaAB_f: " << dcaAB_f << endl;

                        vectoratsA     = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
                        vectoratsB     = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
                        vectorAB       = vectoratsA+vectoratsB;
                        vectorAB       = vectorAB/2.0; // decay vertex

                        //Double_t vex_est_x = vectorAB.x();
                        //Double_t vex_est_y = vectorAB.y();
                        //Double_t vex_est_z = vectorAB.z();

                        StThreeVectorF baseA,dirA,baseB,dirB;
                        baseA = helixA.at(pathA_test);
                        baseB = helixB.at(pathB_test);
                        dirA  = helixA.at(pathA_test-0.01) - helixA.at(pathA_test+0.01);
                        dirB  = helixB.at(pathB_test-0.01) - helixB.at(pathB_test+0.01);

                        StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

                        //Double_t vex_lin_x = vectorAB_lin.x();
                        //Double_t vex_lin_y = vectorAB_lin.y();
                        //Double_t vex_lin_z = vectorAB_lin.z();

                        Double_t dcaAB_lin = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
                        StThreeVectorF vectorABtoPrim_lin = vectorAB_lin - vector_prim; // vector primary vertex to decay vertex
                        Float_t VerdistX_lin = vectorABtoPrim_lin.mag(); // distance between primary vertex and decay verte


                        // calculate the scalar product with the approximated secondary vertex position
                        vectornewA_lin     = helixA.cat(pathA_test); // direction vector at dca for helixA
                        vectornewB_lin     = helixB.cat(pathB_test); // direction vector at dca for helixB
                        vectornewA_lin = MomentumA*vectornewA_lin/vectornewA_lin.mag(); // new momentum vector at decay vertex
                        vectornewB_lin = MomentumB*vectornewB_lin/vectornewB_lin.mag(); // new momentum vector at decay vertex
                        ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),mass_array[ParticleA]);
                        ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),mass_array[ParticleB]);
                        TLorentzVector trackAB_lin      = ltrackA_lin+ltrackB_lin; // mother particle
                        dirY_lin.set(trackAB_lin.Px(),trackAB_lin.Py(),trackAB_lin.Pz());
                        dirY_lin = dirY_lin/dirY_lin.mag();
                        Double_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim_lin/vectorABtoPrim_lin.mag());

                        //cout << "i = " << i << ", j = " << j << ", VerdistX_lin = " << VerdistX_lin << ", dcaAB_lin = " << dcaAB_lin << ", scalarProduct_lin = " << scalarProduct_lin << endl;

                        //if( (VerdistX_lin > 3.0 && dcaAB_lin < 1.5 && fDCA_Helix_out == 1) ||  fDCA_Helix_out == 0 )
                        if( VerdistX_lin > VerdistX_cut && dcaAB_lin < 0.2 && scalarProduct_lin > 0.0 )
                        {
                            Float_t VerdistX;
                            if(fDCA_Helix_out == 1)
                            {
                                VerdistX       = VerdistX_lin;
                                vectornewA     = vectornewA_lin;
                                vectornewB     = vectornewB_lin;
                                dcaAB_f        = dcaAB_lin;
                                vectorAB       = vectorAB_lin;
                                vectorABtoPrim = vectorABtoPrim_lin;
                                vectornewA     = vectornewA_lin;
                                vectornewB     = vectornewB_lin;

                                //fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
                            }
                            else
                            {
                                fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices

                                vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                                vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                                vectorAB       = vectoratsA+vectoratsB;
                                vectorAB       = vectorAB/2.0; // decay vertex

                                vectorABtoPrim = vectorAB - vector_prim; // vector primary vertex to decay vertex
                                VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

                                vectornewA     = helixA.cat(pathA_f); // direction vector at dca for helixA
                                vectornewB     = helixB.cat(pathB_f); // direction vector at dca for helixB

                                vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
                                vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex
                            }


                            Float_t etaA_c   = vectornewA.pseudoRapidity();
                            Float_t etaB_c   = vectornewB.pseudoRapidity();

                            ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),mass_array[ParticleA]);
                            ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),mass_array[ParticleB]);

                            // Missing mass and invariant mass calculations
                            TLorentzVector trackAB      = ltrackA+ltrackB; // mother particle
                            Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                            Float_t MomentumAB          = trackAB.P(); // momentum of mother particle
                            Float_t BetaAB = TMath::Sqrt(1.0/(1.0+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            dirY = dirY/dirY.mag();

                            Double_t scalarProduct = dirY.dot(vectorABtoPrim/vectorABtoPrim.mag());
                            //cout << "scalarProduct = " << scalarProduct << endl;

                            baseY = vectorAB;
                            Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vector_prim);

                            //cout << "VerdistX: " << VerdistX  << ", VerdistY: " << VerdistY << ", dcaAB_f: " << dcaAB_f <<  endl;

                            h_dca_diff          ->Fill(dcaAB_f-dcaAB_lin);
                            h_decay_length_diff ->Fill(VerdistX-VerdistX_lin);

                            Float_t pt          = trackAB.Pt();  // Transverse momentum of mother particle
                            Float_t rap         = trackAB.Rapidity(); // Rapidity of mother particle

                            Float_t phiAB   = dirY.phi();
                            Float_t thetaAB = dirY.theta();

                            if(
                               VerdistX         > VerdistX_cut
                               && VerdistX      < 0.5
                               && dcaAB_f       < 0.009
                               && VerdistY      < VerdistY_cut
                               && scalarProduct > 0.0
                               && InvMassAB     > 1.6
                               && InvMassAB     < 2.1
                              )
                            {
                                //cout << "Mass checked, topology checked --> combination accepted" << endl;
                                //****************** Event plane calculations ***********************************************************
#if 0
                                // Calculate the same vectors as for the event plane -> no V0!
                                Float_t pathA_prim,dcaA_prim,pathB_prim,dcaB_prim;
                                fHelixAtoPointdca(vector_prim,helixA,pathA_prim,dcaA_prim);
                                fHelixAtoPointdca(vector_prim,helixB,pathB_prim,dcaB_prim);

                                StThreeVectorF vectornewA_prim,vectornewB_prim;
                                vectornewA_prim     = helixA.cat(pathA_prim);
                                vectornewB_prim     = helixB.cat(pathB_prim);

                                vectornewA_prim     = MomentumA*vectornewA_prim/vectornewA_prim.mag(); // momentum vector at primary vertex
                                vectornewB_prim     = MomentumB*vectornewB_prim/vectornewB_prim.mag(); // momentum vector at primary vertex

                                Int_t delta_phi_ME = 0;
                                Float_t phi_event_plane                 = -400.0;
                                Float_t phi_event_plane_eta_gap         = -400.0;
                                Float_t delta_phi_ME_AB_weight          = 0.0;
                                Float_t delta_phi_ME_AB_weight_eta_gap  = 0.0;

                                // Calculate the event plane anlges
                                calc_event_plane_angles(2,rap,trackA,trackB,trackA,vectornewA_prim,vectornewB_prim,vectornewA_prim,ME_Flag,SE_ME_Flag,RunIdA,EventVertexXA,EventVertexYA,
                                                        EventVertexZA,RunIdB,EventVertexXB,EventVertexYB,EventVertexZB,
                                                        phi_event_plane,phi_event_plane_eta_gap,delta_phi_ME_AB_weight,delta_phi_ME_AB_weight_eta_gap);

                                // check whether the event planes between event A and B are close to each other
                                if(
                                   fabs(delta_phi_ME_AB_weight) < TMath::DegToRad()*30.0
                                   || (SE_ME_Flag == 0)
                                  )
                                {
                                    delta_phi_ME = 1;
                                }
#endif
                                //******************************************************************************************************


                                //if(delta_phi_ME == 1)
                                if(1)
                                {
#if 0
                                    Float_t p_xA_c   = vectornewA_prim.x();
                                    Float_t p_yA_c   = vectornewA_prim.y();
                                    Float_t p_tA_c   = sqrt(p_xA_c*p_xA_c + p_yA_c*p_yA_c);
                                    Float_t etaA_c   = vectornewA_prim.pseudoRapidity();
                                    Float_t phiA_c   = vectornewA_prim.phi();
                                    Double_t p_t_weightA = 1.0;
                                    if(p_tA_c < 2.0)  p_t_weightA = p_tA_c;
                                    if(p_tA_c >= 2.0) p_t_weightA = 2.0;
                                    Float_t iQxA     = p_t_weightA*TMath::Cos(2.0*phiA_c);
                                    Float_t iQyA     = p_t_weightA*TMath::Sin(2.0*phiA_c);

                                    Float_t p_xB_c   = vectornewB_prim.x();
                                    Float_t p_yB_c   = vectornewB_prim.y();
                                    Float_t p_tB_c   = sqrt(p_xB_c*p_xB_c + p_yB_c*p_yB_c);
                                    Float_t etaB_c   = vectornewB_prim.pseudoRapidity();
                                    Float_t phiB_c   = vectornewB_prim.phi();
                                    Double_t p_t_weightB = 1.0;
                                    if(p_tB_c < 2.0)  p_t_weightB = p_tB_c;
                                    if(p_tB_c >= 2.0) p_t_weightB = 2.0;
                                    Float_t iQxB     = p_t_weightB*TMath::Cos(2.0*phiB_c);
                                    Float_t iQyB     = p_t_weightB*TMath::Sin(2.0*phiB_c);
#endif

                                    Float_t iQxA   = 0.0;
                                    Float_t iQyA   = 0.0;
                                    Float_t iQxB   = 0.0;
                                    Float_t iQyB   = 0.0;
                                    //Float_t etaA_c = 0.0;
                                    //Float_t etaB_c = 0.0;

                                    // Filling Ntuple
                                    D0_track = D0_event.createTrack();
                                    D0_track->setm2A(Mass2A);
                                    D0_track->setm2B(Mass2B);
                                    D0_track->setnsA(nSigmaKaonA);
                                    D0_track->setnsB(nSigmaPionB);
                                    D0_track->setdcaA(dcaA);
                                    D0_track->setdcaB(dcaB);
                                    D0_track->setiQxA(iQxA);
                                    D0_track->setiQyA(iQyA);
                                    D0_track->setiQxB(iQxB);
                                    D0_track->setiQyB(iQyB);
                                    D0_track->setetaA(etaA_c);
                                    D0_track->setetaB(etaB_c);
                                    D0_track->setInvAB(InvMassAB);
                                    D0_track->setpt(pt);
                                    D0_track->setrap(rap);
                                    D0_track->setphi(phiAB);
                                    D0_track->settheta(thetaAB);
                                    D0_track->setqpA(MomentumA*PolarityA);
                                    D0_track->setqpB(MomentumB*PolarityB);
                                    D0_track->setVerdistX(VerdistX);
                                    D0_track->setVerdistY(VerdistY);
                                    D0_track->setdcaAB(dcaAB_f);

                                }
                            }
                        }
                    }
                }
            }
        }
        Tree_D0_v2  ->Fill();
        return 1;
    }
    else return 0;
}




Int_t ThetaPlus_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs], Int_t PID_counter_Array_C[][N_max_PIDs],
                           Int_t PID_Array[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks], Int_t PID_Array_C[][N_max_PIDs][N_max_tracks],
                           StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B, StPicoAlexEvent* picoDst_C,
                           Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t Ana_Num,
                           Int_t SE_ME_Flag)
{
    //  mass = unknown MeV/c2
    //  ctau = unknown cm
    //
    //  Au + Au -> K0S + p + X
    //              |    -> (C)
    //              -> pi+(A) + pi-(B)
    // Event vertex information
    //cout << "Lambda analysis started" << endl;
    StThreeVectorF vector_prim,vector_primB,vectordiff;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist,ref_mult_frac;
    Int_t refMultA,refMultB,RunIdA,RunIdB;

    event_A_ana   = picoDst_A; // pi+
    event_B_ana   = picoDst_B; // pi-
    event_C_ana   = picoDst_C; // proton

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();
    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    refMultB          = refMultA;
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    //Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    //Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    //Float_t vzVpd     = event_A_ana->vzVpd();
    vector_prim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);
    Double_t pion_low_cut        = -0.5; //
    Double_t pion_high_cut       = 0.1;  // 0.05
    Double_t proton_low_cut      = -0.5;  // 0.76
    Double_t proton_high_cut     = 0.1;  // 1.12

    //cout << "refMultA = " << refMultA <<  endl;

    if(
       SE_ME_Flag == 1  // mixed event analysis for Lambdas
       && (fAnalysisNum == 3)
      )
    {
        EventVertexXB  = event_C_ana->primaryVertex().x();
        EventVertexYB  = event_C_ana->primaryVertex().y();
        EventVertexZB  = event_C_ana->primaryVertex().z();
        refMultB       = event_C_ana->refMult();
        RunIdB         = event_C_ana->runId();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex

        ref_mult_frac = fabs(100*2.0*((Float_t)(refMultA - refMultB))/((Float_t)(refMultA + refMultB)));
    }

    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm
    //Float_t ME_ref_mult_frac_cut = 100.0; // 100% difference in reference multiplicity is allowed

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
       //&& ref_mult_frac < ME_ref_mult_frac_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
        //cout << "Triggers checked, Vertex checked..." << endl;
    }

    //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
    //    << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;

    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]   > 0   // pi+
       && PID_counter_Array_B[Ana_Num][ParticleB] > 0   // pi-
       && PID_counter_Array_C[Ana_Num][ParticleC] > 0   // proton
       //&& event->getNVertices() >= 1
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana   ->isMinBias()
       && event_B_ana   ->isMinBias()
       && event_C_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
        //    << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB, helixC;
        StThreeVectorF vectorA, vectorB, vectorC, vectoratsA, vectoratsB, vectorAB, vector_primAB, vectornewA, vectornewB, vectornewA_lin, vectornewB_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY, dirY_lin;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackA_lin, ltrackB_lin, ltrackB2;
        TLorentzVector ltrackA_pip;

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // pi+ candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction

            Float_t Mass2A        = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }

            Float_t MassA       = SquareRoot(Mass2A);
            //Float_t TPCdEdxA    = trackA.dEdx(); // TPC dE/dx
            Float_t nSigmaPA    = trackA.nSigmaPion();
            //Float_t TofA        = trackA.btof();
            if(nHitsPossA <= 0) nHitsPossA = 10000.0;

            Double_t dcaA_cut     = 0.5;
            Double_t dcaB_cut     = 0.5;
            Double_t VerdistX_cut = 2.5;
            Double_t VerdistY_cut = 1.5;
            Int_t    flag_proton  = 0;
            Int_t    flag_pion    = 0;
            Int_t    flag_tof_proton  = 0;
            Int_t    flag_tof_pion    = 0;

            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                flag_tof_proton = 1; // pi+
            }
            else
            {
                flag_tof_proton = 0;
            }

            if(flag_tof_proton == 1 && Mass2A < 0.1 && Mass2A > -0.5)
            {
                dcaA_cut     = 0.3;
                flag_proton  = 1;
            }
            else
            {
                flag_proton = 0;
            }

            //cout << "i = " << i << ", dcaA = " << dcaA << ", fitA = " << nHitsFitA << ", possA = " << nHitsPossA << ", pA =" << MomentumA << endl;

            if(
               dcaA         > dcaA_cut // 0.5
               && nHitsFitA > 14 // 14
               && MomentumA > 0.1
               && MomentumA < 10.0
               && (nHitsFitA/nHitsPossA) > 0.52
              )
            {

                for(Int_t j = 0; j < PID_counter_Array_B[Ana_Num][ParticleB]; j++)  // pi- candidates
                {
                    Int_t trackB_num = PID_Array_B[Ana_Num][ParticleB][j];
                    StPicoAlexTrack trackB    = *picoDst_B->track( trackB_num );
                    vectorB = trackB.origin();
                    if(fAnalysisNum == 3) vectorB = vectorB;
                    helixB = StPhysicalHelixD(trackB.gMom(),vectorB,event_B_ana->bField()*MAGFIELDFACTOR,trackB.charge());
                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    //Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction

                    Float_t Mass2B        = -100.0;
                    // calculate mass2
                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                    }


                    Float_t MassB       = SquareRoot(Mass2B);
                    //Float_t TPCdEdxB    = trackB.dEdx(); // TPC dE/dx
                    Float_t nSigmaPionB = trackB.nSigmaPion();
                    //Float_t TofB        = trackB.btof();
                    if(nHitsPossB <= 0) nHitsPossB = 10000.0;

                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        flag_tof_pion = 1;
                    }
                    else
                    {
                        flag_tof_pion = 0;
                    }

                    if(flag_tof_pion == 1 && Mass2B < 0.1 && Mass2B > -0.5) // time-of-flight is available and particle B is a pion
                    {
                        flag_pion = 1;
                    }
                    else
                    {
                        flag_pion = 0;
                    }
                    if(flag_proton == 0 && flag_pion == 0)
                    {
                        dcaA_cut     = 0.5;
                        dcaB_cut     = 0.5;
                        VerdistX_cut = 3.0;
                        VerdistY_cut = 0.8;
                    }
                    if((flag_proton == 0 && flag_pion == 1) || (flag_proton == 1 && flag_pion == 0))
                    {
                        dcaA_cut     = 0.4;
                        dcaB_cut     = 0.4;
                        VerdistX_cut = 3.0;
                        VerdistY_cut = 0.9;
                    }
                    if(flag_proton == 1 && flag_pion == 1)
                    {
                        dcaA_cut     = 0.3;
                        dcaB_cut     = 0.3;
                        VerdistX_cut = 2.0;
                        VerdistY_cut = 1.5;
                    }

                    //cout << "j = " << j << ", dcaB = " << dcaB << ", fitB = " << nHitsFitB << ", possB = " << nHitsPossB << ", pB =" << MomentumB << endl;

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       //&& dcaB > dcaA
                       && dcaA      > dcaA_cut // 1.0
                       && dcaB      > dcaB_cut // 1.0
                       && nHitsFitB > 14 // 14
                       && MomentumB > 0.1
                       && MomentumB < 10.0
                       && (nHitsFitB/nHitsPossB) > 0.52
                      )
                    {
                        Float_t pathA_f, pathB_f, dcaAB_f;
                        Float_t pathA_test = -1.0;
                        Float_t pathB_test = -1.0;
                        Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);

                        //cout << "pathA_test = " << pathA_test << ", pathB_test = " << pathB_test << endl;

                        vectoratsA     = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
                        vectoratsB     = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
                        vectorAB       = vectoratsA+vectoratsB;
                        vectorAB       = vectorAB/2.0; // decay vertex

                        //Double_t vex_est_x = vectorAB.x();
                        //Double_t vex_est_y = vectorAB.y();
                        //Double_t vex_est_z = vectorAB.z();

                        StThreeVectorF baseA,dirA,baseB,dirB;
                        baseA = helixA.at(pathA_test);
                        baseB = helixB.at(pathB_test);
                        dirA  = helixA.at(pathA_test-2.0) - helixA.at(pathA_test+2.0);
                        dirB  = helixB.at(pathB_test-2.0) - helixB.at(pathB_test+2.0);

                        StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

                        //Double_t vex_lin_x = vectorAB_lin.x();
                        //Double_t vex_lin_y = vectorAB_lin.y();
                        //Double_t vex_lin_z = vectorAB_lin.z();

                        Double_t dcaAB_lin = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
                        StThreeVectorF vectorABtoPrim_lin = vectorAB_lin - vector_prim; // vector primary vertex to decay vertex
                        Float_t VerdistX_lin = vectorABtoPrim_lin.mag(); // distance between primary vertex and decay verte


                        // calculate the scalar product with the approximated secondary vertex position
                        vectornewA_lin     = helixA.cat(pathA_test); // direction vector at dca for helixA
                        vectornewB_lin     = helixB.cat(pathB_test); // direction vector at dca for helixB
                        vectornewA_lin = MomentumA*vectornewA_lin/vectornewA_lin.mag(); // new momentum vector at decay vertex
                        vectornewB_lin = MomentumB*vectornewB_lin/vectornewB_lin.mag(); // new momentum vector at decay vertex
                        ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),0.13957018);
                        ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),0.13957018);
                        TLorentzVector trackAB_lin      = ltrackA_lin+ltrackB_lin; // mother particle
                        dirY_lin.set(trackAB_lin.Px(),trackAB_lin.Py(),trackAB_lin.Pz());
                        dirY_lin = dirY_lin/dirY_lin.mag();
                        Double_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim_lin/vectorABtoPrim_lin.mag());

                        //cout << "i = " << i << ", j = " << j << ", VerdistX_lin = " << VerdistX_lin << ", dcaAB_lin = " << dcaAB_lin << ", scalarProduct_lin = " << scalarProduct_lin << endl;

                        //if( (VerdistX_lin > 3.0 && dcaAB_lin < 1.5 && fDCA_Helix_out == 1) ||  fDCA_Helix_out == 0 )
                        if( VerdistX_lin > VerdistX_cut && dcaAB_lin < 1.5 && scalarProduct_lin > 0.0 )
                        {

                            if(fDCA_Helix_out == 1)
                            {
                                fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
                            }
                            else
                            {
                                fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
                            }

                            vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                            vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                            vectorAB       = vectoratsA+vectoratsB;
                            vectorAB       = vectorAB/2.0; // decay vertex

                            //Double_t vex_x = vectorAB.x();
                            //Double_t vex_y = vectorAB.y();
                            //Double_t vex_z = vectorAB.z();

                            //cout << "Vertex = {" << vex_x << ", " << vex_y << ", " << vex_z << "}" <<
                            //    ", Vertex_est = {" << vex_est_x << ", " << vex_est_y << ", " << vex_est_z << "}"<<
                            //    ", Vertex_lin = {" << vex_lin_x << ", " << vex_lin_y << ", " << vex_lin_z << "}"
                            //    << ", fDCA_Helix_out = " << fDCA_Helix_out << endl;

                            vectorABtoPrim = vectorAB - vector_prim; // vector primary vertex to decay vertex
                            Float_t VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

                            vectornewA     = helixA.cat(pathA_f); // direction vector at dca for helixA
                            vectornewB     = helixB.cat(pathB_f); // direction vector at dca for helixB

                            vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
                            vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex

                            ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.13957018);
                            ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),0.13957018);
                            ltrackA_pip.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.93827203);

                            //Float_t thetaA = vectoratsA.theta();
                            //Float_t thetaB = vectoratsB.theta();

                            // Missing mass and invariant mass calculations
                            TLorentzVector trackAB      = ltrackA+ltrackB; // mother particle
                            TLorentzVector trackAB_K0S  = ltrackA_pip+ltrackB; // mother particle
                            Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                            //Double_t InvMassAB_K0S      = trackAB_K0S.M(); // invariant mass of mother particle
                            Float_t MomentumAB          = trackAB.P(); // momentum of mother particle
                            Float_t BetaAB = TMath::Sqrt(1.0/(1.0+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            dirY = dirY/dirY.mag();

                            Double_t scalarProduct = dirY.dot(vectorABtoPrim/vectorABtoPrim.mag());
                            //cout << "scalarProduct = " << scalarProduct << endl;

                            baseY = vectorAB;
                            Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vector_prim);

                            //Float_t pt          = trackAB.Pt();  // Transverse momentum of mother particle
                            //Float_t rap         = trackAB.Rapidity(); // Rapidity of mother particle

                            //Float_t phiAB   = dirY.phi();
                            //Float_t thetaAB = dirY.theta();

                            // beta correction Lambda ========================================================================================
                            Float_t BetaACorr  = BetaA;
                            Float_t Mass2ACorr = Mass2A;
                            Float_t BetaBCorr  = BetaB;
                            Float_t Mass2BCorr = Mass2B;
                            Float_t MassACorr  = MassA;
                            Float_t MassBCorr  = MassB;
                            Float_t PathAB = vectorABtoPrim.mag(); // Lambda uncharged !!


                            if((trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0) || (trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0))
                            {
                                if( debug_flag )
                                {
                                    cout << "------------------------------------------------- LAMBDA --------------------------------------------------------------------------" << endl;
                                    cout << "PathAB = " << PathAB << "\tInvMassAB = " << InvMassAB << "\tMomentumAB = " << MomentumAB
                                        << "\tBetaAB = " << BetaAB << "\tTofAB = " << PathAB / (BetaAB*29.9792458)<< endl;
                                }
                            }
                            if(flag_tof_proton == 1)
                            {
                                BetaACorr = correctBeta4SingleDecay(trackA,trackAB,helixA,vectorAB,PathAB);
                                Mass2ACorr = MomentumA * MomentumA * (1.0/(BetaACorr*BetaACorr)-1.0);
                                MassACorr = SquareRoot(Mass2ACorr);
                                if( debug_flag ) cout << "A) MassA = " << MassA << "\tMassACorr = " <<  MassACorr << endl;
                            }
                            else Mass2ACorr = -100;

                            if(flag_tof_pion == 1)
                            {
                                BetaBCorr = correctBeta4SingleDecay(trackB,trackAB,helixB,vectorAB,PathAB);
                                Mass2BCorr = MomentumB * MomentumB * (1.0/(BetaBCorr*BetaBCorr)-1.0);
                                MassBCorr = SquareRoot(Mass2BCorr);
                                if( debug_flag ) cout << "B) MassB = " << MassB << "\tMassBCorr = " <<  MassBCorr << endl;
                            }
                            else Mass2BCorr = -100;


                            if(
                               InvMassAB        > 0.48
                               && InvMassAB     < 0.51
                               && VerdistX      > VerdistX_cut  // 3.5
                               && dcaAB_f       < 1.5  // 1.5
                               && VerdistY      < VerdistY_cut  // 1.5
                               && scalarProduct > 0.0
                               && ((flag_tof_proton == 0) || (Mass2ACorr > 100.0) || (Mass2ACorr  > proton_low_cut && Mass2ACorr < proton_high_cut))
                               && ((flag_tof_pion   == 0) || (Mass2BCorr > 100.0) || (Mass2BCorr  > pion_low_cut && Mass2BCorr   < pion_high_cut))
                              )
                            {
                                for(Int_t k = 0; k < PID_counter_Array_C[Ana_Num][ParticleC]; k++)  // proton candidates
                                {
                                    Int_t trackC_num = PID_Array_C[Ana_Num][ParticleC][k];
                                    StPicoAlexTrack trackC    = *picoDst_C->track( trackC_num );
                                    vectorC = event_C_ana->primaryVertex();
                                    vectorC = vectorC + vectordiff; // vectordiff == {0,0,0} for same event, don't do it for Sigma(1385) analysis
                                    helixC = StPhysicalHelixD(trackC.pMom(),vectorC,event_C_ana->bField()*MAGFIELDFACTOR,trackC.charge());

                                    Float_t MomentumC    = trackC.pMom().mag();
                                    Float_t dcaC         = trackC.dca();   // distance of closest approach to primary vertex
                                    Float_t nHitsPossC   = trackC.nHitsMax();
                                    Float_t nHitsFitC    = trackC.nHitsFit();
                                    //Float_t PolarityC    = trackC.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                    Float_t BetaC        = trackC.btofBeta();  // Velocity after time-of-flight reconstruction
                                    Float_t nSigmaProton = trackC.nSigmaProton();

                                    Float_t Mass2C       = -100.0;
                                    // calculate mass2
                                    if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                                    {
                                        Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                                    }

                                    if(
                                       trackA_num != trackC_num // Prevent that a track is used twice
                                       && trackB_num != trackC_num // Prevent that a track is used twice
                                       && dcaC      < 2.0 // 1.0
                                       && nHitsFitC > 14 // 14
                                       && MomentumC > 0.1
                                       && MomentumC < 10.0
                                       && (nHitsFitC/nHitsPossC) > 0.52
                                       && ((Mass2C < 1.5 && Mass2C > 0.45) || Mass2C < -10)
                                       && (nSigmaProton > -2.5 && nSigmaProton < 2.5)
                                      )
                                    {

                                        StThreeVectorF primdirC = helixC.cat(helixC.pathLength(vector_prim)); //
                                        primdirC = MomentumC*primdirC/primdirC.mag();
                                        ltrackC.SetXYZM(primdirC.x(),primdirC.y(),primdirC.z(),0.93827203);

                                        TLorentzVector trackABC      = trackAB+ltrackC; // mother particle
                                        Double_t InvMassABC          = trackABC.M(); // invariant mass of mother particle
                                        //Float_t MomentumABC          = trackABC.P(); // momentum of mother particle
                                        Float_t ptABC                = trackABC.Pt();  // Transverse momentum of mother particle
                                        Float_t rapABC               = trackABC.Rapidity(); // Rapidity of mother particle

                                        // Filling Ntuple
                                        ThetaPlus_NTDataArray[0]      =(Float_t)InvMassAB;
                                        ThetaPlus_NTDataArray[1]      =(Float_t)InvMassABC;
                                        ThetaPlus_NTDataArray[2]      =(Float_t)Mass2ACorr;
                                        ThetaPlus_NTDataArray[3]      =(Float_t)Mass2BCorr;
                                        ThetaPlus_NTDataArray[4]      =(Float_t)Mass2C;
                                        ThetaPlus_NTDataArray[5]      =(Float_t)nSigmaPA;
                                        ThetaPlus_NTDataArray[6]      =(Float_t)nSigmaPionB;
                                        ThetaPlus_NTDataArray[7]      =(Float_t)nSigmaProton;
                                        ThetaPlus_NTDataArray[8]      =(Float_t)MomentumA;
                                        ThetaPlus_NTDataArray[9]      =(Float_t)MomentumB;
                                        ThetaPlus_NTDataArray[10]     =(Float_t)MomentumC;
                                        ThetaPlus_NTDataArray[11]     =(Float_t)dcaA;
                                        ThetaPlus_NTDataArray[12]     =(Float_t)dcaB;
                                        ThetaPlus_NTDataArray[13]     =(Float_t)dcaC;
                                        ThetaPlus_NTDataArray[14]     =(Float_t)refMultA;
                                        ThetaPlus_NTDataArray[15]     =(Float_t)dcaAB_f;
                                        ThetaPlus_NTDataArray[16]     =(Float_t)VerdistX;
                                        ThetaPlus_NTDataArray[17]     =(Float_t)VerdistY;
                                        ThetaPlus_NTDataArray[18]     =(Float_t)ptABC;
                                        ThetaPlus_NTDataArray[19]     =(Float_t)rapABC;
                                        ThetaPlus_NTDataArray[20]     =(Float_t)RunIdA;
                                        ThetaPlus_NTDataArray[21]     =(Float_t)scalarProduct;
                                        ThetaPlus_NTDataArray[22]     =(Float_t)EventVertexXA;
                                        ThetaPlus_NTDataArray[23]     =(Float_t)EventVertexYA;
                                        ThetaPlus_NTDataArray[24]     =(Float_t)EventVertexZA;

                                        ThetaPlus_NT->Fill(ThetaPlus_NTDataArray);
                                        if(ME_Flag == 1)
                                        {
                                            Lambda_ME_counter++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        //if(fAnalysisNum == 2 || fAnalysisNum == 14) Tree_V0_v2  ->Fill();
        return 1;
    }
    else return 0;
}







void RefMult_QA_analysis(Int_t PID_counter_Array[][N_max_PIDs],Int_t PID_Array[][N_max_PIDs][N_max_tracks],StPicoAlexEvent* picoDst_A,Int_t ParticleA,Int_t Ana_Num)
{
    // Event vertex information
    event_A_ana   = picoDst_A;

    Float_t EventVertexX  = event_A_ana->primaryVertex().x();
    Float_t EventVertexY  = event_A_ana->primaryVertex().y();
    Float_t EventVertexZ  = event_A_ana->primaryVertex().z();
    Int_t   refMult       = event_A_ana->refMult();
    Int_t   RunId         = event_A_ana->runId();
    Float_t BG_rate       = event_A_ana->backgroundRate();
    Float_t ZDCx          = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx          = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd         = event_A_ana->vzVpd();
    Int_t   triggWord     = event_A_ana->triggerWord();

    //cout << "vertex = {" << EventVertexX << ", " << EventVertexY << ", " << EventVertexZ << "}, BG_rate = " << BG_rate << ", triggWord = " << triggWord << ", refMult = " << refMult << endl;

    StThreeVectorF vector_prim;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    if(
       event_A_ana   ->isMinBias()
      )
    {
        Int_t Primary_counter           = 0;
        Int_t Primary_tof_match_counter = 0;
        Int_t Global_counter            = 0;
        Float_t mean_apt                = 0.0;  // mean pt of all tracks
        Float_t mean_ppt                = 0.0;  // mean pt of primary tracks
        Float_t mean_px_QA              = 0.0;
        Float_t mean_py_QA              = 0.0;
        Float_t mean_pz_QA              = 0.0;
        Float_t EP_Qx_no_weight_QA      = 0.0;
        Float_t EP_Qy_no_weight_QA      = 0.0;
        Float_t EP_Qx_w_weight_QA       = 0.0;
        Float_t EP_Qy_w_weight_QA       = 0.0;
        Int_t   N_mean_apt              = 0;
        Int_t   N_mean_ppt              = 0;
        Float_t refMult_check           = 0.0;
        Float_t global_check            = 0.0;
        Float_t refMult_checkA          = 0.0;
        Float_t refMult_checkB          = 0.0;
        Float_t refMult_checkC          = 0.0;
        Int_t Qtracks_used_QA           = 0;
        Float_t EP_Qx_eta_pos_QA        = 0.0;
        Float_t EP_Qy_eta_pos_QA        = 0.0;
        Int_t Qtracks_used_eta_pos_QA   = 0;
        Float_t EP_Qx_eta_neg_QA        = 0.0;
        Float_t EP_Qy_eta_neg_QA        = 0.0;
        Int_t Qtracks_used_eta_neg_QA   = 0;


        if(
           PID_counter_Array[Ana_Num][ParticleA] > 0
           && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
          )
        {
            //cout << "RefMult QA with " << PID_counter_Array[Ana_Num][ParticleA] << " tracks" << endl;
            // Loop over all particle combinations
            StPhysicalHelixD helixA;
            TLorentzVector   ltrackA;
            StThreeVectorF   primdirA, vectoratsA;

            //for(Int_t n_polarity = 0; n_polarity < 2; n_polarity++) // + charges, - charges
            {

                for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
                {
                    // Get the tracks and calculate the direction and base vectors
                    Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
                    StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
                    Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle

                    //if(
                    //   (n_polarity == 0 && PolarityA > 0.0) ||
                    //   (n_polarity == 1 && PolarityA < 0.0)
                    //  ) // Check for GPC...
                    {

                        helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

                        Float_t dcaA        = trackA.dca();
                        Float_t nHitsPossA  = trackA.nHitsMax();
                        Float_t nHitsFitA   = trackA.nHitsFit();
                        Float_t MomentumA   = trackA.gMom().mag();
                        Float_t BetaA       = trackA.btofBeta();
                        //Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
                        //Float_t nSigmaKaonA = trackA.nSigmaKaon();



                        primdirA = helixA.cat(helixA.pathLength(vector_prim)); //
                        primdirA = MomentumA*primdirA/primdirA.mag();
                        ltrackA.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),0.493677);
                        Float_t etaA = primdirA.pseudoRapidity();

                        if(
                           fabs(etaA) < 0.5
                           && MomentumA  > 0.1
                           && MomentumA  < 10.0
                           && nHitsFitA  >= 15 // minimum in PicoDsts is 15
                           && (nHitsFitA/nHitsPossA) > 0.52
                          )
                        {
                            if(dcaA < 3.0)
                            {
                                refMult_check += 1.0;
                            }
                            else
                            {
                                global_check  += 1.0;
                            }
                        }

                        if(
                           MomentumA  > 0.1
                           && MomentumA  < 10.0
                           && nHitsFitA  >= 15 // minimum in PicoDsts is 15
                           && (nHitsFitA/nHitsPossA) > 0.52
                          )
                        {
                            if(
                               dcaA < 1.0
                               && fabs(etaA) < 0.5
                              )
                            {
                                refMult_checkA += 1.0;
                            }
                            if(
                               dcaA < 3.0
                               && fabs(etaA) < 0.1
                              )
                            {
                                refMult_checkB += 1.0;
                            }
                            if(
                               dcaA < 1.0
                               && fabs(etaA) < 0.1
                              )
                            {
                                refMult_checkC += 1.0;
                            }
                        }


                        Float_t Mass2A      = -100.0;
                        // calculate mass2
                        if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
                        {
                            Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
                        }

                        Float_t ptA         = ltrackA.Pt();

                        if(
                           nHitsFitA >= 15 // 15
                           && MomentumA > 0.1
                           && MomentumA < 10.0
                           && (nHitsFitA/nHitsPossA) > 0.52
                          )
                        {
                            mean_apt += ptA;
                            N_mean_apt++;

                            if(dcaA < 3.0)
                            {
                                Primary_counter++;
                                mean_ppt += ptA;
                                N_mean_ppt++;
                                if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
                                {
                                    Primary_tof_match_counter++;
                                }
                            }
                            else Global_counter++;
                        }


                        // For re-centering method
                        if(
                           (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
                           //&& fabs(EventVertexZ) < vertex_z_cut
                           && nHitsFitA          > nHitsFitA_EP_cut  // 14
                           && nHitsPossA         > nHitsPossA_EP_cut
                           && MomentumA          > MomentumA_EP_low_cut
                           && MomentumA          < MomentumA_EP_high_cut
                          )
                        {
                            if(
                               ((Float_t)nHitsFitA)/((Float_t)nHitsPossA)    > nHits_ratio_EP_cut
                              )
                            {
                                StThreeVectorF dirA;

                                Float_t pathA = -999.0;
                                Float_t dcaAB = -999.0;

                                fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);
                                vectoratsA = helixA.cat(pathA);

                                Float_t eta = vectoratsA.pseudoRapidity();

                                if(
                                   dcaA          < dcaAB_EP_cut  // 2.0
                                   && fabs(eta)  < eta_EP_cut
                                  )
                                {
                                    dirA = helixA.cat(pathA); // direction vector
                                    Double_t phiA = dirA.phi(); // phiA has values from -pi..pi

                                    Float_t p_x = MomentumA*dirA.x();
                                    Float_t p_y = MomentumA*dirA.y();
                                    Float_t p_z = MomentumA*dirA.z();
                                    Float_t p_t = sqrt(p_x*p_x + p_y*p_y);

                                    Double_t p_t_weight = 1.0;
                                    if(p_t < 2.0)  p_t_weight = p_t;
                                    if(p_t >= 2.0) p_t_weight = 2.0;

                                    mean_px_QA += p_x;
                                    mean_py_QA += p_y;
                                    mean_pz_QA += p_z;

                                    Float_t iQx = TMath::Cos(2.0*phiA);
                                    Float_t iQy = TMath::Sin(2.0*phiA);

                                    EP_Qx_no_weight_QA  += iQx;
                                    EP_Qy_no_weight_QA  += iQy;
                                    EP_Qx_w_weight_QA   += p_t_weight*iQx;
                                    EP_Qy_w_weight_QA   += p_t_weight*iQy;

                                    Qtracks_used_QA++;

                                    // For eta sub-event method
                                    if(
                                       fabs(eta) > eta_gap
                                      )
                                    {
                                        if(eta >= 0.0)
                                        {
                                            EP_Qx_eta_pos_QA += p_t_weight*iQx;
                                            EP_Qy_eta_pos_QA += p_t_weight*iQy;
                                            Qtracks_used_eta_pos_QA++;
                                        }
                                        if(eta < 0.0)
                                        {
                                            EP_Qx_eta_neg_QA += p_t_weight*iQx;
                                            EP_Qy_eta_neg_QA += p_t_weight*iQy;
                                            Qtracks_used_eta_neg_QA++;
                                        }
                                    }

                                    // Phi-weight correction
                                    // determine the file time
                                    TString name;
                                    Long64_t val;
                                    Float_t day = 0.0;
                                    char NoP[50];
                                    sprintf(NoP,"%u",RunId);
                                    name = NoP[2];
                                    sscanf(name.Data(),"%Li",&val);
                                    day = 100.0 * val;
                                    name = NoP[3];
                                    sscanf(name.Data(),"%Li",&val);
                                    day = day + 10.0 * val;
                                    name = NoP[4];
                                    sscanf(name.Data(),"%Li",&val);
                                    day = day + 1.0 * val;

                                    Int_t day_bin     = (Int_t)day;
                                    Int_t z_bin       = (Int_t)((EventVertexZ-phi_corr_z_start)/phi_corr_delta_z);
                                    Int_t eta_bin     = (Int_t)((eta-phi_corr_eta_start)/phi_corr_delta_eta);
                                    Int_t pt_bin      = (Int_t)((p_t-phi_corr_pt_start)/phi_corr_delta_pt);
                                    Int_t pol_bin     = 0;
                                    if(PolarityA > 0) pol_bin = 0;
                                    if(PolarityA < 0) pol_bin = 1;

                                    if(pt_bin >= nPhi_corr_pt_bins) pt_bin = nPhi_corr_pt_bins-1;

                                    //cout << "z = " << EventVertexZ << ", z_bin = " << z_bin << endl;

                                    if(
                                       day_bin    >= 0
                                       && eta_bin >= 0
                                       && z_bin   >= 0
                                       && pt_bin  >= 0
                                       && day_bin < nPhi_corr_days
                                       && eta_bin < nPhi_corr_eta_bins
                                       && pt_bin  < nPhi_corr_pt_bins
                                       && z_bin   < nPhi_corr_z_bins
                                       && fabs(EventVertexZ) < vertex_z_cut
                                      )
                                    {
                                        hPhi_days_use[day_bin] = 1;
                                        hPhi_corr[day_bin][z_bin][eta_bin][pol_bin][pt_bin] ->Fill(phiA);
                                    }


                                }
                            }
                        }

                    }

                }

                //cout << "z = " << EventVertexZ << ", refMult = " << refMult << ", refMult_check = " << refMult_check << endl;

                if(N_mean_apt > 0) mean_apt/=((Double_t)N_mean_apt);
                if(N_mean_ppt > 0) mean_ppt/=((Double_t)N_mean_ppt);

                RefMult_QA_NTDataArray[0]    =(Float_t)refMult;
                RefMult_QA_NTDataArray[1]    =(Float_t)EventVertexX;
                RefMult_QA_NTDataArray[2]    =(Float_t)EventVertexY;
                RefMult_QA_NTDataArray[3]    =(Float_t)EventVertexZ;
                RefMult_QA_NTDataArray[4]    =(Float_t)RunId;
                RefMult_QA_NTDataArray[5]    =(Float_t)Primary_counter;
                RefMult_QA_NTDataArray[6]    =(Float_t)Primary_tof_match_counter;
                RefMult_QA_NTDataArray[7]    =(Float_t)Global_counter;
                RefMult_QA_NTDataArray[8]    =(Float_t)mean_apt;
                RefMult_QA_NTDataArray[9]    =(Float_t)mean_ppt;
                RefMult_QA_NTDataArray[10]   =(Float_t)BG_rate;
                RefMult_QA_NTDataArray[11]   =(Float_t)ZDCx;
                RefMult_QA_NTDataArray[12]   =(Float_t)BBCx;
                RefMult_QA_NTDataArray[13]   =(Float_t)vzVpd;
                RefMult_QA_NTDataArray[14]   =(Float_t)EP_Qx_eta_pos;
                RefMult_QA_NTDataArray[15]   =(Float_t)EP_Qy_eta_pos;
                RefMult_QA_NTDataArray[16]   =(Float_t)EP_Qx_eta_neg;
                RefMult_QA_NTDataArray[17]   =(Float_t)EP_Qy_eta_neg;
                RefMult_QA_NTDataArray[18]   =(Float_t)refMult_check;
                RefMult_QA_NTDataArray[19]   =(Float_t)global_check;
                RefMult_QA_NTDataArray[20]   =(Float_t)refMult_checkA;
                RefMult_QA_NTDataArray[21]   =(Float_t)refMult_checkB;
                RefMult_QA_NTDataArray[22]   =(Float_t)refMult_checkC;
                RefMult_QA_NTDataArray[23]   =(Float_t)triggWord; // triggerWord
                RefMult_QA_NTDataArray[24]   =(Float_t)EP_Qx;
                RefMult_QA_NTDataArray[25]   =(Float_t)EP_Qy;
                RefMult_QA_NTDataArray[26]   =(Float_t)Qtracks_used_eta_pos;
                RefMult_QA_NTDataArray[27]   =(Float_t)Qtracks_used_eta_neg;
                RefMult_QA_NTDataArray[28]   =(Float_t)Qtracks_used;
                RefMult_QA_NTDataArray[29]   =(Float_t)EP_Qx_w_weight_QA;
                RefMult_QA_NTDataArray[30]   =(Float_t)EP_Qy_w_weight_QA;
                RefMult_QA_NTDataArray[31]   =(Float_t)EP_Qx_eta_pos_QA;
                RefMult_QA_NTDataArray[32]   =(Float_t)EP_Qy_eta_pos_QA;
                RefMult_QA_NTDataArray[33]   =(Float_t)Qtracks_used_eta_pos_QA;
                RefMult_QA_NTDataArray[34]   =(Float_t)EP_Qx_eta_neg_QA;
                RefMult_QA_NTDataArray[35]   =(Float_t)EP_Qy_eta_neg_QA;
                RefMult_QA_NTDataArray[36]   =(Float_t)Qtracks_used_eta_neg_QA;
                RefMult_QA_NTDataArray[37]   =(Float_t)Qtracks_used_QA;

                RefMult_QA_NT->Fill(RefMult_QA_NTDataArray);

            } // end of polarity loop

            //cout << "EP_Qx_eta_pos_QA = " << EP_Qx_eta_pos_QA << ", EP_Qy_eta_pos_QA = " <<  EP_Qy_eta_pos_QA
            //    << ", Qtracks_used_eta_pos_QA = " << Qtracks_used_eta_pos_QA << ", refMult = " << refMult << endl;

        }
    }
}


void Ach_analysis(Int_t PID_counter_Array[][N_max_PIDs],Int_t PID_Array[][N_max_PIDs][N_max_tracks],StPicoAlexEvent* picoDst_A,Int_t ParticleA,Int_t Ana_Num)
{
    // Event vertex information
    event_A_ana   = picoDst_A;

    Float_t EventVertexX  = event_A_ana->primaryVertex().x();
    Float_t EventVertexY  = event_A_ana->primaryVertex().y();
    Float_t EventVertexZ  = event_A_ana->primaryVertex().z();
    Int_t   refMult       = event_A_ana->refMult();
    Int_t   RunId         = event_A_ana->runId();
    Float_t BG_rate       = event_A_ana->backgroundRate();
    Float_t ZDCx          = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx          = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd         = event_A_ana->vzVpd();
    Int_t   triggWord     = event_A_ana->triggerWord();

    //cout << "vertex = {" << EventVertexX << ", " << EventVertexY << ", " << EventVertexZ << "}, BG_rate = " << BG_rate << ", triggWord = " << triggWord << ", refMult = " << refMult << endl;

    StThreeVectorF vector_prim;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    const Int_t max_tracks = 2000;
    Double_t pt_tracks[2][max_tracks]; // pt values for both charges for every single track


    // 0: 70%-80%
    // 1: 60%-70%
    // 2: 50%-60%
    // 3: 40%-50%
    // 4: 30%-40%
    // 5: 20%-30%
    // 6: 10%-20%
    // 7: 5%-10%
    // 8: 0%-5%

    if(
       event_A_ana   ->isMinBias()
       && erefMult_bin > 3 && erefMult_bin < 5 // 30%-40% central
      )
    {
        if(
           PID_counter_Array[Ana_Num][ParticleA] > 0
           && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
           && fabs(EventVertexZ) < 30.0
           && fabs(vzVpd-EventVertexZ) < 3.0
          )
        {
            // Loop over all particle combinations
            StPhysicalHelixD helixA;
            TLorentzVector   ltrackA;
            StThreeVectorF   primdirA, vectoratsA;

            Int_t N_pos_neg_charges_reconstructed[2]     = {0,0};
            Int_t N_pos_neg_charges_reconstructed_Ach[2] = {0,0};

            for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
            {
                // Get the tracks and calculate the direction and base vectors
                Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
                StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );

                //helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());
                helixA = StPhysicalHelixD(trackA.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

                Float_t dcaA        = trackA.dca();
                Float_t nHitsPossA  = trackA.nHitsMax();
                Float_t nHitsFitA   = trackA.nHitsFit();
                Float_t MomentumA   = trackA.gMom().mag();
                Float_t BetaA       = trackA.btofBeta();
                Float_t nSigmaKaonA = trackA.nSigmaKaon();
                Float_t nSigmaPionA = trackA.nSigmaPion();
                Float_t nSigmaPA    = trackA.nSigmaProton();
                Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle

                Int_t charge = 0;
                if(PolarityA < 0.0) charge = 1;

                primdirA = helixA.cat(helixA.pathLength(vector_prim)); //
                primdirA = MomentumA*primdirA/primdirA.mag();
                ltrackA.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),0.13957018); // Default a pion
                if(fabs(nSigmaKaonA) < 2.0) ltrackA.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),0.493677);
                if(fabs(nSigmaPA)    < 2.0) ltrackA.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),0.938272013);

                Float_t etaA = primdirA.pseudoRapidity();
                Float_t pT   = ltrackA.Pt();

                // Cuts to determine Ach
                if(
                   fabs(etaA) < 1.0
                   && pT > 0.15
                   && pT < 12.0
                   && !(pT < 0.4 && fabs(nSigmaPA) < 3.0)
                   && dcaA < 1.0
                   && nHitsFitA >= 15
                   && (nHitsPossA/nHitsFitA) > 0.52
                   && N_pos_neg_charges_reconstructed_Ach[charge] < max_tracks
                  )
                {
                    N_pos_neg_charges_reconstructed_Ach[charge]++;
                }

                // Cuts to determine pion pT
                if(
                   fabs(etaA) < 1.0
                   && pT > 0.15
                   && pT < 0.5
                   && dcaA < 1.0
                   && nHitsFitA >= 15
                   && (nHitsPossA/nHitsFitA) > 0.52
                   && N_pos_neg_charges_reconstructed[charge] < max_tracks
                   && fabs(nSigmaPionA) < 2.0
                  )
                {
                    //cout << "pT: " << pT << endl;
                    pt_tracks[charge][N_pos_neg_charges_reconstructed[charge]] = pT;
                    N_pos_neg_charges_reconstructed[charge]++;
                }
            } // end of track loop

            if(N_pos_neg_charges_reconstructed_Ach[0] + N_pos_neg_charges_reconstructed_Ach[1] > 0)
            {
                Double_t Ach  = ((Double_t)(N_pos_neg_charges_reconstructed_Ach[0] - N_pos_neg_charges_reconstructed_Ach[1]))/((Double_t)(N_pos_neg_charges_reconstructed_Ach[0] + N_pos_neg_charges_reconstructed_Ach[1]));
                Double_t Mult = (Double_t)(N_pos_neg_charges_reconstructed_Ach[0] + N_pos_neg_charges_reconstructed_Ach[1]);

                h_Ach           ->Fill(Ach);

                for(Int_t i_pos_neg = 0; i_pos_neg < 2; i_pos_neg++) // loop first over all positive and then over all negative charges
                {
                    for(Int_t i_tracks = 0; i_tracks < N_pos_neg_charges_reconstructed[i_pos_neg]; i_tracks++)
                    {
                        p_pt_Ach[i_pos_neg] ->Fill(Ach,pt_tracks[i_pos_neg][i_tracks]);
                    }
                }

            }

        }
    }
}


Int_t Hadron_v2_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                   Int_t ParticleA, Int_t Ana_Num,Int_t run_events,Int_t event_number)
{
    //cout << "Hadron_v2_analysis started" << endl;
    event_A_ana   = picoDst_A;

    // Event vertex information
    Float_t EventVertexX  = event_A_ana->primaryVertex().x();
    Float_t EventVertexY  = event_A_ana->primaryVertex().y();
    Float_t EventVertexZ  = event_A_ana->primaryVertex().z();
    Int_t   refMult       = event_A_ana->refMult();
    Int_t   RunId         = event_A_ana->runId();
    Float_t ZDCx          = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx          = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd         = event_A_ana->vzVpd();
    StThreeVectorF vector_prim, vector_prim_new;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    if(
       PID_counter_Array[Ana_Num][ParticleA] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && event_A_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " particles *************" << endl;
        //cout << "NVertices = " << event->getNVertices() << ", Vertex = {" << EventVertexX << ", " << EventVertexY << ", " << EventVertexZ << "}" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA;
        StThreeVectorF vectorA, vectorAB, vector_primAB, vectornewA, vectornewB, vector_sum, dirA;
        StThreeVectorF testA, testAB, vectorABtoPrim, baseY, dirY;

        Int_t Qtracks_used_hadron_v2 = 0;

        alex_event.clearTrackList();
        alex_event.setx(EventVertexX);
        alex_event.sety(EventVertexY);
        alex_event.setz(EventVertexZ);
        alex_event.setid(RunId);
        alex_event.setmult(refMult);
        alex_event.setn_prim(n_primaries);
        alex_event.setn_non_prim(n_non_primaries);
        alex_event.setn_tof_prim(n_tofmatch_prim);
        alex_event.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        alex_event.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        alex_event.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        alex_event.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        alex_event.setEP_Qx_ptw(EP_Qx_ptw);
        alex_event.setEP_Qy_ptw(EP_Qy_ptw);
        alex_event.setQtracks_eta_pos(Qtracks_used_eta_pos);
        alex_event.setQtracks_eta_neg(Qtracks_used_eta_neg);
        alex_event.setQtracks_full(Qtracks_used);
        alex_event.setZDCx(ZDCx);
        alex_event.setBBCx(BBCx);
        alex_event.setvzVpd(vzVpd);

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            //cout << "i = " << i << endl;
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t BetaA       = trackA.btofBeta();
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
            Float_t nSigmaKaonA = trackA.nSigmaKaon();
            Float_t nSigmaPionA = trackA.nSigmaPion();
            Float_t nSigmaPA    = trackA.nSigmaProton();
            Float_t nSigmaElA   = trackA.nSigmaElectron();

            //primdirA = helixA.cat(helixA.pathLength(vector_prim)); //
            //primdirA = MomentumA*primdirA/primdirA.mag();
            //ltrackA.SetXYZM(primdirA.x(),primdirA.y(),primdirA.z(),0.493677);
            //Float_t etaA = primdirA.pseudoRapidity();

            Float_t Mass2A      = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }

            if(
               nHitsFitA     >= 15  // 15
               && nHitsPossA > 0
               && MomentumA  > 0.1  // 0.15
               && MomentumA  < 10.0
              )
            {
                if(
                   ((Float_t)nHitsFitA)/((Float_t)nHitsPossA)    > 0.52
                  )
                {
                    Float_t pathA = -999.0;
                    Float_t dcaAB = -999.0;
                    fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);

                    dirA = helixA.cat(pathA); // direction vector

                    Float_t eta = dirA.pseudoRapidity();

                    if(
                       dcaA          < 3.0 // 2.0
                       && fabs(eta)  < 2.0 // 1.0
                      )
                    {

                        Float_t p_x = MomentumA*dirA.x();
                        Float_t p_y = MomentumA*dirA.y();
                        Float_t p_t = sqrt(p_x*p_x + p_y*p_y);


                        Double_t phiA   = dirA.phi();
                        Double_t thetaA = dirA.theta();


                        Float_t iQx = TMath::Cos(2.0*phiA);
                        Float_t iQy = TMath::Sin(2.0*phiA);

                        Double_t phi_w;
                        Double_t total_weight = calc_event_plane_weight(phiA,p_t,eta,RunId,EventVertexZ,(Int_t)PolarityA,phi_w);
                        Double_t total_weight_eta_gap = total_weight;

                        if(!(
                             dcaAB                     < dcaAB_EP_cut
                             && fabs(eta)              < eta_EP_cut
                             && MomentumA              > MomentumA_EP_low_cut
                             && MomentumA              < MomentumA_EP_high_cut
                             && (nHitsFitA/nHitsPossA) > nHits_ratio_EP_cut
                             && nHitsFitA              > nHitsFitA_EP_cut
                             && nHitsPossA             > nHitsPossA_EP_cut
                             && fabs(EventVertexZ)     < vertex_z_cut
                             && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
                            )
                          )
                        {
                            total_weight         = 0.0;  // track was not used in EP calculation
                            total_weight_eta_gap = 0.0;  // track was not used in eta gap EP calculation
                            if(
                               !(fabs(eta) > eta_gap)
                              )
                            {
                                total_weight_eta_gap = 0.0;  // track was not used in eta gap EP calculation
                            }
                        }

                        //cout << "Hadron v2, phiA = " << phiA << ", p_t = " << p_t << ", eta = " << eta << ", RunId = "  << RunId
                        //    << ", z = " << EventVertexZ << ", total_weight = " << total_weight << ", refMult = " << refMult << endl;

                        Float_t EP_Qx_new = EP_Qx;
                        Float_t EP_Qy_new = EP_Qy;

                        Float_t EP_Qx_new_eta_gap = 0;
                        Float_t EP_Qy_new_eta_gap = 0;

                        if(eta < 0.0) // Use positive eta EP for netagive eta particle track
                        {
                            EP_Qx_new_eta_gap = EP_Qx_eta_pos;
                            EP_Qy_new_eta_gap = EP_Qy_eta_pos;
                        }
                        if(eta > 0.0) // Use negative eta EP for positive eta particle track
                        {
                            EP_Qx_new_eta_gap = EP_Qx_eta_neg;
                            EP_Qy_new_eta_gap = EP_Qy_eta_neg;
                        }

                        // Subtract the Q-vector from the event plane vector to avoid auto correlation
                        EP_Qx_new = EP_Qx_new - total_weight*iQx;
                        EP_Qy_new = EP_Qy_new - total_weight*iQy;

                        // Particle is not used for EP calculation by construction!
                        //EP_Qx_new_eta_gap = EP_Qx_new_eta_gap - total_weight_eta_gap*iQx;
                        //EP_Qy_new_eta_gap = EP_Qy_new_eta_gap - total_weight_eta_gap*iQy;


                        // Event plane has values from 0..pi
                        Float_t phi_event_plane         = -100.0;
                        Float_t phi_event_plane_raw     = -100.0;
                        Float_t phi_event_plane_eta_gap = -100.0;
                        phi_event_plane         = calc_phi_event_plane_2nd(EP_Qx_new,EP_Qy_new);
                        phi_event_plane_raw     = calc_phi_event_plane_2nd(EP_Qx,EP_Qy);
                        phi_event_plane_eta_gap = calc_phi_event_plane_2nd(EP_Qx_new_eta_gap,EP_Qy_new_eta_gap);

                        Float_t delta_phi = 0.0;

                        //cout << "phi_event_plane = " << phi_event_plane << ", phiA = " << phiA << endl;
                        // phiA has values from -pi..pi
                        if(phiA > phi_event_plane)
                        {
                            delta_phi = phiA-phi_event_plane;
                        }
                        else
                        {
                            delta_phi = 2.0*TMath::Pi()+(phiA-phi_event_plane);
                        }
                        delta_phi = delta_phi - TMath::Pi(); // delta_phi has values from -pi..pi

                        Qtracks_used_hadron_v2++;


                        //if(fabs(delta_phi) > 1.570 && fabs(delta_phi) < 1.571)
                        //{
                        //    cout << "i = " << i << ", delta_phi = " << delta_phi*TMath::RadToDeg() << ", EP = " << phi_event_plane*TMath::RadToDeg()
                        //        << ", phi = " << phiA*TMath::RadToDeg()
                        //        << ", EP_Qx = " << EP_Qx << ", EP_Qy = " << EP_Qy << ", EP_Qx_new = " << EP_Qx_new << ", EP_Qy_new = " << EP_Qy_new
                        //        << ", qx = " << total_weight*iQx << ", qy = " << total_weight*iQy << endl;
                        //}

                        alex_track = alex_event.createTrack();
                        alex_track->setphi(phiA);
                        alex_track->settheta(thetaA);
                        alex_track->setp_t(p_t);
                        alex_track->seteta(eta);
                        alex_track->setdca(dcaAB);
                        alex_track->setpq(MomentumA*PolarityA);
                        alex_track->setm2(Mass2A);
                        alex_track->setnsK(nSigmaKaonA);
                        alex_track->setnsPi(nSigmaPionA);
                        alex_track->setnsP(nSigmaPA);
                        alex_track->setnsEl(nSigmaElA);
                        alex_track->setdEdx(TPCdEdxA);
                        alex_track->setiQx(iQx);
                        alex_track->setiQy(iQy);
                        alex_track->setphi_EP(phi_event_plane);
                        alex_track->setphi_EP_raw(phi_event_plane_raw);
                        alex_track->setphi_EP_eta_gap(phi_event_plane_eta_gap);

                        //cout << "Tree_hadron_v2 filled, phiA = " << phiA << endl;

                        /*
                        // Track based quantities
                        Hadron_v2_NTDataArray[0]     =(Float_t)phiA;
                        Hadron_v2_NTDataArray[1]     =(Float_t)p_t;
                        Hadron_v2_NTDataArray[2]     =(Float_t)eta;
                        Hadron_v2_NTDataArray[3]     =(Float_t)dcaAB;
                        Hadron_v2_NTDataArray[4]     =(Float_t)MomentumA*PolarityA;
                        Hadron_v2_NTDataArray[5]     =(Float_t)Mass2A;
                        Hadron_v2_NTDataArray[6]     =(Float_t)nSigmaKaonA;
                        Hadron_v2_NTDataArray[7]     =(Float_t)nSigmaPionA;
                        Hadron_v2_NTDataArray[8]     =(Float_t)nSigmaPA;
                        Hadron_v2_NTDataArray[9]     =(Float_t)nSigmaElA;
                        Hadron_v2_NTDataArray[10]    =(Float_t)TPCdEdxA;
                        Hadron_v2_NTDataArray[11]    =(Float_t)iQx;
                        Hadron_v2_NTDataArray[12]    =(Float_t)iQy;
                        Hadron_v2_NTDataArray[13]    =(Float_t)phi_event_plane;
                        Hadron_v2_NTDataArray[14]    =(Float_t)phi_event_plane_raw;
                        Hadron_v2_NTDataArray[25]    =(Float_t)phi_event_plane_eta_gap;

                        // Event based quantities
                        Hadron_v2_NTDataArray[13]    =(Float_t)EventVertexX;
                        Hadron_v2_NTDataArray[14]    =(Float_t)EventVertexY;
                        Hadron_v2_NTDataArray[15]    =(Float_t)EventVertexZ;
                        Hadron_v2_NTDataArray[16]    =(Float_t)RunId;
                        Hadron_v2_NTDataArray[17]    =(Float_t)refMult;
                        Hadron_v2_NTDataArray[21]    =(Float_t)n_primaries;
                        Hadron_v2_NTDataArray[22]    =(Float_t)n_non_primaries;
                        Hadron_v2_NTDataArray[23]    =(Float_t)n_tofmatch_prim;
                        Hadron_v2_NTDataArray[24]    =(Float_t)EP_Qx_eta_pos_ptw;
                        Hadron_v2_NTDataArray[25]    =(Float_t)EP_Qy_eta_pos_ptw;
                        Hadron_v2_NTDataArray[26]    =(Float_t)EP_Qx_eta_neg_ptw;
                        Hadron_v2_NTDataArray[27]    =(Float_t)EP_Qy_eta_neg_ptw;
                        Hadron_v2_NTDataArray[28]    =(Float_t)EP_Qx_ptw;
                        Hadron_v2_NTDataArray[29]    =(Float_t)EP_Qy_ptw;
                        Hadron_v2_NTDataArray[30]    =(Float_t)Qtracks_used_eta_pos;
                        Hadron_v2_NTDataArray[31]    =(Float_t)Qtracks_used_eta_neg;
                        Hadron_v2_NTDataArray[32]    =(Float_t)Qtracks_used;

                        Hadron_v2_NT->Fill(Hadron_v2_NTDataArray);
                        */
                    }

                }
            }
        }

        Tree_hadron_v2  ->Fill();

        return 1;
    }
    else return 0;
}




Int_t Hadron_spectra_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                   Int_t ParticleA, Int_t Ana_Num,Int_t run_events,Int_t event_number)
{
    //cout << "Hadron_v2_analysis started" << endl;
    event_A_ana   = picoDst_A;

    // Event vertex information
    Float_t EventVertexX  = event_A_ana->primaryVertex().x();
    Float_t EventVertexY  = event_A_ana->primaryVertex().y();
    Float_t EventVertexZ  = event_A_ana->primaryVertex().z();
    Int_t   refMult       = event_A_ana->refMult();
    Int_t   RunId         = event_A_ana->runId();
    Float_t ZDCx          = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx          = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd         = event_A_ana->vzVpd();
    StThreeVectorF vector_prim, vector_prim_new;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    if(
       PID_counter_Array[Ana_Num][ParticleA] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && event_A_ana   ->isMinBias()
       && n_tofmatch_prim >= 2
       && erefMult_bin >= 0 && erefMult_bin < 9
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " particles *************" << endl;
        //cout << "NVertices = " << event->getNVertices() << ", Vertex = {" << EventVertexX << ", " << EventVertexY << ", " << EventVertexZ << "}" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA;
        StThreeVectorF vectorA, vectorAB, vector_primAB, vectornewA, vectornewB, vector_sum, dirA;
        StThreeVectorF testA, testAB, vectorABtoPrim, baseY, dirY;

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            //cout << "i = " << i << endl;
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            helixA = StPhysicalHelixD(trackA.pMom(),event_A_ana->primaryVertex(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.pMom().mag();
            Float_t BetaA       = trackA.btofBeta();
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
            Float_t nSigmaKaonA = trackA.nSigmaKaon();
            Float_t nSigmaPionA = trackA.nSigmaPion();
            Float_t nSigmaPA    = trackA.nSigmaProton();
            Float_t nSigmaElA   = trackA.nSigmaElectron();
            Float_t nHitsDedxA  = trackA.nHitsDedx();
            Float_t tofYLocal   = trackA.btofYLocal();
            Float_t tofZLocal   = trackA.btofZLocal();

            Int_t charge_bin = 0;
            if(PolarityA < 0.0) charge_bin = 1;

            Float_t Mass2A      = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }

            if(
               nHitsFitA     > 26  // 15
               && nHitsDedxA > 16
               && nHitsPossA > 0
               && MomentumA  > 0.1  // 0.15
               && MomentumA  < 10.0
              )
            {
                if(
                   ((Float_t)nHitsFitA)/((Float_t)nHitsPossA) > 0.52
                  )
                {
                    Float_t pathA = -999.0;
                    Float_t dcaAB = -999.0;
                    fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);
                    TLorentzVector ltrackA;
                    StThreeVectorF dirA = helixA.cat(pathA); // direction vector
                    dirA *= MomentumA/dirA.mag();
                    ltrackA.SetXYZM(dirA.x(),dirA.y(),dirA.z(),mass_array[ParticleA]);
                    Float_t eta      = dirA.pseudoRapidity();
                    Float_t rapidity = ltrackA.Rapidity();
                    Float_t p_t      = ltrackA.Pt();

                    Int_t pt_bin = -1;
                    for(Int_t r = 0; r < (npt_bins_pt_spectra-1); r++)
                    {
                        if(p_t >= pt_bin_ranges_spectra[0])
                        {
                            if(p_t >= pt_bin_ranges_spectra[r] && p_t < pt_bin_ranges_spectra[r+1])
                            {
                                pt_bin = r;
                                break;
                            }
                        }
                        else break;
                    }

                    if(
                       dcaA              < 3.0 // 2.0
                       && fabs(rapidity) < 0.05 // 1.0
                       && pt_bin >= 0
                       && pt_bin < npt_bins_pt_spectra
                      )
                    {
                        if(p_t < 0.6 && erefMult_bin > 6)
                        {
                            h_ToF_local_YZ_vs_m2[0] ->Fill(Mass2A,tofYLocal);
                            h_ToF_local_YZ_vs_m2[1] ->Fill(Mass2A,tofZLocal);
                        }

                        h_m2_mult_pt[erefMult_bin][pt_bin][charge_bin][0]        ->Fill(Mass2A);

                        if(fabs(tofYLocal) < 1.7)
                        {
                            combPID->setInitValues(p_t,nSigmaPionA*nsigma_scaling_fac,Mass2A);
                            Double_t NewX = combPID->getNewX();
                            Double_t NewY = combPID->getNewY();

                            h_NewY_vs_NewX_mult_pt[erefMult_bin][pt_bin][charge_bin] ->Fill(NewX,NewY);
                            h_m2_mult_pt[erefMult_bin][pt_bin][charge_bin][1]        ->Fill(Mass2A);
                        }
                    }

                }
            }
        }

        return 1;
    }
    else return 0;
}


Int_t Fill_m2_nSigmaP_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                   Int_t ParticleA, Int_t Ana_Num,Int_t run_events,Int_t event_number)
{
    //cout << "Hadron_v2_analysis started" << endl;
    event_A_ana   = picoDst_A;

    // Event vertex information
    Float_t EventVertexX  = event_A_ana->primaryVertex().x();
    Float_t EventVertexY  = event_A_ana->primaryVertex().y();
    Float_t EventVertexZ  = event_A_ana->primaryVertex().z();
    Int_t   refMult       = event_A_ana->refMult();
    Int_t   RunId         = event_A_ana->runId();
    Float_t ZDCx          = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx          = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd         = event_A_ana->vzVpd();
    StThreeVectorF vector_prim, vector_prim_new;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    if(
       PID_counter_Array[Ana_Num][ParticleA] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && event_A_ana   ->isMinBias()
      )
    {
        StPhysicalHelixD helixA;
        StThreeVectorF vectorA, vectorAB, vector_primAB, vectornewA, vectornewB, vector_sum, dirA;
        StThreeVectorF testA, testAB, vectorABtoPrim, baseY, dirY;

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            //cout << "i = " << i << endl;
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t BetaA       = trackA.btofBeta();
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
            Float_t nSigmaKaonA = trackA.nSigmaKaon();
            Float_t nSigmaPionA = trackA.nSigmaPion();
            Float_t nSigmaPA    = trackA.nSigmaProton();
            Float_t nSigmaElA   = trackA.nSigmaElectron();

            Int_t charge = 0;
            if(PolarityA > 0) {charge = 0;}
            else {charge = 1;}

            //primdirA = helixA.cat(helixA.pathLength(vector_prim)); //
            //primdirA = MomentumA*primdirA/primdirA.mag();
            //Float_t etaA = primdirA.pseudoRapidity();

            Float_t Mass2A      = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }

            if(
               nHitsFitA     > 15  // 15
               && nHitsFitA  < 50  // 15
               && nHitsPossA > 0
               && MomentumA  > 0.4  // 0.15
               && MomentumA  < 3.4
              )
            {
                if(
                   ((Float_t)nHitsFitA)/((Float_t)nHitsPossA) > 0.52
                  )
                {
                    Float_t pathA = -999.0;
                    Float_t dcaAB = -999.0;
                    fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);

                    dirA = helixA.cat(pathA); // direction vector
                    dirA *= MomentumA;
                    TLorentzVector ltrackA;
                    ltrackA.SetXYZM(dirA.x(),dirA.y(),dirA.z(),mass_array[ParticleA]);

                    Float_t rapidity = ltrackA.Rapidity();
                    Float_t p_x      = dirA.x();
                    Float_t p_y      = dirA.y();
                    Float_t p_t      = sqrt(p_x*p_x + p_y*p_y);
                    Float_t eta      = dirA.pseudoRapidity();
                    Float_t phiA     = dirA.phi();
                    Float_t thetaA   = dirA.theta();

                    Int_t pt_bin = (Int_t)((p_t-low_pT_m2_nSigmaP)/delta_pT_m2_nSigmaP);

                    if(
                       dcaA   < 2.0 // 2.0
                       && p_t > 0.4
                       && p_t < 3.4
                       && rapidity > -0.7
                       && rapidity < 0.7
                       && pt_bin >= 0
                       && pt_bin < N_2D_m2_nSigmaP_pT_bins
                       && erefMult_bin >= 0
                       && erefMult_bin < N_refmult_bins
                       && Mass2A > -1.0
                       && Mass2A < 10.0
                      )
                    {
                        h2D_m2_nSigmaP[erefMult_bin][pt_bin][charge] ->Fill(nSigmaPA,Mass2A);
                    }

                }
            }
        }

        return 1;
    }
    else return 0;
}


Int_t fill_phi_hist_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                   Int_t ParticleA, Int_t ParticleB, Int_t Ana_Num,Int_t run_events,Int_t event_number)
{
    //cout << "Hadron_v2_analysis started" << endl;
    event_A_ana   = picoDst_A;

    // Event vertex information
    Float_t EventVertexX  = event_A_ana->primaryVertex().x();
    Float_t EventVertexY  = event_A_ana->primaryVertex().y();
    Float_t EventVertexZ  = event_A_ana->primaryVertex().z();
    Int_t   refMult       = event_A_ana->refMult();
    Int_t   RunId         = event_A_ana->runId();
    Int_t   mult_bin      = refmultCorrUtil->getCentralityBin9();
    StThreeVectorF vector_prim, vector_prim_new;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    if(
       PID_counter_Array[Ana_Num][ParticleA]    > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && mult_bin >= 0 && mult_bin < 9
       && event_A_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " particles *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB;
        StThreeVectorF vectorA, dirA, dirB;


        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // particle loop A
        {
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t BetaA       = trackA.btofBeta();
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t nSigmaPA    = trackA.nSigmaProton();

            Float_t pathA = -999.0;
            Float_t dcaAB = -999.0;
            fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);

            dirA = helixA.cat(pathA); // direction vector

            Float_t phiA = dirA.phi();
            Float_t p_xA = MomentumA*dirA.x();
            Float_t p_yA = MomentumA*dirA.y();
            Float_t p_tA = sqrt(p_xA*p_xA + p_yA*p_yA);

            if(phiA < 0.0) phiA = 2.0*TMath::Pi() + phiA;

            Int_t charge_bin = 0;
            if(PolarityA < 0.0) charge_bin = 1;

            if(
               mult_bin >= 0 && mult_bin < N_refmult_bins
               && nHitsFitA  >= 15  // 15
               && nHitsPossA > 0
               && MomentumA  > 0.1  // 0.15
               && MomentumA  < 10.0
              )
            {
                h_phi_vs_track_number[mult_bin][charge_bin] ->Fill((Float_t)i,phiA);
                h_pT_vs_track_number[mult_bin][charge_bin]  ->Fill((Float_t)i,p_tA);
                h_mean_pt_vs_phi[mult_bin][charge_bin]      ->Fill(phiA,p_tA);
            }

        }
    }
}



Int_t pp_correlation_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                   Int_t ParticleA, Int_t ParticleB, Int_t Ana_Num,Int_t run_events,Int_t event_number)
{
    //cout << "Hadron_v2_analysis started" << endl;
    event_A_ana   = picoDst_A;

    Int_t n_pp = 0;
    if(ParticleA == 14 && ParticleB == 14) n_pp = 0;
    if(ParticleA == 15 && ParticleB == 15) n_pp = 1;
    if(ParticleA == 14 && ParticleB == 15) n_pp = 2;

    //cout << "ParticleA = " << ParticleA << ", ParticleB = " << ParticleB << ", n_pp = " << n_pp << ", n = " << PID_counter_Array[Ana_Num][ParticleA] << endl;

    // Event vertex information
    Float_t EventVertexX  = event_A_ana->primaryVertex().x();
    Float_t EventVertexY  = event_A_ana->primaryVertex().y();
    Float_t EventVertexZ  = event_A_ana->primaryVertex().z();
    Int_t   refMult       = event_A_ana->refMult();
    Int_t   RunId         = event_A_ana->runId();
    Int_t   mult_bin      = refmultCorrUtil->getCentralityBin9();
    StThreeVectorF vector_prim, vector_prim_new;
    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    if(
       PID_counter_Array[Ana_Num][ParticleA]    > 0
       && PID_counter_Array[Ana_Num][ParticleB] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && mult_bin >= 0 && mult_bin < 9
       && event_A_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " particles *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB;
        StThreeVectorF vectorA, dirA, dirB;

        Double_t pp_corr_array[N_refmult_bins][n_pt_bins_pp_corr][n_pt_bins_pp_corr][2][n_pp_harm]; // last one is for weights [w/o,with]
        Double_t pp_corr_array_counter[N_refmult_bins][n_pt_bins_pp_corr][n_pt_bins_pp_corr][2];

        for(Int_t n_mult = 0; n_mult < N_refmult_bins; n_mult++)
        {
            for(Int_t n_ptA = 0; n_ptA < n_pt_bins_pp_corr; n_ptA++)
            {
                for(Int_t n_ptB = 0; n_ptB < n_pt_bins_pp_corr; n_ptB++)
                {
                    for(Int_t n_w = 0; n_w < 2; n_w++)
                    {
                        pp_corr_array_counter[n_mult][n_ptA][n_ptB][n_w] = 0.0;
                        for(Int_t n_harm = 0; n_harm < n_pp_harm; n_harm++)
                        {
                            pp_corr_array[n_mult][n_ptA][n_ptB][n_w][n_harm] = 0.0;
                        }
                    }
                }
            }
        }

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // particle loop A
        {
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t BetaA       = trackA.btofBeta();
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t nSigmaPA    = trackA.nSigmaProton();

            Float_t Mass2A      = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }

            if(
               nHitsFitA     >= 15  // 15
               && nHitsPossA > 0
               && MomentumA  > 0.1  // 0.15
               && MomentumA  < 10.0
               && nSigmaPA   < 2.5
               && nSigmaPA   > -2.5
              )
            {
                if(
                   ((Float_t)nHitsFitA)/((Float_t)nHitsPossA)    > 0.52
                  )
                {
                    Float_t pathA = -999.0;
                    Float_t dcaAB = -999.0;
                    fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);

                    dirA = helixA.cat(pathA); // direction vector

                    Float_t etaA = dirA.pseudoRapidity();

                    if(
                       dcaA          < 2.0 // 2.0
                       && fabs(etaA) < 1.0 // 1.0
                      )
                    {

                        Float_t p_xA = MomentumA*dirA.x();
                        Float_t p_yA = MomentumA*dirA.y();
                        Float_t p_tA = sqrt(p_xA*p_xA + p_yA*p_yA);

                        Int_t pt_binA = -1;
                        for(Int_t ipt = 0; ipt < n_pt_bins_pp_corr; ipt++)
                        {
                            if(p_tA > pp_corr_pt_bins[0][ipt] && p_tA <= pp_corr_pt_bins[1][ipt])
                            {
                                pt_binA = ipt;
                                break;
                            }
                        }


                        if(
                           pt_binA >= 0 && pt_binA < n_pt_bins_pp_corr
                           &&
                           (
                            (p_tA < 2.6 && Mass2A > 0.6) ||
                            (p_tA >= 2.6 && p_tA < 3.2 && Mass2A > 0.7)
                           )
                           &&
                           (
                            (p_tA < 1.6 && Mass2A < 1.1) ||
                            (p_tA >= 1.6 && p_tA < 2.4 && Mass2A < 1.3) ||
                            (p_tA >= 2.4 && p_tA < 3.2 && Mass2A < 1.5)
                           )
                          )
                        {

                            Double_t phiA   = dirA.phi();
                            Double_t thetaA = dirA.theta();

                            Double_t phi_wA;
                            Double_t total_weightA = calc_event_plane_weight(phiA,p_tA,etaA,RunId,EventVertexZ,(Int_t)PolarityA,phi_wA);

                            Int_t start_track = 0;
                            if(ParticleA == ParticleB) start_track = i+1;

                            for(Int_t j = start_track; j < PID_counter_Array[Ana_Num][ParticleB]; j++) // particle loop B
                            {
                                Int_t trackB_num = PID_Array[Ana_Num][ParticleB][j];
                                if(trackA_num != trackB_num)
                                {
                                    // Get the tracks and calculate the direction and base vectors
                                    StPicoAlexTrack trackB = *picoDst_A->track( trackB_num );
                                    helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackB.charge());

                                    Float_t dcaB        = trackB.dca();
                                    Float_t nHitsPossB  = trackB.nHitsMax();
                                    Float_t nHitsFitB   = trackB.nHitsFit();
                                    Float_t MomentumB   = trackB.gMom().mag();
                                    Float_t BetaB       = trackB.btofBeta();
                                    Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                    Float_t nSigmaPB    = trackB.nSigmaProton();

                                    Float_t Mass2B      = -100.0;
                                    // calculate mass2
                                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                                    {
                                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                                    }

                                    if(
                                       nHitsFitB     >= 15  // 15
                                       && nHitsPossB > 0
                                       && MomentumB  > 0.1  // 0.15
                                       && MomentumB  < 10.0
                                       && nSigmaPB   < 2.5
                                       && nSigmaPB   > -2.5
                                      )
                                    {
                                        if(
                                           ((Float_t)nHitsFitB)/((Float_t)nHitsPossB)    > 0.52
                                          )
                                        {
                                            Float_t pathB = -999.0;
                                            Float_t dcaAB = -999.0;
                                            fHelixAtoPointdca(vector_prim,helixB,pathB,dcaAB);

                                            dirB = helixB.cat(pathB); // direction vector

                                            Float_t etaB = dirB.pseudoRapidity();

                                            if(
                                               dcaB          < 2.0 // 2.0
                                               && fabs(etaB) < 1.0 // 1.0
                                              )
                                            {

                                                Float_t p_xB = MomentumB*dirB.x();
                                                Float_t p_yB = MomentumB*dirB.y();
                                                Float_t p_tB = sqrt(p_xB*p_xB + p_yB*p_yB);

                                                Int_t pt_binB = -1;
                                                for(Int_t ipt = 0; ipt < n_pt_bins_pp_corr; ipt++)
                                                {
                                                    if(p_tB > pp_corr_pt_bins[0][ipt] && p_tB <= pp_corr_pt_bins[1][ipt])
                                                    {
                                                        pt_binB = ipt;
                                                        break;
                                                    }
                                                }

                                                if(
                                                   pt_binB >= 0 && pt_binB < n_pt_bins_pp_corr
                                                   &&
                                                   (
                                                    (p_tB < 2.6 && Mass2B > 0.6) ||
                                                    (p_tB >= 2.6 && p_tB < 3.2 && Mass2B > 0.7)
                                                   )
                                                   &&
                                                   (
                                                    (p_tB < 1.6 && Mass2B < 1.1) ||
                                                    (p_tB >= 1.6 && p_tB < 2.4 && Mass2B < 1.3) ||
                                                    (p_tB >= 2.4 && p_tB < 3.2 && Mass2B < 1.5)
                                                   )
                                                  )
                                                {

                                                    Double_t phiB   = dirB.phi();

                                                    Double_t phi_wB;
                                                    Double_t total_weightB = calc_event_plane_weight(phiB,p_tB,etaB,RunId,EventVertexZ,(Int_t)PolarityB,phi_wB);

                                                    for(Int_t n_harm = 0; n_harm < n_pp_harm; n_harm++)
                                                    {
                                                        pp_corr_array[mult_bin][pt_binA][pt_binB][0][n_harm] += TMath::Cos(((Double_t)(n_harm+1))*(phiA-phiB));
                                                        pp_corr_array[mult_bin][pt_binA][pt_binB][1][n_harm] += TMath::Cos(((Double_t)(n_harm+1))*(phiA-phiB))*phi_wA*phi_wB;
                                                        //if(n_harm == 1 && ParticleA == 15 && ParticleB == 15) cout << "40 ParticleA = " << ParticleA << ", ParticleB = " << ParticleB << ", i = " << i << ", j = " << j
                                                        //    << ", phiA = " << phiA << ", phiB = " << phiB
                                                        //        << ", 2*delta_phi = " << ((Double_t)(n_harm+1))*(phiA-phiB) << ", cos = " << TMath::Cos(((Double_t)(n_harm+1))*(phiA-phiB))
                                                        //        << ", ptA = " << p_tA << ", ptB = " << p_tB << ", pt_binA = " << pt_binA << ", pt_binB = " << pt_binB
                                                        //        << endl;

                                                        if(n_harm == 1 && ParticleA == 14 && ParticleB == 14 && mult_bin == 4)
                                                        {
                                                            if(pt_binA == 0 && pt_binB == 4)
                                                            {
                                                                h_phi_pp_corr[0][0] ->Fill(phiA);
                                                                h_phi_pp_corr[0][1] ->Fill(phiB);
                                                                h_phi_pp_corr[0][2] ->Fill(phiA-phiB);
                                                                h_cos_phi_pp_corr[0][2] ->Fill(TMath::Cos(((Double_t)(n_harm+1))*(phiA-phiB)));
                                                            }
                                                            if(pt_binA == 4 && pt_binB == 0)
                                                            {
                                                                h_phi_pp_corr[1][0] ->Fill(phiA);
                                                                h_phi_pp_corr[1][1] ->Fill(phiB);
                                                                h_phi_pp_corr[1][2] ->Fill(phiA-phiB);
                                                                h_cos_phi_pp_corr[1][2] ->Fill(TMath::Cos(((Double_t)(n_harm+1))*(phiA-phiB)));
                                                            }
                                                        }

                                                        //if(n_harm == 1 && ParticleA == 14 && ParticleB == 14 && pt_binA == 0 && pt_binB == 4) cout << "04 ParticleA = " << ParticleA << ", ParticleB = " << ParticleB << ", i = " << i << ", j = " << j
                                                        //    << ", phiA = " << phiA << ", phiB = " << phiB
                                                        //    << ", 2*delta_phi = " << ((Double_t)(n_harm+1))*(phiA-phiB) << ", cos = " << TMath::Cos(((Double_t)(n_harm+1))*(phiA-phiB)) << endl;
                                                    }

                                                    pp_corr_array_counter[mult_bin][pt_binA][pt_binB][0] = pp_corr_array_counter[mult_bin][pt_binA][pt_binB][0]+1.0;
                                                    pp_corr_array_counter[mult_bin][pt_binA][pt_binB][1] = pp_corr_array_counter[mult_bin][pt_binA][pt_binB][1]+phi_wA*phi_wB;
                                                }
                                            }

                                        }
                                    }
                                }
                            } // end particle loop B


                        }


                    }

                }
            }
        } // end particle loop A

        //cout << "" << endl;
        //cout << "Calculate bin entries" << endl;
        for(Int_t n_mult = 0; n_mult < N_refmult_bins; n_mult++)
        {
            for(Int_t n_ptA = 0; n_ptA < n_pt_bins_pp_corr; n_ptA++)
            {
                for(Int_t n_ptB = 0; n_ptB < n_pt_bins_pp_corr; n_ptB++)
                {
                    for(Int_t n_w = 0; n_w < 2; n_w++)
                    {
                        for(Int_t n_harm = 0; n_harm < n_pp_harm; n_harm++)
                        {
                            Double_t pp_corr_val   = pp_corr_array[n_mult][n_ptA][n_ptB][n_w][n_harm];
                            Double_t pp_weight_val = pp_corr_array_counter[n_mult][n_ptA][n_ptB][n_w];

                            if(pp_weight_val > 0.0)
                            {
                                Double_t pp_corr_norm_val = pp_corr_val/pp_weight_val; // event-by-event average of correlation
                                p_pp_corr[n_mult][n_w][n_harm][n_pp] ->Fill(n_ptA,n_ptB,pp_corr_norm_val);
                                //if(n_harm == 1 && ParticleA == 15 && ParticleB == 15 && n_w == 0) cout << "ptA = " << n_ptA << ", n_ptB = " << n_ptB
                                //    << ", val = " << pp_corr_val << ", weight = " << pp_weight_val << ", val_norm = " << pp_corr_norm_val << endl;
                            }
                        }
                    }
                }
            }
        }

        return 1;
    }
    else return 0;
}



Int_t Event_Display(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                    Int_t ParticleA, Int_t Ana_Num,Int_t run_events,Int_t event_number)
{
    event_A_ana   = picoDst_A;

    // Event vertex information
    Float_t EventVertexX  = event_A_ana->primaryVertex().x();
    Float_t EventVertexY  = event_A_ana->primaryVertex().y();
    Float_t EventVertexZ  = event_A_ana->primaryVertex().z();

    //x_mean = EventVertexX;
    //y_mean = EventVertexY;
    //z_mean = EventVertexZ;

    //cout << "Vertex = {" << EventVertexX << ", " << EventVertexY << ", " << EventVertexZ << "}" << endl;
    StThreeVectorF vector_prim_old, vector_prim_new,vector_prim_new_secondary;
    vector_prim_old.set(EventVertexX,EventVertexY,EventVertexZ);
    vector_prim_new.set(x_mean,y_mean,z_mean);
    vector_prim_new_secondary.set(x_mean2,y_mean2,z_mean2);
    //Double_t vector_dist  = sqrt((x_mean-x_mean2)*(x_mean-x_mean2)+(y_mean-y_mean2)*(y_mean-y_mean2)+(z_mean-z_mean2)*(z_mean-z_mean2));
    //Double_t vector_dist2 = sqrt((x_mean-EventVertexX)*(x_mean-EventVertexX)+(y_mean-EventVertexY)*(y_mean-EventVertexY)+(z_mean-EventVertexZ)*(z_mean-EventVertexZ));
    Int_t   refMult       = event_A_ana->refMult();
    Int_t   RunId         = event_A_ana->runId();
    //cout << "Event display started, number of tracks = " << PID_counter_Array[Ana_Num][ParticleA] << ", RunId = " << RunId << ", refMult = " << refMult << endl;
    //total_event++;
    //

   // cout << "|v1-v2| = " << vector_dist_prim_second << ", N_prim = " << N_primaries_new << ", N_second = " << N_primaries_new2
   //     << ", N_prim_tof = " << N_tofmatch_new << ", N_second_tof = " << N_tofmatch_new2 << endl;

    if(
       1
       //&& z_mean > 208.5 && z_mean < 210.5
       //&& x_mean > -3.0 && x_mean < 2.2
       //&& y_mean > -4.0 && y_mean < -3.0
       //&& comb_counter_global > 50

       //PID_counter_Array[Ana_Num][ParticleA] > 2000

       //&& refMult > 300
       //&& event_A_ana   ->isMinBias()

       //&& vector_dist_prim_second > 20.0
       //&& N_primaries_new > 30
       //&& N_primaries_new2 > 3
       //&& N_tofmatch_new < N_primaries_new*0.2
       //&& N_tofmatch_new2 == 0

       //&& (x_mean*x_mean + y_mean*y_mean) < 4.0
       //&& fabs(z_mean) < 70.0

       && x_mean > -3.0 && x_mean < 2.2
       && y_mean > -4.0 && y_mean < -3.0
       && z_mean > 208.5 && z_mean < 210.5
       //&& refMult > 10
       //&& N_primaries_new2 > 9
       //&& fabs(tracks_left_to_right_diff) > 0.2
       //&& tracks_left_to_right_ratio > 1.2
       //&& non_prim_to_prim_ratio > 0.7
       //&& vector_dist  > 10.0
       //&& vector_dist2 < 1.0
       //&& x_mean2 > -200.0
       //&& event_A_ana->getNVertices() == 1
       //&& (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < 2.0*2.0
       //&& fabs(EventVertexZ) < 70.0
       //&& (event_A_ana->getTriggerId() == 290001 || event_A_ana->getTriggerId() == 290004)
       //&& goodHLT == 1

      )
    {
        n_poly_track_counter = 0;
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " particles *************" << endl;
        //cout << "NVertices = " << event_A_ana->getNVertices() << ", Vertex = {" << EventVertexX << ", " << EventVertexY << ", " << EventVertexZ << "}" << endl;
        // Loop over all particle combinations


        StPhysicalHelixD helixA;
        StThreeVectorF vectorA, vectoratsA, vectorAB, vector_primAB, vectornewA, vectornewB, vector_sum, dirA, vector_diff;
        StThreeVectorF testA, testAB, vectorABtoPrim, baseY, dirY;

        hEventDCA_new->Reset();
        hEventDCA_old->Reset();

        Int_t flag_high_tower = 0;
        // Pre loop
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            //cout << "i = " << i << endl;
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );

            Float_t BEMC_e0 = trackA.e0();
            if(BEMC_e0 > 10.0)
            {
                flag_high_tower = 1;
                break;
            }
        }

        //if(flag_high_tower)
        {
            cout << "High tower event with refMult = " << refMult << endl;
            for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
            {
                //cout << "i = " << i << endl;
                // Get the tracks and calculate the direction and base vectors
                Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
                StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
                helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());
                //Float_t dcaA        = trackA.dca();
                //Float_t nHitsPossA  = trackA.nHitsMax();
                //Float_t nHitsFitA   = trackA.nHitsFit();
                Float_t MomentumA   = trackA.gMom().mag();
                Float_t Beta        = trackA.btofBeta();

                Float_t nSigmaPion     = trackA.nSigmaPion();
                Float_t nSigmaKaon     = trackA.nSigmaKaon();
                Float_t nSigmaProton   = trackA.nSigmaProton();
                Float_t nSigmaElectron = trackA.nSigmaElectron();

                Float_t Mass2A        = -100.0;
                // calculate mass2
                if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && Beta != 0)
                {
                    Mass2A = MomentumA*MomentumA*(1.0/(Beta*Beta) - 1.0);
                }

                StThreeVectorF tofhitpos =  trackA.btofHisPos();

                Float_t TofHitX     = tofhitpos.x();
                Float_t TofHitY     = tofhitpos.y();
                Float_t TofHitZ     = tofhitpos.z();

                vector_diff = helixA.at(10.0)-helixA.at(0.0);
                //Double_t vector_diff_z = vector_diff.z();

                EDisplay_NTDataArray[0]      =(Float_t)x_mean;         //
                EDisplay_NTDataArray[1]      =(Float_t)y_mean;         //
                EDisplay_NTDataArray[2]      =(Float_t)z_mean;         //
                EDisplay_NTDataArray[3]      =(Float_t)run_events;
                EDisplay_NTDataArray[4]      =(Float_t)EventVertexX;
                EDisplay_NTDataArray[5]      =(Float_t)EventVertexY;
                EDisplay_NTDataArray[6]      =(Float_t)EventVertexZ;
                EDisplay_NTDataArray[7]      =(Float_t)x_mean2;         //
                EDisplay_NTDataArray[8]      =(Float_t)y_mean2;         //
                EDisplay_NTDataArray[9]      =(Float_t)z_mean2;         //
                EDisplay_NTDataArray[10]     =(Float_t)refMult;
                EDisplay_NTDataArray[11]     =(Float_t)N_tofmatch_new;
                EDisplay_NTDataArray[12]     =(Float_t)N_tofmatch_new2;
                EDisplay_NTDataArray[13]     =(Float_t)vector_dist_prim_second;
                EDisplay_NTDataArray[14]     =(Float_t)RunId;
                EDisplay_NTDataArray[15]     =(Float_t)PID_counter_Array[Ana_Num][ParticleA];         //
                EDisplay_NTDataArray[16]     =(Float_t)trackA.gMom().x();         //
                EDisplay_NTDataArray[17]     =(Float_t)trackA.gMom().y();         //
                EDisplay_NTDataArray[18]     =(Float_t)trackA.gMom().z();         //
                EDisplay_NTDataArray[19]     =(Float_t)trackA.origin().x();         //
                EDisplay_NTDataArray[20]     =(Float_t)trackA.origin().y();         //
                EDisplay_NTDataArray[21]     =(Float_t)trackA.origin().z();         //
                EDisplay_NTDataArray[22]     =(Float_t)event_A_ana->bField()*MAGFIELDFACTOR;         //
                EDisplay_NTDataArray[23]     =(Float_t)trackA.charge();         //
                EDisplay_NTDataArray[24]     =(Float_t)trackA.dca();         //
                EDisplay_NTDataArray[25]     =(Float_t)trackA.btofBeta();         //
                EDisplay_NTDataArray[26]     =(Float_t)tofhitpos.x();         //
                EDisplay_NTDataArray[27]     =(Float_t)tofhitpos.y();         //
                EDisplay_NTDataArray[28]     =(Float_t)tofhitpos.z();         //
                EDisplay_NTDataArray[29]     =(Float_t)ED_NT_counter;         //
                EDisplay_NTDataArray[30]     =(Float_t)trackA.e1();         //
                EDisplay_NTDataArray[31]     =(Float_t)trackA.e2();         //
                EDisplay_NTDataArray[32]     =(Float_t)trackA.e3();         //
                EDisplay_NTDataArray[33]     =(Float_t)trackA.e();         //
                EDisplay_NTDataArray[34]     =(Float_t)trackA.e0();         //
                EDisplay_NTDataArray[35]     =(Float_t)trackA.btowId();         //
                EDisplay_NTDataArray[36]     =(Float_t)Mass2A;         //
                EDisplay_NTDataArray[37]     =(Float_t)nSigmaPion;         //
                EDisplay_NTDataArray[38]     =(Float_t)nSigmaKaon;         //
                EDisplay_NTDataArray[39]     =(Float_t)nSigmaProton;         //
                EDisplay_NTDataArray[40]     =(Float_t)nSigmaElectron;         //



                if((fAnalysisNum == 9))
                {
                    EDisplay_NT->Fill(EDisplay_NTDataArray);
                }


                //if(EDisplay_NTDataArray[34] > 0.0) cout << "Good tower hit" << endl;




                if(
                   1 == 2  &&
                   MomentumA > 0.15
                   //&& dcaA > 3.0
                  )
                {

                    testA     = helixA.at(helixA.pathLength(vector_prim_new)); // 3D-vector of helixA point at dca to testB
                    //Float_t dca_calc = (testA-vector_prim_new).mag(); // dca between helixA and helixB
                    Float_t pathA = -999.0;
                    Float_t dcaAB = -999.0;
                    Float_t pathA_secondary = -999.0;
                    Float_t dcaAB_secondary = -999.0;
                    fHelixAtoPointdca(vector_prim_old,helixA,pathA,dcaAB);
                    hEventDCA_old ->Fill(dcaAB);
                    fHelixAtoPointdca(vector_prim_new,helixA,pathA,dcaAB);
                    hEventDCA_new ->Fill(dcaAB);
                    fHelixAtoPointdca(vector_prim_new_secondary,helixA,pathA_secondary,dcaAB_secondary);
                    //if(dca_calc < dcaA)
                    //{
                    //cout << "dca = " << dcaA << ", dca_calc = " << dca_calc << ", dcaAB (Alex) = " << dcaAB
                    //    << ", path helix = " << helixA.pathLength(vector_prim) << ", path (Alex) = " << pathA << endl;
                    //}

                    Int_t n_max_points = 100;
                    Float_t step_size  = 5.0; // path length [cm]
                    Float_t z_max      = 200.0;
                    Float_t z_min      = -200.0;
                    Float_t radius_max = 200.0;

                    pTofHit[n_poly_track_counter] ->Delete();
                    pTofHit[n_poly_track_counter] = new TMarker3DBox();
                    pTofHit[n_poly_track_counter] ->SetSize(1.5,1.5,1.5);
                    pTofHit[n_poly_track_counter] ->SetLineColor(4);
                    pTofHit[n_poly_track_counter] ->SetLineStyle(1);
                    pTofHit[n_poly_track_counter] ->SetLineWidth(1);
                    pTofHit[n_poly_track_counter] ->SetPosition(TofHitX,TofHitY,TofHitZ);

                    pTrack[n_poly_track_counter]   ->Delete();
                    pTrack[n_poly_track_counter]   = new TPolyLine3D();
                    pTrack[n_poly_track_counter]   ->SetLineWidth(1);
                    pTrack[n_poly_track_counter]   ->SetLineColor(2);
                    pTrack[n_poly_track_counter]   ->SetLineStyle(1);
                    if(
                       (Beta > -300.0)
                      )
                    {
                        pTrack[n_poly_track_counter]   ->SetLineColor(2);
                    }
                    //if(
                    //   (x_mean*x_mean + y_mean*y_mean) < 20.0*20.0
                    //   && fabs(z_mean) < 70.0
                    //  )
                    //{
                    //    cout << "Beta = " << Beta << ", MomentumA = " << MomentumA << endl;
                    //}

                    pTrack_extrapolate[n_poly_track_counter]   ->Delete();
                    pTrack_extrapolate[n_poly_track_counter]   = new TPolyLine3D();
                    pTrack_extrapolate[n_poly_track_counter]   ->SetLineWidth(1);
                    pTrack_extrapolate[n_poly_track_counter]   ->SetLineColor(5);
                    pTrack_extrapolate[n_poly_track_counter]   ->SetLineStyle(1);

                    //if(vector_diff_z > 0.0) pTrack_extrapolate[n_poly_track_counter]   ->SetLineColor(3);

                    pTrack_extrapolate_to_Tof[n_poly_track_counter]   ->Delete();
                    pTrack_extrapolate_to_Tof[n_poly_track_counter]   = new TPolyLine3D();
                    pTrack_extrapolate_to_Tof[n_poly_track_counter]   ->SetLineWidth(1);
                    pTrack_extrapolate_to_Tof[n_poly_track_counter]   ->SetLineColor(7);
                    pTrack_extrapolate_to_Tof[n_poly_track_counter]   ->SetLineStyle(1);

                    Float_t path_length_end;

                    // from first measured point to end of TPC
                    for(Int_t j = 0; j < n_max_points; j++)
                    {
                        Float_t path_length = 0.0 + j * step_size;
                        testA     = helixA.at(path_length);
                        Float_t x_val = testA.x();
                        Float_t y_val = testA.y();
                        Float_t z_val = testA.z();
                        if(
                           (x_val*x_val + y_val*y_val) < radius_max*radius_max
                           && z_val < z_max
                           && z_val > z_min
                           && n_poly_track_counter < n_poly_marker_track
                          )
                        {
                            pTrack[n_poly_track_counter] ->SetNextPoint(x_val,y_val,z_val);
                            path_length_end = path_length;
                        }
                        else break;
                    }

                    // from end of TPC tof Tof
                    for(Int_t j = 0; j < n_max_points; j++)
                    {
                        Float_t path_length = path_length_end + j * step_size;
                        testA     = helixA.at(path_length);
                        Float_t x_val = testA.x();
                        Float_t y_val = testA.y();
                        Float_t z_val = testA.z();
                        if(
                           (x_val*x_val + y_val*y_val) < 217.0*217.0
                           //&& ( z_val < z_max && z_val > z_min)
                           && (TofHitX*TofHitX + TofHitY*TofHitY) > 10.0
                           && n_poly_track_counter < n_poly_marker_track
                          )
                        {
                            pTrack_extrapolate_to_Tof[n_poly_track_counter] ->SetNextPoint(x_val,y_val,z_val);
                        }
                        else break;
                    }

                    // from dca to primary vertex to first measured point
                    if(dcaAB_secondary < dcaAB) pathA = pathA_secondary;
                    for(Int_t j = 0; j < n_max_points; j++)
                    {
                        Float_t path_length = pathA + j * step_size;
                        testA     = helixA.at(path_length);
                        Float_t x_val = testA.x();
                        Float_t y_val = testA.y();
                        Float_t z_val = testA.z();
                        if(
                           (x_val*x_val + y_val*y_val) < radius_max*radius_max
                           && z_val < z_max
                           && z_val > z_min
                           && n_poly_track_counter < n_poly_marker_track
                           && path_length < 0.0
                          )
                        {
                            pTrack_extrapolate[n_poly_track_counter] ->SetNextPoint(x_val,y_val,z_val);
                        }
                        else break;
                    }

                    n_poly_track_counter++;
                }
            }

            cout << "Event counter = " << ED_NT_counter << ", event number = "  << event_number
                << ", n_poly_track_counter = " << n_poly_track_counter
                << ", new vertex = {" << x_mean << ", " << y_mean << ", " << z_mean << "}" << endl;

            event_counter++;
            ED_NT_counter++;
            //****************************** Event display *********************************************************
            if(
               1 == 2 &&
               ED_counter < n_cEventDisplay_array
              )
            {
                cout << "Event display counter = " << ED_counter << ", event number = "  << event_number
                    << ", n_poly_track_counter = " << n_poly_track_counter
                    << ", new vertex = {" << x_mean << ", " << y_mean << ", " << z_mean << "}" << endl;
                pPrimaryVertex->SetPosition(x_mean,y_mean,z_mean);
                pPrimaryVertex2->SetPosition(x_mean2,y_mean2,z_mean2);
                pPrimaryVertex_old->SetPosition(EventVertexX,EventVertexY,EventVertexZ);
                cEventDisplay_array[ED_counter]->cd(0);
                for(Int_t r = 0; r < 4; r++)
                {
                    TPC_endcaps[r] ->DrawClone();
                }
                for(Int_t g = 0; g < n_poly_track_counter; g++)
                {
                    pTrack[g]                    -> DrawClone();
                    pTrack_extrapolate[g]        -> DrawClone();
                    pTrack_extrapolate_to_Tof[g] -> DrawClone();
                    pTofHit[g]                   -> DrawClone();
                }
                BeamLine        ->DrawClone();
                XAxis           ->DrawClone();
                YAxis           ->DrawClone();
                ZAxis           ->DrawClone();
                pPrimaryVertex  ->DrawClone();
                if(z_mean2 > -900) pPrimaryVertex2 ->DrawClone();
                pPrimaryVertex_old ->DrawClone();
                cEventDisplay_array[ED_counter]->Update();

                cEventDCA_array[ED_counter]->cd(0);
                hEventDCA_new->SetLineColor(1);
                hEventDCA_old->SetLineColor(2);
                hEventDCA_old->DrawCopy("h");
                hEventDCA_new->DrawCopy("same h");

                ED_counter++;
            }
            //******************************************************************************************************
        }


        return 1;
    }
    else return 0;
}



Int_t Vertex_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                      Int_t ParticleA, Int_t ParticleB, Int_t Ana_Num,Int_t run_events,Int_t event_number,
                      Float_t& x_mean, Float_t& y_mean, Float_t& z_mean, Int_t Save_flag)
{
    // V1.2 09/19/2011 changed to new pico structure
    // V1.1
    event_A_ana   = picoDst_A;

    //********************************************************************************************
    //Vertexer parameters
    Int_t mix_tracks_min = 10;    // minimum number of cross combinations
    Float_t min_nhitsfit = 8;   // minimum number of TPC fit points for each track
    Float_t min_momentum = 0.1;  // minimum momentum for each track
    //********************************************************************************************

    Int_t mix_tracks     = 20;   // number of cross combinations -> adaptive

    N_tofmatch_new        = 0;
    N_tofmatch_new2       = 0;
    N_tofmatch_old        = 0;
    N_primaries_old       = 0;
    N_primaries_new       = 0;
    N_primaries_new2      = 0;

    // Event vertex information
    Float_t EventVertexX  = event_A_ana->primaryVertex().x();
    Float_t EventVertexY  = event_A_ana->primaryVertex().y();
    Float_t EventVertexZ  = event_A_ana->primaryVertex().z();

    Int_t   refMult       = event_A_ana->refMult();
    Int_t   RunId         = event_A_ana->runId();
    Int_t   trigger_word  = event_A_ana->triggerWord();
    Int_t   EventId       = event_A_ana->eventId();

    Float_t BBC           = 0.0;
    Float_t VPD           = 0.0;

    //Int_t   flag_counter  = 0; // counts the number of tracks which belong to the primary vertex

    hEvent_MeanX->Reset();
    hEvent_MeanY->Reset();
    hEvent_MeanZ->Reset();


    total_event++;
    //

    if(
       1 == 1
       // && event_A_ana   ->isMinBias()
       //refMult > 50
       //&& refMult <= 50
       //&& tracks_left_to_right_ratio > 1.2
       //&& fabs(tracks_left_to_right_diff) > 0.2
       //&& non_prim_to_prim_ratio > 0.7
       //PID_counter_Array[Ana_Num][ParticleA] > 2000  //
       //(EventVertexX*EventVertexX + EventVertexY*EventVertexY) < 2.0*2.0
       //&& fabs(EventVertexZ) < 70.0
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " particles, runId = " << RunId << "  *************" << endl;
        // Loop over all particle combinations
        StPhysicalHelixD helixA, helixB;
        StThreeVectorD vectorA, vectorB, vectoratsA, vectoratsB, vectoratsA2, vectorAB, vector_primAB, vectornewA, vectornewB, vector_sum;
        StThreeVectorD testA, testB, testAB, vectorABtoPrim, baseY, dirY;
        vector_sum.set(0.0,0.0,0.0);
        dirY.set(0.0,0.0,1.0);
        baseY.set(0.0,0.0,0.0);
        Int_t comb_counter = 0;
        const Int_t max_tracks = PID_counter_Array[Ana_Num][ParticleA];
        

        //*********************************************************************************
        // set mix tracks
        Int_t max_mix_tracks = max_tracks - 1; // Maximum number of tracks which can be combined with each track without double counting
        if(max_tracks < 40)                     mix_tracks = max_tracks;
        //if(max_tracks >= 20 && max_tracks < 40) mix_tracks = (Int_t)(100.0/((Float_t)max_tracks));
        if(max_tracks >= 40)                    mix_tracks = mix_tracks_min;

        if(mix_tracks > max_mix_tracks) mix_tracks = max_mix_tracks;
        //*********************************************************************************

        const Int_t N_xyz_table = (Int_t)(max_tracks*mix_tracks);
        Float_t x_table[N_xyz_table][2]; // [value, weight]
        Float_t y_table[N_xyz_table][2];
        Float_t z_table[N_xyz_table][2];
        Int_t   track_num_array[max_tracks]; // this array is used for the random track combinations
        Int_t   track_ran_num_array[max_tracks]; // randomized order of the tracks

        // Tracks are ordered in azimuthal space --> randomize track number for later combinations
        TRandom ran_number;
        //Int_t seed_numberB = seed_number + max_tracks; // changes from event to event
        ran_number.SetSeed((UInt_t)seed_number);

        for(Int_t m = 0; m < max_tracks; m++)
        {
            track_ran_num_array[m] = m;
        }

        for(Int_t m = 0; m < max_tracks; m++)
        {
            Int_t track_ran_position = ran_number.Integer(max_tracks-m);
            Int_t save_num = track_ran_num_array[track_ran_position];
            track_ran_num_array[track_ran_position] = track_ran_num_array[max_tracks-m-1];
            track_ran_num_array[max_tracks-m-1] = save_num;
        }

        const Int_t tracks_max_mix = max_tracks + mix_tracks;
        Int_t track_num_array_max_mix[tracks_max_mix];
        for(Int_t m = 0; m < max_tracks; m++)
        {
            track_num_array[m] = track_ran_num_array[m]; // first loop: direct mapping
            track_num_array_max_mix[m] = track_num_array[m]; // the first block of the array is identical to track_num_array
        }
        for(Int_t m = 0; m < mix_tracks; m++)
        {
            track_num_array_max_mix[max_tracks+m] = track_num_array[m]; // the additional mix_tracks block is a copy of the first entries from track_num_array
        }

        //if(max_tracks < 40) cout << "Start vertex analysis, refMult = " << refMult << ", max_tracks = " << max_tracks << endl;


        //cout << "Start vertex analysis: ED_counter = " << ED_counter << ", max_tracks = " << max_tracks << ", mix_tracks = " << mix_tracks << endl;

        // Get the tracks and calculate the direction and base vectors
        for(Int_t i = 0; i < max_tracks; i++) //
        {
            Int_t track_Ai         = track_num_array[i];
            Int_t trackA_num       = PID_Array[Ana_Num][ParticleA][track_Ai];

            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            //Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            //Float_t BetaA       = trackA.btofBeta();
            //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            //Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
            //Float_t nSigmaKaonA = trackA.nSigmaKaon();
            //Float_t nSigmaPionA = trackA.nSigmaPion();
            //Float_t nSigmaPA    = trackA.nSigmaProton();


            //cout << "i = " << i << ", max_tracks = "  << max_tracks << ", nHitsFitA = " << nHitsFitA << ", nHitsPossA = "
            //    << nHitsPossA << ", MomentumA = " << MomentumA << endl;

            if(
               nHitsFitA     > min_nhitsfit
               && nHitsPossA > 0
               && MomentumA  > min_momentum
              )
            {
                helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

                //for(Int_t m = 0; m < max_tracks; m++)
                //{
                //    track_num_array[m] = track_ran_num_array[m]; // first loop: direct mapping
                //}

                for(Int_t j = 0; j < mix_tracks; j++)  // loop over the cross combination
                {
                    //Int_t track_ran_number = 0;
                    //if(max_tracks-i-j-1 >= 0) track_ran_number = ran_number.Integer(max_tracks-i-j-1);  // will generate a random integer from 0..(max_tracks-2-i-j)
                    //if(max_tracks < 40) cout << "i = " << i << ", j = " << j  << ", max_tracks = "  << max_tracks << ", max_tracks-i-j-1 = " << max_tracks-i-j-1 << ", track_ran_number = " << track_ran_number << endl;
                    //track_ran_number += i+j+1; // exclude the first block of the track_num_array which was already used
                    //cout << "track_ran_numberB = " << track_ran_number << endl;

                    if(i+j+1 < max_tracks+mix_tracks)
                    {

                        Int_t track_Bi = track_num_array_max_mix[i+j+1]; // has a length of (max_tracks + mix_tracks)
                        //cout << "track_ran_number_array = " << track_ran_number_array << endl;

                        //track_num_array[track_ran_number] = i+j+1; // track_ran_number was used, this entry is now exchanged by i+j+1

                        //for(Int_t m = 0; m < max_tracks; m++)
                        //{
                        //    cout << "m = " << m << ", track_num_array = " <<  track_num_array[m] << endl;
                        //}

                        Int_t trackB_num = PID_Array[Ana_Num][ParticleB][track_Bi];

                        if(
                           trackA_num != trackB_num
                          )
                        {
                            //if(max_tracks < 40) cout << "trackA_num = " << trackA_num << ", trackB_num = " << trackB_num << endl;

                            StPicoAlexTrack trackB = *picoDst_A->track( trackB_num );
                            //Float_t dcaB        = trackB.dca();
                            Float_t nHitsPossB  = trackB.nHitsMax();
                            Float_t nHitsFitB   = trackB.nHitsFit();
                            Float_t MomentumB   = trackB.gMom().mag();
                            //Float_t BetaB       = trackB.btofBeta();
                            //Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                            //Float_t TPCdEdxB    = trackB.dEdx(); // Combined inner and outer MDC dE/dx
                            //Float_t nSigmaKaonB = trackB.nSigmaKaon();
                            //Float_t nSigmaPionB = trackB.nSigmaPion();
                            //Float_t nSigmaPB    = trackB.nSigmaProton();

                            if(
                               nHitsFitB     > min_nhitsfit
                               && nHitsPossB > 0
                               && MomentumB  > min_momentum
                              )
                            {
                                helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackB.charge());

                                if(
                                   1
                                   //((Float_t)nHitsFitA)/((Float_t)nHitsPossA)    > 0.55
                                   //&& ((Float_t)nHitsFitB)/((Float_t)nHitsPossB) > 0.55
                                  )
                                {
                                    //Float_t pathA2_f, pathB2_f,dcaBB2_f,dcaAB2_f;
                                    //fHelixAtoLinedca(dirY,baseY,helixA,pathB2_f,pathA2_f,dcaAB2_f); // calculates dca of helix to line
                                    //fHelixAtoLinedca(dirY,baseY,helixB,pathB2_f,pathA2_f,dcaBB2_f); // calculates dca of helix to line

                                    Float_t pathA_f, pathB_f, dcaAB_f;
                                    Float_t pathA_test = 0.0;
                                    Float_t pathB_test = 0.0;
                                    Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);
                                    Double_t DCA_helixAB_est = 0.0;

                                    vectoratsA     = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
                                    vectoratsB     = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
                                    vectorAB       = vectoratsA+vectoratsB;
                                    DCA_helixAB_est = (vectoratsA-vectoratsB).mag();
                                    vectorAB       = vectorAB/2.0; // decay vertex

                                    //Double_t vex_est_x = vectorAB.x();
                                    //Double_t vex_est_y = vectorAB.y();
                                    //Double_t vex_est_z = vectorAB.z();

                                    if(DCA_helixAB_est < 5.0)
                                    {

                                        if(fDCA_Helix_out == 1)
                                        {
                                            fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
                                        }
                                        else
                                        {
                                            fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
                                        }

                                        vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                                        vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                                        vectorAB       = vectoratsA+vectoratsB;
                                        vectorAB       = vectorAB/2.0; // dca vertex

                                        //cout << "i = " << i << ", j = " << j << ", track_Ai = " << track_Ai << ", track_Bi = " << track_Bi << ", nHitsFitA = "
                                        //    << nHitsFitA << ", nHitsFitB = " << nHitsFitB  << ", pA = " << MomentumA << ", pB = " << MomentumB
                                        //    << ", DCA_helixAB_est = " << DCA_helixAB_est << ", dcaAB = " << dcaAB_f << endl;

                                        if(
                                           dcaAB_f < 3.0
                                           //&& dcaAB2_f < 20.0
                                           //&& dcaBB2_f < 20.0
                                           //&& (vectorAB.x()*vectorAB.x() + vectorAB.y()*vectorAB.y()) < 30.0
                                           //&& fabs(vectorAB.x()) < 60.0
                                           //&& fabs(vectorAB.y()) < 60.0
                                           //&& fabs(vectorAB.z()) < 200.0
                                          )
                                        {
                                            //cout << "i = " << i << ", j = " << j << ", comb_counter = " << comb_counter << ", N_xyz_table = " << N_xyz_table << ", max_tracks = " << max_tracks << ", mix_tracks = " << mix_tracks << endl;
                                            // vertex values
                                            x_table[comb_counter][0] = vectorAB.x();
                                            y_table[comb_counter][0] = vectorAB.y();
                                            z_table[comb_counter][0] = vectorAB.z();

                                            // weights
                                            //x_table[comb_counter][1] = (nHitsFitA/nHitsPossA)*(nHitsFitB/nHitsPossB);
                                            //y_table[comb_counter][1] = (nHitsFitA/nHitsPossA)*(nHitsFitB/nHitsPossB);
                                            //z_table[comb_counter][1] = (nHitsFitA/nHitsPossA)*(nHitsFitB/nHitsPossB);

                                            x_table[comb_counter][1] = 1.0;
                                            y_table[comb_counter][1] = 1.0;
                                            z_table[comb_counter][1] = 1.0;

                                            hx_vertex_distr->Fill(vectorAB.x());
                                            hy_vertex_distr->Fill(vectorAB.y());
                                            hz_vertex_distr->Fill(vectorAB.z());

                                            hEvent_MeanX ->Fill(vectorAB.x());
                                            hEvent_MeanY ->Fill(vectorAB.y());
                                            hEvent_MeanZ ->Fill(vectorAB.z());

                                            comb_counter++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        x_mean  = -999.0;
        y_mean  = -999.0;
        z_mean  = -999.0;
        x_mean2 = -999.0;
        y_mean2 = -999.0;
        z_mean2 = -999.0;

        Float_t tracks_used = (Float_t)run_events;
        if(run_events > max_tracks) tracks_used = (Float_t)max_tracks;
        Float_t max_comb_counter = ((tracks_used*tracks_used)-tracks_used)/2.0;
        Float_t fraction_comb_counter = 100.0*comb_counter/max_comb_counter;

        if(
           comb_counter > 1
           //&& fraction_comb_counter > 20.0
          )  // 40
        {
            x_mean = CalcTruncMean(x_table,1.0,0.3,comb_counter);
            y_mean = CalcTruncMean(y_table,1.0,0.3,comb_counter);
            z_mean = CalcTruncMean(z_table,1.0,0.3,comb_counter);

            if(
               ED_counter < n_cEventDisplay_array
              )
            {
                //cout << "Plot MeanXYZ for ED_counter = " << ED_counter << ", comb_counter = " << comb_counter << ", refMult = " << refMult
                //    << ", Vertex_alex = {" << x_mean << ", " << y_mean << ", " << z_mean << "}" << endl;

                cEvent_MeanX_array[ED_counter]->cd(0);
                hEvent_MeanX->DrawCopy("h");
                PlotLine(x_mean,x_mean,0.0,10,2,2,1); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle

                cEvent_MeanY_array[ED_counter]->cd(0);
                hEvent_MeanY->DrawCopy("h");
                PlotLine(y_mean,y_mean,0.0,10,2,2,1); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle

                cEvent_MeanZ_array[ED_counter]->cd(0);
                hEvent_MeanZ->DrawCopy("h");
                PlotLine(z_mean,z_mean,0.0,10,2,2,1); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle

                ED_counter++;
            }

            // Revome all tracks with a dca < 3 cm to the primary vertex
            Int_t non_primary_tracks = 0; // counts the number of tracks which are NOT close to the primary vertex
            Int_t track_num_array_secondary[max_tracks]; // this array is used for the random track combinations (secondary_vertex)
            StThreeVectorD prim_vertex;
            prim_vertex.set(x_mean,y_mean,z_mean);
            for(Int_t i = 0; i < max_tracks; i++) //
            {
                Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];

                StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
                //Float_t dcaA        = trackA.dca();
                Float_t nHitsPossA  = trackA.nHitsMax();
                Float_t nHitsFitA   = trackA.nHitsFit();
                Float_t MomentumA   = trackA.gMom().mag();
                Float_t BetaA       = trackA.btofBeta();
                //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                //Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
                //Float_t nSigmaKaonA = trackA.nSigmaKaon();
                //Float_t nSigmaPionA = trackA.nSigmaPion();
                //Float_t nSigmaPA    = trackA.nSigmaProton();

                if(
                   nHitsFitA     > min_nhitsfit
                   && nHitsPossA > 0
                   && MomentumA  > min_momentum
                  )
                {
                    helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());

                    Float_t pathA_prim_dca = -999.0;
                    Float_t dcaAB_prim_dca = -999.0;
                    fHelixAtoPointdca(prim_vertex,helixA,pathA_prim_dca,dcaAB_prim_dca); // new helix to point dca calculation

                    if(dcaAB_prim_dca > 3.0) // non primaries
                    {
                        track_num_array_secondary[non_primary_tracks] = i;
                        non_primary_tracks++;
                        N_primaries_new2++;
                        if(BetaA > -500) N_tofmatch_new2++;
                    }
                    else // primary
                    {
                        N_primaries_new++;
                        if(BetaA > -500) N_tofmatch_new++;
                    }
                }
            }

            //cout << "Primary vertex = {" << x_mean << ", " << y_mean << ", " << z_mean << "}" << ", non_primary_tracks = " << non_primary_tracks << endl;

            // Loop for secondary vertex
            const Int_t N_xyz_table2 = (Int_t)(non_primary_tracks*mix_tracks);
            Float_t x_table2[N_xyz_table2][2]; // [value, weight]
            Float_t y_table2[N_xyz_table2][2];
            Float_t z_table2[N_xyz_table2][2];

            Int_t comb_counter2 = 0;

            // Secondary vertex
            for(Int_t i = 0; i < non_primary_tracks; i++) //
            {
                Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
                StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
                //Float_t dcaA        = trackA.dca();
                Float_t nHitsPossA  = trackA.nHitsMax();
                Float_t nHitsFitA   = trackA.nHitsFit();
                Float_t MomentumA   = trackA.gMom().mag();
                //Float_t BetaA       = trackA.btofBeta();
                //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                //Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
                //Float_t nSigmaKaonA = trackA.nSigmaKaon();
                //Float_t nSigmaPionA = trackA.nSigmaPion();
                //Float_t nSigmaPA    = trackA.nSigmaProton();

                if(
                   nHitsFitA     > min_nhitsfit
                   && nHitsPossA > 0
                   && MomentumA  > min_momentum
                  )
                {
                    helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());
                    for(Int_t m = 0; m < non_primary_tracks; m++)
                    {
                        track_num_array[m] = track_num_array_secondary[m]; // first loop: direct mapping m = 0 -> track_num_array_secondary[m],...
                    }

                    for(Int_t j = 0; j < mix_tracks; j++)  //
                    {
                        Int_t track_ran_number = 0;
                        if(non_primary_tracks-i-j-1 >= 0)
                        {
                            track_ran_number = ran_number.Integer(non_primary_tracks-i-j-1);  // will generate a random integer from 0..(max_tracks-2-i-j)

                        }
                        //cout << "i = " << i << ", j = " << j  << ",non_primary_tracks  = "  << non_primary_tracks << ", non_primary_tracks-i-j-1 = " << non_primary_tracks-i-j-1 << ", track_ran_number = " << track_ran_number << endl;
                        track_ran_number += i+j+1; // exclude the first block of the track_num_array which was already used
                        //cout << "track_ran_numberB = " << track_ran_number << endl;

                        if(track_ran_number < non_primary_tracks && (non_primary_tracks-i-j-1) > 0)
                        {
                            Int_t track_ran_number_array = track_num_array[track_ran_number];
                            //cout << "track_ran_number_array = " << track_ran_number_array << endl;
                            track_num_array[track_ran_number] = i+j+1; // track_ran_number was used, this entry is now exchanged by i+j+1

                            //for(Int_t m = 0; m < max_tracks; m++)
                            //{
                            //    cout << "m = " << m << ", track_num_array = " <<  track_num_array[m] << endl;
                            //}
                            Int_t trackB_num = PID_Array[Ana_Num][ParticleB][track_ran_number_array];
                            if(
                               trackA_num != trackB_num
                              )
                            {
                                StPicoAlexTrack trackB = *picoDst_A->track( trackB_num );
                                //Float_t dcaB        = trackB.dca();
                                Float_t nHitsPossB  = trackB.nHitsMax();
                                Float_t nHitsFitB   = trackB.nHitsFit();
                                Float_t MomentumB   = trackB.gMom().mag();
                                //Float_t BetaB       = trackB.btofBeta();
                                //Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                //Float_t TPCdEdxB    = trackB.dEdx(); // Combined inner and outer MDC dE/dx
                                //Float_t nSigmaKaonB = trackB.nSigmaKaon();
                                //Float_t nSigmaPionB = trackB.nSigmaPion();
                                //Float_t nSigmaPB    = trackB.nSigmaProton();

                                if(
                                   nHitsFitB     > min_nhitsfit
                                   && nHitsPossB > 0
                                   && MomentumB  > min_momentum
                                  )
                                {
                                    helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackB.charge());

                                    //cout << "i = " << i << ", j = " << j << ", nHitsFitA = "
                                    //    << nHitsFitA << ", nHitsFitB = " << nHitsFitB  << ", pA = " << MomentumA << ", pB = " << MomentumB << endl;

                                    if(
                                       1 == 1
                                       //((Float_t)nHitsFitA)/((Float_t)nHitsPossA)    > 0.55
                                       //&& ((Float_t)nHitsFitB)/((Float_t)nHitsPossB) > 0.55
                                      )
                                    {
                                        //Float_t pathA2_f, pathB2_f,dcaBB2_f,dcaAB2_f;
                                        //fHelixAtoLinedca(dirY,baseY,helixA,pathB2_f,pathA2_f,dcaAB2_f); // calculates dca of helix to line
                                        //fHelixAtoLinedca(dirY,baseY,helixB,pathB2_f,pathA2_f,dcaBB2_f); // calculates dca of helix to line

                                        Float_t pathA_f, pathB_f, dcaAB_f;
                                        fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
                                        vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                                        vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                                        vectorAB       = vectoratsA+vectoratsB;
                                        vectorAB       = vectorAB/2.0; // dca vertex
                                        if(
                                           dcaAB_f < 3.0
                                           //&& dcaAB2_f < 20.0
                                           //&& dcaBB2_f < 20.0
                                           //&& (vectorAB.x()*vectorAB.x() + vectorAB.y()*vectorAB.y()) < 30.0
                                           //&& fabs(vectorAB.x()) < 60.0
                                           //&& fabs(vectorAB.y()) < 60.0
                                           //&& fabs(vectorAB.z()) < 200.0
                                          )
                                        {
                                            // vertex values
                                            x_table2[comb_counter2][0] = vectorAB.x();
                                            y_table2[comb_counter2][0] = vectorAB.y();
                                            z_table2[comb_counter2][0] = vectorAB.z();

                                            // weights
                                            x_table2[comb_counter2][1] = (nHitsFitA/nHitsPossA)*(nHitsFitB/nHitsPossB);
                                            y_table2[comb_counter2][1] = (nHitsFitA/nHitsPossA)*(nHitsFitB/nHitsPossB);
                                            z_table2[comb_counter2][1] = (nHitsFitA/nHitsPossA)*(nHitsFitB/nHitsPossB);

                                            x_table2[comb_counter2][1] = 1.0;
                                            y_table2[comb_counter2][1] = 1.0;
                                            z_table2[comb_counter2][1] = 1.0;

                                            comb_counter2++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if(
               comb_counter2 > 1
               //&& fraction_comb_counter > 20.0
              )  // 40
            {
                x_mean2 = CalcTruncMean(x_table2,1.0,0.3,comb_counter2);
                y_mean2 = CalcTruncMean(y_table2,1.0,0.3,comb_counter2);
                z_mean2 = CalcTruncMean(z_table2,1.0,0.3,comb_counter2);
            }

            //cout << "Secondary vertex = {" << x_mean2 << ", " << y_mean2 << ", " << z_mean2 << "}" << endl;


            // Calculating the number of primaries
            Float_t eta_sum_old   = 0;
            //Int_t eta_counter_old = 0;
            Float_t eta_sum_new   = 0;
            //Int_t eta_counter_new = 0;

            vector_dist_prim_second  = sqrt((x_mean-x_mean2)*(x_mean-x_mean2)+(y_mean-y_mean2)*(y_mean-y_mean2)+(z_mean-z_mean2)*(z_mean-z_mean2));

            Vertex_NTDataArray[0]     =(Float_t)x_mean;         //
            Vertex_NTDataArray[1]     =(Float_t)y_mean;         //
            Vertex_NTDataArray[2]     =(Float_t)z_mean;         //
            Vertex_NTDataArray[3]     =(Float_t)run_events;
            Vertex_NTDataArray[4]     =(Float_t)EventVertexX;
            Vertex_NTDataArray[5]     =(Float_t)EventVertexY;
            Vertex_NTDataArray[6]     =(Float_t)EventVertexZ;
            Vertex_NTDataArray[7]     =(Float_t)comb_counter;
            Vertex_NTDataArray[8]     =(Float_t)0;
            Vertex_NTDataArray[9]     =(Float_t)BBC;
            Vertex_NTDataArray[10]    =(Float_t)VPD;
            Vertex_NTDataArray[11]    =(Float_t)N_primaries_old;
            Vertex_NTDataArray[12]    =(Float_t)N_primaries_new;
            Vertex_NTDataArray[13]    =(Float_t)fraction_comb_counter;
            Vertex_NTDataArray[14]    =(Float_t)eta_sum_old;
            Vertex_NTDataArray[15]    =(Float_t)eta_sum_new;
            Vertex_NTDataArray[16]    =(Float_t)x_mean2;         //
            Vertex_NTDataArray[17]    =(Float_t)y_mean2;         //
            Vertex_NTDataArray[18]    =(Float_t)z_mean2;         //
            Vertex_NTDataArray[19]    =(Float_t)N_primaries_new2;
            Vertex_NTDataArray[20]    =(Float_t)refMult;
            Vertex_NTDataArray[21]    =(Float_t)N_tofmatch_new;
            Vertex_NTDataArray[22]    =(Float_t)N_tofmatch_new2;
            Vertex_NTDataArray[23]    =(Float_t)vector_dist_prim_second;
            Vertex_NTDataArray[24]    =(Float_t)RunId;
            Vertex_NTDataArray[25]    =(Float_t)N_tofmatch_old;
            Vertex_NTDataArray[26]    =(Float_t)trigger_word;
            Vertex_NTDataArray[27]    =(Float_t)x_meanA;         //
            Vertex_NTDataArray[28]    =(Float_t)y_meanA;         //
            Vertex_NTDataArray[29]    =(Float_t)z_meanA;         //
            Vertex_NTDataArray[30]    =(Float_t)x_meanB;         //
            Vertex_NTDataArray[31]    =(Float_t)y_meanB;         //
            Vertex_NTDataArray[32]    =(Float_t)z_meanB;         //
            Vertex_NTDataArray[33]    =(Float_t)EventId;         //

            comb_counter_global = comb_counter;

            if((fAnalysisNum == 9) && Save_flag == 1)
            {
                Vertex_NT->Fill(Vertex_NTDataArray);
                vertex_counter++;
            }

            if(
               (x_mean*x_mean + y_mean*y_mean) < 4.0
               && fabs(z_mean) < 70.0
              )
            {
                good_event++;
                //cout << "event_number = " << event_number << ", run_events = " << run_events << ", comb_counter = "
                //    << comb_counter << ", fraction used = " << fraction_comb_counter
                //    << "%, vertex = {" << x_mean << ", " << y_mean << ", " << z_mean << "}" << endl;
            }
        }
        else
        {
            bad_event++;
        }

        //Float_t new_old_diff = (((x_mean-EventVertexX)*(x_mean-EventVertexX))+((y_mean-EventVertexY)*(y_mean-EventVertexY))+((z_mean-EventVertexZ)*(z_mean-EventVertexZ)));

        return 1;
    }
    else return 0;
}


/*
Int_t K_Pi_three_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs], Int_t PID_counter_Array_D[][N_max_PIDs],
                          Int_t PID_Array[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks], Int_t PID_Array_D[][N_max_PIDs][N_max_tracks],
                          StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B, StPicoAlexEvent* picoDst_D,
                          Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t ParticleD,
                          Int_t Ana_Num, Int_t run_events, Int_t run_events_B, Int_t event_number, Int_t SE_ME_Flag, Int_t eAnalysisNum)
{
    // mass =  493.677 MeV/c2
    // ctau = 2.461 cm
    //
    // Au + Au -> K+/- + N + N
    //             |
    //             -> pi+/- + pi+/- + pi-/+

    dummy_counter_A++;
    // Event vertex information
    StThreeVectorD vector_prim,vector_primB,vectordiff, vectornewA_lin, vectornewB_lin;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB,RunIdA,RunIdB;

    event_A_ana   = picoDst_A;
    event_B_ana   = picoDst_B;

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();

    Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_A_ana->vzVpd();

    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    vector_prim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        RunIdB         = event_B_ana->runId();
        vector_primB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vector_prim - vector_primB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex

    }

    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing
    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
    }

    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]   > 0   // p(A)/bar(p(A))
       && PID_counter_Array[Ana_Num][ParticleB]   > 0   // pi-/+(B)
       && (PID_counter_Array_B[Ana_Num][ParticleC] > 0 || PID_counter_Array_B[Ana_Num][ParticleD] > 0)   //
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana->isMinBias()
       && event_B_ana->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
        //   << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;
        // Loop over all particle combinations

        Double_t dcaA_cut     = 2.0; // 2.4 all in cm, optimized cuts by Kamila
        Double_t dcaB_cut     = 2.0; // 2.4
        Double_t dcaC_cut     = 2.0;  // 2.4
        Double_t VerdistX_cut = 18.0; // 20.0

        Double_t dca_Kaon_cut = 6.0; // 4.9
        Double_t dcaAB_f_cut  = 2.0; // 1
        Double_t dcaABC_cut   = 2.0; // 1.2

        StPhysicalHelixD helixA, helixB;
        StThreeVectorD vectorA, vectorB, vectoratsA, vectoratsB, vectorAB, vectorAB_est, vector_primAB, vectornewA, vectornewB, dirY_lin;
        StThreeVectorD testA, testB, testAB, vectorABtoPrim, baseY, dirY;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackB2, ltrackD, ltrackA_lin, ltrackB_lin, ltrackB2_pi;


        // Fill event information for K+/-
        alexV0_event_A.clearTrackList();
        alexV0_event_A.setx(EventVertexXA);
        alexV0_event_A.sety(EventVertexYA);
        alexV0_event_A.setz(EventVertexZA);
        alexV0_event_A.setid(RunIdA);
        alexV0_event_A.setmult(refMultA);
        alexV0_event_A.setn_prim(n_primaries);
        alexV0_event_A.setn_non_prim(n_non_primaries);
        alexV0_event_A.setn_tof_prim(n_tofmatch_prim);
        alexV0_event_A.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        alexV0_event_A.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        alexV0_event_A.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        alexV0_event_A.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        alexV0_event_A.setEP_Qx_ptw(EP_Qx_ptw);
        alexV0_event_A.setEP_Qy_ptw(EP_Qy_ptw);
        alexV0_event_A.setQtracks_eta_pos(Qtracks_used_eta_pos);
        alexV0_event_A.setQtracks_eta_neg(Qtracks_used_eta_neg);
        alexV0_event_A.setQtracks_full(Qtracks_used);
        alexV0_event_A.setZDCx(ZDCx);
        alexV0_event_A.setBBCx(BBCx);
        alexV0_event_A.setvzVpd(vzVpd);


        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // p+/- candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t BetaA       = trackA.btofBeta();
            //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            //Float_t TPCdEdxA    = trackA.dEdx(); // Combined inner and outer MDC dE/dx
            //Float_t nSigmaKaonA = trackA.nSigmaKaon();
            Float_t nSigmaPionA = trackA.nSigmaPion();
            //Float_t nSigmaPA    = trackA.nSigmaProton();

            if(
               dcaA          > dcaA_cut  // 15.0
               && nHitsFitA  > 10
               && nHitsPossA > 0.0
               && MomentumA  > 0.1
               && MomentumA  < 10.0
               && (nHitsFitA/nHitsPossA) > 0.52
               && nSigmaPionA > -3.0
               && nSigmaPionA < 3.0
              )
            {
                helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackA.charge());


                for(Int_t j = i+1; j < PID_counter_Array[Ana_Num][ParticleB]; j++)  // pi+/- candidates
                {
                    Int_t trackB_num = PID_Array[Ana_Num][ParticleB][j];
                    StPicoAlexTrack trackB = *picoDst_A->track( trackB_num );

                    Float_t dcaB        = trackB.dca();
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t BetaB       = trackB.btofBeta();
                    //Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    //Float_t TPCdEdxB    = trackB.dEdx(); // Combined inner and outer MDC dE/dx
                    //Float_t nSigmaKaonB = trackB.nSigmaKaon();
                    Float_t nSigmaPionB = trackB.nSigmaPion();
                    //Float_t nSigmaPB    = trackB.nSigmaProton();

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       && dcaB       > dcaB_cut // 15.0
                       && nHitsFitB  > 10
                       && nHitsPossB > 0.0
                       && MomentumB  > 0.1
                       && MomentumB  < 10.0
                       && (nHitsFitB/nHitsPossB) > 0.52
                       && nSigmaPionB > -3.0
                       && nSigmaPionB < 3.0
                      )
                    {
                        helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackB.charge());

                        //cout << "i = " << i << ", j = " << j << endl;
                        Float_t pathA_f, pathB_f, dcaAB_f;
                        Float_t pathA_test = 0.0;
                        Float_t pathB_test = 0.0;
                        Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);

                        vectoratsA     = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
                        vectoratsB     = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
                        vectorAB       = vectoratsA+vectoratsB;
                        vectorAB       = vectorAB/2.0; // decay vertex

                        //Double_t vex_est_x = vectorAB.x();
                        //Double_t vex_est_y = vectorAB.y();
                        //Double_t vex_est_z = vectorAB.z();

                        StThreeVectorF baseA,dirA,baseB,dirB;
                        baseA = helixA.at(pathA_test);
                        baseB = helixB.at(pathB_test);
                        dirA  = helixA.at(pathA_test-2.0) - helixA.at(pathA_test+2.0);
                        dirB  = helixB.at(pathB_test-2.0) - helixB.at(pathB_test+2.0);

                        StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

                        //Double_t vex_lin_x = vectorAB_lin.x();
                        //Double_t vex_lin_y = vectorAB_lin.y();
                        //Double_t vex_lin_z = vectorAB_lin.z();

                        Double_t dcaAB_lin = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
                        StThreeVectorD vectorABtoPrim_lin = vectorAB_lin - vector_prim; // vector primary vertex to decay vertex
                        Float_t VerdistX_lin = vectorABtoPrim_lin.mag(); // distance between primary vertex and decay verte


                        // calculate the scalar product with the approximated secondary vertex position
                        vectornewA_lin     = helixA.cat(pathA_test); // direction vector at dca for helixA
                        vectornewB_lin     = helixB.cat(pathB_test); // direction vector at dca for helixB
                        vectornewA_lin = MomentumA*vectornewA_lin/vectornewA_lin.mag(); // new momentum vector at decay vertex
                        vectornewB_lin = MomentumB*vectornewB_lin/vectornewB_lin.mag(); // new momentum vector at decay vertex
                        ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),0.13957018);
                        ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),0.13957018);
                        TLorentzVector trackAB_lin      = ltrackA_lin+ltrackB_lin; // mother particle
                        dirY_lin.set(trackAB_lin.Px(),trackAB_lin.Py(),trackAB_lin.Pz());
                        dirY_lin = dirY_lin/dirY_lin.mag();
                        Double_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim_lin/vectorABtoPrim_lin.mag());


                        if( VerdistX_lin > 0.7*VerdistX_cut && dcaAB_lin < 1.5*dcaAB_f_cut && scalarProduct_lin > 0.0 )
                        {
                            if(fDCA_Helix_out == 1)
                            {
                                fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
                            }
                            else
                            {
                                fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
                            }

                            vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                            vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                            vectorAB       = vectoratsA+vectoratsB;
                            vectorAB       = vectorAB/2.0; // decay vertex
                            vectorABtoPrim = vectorAB - vector_prim; // vector primary vertex to decay vertex
                            Float_t VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

                            vectornewA     = helixA.cat(pathA_f); // direction vector at dca for helixA
                            vectornewB     = helixB.cat(pathB_f); // direction vector at dca for helixB

                            vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
                            vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex

                            ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.13957018);
                            ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),0.13957018);

                            // Missing mass and invariant mass calculations
                            TLorentzVector trackAB      = ltrackA+ltrackB; // mother particle

                            Float_t Mass2A        = -100.0;
                            // calculate mass2
                            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
                            {
                                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
                            }
                            //Float_t MassA       = SquareRoot(Mass2A);

                            Float_t Mass2B        = -100.0;
                            // calculate mass2
                            if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                            {
                                Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                            }
                            //Float_t MassB       = SquareRoot(Mass2B);

                            if(
                               VerdistX          > VerdistX_cut // 60.0
                               && dcaAB_f        < dcaAB_f_cut  // 3.0
                              )
                            {

                                for(Int_t k = 0; k < PID_counter_Array_B[Ana_Num][ParticleC]; k++)  // p-/+ candidates
                                {
                                    Int_t trackC_num = PID_Array_B[Ana_Num][ParticleC][k];
                                    if(
                                       ((trackC_num != trackA_num)
                                       && (trackC_num != trackB_num))
                                       || (SE_ME_Flag == 1)
                                      )
                                    {
                                        StPhysicalHelixD helixC;
                                        StThreeVectorD vectorC, vectoratsC, vectoratsA2, vectorA2C, vectorA2CtoPrim, vectornewC, dirY2, baseY2;

                                        StPicoAlexTrack trackC = *picoDst_B->track( trackC_num );

                                        Float_t dcaC        = trackC.dca();
                                        Float_t nHitsPossC  = trackC.nHitsMax();
                                        Float_t nHitsFitC   = trackC.nHitsFit();
                                        Float_t MomentumC   = trackC.gMom().mag();
                                        Float_t BetaC       = trackC.btofBeta();
                                        //Float_t PolarityC   = trackC.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                        //Float_t TPCdEdxC    = trackC.dEdx(); // Combined inner and outer MDC dE/dx
                                        //Float_t nSigmaKaonC = trackC.nSigmaKaon();
                                        Float_t nSigmaPionC = trackC.nSigmaPion();
                                        //Float_t nSigmaPC    = trackC.nSigmaProton();

                                        if(
                                           dcaC          > dcaC_cut  // 15.0
                                           && nHitsFitC  > 10
                                           && nHitsPossC > 0
                                           && MomentumC  > 0.1
                                           && MomentumC  < 10.0
                                           && (nHitsFitC/nHitsPossC) > 0.52
                                           && nSigmaPionC > -3.0
                                           && nSigmaPionC < 3.0
                                          )
                                        {
                                            helixC = StPhysicalHelixD(trackC.gMom(),trackC.origin()+ vectordiff,event_B_ana->bField()*MAGFIELDFACTOR,trackC.charge());

                                            Float_t Mass2C        = -100.0;
                                            // calculate mass2
                                            if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                                            {
                                                Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                                            }
                                            //Float_t MassC       = SquareRoot(Mass2C);


                                            Float_t dcaABC  = -999.0;
                                            Float_t pathC_f = 0.0;
                                            fHelixAtoPointdca(vectorAB,helixC,pathC_f,dcaABC);
                                            vectoratsC = helixC.at(pathC_f);

                                            StThreeVectorD vectorABC;
                                            vectorABC = (vectorAB + vectoratsC)/2.0;

                                            if(dcaABC < dcaABC_cut)
                                            {
                                                Float_t pathA_g, dcaA_g, pathB_g, dcaB_g, pathC_g, dcaC_g;
                                                fHelixAtoPointdca(vectorABC,helixA,pathA_g,dcaA_g);
                                                fHelixAtoPointdca(vectorABC,helixB,pathB_g,dcaB_g);
                                                fHelixAtoPointdca(vectorABC,helixC,pathC_g,dcaC_g);

                                                vectornewA     = helixA.cat(pathA_g); // direction vector at dca for helixA
                                                vectornewB     = helixB.cat(pathB_g); // direction vector at dca for helixB
                                                vectornewC     = helixC.cat(pathB_g); // direction vector at dca for helixB

                                                vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
                                                vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex
                                                vectornewC = MomentumC*vectornewC/vectornewC.mag(); // new momentum vector at decay vertex

                                                ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.13957018);
                                                ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),0.13957018);
                                                ltrackC.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.13957018);

                                                TLorentzVector ltrackABC = ltrackA + ltrackB + ltrackC; // mother particle
                                                Double_t InvMassABC      = ltrackABC.M(); // invariant mass of mother particle
                                                //Double_t MomentumABC     = ltrackABC.P(); // momentum of mother particle
                                                Float_t  pt2             = ltrackABC.Pt();  // Transverse momentum of mother particle
                                                Float_t  rap2            = ltrackABC.Rapidity(); // Rapidity of mother particle

                                                dirY2.set(ltrackABC.Px(),ltrackABC.Py(),ltrackABC.Pz()); // direction vector of Xi
                                                dirY2 = dirY2/dirY2.mag();
                                                StThreeVectorD vectorA2CtoPrim   = vectorABC - vector_prim; // vector primary vertex to decay vertex
                                                Double_t scalarProduct2 = dirY2.dot(vectorA2CtoPrim/vectorA2CtoPrim.mag());
                                                Float_t phiABC   = dirY2.phi();
                                                Float_t thetaABC = dirY2.theta();


                                                 // Calculate the helix of the Kaon 
                                                 StPhysicalHelixD helix_Kaon;
                                                 StThreeVectorD origin_Kaon;   // decay vertex of Kaon
                                                 origin_Kaon.set(vectorABC.x(),vectorABC.y(),vectorABC.z());
                                                 Float_t h_Kaon      = helixA.h(); //
                                                 Float_t phase_Kaon  = TMath::ATan2(-1.0*h_Kaon*ltrackABC.Px(),h_Kaon*ltrackABC.Py());
                                                 Float_t dip_Kaon    = TMath::ATan2(ltrackABC.Pz(),pt2);  // correct
                                                 Float_t curv_Kaon   = 1.0;
                                                 if(pt2 != 0.0)
                                                 {
                                                     curv_Kaon = curv_to_invpt_ratio/pt2;
                                                 }

                                                 helix_Kaon.setParameters(curv_Kaon,dip_Kaon,phase_Kaon,origin_Kaon,h_Kaon);
                                                 Float_t path_Kaon  = -999.0;
                                                 Float_t dca_Kaon   = -999.0;
                                                 fHelixAtoPointdca(vector_prim,helix_Kaon,path_Kaon,dca_Kaon);

                                                 Float_t path_Kaon2  = -999.0;
                                                 Float_t dca_Kaon2   = -999.0;
                                                 fHelixAtoPointdca(vectorABC,helix_Kaon,path_Kaon2,dca_Kaon2);
                                                 Double_t decay_length = path_Kaon2-path_Kaon;
                                                 // dca_omega is almost identical to VerdistY2


                                                if(
                                                   InvMassABC        < 0.7 //
                                                   && dca_Kaon       < dca_Kaon_cut // 6.0
                                                  )
                                                {
                                                    Float_t dip_Kaon   = TMath::ATan2(ltrackABC.Pz(),pt2);  // correct

                                                    //------------------------------------------------------------------------------------
                                                    // beta tof correction
                                                    Float_t Mass2ACorr = -999.;
                                                    StLorentzVectorD SttrackABC(ltrackABC.X(),ltrackABC.Y(),ltrackABC.Z(),ltrackABC.E());
                                                    if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 )
                                                    { // pi+/-
                                                        Float_t BetaAOrig = trackA.btof();
                                                        Float_t BetaACorr;
                                                        StThreeVectorD vtofA = trackA.btofHisPos();
                                                        v0tofcorr->setVectors3D(vector_prim)(vectorABC)(vtofA);
                                                        v0tofcorr->setMotherTracks(SttrackABC);
                                                        v0tofcorr->correctBeta(helixA,BetaAOrig,BetaACorr);
                                                        v0tofcorr->clearContainers();
                                                        //cout << "BetaAOrig = " << BetaAOrig << "   BetaACorr = " << BetaACorr << endl;
                                                        Mass2ACorr = MomentumA*MomentumA*(1.0/(BetaACorr*BetaACorr) - 1.0);
                                                        //cout << "Mass2AOrig = " << Mass2A << "   Mass2ACorr = " << Mass2ACorr << endl;
                                                    }

                                                    Float_t Mass2BCorr = -999.;
                                                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 )
                                                    { // pi+/-
                                                        Float_t BetaBOrig = trackB.btof();
                                                        Float_t BetaBCorr;
                                                        StThreeVectorD vtofB = trackB.btofHisPos();
                                                        v0tofcorr->setVectors3D(vector_prim)(vectorABC)(vtofB);
                                                        v0tofcorr->setMotherTracks(SttrackABC);
                                                        v0tofcorr->correctBeta(helixB,BetaBOrig,BetaBCorr);
                                                        v0tofcorr->clearContainers();
                                                        //cout << "BetaAOrig = " << BetaAOrig << "   BetaACorr = " << BetaACorr << endl;
                                                        Mass2BCorr = MomentumB*MomentumB*(1.0/(BetaBCorr*BetaBCorr) - 1.0);
                                                        //cout << "Mass2AOrig = " << Mass2A << "   Mass2ACorr = " << Mass2ACorr << endl;
                                                    }

                                                    Float_t Mass2CCorr = -999.;
                                                    if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 )
                                                    { // pi+/-
                                                        Float_t BetaCOrig = trackC.btof();
                                                        Float_t BetaCCorr;
                                                        StThreeVectorD vtofC = trackC.btofHisPos();
                                                        v0tofcorr->setVectors3D(vector_prim)(vectorABC)(vtofC);
                                                        v0tofcorr->setMotherTracks(SttrackABC);
                                                        v0tofcorr->correctBeta(helixC,BetaCOrig,BetaCCorr);
                                                        v0tofcorr->clearContainers();
                                                        //cout << "BetaAOrig = " << BetaAOrig << "   BetaACorr = " << BetaACorr << endl;
                                                        Mass2CCorr = MomentumC*MomentumC*(1.0/(BetaCCorr*BetaCCorr) - 1.0);
                                                        //cout << "Mass2AOrig = " << Mass2A << "   Mass2ACorr = " << Mass2ACorr << endl;
                                                    }
                                                    //------------------------------------------------------------------------------------



                                                    // Event plane calculations 
                                                    // Calculate the same vectors as for the event plane -> no V0!
                                                    Float_t pathA_prim,dcaA_prim,pathB_prim,dcaB_prim,pathC_prim,dcaC_prim;
                                                    fHelixAtoPointdca(vector_prim,helixA,pathA_prim,dcaA_prim);
                                                    fHelixAtoPointdca(vector_prim,helixB,pathB_prim,dcaB_prim);
                                                    fHelixAtoPointdca(vector_prim,helixC,pathC_prim,dcaC_prim);

                                                    StThreeVectorF vectornewA_prim,vectornewB_prim,vectornewC_prim;
                                                    vectornewA_prim     = helixA.cat(pathA_prim);
                                                    vectornewB_prim     = helixB.cat(pathB_prim);
                                                    vectornewC_prim     = helixC.cat(pathC_prim);

                                                    vectornewA_prim     = MomentumA*vectornewA_prim/vectornewA_prim.mag(); // momentum vector at primary vertex
                                                    vectornewB_prim     = MomentumB*vectornewB_prim/vectornewB_prim.mag(); // momentum vector at primary vertex
                                                    vectornewC_prim     = MomentumC*vectornewC_prim/vectornewC_prim.mag(); // momentum vector at primary vertex

                                                    Float_t phi_event_plane                = -400.0;
                                                    Float_t phi_event_plane_eta_gap        = -400.0;
                                                    Float_t delta_phi_ME_AB_weight         = 0.0;
                                                    Float_t delta_phi_ME_AB_weight_eta_gap = 0.0;

                                                    // Calculate the event plane anlges
                                                    calc_event_plane_angles(3,rap2,trackA,trackB,trackC,vectornewA_prim,vectornewB_prim,vectornewC_prim,ME_Flag,SE_ME_Flag,RunIdA,EventVertexXA,EventVertexYA,
                                                                            EventVertexZA,RunIdB,EventVertexXB,EventVertexYB,EventVertexZB,
                                                                            phi_event_plane,phi_event_plane_eta_gap,delta_phi_ME_AB_weight,delta_phi_ME_AB_weight_eta_gap);


                                                    Int_t delta_phi_ME = 0;
                                                    // check whether the event planes between event A and B are close to each other
                                                    if(
                                                       fabs(delta_phi_ME_AB_weight) < TMath::DegToRad()*30.0
                                                       || (SE_ME_Flag == 0)
                                                      )
                                                    {
                                                        delta_phi_ME = 1;
                                                    }


                                                    delta_phi_ME = 1; // ignore the delta_phi_ME

                                                    if(delta_phi_ME == 1)
                                                    {
                                                        Float_t p_xA_c   = vectornewA_prim.x();
                                                        Float_t p_yA_c   = vectornewA_prim.y();
                                                        Float_t p_tA_c   = sqrt(p_xA_c*p_xA_c + p_yA_c*p_yA_c);
                                                        Float_t etaA_c   = vectornewA_prim.pseudoRapidity();
                                                        Float_t phiA_c   = vectornewA_prim.phi();
                                                        Double_t p_t_weightA = 1.0;
                                                        if(p_tA_c < 2.0)  p_t_weightA = p_tA_c;
                                                        if(p_tA_c >= 2.0) p_t_weightA = 2.0;
                                                        Float_t iQxA     = p_t_weightA*TMath::Cos(2.0*phiA_c);
                                                        Float_t iQyA     = p_t_weightA*TMath::Sin(2.0*phiA_c);

                                                        Float_t p_xB_c   = vectornewB_prim.x();
                                                        Float_t p_yB_c   = vectornewB_prim.y();
                                                        Float_t p_tB_c   = sqrt(p_xB_c*p_xB_c + p_yB_c*p_yB_c);
                                                        Float_t etaB_c   = vectornewB_prim.pseudoRapidity();
                                                        Float_t phiB_c   = vectornewB_prim.phi();
                                                        Double_t p_t_weightB = 1.0;
                                                        if(p_tB_c < 2.0)  p_t_weightB = p_tB_c;
                                                        if(p_tB_c >= 2.0) p_t_weightB = 2.0;
                                                        Float_t iQxB     = p_t_weightB*TMath::Cos(2.0*phiB_c);
                                                        Float_t iQyB     = p_t_weightB*TMath::Sin(2.0*phiB_c);

                                                        Float_t p_xC_c   = vectornewC_prim.x();
                                                        Float_t p_yC_c   = vectornewC_prim.y();
                                                        Float_t p_tC_c   = sqrt(p_xC_c*p_xC_c + p_yC_c*p_yC_c);
                                                        Float_t etaC_c   = vectornewC_prim.pseudoRapidity();
                                                        Float_t phiC_c   = vectornewC_prim.phi();
                                                        Double_t p_t_weightC = 1.0;
                                                        if(p_tC_c < 2.0)  p_t_weightC = p_tC_c;
                                                        if(p_tC_c >= 2.0) p_t_weightC = 2.0;
                                                        Float_t iQxC     = p_t_weightC*TMath::Cos(2.0*phiC_c);
                                                        Float_t iQyC     = p_t_weightC*TMath::Sin(2.0*phiC_c);




                                                        alexV0_track_A = alexV0_event_A.createTrack();
                                                        alexV0_track_A->setm2A(Mass2ACorr); //
                                                        alexV0_track_A->setm2B(Mass2BCorr); //
                                                        alexV0_track_A->setm2C(Mass2CCorr); //
                                                        alexV0_track_A->setnsA(nSigmaPionA);//
                                                        alexV0_track_A->setnsB(nSigmaPionB);//
                                                        alexV0_track_A->setnsC(nSigmaPionC);//
                                                        alexV0_track_A->setdcaA(dcaA);  //
                                                        alexV0_track_A->setdcaB(dcaB);  //
                                                        alexV0_track_A->setdcaC(dcaC);  //
                                                        alexV0_track_A->setiQxA(iQxA);  //
                                                        alexV0_track_A->setiQyA(iQyA);  //
                                                        alexV0_track_A->setiQxB(iQxB);  //
                                                        alexV0_track_A->setiQyB(iQyB);  //
                                                        alexV0_track_A->setiQxC(iQxC);  //
                                                        alexV0_track_A->setiQyC(iQyC);  //
                                                        alexV0_track_A->setetaA(etaA_c);  //
                                                        alexV0_track_A->setetaB(etaB_c);  //
                                                        alexV0_track_A->setetaC(etaC_c);  //
                                                        alexV0_track_A->setInvAB(trackC.charge());  //
                                                        alexV0_track_A->setInvABC(InvMassABC);      //
                                                        alexV0_track_A->setInvAB_miss(0);           //
                                                        alexV0_track_A->setInvABC_miss(0);          //
                                                        alexV0_track_A->setdcaAB(dcaAB_f);          //
                                                        alexV0_track_A->setdcaBC(0);                //
                                                        alexV0_track_A->setdcaABC(dip_Kaon);        //
                                                        alexV0_track_A->setVerdistX(VerdistX);      //
                                                        alexV0_track_A->setVerdistY(dca_Kaon);      //
                                                        alexV0_track_A->setVerdistX2(decay_length); //
                                                        alexV0_track_A->setVerdistY2(0);            //
                                                        alexV0_track_A->setpt(pt2);                 //
                                                        alexV0_track_A->setrap(rap2);               //
                                                        alexV0_track_A->setphi(phiABC);             //
                                                        alexV0_track_A->settheta(thetaABC);         //
                                                        alexV0_track_A->setPsi_ep(phi_event_plane); //
                                                        alexV0_track_A->setPsi_ep_eta(phi_event_plane_eta_gap); //
                                                        alexV0_track_A->setPsi_diff_ME(delta_phi_ME_AB_weight); //
                                                        alexV0_track_A->setscal_prod(0); //
                                                        alexV0_track_A->setscal_prod2(scalarProduct2); //

#if 0
                                                        // Filling Ntuple
                                                        Kaon_NTDataArray[0]     =(Float_t)VerdistX;
                                                        Kaon_NTDataArray[1]     =(Float_t)InvMassABC;
                                                        Kaon_NTDataArray[2]     =(Float_t)MomentumABC;
                                                        Kaon_NTDataArray[3]     =(Float_t)dcaABC;
                                                        Kaon_NTDataArray[4]     =(Float_t)dcaAB_f;
                                                        Kaon_NTDataArray[5]     =(Float_t)pt2;
                                                        Kaon_NTDataArray[6]     =(Float_t)rap2;
                                                        Kaon_NTDataArray[7]     =(Float_t)dcaA;
                                                        Kaon_NTDataArray[8]     =(Float_t)dcaB;
                                                        Kaon_NTDataArray[9]     =(Float_t)dcaC;
                                                        Kaon_NTDataArray[10]    =(Float_t)MomentumA;
                                                        Kaon_NTDataArray[11]    =(Float_t)MomentumB;
                                                        Kaon_NTDataArray[12]    =(Float_t)MomentumC;
                                                        Kaon_NTDataArray[13]    =(Float_t)nSigmaPionA;
                                                        Kaon_NTDataArray[14]    =(Float_t)nSigmaPionB;
                                                        Kaon_NTDataArray[15]    =(Float_t)nSigmaPionC;
                                                        Kaon_NTDataArray[16]    =(Float_t)dca_Kaon;
                                                        Kaon_NTDataArray[17]    =(Float_t)vectorABC.x();
                                                        Kaon_NTDataArray[18]    =(Float_t)vectorABC.y();
                                                        Kaon_NTDataArray[19]    =(Float_t)vectorABC.z();
                                                        Kaon_NTDataArray[20]    =(Float_t)decay_length;

                                                        Kaon_NT->Fill(Kaon_NTDataArray);
#endif
                                                        dummy_counter++;
                                                        dummy_counter_loop++;
                                                    }

                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        Tree_KaonV0_v2  ->Fill();

        return 1;
    }
    else return 0;
}
*/


Int_t Merged_Omega_Xi_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs],
                  Int_t PID_Array[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                  StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                  Int_t ParticleA, Int_t ParticleB,Int_t Ana_Num,Int_t SE_ME_Flag)
{
    // Omega
    // mass = 1672.45 MeV/c2
    // ctau = 2.461 cm
    //
    // Au + Au -> Omega-/+ + N + N
    //             |
    //             -> Lambda(anti-) + K-/+(C)
    //                  |
    //                  -> p(A) + pi-(B)

    // Xi
    // mass = 1321.31 MeV/c2
    // ctau = 4.91 cm
    //
    // Au + Au -> Xi-/+ + N + N
    //             |
    //             -> Lambda(anti-) + pi-/+(C)
    //                  |
    //                  -> p(A) + pi-(B)

    dummy_counter_A++;
    // Event vertex information
    StThreeVectorF vectorprim,vectorprimB,vectordiff, vectornewA_lin, vectornewB_lin;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB,RunIdA,RunIdB;


    event_A_ana       = picoDst_A;
    event_B_ana       = picoDst_A;
    event_C_ana       = picoDst_B;

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();

    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    refMultB          = refMultA;
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_A_ana->vzVpd();

    vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    Int_t flag_tof_particleA = 0;
    Int_t flag_tof_particleB = 0;
    Int_t flag_tof_particleC = 0;

    Int_t ParticleC = 0;

    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        EventVertexXB  = event_C_ana->primaryVertex().x();
        EventVertexYB  = event_C_ana->primaryVertex().y();
        EventVertexZB  = event_C_ana->primaryVertex().z();
        refMultB       = event_C_ana->refMult();
        RunIdB         = event_C_ana->runId();
        vectorprimB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vectorprim - vectorprimB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }


    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
        //cout << "Triggers checked, Vertex checked..." << endl;
    }

    //
    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]    > 0   // p(A)/bar(p(A))
       && PID_counter_Array[Ana_Num][ParticleB]    > 0   // pi-/+(B)
       && ((ParticleA == 14 && (PID_counter_Array_B[Ana_Num][9] > 0 || PID_counter_Array_B[Ana_Num][12] > 0)) ||
           (ParticleA == 15 && (PID_counter_Array_B[Ana_Num][8] > 0 || PID_counter_Array_B[Ana_Num][11] > 0))
          )
       //&& event->getNVertices() >= 1
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana   ->isMinBias()
       && event_C_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
        //   << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;
        // Loop over all particle combinations

        Double_t dcaA_cut        = 0.0;
        Double_t dcaB_cut        = 0.9;
        Double_t dcaC_cut        = 0.5;
        Double_t VerdistX_cut    = 3.5;
        Double_t InvMassABC_cut  = 2.5;
        Double_t dca_mother_cut  = 0.7;
        Double_t dcaBC_cut       = 1.0;
        Double_t VerdistX2_cut   = 2.5;


        StPhysicalHelixD helixA, helixB;
        StThreeVectorF vectorA, vectorB, vectoratsA, vectoratsB, vectorAB, vectorAB_est, vectorprimAB, vectornewA, vectornewB, dirY_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackC_pi, ltrackB2, ltrackD, ltrackA_lin, ltrackB_lin, ltrackB2_pi, ltrackA_pip;

        // Fill event information for Omega-/+
        alexV0_event_A.clearTrackList();
        alexV0_event_A.setx(EventVertexXA);
        alexV0_event_A.sety(EventVertexYA);
        alexV0_event_A.setz(EventVertexZA);
        alexV0_event_A.setid(RunIdA);
        alexV0_event_A.setmult(refMultA);
        alexV0_event_A.setn_prim(n_primaries);
        alexV0_event_A.setn_non_prim(n_non_primaries);
        alexV0_event_A.setn_tof_prim(n_tofmatch_prim);
        alexV0_event_A.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        alexV0_event_A.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        alexV0_event_A.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        alexV0_event_A.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        alexV0_event_A.setEP_Qx_ptw(EP_Qx_ptw);
        alexV0_event_A.setEP_Qy_ptw(EP_Qy_ptw);
        alexV0_event_A.setQtracks_eta_pos(Qtracks_used_eta_pos);
        alexV0_event_A.setQtracks_eta_neg(Qtracks_used_eta_neg);
        alexV0_event_A.setQtracks_full(Qtracks_used);
        alexV0_event_A.setZDCx(ZDCx);
        alexV0_event_A.setBBCx(BBCx);
        alexV0_event_A.setvzVpd(vzVpd);

        // Fill event information for Xi-/+
        alexV0_event_B.clearTrackList();
        alexV0_event_B.setx(EventVertexXA);
        alexV0_event_B.sety(EventVertexYA);
        alexV0_event_B.setz(EventVertexZA);
        alexV0_event_B.setid(RunIdA);
        alexV0_event_B.setmult(refMultA);
        alexV0_event_B.setn_prim(n_primaries);
        alexV0_event_B.setn_non_prim(n_non_primaries);
        alexV0_event_B.setn_tof_prim(n_tofmatch_prim);
        alexV0_event_B.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        alexV0_event_B.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        alexV0_event_B.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        alexV0_event_B.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        alexV0_event_B.setEP_Qx_ptw(EP_Qx_ptw);
        alexV0_event_B.setEP_Qy_ptw(EP_Qy_ptw);
        alexV0_event_B.setQtracks_eta_pos(Qtracks_used_eta_pos);
        alexV0_event_B.setQtracks_eta_neg(Qtracks_used_eta_neg);
        alexV0_event_B.setQtracks_full(Qtracks_used);
        alexV0_event_B.setZDCx(ZDCx);
        alexV0_event_B.setBBCx(BBCx);
        alexV0_event_B.setvzVpd(vzVpd);

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) // p candidates
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];

            StPicoAlexTrack trackA  = *event_A_ana->track( trackA_num );
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t nSigmaPA    = trackA.nSigmaProton();
            //Float_t TofA        = trackA.btof();
            Float_t Mass2A      = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                flag_tof_particleA = 1;
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }
            else
            {
                flag_tof_particleA = 0;
            }
            Float_t MassA       = SquareRoot(Mass2A);

            if(
               dcaA          > dcaA_cut  // 0.15
               && nHitsFitA  > 14
               && nHitsPossA > 0.0
               && MomentumA  > 0.1
               && MomentumA  < 10.0
               && (nHitsFitA/nHitsPossA) > 0.52
               && (flag_tof_particleA == 0 || (flag_tof_particleA == 1 && Mass2A > 0.4 && Mass2A < 1.5)) // proton mass
              )
            {
                helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());

                for(Int_t j = 0; j < PID_counter_Array[Ana_Num][ParticleB]; j++)  // pi- candidates
                {
                    Int_t trackB_num = PID_Array[Ana_Num][ParticleB][j];
                    StPicoAlexTrack trackB  = *event_B_ana->track( trackB_num ); // take again event A
                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    //Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                    Float_t nSigmaPiB   = trackB.nSigmaPion();
                    //Float_t TofB        = trackB.btof();
                    Float_t Mass2B      = -100.0;
                    // calculate mass2
                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        flag_tof_particleB = 1;
                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                    }
                    else
                    {
                        flag_tof_particleB = 0;
                    }
                    Float_t MassB       = SquareRoot(Mass2B);

                    //if( flag_tof_particleA == 0 && flag_tof_particleB == 0) // apply more strict cuts if there is no TOF information for the two particles
                    //{
                    //    dcaA_cut  = 0.1;
                    //    dcaB_cut  = 1.5;
                    //    dcaC_cut  = 0.7;
                    //}

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       && dcaB       > dcaB_cut // 1.3
                       && dcaA       > dcaA_cut //
                       && nHitsFitB  > 14
                       && nHitsPossB > 0.0
                       && MomentumB  > 0.1
                       && MomentumB  < 10.0
                       && (nHitsFitB/nHitsPossB) > 0.52
                       && (flag_tof_particleB == 0 || (flag_tof_particleB == 1 && Mass2B > -0.1 && Mass2B < 0.1)) // pion mass
                      )
                    {
                        helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin(),event_B_ana->bField()*MAGFIELDFACTOR, trackB.charge());

                        //cout << "i = " << i << ", j = " << j << endl;
                        Float_t pathA_f, pathB_f, dcaAB_f;
                        Float_t pathA_test = 0.0;
                        Float_t pathB_test = 0.0;
                        Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);

                        vectoratsA     = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
                        vectoratsB     = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
                        vectorAB       = vectoratsA+vectoratsB;
                        vectorAB       = vectorAB/2.0; // decay vertex

                        //Double_t vex_est_x = vectorAB.x();
                        //Double_t vex_est_y = vectorAB.y();
                        //Double_t vex_est_z = vectorAB.z();

                        StThreeVectorF baseA,dirA,baseB,dirB;
                        baseA = helixA.at(pathA_test);
                        baseB = helixB.at(pathB_test);
                        dirA  = helixA.at(pathA_test-2.0) - helixA.at(pathA_test+2.0);
                        dirB  = helixB.at(pathB_test-2.0) - helixB.at(pathB_test+2.0);

                        StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

                        //Double_t vex_lin_x = vectorAB_lin.x();
                        //Double_t vex_lin_y = vectorAB_lin.y();
                        //Double_t vex_lin_z = vectorAB_lin.z();

                        Double_t dcaAB_lin = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
                        StThreeVectorF vectorABtoPrim_lin = vectorAB_lin - vectorprim; // vector primary vertex to decay vertex
                        Float_t VerdistX_lin = vectorABtoPrim_lin.mag(); // distance between primary vertex and decay verte


                        // calculate the scalar product with the approximated secondary vertex position
                        vectornewA_lin     = helixA.cat(pathA_test); // direction vector at dca for helixA
                        vectornewB_lin     = helixB.cat(pathB_test); // direction vector at dca for helixB
                        vectornewA_lin = MomentumA*vectornewA_lin/vectornewA_lin.mag(); // new momentum vector at decay vertex
                        vectornewB_lin = MomentumB*vectornewB_lin/vectornewB_lin.mag(); // new momentum vector at decay vertex
                        ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),0.93827203);
                        ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),0.13957018);
                        TLorentzVector trackAB_lin      = ltrackA_lin+ltrackB_lin; // mother particle
                        dirY_lin.set(trackAB_lin.Px(),trackAB_lin.Py(),trackAB_lin.Pz());
                        dirY_lin = dirY_lin/dirY_lin.mag();
                        Double_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim_lin/vectorABtoPrim_lin.mag());


                        if( VerdistX_lin > 3.5 && dcaAB_lin < 1.5 && scalarProduct_lin > 0.0 )
                        {
                            if(fDCA_Helix_out == 1)
                            {
                                fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
                            }
                            else
                            {
                                fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
                            }

                            vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                            vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                            vectorAB       = vectoratsA+vectoratsB;
                            vectorAB       = vectorAB/2.0; // decay vertex
                            vectorABtoPrim = vectorAB - vectorprim; // vector primary vertex to decay vertex
                            Float_t VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

                            vectornewA     = helixA.cat(pathA_f); // direction vector at dca for helixA
                            vectornewB     = helixB.cat(pathB_f); // direction vector at dca for helixB

                            vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
                            vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex

                            ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.93827203);
                            ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),0.13957018);
                            ltrackA_pip.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.13957018);

                            TLorentzVector trackAB_K0S  = ltrackA_pip+ltrackB; // mother particle
                            Double_t InvMassAB_K0S      = trackAB_K0S.M(); // invariant mass of mother particle

                            // Missing mass and invariant mass calculations
                            TLorentzVector trackAB      = ltrackA+ltrackB; // mother particle
                            Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                            Float_t MomentumAB          = trackAB.P(); // momentum of mother particle
                            Float_t BetaAB = TMath::Sqrt(1./(1+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            dirY = dirY/dirY.mag();
                            Double_t scalarProduct = dirY.dot(vectorABtoPrim/vectorABtoPrim.mag());

                            baseY = vectorAB;
                            Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vectorprim);


                            // Apply Lambda cuts
                            if(
                               InvMassAB         > 1.11
                               && InvMassAB      < 1.12
                               && VerdistX       > VerdistX_cut // 4.0
                               && dcaAB_f        < 1.0  // 1.0
                               && VerdistY       > 0.2  // 0.2
                               && scalarProduct  > 0.0
                              )
                            {
                                for(Int_t par_Xi_Omega = 0; par_Xi_Omega < 2; par_Xi_Omega++) // First loop is Xi, second is Omega
                                {
                                    if(ParticleA == 14 && ParticleB == 9 && par_Xi_Omega == 0) ParticleC = 9;  // Xi-
                                    if(ParticleA == 14 && ParticleB == 9 && par_Xi_Omega == 1) ParticleC = 12; // Omega-
                                    if(ParticleA == 15 && ParticleB == 8 && par_Xi_Omega == 0) ParticleC = 8;  // Xi+
                                    if(ParticleA == 15 && ParticleB == 8 && par_Xi_Omega == 1) ParticleC = 11; // Omega+
                                    //****************************************************************************************************************
                                    // Omega or Xi loop
                                    for(Int_t k = 0; k < PID_counter_Array_B[Ana_Num][ParticleC]; k++)  // K-/+, pi-/+ candidates
                                    {
                                        Int_t trackC_num = PID_Array_B[Ana_Num][ParticleC][k];
                                        if(
                                           (
                                            trackC_num != trackB_num
                                            && trackC_num != trackA_num)
                                           || (SE_ME_Flag == 1)
                                          )
                                        {
                                            StPhysicalHelixD helixC;
                                            StThreeVectorF vectorC, vectoratsC, vectoratsA2, vectorA2C, vectorA2CtoPrim, vectornewC, dirY2, baseY2;

                                            StPicoAlexTrack trackC  = *event_C_ana->track( trackC_num ); // take now event B
                                            Float_t MomentumC   = trackC.gMom().mag();
                                            Float_t dcaC        = trackC.dca();   // distance of closest approach to primary vertex
                                            Float_t nHitsPossC  = trackC.nHitsMax();
                                            Float_t nHitsFitC   = trackC.nHitsFit();
                                            //Float_t PolarityC   = trackC.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                            Float_t BetaC       = trackC.btofBeta();  // Velocity after time-of-flight reconstruction
                                            Float_t nSigmaC     = -100.0;
                                            if(ParticleC == 8 || ParticleC == 9)   nSigmaC = trackC.nSigmaPion(); // for Xi analysis
                                            if(ParticleC == 11 || ParticleC == 12) nSigmaC = trackC.nSigmaKaon(); // for Omega analysis
                                            //Float_t TofC        = trackC.btof();
                                            //Float_t TPCdEdxC    = trackC.dEdx(); // TPC dE/dx
                                            Float_t Mass2C      = -100.0;
                                            // calculate mass2
                                            if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                                            {
                                                flag_tof_particleC = 1;
                                                Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                                            }
                                            else
                                            {
                                                flag_tof_particleC = 0;
                                            }
                                            Float_t MassC       = SquareRoot(Mass2C);

                                            if(
                                               dcaC          > dcaC_cut  // 0.4
                                               && nHitsFitC  > 14
                                               && nHitsPossC > 0
                                               && MomentumC  > 0.1
                                               && MomentumC  < 10.0
                                               && (nHitsFitC/nHitsPossC) > 0.52
                                               && (
                                                   flag_tof_particleC == 0 ||
                                                   ((ParticleC == 8  || ParticleC == 9)  && (flag_tof_particleC == 1 && Mass2C > -0.1 && Mass2C < 0.1)) ||
                                                   ((ParticleC == 11 || ParticleC == 12) && (flag_tof_particleC == 1 && Mass2C > 0.125 && Mass2C < 0.36))
                                                  )
                                              )
                                            {
                                                helixC = StPhysicalHelixD(trackC.gMom(),trackC.origin()+vectordiff,event_C_ana->bField()*MAGFIELDFACTOR,trackC.charge());

                                                Float_t pathA2_f, pathC_f,dcaBC_f;
                                                fHelixAtoLinedca(dirY,baseY,helixC,pathC_f,pathA2_f,dcaBC_f); // calculates dca of second pi- helix to Lambda track

                                                vectoratsA2       = baseY+pathA2_f*dirY;  // space vector of Lambda at dca to helixC
                                                vectoratsC        = helixC.at(pathC_f);  // space vector of helixC at dca to Lambda
                                                vectorA2C         = vectoratsA2+vectoratsC;
                                                vectorA2C         = vectorA2C/2.0; // decay vertex
                                                vectorA2CtoPrim   = vectorA2C - vectorprim; // vector primary vertex to decay vertex
                                                Float_t VerdistX2 = vectorA2CtoPrim.mag(); // distance between primary vertex and decay vertex

                                                vectornewC       = helixC.cat(pathC_f); // direction vector at dca for helixC
                                                vectornewC       = MomentumC*vectornewC/vectornewC.mag(); // new momentum vector at decay vertex

                                                if(ParticleC == 8 || ParticleC == 9) // Xi analysis
                                                {
                                                    ltrackC.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.13957018);
                                                    ltrackC_pi.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.493677);
                                                }
                                                if(ParticleC == 11 || ParticleC == 12) // Omega analysis
                                                {
                                                    ltrackC.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.493677);
                                                    ltrackC_pi.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.13957018);
                                                }

                                                // Missing mass and invariant mass calculations
                                                TLorentzVector trackAB_Lambda; // mother particle with nominal Lambda mass
                                                trackAB_Lambda.SetXYZM(trackAB.Px(),trackAB.Py(),trackAB.Pz(),1.115683);
                                                TLorentzVector trackLambdaC    = trackAB_Lambda+ltrackC; // mother particle
                                                TLorentzVector trackLambdaC_pi = trackAB_Lambda+ltrackC_pi; // mother particle

                                                Double_t InvMassABC      = trackLambdaC.M(); // invariant mass of mother particle
                                                Double_t InvMassABC_Xi   = trackLambdaC_pi.M(); // invariant mass of mother particle
                                                Double_t MomentumABC     = trackLambdaC.P(); // momentum of mother particle
                                                Float_t BetaABC = TMath::Sqrt(1./(1+(InvMassABC/MomentumABC)*(InvMassABC/MomentumABC)));

                                                dirY2.set(trackLambdaC.Px(),trackLambdaC.Py(),trackLambdaC.Pz()); // direction vector of Xi
                                                dirY2 = dirY2/dirY2.mag();
                                                Double_t scalarProduct2 = dirY2.dot(vectorA2CtoPrim/vectorA2CtoPrim.mag());

                                                baseY2 = vectorA2C;
                                                //Double_t  VerdistY2  = calculateMinimumDistanceStraightToPoint(baseY2,dirY2,vectorprim);

                                                Float_t pt2          = trackLambdaC.Pt();  // Transverse momentum of mother particle
                                                Float_t rap2         = trackLambdaC.Rapidity(); // Rapidity of mother particle

                                                Float_t phiABC   = dirY2.phi();
                                                Float_t thetaABC = dirY2.theta();


                                                //********************** Calculate the helix of the Omega/Xi ***************************
                                                StPhysicalHelixD helix_omega;
                                                StThreeVectorF origin_omega;   // decay vertex of Omega
                                                origin_omega.set(vectorA2C.x(),vectorA2C.y(),vectorA2C.z());
                                                //Float_t h_omega     = trackC.geth(); // +1 for K+ and -1 for K-
                                                Float_t h_omega;
                                                if(ParticleC == 12 || ParticleC == 9) h_omega = -1.0;
                                                if(ParticleC == 11 || ParticleC == 8) h_omega = 1.0;
                                                Float_t phase_omega = TMath::ATan2(-1.0*h_omega*trackLambdaC.Px(),h_omega*trackLambdaC.Py());
                                                Float_t dip_omega   = TMath::ATan2(trackLambdaC.Pz(),pt2);  // correct
                                                Float_t curv_omega  = 1.0;
                                                if(pt2 != 0.0)
                                                {
                                                    curv_omega = curv_to_invpt_ratio/pt2;
                                                }

                                                helix_omega.setParameters(curv_omega,dip_omega,phase_omega,origin_omega,h_omega);
                                                Float_t path_omega = -999.0;
                                                Float_t dca_omega  = -999.0;
                                                fHelixAtoPointdca(vectorprim,helix_omega,path_omega,dca_omega);
                                                //**************************************************************************************

                                                if(
                                                   InvMassABC        < InvMassABC_cut // 2.5
                                                   && dca_omega      < dca_mother_cut  // 0.5
                                                   && dcaBC_f        < dcaBC_cut  // 0.8
                                                   && VerdistX2      > VerdistX2_cut  // 2.5
                                                   && scalarProduct2 > 0.0
                                                  )
                                                {

                                                    //****************************** BETA CORRECTION OMEGA ******************************
                                                    Float_t BetaACorr   = BetaA;
                                                    Float_t Mass2ACorr  = Mass2A;
                                                    Float_t BetaBCorr   = BetaB;
                                                    Float_t Mass2BCorr  = Mass2B;
                                                    Float_t BetaCCorr   = BetaC;
                                                    Float_t Mass2CCorr  = Mass2C;
                                                    Float_t MassACorr   = MassA;
                                                    Float_t MassBCorr   = MassB;
                                                    Float_t MassCCorr   = MassC;
                                                    Float_t PathAB      = (vectorAB-vectorA2C).mag(); // lambda uncharged
                                                    Float_t PathABC     = vectorA2CtoPrim.mag(); // Xi charged, but curvature negligible

                                                    if( flag_tof_particleA == 1 || flag_tof_particleB == 1 || flag_tof_particleC == 1 )
                                                    {
                                                        if( debug_flag )
                                                        {
                                                            cout << "----------------------------------------------- OMEGA -----------------------------------------------------------------------" << endl;
                                                            cout << "PathAB = " << PathAB << "\tInvMassAB = " << InvMassAB << "\tMomentumAB = " << MomentumAB
                                                                << "\tBetaAB = " << BetaAB << "\tTofAB = " << PathAB / (BetaAB*clight) << endl;
                                                            cout << "PathABC = " << PathABC << "\tInvMassABC = " << InvMassABC << "\tMomentumABC = " << MomentumABC
                                                                << "\tBetaABC = " << BetaABC << "\tTofABC = " << PathABC / (BetaABC*clight) << endl;
                                                        }
                                                    }
                                                    if( flag_tof_particleA == 1 )
                                                    {
                                                        // proton
                                                        BetaACorr = correctBeta4DoubleDecay(trackA,trackAB_Lambda,trackLambdaC,helixA,vectorAB,PathAB,PathABC);
                                                        Mass2ACorr = MomentumA * MomentumA * (1./(BetaACorr*BetaACorr)-1.);
                                                        MassACorr = SquareRoot(Mass2ACorr);
                                                        if ( debug_flag ) cout << "A) MassA = " << MassA << "\tMassACorr = " <<  MassACorr << endl;
                                                    }
                                                    if( flag_tof_particleB == 1 )
                                                    {
                                                        // pi-
                                                        BetaBCorr = correctBeta4DoubleDecay(trackB,trackAB_Lambda,trackLambdaC,helixB,vectorAB,PathAB,PathABC);
                                                        Mass2BCorr = MomentumB * MomentumB * (1./(BetaBCorr*BetaBCorr)-1.);
                                                        MassBCorr = SquareRoot(Mass2BCorr);
                                                        if ( debug_flag ) cout << "B) MassB = " << MassB << "\tMassBCorr = " <<  MassBCorr << endl;
                                                    }
                                                    if( flag_tof_particleC == 1 )
                                                    {
                                                        // pi- (C)
                                                        BetaCCorr = correctBeta4SingleDecay(trackC,trackLambdaC,helixC,vectorA2C,PathABC);
                                                        Mass2CCorr = MomentumC * MomentumC * (1./(BetaCCorr*BetaCCorr)-1.);
                                                        MassCCorr = SquareRoot(Mass2CCorr);
                                                        if ( debug_flag ) cout << "C) MassC = " << MassC << "\tMassCCorr = " <<  MassCCorr << endl;
                                                    }


                                                    if(flag_tof_particleA == 0) Mass2ACorr  = -100;
                                                    if(flag_tof_particleB == 0) Mass2BCorr  = -100;
                                                    if(flag_tof_particleC == 0) Mass2CCorr  = -100;

                                                    //****************************************************************************************


                                                    //****************** Event plane calculations ***********************************************************
                                                    // Calculate the same vectors as for the event plane -> no V0!
                                                    Float_t pathA_prim,dcaA_prim,pathB_prim,dcaB_prim,pathC_prim,dcaC_prim;
                                                    fHelixAtoPointdca(vectorprim,helixA,pathA_prim,dcaA_prim);
                                                    fHelixAtoPointdca(vectorprim,helixB,pathB_prim,dcaB_prim);
                                                    fHelixAtoPointdca(vectorprim,helixC,pathC_prim,dcaC_prim);

                                                    StThreeVectorF vectornewA_prim,vectornewB_prim,vectornewC_prim;
                                                    vectornewA_prim     = helixA.cat(pathA_prim);
                                                    vectornewB_prim     = helixB.cat(pathB_prim);
                                                    vectornewC_prim     = helixC.cat(pathC_prim);

                                                    vectornewA_prim     = MomentumA*vectornewA_prim/vectornewA_prim.mag(); // momentum vector at primary vertex
                                                    vectornewB_prim     = MomentumB*vectornewB_prim/vectornewB_prim.mag(); // momentum vector at primary vertex
                                                    vectornewC_prim     = MomentumC*vectornewC_prim/vectornewC_prim.mag(); // momentum vector at primary vertex

                                                    Float_t phi_event_plane                = -400.0;
                                                    Float_t phi_event_plane_eta_gap        = -400.0;
                                                    Float_t delta_phi_ME_AB_weight         = 0.0;
                                                    Float_t delta_phi_ME_AB_weight_eta_gap = 0.0;

                                                    // Calculate the event plane anlges
                                                    calc_event_plane_angles(3,rap2,trackA,trackB,trackC,vectornewA_prim,vectornewB_prim,vectornewC_prim,ME_Flag,SE_ME_Flag,RunIdA,EventVertexXA,EventVertexYA,
                                                                            EventVertexZA,RunIdB,EventVertexXB,EventVertexYB,EventVertexZB,
                                                                            phi_event_plane,phi_event_plane_eta_gap,delta_phi_ME_AB_weight,delta_phi_ME_AB_weight_eta_gap);


                                                    Int_t delta_phi_ME = 0;
                                                    // check whether the event planes between event A and B are close to each other
                                                    if(
                                                       fabs(delta_phi_ME_AB_weight) < TMath::DegToRad()*30.0
                                                       || (SE_ME_Flag == 0)
                                                      )
                                                    {
                                                        delta_phi_ME = 1;
                                                    }


                                                    //******************************************************************************************************
                                                    delta_phi_ME = 1; // ignore the delta_phi_ME

                                                    if(delta_phi_ME == 1)
                                                    {
                                                        Float_t p_xA_c   = vectornewA_prim.x();
                                                        Float_t p_yA_c   = vectornewA_prim.y();
                                                        Float_t p_tA_c   = sqrt(p_xA_c*p_xA_c + p_yA_c*p_yA_c);
                                                        Float_t etaA_c   = vectornewA_prim.pseudoRapidity();
                                                        Float_t phiA_c   = vectornewA_prim.phi();
                                                        Double_t p_t_weightA = 1.0;
                                                        if(p_tA_c < 2.0)  p_t_weightA = p_tA_c;
                                                        if(p_tA_c >= 2.0) p_t_weightA = 2.0;
                                                        Float_t iQxA     = p_t_weightA*TMath::Cos(2.0*phiA_c);
                                                        Float_t iQyA     = p_t_weightA*TMath::Sin(2.0*phiA_c);

                                                        Float_t p_xB_c   = vectornewB_prim.x();
                                                        Float_t p_yB_c   = vectornewB_prim.y();
                                                        Float_t p_tB_c   = sqrt(p_xB_c*p_xB_c + p_yB_c*p_yB_c);
                                                        Float_t etaB_c   = vectornewB_prim.pseudoRapidity();
                                                        Float_t phiB_c   = vectornewB_prim.phi();
                                                        Double_t p_t_weightB = 1.0;
                                                        if(p_tB_c < 2.0)  p_t_weightB = p_tB_c;
                                                        if(p_tB_c >= 2.0) p_t_weightB = 2.0;
                                                        Float_t iQxB     = p_t_weightB*TMath::Cos(2.0*phiB_c);
                                                        Float_t iQyB     = p_t_weightB*TMath::Sin(2.0*phiB_c);

                                                        Float_t p_xC_c   = vectornewC_prim.x();
                                                        Float_t p_yC_c   = vectornewC_prim.y();
                                                        Float_t p_tC_c   = sqrt(p_xC_c*p_xC_c + p_yC_c*p_yC_c);
                                                        Float_t etaC_c   = vectornewC_prim.pseudoRapidity();
                                                        Float_t phiC_c   = vectornewC_prim.phi();
                                                        Double_t p_t_weightC = 1.0;
                                                        if(p_tC_c < 2.0)  p_t_weightC = p_tC_c;
                                                        if(p_tC_c >= 2.0) p_t_weightC = 2.0;
                                                        Float_t iQxC     = p_t_weightC*TMath::Cos(2.0*phiC_c);
                                                        Float_t iQyC     = p_t_weightC*TMath::Sin(2.0*phiC_c);


                                                        if(ParticleC == 11 || ParticleC == 12) // K+/- --> Omega-/+ analysis
                                                        {
                                                            alexV0_track_A = alexV0_event_A.createTrack();
                                                            alexV0_track_A->setm2A(Mass2ACorr);
                                                            alexV0_track_A->setm2B(Mass2BCorr);
                                                            alexV0_track_A->setm2C(Mass2CCorr);
                                                            alexV0_track_A->setnsA(nSigmaPA);
                                                            alexV0_track_A->setnsB(nSigmaPiB);
                                                            alexV0_track_A->setnsC(nSigmaC);
                                                            alexV0_track_A->setdcaA(dcaA);
                                                            alexV0_track_A->setdcaB(dcaB);
                                                            alexV0_track_A->setdcaC(dcaC);
                                                            alexV0_track_A->setiQxA(iQxA);
                                                            alexV0_track_A->setiQyA(iQyA);
                                                            alexV0_track_A->setiQxB(iQxB);
                                                            alexV0_track_A->setiQyB(iQyB);
                                                            alexV0_track_A->setiQxC(iQxC);
                                                            alexV0_track_A->setiQyC(iQyC);
                                                            alexV0_track_A->setetaA(etaA_c);
                                                            alexV0_track_A->setetaB(etaB_c);
                                                            alexV0_track_A->setetaC(etaC_c);
                                                            alexV0_track_A->setInvAB(InvMassAB);
                                                            alexV0_track_A->setInvABC(InvMassABC);
                                                            alexV0_track_A->setInvAB_miss(InvMassAB_K0S);
                                                            alexV0_track_A->setInvABC_miss(InvMassABC_Xi);
                                                            alexV0_track_A->setdcaAB(dcaAB_f);
                                                            alexV0_track_A->setdcaBC(dcaBC_f);
                                                            alexV0_track_A->setdcaABC(dip_omega);
                                                            alexV0_track_A->setVerdistX(VerdistX);
                                                            alexV0_track_A->setVerdistY(VerdistY);
                                                            alexV0_track_A->setVerdistX2(VerdistX2);
                                                            alexV0_track_A->setVerdistY2(dca_omega);
                                                            alexV0_track_A->setpt(pt2);
                                                            alexV0_track_A->setrap(rap2);
                                                            alexV0_track_A->setphi(phiABC);
                                                            alexV0_track_A->settheta(thetaABC);
                                                            alexV0_track_A->setPsi_ep(phi_event_plane);
                                                            alexV0_track_A->setPsi_ep_eta(phi_event_plane_eta_gap);
                                                            alexV0_track_A->setPsi_diff_ME(delta_phi_ME_AB_weight);
                                                            alexV0_track_A->setscal_prod(scalarProduct);
                                                            alexV0_track_A->setscal_prod2(scalarProduct2);
                                                        }

                                                        if(ParticleC == 8 || ParticleC == 9) // pi+/- --> Xi-/+ analysis
                                                        {
                                                            alexV0_track_B = alexV0_event_B.createTrack();
                                                            alexV0_track_B->setm2A(Mass2ACorr);
                                                            alexV0_track_B->setm2B(Mass2BCorr);
                                                            alexV0_track_B->setm2C(Mass2CCorr);
                                                            alexV0_track_B->setnsA(nSigmaPA);
                                                            alexV0_track_B->setnsB(nSigmaPiB);
                                                            alexV0_track_B->setnsC(nSigmaC);
                                                            alexV0_track_B->setdcaA(dcaA);
                                                            alexV0_track_B->setdcaB(dcaB);
                                                            alexV0_track_B->setdcaC(dcaC);
                                                            alexV0_track_B->setiQxA(iQxA);
                                                            alexV0_track_B->setiQyA(iQyA);
                                                            alexV0_track_B->setiQxB(iQxB);
                                                            alexV0_track_B->setiQyB(iQyB);
                                                            alexV0_track_B->setiQxC(iQxC);
                                                            alexV0_track_B->setiQyC(iQyC);
                                                            alexV0_track_B->setetaA(etaA_c);
                                                            alexV0_track_B->setetaB(etaB_c);
                                                            alexV0_track_B->setetaC(etaC_c);
                                                            alexV0_track_B->setInvAB(InvMassAB);
                                                            alexV0_track_B->setInvABC(InvMassABC);
                                                            alexV0_track_B->setInvAB_miss(InvMassAB_K0S);
                                                            alexV0_track_B->setInvABC_miss(InvMassABC_Xi);
                                                            alexV0_track_B->setdcaAB(dcaAB_f);
                                                            alexV0_track_B->setdcaBC(dcaBC_f);
                                                            alexV0_track_B->setdcaABC(dip_omega);
                                                            alexV0_track_B->setVerdistX(VerdistX);
                                                            alexV0_track_B->setVerdistY(VerdistY);
                                                            alexV0_track_B->setVerdistX2(VerdistX2);
                                                            alexV0_track_B->setVerdistY2(dca_omega);
                                                            alexV0_track_B->setpt(pt2);
                                                            alexV0_track_B->setrap(rap2);
                                                            alexV0_track_B->setphi(phiABC);
                                                            alexV0_track_B->settheta(thetaABC);
                                                            alexV0_track_B->setPsi_ep(phi_event_plane);
                                                            alexV0_track_B->setPsi_ep_eta(phi_event_plane_eta_gap);
                                                            alexV0_track_B->setPsi_diff_ME(delta_phi_ME_AB_weight);
                                                            alexV0_track_B->setscal_prod(scalarProduct);
                                                            alexV0_track_B->setscal_prod2(scalarProduct2);
                                                        }

                                                        dummy_counter++;
                                                        dummy_counter_loop++;

                                                    }
                                                }
                                            }
                                        }
                                    }

                                    //****************************************************************************************************************
                                }
                            }
                        }
                    }
                }
            }
        }
        if(ParticleA == 14 && ParticleB == 9)
        {
            Tree_OmegaMV0_v2  ->Fill();
            Tree_XiMV0_v2     ->Fill();
        }
        if(ParticleA == 15 && ParticleB == 8)
        {
            Tree_OmegaPV0_v2  ->Fill();
            Tree_XiPV0_v2     ->Fill();
        }
        return 1;
    }
    else return 0;
}



Int_t LambdaCPlus_V0_analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs],
                  Int_t PID_Array[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                  StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                  Int_t ParticleA, Int_t ParticleB, Int_t ParticleC,Int_t Ana_Num,Int_t SE_ME_Flag)
{
    // LambdaC+
    // mass = 2286.46 MeV/c2
    // ctau = 59.9 mum
    //
    // Au + Au -> LambdaC+ + N + N
    //             |
    //             -> Lambda + pi+
    //                  |
    //                  -> p + pi-


    // Au + Au -> LambdaC+ + N + N
    //             |
    //             -> K0S + pi+
    //                 |
    //                 -> pi+ + pi-

    dummy_counter_A++;
    // Event vertex information
    StThreeVectorF vectorprim,vectorprimB,vectordiff, vectornewA_lin, vectornewB_lin;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB,RunIdA,RunIdB;


    event_A_ana       = picoDst_A;
    event_B_ana       = picoDst_A;
    event_C_ana       = picoDst_B;

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();

    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    refMultB          = refMultA;
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_A_ana->vzVpd();

    vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    Int_t flag_tof_particleA = 0;
    Int_t flag_tof_particleB = 0;
    Int_t flag_tof_particleC = 0;


    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        EventVertexXB  = event_C_ana->primaryVertex().x();
        EventVertexYB  = event_C_ana->primaryVertex().y();
        EventVertexZB  = event_C_ana->primaryVertex().z();
        refMultB       = event_C_ana->refMult();
        RunIdB         = event_C_ana->runId();
        vectorprimB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vectorprim - vectorprimB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }


    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
        //cout << "Triggers checked, Vertex checked..." << endl;
    }

    //
    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array[Ana_Num][ParticleA]      > 0
       && PID_counter_Array[Ana_Num][ParticleB]      > 0
       && PID_counter_Array_B[Ana_Num][ParticleC]    > 0
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana   ->isMinBias()
       && event_C_ana   ->isMinBias()
      )
    {
        // Loop over all particle combinations
        Double_t dcaA_cut           = 0.1;
        Double_t dcaB_cut           = 0.7;
        Double_t InvMassAB_low_cut  = 1.112;
        Double_t InvMassAB_high_cut = 1.1184;
        Double_t mass_mother        = 1.115683;
        Int_t flag_mother_particle  = 0; // 0 = Lambda, 1 = K0S
        if(ParticleA == 14 && ParticleB == 9)
        {
            dcaA_cut           = 0.1;
            dcaB_cut           = 0.7;
            InvMassAB_low_cut  = 1.112;
            InvMassAB_high_cut = 1.1184;
            mass_mother        = 1.115683;
            flag_mother_particle  = 0;
        }
        if(ParticleA == 8 && ParticleB == 9)
        {
            dcaA_cut           = 0.7;
            dcaB_cut           = 0.7;
            InvMassAB_low_cut  = 0.487;
            InvMassAB_high_cut = 0.506;
            mass_mother        = 0.497648;
            flag_mother_particle  = 1;
        }


        StPhysicalHelixD helixA, helixB;
        StThreeVectorF vectorA, vectorB, vectoratsA, vectoratsB, vectorAB, vectorAB_est, vectorprimAB, vectornewA, vectornewB, dirY_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackC_pi, ltrackB2, ltrackD, ltrackA_lin, ltrackB_lin, ltrackB2_pi, ltrackA_pip;

        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++)
        {
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];

            StPicoAlexTrack trackA  = *event_A_ana->track( trackA_num );
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t nSigmaA     = trackA.nSigmaProton();
            Float_t TofA        = trackA.btof();
            if(ParticleA == 8  || ParticleA == 9)  nSigmaA = trackA.nSigmaPion();
            if(ParticleA == 11 || ParticleA == 12) nSigmaA = trackA.nSigmaKaon();
            if(ParticleA == 14 || ParticleA == 15) nSigmaA = trackA.nSigmaProton();

            Float_t Mass2A      = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                flag_tof_particleA = 1;
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }
            else
            {
                flag_tof_particleA = 0;
            }
            Float_t MassA       = SquareRoot(Mass2A);

            if(
               dcaA          > dcaA_cut  // 0.15
               && nHitsFitA  > 14
               && nHitsPossA > 0.0
               && MomentumA  > 0.1
               && MomentumA  < 10.0
               && (nHitsFitA/nHitsPossA) > 0.52
               && (flag_tof_particleA == 0 || (flag_tof_particleA == 1 && Mass2A > Mass2_low_cut[ParticleA] && Mass2A < Mass2_high_cut[ParticleA]))
               // && fabs(nSigmaA) < 2.5
              )
            {
                helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());

                for(Int_t j = 0; j < PID_counter_Array[Ana_Num][ParticleB]; j++)
                {
                    Int_t trackB_num = PID_Array[Ana_Num][ParticleB][j];
                    StPicoAlexTrack trackB  = *event_B_ana->track( trackB_num ); // take again event A
                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                    Float_t nSigmaB     = trackB.nSigmaPion();
                    Float_t TofB        = trackB.btof();
                    if(ParticleB == 8  || ParticleB == 9)  nSigmaB = trackB.nSigmaPion();
                    if(ParticleB == 11 || ParticleB == 12) nSigmaB = trackB.nSigmaKaon();
                    if(ParticleB == 14 || ParticleB == 15) nSigmaB = trackB.nSigmaProton();

                    Float_t Mass2B      = -100.0;
                    // calculate mass2
                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        flag_tof_particleB = 1;
                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                    }
                    else
                    {
                        flag_tof_particleB = 0;
                    }
                    Float_t MassB       = SquareRoot(Mass2B);

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       && dcaB       > dcaB_cut
                       && nHitsFitB  > 14
                       && nHitsPossB > 0.0
                       && MomentumB  > 0.1
                       && MomentumB  < 10.0
                       && (nHitsFitB/nHitsPossB) > 0.52
                       && (flag_tof_particleB == 0 || (flag_tof_particleB == 1 && Mass2B > Mass2_low_cut[ParticleB] && Mass2B < Mass2_high_cut[ParticleB]))
                      )
                    {
                        helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin(),event_B_ana->bField()*MAGFIELDFACTOR, trackB.charge());

                        //cout << "i = " << i << ", j = " << j << endl;
                        Float_t pathA_f, pathB_f, dcaAB_f;
                        Float_t pathA_test = 0.0;
                        Float_t pathB_test = 0.0;
                        Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);

                        vectoratsA     = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
                        vectoratsB     = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
                        vectorAB       = vectoratsA+vectoratsB;
                        vectorAB       = vectorAB/2.0; // decay vertex


                        StThreeVectorF baseA,dirA,baseB,dirB;
                        baseA = helixA.at(pathA_test);
                        baseB = helixB.at(pathB_test);
                        dirA  = helixA.at(pathA_test-2.0) - helixA.at(pathA_test+2.0);
                        dirB  = helixB.at(pathB_test-2.0) - helixB.at(pathB_test+2.0);

                        StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

                        Double_t dcaAB_lin = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
                        StThreeVectorF vectorABtoPrim_lin = vectorAB_lin - vectorprim; // vector primary vertex to decay vertex
                        Float_t VerdistX_lin = vectorABtoPrim_lin.mag(); // distance between primary vertex and decay verte


                        // calculate the scalar product with the approximated secondary vertex position
                        vectornewA_lin     = helixA.cat(pathA_test); // direction vector at dca for helixA
                        vectornewB_lin     = helixB.cat(pathB_test); // direction vector at dca for helixB
                        vectornewA_lin = MomentumA*vectornewA_lin/vectornewA_lin.mag(); // new momentum vector at decay vertex
                        vectornewB_lin = MomentumB*vectornewB_lin/vectornewB_lin.mag(); // new momentum vector at decay vertex
                        ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),mass_array[ParticleA]);
                        ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),mass_array[ParticleB]);
                        TLorentzVector trackAB_lin      = ltrackA_lin+ltrackB_lin; // mother particle
                        dirY_lin.set(trackAB_lin.Px(),trackAB_lin.Py(),trackAB_lin.Pz());
                        dirY_lin = dirY_lin/dirY_lin.mag();
                        Double_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim_lin/vectorABtoPrim_lin.mag());


                        if( VerdistX_lin > 1.8 && dcaAB_lin < 2.0 && scalarProduct_lin > 0.0 )
                        {
                            if(fDCA_Helix_out == 1)
                            {
                                fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
                            }
                            else
                            {
                                fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
                            }

                            vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                            vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                            vectorAB       = vectoratsA+vectoratsB;
                            vectorAB       = vectorAB/2.0; // decay vertex
                            vectorABtoPrim = vectorAB - vectorprim; // vector primary vertex to decay vertex
                            Float_t VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

                            vectornewA     = helixA.cat(pathA_f); // direction vector at dca for helixA
                            vectornewB     = helixB.cat(pathB_f); // direction vector at dca for helixB

                            vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
                            vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex

                            ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),mass_array[ParticleA]);
                            ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),mass_array[ParticleB]);

                            // Invariant mass calculations
                            TLorentzVector trackAB      = ltrackA+ltrackB; // mother particle
                            Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                            Float_t MomentumAB          = trackAB.P(); // momentum of mother particle
                            Float_t BetaAB = TMath::Sqrt(1./(1+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            dirY = dirY/dirY.mag();
                            Double_t scalarProduct = dirY.dot(vectorABtoPrim/vectorABtoPrim.mag());

                            baseY = vectorAB;
                            Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vectorprim);



                            //-----------------------------------------------------------
                            // Apply Patrick's time-of-flight correction for V0 particles
                            StThreeVectorD tofhitA = trackA.btofHisPos();
                            StThreeVectorD tofhitB = trackB.btofHisPos();

                            StThreeVectorD momGlob(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            StLorentzVectorD SttrackAB(momGlob,momGlob.massHypothesis(mass_mother));

                            Float_t BetaA_corr  = BetaA;
                            Float_t BetaB_corr  = BetaB;

                            Float_t Mass2A_corr = Mass2A;
                            Float_t Mass2B_corr = Mass2B;


                            if(trackA.btofMatchFlag() > 0)
                            {
                                V0_tof_corr->setVectors3D(vectorprim)(vectorAB)(tofhitA);
                                V0_tof_corr->setMotherTracks(SttrackAB);
                                V0_tof_corr->correctBeta(helixA,TofA,BetaA_corr);
                                //Mass2A_corr = MomentumA*MomentumA*(1.0/(BetaA_corr*BetaA_corr) - 1.0);
                                Mass2A_corr = V0_tof_corr->calcM2(MomentumA,BetaA_corr);
                                V0_tof_corr->clearContainers();
                            }
                            if(trackB.btofMatchFlag() > 0)
                            {
                                V0_tof_corr->setVectors3D(vectorprim)(vectorAB)(tofhitB);
                                V0_tof_corr->setMotherTracks(SttrackAB);
                                V0_tof_corr->correctBeta(helixB,TofB,BetaB_corr);
                                //Mass2B_corr = MomentumB*MomentumB*(1.0/(BetaB_corr*BetaB_corr) - 1.0);
                                Mass2B_corr = V0_tof_corr->calcM2(MomentumB,BetaB_corr);
                                V0_tof_corr->clearContainers();
                            }
                            //-----------------------------------------------------------


                            // Apply Lambda or K0S cuts
                            if(
                               dcaAB_f           < 1.3
                               && scalarProduct  > 0.0
                               &&
                               (
                                (
                                 flag_mother_particle == 0 &&
                                 (
                                  ( // no particle has a mass
                                   (trackA.btofMatchFlag() <= 0 && trackB.btofMatchFlag() <= 0)
                                   && dcaA          > 0.6
                                   && dcaB          > 1.7
                                   && VerdistX      > 4.0
                                   && VerdistY      < 0.75
                                  )
                                  ||
                                  ( // only pion has a mass
                                   (trackA.btofMatchFlag() <= 0 && trackB.btofMatchFlag() > 0)
                                   && dcaA          > 0.5
                                   && dcaB          > 1.5
                                   && VerdistX      > 3.5
                                   && VerdistY      < 0.75
                                  )
                                  ||
                                  ( // only proton has a mass
                                   (trackA.btofMatchFlag() > 0 && trackB.btofMatchFlag() <= 0)
                                   && dcaA          > 0.15
                                   && dcaB          > 0.8
                                   && VerdistX      > 2.5
                                   && VerdistY      < 1.2
                                  )
                                  ||
                                  ( // both particles have a mass
                                   (trackA.btofMatchFlag() > 0 && trackB.btofMatchFlag() > 0)
                                   && dcaA          > 0.1
                                   && dcaB          > 0.7
                                   && VerdistX      > 2.0
                                   && VerdistY      < 1.3
                                  )
                                 ) // end Lambda cuts
                                )
                                ||
                                (
                                 flag_mother_particle == 1 &&
                                 (
                                  dcaA             > 0.7
                                  && dcaB          > 0.7
                                  && VerdistX      > 3.0
                                  && VerdistY      < 0.8
                                 )
                                )// end K0S cuts
                               )
                              )
                            {
                                //-----------------------------------------------------------
                                // fill Lambda invariant mass spectra

                                h_massA_LambdaCplus      ->Fill(Mass2A);
                                h_massB_LambdaCplus      ->Fill(Mass2B);
                                h_massA_corr_LambdaCplus ->Fill(Mass2A_corr);
                                h_massB_corr_LambdaCplus ->Fill(Mass2B_corr);

                                h_InvMass_AB_LambdaCplus[0]->Fill(InvMassAB);
                                if(
                                   (trackA.btofMatchFlag() <= 0 || (trackA.btofMatchFlag() > 0 && Mass2A_corr > Mass2_low_cut[ParticleA] && Mass2A_corr < Mass2_high_cut[ParticleA])) &&
                                   (trackB.btofMatchFlag() <= 0 || (trackB.btofMatchFlag() > 0 && Mass2B_corr > Mass2_low_cut[ParticleB] && Mass2B_corr < Mass2_high_cut[ParticleB]))
                                  )
                                {
                                    h_InvMass_AB_LambdaCplus[1]->Fill(InvMassAB);
                                }
                                if(
                                   (flag_tof_particleA == 0 || (flag_tof_particleA == 1 && Mass2A > Mass2_low_cut[ParticleA] && Mass2A < Mass2_high_cut[ParticleA])) &&
                                   (flag_tof_particleB == 0 || (flag_tof_particleB == 1 && Mass2B > Mass2_low_cut[ParticleB] && Mass2B < Mass2_high_cut[ParticleB]))
                                  )
                                {
                                    h_InvMass_AB_LambdaCplus[2]->Fill(InvMassAB);
                                }
                                if(
                                   ((trackA.btofMatchFlag() > 0 && Mass2A_corr > Mass2_low_cut[ParticleA] && Mass2A_corr < Mass2_high_cut[ParticleA])) &&
                                   ((trackB.btofMatchFlag() > 0 && Mass2B_corr > Mass2_low_cut[ParticleB] && Mass2B_corr < Mass2_high_cut[ParticleB]))
                                  )
                                {
                                    h_InvMass_AB_LambdaCplus[3]->Fill(InvMassAB);
                                }

                                if(Mass2A > 0.7 && Mass2A < 1.1 && Mass2B > -0.03 && Mass2B < 0.05) h_InvMass_AB_LambdaCplus[4]->Fill(InvMassAB);
                                if(Mass2A_corr > 0.7 && Mass2A_corr < 1.1 && Mass2B_corr > -0.03 && Mass2B_corr < 0.05) h_InvMass_AB_LambdaCplus[5]->Fill(InvMassAB);

                                if(
                                   (trackA.btofMatchFlag() <= 0 || (trackA.btofMatchFlag() > 0 && Mass2A_corr > Mass2_low_cut_strict[ParticleA] && Mass2A_corr < Mass2_high_cut_strict[ParticleA])) &&
                                   (trackB.btofMatchFlag() <= 0 || (trackB.btofMatchFlag() > 0 && Mass2B_corr > Mass2_low_cut_strict[ParticleB] && Mass2B_corr < Mass2_high_cut_strict[ParticleB]))
                                  )
                                {
                                    h_InvMass_AB_LambdaCplus[6]->Fill(InvMassAB);
                                }

                                if(
                                   (flag_tof_particleA == 0 || (flag_tof_particleA == 1 && Mass2A > Mass2_low_cut_strict[ParticleA] && Mass2A < Mass2_high_cut_strict[ParticleA])) &&
                                   (flag_tof_particleB == 0 || (flag_tof_particleB == 1 && Mass2B > Mass2_low_cut_strict[ParticleB] && Mass2B < Mass2_high_cut_strict[ParticleB]))
                                  )
                                {
                                    h_InvMass_AB_LambdaCplus[7]->Fill(InvMassAB);
                                }
                                //-----------------------------------------------------------



                                // Apply invariant mass cuts and strict mass2 cuts on tof-corrected mass2 values
                                if(
                                   InvMassAB      > InvMassAB_low_cut  &&
                                   InvMassAB      < InvMassAB_high_cut &&
                                   (trackA.btofMatchFlag() <= 0 || (trackA.btofMatchFlag() > 0 && Mass2A_corr > Mass2_low_cut_strict[ParticleA] && Mass2A_corr < Mass2_high_cut_strict[ParticleA])) &&
                                   (trackB.btofMatchFlag() <= 0 || (trackB.btofMatchFlag() > 0 && Mass2B_corr > Mass2_low_cut_strict[ParticleB] && Mass2B_corr < Mass2_high_cut_strict[ParticleB]))
                                  )
                                {
                                    for(Int_t k = 0; k < PID_counter_Array_B[Ana_Num][ParticleC]; k++)
                                    {
                                        Int_t trackC_num = PID_Array_B[Ana_Num][ParticleC][k];
                                        if(
                                           (
                                            trackC_num != trackB_num
                                            && trackC_num != trackA_num)
                                           || (SE_ME_Flag == 1)
                                          )
                                        {
                                            StPhysicalHelixD helixC;
                                            StThreeVectorF vectorC, vectoratsC, vectoratsA2, vectorA2C, vectorA2CtoPrim, vectornewC, dirY2, baseY2;

                                            StPicoAlexTrack trackC  = *event_C_ana->track( trackC_num ); // take now event B
                                            Float_t MomentumC   = trackC.pMom().mag();
                                            Float_t dcaC        = trackC.dca();   // distance of closest approach to primary vertex
                                            Float_t nHitsPossC  = trackC.nHitsMax();
                                            Float_t nHitsFitC   = trackC.nHitsFit();
                                            Float_t BetaC       = trackC.btofBeta();  // Velocity after time-of-flight reconstruction
                                            Float_t nSigmaC     = trackC.nSigmaProton();
                                            if(ParticleC == 8  || ParticleC == 9)  nSigmaC = trackC.nSigmaPion();
                                            if(ParticleC == 11 || ParticleC == 12) nSigmaC = trackC.nSigmaKaon();
                                            if(ParticleC == 14 || ParticleC == 15) nSigmaC = trackC.nSigmaProton();

                                            Float_t Mass2C      = -100.0;
                                            // calculate mass2
                                            if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                                            {
                                                flag_tof_particleC = 1;
                                                Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                                            }
                                            else
                                            {
                                                flag_tof_particleC = 0;
                                            }
                                            Float_t MassC       = SquareRoot(Mass2C);

                                            if(
                                               dcaC          < 1.3
                                               && nHitsFitC  > 14
                                               && nHitsPossC > 0
                                               && MomentumC  > 0.1
                                               && MomentumC  < 10.0
                                               && (nHitsFitC/nHitsPossC) > 0.52
                                               && (flag_tof_particleC == 0 || (flag_tof_particleC == 1 && Mass2C > Mass2_low_cut[ParticleC] && Mass2C < Mass2_high_cut[ParticleC]))
                                              )
                                            {
                                                helixC = StPhysicalHelixD(trackC.pMom(),event_C_ana->primaryVertex()+vectordiff,event_C_ana->bField()*MAGFIELDFACTOR,trackC.charge());

                                                StThreeVectorF primdirC = helixC.cat(helixC.pathLength(vectorprim));
                                                primdirC = MomentumC*primdirC/primdirC.mag();
                                                ltrackC.SetXYZM(primdirC.x(),primdirC.y(),primdirC.z(),mass_array[ParticleC]);

                                                // Missing mass and invariant mass calculations
                                                TLorentzVector trackAB_nominal; // mother particle with nominal mass
                                                trackAB_nominal.SetXYZM(trackAB.Px(),trackAB.Py(),trackAB.Pz(),mass_mother);
                                                TLorentzVector trackABC = trackAB_nominal+ltrackC;

                                                Double_t InvMassABC      = trackABC.M(); // invariant mass of mother particle
                                                Double_t MomentumABC     = trackABC.P(); // momentum of mother particle
                                                Float_t BetaABC = TMath::Sqrt(1./(1+(InvMassABC/MomentumABC)*(InvMassABC/MomentumABC)));

                                                dirY2.set(trackABC.Px(),trackABC.Py(),trackABC.Pz()); // direction vector of Xi
                                                dirY2 = dirY2/dirY2.mag();

                                                Float_t pt2          = trackABC.Pt();  // Transverse momentum of mother particle
                                                Float_t rap2         = trackABC.Rapidity(); // Rapidity of mother particle

                                                Float_t phiABC   = dirY2.phi();
                                                Float_t thetaABC = dirY2.theta();


                                                if(
                                                   1 == 1
                                                  )
                                                {
                                                    h_InvMass_ABC_LambdaCplus[0] ->Fill(InvMassABC);
                                                    h_massC_LambdaCplus          ->Fill(Mass2C);
                                                }
                                                if(
                                                   fabs(nSigmaA)    < 2.5
                                                   && fabs(nSigmaB) < 2.5
                                                   && fabs(nSigmaC) < 2.5
                                                  )
                                                {
                                                    h_InvMass_ABC_LambdaCplus[1]->Fill(InvMassABC);
                                                }
                                                if(
                                                   fabs(nSigmaA)    < 2.5
                                                   && fabs(nSigmaB) < 2.5
                                                   && fabs(nSigmaC) < 2.5
                                                   && trackA.btofMatchFlag() > 0
                                                   && trackB.btofMatchFlag() > 0
                                                   && trackC.btofMatchFlag() > 0
                                                  )
                                                {
                                                    h_InvMass_ABC_LambdaCplus[2]->Fill(InvMassABC);
                                                }
                                                if(erefMult_bin >= 0 && erefMult_bin < 9)
                                                {
                                                    h_InvMass_ABC_LambdaCplus[3+erefMult_bin]->Fill(InvMassABC);
                                                    if(pt2 > 0.0 && pt2 < 1.0)
                                                    {
                                                        h_InvMass_ABC_LambdaCplus[12+erefMult_bin]->Fill(InvMassABC);
                                                    }
                                                    if(pt2 >= 1.0 && pt2 < 2.0)
                                                    {
                                                        h_InvMass_ABC_LambdaCplus[21+erefMult_bin]->Fill(InvMassABC);
                                                    }
                                                    if(pt2 >= 2.0 && pt2 < 3.0)
                                                    {
                                                        h_InvMass_ABC_LambdaCplus[30+erefMult_bin]->Fill(InvMassABC);
                                                    }
                                                    if(pt2 >= 3.0 && pt2 < 20.0)
                                                    {
                                                        h_InvMass_ABC_LambdaCplus[39+erefMult_bin]->Fill(InvMassABC);
                                                    }

                                                    if(
                                                       ((MomentumA <= 1.2) || (MomentumA > 1.2 && fabs(nSigmaA) < 2.5)) &&
                                                       ((MomentumB <= 0.8) || (MomentumB > 0.8 && fabs(nSigmaB) < 2.5)) &&
                                                       ((MomentumC <= 0.8) || (MomentumC > 0.8 && trackC.btofMatchFlag() > 0 && fabs(nSigmaC) < 2.5))
                                                      )
                                                    {
                                                        h_InvMass_ABC_LambdaCplus[48+erefMult_bin]->Fill(InvMassABC);
                                                        if(pt2 > 0.0 && pt2 < 1.0)
                                                        {
                                                            h_InvMass_ABC_LambdaCplus[48+9+erefMult_bin]->Fill(InvMassABC);
                                                        }
                                                        if(pt2 >= 1.0 && pt2 < 2.0)
                                                        {
                                                            h_InvMass_ABC_LambdaCplus[48+18+erefMult_bin]->Fill(InvMassABC);
                                                        }
                                                        if(pt2 >= 2.0 && pt2 < 3.0)
                                                        {
                                                            h_InvMass_ABC_LambdaCplus[48+27+erefMult_bin]->Fill(InvMassABC);
                                                        }
                                                        if(pt2 >= 3.0 && pt2 < 20.0)
                                                        {
                                                            h_InvMass_ABC_LambdaCplus[48+36+erefMult_bin]->Fill(InvMassABC);
                                                        }
                                                    }


                                                }

                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return 1;
    }
    else return 0;
}






Int_t Omega2250_analysis(Int_t PID_counter_Array_A[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs],
                  Int_t PID_Array_A[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                  StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                  Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t ParticleD, Int_t ParticleE, Int_t Ana_Num,Int_t SE_ME_Flag)
{
    // Omega(2250)
    // mass = 2252 MeV/c2
    // width = 55 +/- 18 MeV
    //
    // Au + Au -> Omega(2250)- + N + N
    //             |
    //             -> Xi- + pi+(D) + K-(E)
    //                |
    //                -> Lambda + pi-(C)
    //                  |
    //                  -> p(A) + pi-(B)
    //
    // (A,B,C) -> eventA
    // (D,E)   -> eventB

    dummy_counter_A++;
    // Event vertex information
    StThreeVectorF vectorprim,vectorprimB,vectordiff, vectornewA_lin, vectornewB_lin;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB,RunIdA,RunIdB;


    event_A_ana       = picoDst_A;
    event_B_ana       = picoDst_B;

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();

    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    refMultB          = refMultA;
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_A_ana->vzVpd();

    vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    Int_t flag_tof_particleA = 0;
    Int_t flag_tof_particleB = 0;
    Int_t flag_tof_particleC = 0;
    Int_t flag_tof_particleD = 0;
    Int_t flag_tof_particleE = 0;


    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        RunIdB         = event_B_ana->runId();
        vectorprimB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vectorprim - vectorprimB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }


    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
        //cout << "Triggers checked, Vertex checked..." << endl;
    }

    //
    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array_A[Ana_Num][ParticleA]    > 0   // p(A)
       && PID_counter_Array_A[Ana_Num][ParticleB]    > 0   // pi-
       && PID_counter_Array_A[Ana_Num][ParticleC]    > 0   // pi-
       && PID_counter_Array_B[Ana_Num][ParticleD]    > 0   // pi+
       && PID_counter_Array_B[Ana_Num][ParticleE]    > 0   // K-
       //&& event->getNVertices() >= 1
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana   ->isMinBias()
       && event_B_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
        //   << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;
        // Loop over all particle combinations

        Double_t dcaA_cut            = 0.1;
        Double_t dcaB_cut            = 0.9;
        Double_t dcaC_cut            = 0.5;
        Double_t dcaD_cut            = 1.0;
        Double_t dcaE_cut            = 1.0;
        Double_t VerdistX_cut        = 3.5;
        Double_t InvMassABC_cut      = 1.335;
        Double_t InvMassABC_low_cut  = 1.31;
        Double_t dca_mother_cut      = 0.7;
        Double_t dcaBC_cut           = 1.0;
        Double_t VerdistX2_cut       = 2.5;


        StPhysicalHelixD helixA, helixB, helixC, helixD, helixE;
        StThreeVectorF vectorA, vectorB, vectoratsA, vectoratsB, vectorAB, vectorAB_est, vectorprimAB, vectornewA, vectornewB, vectornewD, vectornewE, dirY_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackD, ltrackE, ltrackC_pi, ltrackB2, ltrackA_lin, ltrackB_lin, ltrackB2_pi, ltrackA_pip;

        /*
        // Fill event information for Xi-/+
        alexV0_event_A.clearTrackList();
        alexV0_event_A.setx(EventVertexXA);
        alexV0_event_A.sety(EventVertexYA);
        alexV0_event_A.setz(EventVertexZA);
        alexV0_event_A.setid(RunIdA);
        alexV0_event_A.setmult(refMultA);
        alexV0_event_A.setn_prim(n_primaries);
        alexV0_event_A.setn_non_prim(n_non_primaries);
        alexV0_event_A.setn_tof_prim(n_tofmatch_prim);
        alexV0_event_A.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        alexV0_event_A.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        alexV0_event_A.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        alexV0_event_A.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        alexV0_event_A.setEP_Qx_ptw(EP_Qx_ptw);
        alexV0_event_A.setEP_Qy_ptw(EP_Qy_ptw);
        alexV0_event_A.setQtracks_eta_pos(Qtracks_used_eta_pos);
        alexV0_event_A.setQtracks_eta_neg(Qtracks_used_eta_neg);
        alexV0_event_A.setQtracks_full(Qtracks_used);
        alexV0_event_A.setZDCx(ZDCx);
        alexV0_event_A.setBBCx(BBCx);
        alexV0_event_A.setvzVpd(vzVpd);
        */

        if(0)
        {
            cout << "Omega event accepted, A = " << PID_counter_Array_A[Ana_Num][ParticleA] <<
                ", B = " << PID_counter_Array_A[Ana_Num][ParticleB] <<
                ", C = " << PID_counter_Array_A[Ana_Num][ParticleC] <<
                ", D = " << PID_counter_Array_B[Ana_Num][ParticleD] <<
                ", E = " << PID_counter_Array_B[Ana_Num][ParticleE] <<
                ", erefMult_bin = " << erefMult_bin <<
                endl;
        }

        for(Int_t i = 0; i < PID_counter_Array_A[Ana_Num][ParticleA]; i++) // p candidates
        {
            Int_t trackA_num = PID_Array_A[Ana_Num][ParticleA][i];

            StPicoAlexTrack trackA  = *event_A_ana->track( trackA_num );
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t nSigmaPA    = trackA.nSigmaProton();
            //Float_t TofA        = trackA.btof();
            Float_t Mass2A      = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                flag_tof_particleA = 1;
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }
            else
            {
                flag_tof_particleA = 0;
            }
            Float_t MassA       = SquareRoot(Mass2A);

            if(
               dcaA          > dcaA_cut  // 0.15
               && nHitsFitA  > 14
               && nHitsPossA > 0.0
               && MomentumA  > 0.1
               && MomentumA  < 10.0
               && (nHitsFitA/nHitsPossA) > 0.52
               && (flag_tof_particleA == 0 || (flag_tof_particleA == 1 && Mass2A > 0.4 && Mass2A < 1.5)) // proton mass
              )
            {
                helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());

                for(Int_t j = 0; j < PID_counter_Array_A[Ana_Num][ParticleB]; j++)  // pi- candidates
                {
                    Int_t trackB_num = PID_Array_A[Ana_Num][ParticleB][j];
                    StPicoAlexTrack trackB  = *event_A_ana->track( trackB_num ); // take again event A
                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    //Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                    Float_t nSigmaPiB   = trackB.nSigmaPion();
                    //Float_t TofB        = trackB.btof();
                    Float_t Mass2B      = -100.0;
                    // calculate mass2
                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        flag_tof_particleB = 1;
                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                    }
                    else
                    {
                        flag_tof_particleB = 0;
                    }
                    Float_t MassB       = SquareRoot(Mass2B);

                    //if( flag_tof_particleA == 0 && flag_tof_particleB == 0) // apply more strict cuts if there is no TOF information for the two particles
                    //{
                    //    dcaA_cut  = 0.1;
                    //    dcaB_cut  = 1.5;
                    //    dcaC_cut  = 0.7;
                    //}

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       && dcaB       > dcaB_cut // 1.3
                       && nHitsFitB  > 14
                       && nHitsPossB > 0.0
                       && MomentumB  > 0.1
                       && MomentumB  < 10.0
                       && (nHitsFitB/nHitsPossB) > 0.52
                       && (flag_tof_particleB == 0 || (flag_tof_particleB == 1 && Mass2B > -0.1 && Mass2B < 0.1)) // pion mass
                      )
                    {
                        helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackB.charge());

                        //cout << "i = " << i << ", j = " << j << endl;
                        Float_t pathA_f, pathB_f, dcaAB_f;
                        Float_t pathA_test = 0.0;
                        Float_t pathB_test = 0.0;
                        Int_t fDCA_Helix_out = fDCA_Helix_Estimate(helixA,helixB,pathA_test,pathB_test,dcaAB_f);

                        vectoratsA     = helixA.at(pathA_test);  // space vector of helixA at dca to helixB
                        vectoratsB     = helixB.at(pathB_test);  // space vector of helixB at dca to helixA
                        vectorAB       = vectoratsA+vectoratsB;
                        vectorAB       = vectorAB/2.0; // decay vertex

                        //Double_t vex_est_x = vectorAB.x();
                        //Double_t vex_est_y = vectorAB.y();
                        //Double_t vex_est_z = vectorAB.z();

                        StThreeVectorF baseA,dirA,baseB,dirB;
                        baseA = helixA.at(pathA_test);
                        baseB = helixB.at(pathB_test);
                        dirA  = helixA.at(pathA_test-2.0) - helixA.at(pathA_test+2.0);
                        dirB  = helixB.at(pathB_test-2.0) - helixB.at(pathB_test+2.0);

                        StThreeVectorF vectorAB_lin  = calcVertexAnalytical(baseA,dirA,baseB,dirB); // vertex of the two tracks

                        //Double_t vex_lin_x = vectorAB_lin.x();
                        //Double_t vex_lin_y = vectorAB_lin.y();
                        //Double_t vex_lin_z = vectorAB_lin.z();

                        Double_t dcaAB_lin = calculateMinimumDistance(baseA,dirA,baseB,dirB);       // minimum distance between the two tracks
                        StThreeVectorF vectorABtoPrim_lin = vectorAB_lin - vectorprim; // vector primary vertex to decay vertex
                        Float_t VerdistX_lin = vectorABtoPrim_lin.mag(); // distance between primary vertex and decay verte


                        // calculate the scalar product with the approximated secondary vertex position
                        vectornewA_lin     = helixA.cat(pathA_test); // direction vector at dca for helixA
                        vectornewB_lin     = helixB.cat(pathB_test); // direction vector at dca for helixB
                        vectornewA_lin = MomentumA*vectornewA_lin/vectornewA_lin.mag(); // new momentum vector at decay vertex
                        vectornewB_lin = MomentumB*vectornewB_lin/vectornewB_lin.mag(); // new momentum vector at decay vertex
                        ltrackA_lin.SetXYZM(vectornewA_lin.x(),vectornewA_lin.y(),vectornewA_lin.z(),0.93827203);
                        ltrackB_lin.SetXYZM(vectornewB_lin.x(),vectornewB_lin.y(),vectornewB_lin.z(),0.13957018);
                        TLorentzVector trackAB_lin      = ltrackA_lin+ltrackB_lin; // mother particle
                        dirY_lin.set(trackAB_lin.Px(),trackAB_lin.Py(),trackAB_lin.Pz());
                        dirY_lin = dirY_lin/dirY_lin.mag();
                        Double_t scalarProduct_lin = dirY_lin.dot(vectorABtoPrim_lin/vectorABtoPrim_lin.mag());


                        if( VerdistX_lin > 3.5 && dcaAB_lin < 1.5 && scalarProduct_lin > 0.0 )
                        {
                            if(fDCA_Helix_out == 1)
                            {
                                fHelixABdca_start_params(helixA,helixB,pathA_f,pathB_f,dcaAB_f,pathA_test,pathB_test); // calculate dca between two helices
                            }
                            else
                            {
                                fHelixABdca(helixA,helixB,pathA_f,pathB_f,dcaAB_f); // calculate dca between two helices
                            }

                            vectoratsA     = helixA.at(pathA_f);  // space vector of helixA at dca to helixB
                            vectoratsB     = helixB.at(pathB_f);  // space vector of helixB at dca to helixA
                            vectorAB       = vectoratsA+vectoratsB;
                            vectorAB       = vectorAB/2.0; // decay vertex
                            vectorABtoPrim = vectorAB - vectorprim; // vector primary vertex to decay vertex
                            Float_t VerdistX = vectorABtoPrim.mag(); // distance between primary vertex and decay vertex

                            vectornewA     = helixA.cat(pathA_f); // direction vector at dca for helixA
                            vectornewB     = helixB.cat(pathB_f); // direction vector at dca for helixB

                            vectornewA = MomentumA*vectornewA/vectornewA.mag(); // new momentum vector at decay vertex
                            vectornewB = MomentumB*vectornewB/vectornewB.mag(); // new momentum vector at decay vertex

                            ltrackA.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.93827203);
                            ltrackB.SetXYZM(vectornewB.x(),vectornewB.y(),vectornewB.z(),0.13957018);
                            ltrackA_pip.SetXYZM(vectornewA.x(),vectornewA.y(),vectornewA.z(),0.13957018);

                            TLorentzVector trackAB_K0S  = ltrackA_pip+ltrackB; // mother particle
                            Double_t InvMassAB_K0S      = trackAB_K0S.M(); // invariant mass of mother particle

                            // Missing mass and invariant mass calculations
                            TLorentzVector trackAB      = ltrackA+ltrackB; // mother particle
                            Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                            Float_t MomentumAB          = trackAB.P(); // momentum of mother particle
                            Float_t BetaAB = TMath::Sqrt(1./(1+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            dirY = dirY/dirY.mag();
                            Double_t scalarProduct = dirY.dot(vectorABtoPrim/vectorABtoPrim.mag());

                            baseY = vectorAB;
                            Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vectorprim);


                            // Apply Lambda cuts
                            if(
                               InvMassAB         > 1.11
                               && InvMassAB      < 1.12
                               && VerdistX       > VerdistX_cut // 4.0
                               && dcaAB_f        < 1.0  // 1.0
                               && VerdistY       > 0.2  // 0.2
                               && scalarProduct  > 0.0
                              )
                            {
                                //cout << "Lambda accepted" << endl;
                                // stop loop at 1 -> can be changed in the future
                                for(Int_t par_Xi_Omega = 0; par_Xi_Omega < 1; par_Xi_Omega++) // First loop is Xi, second is Omega
                                {
                                    //if(ParticleA == 14 && ParticleB == 9 && par_Xi_Omega == 0) ParticleC = 9;  // Xi-
                                    //if(ParticleA == 14 && ParticleB == 9 && par_Xi_Omega == 1) ParticleC = 12; // Omega-
                                    //if(ParticleA == 15 && ParticleB == 8 && par_Xi_Omega == 0) ParticleC = 8;  // Xi+
                                    //if(ParticleA == 15 && ParticleB == 8 && par_Xi_Omega == 1) ParticleC = 11; // Omega+
                                    //****************************************************************************************************************
                                    // Omega or Xi loop
                                    for(Int_t k = 0; k < PID_counter_Array_A[Ana_Num][ParticleC]; k++)  // pi-
                                    {
                                        Int_t trackC_num = PID_Array_A[Ana_Num][ParticleC][k];
                                        if(
                                           trackC_num != trackB_num
                                           && trackC_num != trackA_num
                                          )
                                        {
                                            StThreeVectorF vectorC, vectoratsC, vectoratsA2, vectorA2C, vectorA2CtoPrim, vectornewC, dirY2, baseY2;

                                            StPicoAlexTrack trackC  = *event_A_ana->track( trackC_num ); // take now event B
                                            Float_t MomentumC   = trackC.gMom().mag();
                                            Float_t dcaC        = trackC.dca();   // distance of closest approach to primary vertex
                                            Float_t nHitsPossC  = trackC.nHitsMax();
                                            Float_t nHitsFitC   = trackC.nHitsFit();
                                            //Float_t PolarityC   = trackC.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                            Float_t BetaC       = trackC.btofBeta();  // Velocity after time-of-flight reconstruction
                                            Float_t nSigmaC     = -100.0;
                                            if(ParticleC == 8 || ParticleC == 9)   nSigmaC = trackC.nSigmaPion(); // for Xi analysis
                                            if(ParticleC == 11 || ParticleC == 12) nSigmaC = trackC.nSigmaKaon(); // for Omega analysis
                                            //Float_t TofC        = trackC.btof();
                                            //Float_t TPCdEdxC    = trackC.dEdx(); // TPC dE/dx
                                            Float_t Mass2C      = -100.0;
                                            // calculate mass2
                                            if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                                            {
                                                flag_tof_particleC = 1;
                                                Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                                            }
                                            else
                                            {
                                                flag_tof_particleC = 0;
                                            }
                                            Float_t MassC       = SquareRoot(Mass2C);

                                            if(
                                               dcaC          > dcaC_cut  // 0.4
                                               && nHitsFitC  > 14
                                               && nHitsPossC > 0
                                               && MomentumC  > 0.1
                                               && MomentumC  < 10.0
                                               && (nHitsFitC/nHitsPossC) > 0.52
                                               && (
                                                   flag_tof_particleC == 0 ||
                                                   ((ParticleC == 8  || ParticleC == 9)  && (flag_tof_particleC == 1 && Mass2C > -0.1 && Mass2C < 0.1)) ||
                                                   ((ParticleC == 11 || ParticleC == 12) && (flag_tof_particleC == 1 && Mass2C > 0.125 && Mass2C < 0.36))
                                                  )
                                              )
                                            {
                                                helixC = StPhysicalHelixD(trackC.gMom(),trackC.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackC.charge());

                                                Float_t pathA2_f, pathC_f,dcaBC_f;
                                                fHelixAtoLinedca(dirY,baseY,helixC,pathC_f,pathA2_f,dcaBC_f); // calculates dca of second pi- helix to Lambda track

                                                vectoratsA2       = baseY+pathA2_f*dirY;  // space vector of Lambda at dca to helixC
                                                vectoratsC        = helixC.at(pathC_f);  // space vector of helixC at dca to Lambda
                                                vectorA2C         = vectoratsA2+vectoratsC;
                                                vectorA2C         = vectorA2C/2.0; // decay vertex
                                                vectorA2CtoPrim   = vectorA2C - vectorprim; // vector primary vertex to decay vertex
                                                Float_t VerdistX2 = vectorA2CtoPrim.mag(); // distance between primary vertex and decay vertex

                                                vectornewC       = helixC.cat(pathC_f); // direction vector at dca for helixC
                                                vectornewC       = MomentumC*vectornewC/vectornewC.mag(); // new momentum vector at decay vertex

                                                if(ParticleC == 8 || ParticleC == 9) // Xi analysis
                                                {
                                                    ltrackC.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.13957018);
                                                    ltrackC_pi.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.493677);
                                                }
                                                if(ParticleC == 11 || ParticleC == 12) // Omega analysis
                                                {
                                                    ltrackC.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.493677);
                                                    ltrackC_pi.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.13957018);
                                                }

                                                // Missing mass and invariant mass calculations
                                                TLorentzVector trackAB_Lambda; // mother particle with nominal Lambda mass
                                                trackAB_Lambda.SetXYZM(trackAB.Px(),trackAB.Py(),trackAB.Pz(),1.115683);
                                                TLorentzVector trackLambdaC    = trackAB_Lambda+ltrackC; // mother particle
                                                TLorentzVector trackLambdaC_pi = trackAB_Lambda+ltrackC_pi; // mother particle

                                                Double_t InvMassABC      = trackLambdaC.M(); // invariant mass of mother particle
                                                Double_t InvMassABC_Xi   = trackLambdaC_pi.M(); // invariant mass of mother particle
                                                Double_t MomentumABC     = trackLambdaC.P(); // momentum of mother particle
                                                Float_t BetaABC = TMath::Sqrt(1./(1+(InvMassABC/MomentumABC)*(InvMassABC/MomentumABC)));

                                                dirY2.set(trackLambdaC.Px(),trackLambdaC.Py(),trackLambdaC.Pz()); // direction vector of Xi
                                                dirY2 = dirY2/dirY2.mag();
                                                Double_t scalarProduct2 = dirY2.dot(vectorA2CtoPrim/vectorA2CtoPrim.mag());

                                                baseY2 = vectorA2C;
                                                //Double_t  VerdistY2  = calculateMinimumDistanceStraightToPoint(baseY2,dirY2,vectorprim);

                                                Float_t pt2          = trackLambdaC.Pt();  // Transverse momentum of mother particle
                                                Float_t rap2         = trackLambdaC.Rapidity(); // Rapidity of mother particle

                                                Float_t phiABC   = dirY2.phi();
                                                Float_t thetaABC = dirY2.theta();


                                                //********************** Calculate the helix of the Omega/Xi ***************************
                                                StPhysicalHelixD helix_omega;
                                                StThreeVectorF origin_omega;   // decay vertex of Omega
                                                origin_omega.set(vectorA2C.x(),vectorA2C.y(),vectorA2C.z());
                                                //Float_t h_omega     = trackC.geth(); // +1 for K+ and -1 for K-
                                                Float_t h_omega;
                                                if(ParticleC == 12 || ParticleC == 9) h_omega = -1.0;
                                                if(ParticleC == 11 || ParticleC == 8) h_omega = 1.0;
                                                Float_t phase_omega = TMath::ATan2(-1.0*h_omega*trackLambdaC.Px(),h_omega*trackLambdaC.Py());
                                                Float_t dip_omega   = TMath::ATan2(trackLambdaC.Pz(),pt2);  // correct
                                                Float_t curv_omega  = 1.0;
                                                if(pt2 != 0.0)
                                                {
                                                    curv_omega = curv_to_invpt_ratio/pt2;
                                                }

                                                helix_omega.setParameters(curv_omega,dip_omega,phase_omega,origin_omega,h_omega);
                                                Float_t path_omega = -999.0;
                                                Float_t dca_omega  = -999.0;
                                                fHelixAtoPointdca(vectorprim,helix_omega,path_omega,dca_omega);
                                                //**************************************************************************************

                                                if(
                                                   InvMassABC        < InvMassABC_cut // 2.5
                                                   && InvMassABC     > InvMassABC_low_cut // 2.5
                                                   && dca_omega      < dca_mother_cut  // 0.5
                                                   && dcaBC_f        < dcaBC_cut  // 0.8
                                                   && VerdistX2      > VerdistX2_cut  // 2.5
                                                   && scalarProduct2 > 0.0
                                                  )
                                                {

                                                    //cout << "Xi accepted" << endl;

                                                    //****************************** BETA CORRECTION OMEGA ******************************
                                                    Float_t BetaACorr   = BetaA;
                                                    Float_t Mass2ACorr  = Mass2A;
                                                    Float_t BetaBCorr   = BetaB;
                                                    Float_t Mass2BCorr  = Mass2B;
                                                    Float_t BetaCCorr   = BetaC;
                                                    Float_t Mass2CCorr  = Mass2C;
                                                    Float_t MassACorr   = MassA;
                                                    Float_t MassBCorr   = MassB;
                                                    Float_t MassCCorr   = MassC;
                                                    Float_t PathAB      = (vectorAB-vectorA2C).mag(); // lambda uncharged
                                                    Float_t PathABC     = vectorA2CtoPrim.mag(); // Xi charged, but curvature negligible

                                                    if( flag_tof_particleA == 1 || flag_tof_particleB == 1 || flag_tof_particleC == 1 )
                                                    {
                                                        if( debug_flag )
                                                        {
                                                            cout << "----------------------------------------------- OMEGA -----------------------------------------------------------------------" << endl;
                                                            cout << "PathAB = " << PathAB << "\tInvMassAB = " << InvMassAB << "\tMomentumAB = " << MomentumAB
                                                                << "\tBetaAB = " << BetaAB << "\tTofAB = " << PathAB / (BetaAB*clight) << endl;
                                                            cout << "PathABC = " << PathABC << "\tInvMassABC = " << InvMassABC << "\tMomentumABC = " << MomentumABC
                                                                << "\tBetaABC = " << BetaABC << "\tTofABC = " << PathABC / (BetaABC*clight) << endl;
                                                        }
                                                    }
                                                    if( flag_tof_particleA == 1 )
                                                    {
                                                        // proton
                                                        BetaACorr = correctBeta4DoubleDecay(trackA,trackAB_Lambda,trackLambdaC,helixA,vectorAB,PathAB,PathABC);
                                                        Mass2ACorr = MomentumA * MomentumA * (1./(BetaACorr*BetaACorr)-1.);
                                                        MassACorr = SquareRoot(Mass2ACorr);
                                                        if ( debug_flag ) cout << "A) MassA = " << MassA << "\tMassACorr = " <<  MassACorr << endl;
                                                    }
                                                    if( flag_tof_particleB == 1 )
                                                    {
                                                        // pi-
                                                        BetaBCorr = correctBeta4DoubleDecay(trackB,trackAB_Lambda,trackLambdaC,helixB,vectorAB,PathAB,PathABC);
                                                        Mass2BCorr = MomentumB * MomentumB * (1./(BetaBCorr*BetaBCorr)-1.);
                                                        MassBCorr = SquareRoot(Mass2BCorr);
                                                        if ( debug_flag ) cout << "B) MassB = " << MassB << "\tMassBCorr = " <<  MassBCorr << endl;
                                                    }
                                                    if( flag_tof_particleC == 1 )
                                                    {
                                                        // pi- (C)
                                                        BetaCCorr = correctBeta4SingleDecay(trackC,trackLambdaC,helixC,vectorA2C,PathABC);
                                                        Mass2CCorr = MomentumC * MomentumC * (1./(BetaCCorr*BetaCCorr)-1.);
                                                        MassCCorr = SquareRoot(Mass2CCorr);
                                                        if ( debug_flag ) cout << "C) MassC = " << MassC << "\tMassCCorr = " <<  MassCCorr << endl;
                                                    }


                                                    if(flag_tof_particleA == 0) Mass2ACorr  = -100;
                                                    if(flag_tof_particleB == 0) Mass2BCorr  = -100;
                                                    if(flag_tof_particleC == 0) Mass2CCorr  = -100;
                                                    //****************************************************************************************



                                                    for(Int_t l = 0; l < PID_counter_Array_B[Ana_Num][ParticleD]; l++) // pi+ candidates
                                                    {
                                                        Int_t trackD_num = PID_Array_B[Ana_Num][ParticleD][l];

                                                        if(
                                                           (
                                                            trackD_num != trackA_num
                                                            && trackD_num != trackB_num
                                                            && trackD_num != trackC_num
                                                           )
                                                           || (SE_ME_Flag == 1)
                                                          )
                                                        {
                                                            StPicoAlexTrack trackD  = *event_B_ana->track( trackD_num );
                                                            Float_t MomentumD   = trackD.gMom().mag();
                                                            Float_t dcaD        = trackD.dca();   // distance of closest approach to primary vertex
                                                            Float_t nHitsPossD  = trackD.nHitsMax();
                                                            Float_t nHitsFitD   = trackD.nHitsFit();
                                                            //Float_t PolarityD   = trackD.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                                            Float_t BetaD       = trackD.btofBeta();  // Velocity after time-of-flight reconstruction
                                                            Float_t nSigmaPD    = trackD.nSigmaPion();
                                                            //Float_t TofD        = trackD.btof();
                                                            Float_t Mass2D      = -100.0;
                                                            // calculate mass2
                                                            if(trackD.btofMatchFlag() > 0 && trackD.btof() != 0 && BetaD != 0)
                                                            {
                                                                flag_tof_particleD = 1;
                                                                Mass2D = MomentumD*MomentumD*(1.0/(BetaD*BetaD) - 1.0);
                                                            }
                                                            else
                                                            {
                                                                flag_tof_particleD = 0;
                                                            }
                                                            Float_t MassD       = SquareRoot(Mass2D);

                                                            if(
                                                               dcaD          < dcaD_cut  // 0.15
                                                               && nHitsFitD  > 14
                                                               && nHitsPossD > 0.0
                                                               && MomentumD  > 0.1
                                                               && MomentumD  < 10.0
                                                               && (nHitsFitD/nHitsPossD) > 0.52
                                                               && flag_tof_particleD == 1
                                                               && nSigmaPD < 2.5
                                                               && nSigmaPD > -2.5
                                                               && (flag_tof_particleD == 0 || (flag_tof_particleD == 1 && Mass2D > -0.1 && Mass2D < 0.1)) // pion mass
                                                              )
                                                            {
                                                                helixD = StPhysicalHelixD(trackD.gMom(),trackD.origin()+vectordiff,event_B_ana->bField()*MAGFIELDFACTOR, trackD.charge());



                                                                for(Int_t m = 0; m < PID_counter_Array_B[Ana_Num][ParticleE]; m++) // K- candidates
                                                                {
                                                                    Int_t trackE_num = PID_Array_B[Ana_Num][ParticleE][m];

                                                                    if(
                                                                       (
                                                                        trackE_num != trackA_num
                                                                        && trackE_num != trackB_num
                                                                        && trackE_num != trackC_num
                                                                       )
                                                                       || (SE_ME_Flag == 1)
                                                                      )
                                                                    {
                                                                        StPicoAlexTrack trackE  = *event_B_ana->track( trackE_num );
                                                                        Float_t MomentumE   = trackE.gMom().mag();
                                                                        Float_t dcaE        = trackE.dca();   // distance of closest approach to primary vertex
                                                                        Float_t nHitsPossE  = trackE.nHitsMax();
                                                                        Float_t nHitsFitE   = trackE.nHitsFit();
                                                                        //Float_t PolarityE   = trackE.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                                                        Float_t BetaE       = trackE.btofBeta();  // Velocity after time-of-flight reconstruction
                                                                        Float_t nSigmaPE    = trackE.nSigmaKaon();
                                                                        //Float_t TofE        = trackE.btof();
                                                                        Float_t Mass2E      = -100.0;
                                                                        // calculate mass2
                                                                        if(trackE.btofMatchFlag() > 0 && trackE.btof() != 0 && BetaE != 0)
                                                                        {
                                                                            flag_tof_particleE = 1;
                                                                            Mass2E = MomentumE*MomentumE*(1.0/(BetaE*BetaE) - 1.0);
                                                                        }
                                                                        else
                                                                        {
                                                                            flag_tof_particleE = 0;
                                                                        }
                                                                        Float_t MassE       = SquareRoot(Mass2E);

                                                                        if(
                                                                           trackE_num != trackD_num // not a part of mixed event exclusion!
                                                                           && dcaE       < dcaE_cut
                                                                           && nHitsFitE  > 14
                                                                           && nHitsPossE > 0.0
                                                                           && MomentumE  > 0.1
                                                                           && MomentumE  < 10.0
                                                                           && (nHitsFitE/nHitsPossE) > 0.52
                                                                           && flag_tof_particleE == 1
                                                                           && nSigmaPE < 2.5
                                                                           && nSigmaPE > -2.5
                                                                           && (flag_tof_particleE == 0 || (flag_tof_particleE == 1 && Mass2E > 0.125 && Mass2E < 0.36)) // kaon mass
                                                                          )
                                                                        {
                                                                            helixE = StPhysicalHelixD(trackE.gMom(),trackE.origin()+vectordiff,event_B_ana->bField()*MAGFIELDFACTOR, trackE.charge());

                                                                            Float_t pathD_f, pathABC_f, dcaD_ABC_f, path_ABC_test, dca_omega_D_f;
                                                                            Float_t path_omega_test = 0.0;
                                                                            Float_t pathD_test      = 0.0;
                                                                            Int_t fDCA_Helix_out_D_ABC = fDCA_Helix_Estimate(helix_omega,helixD,path_ABC_test,pathD_test,dca_omega_D_f);

                                                                            if(fDCA_Helix_out_D_ABC == 1)
                                                                            {
                                                                                fHelixABdca_start_params(helix_omega,helixD,pathABC_f,pathD_f,dcaD_ABC_f,path_omega_test,pathD_test); // calculate dca between two helices
                                                                            }
                                                                            else
                                                                            {
                                                                                fHelixABdca(helix_omega,helixD,pathABC_f,pathD_f,dcaD_ABC_f); // calculate dca between two helices
                                                                            }

                                                                            StThreeVectorF vectoratsD       = helixD.at(pathD_f);  // space vector of helixA at dca to helixB
                                                                            StThreeVectorF vectoratsD_ABC   = helix_omega.at(pathABC_f);  // space vector of helixB at dca to helixA
                                                                            StThreeVectorF vectorD_ABC      = vectoratsD+vectoratsD_ABC;
                                                                            vectorD_ABC      = vectorD_ABC/2.0; // decay vertex



                                                                            Float_t pathE_f, dcaE_ABC_f, dca_omega_E_f;
                                                                            Float_t pathE_test      = 0.0;
                                                                            Int_t fDCA_Helix_out_E_ABC = fDCA_Helix_Estimate(helix_omega,helixE,path_ABC_test,pathE_test,dca_omega_E_f);

                                                                            if(fDCA_Helix_out_E_ABC == 1)
                                                                            {
                                                                                fHelixABdca_start_params(helix_omega,helixE,pathABC_f,pathE_f,dcaE_ABC_f,path_omega_test,pathE_test); // calculate dca between two helices
                                                                            }
                                                                            else
                                                                            {
                                                                                fHelixABdca(helix_omega,helixE,pathABC_f,pathE_f,dcaE_ABC_f); // calculate dca between two helices
                                                                            }

                                                                            StThreeVectorF vectoratsE       = helixE.at(pathE_f);  // space vector of helixA at dca to helixB
                                                                            StThreeVectorF vectoratsE_ABC   = helix_omega.at(pathABC_f);  // space vector of helixB at dca to helixA
                                                                            StThreeVectorF vectorE_ABC      = vectoratsE+vectoratsE_ABC;
                                                                            vectorE_ABC      = vectorE_ABC/2.0; // decay vertex

                                                                            StThreeVectorF vectorDE_ABC     = vectorD_ABC+vectorE_ABC;
                                                                            vectorDE_ABC = vectorDE_ABC/2.0;
                                                                            StThreeVectorF vectorDE_Diff    = vectorD_ABC-vectorE_ABC;
                                                                            Float_t mag_vectorDE_Diff = vectorDE_Diff.mag();

                                                                            StThreeVectorF vectorDE_ABC_toPrim     = vectorDE_ABC - vectorprim; // vector primary vertex to decay vertex
                                                                            Float_t VerdistX_DE_ABC = vectorDE_ABC_toPrim.mag(); // distance between primary vertex and decay vertex

                                                                            //cout << "Before last topology cut, VerdistX_DE_ABC = " << VerdistX_DE_ABC
                                                                            //    << ", mag_vectorDE_Diff = " << mag_vectorDE_Diff << ", dcaD = " << dcaD << ", dcaE = " << dcaE << endl;


                                                                            if(
                                                                               VerdistX_DE_ABC < 2.0
                                                                               && mag_vectorDE_Diff < 1.0
                                                                              )
                                                                            {
                                                                                vectornewD       = helixD.cat(pathD_f); // direction vector at dca for helixC
                                                                                vectornewD       = MomentumD*vectornewD/vectornewD.mag(); // new momentum vector at decay vertex

                                                                                vectornewE       = helixE.cat(pathE_f); // direction vector at dca for helixC
                                                                                vectornewE       = MomentumE*vectornewE/vectornewE.mag(); // new momentum vector at decay vertex

                                                                                ltrackD.SetXYZM(vectornewD.x(),vectornewD.y(),vectornewD.z(),0.13957018);
                                                                                ltrackE.SetXYZM(vectornewE.x(),vectornewE.y(),vectornewE.z(),0.493677);

                                                                                TLorentzVector trackABC_Xi_nominal; // mother particle with nominal Xi mass
                                                                                trackABC_Xi_nominal.SetXYZM(trackLambdaC.Px(),trackLambdaC.Py(),trackLambdaC.Pz(),1.32131);
                                                                                TLorentzVector ltrackABCDE    = trackABC_Xi_nominal+ltrackD; // mother particle
                                                                                ltrackABCDE += ltrackE;

                                                                                Double_t InvMassABCDE      = ltrackABCDE.M(); // invariant mass of mother particle
                                                                                Double_t PMassABCDE        = ltrackABCDE.P();

                                                                                //cout << "InvMassABCDE = " << InvMassABCDE << endl;

                                                                                if(erefMult_bin >= 0 && erefMult_bin < 9)
                                                                                {
                                                                                    h2D_Omega2250_m_vs_p_cent[erefMult_bin] ->Fill(PMassABCDE,InvMassABCDE);
                                                                                    h1D_Omega2250_m_cent[erefMult_bin]      ->Fill(InvMassABCDE);

                                                                                    //if(erefMult_bin == 0)
                                                                                    //{
                                                                                    //    dummy_counter_array[erefMult_bin]++;
                                                                                    //    cout << "------ Entry: dummy_counter = " << dummy_counter_array[erefMult_bin] << ", hist entries = "
                                                                                    //        << h1D_Omega2250_m_cent[erefMult_bin]->GetEntries() << endl;
                                                                                    //}

                                                                                }


                                                                                /*
                                                                                 alexV0_track_A = alexV0_event_A.createTrack();
                                                                                 alexV0_track_A->setm2A(Mass2ACorr);
                                                                                 alexV0_track_A->setm2B(Mass2BCorr);
                                                                                 alexV0_track_A->setm2C(Mass2CCorr);
                                                                                 alexV0_track_A->setnsA(nSigmaPA);
                                                                                 alexV0_track_A->setnsB(nSigmaPiB);
                                                                                 alexV0_track_A->setnsC(nSigmaC);
                                                                                 alexV0_track_A->setdcaA(dcaA);
                                                                                 alexV0_track_A->setdcaB(dcaB);
                                                                                 alexV0_track_A->setdcaC(dcaC);
                                                                                 alexV0_track_A->setiQxA(iQxA);
                                                                                 alexV0_track_A->setiQyA(iQyA);
                                                                                 alexV0_track_A->setiQxB(iQxB);
                                                                                 alexV0_track_A->setiQyB(iQyB);
                                                                                 alexV0_track_A->setiQxC(iQxC);
                                                                                 alexV0_track_A->setiQyC(iQyC);
                                                                                 alexV0_track_A->setetaA(etaA_c);
                                                                                 alexV0_track_A->setetaB(etaB_c);
                                                                                 alexV0_track_A->setetaC(etaC_c);
                                                                                 alexV0_track_A->setInvAB(InvMassAB);
                                                                                 alexV0_track_A->setInvABC(InvMassABC);
                                                                                 alexV0_track_A->setInvAB_miss(InvMassAB_K0S);
                                                                                 alexV0_track_A->setInvABC_miss(InvMassABC_Xi);
                                                                                 alexV0_track_A->setdcaAB(dcaAB_f);
                                                                                 alexV0_track_A->setdcaBC(dcaBC_f);
                                                                                 alexV0_track_A->setdcaABC(dip_omega);
                                                                                 alexV0_track_A->setVerdistX(VerdistX);
                                                                                 alexV0_track_A->setVerdistY(VerdistY);
                                                                                 alexV0_track_A->setVerdistX2(VerdistX2);
                                                                                 alexV0_track_A->setVerdistY2(dca_omega);
                                                                                 alexV0_track_A->setpt(pt2);
                                                                                 alexV0_track_A->setrap(rap2);
                                                                                 alexV0_track_A->setphi(phiABC);
                                                                                 alexV0_track_A->settheta(thetaABC);
                                                                                 alexV0_track_A->setPsi_ep(phi_event_plane);
                                                                                 alexV0_track_A->setPsi_ep_eta(phi_event_plane_eta_gap);
                                                                                 alexV0_track_A->setPsi_diff_ME(delta_phi_ME_AB_weight);
                                                                                 alexV0_track_A->setscal_prod(scalarProduct);
                                                                                 alexV0_track_A->setscal_prod2(scalarProduct2);
                                                                                 */

                                                                                dummy_counter_loop++;
                                                                            }

                                                                        }
                                                                    }
                                                                } // End of m-loop (K-)

                                                            }
                                                        }
                                                    } // End of l-loop (pi+)

                                                }
                                            }
                                        }
                                    } // End of k-loop (2nd pi-)
                                }
                            }
                        }
                    }
                } // End of j-loop (1st pi-)

            }
        } // End of i-loop (p)

        if(ParticleA == 14 && ParticleB == 9)
        {
            //Tree_OmegaMV0_v2  ->Fill();
            //Tree_XiMV0_v2     ->Fill();
        }
        if(ParticleA == 15 && ParticleB == 8)
        {
            //Tree_OmegaPV0_v2  ->Fill();
            //Tree_XiPV0_v2     ->Fill();
        }
        return 1;
    }
    else return 0;
}



Int_t d4s_analysis(Int_t PID_counter_Array_A[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs],
                  Int_t PID_Array_A[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                  StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                  Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t ParticleD, Int_t ParticleE, Int_t Ana_Num,Int_t SE_ME_Flag, Int_t rot_background_flag)
{
    // d4s(????)
    // mass = ???? MeV/c2
    // width = ?? +/- ?? MeV
    //
    // Au + Au -> d4s(????)- + N + N
    //             |
    //             -> Xi-                  +     Lambda         +        pi0
    //                |                            |                      |
    //                -> Lambda + pi-(C)           -> p(D) + pi-(E)       -> not measurable
    //                  |
    //                  -> p(A) + pi-(B)
    //
    // (A,B,C) -> eventA
    // (D,E)   -> eventB

    dummy_counter_A++;
    // Event vertex information
    StThreeVectorF vectorprim,vectorprimB,vectordiff, vectornewA_lin, vectornewB_lin;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB,RunIdA,RunIdB;


    event_A_ana       = picoDst_A;
    event_B_ana       = picoDst_B;

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();

    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    refMultB          = refMultA;
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_A_ana->vzVpd();

    vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    Int_t flag_tof_particleA = 0;
    Int_t flag_tof_particleB = 0;
    Int_t flag_tof_particleC = 0;
    Int_t flag_tof_particleD = 0;
    Int_t flag_tof_particleE = 0;


    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        RunIdB         = event_B_ana->runId();
        vectorprimB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vectorprim - vectorprimB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }


    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
        //cout << "Triggers checked, Vertex checked..." << endl;
    }

    //
    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array_A[Ana_Num][ParticleA]    > 0   // p(A)
       && PID_counter_Array_A[Ana_Num][ParticleB]    > 0   // pi-(B)
       && PID_counter_Array_A[Ana_Num][ParticleC]    > 0   // pi-(C)
       && PID_counter_Array_B[Ana_Num][ParticleD]    > 0   // p(D)
       && PID_counter_Array_B[Ana_Num][ParticleE]    > 0   // pi-(E)
       //&& event->getNVertices() >= 1
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana   ->isMinBias()
       && event_B_ana   ->isMinBias()
      )
    {
        //cout << "************** START of event with " << PID_counter_Array[Ana_Num][ParticleA] << " protons and "
        //   << PID_counter_Array_B[Ana_Num][ParticleB] << " pions, refMultA = " << refMultA << ", refMultB = " << refMultB << " ***********" << endl;
        // Loop over all particle combinations

        Float_t dcaA_cut            = 0.1;
        Float_t dcaB_cut            = 0.9;
        Float_t dcaC_cut            = 0.5;
        Float_t dcaD_cut            = 0.1; // 1.0
        Float_t dcaE_cut            = 0.9;
        Float_t VerdistX_cut        = 3.5;
        Float_t VerdistX_DE_cut     = 3.5;
        Float_t dcaAB_cut           = 2.0;
        Float_t dcaDE_cut           = 2.0;
        Float_t InvMassABC_cut      = 1.3273;
        Float_t InvMassABC_low_cut  = 1.3157;
        Float_t dca_mother_cut      = 3.0;
        Float_t dcaBC_cut           = 1.0;
        Float_t VerdistX2_cut       = 2.5;


        StPhysicalHelixD helixA, helixB, helixC, helixD, helixE;
        StThreeVectorF vectorA, vectorB, vectoratsA, vectoratsB, vectorAB, vectorDE, vectorAB_est, vectorprimAB, vectornewA, vectornewB, vectornewD, vectornewE, dirY_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY, dirY_DE;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackD, ltrackE, ltrackC_pi, ltrackB2, ltrackA_lin, ltrackB_lin, ltrackB2_pi, ltrackA_pip;

        /*
        // Fill event information for Xi-/+
        alexV0_event_A.clearTrackList();
        alexV0_event_A.setx(EventVertexXA);
        alexV0_event_A.sety(EventVertexYA);
        alexV0_event_A.setz(EventVertexZA);
        alexV0_event_A.setid(RunIdA);
        alexV0_event_A.setmult(refMultA);
        alexV0_event_A.setn_prim(n_primaries);
        alexV0_event_A.setn_non_prim(n_non_primaries);
        alexV0_event_A.setn_tof_prim(n_tofmatch_prim);
        alexV0_event_A.setEP_Qx_eta_pos_ptw(EP_Qx_eta_pos_ptw);
        alexV0_event_A.setEP_Qy_eta_pos_ptw(EP_Qy_eta_pos_ptw);
        alexV0_event_A.setEP_Qx_eta_neg_ptw(EP_Qx_eta_neg_ptw);
        alexV0_event_A.setEP_Qy_eta_neg_ptw(EP_Qy_eta_neg_ptw);
        alexV0_event_A.setEP_Qx_ptw(EP_Qx_ptw);
        alexV0_event_A.setEP_Qy_ptw(EP_Qy_ptw);
        alexV0_event_A.setQtracks_eta_pos(Qtracks_used_eta_pos);
        alexV0_event_A.setQtracks_eta_neg(Qtracks_used_eta_neg);
        alexV0_event_A.setQtracks_full(Qtracks_used);
        alexV0_event_A.setZDCx(ZDCx);
        alexV0_event_A.setBBCx(BBCx);
        alexV0_event_A.setvzVpd(vzVpd);
        */

        if(0)
        {
            cout << "d4s event accepted, A = " << PID_counter_Array_A[Ana_Num][ParticleA] <<
                ", B = " << PID_counter_Array_A[Ana_Num][ParticleB] <<
                ", C = " << PID_counter_Array_A[Ana_Num][ParticleC] <<
                ", D = " << PID_counter_Array_B[Ana_Num][ParticleD] <<
                ", E = " << PID_counter_Array_B[Ana_Num][ParticleE] <<
                ", erefMult_bin = " << erefMult_bin <<
                endl;
        }

        for(Int_t i = 0; i < PID_counter_Array_A[Ana_Num][ParticleA]; i++) // p candidates
        {
            Int_t trackA_num = PID_Array_A[Ana_Num][ParticleA][i];

            StPicoAlexTrack trackA  = *event_A_ana->track( trackA_num );
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            //Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t nSigmaPA    = trackA.nSigmaProton();
            //Float_t TofA        = trackA.btof();
            Float_t Mass2A      = -100.0;
            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                flag_tof_particleA = 1;
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }
            else
            {
                flag_tof_particleA = 0;
            }
            Float_t MassA       = SquareRoot(Mass2A);

            if(
               dcaA          > dcaA_cut  // 0.15
               && nHitsFitA  > 14
               && nHitsPossA > 0.0
               && MomentumA  > 0.1
               && MomentumA  < 10.0
               && (nHitsFitA/nHitsPossA) > 0.52
               && (flag_tof_particleA == 0 || (flag_tof_particleA == 1 && Mass2A > 0.4 && Mass2A < 1.5)) // proton mass
              )
            {
                helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());

                for(Int_t j = 0; j < PID_counter_Array_A[Ana_Num][ParticleB]; j++)  // pi- candidates
                {
                    Int_t trackB_num = PID_Array_A[Ana_Num][ParticleB][j];
                    StPicoAlexTrack trackB  = *event_A_ana->track( trackB_num ); // take again event A
                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    //Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                    Float_t nSigmaPiB   = trackB.nSigmaPion();
                    //Float_t TofB        = trackB.btof();
                    Float_t Mass2B      = -100.0;
                    // calculate mass2
                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        flag_tof_particleB = 1;
                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                    }
                    else
                    {
                        flag_tof_particleB = 0;
                    }
                    Float_t MassB       = SquareRoot(Mass2B);

                    //if( flag_tof_particleA == 0 && flag_tof_particleB == 0) // apply more strict cuts if there is no TOF information for the two particles
                    //{
                    //    dcaA_cut  = 0.1;
                    //    dcaB_cut  = 1.5;
                    //    dcaC_cut  = 0.7;
                    //}

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       && dcaB       > dcaB_cut // 1.3
                       && nHitsFitB  > 14
                       && nHitsPossB > 0.0
                       && MomentumB  > 0.1
                       && MomentumB  < 10.0
                       && (nHitsFitB/nHitsPossB) > 0.52
                       && (flag_tof_particleB == 0 || (flag_tof_particleB == 1 && Mass2B > -0.1 && Mass2B < 0.1)) // pion mass
                      )
                    {
                        helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackB.charge());

                        TLorentzVector trackAB, ltrackA, ltrackB;
                        Float_t VerdistX, dcaAB;
                        Int_t decay_prop_AB = calc_Decay_Properties(helixA,helixB,mass_array[ParticleA],mass_array[ParticleB],MomentumA,MomentumB,
                                                                 vectorprim,VerdistX_cut,dcaAB_cut,vectorAB,ltrackA,ltrackB,trackAB,VerdistX,dcaAB);


                        //cout << "Lambda 1 reconstructed" << endl;

                        if( decay_prop_AB == 1 )
                        {
                            Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                            Float_t  MomentumAB         = trackAB.P(); // momentum of mother particle
                            Float_t  BetaAB = TMath::Sqrt(1./(1+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            dirY *= 1.0/dirY.mag();

                            baseY = vectorAB;
                            Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vectorprim);

                            if(
                               VerdistX > VerdistX_cut
                               && dcaAB < dcaAB_cut
                               && VerdistY > 0.2
                              )
                            {
                                h1D_Lambda->Fill(InvMassAB);
                                //cout << "InvMassAB = " << InvMassAB << endl;
                            }

                            // Apply Lambda cuts
                            if(
                               InvMassAB         > 1.11
                               && InvMassAB      < 1.12
                               && VerdistX       > VerdistX_cut // 4.0
                               && dcaAB          < 1.0  // 1.0
                               && VerdistY       > 0.2  // 0.2
                              )
                            {
                                //cout << "Lambda 1 accepted" << endl;
                                for(Int_t k = 0; k < PID_counter_Array_A[Ana_Num][ParticleC]; k++)  // pi-
                                {
                                    Int_t trackC_num = PID_Array_A[Ana_Num][ParticleC][k];
                                    if(
                                       trackC_num != trackB_num
                                       && trackC_num != trackA_num
                                      )
                                    {
                                        StThreeVectorF vectorC, vectoratsC, vectoratsA2, vectorA2C, vectorA2CtoPrim, vectornewC, dirY2, baseY2;

                                        StPicoAlexTrack trackC  = *event_A_ana->track( trackC_num ); // take now event B
                                        Float_t MomentumC   = trackC.gMom().mag();
                                        Float_t dcaC        = trackC.dca();   // distance of closest approach to primary vertex
                                        Float_t nHitsPossC  = trackC.nHitsMax();
                                        Float_t nHitsFitC   = trackC.nHitsFit();
                                        //Float_t PolarityC   = trackC.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                        Float_t BetaC       = trackC.btofBeta();  // Velocity after time-of-flight reconstruction
                                        Float_t nSigmaC     = -100.0;
                                        if(ParticleC == 8 || ParticleC == 9)   nSigmaC = trackC.nSigmaPion(); // for Xi analysis
                                        if(ParticleC == 11 || ParticleC == 12) nSigmaC = trackC.nSigmaKaon(); // for Omega analysis
                                        //Float_t TofC        = trackC.btof();
                                        //Float_t TPCdEdxC    = trackC.dEdx(); // TPC dE/dx
                                        Float_t Mass2C      = -100.0;
                                        // calculate mass2
                                        if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                                        {
                                            flag_tof_particleC = 1;
                                            Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                                        }
                                        else
                                        {
                                            flag_tof_particleC = 0;
                                        }
                                        Float_t MassC       = SquareRoot(Mass2C);

                                        if(
                                           dcaC          > dcaC_cut  // 0.4
                                           && nHitsFitC  > 14
                                           && nHitsPossC > 0
                                           && MomentumC  > 0.1
                                           && MomentumC  < 10.0
                                           && (nHitsFitC/nHitsPossC) > 0.52
                                           && (
                                               flag_tof_particleC == 0 ||
                                               ((ParticleC == 8  || ParticleC == 9)  && (flag_tof_particleC == 1 && Mass2C > -0.1 && Mass2C < 0.1)) ||
                                               ((ParticleC == 11 || ParticleC == 12) && (flag_tof_particleC == 1 && Mass2C > 0.125 && Mass2C < 0.36))
                                              )
                                          )
                                        {
                                            helixC = StPhysicalHelixD(trackC.gMom(),trackC.origin(),event_A_ana->bField()*MAGFIELDFACTOR,trackC.charge());

                                            Float_t pathA2_f, pathC_f,dcaBC_f;
                                            fHelixAtoLinedca(dirY,baseY,helixC,pathC_f,pathA2_f,dcaBC_f); // calculates dca of second pi- helix to Lambda track

                                            vectoratsA2       = baseY+pathA2_f*dirY;  // space vector of Lambda at dca to helixC
                                            vectoratsC        = helixC.at(pathC_f);  // space vector of helixC at dca to Lambda
                                            vectorA2C         = vectoratsA2+vectoratsC;
                                            vectorA2C         = vectorA2C/2.0; // decay vertex
                                            vectorA2CtoPrim   = vectorA2C - vectorprim; // vector primary vertex to decay vertex
                                            Float_t VerdistX2 = vectorA2CtoPrim.mag(); // distance between primary vertex and decay vertex

                                            vectornewC       = helixC.cat(pathC_f); // direction vector at dca for helixC
                                            vectornewC       = MomentumC*vectornewC/vectornewC.mag(); // new momentum vector at decay vertex

                                            if(ParticleC == 8 || ParticleC == 9) // Xi analysis
                                            {
                                                ltrackC.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.13957018);
                                                ltrackC_pi.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.493677);
                                            }
                                            if(ParticleC == 11 || ParticleC == 12) // Omega analysis
                                            {
                                                ltrackC.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.493677);
                                                ltrackC_pi.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),0.13957018);
                                            }

                                            // Missing mass and invariant mass calculations
                                            TLorentzVector trackAB_Lambda; // mother particle with nominal Lambda mass
                                            trackAB_Lambda.SetXYZM(trackAB.Px(),trackAB.Py(),trackAB.Pz(),1.115683);
                                            TLorentzVector trackLambdaC    = trackAB_Lambda+ltrackC; // mother particle
                                            TLorentzVector trackLambdaC_pi = trackAB_Lambda+ltrackC_pi; // mother particle

                                            Double_t InvMassABC      = trackLambdaC.M(); // invariant mass of mother particle
                                            Double_t InvMassABC_Xi   = trackLambdaC_pi.M(); // invariant mass of mother particle
                                            Double_t MomentumABC     = trackLambdaC.P(); // momentum of mother particle
                                            Float_t BetaABC = TMath::Sqrt(1./(1+(InvMassABC/MomentumABC)*(InvMassABC/MomentumABC)));

                                            dirY2.set(trackLambdaC.Px(),trackLambdaC.Py(),trackLambdaC.Pz()); // direction vector of Xi
                                            dirY2 = dirY2/dirY2.mag();
                                            Double_t scalarProduct2 = dirY2.dot(vectorA2CtoPrim/vectorA2CtoPrim.mag());

                                            baseY2 = vectorA2C;
                                            //Double_t  VerdistY2  = calculateMinimumDistanceStraightToPoint(baseY2,dirY2,vectorprim);

                                            Float_t pt2          = trackLambdaC.Pt();  // Transverse momentum of mother particle
                                            Float_t rap2         = trackLambdaC.Rapidity(); // Rapidity of mother particle

                                            Float_t phiABC   = dirY2.phi();
                                            Float_t thetaABC = dirY2.theta();


                                            //********************** Calculate the helix of the Omega/Xi ***************************
                                            StPhysicalHelixD helix_omega;
                                            StThreeVectorF origin_omega;   // decay vertex of Omega
                                            origin_omega.set(vectorA2C.x(),vectorA2C.y(),vectorA2C.z());
                                            //Float_t h_omega     = trackC.geth(); // +1 for K+ and -1 for K-
                                            Float_t h_omega;
                                            if(ParticleC == 12 || ParticleC == 9) h_omega = -1.0;
                                            if(ParticleC == 11 || ParticleC == 8) h_omega = 1.0;
                                            Float_t phase_omega = TMath::ATan2(-1.0*h_omega*trackLambdaC.Px(),h_omega*trackLambdaC.Py());
                                            Float_t dip_omega   = TMath::ATan2(trackLambdaC.Pz(),pt2);  // correct
                                            Float_t curv_omega  = 1.0;
                                            if(pt2 != 0.0)
                                            {
                                                curv_omega = curv_to_invpt_ratio/pt2;
                                            }

                                            helix_omega.setParameters(curv_omega,dip_omega,phase_omega,origin_omega,h_omega);
                                            Float_t path_omega = -999.0;
                                            Float_t dca_omega  = -999.0;
                                            fHelixAtoPointdca(vectorprim,helix_omega,path_omega,dca_omega);
                                            //**************************************************************************************

                                            if(
                                               dca_omega         < dca_mother_cut  // 0.5
                                               && dcaBC_f        < dcaBC_cut  // 0.8
                                               && VerdistX2      > VerdistX2_cut  // 2.5
                                               && scalarProduct2 > 0.0
                                              )
                                            {
                                                h1D_Xi->Fill(InvMassABC);
                                            }

                                            if(
                                               InvMassABC        < InvMassABC_cut // 2.5
                                               && InvMassABC     > InvMassABC_low_cut // 2.5
                                               && dca_omega      < dca_mother_cut  // 0.5
                                               && dcaBC_f        < dcaBC_cut  // 0.8
                                               && VerdistX2      > VerdistX2_cut  // 2.5
                                               && scalarProduct2 > 0.0
                                              )
                                            {

                                                cout << "Xi accepted" << endl;

                                                //****************************** BETA CORRECTION OMEGA ******************************
                                                Float_t BetaACorr   = BetaA;
                                                Float_t Mass2ACorr  = Mass2A;
                                                Float_t BetaBCorr   = BetaB;
                                                Float_t Mass2BCorr  = Mass2B;
                                                Float_t BetaCCorr   = BetaC;
                                                Float_t Mass2CCorr  = Mass2C;
                                                Float_t MassACorr   = MassA;
                                                Float_t MassBCorr   = MassB;
                                                Float_t MassCCorr   = MassC;
                                                Float_t PathAB      = (vectorAB-vectorA2C).mag(); // lambda uncharged
                                                Float_t PathABC     = vectorA2CtoPrim.mag(); // Xi charged, but curvature negligible

                                                if( flag_tof_particleA == 1 || flag_tof_particleB == 1 || flag_tof_particleC == 1 )
                                                {
                                                    if( debug_flag )
                                                    {
                                                        cout << "----------------------------------------------- OMEGA -----------------------------------------------------------------------" << endl;
                                                        cout << "PathAB = " << PathAB << "\tInvMassAB = " << InvMassAB << "\tMomentumAB = " << MomentumAB
                                                            << "\tBetaAB = " << BetaAB << "\tTofAB = " << PathAB / (BetaAB*clight) << endl;
                                                        cout << "PathABC = " << PathABC << "\tInvMassABC = " << InvMassABC << "\tMomentumABC = " << MomentumABC
                                                            << "\tBetaABC = " << BetaABC << "\tTofABC = " << PathABC / (BetaABC*clight) << endl;
                                                    }
                                                }
                                                if( flag_tof_particleA == 1 )
                                                {
                                                    // proton
                                                    BetaACorr = correctBeta4DoubleDecay(trackA,trackAB_Lambda,trackLambdaC,helixA,vectorAB,PathAB,PathABC);
                                                    Mass2ACorr = MomentumA * MomentumA * (1./(BetaACorr*BetaACorr)-1.);
                                                    MassACorr = SquareRoot(Mass2ACorr);
                                                    if ( debug_flag ) cout << "A) MassA = " << MassA << "\tMassACorr = " <<  MassACorr << endl;
                                                }
                                                if( flag_tof_particleB == 1 )
                                                {
                                                    // pi-
                                                    BetaBCorr = correctBeta4DoubleDecay(trackB,trackAB_Lambda,trackLambdaC,helixB,vectorAB,PathAB,PathABC);
                                                    Mass2BCorr = MomentumB * MomentumB * (1./(BetaBCorr*BetaBCorr)-1.);
                                                    MassBCorr = SquareRoot(Mass2BCorr);
                                                    if ( debug_flag ) cout << "B) MassB = " << MassB << "\tMassBCorr = " <<  MassBCorr << endl;
                                                }
                                                if( flag_tof_particleC == 1 )
                                                {
                                                    // pi- (C)
                                                    BetaCCorr = correctBeta4SingleDecay(trackC,trackLambdaC,helixC,vectorA2C,PathABC);
                                                    Mass2CCorr = MomentumC * MomentumC * (1./(BetaCCorr*BetaCCorr)-1.);
                                                    MassCCorr = SquareRoot(Mass2CCorr);
                                                    if ( debug_flag ) cout << "C) MassC = " << MassC << "\tMassCCorr = " <<  MassCCorr << endl;
                                                }


                                                if(flag_tof_particleA == 0) Mass2ACorr  = -100;
                                                if(flag_tof_particleB == 0) Mass2BCorr  = -100;
                                                if(flag_tof_particleC == 0) Mass2CCorr  = -100;
                                                //****************************************************************************************



                                                for(Int_t l = 0; l < PID_counter_Array_B[Ana_Num][ParticleD]; l++) // p candidates
                                                {
                                                    Int_t trackD_num = PID_Array_B[Ana_Num][ParticleD][l];

                                                    if(
                                                       (
                                                        trackD_num != trackA_num
                                                        && trackD_num != trackB_num
                                                        && trackD_num != trackC_num
                                                       )
                                                       || (SE_ME_Flag == 1)
                                                      )
                                                    {
                                                        StPicoAlexTrack trackD  = *event_B_ana->track( trackD_num );
                                                        Float_t MomentumD   = trackD.gMom().mag();
                                                        Float_t dcaD        = trackD.dca();   // distance of closest approach to primary vertex
                                                        Float_t nHitsPossD  = trackD.nHitsMax();
                                                        Float_t nHitsFitD   = trackD.nHitsFit();
                                                        //Float_t PolarityD   = trackD.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                                        Float_t BetaD       = trackD.btofBeta();  // Velocity after time-of-flight reconstruction
                                                        Float_t nSigmaPD    = trackD.nSigmaProton();
                                                        //Float_t TofD        = trackD.btof();
                                                        Float_t Mass2D      = -100.0;
                                                        // calculate mass2
                                                        if(trackD.btofMatchFlag() > 0 && trackD.btof() != 0 && BetaD != 0)
                                                        {
                                                            flag_tof_particleD = 1;
                                                            Mass2D = MomentumD*MomentumD*(1.0/(BetaD*BetaD) - 1.0);
                                                        }
                                                        else
                                                        {
                                                            flag_tof_particleD = 0;
                                                        }
                                                        Float_t MassD       = SquareRoot(Mass2D);

                                                        if(
                                                           dcaD          > dcaD_cut  // 0.15
                                                           && nHitsFitD  > 14
                                                           && nHitsPossD > 0.0
                                                           && MomentumD  > 0.1
                                                           && MomentumD  < 10.0
                                                           && (nHitsFitD/nHitsPossD) > 0.52
                                                           //&& flag_tof_particleD == 1
                                                           && nSigmaPD < 3.0
                                                           && nSigmaPD > -3.0
                                                           && (flag_tof_particleD == 0 || (flag_tof_particleD == 1 && Mass2D > 0.4 && Mass2D < 1.5)) // proton mass
                                                          )
                                                        {
                                                            helixD = StPhysicalHelixD(trackD.gMom(),trackD.origin()+vectordiff,event_B_ana->bField()*MAGFIELDFACTOR, trackD.charge());



                                                            for(Int_t m = 0; m < PID_counter_Array_B[Ana_Num][ParticleE]; m++) // pi- candidates
                                                            {
                                                                Int_t trackE_num = PID_Array_B[Ana_Num][ParticleE][m];

                                                                if(
                                                                   (
                                                                    trackE_num != trackA_num
                                                                    && trackE_num != trackB_num
                                                                    && trackE_num != trackC_num
                                                                   )
                                                                   || (SE_ME_Flag == 1)
                                                                  )
                                                                {
                                                                    StPicoAlexTrack trackE  = *event_B_ana->track( trackE_num );
                                                                    Float_t MomentumE   = trackE.gMom().mag();
                                                                    Float_t dcaE        = trackE.dca();   // distance of closest approach to primary vertex
                                                                    Float_t nHitsPossE  = trackE.nHitsMax();
                                                                    Float_t nHitsFitE   = trackE.nHitsFit();
                                                                    //Float_t PolarityE   = trackE.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                                                                    Float_t BetaE       = trackE.btofBeta();  // Velocity after time-of-flight reconstruction
                                                                    Float_t nSigmaPE    = trackE.nSigmaPion();
                                                                    //Float_t TofE        = trackE.btof();
                                                                    Float_t Mass2E      = -100.0;
                                                                    // calculate mass2
                                                                    if(trackE.btofMatchFlag() > 0 && trackE.btof() != 0 && BetaE != 0)
                                                                    {
                                                                        flag_tof_particleE = 1;
                                                                        Mass2E = MomentumE*MomentumE*(1.0/(BetaE*BetaE) - 1.0);
                                                                    }
                                                                    else
                                                                    {
                                                                        flag_tof_particleE = 0;
                                                                    }
                                                                    Float_t MassE       = SquareRoot(Mass2E);

                                                                    if(
                                                                       trackE_num != trackD_num // not a part of mixed event exclusion!
                                                                       && dcaE       > dcaE_cut
                                                                       && nHitsFitE  > 14
                                                                       && nHitsPossE > 0.0
                                                                       && MomentumE  > 0.1
                                                                       && MomentumE  < 10.0
                                                                       && (nHitsFitE/nHitsPossE) > 0.52
                                                                       //&& flag_tof_particleE == 1
                                                                       && nSigmaPE < 3.0
                                                                       && nSigmaPE > -3.0
                                                                       && (flag_tof_particleE == 0 || (flag_tof_particleE == 1 && Mass2E > -0.1 && Mass2E < 0.1)) // pion mass
                                                                      )
                                                                    {
                                                                        helixE = StPhysicalHelixD(trackE.gMom(),trackE.origin()+vectordiff,event_B_ana->bField()*MAGFIELDFACTOR, trackE.charge());


                                                                        TLorentzVector trackDE, ltrackD, ltrackE;
                                                                        Float_t VerdistX_DE, dcaDE;
                                                                        Int_t decay_prop_DE = calc_Decay_Properties(helixD,helixE,mass_array[ParticleD],mass_array[ParticleE],MomentumD,MomentumE,
                                                                                                                    vectorprim,VerdistX_DE_cut,dcaDE_cut,vectorDE,ltrackD,ltrackE,trackDE,VerdistX_DE,dcaDE);

                                                                        //cout << "Lambda 2 reconstructed" << endl;

                                                                        if(rot_background_flag == 1) trackDE.SetPxPyPzE(-1.0*trackDE.Px(),-1.0*trackDE.Py(),-1.0*trackDE.Pz(),trackDE.E());

                                                                        if(
                                                                           decay_prop_DE == 1
                                                                           && VerdistX_DE > VerdistX_DE_cut
                                                                           && dcaDE < dcaDE_cut
                                                                          )
                                                                        {
                                                                            Double_t InvMassDE          = trackDE.M(); // invariant mass of mother particle
                                                                            Float_t  MomentumDE         = trackDE.P(); // momentum of mother particle
                                                                            Float_t  BetaDE = TMath::Sqrt(1./(1+(InvMassDE/MomentumDE)*(InvMassDE/MomentumDE)));

                                                                            StThreeVectorF dir_Lambda_DE;
                                                                            dir_Lambda_DE.set(trackDE.Px(),trackDE.Py(),trackDE.Pz());
                                                                            dir_Lambda_DE *= 1.0/dir_Lambda_DE.mag();

                                                                            Double_t  VerdistY_DE  = calculateMinimumDistanceStraightToPoint(vectorDE,dir_Lambda_DE,vectorprim);

                                                                            if(
                                                                               VerdistY_DE > 0.2
                                                                              )
                                                                            {
                                                                                h1D_Lambda_B->Fill(InvMassDE);
                                                                            }

                                                                            if(VerdistY_DE  > 0.2
                                                                               && InvMassDE > 1.11
                                                                               && InvMassDE < 1.12
                                                                              )
                                                                            {
                                                                                Float_t path_helix_Xi_to_DE, path_line_Xi_to_DE, dca_Xi_to_DE;
                                                                                fHelixAtoLinedca(dir_Lambda_DE,vectorDE,helix_omega,path_helix_Xi_to_DE,path_line_Xi_to_DE,dca_Xi_to_DE); // calculates dca of second pi- helix to Lambda track

                                                                                StThreeVectorF vector_line_Xi_to_DE, vector_helix_Xi_to_DE, vector_Xi_Lambda;
                                                                                vector_line_Xi_to_DE  = vectorDE+path_line_Xi_to_DE*dir_Lambda_DE;  // space vector of Lambda at dca to helixC
                                                                                vector_helix_Xi_to_DE = helix_omega.at(path_helix_Xi_to_DE);  // space vector of helixC at dca to Lambda
                                                                                vector_Xi_Lambda      = vector_line_Xi_to_DE + vector_helix_Xi_to_DE;
                                                                                vector_Xi_Lambda *= 0.5; // Decay vector of Xi and second Lambda

                                                                                TLorentzVector trackABCDE;
                                                                                trackABCDE = trackLambdaC + trackDE;

                                                                                Float_t VerdistX0 = (vector_Xi_Lambda - vectorprim).mag();

                                                                                Float_t PtMassABCDE   = trackABCDE.Pt();
                                                                                Float_t PMassABCDE    = trackABCDE.P();
                                                                                Float_t InvMassABCDE  = trackABCDE.M();

                                                                                StThreeVectorF dir_ABCDE;
                                                                                dir_ABCDE.set(trackABCDE.Px(),trackABCDE.Py(),trackABCDE.Pz());
                                                                                Double_t  VerdistY_ABCDE  = calculateMinimumDistanceStraightToPoint(vector_Xi_Lambda,dir_ABCDE,vectorprim);

                                                                                //********************** Calculate the helix of the d4s ***************************
                                                                                StPhysicalHelixD helix_d4s;
                                                                                StThreeVectorF origin_d4s;   // decay vertex of d4s
                                                                                origin_d4s.set(vector_Xi_Lambda.x(),vector_Xi_Lambda.y(),vector_Xi_Lambda.z());
                                                                                Float_t h_d4s = -1.0;
                                                                                Float_t phase_d4s = TMath::ATan2(-1.0*h_d4s*trackABCDE.Px(),h_d4s*trackABCDE.Py());
                                                                                Float_t dip_d4s   = TMath::ATan2(trackABCDE.Pz(),PtMassABCDE);  // correct
                                                                                Float_t curv_d4s  = 1.0;
                                                                                if(PtMassABCDE != 0.0)
                                                                                {
                                                                                    curv_d4s = curv_to_invpt_ratio/PtMassABCDE;
                                                                                }

                                                                                helix_d4s.setParameters(curv_d4s,dip_d4s,phase_d4s,origin_d4s,h_d4s);
                                                                                Float_t path_d4s = -999.0;
                                                                                Float_t dca_d4s  = -999.0;
                                                                                fHelixAtoPointdca(vectorprim,helix_d4s,path_d4s,dca_d4s);

                                                                                if(
                                                                                   dcaA             > 0.3
                                                                                   && dcaB          > 1.2
                                                                                   && dcaC          > 0.5
                                                                                   && dcaD          > 1.0
                                                                                   && dcaE          > 1.2
                                                                                   && dcaDE         < 1.0
                                                                                   && dcaAB         < 1.0
                                                                                   && dca_Xi_to_DE  < 1.0
                                                                                   && dcaBC_f       < 1.0
                                                                                   && VerdistX      > 3.5
                                                                                   && VerdistX      > (VerdistX2 + 1.0)
                                                                                   && VerdistX_DE   > 4.5
                                                                                   && VerdistX2     > 1.5
                                                                                   && VerdistX0     > 4.0
                                                                                   && VerdistY      > 0.4
                                                                                   && VerdistY_DE   > 0.4
                                                                                   && dca_omega     > 0.4
                                                                                   && InvMassAB     < 1.11848
                                                                                   && InvMassAB     > 1.11231
                                                                                   && InvMassDE     < 1.11848
                                                                                   && InvMassDE     > 1.11231
                                                                                   && InvMassABC    < 1.32497
                                                                                   && InvMassABC    > 1.3183
                                                                                  )
                                                                                {
                                                                                    cout << "VerdistY_ABCDE = " << VerdistY_ABCDE << ", dca_d4s = " << dca_d4s << endl;
                                                                                }


                                                                                if(
                                                                                   VerdistX0 > 0.0
                                                                                  )
                                                                                {

                                                                                    if(erefMult_bin >= 0 && erefMult_bin < 9)
                                                                                    {
                                                                                        h2D_Omega2250_m_vs_p_cent[erefMult_bin] ->Fill(PMassABCDE,InvMassABCDE);
                                                                                        h1D_Omega2250_m_cent[erefMult_bin]      ->Fill(InvMassABCDE);


                                                                                        d4s_NTDataArray[0]     =(Float_t)InvMassABCDE;
                                                                                        d4s_NTDataArray[1]     =(Float_t)InvMassAB;
                                                                                        d4s_NTDataArray[2]     =(Float_t)InvMassABC;
                                                                                        d4s_NTDataArray[3]     =(Float_t)InvMassDE;
                                                                                        d4s_NTDataArray[4]     =(Float_t)VerdistX;
                                                                                        d4s_NTDataArray[5]     =(Float_t)VerdistX2;
                                                                                        d4s_NTDataArray[6]     =(Float_t)VerdistX_DE;
                                                                                        d4s_NTDataArray[7]     =(Float_t)VerdistX0;
                                                                                        d4s_NTDataArray[8]     =(Float_t)refMultA;
                                                                                        d4s_NTDataArray[9]     =(Float_t)PMassABCDE;
                                                                                        d4s_NTDataArray[10]    =(Float_t)VerdistY;
                                                                                        d4s_NTDataArray[11]    =(Float_t)dca_omega;
                                                                                        d4s_NTDataArray[12]    =(Float_t)VerdistY_DE;
                                                                                        d4s_NTDataArray[13]    =(Float_t)dcaA;
                                                                                        d4s_NTDataArray[14]    =(Float_t)dcaB;
                                                                                        d4s_NTDataArray[15]    =(Float_t)dcaC;
                                                                                        d4s_NTDataArray[16]    =(Float_t)dcaD;
                                                                                        d4s_NTDataArray[17]    =(Float_t)dcaE;
                                                                                        d4s_NTDataArray[18]    =(Float_t)dca_Xi_to_DE;
                                                                                        d4s_NTDataArray[19]    =(Float_t)dcaDE;
                                                                                        d4s_NTDataArray[20]    =(Float_t)dcaAB;
                                                                                        d4s_NTDataArray[21]    =(Float_t)dcaBC_f;
                                                                                        d4s_NTDataArray[22]    =(Float_t)dca_d4s;
                                                                                        d4s_NTDataArray[23]    =(Float_t)refMultB;
                                                                                        // dca_Xi_to_DE

                                                                                        d4s_NT->Fill(d4s_NTDataArray);

                                                                                        cout << "refMultA = " << refMultA << ", refMultB = " << refMultB << endl;
                                                                                    }


                                                                                    /*
                                                                                     alexV0_track_A = alexV0_event_A.createTrack();
                                                                                     alexV0_track_A->setm2A(Mass2ACorr);
                                                                                     alexV0_track_A->setm2B(Mass2BCorr);
                                                                                     alexV0_track_A->setm2C(Mass2CCorr);
                                                                                     alexV0_track_A->setnsA(nSigmaPA);
                                                                                     alexV0_track_A->setnsB(nSigmaPiB);
                                                                                     alexV0_track_A->setnsC(nSigmaC);
                                                                                     alexV0_track_A->setdcaA(dcaA);
                                                                                     alexV0_track_A->setdcaB(dcaB);
                                                                                     alexV0_track_A->setdcaC(dcaC);
                                                                                     alexV0_track_A->setiQxA(iQxA);
                                                                                     alexV0_track_A->setiQyA(iQyA);
                                                                                     alexV0_track_A->setiQxB(iQxB);
                                                                                     alexV0_track_A->setiQyB(iQyB);
                                                                                     alexV0_track_A->setiQxC(iQxC);
                                                                                     alexV0_track_A->setiQyC(iQyC);
                                                                                     alexV0_track_A->setetaA(etaA_c);
                                                                                     alexV0_track_A->setetaB(etaB_c);
                                                                                     alexV0_track_A->setetaC(etaC_c);
                                                                                     alexV0_track_A->setInvAB(InvMassAB);
                                                                                     alexV0_track_A->setInvABC(InvMassABC);
                                                                                     alexV0_track_A->setInvAB_miss(InvMassAB_K0S);
                                                                                     alexV0_track_A->setInvABC_miss(InvMassABC_Xi);
                                                                                     alexV0_track_A->setdcaAB(dcaAB_f);
                                                                                     alexV0_track_A->setdcaBC(dcaBC_f);
                                                                                     alexV0_track_A->setdcaABC(dip_omega);
                                                                                     alexV0_track_A->setVerdistX(VerdistX);
                                                                                     alexV0_track_A->setVerdistY(VerdistY);
                                                                                     alexV0_track_A->setVerdistX2(VerdistX2);
                                                                                     alexV0_track_A->setVerdistY2(dca_omega);
                                                                                     alexV0_track_A->setpt(pt2);
                                                                                     alexV0_track_A->setrap(rap2);
                                                                                     alexV0_track_A->setphi(phiABC);
                                                                                     alexV0_track_A->settheta(thetaABC);
                                                                                     alexV0_track_A->setPsi_ep(phi_event_plane);
                                                                                     alexV0_track_A->setPsi_ep_eta(phi_event_plane_eta_gap);
                                                                                     alexV0_track_A->setPsi_diff_ME(delta_phi_ME_AB_weight);
                                                                                     alexV0_track_A->setscal_prod(scalarProduct);
                                                                                     alexV0_track_A->setscal_prod2(scalarProduct2);
                                                                                     */

                                                                                    dummy_counter_loop++;
                                                                                }
                                                                            }
                                                                        }

                                                                    }
                                                                }
                                                            } // End of m-loop (K-)

                                                        }
                                                    }
                                                } // End of l-loop (pi+)

                                            }
                                        }
                                    }
                                } // End of k-loop (2nd pi-)
                            }
                        }
                    }
                } // End of j-loop (1st pi-)

            }
        } // End of i-loop (p)

        if(ParticleA == 14 && ParticleB == 9)
        {
            //Tree_OmegaMV0_v2  ->Fill();
            //Tree_XiMV0_v2     ->Fill();
        }
        if(ParticleA == 15 && ParticleB == 8)
        {
            //Tree_OmegaPV0_v2  ->Fill();
            //Tree_XiPV0_v2     ->Fill();
        }
        return 1;
    }
    else return 0;
}



//----------------------------------------------------------------------------------------------------------------------------------------------
Int_t d4s_analysis_new(Int_t PID_counter_Array_A[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs],
                       Int_t PID_Array_A[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                       StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                       Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t ParticleD, Int_t ParticleE,
                       Int_t Ana_Num,Int_t SE_ME_Flag, Int_t rot_background_flag,
                       Int_t SE_ME_flagC
                      )
{
    // d4s(????)
    // mass = ???? MeV/c2
    // width = ?? +/- ?? MeV
    //
    // Au + Au -> d4s(????)- + N + N
    //             |
    //             -> Xi-                  +     Lambda         +        pi0
    //                |                            |                      |
    //                -> Lambda + pi-(C)           -> p(D) + pi-(E)       -> not measurable
    //                  |
    //                  -> p(A) + pi-(B)
    //
    // (A,B,C) -> eventA
    // (D,E)   -> eventB


    // Event vertex information
    StThreeVectorF vectorprim,vectorprimB,vectordiff, vectornewA_lin, vectornewB_lin;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB,RunIdA,RunIdB;

    event_SE_ME_ana[0] = picoDst_A;
    event_SE_ME_ana[1] = picoDst_B;
    Int_t PID_counter_Array_SE_ME[2][N_Ana][N_max_PIDs]; // counter array for PID_Array
    Int_t PID_Array_SE_ME[2][N_Ana][N_max_PIDs][N_max_tracks]; // [Event A, Event B][Ana][PID][number] stores the track numbers for each PID
    for(Int_t i_ana = 0; i_ana < N_Ana; i_ana++)
    {
        for(Int_t i_pid = 0; i_pid < N_max_PIDs; i_pid++)
        {
            PID_counter_Array_SE_ME[0][i_ana][i_pid] = PID_counter_Array_A[i_ana][i_pid];
            PID_counter_Array_SE_ME[1][i_ana][i_pid] = PID_counter_Array_B[i_ana][i_pid];
            for(Int_t i_track = 0; i_track < PID_counter_Array_SE_ME[0][i_ana][i_pid]; i_track++)
            {
                PID_Array_SE_ME[0][i_ana][i_pid][i_track] = PID_Array_A[i_ana][i_pid][i_track];
            }
            for(Int_t i_track = 0; i_track < PID_counter_Array_SE_ME[1][i_ana][i_pid]; i_track++)
            {
                PID_Array_SE_ME[1][i_ana][i_pid][i_track] = PID_Array_B[i_ana][i_pid][i_track];
            }
        }
    }

    //event_A_ana       = picoDst_A;
    //event_B_ana       = picoDst_B;

    EventVertexXA     = event_SE_ME_ana[0]->primaryVertex().x();
    EventVertexYA     = event_SE_ME_ana[0]->primaryVertex().y();
    EventVertexZA     = event_SE_ME_ana[0]->primaryVertex().z();

    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_SE_ME_ana[0]->refMult();
    refMultB          = refMultA;
    RunIdA            = event_SE_ME_ana[0]->runId();
    RunIdB            = RunIdA;
    Float_t ZDCx      = event_SE_ME_ana[0]->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx      = event_SE_ME_ana[0]->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_SE_ME_ana[0]->vzVpd();

    vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    Int_t flag_tof_particleA = 0;
    Int_t flag_tof_particleB = 0;
    Int_t flag_tof_particleC = 0;
    Int_t flag_tof_particleD = 0;
    Int_t flag_tof_particleE = 0;

    Int_t d4s_reco_counter = 0;

    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        EventVertexXB  = event_SE_ME_ana[1]->primaryVertex().x();
        EventVertexYB  = event_SE_ME_ana[1]->primaryVertex().y();
        EventVertexZB  = event_SE_ME_ana[1]->primaryVertex().z();
        refMultB       = event_SE_ME_ana[1]->refMult();
        RunIdB         = event_SE_ME_ana[1]->runId();
        vectorprimB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vectorprim - vectorprimB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }


    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
    }

    //
    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       &&
       (
        (
         SE_ME_Flag == 0
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleA]    > 0   // p(A)
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleB]    > 0   // pi-(B)
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleC]    > 0   // pi-(C)
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleD]    > 0   // p(D)
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleE]    > 0   // pi-(E)
        )
        ||
        (
         SE_ME_Flag == 1
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleA]    > 0   // p(A)
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleB]    > 0   // pi-(B)
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleC]    > 0   // pi-(C)
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleD]    > 0   // p(D)
         && PID_counter_Array_SE_ME[0][Ana_Num][ParticleE]    > 0   // pi-(E)
         && PID_counter_Array_SE_ME[1][Ana_Num][ParticleA]    > 0   // p(A)
         && PID_counter_Array_SE_ME[1][Ana_Num][ParticleB]    > 0   // pi-(B)
         && PID_counter_Array_SE_ME[1][Ana_Num][ParticleC]    > 0   // pi-(C)
         && PID_counter_Array_SE_ME[1][Ana_Num][ParticleD]    > 0   // p(D)
         && PID_counter_Array_SE_ME[1][Ana_Num][ParticleE]    > 0   // pi-(E)
        )
       )
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_SE_ME_ana[0]   ->isMinBias()
       && event_SE_ME_ana[1]   ->isMinBias()
      )
    {
        // Initial cuts -- d4s new
        Float_t dcaA_cut            = 0.2;
        Float_t dcaB_cut            = 0.9;
        Float_t dcaC_cut            = 0.5;
        Float_t dcaA_upper_cut      = 20.0;
        Float_t dcaB_upper_cut      = 20.0;
        Float_t dcaC_upper_cut      = 20.0;
        Float_t dcaAB_cut           = 2.0;
        Float_t InvMassABC_cut      = 1.3273;
        Float_t InvMassABC_low_cut  = (1.3157-50.0);
        Float_t dcaABC_cut          = 2.0;

        Float_t VerdistX_AB_cut     = 3.5;
        Float_t VerdistY_AB_cut     = 0.2;
        Float_t VerdistY_AB_upper_cut = 10.0;
        Float_t VerdistX_ABC_cut    = 2.0;
        Float_t VerdistY_ABC_cut    = 10.0; // 3.0

        Float_t Lambda_m_upper_cut  = 1.14; // 1.12
        Float_t Lambda_m_lower_cut  = 1.09; // 1.11

        Float_t proton_m2_upper_cut = 1.5;
        Float_t proton_m2_lower_cut = 0.4;
        Float_t pion_m2_upper_cut   = 0.1;
        Float_t pion_m2_lower_cut   = -0.1;
        Float_t kaon_m2_upper_cut   = 0.4;
        Float_t kaon_m2_lower_cut   = 0.1;

        Float_t nHitsFit_cut        = 14.0;
        Float_t nHitsPoss_cut       = 0.0;
        Float_t Momentum_upper_cut  = 10.0;
        Float_t Momentum_lower_cut  = 0.15;
        Float_t nHits_ratio_cut     = 0.52;

        Float_t InvMassAB_use       = 1.115683; // Lambda mass
        Float_t InvMassABC_use      = 1.32131;  // Xi+/- mass
        Float_t InvMassAB_upper_cut = Lambda_m_upper_cut;
        Float_t InvMassAB_lower_cut = Lambda_m_lower_cut;

        if(ParticleC == 11 || ParticleC == 12)
        {
            InvMassABC_cut          = 1.69;
            InvMassABC_low_cut      = 1.64;
            InvMassAB_use           = 1.115683; // Lambda mass
            InvMassABC_use          = 1.67245; // Omega mass
            InvMassAB_upper_cut     = Lambda_m_upper_cut;
            InvMassAB_lower_cut     = Lambda_m_lower_cut;
        }

        if(
           (ParticleA == 8 && ParticleB == 9) ||
           (ParticleA == 9 && ParticleB == 8)
          )
        {
            dcaA_cut                = 0.65;
            dcaB_cut                = 0.65;
            VerdistX_AB_cut         = 3.0;
            VerdistY_AB_cut         = 0.2;
            VerdistY_AB_upper_cut   = 2.5;
            VerdistX_ABC_cut        = 2.5;
            InvMassAB_use           = 0.497648; // K0S mass
            InvMassABC_use          = 0.9; // some mass
            InvMassABC_cut          = InvMassABC_use+0.01;
            InvMassABC_low_cut      = InvMassABC_use-0.01;
            InvMassAB_upper_cut     = 0.52;
            InvMassAB_lower_cut     = 0.48;
            dcaAB_cut               = 1.0;
            dcaABC_cut              = 1.5;
        }
        Float_t InvMassDE_use       = InvMassAB_use;


        StPhysicalHelixD helixA, helixB, helixC, helixD, helixE, helix_ABC;
        StThreeVectorF vectorA, vectorB, vectoratsA, vectoratsB, vectorAB, vectorDE, vectorABC, vectorAB_est, vectorprimAB, vectornewA, vectornewB, vectornewD, vectornewE, dirY_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY, dirY_DE, tv3_perpAB, tv3_perpDE;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackD, ltrackE, ltrackC_pi, ltrackB2, ltrackA_lin, ltrackB_lin, ltrackB2_pi, ltrackA_pip;
        TLorentzVector trackAB, trackDE, trackABC, trackAB_test;
        StThreeVectorF Stv3_perpAB, Stv3_perpDE;


        // calculate event plane angle
        Float_t phi_event_plane_d4s = calc_phi_event_plane_2nd(EP_Qx,EP_Qy);


        if(0)
        {
            cout << "d4s event accepted, A = " << PID_counter_Array_SE_ME[0][Ana_Num][ParticleA] <<
                ", B = " << PID_counter_Array_SE_ME[0][Ana_Num][ParticleB] <<
                ", C = " << PID_counter_Array_SE_ME[0][Ana_Num][ParticleC] <<
                ", D = " << PID_counter_Array_SE_ME[SE_ME_Flag][Ana_Num][ParticleD] <<
                ", E = " << PID_counter_Array_SE_ME[SE_ME_Flag][Ana_Num][ParticleE] <<
                ", erefMult_bin = " << erefMult_bin <<
                endl;
        }

        const Int_t N_d4s = 200;
        const Int_t N_max_track = 5000;
        Int_t Index_LambdaXi[N_d4s];
        Int_t Index_Lambda[N_d4s];
        Int_t Index_Xi[N_d4s];
        Int_t Index_Track[N_max_track];

        Int_t Index_LambdaXi_reco[N_d4s];
        Int_t Index_Lambda_reco[N_d4s];
        Int_t Index_Xi_reco[N_d4s];
        Int_t Index_Track_reco[N_max_track];

        Int_t Index_array_Lambda[N_d4s];
        Int_t Index_array_Xi[N_d4s];
        Int_t Index_array_track_evA[N_max_track];
        Int_t Index_array_track_evB[N_max_track];
        Int_t Index_array_total_track_evA[N_max_track];
        Int_t Index_array_total_track_evB[N_max_track];

        for(Int_t i = 0; i < N_d4s; i++)
        {
            Index_LambdaXi[i] = -1;
            Index_Lambda[i]   = -1;
            Index_Xi[i]       = -1;

            Index_LambdaXi_reco[i] = -1;
            Index_Lambda_reco[i]   = -1;
            Index_Xi_reco[i]       = -1;

            Index_array_Lambda[i]   = -1;
            Index_array_Xi[i]       = -1;
        }
        for(Int_t i = 0; i < N_max_track; i++)
        {
            Index_Track[i]       = -1;
            Index_Track_reco[i]  = -1;
            Index_array_track_evA[i] = -1;
            Index_array_track_evB[i] = -1;
            Index_array_total_track_evA[i] = -1;
            Index_array_total_track_evB[i] = -1;
        }
        Int_t Lambda_all_reco_counter = 0;
        Int_t LambdaXi_reco_counter   = 0;
        Int_t Lambda_reco_counter     = 0;
        Int_t Xi_reco_counter         = 0;
        Int_t Track_reco_counter      = 0;

        Float_t mass_cut_use_upper    = 0.0;
        Float_t mass_cut_use_lower    = 0.0;

        //------------------------------------------------------------------------------------------------------------
        // Lambda/K0S reconstruction
        Int_t Lambda_reco_counter_SE_ME = 0;
        Int_t N_SE_ME = 1; // for same event only one Lambda loop will be performed
        // --> if(SE_ME_Flag == 1) N_SE_ME = 4; // for mixed event Lambdas have to be reconstructed from event A and B
        if(SE_ME_Flag == 1) N_SE_ME = 2; //

        // setup all possible combinations for mixing particle A and particle B from events A and B
        // --> Int_t particle_A_SE_ME[4] = {0,1,0,1};
        // --> Int_t particle_B_SE_ME[4] = {0,1,1,0};
        Int_t particle_A_SE_ME[2] = {0,1};
        Int_t particle_B_SE_ME[2] = {0,1};

        StThreeVectorF vectorprim_add[2];
        vectorprim_add[0].set(0.0,0.0,0.0);
        vectorprim_add[1] = vectordiff;

        for(Int_t iSE_ME = 0; iSE_ME < N_SE_ME; iSE_ME++)
        {
            for(Int_t i = 0; i < PID_counter_Array_SE_ME[particle_A_SE_ME[iSE_ME]][Ana_Num][ParticleA]; i++) // p candidates
            {
                Int_t trackA_num;
                StPicoAlexTrack trackA;
                trackA_num = PID_Array_SE_ME[particle_A_SE_ME[iSE_ME]][Ana_Num][ParticleA][i];
                trackA     = *event_SE_ME_ana[particle_A_SE_ME[iSE_ME]]->track( trackA_num );

                Float_t MomentumA   = trackA.gMom().mag();
                Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
                Float_t nHitsPossA  = trackA.nHitsMax();
                Float_t nHitsFitA   = trackA.nHitsFit();
                Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
                Float_t nSigmaPA    = trackA.nSigmaProton();
                Float_t nSigmaAPi   = trackA.nSigmaPion();
                Float_t nSigmaAK    = trackA.nSigmaKaon();
                Float_t nSigmaAP    = trackA.nSigmaProton();
                if(ParticleA == 8  || ParticleA == 9)
                {
                    nSigmaPA = trackA.nSigmaPion();
                    mass_cut_use_upper = pion_m2_upper_cut;
                    mass_cut_use_lower = pion_m2_lower_cut;
                }
                if(ParticleA == 11 || ParticleA == 12)
                {
                    nSigmaPA = trackA.nSigmaKaon();
                    mass_cut_use_upper = kaon_m2_upper_cut;
                    mass_cut_use_lower = kaon_m2_lower_cut;
                }
                if(ParticleA == 14 || ParticleA == 15)
                {
                    nSigmaPA = trackA.nSigmaProton();
                    mass_cut_use_upper = proton_m2_upper_cut;
                    mass_cut_use_lower = proton_m2_lower_cut;
                }
                Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                Float_t Mass2A      = -100.0;

                // calculate mass2
                if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
                {
                    flag_tof_particleA = 1;
                    Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
                }
                else
                {
                    flag_tof_particleA = 0;
                }
                Float_t MassA       = SquareRoot(Mass2A);

                if(
                   dcaA          > dcaA_cut
                   && dcaA       < dcaA_upper_cut
                   && nHitsFitA  > nHitsFit_cut
                   && nHitsPossA > nHitsPoss_cut
                   && MomentumA  > Momentum_lower_cut
                   && MomentumA  < Momentum_upper_cut
                   && (nHitsFitA/nHitsPossA) > nHits_ratio_cut
                   && (flag_tof_particleA == 0 || (flag_tof_particleA == 1 && Mass2A > mass_cut_use_lower && Mass2A < mass_cut_use_upper)) // proton mass
                  )
                {
                    helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin()+vectorprim_add[particle_A_SE_ME[iSE_ME]],event_SE_ME_ana[particle_A_SE_ME[iSE_ME]]->bField()*MAGFIELDFACTOR, trackA.charge());


                    for(Int_t j = 0; j < PID_counter_Array_SE_ME[particle_B_SE_ME[iSE_ME]][Ana_Num][ParticleB]; j++) // pi- candidates
                    {
                        Int_t trackB_num;
                        StPicoAlexTrack trackB;
                        trackB_num = PID_Array_SE_ME[particle_B_SE_ME[iSE_ME]][Ana_Num][ParticleB][j];
                        trackB     = *event_SE_ME_ana[particle_B_SE_ME[iSE_ME]]->track( trackB_num );

                        Float_t MomentumB   = trackB.gMom().mag();
                        Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                        Float_t nHitsPossB  = trackB.nHitsMax();
                        Float_t nHitsFitB   = trackB.nHitsFit();
                        Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                        Float_t nSigmaPiB   = trackB.nSigmaPion();
                        Float_t nSigmaBPi   = trackB.nSigmaPion();
                        Float_t nSigmaBK    = trackB.nSigmaKaon();
                        Float_t nSigmaBP    = trackB.nSigmaProton();
                        if(ParticleB == 8  || ParticleB == 9)
                        {
                            nSigmaPiB = trackB.nSigmaPion();
                            mass_cut_use_upper = pion_m2_upper_cut;
                            mass_cut_use_lower = pion_m2_lower_cut;
                        }
                        if(ParticleB == 11 || ParticleB == 12)
                        {
                            nSigmaPiB = trackB.nSigmaKaon();
                            mass_cut_use_upper = kaon_m2_upper_cut;
                            mass_cut_use_lower = kaon_m2_lower_cut;
                        }
                        if(ParticleB == 14 || ParticleB == 15)
                        {
                            nSigmaPiB = trackB.nSigmaProton();
                            mass_cut_use_upper = proton_m2_upper_cut;
                            mass_cut_use_lower = proton_m2_lower_cut;
                        }

                        Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                        Float_t Mass2B      = -100.0;

                        // calculate mass2
                        if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                        {
                            flag_tof_particleB = 1;
                            Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                        }
                        else
                        {
                            flag_tof_particleB = 0;
                        }
                        Float_t MassB       = SquareRoot(Mass2B);

                        if(
                           (
                            (particle_A_SE_ME[iSE_ME] == particle_B_SE_ME[iSE_ME] && trackA_num != trackB_num)
                            || (particle_A_SE_ME[iSE_ME] != particle_B_SE_ME[iSE_ME])
                           ) // Prevent that a track is used twice
                           && dcaB       > dcaB_cut
                           && dcaB       < dcaB_upper_cut
                           && nHitsFitB  > nHitsFit_cut
                           && nHitsPossB > nHitsPoss_cut
                           && MomentumB  > Momentum_lower_cut
                           && MomentumB  < Momentum_upper_cut
                           && (nHitsFitB/nHitsPossB) > nHits_ratio_cut
                           && (flag_tof_particleB == 0 || (flag_tof_particleB == 1 && Mass2B > mass_cut_use_lower && Mass2B < mass_cut_use_upper)) // pion mass
                          )
                        {
                            helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin()+vectorprim_add[particle_B_SE_ME[iSE_ME]],event_SE_ME_ana[particle_B_SE_ME[iSE_ME]]->bField()*MAGFIELDFACTOR, trackB.charge());

                            Float_t VerdistX, dcaAB;
                            Int_t decay_prop_AB = calc_Decay_Properties(helixA,helixB,mass_array[ParticleA],mass_array[ParticleB],MomentumA,
                                                                        MomentumB,vectorprim,VerdistX_AB_cut,dcaAB_cut,vectorAB,
                                                                        ltrackA,ltrackB,trackAB,VerdistX,dcaAB);


                            // calculate decay plane angle of Lambda candidate
                            TVector3 tv3_trackA   = ltrackA.Vect();
                            TVector3 tv3_trackB   = ltrackB.Vect();
                            TVector3 tv3_perpAB   = tv3_trackA.Cross(tv3_trackB);
                            Stv3_perpAB.set(tv3_perpAB.X(),tv3_perpAB.Y(),tv3_perpAB.Z());


                            if( decay_prop_AB == 1 )
                            {
                                Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                                Float_t  MomentumAB         = trackAB.P(); // momentum of mother particle
                                Float_t  BetaAB = TMath::Sqrt(1./(1+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                                dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                                dirY *= 1.0/dirY.mag();

                                baseY = vectorAB;
                                Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vectorprim);

                                if(0)
                                {
                                    cout << "VerdistY = " << VerdistY << endl;
                                }
                                // Apply Lambda cuts -- new
                                if(
                                   InvMassAB         > InvMassAB_lower_cut
                                   && InvMassAB      < InvMassAB_upper_cut
                                   && dcaAB          < dcaAB_cut
                                   && VerdistX       > VerdistX_AB_cut
                                   && VerdistY       > VerdistY_AB_cut
                                   && VerdistY       < VerdistY_AB_upper_cut
                                  )
                                {
                                    if(Lambda_reco_counter_SE_ME < N_d4s_Lambda)
                                    {
                                        if(0)
                                        {
                                            cout << "iSE_ME = " << iSE_ME << ", A_SE_ME = " << particle_A_SE_ME[iSE_ME] << ", B_SE_ME = " << particle_B_SE_ME[iSE_ME]
                                                << ", A_mult = " << event_SE_ME_ana[particle_A_SE_ME[iSE_ME]]->refMult() << ", B_mult = " << event_SE_ME_ana[particle_B_SE_ME[iSE_ME]]->refMult() << endl;
                                        }

                                        d4s_Lambda_helix[0][Lambda_reco_counter_SE_ME]      = helixA;
                                        d4s_Lambda_helix[1][Lambda_reco_counter_SE_ME]      = helixB;
                                        d4s_Lambda_vertex[0][Lambda_reco_counter_SE_ME]     = vectorAB;
                                        d4s_Lambda_vertex[1][Lambda_reco_counter_SE_ME]     = Stv3_perpAB;
                                        d4s_Lambda_TLV[0][Lambda_reco_counter_SE_ME]        = ltrackA;
                                        d4s_Lambda_TLV[1][Lambda_reco_counter_SE_ME]        = ltrackB;
                                        d4s_Lambda_TLV[2][Lambda_reco_counter_SE_ME]        = trackAB;
                                        d4s_Lambda_properties[0][Lambda_reco_counter_SE_ME] = VerdistX;
                                        d4s_Lambda_properties[1][Lambda_reco_counter_SE_ME] = VerdistY;
                                        d4s_Lambda_properties[2][Lambda_reco_counter_SE_ME] = dcaAB;
                                        d4s_Lambda_properties[3][Lambda_reco_counter_SE_ME] = iSE_ME;
                                        d4s_Lambda_daughters[0][Lambda_reco_counter_SE_ME]  = nSigmaPA;
                                        d4s_Lambda_daughters[1][Lambda_reco_counter_SE_ME]  = nSigmaPiB;
                                        d4s_Lambda_daughters[2][Lambda_reco_counter_SE_ME]  = Mass2A;
                                        d4s_Lambda_daughters[3][Lambda_reco_counter_SE_ME]  = Mass2B;
                                        d4s_Lambda_daughters[4][Lambda_reco_counter_SE_ME]  = trackA_num;
                                        d4s_Lambda_daughters[5][Lambda_reco_counter_SE_ME]  = trackB_num;
                                        d4s_Lambda_daughters[6][Lambda_reco_counter_SE_ME]  = dcaA;
                                        d4s_Lambda_daughters[7][Lambda_reco_counter_SE_ME]  = dcaB;
                                        d4s_Lambda_daughters[8][Lambda_reco_counter_SE_ME]  = PolarityA*MomentumA;
                                        d4s_Lambda_daughters[9][Lambda_reco_counter_SE_ME]  = PolarityB*MomentumB;
                                        d4s_Lambda_daughters[10][Lambda_reco_counter_SE_ME] = nHitsFitA;
                                        d4s_Lambda_daughters[11][Lambda_reco_counter_SE_ME] = nHitsFitB;
                                        d4s_Lambda_daughters[12][Lambda_reco_counter_SE_ME] = nSigmaAPi;
                                        d4s_Lambda_daughters[13][Lambda_reco_counter_SE_ME] = nSigmaAK;
                                        d4s_Lambda_daughters[14][Lambda_reco_counter_SE_ME] = nSigmaAP;
                                        d4s_Lambda_daughters[15][Lambda_reco_counter_SE_ME] = nSigmaBPi;
                                        d4s_Lambda_daughters[16][Lambda_reco_counter_SE_ME] = nSigmaBK;
                                        d4s_Lambda_daughters[17][Lambda_reco_counter_SE_ME] = nSigmaBP;

                                        if(0)
                                        {
                                            cout << "trackA_num = " << trackA_num << ", trackB_num = " << trackB_num << ", M = " << trackAB.M() << endl;
                                        }
                                        //cout << "iLambdaXi = " << Lambda_reco_counter_SE_ME << ", trackAB = {" << trackAB.Px() << ", " << trackAB.Py() << ", " << trackAB.Pz() << endl;
                                        //trackAB_test         = *d4s_Lambda_TLV[0][Lambda_reco_counter_SE_ME];
                                        //cout << "iLambdaXi = " << Lambda_reco_counter_SE_ME << ", trackAB_test = {" << trackAB_test.Px() << ", " << trackAB_test.Py() << ", " << trackAB_test.Pz() << endl;
                                        h_d4s_InvMassAB[iSE_ME]->Fill(trackAB.M());

                                        Lambda_reco_counter_SE_ME++;
                                    }
                                    else
                                    {
                                        //cout << "Warning: maximum number of Lambda candidates reached" << endl;
                                    }
                                    //cout << "Lambda accepted: " << Lambda_reco_counter_SE_ME[iSE_ME] << ", InvMassAB = " << InvMassAB << ", VerdistX = " << VerdistX
                                    //    << ", VerdistY = " << VerdistY << ", dcaAB = " << dcaAB << ", i = " << i << ", j = " << j << endl;
                                    //cout << "pA = " << MomentumA << ", dcaA = " << dcaA << ", m2A = " << Mass2A << endl;
                                    //cout << "pB = " << MomentumB << ", dcaB = " << dcaB << ", m2B = " << Mass2B << endl;
                                }
                            }
                        }
                    }
                }
            }
        }

        //cout << "N_Lambdas_SE = " << Lambda_reco_counter_SE_ME[0] << endl;
        // Lambda reconstruction done
        //------------------------------------------------------------------------------------------------------------



        //------------------------------------------------------------------------------------------------------------
        // Xi reconstruction
        Int_t Xi_reco_counter_loop = 0;
        for(Int_t iLambda = 0; iLambda < Lambda_reco_counter_SE_ME; iLambda++) // Lambda loop
        {
            // Get Lambda properties
            vectorAB              = d4s_Lambda_vertex[0][iLambda];
            trackAB               = d4s_Lambda_TLV[2][iLambda];
            Float_t VerdistX_AB   = d4s_Lambda_properties[0][iLambda];
            Float_t VerdistY_AB   = d4s_Lambda_properties[1][iLambda];
            Float_t dcaAB         = d4s_Lambda_properties[2][iLambda];
            Int_t   iSE_ME_Lambda = d4s_Lambda_properties[3][iLambda];
            Float_t nSigmaPA      = d4s_Lambda_daughters[0][iLambda];
            Float_t nSigmaPiB     = d4s_Lambda_daughters[1][iLambda];
            Float_t Mass2A        = d4s_Lambda_daughters[2][iLambda];
            Float_t Mass2B        = d4s_Lambda_daughters[3][iLambda];
            Int_t trackA_num      = (Int_t)d4s_Lambda_daughters[4][iLambda];
            Int_t trackB_num      = (Int_t)d4s_Lambda_daughters[5][iLambda];
            Float_t dcaA          = d4s_Lambda_daughters[6][iLambda];
            Float_t dcaB          = d4s_Lambda_daughters[7][iLambda];
            Float_t qpA           = d4s_Lambda_daughters[8][iLambda];
            Float_t qPB           = d4s_Lambda_daughters[9][iLambda];
            //cout << "iLambda = " << iLambda << ", trackAB = {" << trackAB.Px() << ", " << trackAB.Py() << ", " << trackAB.Pz() << endl;


            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
            dirY *= 1.0/dirY.mag();
            baseY = vectorAB;

            for(Int_t k = 0; k < PID_counter_Array_SE_ME[SE_ME_flagC][Ana_Num][ParticleC]; k++)  // pi-
            {
                Int_t trackC_num = PID_Array_SE_ME[SE_ME_flagC][Ana_Num][ParticleC][k];
                if(
                   iSE_ME_Lambda == SE_ME_flagC &&
                   !(SE_ME_flagC == particle_A_SE_ME[iSE_ME_Lambda] && trackC_num == trackA_num) &&
                   !(SE_ME_flagC == particle_B_SE_ME[iSE_ME_Lambda] && trackC_num == trackB_num)
                  )
                {
                    StThreeVectorF vectorC, vectoratsC, vectoratsA2, vectorA2C, vectorA2CtoPrim, vectornewC, dirY2, baseY2;

                    StPicoAlexTrack trackC  = *event_SE_ME_ana[SE_ME_flagC]->track( trackC_num );
                    Float_t MomentumC   = trackC.gMom().mag();
                    Float_t dcaC        = trackC.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossC  = trackC.nHitsMax();
                    Float_t nHitsFitC   = trackC.nHitsFit();
                    Float_t BetaC       = trackC.btofBeta();  // Velocity after time-of-flight reconstruction
                    Float_t nSigmaC     = trackC.nSigmaPion(); // for Xi analysis
                    Float_t nSigmaCPi   = trackC.nSigmaPion();
                    Float_t nSigmaCK    = trackC.nSigmaKaon();
                    Float_t nSigmaCP    = trackC.nSigmaProton();
                    if(ParticleC == 8  || ParticleC == 9)
                    {
                        nSigmaC = trackC.nSigmaPion();
                        mass_cut_use_upper = pion_m2_upper_cut;
                        mass_cut_use_lower = pion_m2_lower_cut;
                    }
                    if(ParticleC == 11 || ParticleC == 12)
                    {
                        nSigmaC = trackC.nSigmaKaon();
                        mass_cut_use_upper = kaon_m2_upper_cut;
                        mass_cut_use_lower = kaon_m2_lower_cut;
                    }
                    if(ParticleC == 14 || ParticleC == 15)
                    {
                        nSigmaC = trackC.nSigmaProton();
                        mass_cut_use_upper = proton_m2_upper_cut;
                        mass_cut_use_lower = proton_m2_lower_cut;
                    }

                    Float_t PolarityC   = trackC.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    Float_t Mass2C      = -100.0;
                    // calculate mass2
                    if(trackC.btofMatchFlag() > 0 && trackC.btof() != 0 && BetaC != 0)
                    {
                        flag_tof_particleC = 1;
                        Mass2C = MomentumC*MomentumC*(1.0/(BetaC*BetaC) - 1.0);
                    }
                    else
                    {
                        flag_tof_particleC = 0;
                    }
                    Float_t MassC       = SquareRoot(Mass2C);

                    if(
                       dcaC          > dcaC_cut
                       && dcaC       < dcaC_upper_cut
                       && nHitsFitC  > nHitsFit_cut
                       && nHitsPossC > nHitsPoss_cut
                       && MomentumC  > Momentum_lower_cut
                       && MomentumC  < Momentum_upper_cut
                       && (nHitsFitC/nHitsPossC) > nHits_ratio_cut
                       && (flag_tof_particleC == 0 || (flag_tof_particleC == 1 && Mass2C > mass_cut_use_lower && Mass2C < mass_cut_use_upper))
                      )
                    {
                        helixC = StPhysicalHelixD(trackC.gMom(),trackC.origin()+vectorprim_add[SE_ME_flagC],event_SE_ME_ana[SE_ME_flagC]->bField()*MAGFIELDFACTOR,trackC.charge());

                        Float_t pathA2_f, pathC_f,dcaBC_f;
                        fHelixAtoLinedca(dirY,baseY,helixC,pathC_f,pathA2_f,dcaBC_f); // calculates dca of second pi- helix to Lambda track

                        vectoratsA2       = baseY+pathA2_f*dirY;  // space vector of Lambda at dca to helixC
                        vectoratsC        = helixC.at(pathC_f);  // space vector of helixC at dca to Lambda
                        vectorA2C         = vectoratsA2+vectoratsC;
                        vectorA2C         = vectorA2C/2.0; // decay vertex
                        vectorA2CtoPrim   = vectorA2C - vectorprim; // vector primary vertex to decay vertex
                        Float_t VerdistX_ABC = vectorA2CtoPrim.mag(); // distance between primary vertex and decay vertex

                        vectornewC       = helixC.cat(pathC_f); // direction vector at dca for helixC
                        vectornewC       = MomentumC*vectornewC/vectornewC.mag(); // new momentum vector at decay vertex
                        ltrackC.SetXYZM(vectornewC.x(),vectornewC.y(),vectornewC.z(),mass_array[ParticleC]);

                        TLorentzVector trackAB_Lambda; // mother particle with nominal mother particle mass
                        trackAB_Lambda.SetXYZM(trackAB.Px(),trackAB.Py(),trackAB.Pz(),InvMassAB_use);
                        TLorentzVector trackLambdaC    = trackAB_Lambda+ltrackC; // mother particle

                        Double_t InvMassABC      = trackLambdaC.M(); // invariant mass of mother particle
                        Double_t MomentumABC     = trackLambdaC.P(); // momentum of mother particle
                        Float_t BetaABC = TMath::Sqrt(1./(1+(InvMassABC/MomentumABC)*(InvMassABC/MomentumABC)));

                        dirY2.set(trackLambdaC.Px(),trackLambdaC.Py(),trackLambdaC.Pz()); // direction vector of Xi
                        dirY2 = dirY2/dirY2.mag();
                        Double_t scalarProduct2 = dirY2.dot(vectorA2CtoPrim/vectorA2CtoPrim.mag());

                        baseY2 = vectorA2C;
                        //Double_t  VerdistY2  = calculateMinimumDistanceStraightToPoint(baseY2,dirY2,vectorprim);

                        Float_t pt2          = trackLambdaC.Pt();  // Transverse momentum of mother particle
                        Float_t rap2         = trackLambdaC.Rapidity(); // Rapidity of mother particle

                        Float_t phiABC   = dirY2.phi();
                        Float_t thetaABC = dirY2.theta();


                        //********************** Calculate the helix of the Omega/Xi ***************************
                        StThreeVectorF origin_ABC;   // decay vertex
                        origin_ABC.set(vectorA2C.x(),vectorA2C.y(),vectorA2C.z());
                        //Float_t h_omega     = trackC.geth(); // +1 for K+ and -1 for K-
                        Float_t h_ABC;
                        if(ParticleC == 8 || ParticleC == 11 || ParticleC == 14) h_ABC = 1.0;
                        if(ParticleC == 9 || ParticleC == 12 || ParticleC == 15) h_ABC = -1.0;
                        Float_t phase_ABC = TMath::ATan2(-1.0*h_ABC*trackLambdaC.Px(),h_ABC*trackLambdaC.Py());
                        Float_t dip_ABC   = TMath::ATan2(trackLambdaC.Pz(),pt2);  // correct
                        Float_t curv_ABC  = 1.0;
                        if(pt2 != 0.0)
                        {
                            curv_ABC = curv_to_invpt_ratio/pt2;
                        }

                        helix_ABC.setParameters(curv_ABC,dip_ABC,phase_ABC,origin_ABC,h_ABC);
                        Float_t path_ABC      = -999.0;
                        Float_t VerdistY_ABC  = -999.0;
                        fHelixAtoPointdca(vectorprim,helix_ABC,path_ABC,VerdistY_ABC);
                        //**************************************************************************************

                        if(
                           InvMassABC        < InvMassABC_cut
                           && InvMassABC     > InvMassABC_low_cut
                           && dcaBC_f        < dcaABC_cut
                           && VerdistX_ABC   > VerdistX_ABC_cut
                           && VerdistY_ABC   < VerdistY_ABC_cut
                           && scalarProduct2 > 0.0
                          )
                        {
                            if(Xi_reco_counter_loop < N_d4s_Lambda)
                            {
                                d4s_Xi_helix[0][Xi_reco_counter_loop]      = helixC;
                                d4s_Xi_helix[1][Xi_reco_counter_loop]      = helix_ABC;
                                d4s_Xi_vertex[Xi_reco_counter_loop]        = vectorA2C;
                                d4s_Xi_TLV[0][Xi_reco_counter_loop]        = ltrackC;
                                d4s_Xi_TLV[1][Xi_reco_counter_loop]        = trackLambdaC;
                                d4s_Xi_properties[0][Xi_reco_counter_loop] = VerdistX_ABC;
                                d4s_Xi_properties[1][Xi_reco_counter_loop] = VerdistY_ABC;
                                d4s_Xi_properties[2][Xi_reco_counter_loop] = dcaBC_f;
                                d4s_Xi_properties[3][Xi_reco_counter_loop] = iLambda;
                                d4s_Xi_daughters[0][Xi_reco_counter_loop]  = nSigmaC;
                                d4s_Xi_daughters[1][Xi_reco_counter_loop]  = Mass2C;
                                d4s_Xi_daughters[2][Xi_reco_counter_loop]  = trackC_num;
                                d4s_Xi_daughters[3][Xi_reco_counter_loop]  = dcaC;
                                d4s_Xi_daughters[4][Xi_reco_counter_loop]  = PolarityC*MomentumC;
                                d4s_Xi_daughters[5][Xi_reco_counter_loop]  = nHitsFitC;
                                d4s_Xi_daughters[6][Xi_reco_counter_loop]  = nSigmaCPi;
                                d4s_Xi_daughters[7][Xi_reco_counter_loop]  = nSigmaCK;
                                d4s_Xi_daughters[8][Xi_reco_counter_loop]  = nSigmaCP;

                                Xi_reco_counter_loop++;
                            }
                            else
                            {
                                //cout << "Warning: maximum number of Xi candidates reached: number of Lambda/K0S candidates = " << Lambda_reco_counter_SE_ME[0] << endl;
                            }
                            //cout << "Xi accepted" << endl;
                        }
                    }
                }
            }
        }
        //cout << "Xi_reco_counter_loop = " << Xi_reco_counter_loop << endl;
        //------------------------------------------------------------------------------------------------------------



        //------------------------------------------------------------------------------------------------------------
        // d4s reconstruction


        //---------------------------------------------
        Int_t Lambda_use_counter        = 0;
        Int_t Track_use_counter_evAB[2] = {0,0}; // for event A and event B
        Int_t Track_use_counter_total   = 0;
        //---------------------------------------------


        for(Int_t iXi = 0; iXi < Xi_reco_counter_loop; iXi++) // Xi loop
        {
            // Get Xi properties
            vectorABC            = d4s_Xi_vertex[iXi];
            helixC               = d4s_Xi_helix[0][iXi];
            helix_ABC            = d4s_Xi_helix[1][iXi];
            ltrackC              = d4s_Xi_TLV[0][iXi];
            trackABC             = d4s_Xi_TLV[1][iXi];
            Float_t VerdistX_ABC = d4s_Xi_properties[0][iXi];
            Float_t VerdistY_ABC = d4s_Xi_properties[1][iXi];
            Float_t dcaABC       = d4s_Xi_properties[2][iXi];
            Int_t iLambdaXi      = (Int_t)d4s_Xi_properties[3][iXi];
            Float_t nSigmaC      = d4s_Xi_daughters[0][iXi];
            Float_t Mass2C       = d4s_Xi_daughters[1][iXi];
            Int_t trackC_num     = (Int_t)d4s_Xi_daughters[2][iXi];
            Float_t dcaC         = d4s_Xi_daughters[3][iXi];
            Float_t qpC          = d4s_Xi_daughters[4][iXi];
            Int_t nhitsFitC      = (Int_t)d4s_Xi_daughters[5][iXi];
            Float_t nSigmaCPi    = d4s_Xi_daughters[6][iXi];
            Float_t nSigmaCK     = d4s_Xi_daughters[7][iXi];
            Float_t nSigmaCP     = d4s_Xi_daughters[8][iXi];
            //cout << "iXi = " << iXi << ", trackABC = {" << trackABC.Px() << ", " << trackABC.Py() << ", " << trackABC.Pz() << endl;

            // Get Lambda (Xi) properties
            vectorAB                = d4s_Lambda_vertex[0][iLambdaXi];
            Stv3_perpAB             = d4s_Lambda_vertex[1][iLambdaXi];
            helixA                  = d4s_Lambda_helix[0][iLambdaXi];
            helixB                  = d4s_Lambda_helix[1][iLambdaXi];
            ltrackA                 = d4s_Lambda_TLV[0][iLambdaXi];
            ltrackB                 = d4s_Lambda_TLV[1][iLambdaXi];
            trackAB                 = d4s_Lambda_TLV[2][iLambdaXi];
            Float_t VerdistX_AB     = d4s_Lambda_properties[0][iLambdaXi];
            Float_t VerdistY_AB     = d4s_Lambda_properties[1][iLambdaXi];
            Float_t dcaAB           = d4s_Lambda_properties[2][iLambdaXi];
            Int_t   iSE_ME_LambdaXi = d4s_Lambda_properties[3][iLambdaXi];
            Float_t nSigmaPA        = d4s_Lambda_daughters[0][iLambdaXi];
            Float_t nSigmaPiB       = d4s_Lambda_daughters[1][iLambdaXi];
            Float_t Mass2A          = d4s_Lambda_daughters[2][iLambdaXi];
            Float_t Mass2B          = d4s_Lambda_daughters[3][iLambdaXi];
            Int_t trackA_num        = (Int_t)d4s_Lambda_daughters[4][iLambdaXi];
            Int_t trackB_num        = (Int_t)d4s_Lambda_daughters[5][iLambdaXi];
            Float_t dcaA            = d4s_Lambda_daughters[6][iLambdaXi];
            Float_t dcaB            = d4s_Lambda_daughters[7][iLambdaXi];
            Float_t qpA             = d4s_Lambda_daughters[8][iLambdaXi];
            Float_t qpB             = d4s_Lambda_daughters[9][iLambdaXi];
            Int_t nhitsFitA         = (Int_t)d4s_Lambda_daughters[10][iLambdaXi];
            Int_t nhitsFitB         = (Int_t)d4s_Lambda_daughters[11][iLambdaXi];
            Float_t nSigmaAPi       = d4s_Lambda_daughters[12][iLambdaXi];
            Float_t nSigmaAK        = d4s_Lambda_daughters[13][iLambdaXi];
            Float_t nSigmaAP        = d4s_Lambda_daughters[14][iLambdaXi];
            Float_t nSigmaBPi       = d4s_Lambda_daughters[15][iLambdaXi];
            Float_t nSigmaBK        = d4s_Lambda_daughters[16][iLambdaXi];
            Float_t nSigmaBP        = d4s_Lambda_daughters[17][iLambdaXi];
            //cout << "iLambdaXi = " << iLambdaXi << ", trackAB = {" << trackAB.Px() << ", " << trackAB.Py() << ", " << trackAB.Pz() << endl;

            for(Int_t iLambda = 0; iLambda < Lambda_reco_counter_SE_ME; iLambda++) // Lambda loop
            {
                // Get Lambda properties
                vectorDE              = d4s_Lambda_vertex[0][iLambda];
                Stv3_perpDE           = d4s_Lambda_vertex[1][iLambda];
                helixD                = d4s_Lambda_helix[0][iLambda];
                helixE                = d4s_Lambda_helix[1][iLambda];
                ltrackD               = d4s_Lambda_TLV[0][iLambda];
                ltrackE               = d4s_Lambda_TLV[1][iLambda];
                trackDE               = d4s_Lambda_TLV[2][iLambda];
                Float_t VerdistX_DE   = d4s_Lambda_properties[0][iLambda];
                Float_t VerdistY_DE   = d4s_Lambda_properties[1][iLambda];
                Float_t dcaDE         = d4s_Lambda_properties[2][iLambda];
                Int_t   iSE_ME_Lambda = d4s_Lambda_properties[3][iLambda];
                Float_t nSigmaPD      = d4s_Lambda_daughters[0][iLambda];
                Float_t nSigmaPiE     = d4s_Lambda_daughters[1][iLambda];
                Float_t Mass2D        = d4s_Lambda_daughters[2][iLambda];
                Float_t Mass2E        = d4s_Lambda_daughters[3][iLambda];
                Int_t trackD_num      = (Int_t)d4s_Lambda_daughters[4][iLambda];
                Int_t trackE_num      = (Int_t)d4s_Lambda_daughters[5][iLambda];
                Float_t dcaD          = d4s_Lambda_daughters[6][iLambda];
                Float_t dcaE          = d4s_Lambda_daughters[7][iLambda];
                Float_t qpD           = d4s_Lambda_daughters[8][iLambda];
                Float_t qpE           = d4s_Lambda_daughters[9][iLambda];
                Int_t nhitsFitD       = (Int_t)d4s_Lambda_daughters[10][iLambda];
                Int_t nhitsFitE       = (Int_t)d4s_Lambda_daughters[11][iLambda];
                Float_t nSigmaDPi     = d4s_Lambda_daughters[12][iLambda];
                Float_t nSigmaDK      = d4s_Lambda_daughters[13][iLambda];
                Float_t nSigmaDP      = d4s_Lambda_daughters[14][iLambda];
                Float_t nSigmaEPi     = d4s_Lambda_daughters[15][iLambda];
                Float_t nSigmaEK      = d4s_Lambda_daughters[16][iLambda];
                Float_t nSigmaEP      = d4s_Lambda_daughters[17][iLambda];
                //cout << "iLambda = " << iLambda << ", trackDE = {" << trackDE.Px() << ", " << trackDE.Py() << ", " << trackDE.Pz() << endl;

                if(
                   /*
                   !(SE_ME_Flag == 0 && iLambda == iLambdaXi) && // same event: both Lambdas have to be different
                   !(SE_ME_Flag == 0 && (trackA_num == trackC_num || trackA_num == trackD_num || trackA_num == trackE_num || trackB_num == trackC_num || trackB_num == trackD_num
                                        || trackB_num == trackE_num || trackC_num == trackD_num || trackC_num == trackE_num)) &&

                   !(SE_ME_Flag == 1 && iSE_ME_Lambda == 0 && iSE_ME_Lambda == 0) && // (pA,piA) + (pA,piA)
                   !(SE_ME_Flag == 1 && iSE_ME_Lambda == 0 && iSE_ME_Lambda == 2) && // (pA,piA) + (pA,piB)
                   !(SE_ME_Flag == 1 && iSE_ME_Lambda == 0 && iSE_ME_Lambda == 3) && // (pA,piA) + (pB,piA)
                   !(SE_ME_Flag == 1 && iSE_ME_Lambda == 1 && iSE_ME_Lambda == 1) && // (pB,piB) + (pB,piB)
                   !(SE_ME_Flag == 1 && iSE_ME_Lambda == 1 && iSE_ME_Lambda == 2) && // (pB,piB) + (pA,piB)
                   !(SE_ME_Flag == 1 && iSE_ME_Lambda == 1 && iSE_ME_Lambda == 3) && // (pB,piB) + (pB,piA)

                   !(SE_ME_Flag == 1 && iSE_ME_Lambda == 0 && iSE_ME_Lambda == 1 && (trackA_num == trackC_num || trackB_num == trackC_num)) && // (pA,piA) + (pB,piB)
                   !(SE_ME_Flag == 1 && iSE_ME_Lambda == 2 && iSE_ME_Lambda == 2 && (trackA_num == trackC_num || trackA_num == trackD_num || trackB_num == trackE_num || trackC_num == trackD_num)) // (pA,piB) + (pA,piB)
                   */

                   iLambda != iLambdaXi && // only one Lambda list, even for mixed event
                   (
                    (
                     SE_ME_Flag == 0 &&
                     (
                      trackD_num != trackA_num &&
                      trackE_num != trackB_num &&
                      trackE_num != trackC_num
                     )
                    ) ||
                    (
                     SE_ME_Flag == 1 &&
                     iSE_ME_LambdaXi != iSE_ME_Lambda
                    )
                   )

                   //(trackE_num != trackC_num
                   // && trackE_num != trackB_num
                   // && trackD_num != trackA_num)
                   //|| (SE_ME_Lambda == 1)
                  )
                {
                    StThreeVectorF dir_Lambda_DE;
                    dir_Lambda_DE.set(trackDE.Px(),trackDE.Py(),trackDE.Pz());
                    dir_Lambda_DE *= 1.0/dir_Lambda_DE.mag();

                    Float_t path_helix_Xi_to_DE, path_line_Xi_to_DE, dca_Xi_to_DE;
                    fHelixAtoLinedca(dir_Lambda_DE,vectorDE,helix_ABC,path_helix_Xi_to_DE,path_line_Xi_to_DE,dca_Xi_to_DE);

                    StThreeVectorF vector_line_Xi_to_DE, vector_helix_Xi_to_DE, dir_helix_Xi_to_DE, vector_Xi_Lambda, vector_prim_to_d4s;
                    vector_line_Xi_to_DE  = vectorDE+path_line_Xi_to_DE*dir_Lambda_DE;  //
                    vector_helix_Xi_to_DE = helix_ABC.at(path_helix_Xi_to_DE);  // space vector of helix_ABC at dca to Lambda
                    dir_helix_Xi_to_DE    = helix_ABC.cat(path_helix_Xi_to_DE);  // direction vector of helix_ABC at dca to Lambda
                    dir_helix_Xi_to_DE /= dir_helix_Xi_to_DE.mag();
                    dir_helix_Xi_to_DE *= trackABC.P(); // momentum vector of Xi at dca to Lambda
                    vector_Xi_Lambda      = vector_line_Xi_to_DE + vector_helix_Xi_to_DE;
                    vector_Xi_Lambda *= 0.5; // Decay vector of Xi and second Lambda
                    vector_prim_to_d4s = vector_Xi_Lambda - vectorprim;
                    Float_t VerdistX_ABCDE = vector_prim_to_d4s.mag();

                    TLorentzVector trackDE_Lambda; // mother particle with nominal Lambda mass
                    trackDE_Lambda.SetXYZM(trackDE.Px(),trackDE.Py(),trackDE.Pz(),InvMassDE_use);

                    TLorentzVector trackABC_Xi; // mother particle with nominal Xi mass
                    trackABC_Xi.SetXYZM(dir_helix_Xi_to_DE.x(),dir_helix_Xi_to_DE.y(),dir_helix_Xi_to_DE.z(),InvMassABC_use);
                    TLorentzVector trackABC_Xi_raw; // mother particle with calculated invariant mass
                    trackABC_Xi_raw.SetXYZM(dir_helix_Xi_to_DE.x(),dir_helix_Xi_to_DE.y(),dir_helix_Xi_to_DE.z(),trackABC.M());

                    TLorentzVector trackABCDE;
                    trackABCDE = trackABC_Xi + trackDE_Lambda;

                    Float_t PtMassABCDE   = trackABCDE.Pt();
                    Float_t PMassABCDE    = trackABCDE.P();
                    Float_t InvMassABCDE  = trackABCDE.M();

                    StThreeVectorF dir_ABCDE;
                    dir_ABCDE.set(trackABCDE.Px(),trackABCDE.Py(),trackABCDE.Pz());
                    Double_t  VerdistY_ABCDE  = calculateMinimumDistanceStraightToPoint(vector_Xi_Lambda,dir_ABCDE,vectorprim);

                    //********************** Calculate the helix of the d4s ***************************
                    StPhysicalHelixD helix_d4s;
                    StThreeVectorF origin_d4s;   // decay vertex of d4s
                    origin_d4s.set(vector_Xi_Lambda.x(),vector_Xi_Lambda.y(),vector_Xi_Lambda.z());
                    Float_t h_d4s = -1.0;
                    if(ParticleC == 8 || ParticleC == 11 || ParticleC == 14) h_d4s = 1.0;
                    if(ParticleC == 9 || ParticleC == 12 || ParticleC == 15) h_d4s = -1.0;
                    Float_t phase_d4s = TMath::ATan2(-1.0*h_d4s*trackABCDE.Px(),h_d4s*trackABCDE.Py());
                    Float_t dip_d4s   = TMath::ATan2(trackABCDE.Pz(),PtMassABCDE);  // correct
                    Float_t curv_d4s  = 1.0;
                    if(PtMassABCDE != 0.0)
                    {
                        curv_d4s = curv_to_invpt_ratio/PtMassABCDE;
                    }

                    helix_d4s.setParameters(curv_d4s,dip_d4s,phase_d4s,origin_d4s,h_d4s);
                    Float_t path_d4s = -999.0;
                    Float_t dca_d4s  = -999.0;
                    fHelixAtoPointdca(vectorprim,helix_d4s,path_d4s,dca_d4s);

                    dir_ABCDE /= dir_ABCDE.mag();
                    Double_t scalarProductABCDE = dir_ABCDE.dot(vector_prim_to_d4s/vector_prim_to_d4s.mag());

                    // calculate angle between Lambda decay planes
                    Double_t angle_Lambda = Stv3_perpAB.angle(Stv3_perpDE);

                    if(d4s_reco_counter < (N_d4s-1))
                    {
                        if(d4s_reco_counter == 0)
                        {
                            // Fill event information for d4s
                            d4s_event.cleard4sList();
                            d4s_event.clearXiList();
                            d4s_event.clearLambdaList();
                            d4s_event.clearTrackList();
                            d4s_event.setx(EventVertexXA);
                            d4s_event.sety(EventVertexYA);
                            d4s_event.setz(EventVertexZA);
                            d4s_event.setid(RunIdA);
                            d4s_event.setmult(refMultA);
                            d4s_event.setn_prim(n_primaries);
                            d4s_event.setn_non_prim(n_non_primaries);
                            d4s_event.setn_tof_prim(n_tofmatch_prim);
                            d4s_event.setSE_ME_flag(SE_ME_Flag);
                            d4s_event.setZDCx(ZDCx);
                            d4s_event.setBBCx(BBCx);
                            d4s_event.setvzVpd(vzVpd);
                            d4s_event.setPsi2(phi_event_plane_d4s);
                            d4s_event.setBfield(event_SE_ME_ana[0]->bField());
                        }


                        Int_t Index_Tracks_AB[5]   = {trackA_num,trackB_num,trackC_num,trackD_num,trackE_num};
                        Int_t flag_SE_ME_Tracks[5] = {iSE_ME_LambdaXi,iSE_ME_LambdaXi,SE_ME_flagC,iSE_ME_Lambda,iSE_ME_Lambda};

                        TLorentzVector tlv_array_tracks[5]     = {ltrackA,ltrackB,ltrackC,ltrackD,ltrackE};
                        StPhysicalHelixD helix_array_tracks[5] = {helixA,helixB,helixC,helixD,helixE};
                        Float_t dca_array_tracks[5]            = {dcaA,dcaB,dcaC,dcaD,dcaE};
                        Float_t m2_array_track[5]              = {Mass2A,Mass2B,Mass2C,Mass2D,Mass2E};
                        Float_t nSigmaPi_array_tracks[5]       = {nSigmaAPi,nSigmaBPi,nSigmaCPi,nSigmaDPi,nSigmaEPi};
                        Float_t nSigmaK_array_tracks[5]        = {nSigmaAK,nSigmaBK,nSigmaCK,nSigmaDK,nSigmaEK};
                        Float_t nSigmaP_array_tracks[5]        = {nSigmaAP,nSigmaBP,nSigmaCP,nSigmaDP,nSigmaEP};
                        Float_t qp_array_tracks[5]             = {qpA,qpB,qpC,qpD,qpE};
                        Int_t   nhits_array_tracks[5]          = {nhitsFitA,nhitsFitB,nhitsFitC,nhitsFitD,nhitsFitE};
                        Int_t   trackId_array_tracks[5]        = {trackA_num,trackB_num,trackC_num,trackD_num,trackE_num};

                        Int_t Index_track_event_array[5] = {-1,-1,-1,-1,-1};
                        for(Int_t iTrack_AB = 0; iTrack_AB < 5; iTrack_AB++)
                        {
                            Int_t count_Index_Track_evAB[2] = {-1,-1}; // for event A and event B
                            if(flag_SE_ME_Tracks[iTrack_AB] == 0)
                            {
                                count_Index_Track_evAB[flag_SE_ME_Tracks[iTrack_AB]] = std::count(Index_array_track_evA,Index_array_track_evA+Track_use_counter_evAB[flag_SE_ME_Tracks[iTrack_AB]]+1,Index_Tracks_AB[iTrack_AB]);
                            }
                            if(flag_SE_ME_Tracks[iTrack_AB] == 1)
                            {
                                count_Index_Track_evAB[flag_SE_ME_Tracks[iTrack_AB]] = std::count(Index_array_track_evB,Index_array_track_evB+Track_use_counter_evAB[flag_SE_ME_Tracks[iTrack_AB]]+1,Index_Tracks_AB[iTrack_AB]);
                            }

                            if(count_Index_Track_evAB[0] <= 0 && count_Index_Track_evAB[1] <= 0)
                            {
                                d4s_Track = d4s_event.createTrack();
                                d4s_Track ->set_THelix_Track(helix_array_tracks[iTrack_AB]);
                                d4s_Track ->set_dca_to_prim(dca_array_tracks[iTrack_AB]);
                                d4s_Track ->set_m2(m2_array_track[iTrack_AB]);
                                d4s_Track ->set_Pi_nSigma(nSigmaPi_array_tracks[iTrack_AB]);
                                d4s_Track ->set_K_nSigma(nSigmaK_array_tracks[iTrack_AB]);
                                d4s_Track ->set_P_nSigma(nSigmaP_array_tracks[iTrack_AB]);
                                d4s_Track ->set_qp(qp_array_tracks[iTrack_AB]);
                                d4s_Track ->set_nhits(nhits_array_tracks[iTrack_AB]);
                                d4s_Track ->set_trackId(trackId_array_tracks[iTrack_AB]);
                                d4s_Track ->set_track_SE_ME(flag_SE_ME_Tracks[iTrack_AB]);

                                if(flag_SE_ME_Tracks[iTrack_AB] == 0)
                                {
                                    Index_array_track_evA[Track_use_counter_evAB[flag_SE_ME_Tracks[iTrack_AB]]] = Index_Tracks_AB[iTrack_AB];
                                    Index_array_total_track_evA[Track_use_counter_evAB[flag_SE_ME_Tracks[iTrack_AB]]] = Track_use_counter_total; // Index of list which is finally saved to file
                                }
                                if(flag_SE_ME_Tracks[iTrack_AB] == 1)
                                {
                                    Index_array_track_evB[Track_use_counter_evAB[flag_SE_ME_Tracks[iTrack_AB]]] = Index_Tracks_AB[iTrack_AB];
                                    Index_array_total_track_evB[Track_use_counter_evAB[flag_SE_ME_Tracks[iTrack_AB]]] = Track_use_counter_total; // Index of list which is finally saved to file
                                }
                                Track_use_counter_evAB[flag_SE_ME_Tracks[iTrack_AB]]++;
                                Track_use_counter_total++;
                            }

                            Int_t *Track_Index;
                            Int_t Track_global_Index = 0;
                            if(flag_SE_ME_Tracks[iTrack_AB] == 0)
                            {
                                Track_Index = std::find(Index_array_track_evA,Index_array_track_evA+Track_use_counter_evAB[flag_SE_ME_Tracks[iTrack_AB]]+1,Index_Tracks_AB[iTrack_AB]);
                                Track_global_Index = Index_array_total_track_evA[Track_Index - &Index_array_track_evA[0]];
                            }
                            if(flag_SE_ME_Tracks[iTrack_AB] == 1)
                            {
                                Track_Index = std::find(Index_array_track_evB,Index_array_track_evB+Track_use_counter_evAB[flag_SE_ME_Tracks[iTrack_AB]]+1,Index_Tracks_AB[iTrack_AB]);
                                Track_global_Index = Index_array_total_track_evB[Track_Index - &Index_array_track_evB[0]];
                            }
                            Index_track_event_array[iTrack_AB] = Track_global_Index;

                        }


                        //cout << "" << endl;

                        Int_t Index_Lambda_AB[2] = {iLambdaXi,iLambda};
                        Int_t count_Index_Lambda = -1;

                        Float_t VDY_Lambda_array[2]        = {VerdistY_AB,VerdistY_DE};
                        Float_t VDX_Lambda_array[2]        = {VerdistX_AB,VerdistX_DE};
                        Float_t dca_Lambda_array[2]        = {dcaAB,dcaDE};
                        TLorentzVector TLV_Lambda_array[2] = {trackAB,trackDE};
                        StThreeVectorF TV_Lambda_array[2]  = {vectorAB,vectorDE};
                        Int_t IndexPi_Lambda_array[2]      = {Index_track_event_array[1],Index_track_event_array[4]};
                        Int_t IndexP_Lambda_array[2]       = {Index_track_event_array[0],Index_track_event_array[3]};
                        StThreeVectorF TV_Lambda_plane[2]  = {Stv3_perpAB,Stv3_perpDE};

                        for(Int_t iLambda_AB = 0; iLambda_AB < 2; iLambda_AB++)
                        {
                            count_Index_Lambda = std::count(Index_array_Lambda,Index_array_Lambda+Lambda_use_counter+1,Index_Lambda_AB[iLambda_AB]);
                            if(count_Index_Lambda == 0)
                            {
                                d4s_Lambda = d4s_event.createLambda();
                                d4s_Lambda ->set_dca_to_prim(VDY_Lambda_array[iLambda_AB]);
                                d4s_Lambda ->set_vertex_dist_to_prim(VDX_Lambda_array[iLambda_AB]);
                                d4s_Lambda ->set_dca_daughters(dca_Lambda_array[iLambda_AB]);
                                d4s_Lambda ->set_TLV_Lambda(TLV_Lambda_array[iLambda_AB]);
                                d4s_Lambda ->set_TV3_Lambda(TV_Lambda_array[iLambda_AB]);
                                d4s_Lambda ->set_IndexPi(IndexPi_Lambda_array[iLambda_AB]);
                                d4s_Lambda ->set_IndexP(IndexP_Lambda_array[iLambda_AB]);
                                d4s_Lambda ->set_TV3_Lambda_plane(TV_Lambda_plane[iLambda_AB]);

                                Index_array_Lambda[Lambda_use_counter] = Index_Lambda_AB[iLambda_AB];
                                Lambda_use_counter++;
                            }
                        }

                        Int_t Lambda_Index_reco[2] = {-1,-1};
                        for(Int_t iLambda_AB = 0; iLambda_AB < 2; iLambda_AB++)
                        {
                            Int_t *Lambda_Index;
                            Lambda_Index = std::find(Index_array_Lambda,Index_array_Lambda+Lambda_use_counter+1,Index_Lambda_AB[iLambda_AB]);
                            Lambda_Index_reco[iLambda_AB] = Lambda_Index - &Index_array_Lambda[0];
                            if( Lambda_Index_reco[iLambda_AB] == Lambda_use_counter+1) // not found in list
                            {
                                cout << "ERROR: Lambda_Index not found!" << endl;
                            }
                        }

                        Int_t iXi_count = std::count(Index_Xi,Index_Xi+Xi_reco_counter+1,iXi);
                        if(iXi_count == 0)
                        {
                            d4s_Xi = d4s_event.createXi();
                            d4s_Xi ->set_dca_to_prim(VerdistY_ABC);
                            d4s_Xi ->set_vertex_dist_to_prim(VerdistX_ABC);
                            d4s_Xi ->set_dca_daughters(dcaABC);
                            d4s_Xi ->set_IndexLambdaXi(Lambda_Index_reco[0]);
                            d4s_Xi ->set_IndexXiPi(Index_track_event_array[2]);
                            d4s_Xi ->set_TV3_Xi(vectorABC);
                            d4s_Xi ->set_THelix_Xi(helix_ABC);

                            Index_Xi[Xi_reco_counter]       = iXi;
                            Index_Xi_reco[Xi_reco_counter]  = Xi_reco_counter;
                            Xi_reco_counter++;
                        }

                        Int_t Xi_Index_saved = -1;
                        Int_t *pXi_Index;
                        Int_t pXi_Index_reco;
                        pXi_Index = std::find(Index_Xi,Index_Xi+Xi_reco_counter+1,iXi);
                        pXi_Index_reco = pXi_Index - &Index_Xi[0];
                        if( pXi_Index_reco == Xi_reco_counter+1) // not found in list
                        {
                            cout << "ERROR: pXi_Index not found!" << endl;
                        }
                        else
                        {
                            Xi_Index_saved = Index_Xi_reco[pXi_Index_reco];
                        }

                        d4s_Mother = d4s_event.created4s();
                        d4s_Mother ->set_dca_to_prim(dca_d4s);
                        d4s_Mother ->set_vertex_dist_to_prim(VerdistX_ABCDE);
                        d4s_Mother ->set_dca_daughters(dca_Xi_to_DE);
                        d4s_Mother ->set_scalar_product(scalarProductABCDE);
                        d4s_Mother ->set_IndexLambda(Lambda_Index_reco[1]);
                        d4s_Mother ->set_IndexXi(Xi_Index_saved);
                        d4s_Mother ->set_TLV_d4s(trackABCDE);
                        d4s_Mother ->set_TLV_Xi(trackABC_Xi_raw);
                        d4s_Mother ->set_TV3_d4s(vector_Xi_Lambda);
                        d4s_Mother ->set_angle_Lambdas((Float_t)angle_Lambda);
                        d4s_Mother ->set_THelix_d4s(helix_d4s);

                        if(0)
                        {
                            cout << "d4s_reco_counter = " << d4s_reco_counter
                                << ", iXi(evt) = " << iXi
                                << ", Xi_idx(file) = " << Xi_Index_saved
                                << ", iLambdaXi(evt) = " << iLambdaXi
                                << ", LambdaXi_idx(file) = " << Lambda_Index_reco[0]
                                << ", iLambda(evt) = " << iLambda
                                << ", Lambda_idx(file) = " << Lambda_Index_reco[1]
                                << endl;
                            TString track_label[5] = {"A","B","C","D","E"};
                            for(Int_t y = 0; y < 5; y++)
                            {
                                cout << "track" << track_label[y] << "index(evt) = " << Index_Tracks_AB[y]
                                    << ", SE_ME = " << flag_SE_ME_Tracks[y]
                                    << ", counter = " << Track_use_counter_evAB[flag_SE_ME_Tracks[y]]
                                    << ", total = " << Track_use_counter_total
                                    << ", index(file) = " << Index_track_event_array[y]
                                    << endl;
                            }
                            //cout << "pXi_Index = " << *pXi_Index << ", at pos " << pXi_Index_reco
                            //    << ", pLambdaXi_Index = " << *pLambdaXi_Index << ", at pos = " << pLambdaXi_Index_reco << ", sav pos = " << LambdaXi_Index_saved
                            //    << ", pLambda_Index = " << *pLambda_Index << ", at pos = " << pLambda_Index_reco << ", sav pos = " << Lambda_Index_saved
                            //    << endl;

                            //cout << "d4s_reco_counter = " << d4s_reco_counter << ", scalarProductABCDE = " << scalarProductABCDE <<
                            //    ", d4s_vec = {" << dir_ABCDE.x() << ", " << dir_ABCDE.y() << ", " << dir_ABCDE.z() << "}, iXi = " <<
                            //    ", vert_vec = {" << vector_prim_to_d4s.x() << ", " << vector_prim_to_d4s.y() << ", " << vector_prim_to_d4s.z() << "}, iXi = " << iXi << ", iLambdaXi = " << iLambdaXi << ", iLambda = " << iLambda << endl;
                            //cout << "ABC_vec = {" << trackABC_Xi.Px() << ", " << trackABC_Xi.Py() << ", " << trackABC_Xi.Pz() << endl;
                            //cout << "DE_vec = {" << trackDE_Lambda.Px() << ", " << trackDE_Lambda.Py() << ", " << trackDE_Lambda.Pz() << endl;
                        }

                        d4s_reco_counter++;
                    }
                    else
                    {
                        //cout << "Warning: maximum number of d4s candidates reached: N_Xi = " << Xi_reco_counter_loop << ", N_L = " << Lambda_reco_counter_SE_ME[SE_ME_Lambda] << endl;
                    }
                }
            }
        }
        //------------------------------------------------------------------------------------------------------------



    }
    if(d4s_reco_counter > 0)
    {
        Tree_d4s  ->Fill();
    }

}
//----------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------
Int_t Lambda_mixed_analysis(Int_t PID_counter_Array_A[][N_max_PIDs], Int_t PID_counter_Array_B[][N_max_PIDs],
                  Int_t PID_Array_A[][N_max_PIDs][N_max_tracks], Int_t PID_Array_B[][N_max_PIDs][N_max_tracks],
                  StPicoAlexEvent* picoDst_A, StPicoAlexEvent* picoDst_B,
                  Int_t ParticleA, Int_t ParticleB, Int_t ParticleC, Int_t ParticleD, Int_t ParticleE, Int_t Ana_Num,Int_t SE_ME_Flag, Int_t rot_background_flag)
{

    // Event vertex information
    StThreeVectorF vectorprim,vectorprimB,vectordiff, vectornewA_lin, vectornewB_lin;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA,EventVertexXB,EventVertexYB,EventVertexZB,vertexAB_dist;
    Int_t refMultA,refMultB,RunIdA,RunIdB;

    event_A_ana       = picoDst_A;
    event_B_ana       = picoDst_B;

    EventVertexXA     = event_A_ana->primaryVertex().x();
    EventVertexYA     = event_A_ana->primaryVertex().y();
    EventVertexZA     = event_A_ana->primaryVertex().z();

    EventVertexXB     = EventVertexXA;
    EventVertexYB     = EventVertexYA;
    EventVertexZB     = EventVertexZA;
    refMultA          = event_A_ana->refMult();
    refMultB          = refMultA;
    RunIdA            = event_A_ana->runId();
    RunIdB            = RunIdA;
    Float_t ZDCx      = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx      = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd     = event_A_ana->vzVpd();

    vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);
    vectordiff.set(0.0,0.0,0.0);

    Int_t flag_tof_particleA = 0;
    Int_t flag_tof_particleB = 0;
    Int_t flag_tof_particleC = 0;
    Int_t flag_tof_particleD = 0;
    Int_t flag_tof_particleE = 0;

    Int_t Lambda_reco_counter = 0;

    if(
       SE_ME_Flag == 1  // mixed event analysis
      )
    {
        EventVertexXB  = event_B_ana->primaryVertex().x();
        EventVertexYB  = event_B_ana->primaryVertex().y();
        EventVertexZB  = event_B_ana->primaryVertex().z();
        refMultB       = event_B_ana->refMult();
        RunIdB         = event_B_ana->runId();
        vectorprimB.set(EventVertexXB,EventVertexYB,EventVertexZB);

        vectordiff     = (vectorprim - vectorprimB);
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex
    }


    Float_t radius_cut           = 2.0*2.0; // 2.0 cm radius cut for good events
    Float_t z_axis_cut           = vertex_z_cut;    // 70.0 cm
    Float_t ME_vertex_dist_cut   = 10.0; // 3.0 cm

    Int_t ME_Flag = 0;  // 0 == not accepted for mixing, 1 == accepted for mixing

    if(
       SE_ME_Flag == 1
       && vertexAB_dist < ME_vertex_dist_cut
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && (EventVertexXB*EventVertexXB + EventVertexYB*EventVertexYB) < radius_cut
       && fabs(EventVertexZB) < z_axis_cut
      )
    {
        ME_Flag = 1; // ok for mixed event analysis
    }

    //
    if(
       ((ME_Flag == 1 && SE_ME_Flag == 1) // mixed event analysis was selected
        || (SE_ME_Flag == 0)) // same event analysis was selected
       && PID_counter_Array_A[Ana_Num][ParticleA]    > 0   // p(A)
       && PID_counter_Array_A[Ana_Num][ParticleB]    > 0   // pi-(B)
       && PID_counter_Array_A[Ana_Num][ParticleC]    > 0   // pi-(C)
       && PID_counter_Array_B[Ana_Num][ParticleD]    > 0   // p(D)
       && PID_counter_Array_B[Ana_Num][ParticleE]    > 0   // pi-(E)
       && (EventVertexXA*EventVertexXA + EventVertexYA*EventVertexYA) < radius_cut
       && fabs(EventVertexZA) < z_axis_cut
       && event_A_ana   ->isMinBias()
       && event_B_ana   ->isMinBias()
      )
    {
        // Initial cuts
        Float_t dcaA_cut            = 0.2;
        Float_t dcaB_cut            = 0.9;
        Float_t dcaA_upper_cut      = 20.0;
        Float_t dcaB_upper_cut      = 20.0;
        Float_t dcaAB_cut           = 1.0;

        Float_t VerdistX_AB_cut     = 3.5;
        Float_t VerdistY_AB_cut     = 0.2;
        Float_t VerdistY_AB_upper_cut = 10.0;

        Float_t Lambda_m_upper_cut  = 1.13; // 1.12
        Float_t Lambda_m_lower_cut  = 1.1; // 1.11

        Float_t proton_m2_upper_cut = 1.5;
        Float_t proton_m2_lower_cut = 0.4;
        Float_t pion_m2_upper_cut   = 0.1;
        Float_t pion_m2_lower_cut   = -0.1;
        Float_t kaon_m2_upper_cut   = 0.4;
        Float_t kaon_m2_lower_cut   = 0.1;

        Float_t nHitsFit_cut        = 14.0;
        Float_t nHitsPoss_cut       = 0.0;
        Float_t Momentum_upper_cut  = 10.0;
        Float_t Momentum_lower_cut  = 0.15;
        Float_t nHits_ratio_cut     = 0.52;

        Float_t InvMassAB_use       = 1.115683; // Lambda mass
        Float_t InvMassABC_use      = 1.32131;  // Xi+/- mass
        Float_t InvMassAB_upper_cut = Lambda_m_upper_cut;
        Float_t InvMassAB_lower_cut = Lambda_m_lower_cut;

        Float_t InvMassDE_use       = InvMassAB_use;


        StPhysicalHelixD helixA, helixB, helixC, helixD, helixE, helix_ABC;
        StThreeVectorF vectorA, vectorB, vectoratsA, vectoratsB, vectorAB, vectorDE, vectorABC, vectorAB_est, vectorprimAB, vectornewA, vectornewB, vectornewD, vectornewE, dirY_lin;
        StThreeVectorF testA, testB, testAB, vectorABtoPrim, baseY, dirY, dirY_DE;
        TLorentzVector ltrackA, ltrackB, ltrackC, ltrackD, ltrackE, ltrackC_pi, ltrackB2, ltrackA_lin, ltrackB_lin, ltrackB2_pi, ltrackA_pip;
        TLorentzVector trackAB, trackDE, trackABC, trackAB_test;


        if(0)
        {
            cout << "d4s event accepted, A = " << PID_counter_Array_A[Ana_Num][ParticleA] <<
                ", B = " << PID_counter_Array_A[Ana_Num][ParticleB] <<
                ", C = " << PID_counter_Array_A[Ana_Num][ParticleC] <<
                ", D = " << PID_counter_Array_B[Ana_Num][ParticleD] <<
                ", E = " << PID_counter_Array_B[Ana_Num][ParticleE] <<
                ", erefMult_bin = " << erefMult_bin <<
                endl;
        }


        Float_t mass_cut_use_upper    = 0.0;
        Float_t mass_cut_use_lower    = 0.0;

        //------------------------------------------------------------------------------------------------------------
        // Lambda/K0S reconstruction
        StThreeVectorF vectorprim_add;
        vectorprim_add.set(0.0,0.0,0.0);
        Float_t B_Field = event_A_ana->bField();

        Int_t N_PID_counter_Array_protons = PID_counter_Array_A[Ana_Num][ParticleA];
        for(Int_t i = 0; i < N_PID_counter_Array_protons; i++) // p candidates
        {
            Int_t trackA_num;
            StPicoAlexTrack trackA;
            Int_t ParticleAD_use = ParticleA;
            trackA_num = PID_Array_A[Ana_Num][ParticleA][i];
            trackA     = *event_A_ana->track( trackA_num );
            ParticleAD_use = ParticleA;

            Float_t MomentumA   = trackA.gMom().mag();
            Float_t dcaA        = trackA.dca();   // distance of closest approach to primary vertex
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t BetaA       = trackA.btofBeta();  // Velocity after time-of-flight reconstruction
            Float_t nSigmaPA    = trackA.nSigmaProton();
            if(ParticleAD_use == 8  || ParticleAD_use == 9)
            {
                nSigmaPA = trackA.nSigmaPion();
                mass_cut_use_upper = pion_m2_upper_cut;
                mass_cut_use_lower = pion_m2_lower_cut;
            }
            if(ParticleAD_use == 11 || ParticleAD_use == 12)
            {
                nSigmaPA = trackA.nSigmaKaon();
                mass_cut_use_upper = kaon_m2_upper_cut;
                mass_cut_use_lower = kaon_m2_lower_cut;
            }
            if(ParticleAD_use == 14 || ParticleAD_use == 15)
            {
                nSigmaPA = trackA.nSigmaProton();
                mass_cut_use_upper = proton_m2_upper_cut;
                mass_cut_use_lower = proton_m2_lower_cut;
            }
            Float_t PolarityA   = trackA.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
            Float_t Mass2A      = -100.0;

            // calculate mass2
            if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && BetaA != 0)
            {
                flag_tof_particleA = 1;
                Mass2A = MomentumA*MomentumA*(1.0/(BetaA*BetaA) - 1.0);
            }
            else
            {
                flag_tof_particleA = 0;
            }
            Float_t MassA       = SquareRoot(Mass2A);

            if(
               dcaA          > dcaA_cut
               && dcaA       < dcaA_upper_cut
               && nHitsFitA  > nHitsFit_cut
               && nHitsPossA > nHitsPoss_cut
               && MomentumA  > Momentum_lower_cut
               && MomentumA  < Momentum_upper_cut
               && (nHitsFitA/nHitsPossA) > nHits_ratio_cut
               && (flag_tof_particleA == 0 || (flag_tof_particleA == 1 && Mass2A > mass_cut_use_lower && Mass2A < mass_cut_use_upper)) // proton mass
              )
            {
                helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin()+vectorprim_add,B_Field*MAGFIELDFACTOR, trackA.charge());

                //cout << "B_Field*MAGFIELDFACTOR = " << B_Field*MAGFIELDFACTOR << endl;

                Int_t N_PID_counter_Array_pions = PID_counter_Array_A[Ana_Num][ParticleB];
                for(Int_t j = 0; j < N_PID_counter_Array_pions; j++)  // pi- candidates
                {
                    Int_t trackB_num;
                    StPicoAlexTrack trackB;
                    Int_t ParticleBE_use = ParticleB;
                    trackB_num = PID_Array_A[Ana_Num][ParticleB][j];
                    trackB     = *event_A_ana->track( trackB_num );
                    ParticleBE_use = ParticleB;

                    Float_t MomentumB   = trackB.gMom().mag();
                    Float_t dcaB        = trackB.dca();   // distance of closest approach to primary vertex
                    Float_t nHitsPossB  = trackB.nHitsMax();
                    Float_t nHitsFitB   = trackB.nHitsFit();
                    Float_t BetaB       = trackB.btofBeta();  // Velocity after time-of-flight reconstruction
                    Float_t nSigmaPiB   = trackB.nSigmaPion();
                    if(ParticleBE_use == 8  || ParticleBE_use == 9)
                    {
                        nSigmaPiB = trackB.nSigmaPion();
                        mass_cut_use_upper = pion_m2_upper_cut;
                        mass_cut_use_lower = pion_m2_lower_cut;
                    }
                    if(ParticleBE_use == 11 || ParticleBE_use == 12)
                    {
                        nSigmaPiB = trackB.nSigmaKaon();
                        mass_cut_use_upper = kaon_m2_upper_cut;
                        mass_cut_use_lower = kaon_m2_lower_cut;
                    }
                    if(ParticleBE_use == 14 || ParticleBE_use == 15)
                    {
                        nSigmaPiB = trackB.nSigmaProton();
                        mass_cut_use_upper = proton_m2_upper_cut;
                        mass_cut_use_lower = proton_m2_lower_cut;
                    }

                    Float_t PolarityB   = trackB.charge(); // > 0 . positive charged particle, < 0 -> negative charged particle
                    Float_t Mass2B      = -100.0;

                    // calculate mass2
                    if(trackB.btofMatchFlag() > 0 && trackB.btof() != 0 && BetaB != 0)
                    {
                        flag_tof_particleB = 1;
                        Mass2B = MomentumB*MomentumB*(1.0/(BetaB*BetaB) - 1.0);
                    }
                    else
                    {
                        flag_tof_particleB = 0;
                    }
                    Float_t MassB       = SquareRoot(Mass2B);

                    if(
                       trackA_num != trackB_num // Prevent that a track is used twice
                       && dcaB       > dcaB_cut
                       && dcaB       < dcaB_upper_cut
                       && nHitsFitB  > nHitsFit_cut
                       && nHitsPossB > nHitsPoss_cut
                       && MomentumB  > Momentum_lower_cut
                       && MomentumB  < Momentum_upper_cut
                       && (nHitsFitB/nHitsPossB) > nHits_ratio_cut
                       && (flag_tof_particleB == 0 || (flag_tof_particleB == 1 && Mass2B > mass_cut_use_lower && Mass2B < mass_cut_use_upper)) // pion mass
                      )
                    {
                        helixB = StPhysicalHelixD(trackB.gMom(),trackB.origin()+vectorprim_add,B_Field*MAGFIELDFACTOR, trackB.charge());

                        Float_t VerdistX, dcaAB;
                        Int_t decay_prop_AB = calc_Decay_Properties(helixA,helixB,mass_array[ParticleA],mass_array[ParticleB],MomentumA,MomentumB,
                                                                    vectorprim,VerdistX_AB_cut,dcaAB_cut,vectorAB,ltrackA,ltrackB,trackAB,VerdistX,dcaAB);


                        if( decay_prop_AB == 1 )
                        {
                            Double_t InvMassAB          = trackAB.M(); // invariant mass of mother particle
                            Float_t  MomentumAB         = trackAB.P(); // momentum of mother particle
                            Float_t  BetaAB = TMath::Sqrt(1./(1+(InvMassAB/MomentumAB)*(InvMassAB/MomentumAB)));

                            dirY.set(trackAB.Px(),trackAB.Py(),trackAB.Pz());
                            dirY *= 1.0/dirY.mag();

                            baseY = vectorAB;
                            Double_t  VerdistY  = calculateMinimumDistanceStraightToPoint(vectorAB,dirY,vectorprim);

                            // Apply Lambda cuts
                            if(
                               InvMassAB         > InvMassAB_lower_cut
                               && InvMassAB      < InvMassAB_upper_cut
                               && dcaAB          < dcaAB_cut
                               && VerdistX       > VerdistX_AB_cut
                               && VerdistY       < VerdistY_AB_upper_cut
                              )
                            {

                                if(Lambda_reco_counter == 0)
                                {
                                    // Fill event information for d4s
                                    Lambda_event.clearLambdaList();
                                    Lambda_event.setx(EventVertexXA);
                                    Lambda_event.sety(EventVertexYA);
                                    Lambda_event.setz(EventVertexZA);
                                    Lambda_event.setid(RunIdA);
                                    Lambda_event.setmult(refMultA);
                                    Lambda_event.setn_prim(n_primaries);
                                    Lambda_event.setn_non_prim(n_non_primaries);
                                    Lambda_event.setn_tof_prim(n_tofmatch_prim);
                                    Lambda_event.setSE_ME_flag(SE_ME_Flag);
                                    Lambda_event.setZDCx(ZDCx);
                                    Lambda_event.setBBCx(BBCx);
                                    Lambda_event.setvzVpd(vzVpd);
                                }

                                Lambda_Lambda = Lambda_event.createLambda();
                                Lambda_Lambda ->set_P_dca(dcaA);
                                Lambda_Lambda ->set_P_m2(Mass2A);
                                Lambda_Lambda ->set_P_nSigma(nSigmaPA);
                                Lambda_Lambda ->set_P_qp(PolarityA*MomentumA);
                                Lambda_Lambda ->set_P_Id(trackA_num);
                                Lambda_Lambda ->set_Pi_dca(dcaB);
                                Lambda_Lambda ->set_Pi_m2(Mass2B);
                                Lambda_Lambda ->set_Pi_nSigma(nSigmaPiB);
                                Lambda_Lambda ->set_Pi_qp(PolarityB*MomentumB);
                                Lambda_Lambda ->set_Pi_Id(trackB_num);
                                Lambda_Lambda ->set_dca_to_prim(VerdistY);
                                Lambda_Lambda ->set_vertex_dist_to_prim(VerdistX);
                                Lambda_Lambda ->set_dca_daughters(dcaAB);
                                Lambda_Lambda ->set_TLV_Lambda(trackAB);
                                Lambda_Lambda ->set_TV3_Lambda(vectorAB);
                                Lambda_Lambda ->set_TV3_proton_gMom(trackA.gMom());
                                Lambda_Lambda ->set_TV3_proton_origin(trackA.origin());
                                Lambda_Lambda ->set_TV3_pion_gMom(trackB.gMom());
                                Lambda_Lambda ->set_TV3_pion_origin(trackB.origin());


                                //cout << "iLambdaXi = " << Lambda_reco_counter_SE_ME[iSE_ME] << ", trackAB = {" << trackAB.Px() << ", " << trackAB.Py() << ", " << trackAB.Pz() << endl;
                                //trackAB_test         = *d4s_Lambda_TLV[0][Lambda_reco_counter_SE_ME[iSE_ME]];
                                //cout << "iLambdaXi = " << Lambda_reco_counter_SE_ME[iSE_ME] << ", trackAB_test = {" << trackAB_test.Px() << ", " << trackAB_test.Py() << ", " << trackAB_test.Pz() << endl;
                                //h_d4s_InvMassAB->Fill(trackAB.M());
                                Lambda_reco_counter++;
                            }
                        }
                    }
                }
            }
        }

        //cout << "N_Lambdas_SE = " << Lambda_reco_counter_SE_ME[0] << endl;
        // Lambda reconstruction done
        //------------------------------------------------------------------------------------------------------------
    }
    if(Lambda_reco_counter > 0)
    {
        Tree_Lambda  ->Fill();
    }
}
//----------------------------------------------------------------------------------------------------------------------------------------------


/*
//----------------------------------------------------------------------------------------------------------------------------------------------
Int_t Jet_Analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                     Int_t ParticleA, Int_t Ana_Num, Int_t ran_tracks)
{
    if(0)
    {
        cout << "" << endl;
        cout << "-------------------------------------------------" << endl;
        cout << "New event with " << PID_counter_Array[Ana_Num][ParticleA] << " tracks" << endl;
    }
    // Event vertex information
    event_A_ana            = picoDst_A;
    StThreeVectorF pVertex = event_A_ana   ->primaryVertex();
    Float_t EventVertexX   = pVertex.x();
    Float_t EventVertexY   = pVertex.y();
    Float_t EventVertexZ   = pVertex.z();
    Int_t   refMult        = event_A_ana->refMult();
    Int_t   RunId          = event_A_ana->runId();
    Float_t ZDCx           = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx           = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd          = event_A_ana->vzVpd();
    Int_t   eventId        = event_A_ana->eventId();

    //cout << "refMult = " << refMult << ", erefMult_bin = " << erefMult_bin << endl;

    StThreeVectorF vector_prim, vectoratsA, dirA;
    StPhysicalHelixD helixA;

    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    TRandom ran_gen;
    ran_gen.SetSeed(0);

    // choose a jet definition
    Double_t jet_R = 0.4;
    JetDefinition jet_def(antikt_algorithm, jet_R);

    if(
       PID_counter_Array[Ana_Num][ParticleA] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && fabs(EventVertexZ) < vertex_z_cut
       && event_A_ana->isMinBias()
       && erefMult_bin >= 8 // 0%-5% central
      )
    {
        //----------------------------------------------------------------------------------------------------
        // Loop over all particles
        vector<PseudoJet> particles;
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            //cout << "i = " << i << endl;
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];

            //StPicoAlexTrack trackA = *event->track( trackA_num );
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            StThreeVectorF vectorA = trackA.origin();

            // Requires: momentum, origin, signed Magnetic Field
            //           and Charge of particle (+/- 1)
            StThreeVectorF vec_glob_mom = trackA.gMom();
            if(ran_tracks == 1)
            {
                Double_t ran_angle = ran_gen.Rndm()*TMath::Pi()*2.0;
                vec_glob_mom.rotateZ(ran_angle);
            }
            helixA = StPhysicalHelixD(vec_glob_mom,trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
            //helixA.setParameters(trackA.getcurvature(),trackA.getdipAngle(),trackA.getphase(),vectorA,trackA.geth());

            Float_t dcaA         = trackA.dca();
            Float_t nHitsPossA   = trackA.nHitsMax();
            Float_t nHitsFitA    = trackA.nHitsFit();
            Float_t MomentumA    = trackA.gMom().mag();
            Float_t PolarityA    = trackA.charge(); // > 0 -> positive charged particle, < 0 -> negative charged particle
            //Float_t BetaA       = trackA.btofBeta();
            Float_t nSigmaPion   = trackA.nSigmaPion();
            Float_t nSigmaKaon   = trackA.nSigmaKaon();
            Float_t nSigmaProton = trackA.nSigmaProton();

            Int_t charge = 0;
            if(PolarityA < 0) charge = 1;

            if(
               nHitsFitA     > nHitsFitA_EP_cut  // 14
               && nHitsPossA > nHitsPossA_EP_cut
               && MomentumA  > MomentumA_EP_low_cut
               && MomentumA  < MomentumA_EP_high_cut
               && ((Float_t)nHitsFitA)/((Float_t)nHitsPossA) > nHits_ratio_EP_cut
               && dcaA < 2.0
              )
            {
                Float_t pathA = -999.0;
                Float_t dcaAB = -999.0;

                fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);

                vectoratsA  = helixA.cat(pathA);
                Float_t eta = vectoratsA.pseudoRapidity();
                dirA = helixA.cat(pathA); // direction vector
                Double_t phiA = dirA.phi(); // phiA has values from -pi..pi

                Float_t p_x = MomentumA*dirA.x();
                Float_t p_y = MomentumA*dirA.y();
                Float_t p_z = MomentumA*dirA.z();
                Float_t p_t = sqrt(p_x*p_x + p_y*p_y);

                // TO DO: Add tof info for PID
                Float_t mass2 = mass_array[8]; // pion mass as default value
                if(fabs(nSigmaProton) < 3.0) mass2 = mass_array[14]; // proton
                if(fabs(nSigmaKaon)   < 3.0) mass2 = mass_array[11]; // kaon
                if(fabs(nSigmaPion)   < 3.0) mass2 = mass_array[8];  // pion
                mass2 *= mass2; // squared mass

                Float_t E = TMath::Sqrt(p_x*p_x + p_y*p_y + p_z*p_z + mass2);
                particles.push_back( PseudoJet(p_x,p_y,p_z,E) );
            }
        }

        // jet area definition
        Double_t ghost_maxrap = 1.0;
        GhostedAreaSpec area_spec(ghost_maxrap);
        AreaDefinition area_def(active_area, area_spec);

        //AreaDefinition area_def(active_area, GhostedAreaSpec(ghost_maxrap));
        //ClusterSequenceArea clust_seq(particles, jet_def, area_def);
        ClusterSequenceArea clust_seq_hard(particles, jet_def, area_def);

        // run the clustering, extract the jets
        //ClusterSequence clust_seq(particles, jet_def);
        double ptmin = 0.2;
        vector<PseudoJet> jets = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
        //vector<PseudoJet> jets = sorted_by_pt(clust_seq.inclusive_jets());

        // print out some info
        //cout << "Clustered with " << jet_def.description() << endl;



        // background estimation
        //cout << "Define JetDefinition" << endl;
        JetDefinition jet_def_bkgd(kt_algorithm, 0.4);
        //cout << "Define AreaDefinition" << endl;
        AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap));
        //cout << "Define selector" << endl;
        Selector selector = SelectorAbsRapMax(1.0) * (!SelectorNHardest(2));


        //cout << "Define JetMedianBackgroundEstimator" << endl;
        JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
        //cout << "Define Subtractor" << endl;
        Subtractor subtractor(&bkgd_estimator);
        //cout << "Define bkgd_estimator" << endl;
        bkgd_estimator.set_particles(particles);

        //cout << "Calculate jet_rho and jet_sigma" << endl;
        Double_t jet_rho   = bkgd_estimator.rho();
        //cout << "jet_sigma" << endl;
        Double_t jet_sigma = bkgd_estimator.sigma();

        //cout << "jet_rho = " << jet_rho << ", jet_sigma = " << jet_sigma << endl;

        //cout << "Jets above " << ptmin << " GeV in jets (" << jets.size() << " particles)" << endl;
        //cout << "---------------------------------------\n";
        //printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt", "area");
        //for (unsigned int i = 0; i < jets.size(); i++) {
        //    printf("%5u %15.8f %15.8f %15.8f %15.8f\n", i,
        //           jets[i].rap(), jets[i].phi(), jets[i].perp(),
        //           jets[i].area());
        //}
        //cout << endl;

        h_jet_per_event[0] ->Fill(jets.size());

        // get the subtracted jets
        if(jet_rho > 0.0)
        {
            h_jet_per_event[1] ->Fill(jets.size());
            for(Int_t i = 0; i < jets.size(); i++)
            {
                Float_t jet_pt   = jets[i].perp();
                Float_t jet_area = jets[i].area();
                Float_t jet_pt_sub = jets[i].perp() - jet_rho*jet_area;
                h_jet_pt_sub[0]   ->Fill(jet_pt_sub);
                h_jet_pt[0]       ->Fill(jet_pt);
                h_jet_area[0]     ->Fill(jet_area);
                if(jet_area > 0.15)
                {
                    h_jet_pt_sub[1]   ->Fill(jet_pt_sub);
                    h_jet_pt[1]       ->Fill(jet_pt);
                    h_jet_area[1]     ->Fill(jet_area);
                }
                if(jet_area > 0.3)
                {
                    h_jet_pt_sub[2]   ->Fill(jet_pt_sub);
                    h_jet_pt[2]       ->Fill(jet_pt);
                    h_jet_area[2]     ->Fill(jet_area);
                }
                if(jet_area > 0.4)
                {
                    h_jet_pt_sub[3]   ->Fill(jet_pt_sub);
                    h_jet_pt[3]       ->Fill(jet_pt);
                    h_jet_area[3]     ->Fill(jet_area);
                }
                if(jet_area > 0.5)
                {
                    h_jet_pt_sub[4]   ->Fill(jet_pt_sub);
                    h_jet_pt[4]       ->Fill(jet_pt);
                    h_jet_area[4]     ->Fill(jet_area);
                }
            }
        }

        // print the jets
        if(0)
        {
            cout << " pt y phi" << endl;
            for (unsigned i = 0; i < jets.size(); i++)
            {
                cout << "jet " << i << ": "<< jets[i].perp() << " "
                    << jets[i].rap() << " " << jets[i].phi() << endl;
                vector<PseudoJet> constituents = jets[i].constituents();
                for (unsigned j = 0; j < constituents.size(); j++)
                {
                    cout << " constituent " << j << " pt: " << constituents[j].perp() << endl;
                }
            }
        }

    }
}
//----------------------------------------------------------------------------------------------------------------------------------------------
*/


//----------------------------------------------------------------------------------------------------------------------------------------------
Int_t FillJet_Analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                     Int_t ParticleA, Int_t Ana_Num, Int_t SE_ME_Flag)
{
    if(0)
    {
        cout << "" << endl;
        cout << "-------------------------------------------------" << endl;
        cout << "New event with " << PID_counter_Array[Ana_Num][ParticleA] << " tracks" << endl;
    }
    // Event vertex information
    event_A_ana            = picoDst_A;
    StThreeVectorF pVertex = event_A_ana   ->primaryVertex();
    Float_t EventVertexX   = pVertex.x();
    Float_t EventVertexY   = pVertex.y();
    Float_t EventVertexZ   = pVertex.z();
    Int_t   refMult        = event_A_ana->refMult();
    Int_t   RunId          = event_A_ana->runId();
    Float_t ZDCx           = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx           = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd          = event_A_ana->vzVpd();
    Int_t   eventId        = event_A_ana->eventId();

    //cout << "refMult = " << refMult << ", erefMult_bin = " << erefMult_bin << endl;

    StThreeVectorF vector_prim, vectoratsA[2], dirA[2]; // [primary,global]
    TLorentzVector TLV_track[2]; // [primary,global]
    StPhysicalHelixD helixA[2];  // [primary,global]

    vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);

    // 0: 70%-80%
    // 1: 60%-70%
    // 2: 50%-60%
    // 3: 40%-50%
    // 4: 30%-40%
    // 5: 20%-30%
    // 6: 10%-20%
    // 7: 5%-10%
    // 8: 0%-5%

    if(
       PID_counter_Array[Ana_Num][ParticleA] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && fabs(EventVertexZ) < vertex_z_cut
       //&& event_A_ana->isMinBias()
       //&& erefMult_bin >= 7 // 0%-10% central
       //&& erefMult_bin >= 5 // 0%-30% central
       // && erefMult_bin > 1 && erefMult_bin < 7 // 10%-60% central
       //&& erefMult_bin <= 1  // 60%-80% central
       //&& erefMult_bin <= 2  // 50%-80% central
       //&& erefMult_bin > 2 && erefMult_bin < 6 // 20%-50% central
       //&& erefMult_bin > 2 && erefMult_bin < 5 // 30%-50% central
      )
    {
        // calculate event plane angle
        //Float_t phi_event_plane = calc_phi_event_plane_2nd(EP_Qx,EP_Qy);

        TVector2 QvecEtaPos, QvecEtaNeg;
        QvecEtaPos.Set(EP_Qx_eta_pos,EP_Qy_eta_pos);
        QvecEtaNeg.Set(EP_Qx_eta_neg,EP_Qy_eta_neg);

        // Fill event information for d4s
        JetTrackEvent.clearParticleList();
        JetTrackEvent.setx(EventVertexX);
        JetTrackEvent.sety(EventVertexY);
        JetTrackEvent.setz(EventVertexZ);
        JetTrackEvent.setid(RunId);
        JetTrackEvent.setmult(refMult);
        JetTrackEvent.setn_prim(n_primaries);
        JetTrackEvent.setn_non_prim(n_non_primaries);
        JetTrackEvent.setn_tof_prim(n_tofmatch_prim);
        JetTrackEvent.setSE_ME_flag(SE_ME_Flag);
        JetTrackEvent.setZDCx(ZDCx);
        JetTrackEvent.setBBCx(BBCx);
        JetTrackEvent.setvzVpd(vzVpd);
        JetTrackEvent.setQvecEtaPos(QvecEtaPos);
        JetTrackEvent.setQvecEtaNeg(QvecEtaNeg);
        JetTrackEvent.setcent9(erefMult_bin);

        //QvecEtaPos.Print();
        //QvecEtaNeg.Print();

        //----------------------------------------------------------------------------------------------------
        // Loop over all particles
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            //cout << "i = " << i << endl;
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];

           //StPicoAlexTrack trackA = *event->track( trackA_num );
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            StThreeVectorF vectorA = trackA.origin();

            // Requires: momentum, origin, signed Magnetic Field
            //           and Charge of particle (+/- 1)
            StThreeVectorF mom_vec[2];
            mom_vec[0] = trackA.pMom();
            mom_vec[1] = trackA.gMom();

            helixA[0] = StPhysicalHelixD(mom_vec[0],vector_prim,event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());      // primary
            helixA[1] = StPhysicalHelixD(mom_vec[1],trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());  // global

            Float_t MomentumA[2];
            Float_t dcaA           = trackA.dca();
            Float_t nHitsPossA     = trackA.nHitsMax();
            Float_t nHitsFitA      = trackA.nHitsFit();
            MomentumA[0]           = trackA.pMom().mag();
            MomentumA[1]           = trackA.gMom().mag();
            Float_t PolarityA      = trackA.charge(); // > 0 -> positive charged particle, < 0 -> negative charged particle
            //Float_t BetaA        = trackA.btofBeta();
            Float_t nSigmaPion     = trackA.nSigmaPion();
            Float_t nSigmaKaon     = trackA.nSigmaKaon();
            Float_t nSigmaProton   = trackA.nSigmaProton();
            Float_t Beta           = trackA.btofBeta();  // Velocity after time-of-flight reconstruction

            Int_t charge = 0;
            if(PolarityA < 0) charge = 1;

            // brandon
            if(
               nHitsFitA          > nHitsFitA_EP_cut  // 14
               //nHitsFitA        >= 20
               //&& nHitsPossA      > nHitsPossA_EP_cut
               && MomentumA[1]    > MomentumA_EP_low_cut
               && MomentumA[1]    < 500.0
               && ((Float_t)nHitsFitA)/((Float_t)nHitsPossA) > nHits_ratio_EP_cut
               && dcaA < 2.0
              )
            {
                std::cout << "\n nHistPossA = " << nHitsPossA;
                Float_t pathA[2] = {-999.0,-999.0};
                Float_t dcaAB[2] = {-999.0,-999.0};
                Float_t eta[2]   = {-999.0,-999.0};
                Float_t phi[2]   = {-999.0,-999.0};
                Float_t mass2[2] = {-100.0,-100.0};

                for(Int_t i_prim_glob = 0; i_prim_glob < 2; i_prim_glob++) // calculate track properties for primaries and globals
                {
                    fHelixAtoPointdca(vector_prim,helixA[i_prim_glob],pathA[i_prim_glob],dcaAB[i_prim_glob]);
                    dirA[i_prim_glob]        = helixA[i_prim_glob].cat(pathA[i_prim_glob]); // direction vector
                    dirA[i_prim_glob]       /= dirA[i_prim_glob].mag();
                    //eta[i_prim_glob]         = dirA[i_prim_glob].pseudoRapidity();
                    //phi[i_prim_glob]         = dirA[i_prim_glob].phi(); // phiA has values from -pi..pi

                    Float_t p_x = MomentumA[i_prim_glob]*dirA[i_prim_glob].x();
                    Float_t p_y = MomentumA[i_prim_glob]*dirA[i_prim_glob].y();
                    Float_t p_z = MomentumA[i_prim_glob]*dirA[i_prim_glob].z();
                    Float_t p_t = sqrt(p_x*p_x + p_y*p_y);

                    // calculate mass2
                    if(trackA.btofMatchFlag() > 0 && trackA.btof() != 0 && Beta != 0)
                    {
                        mass2[i_prim_glob] = MomentumA[i_prim_glob]*MomentumA[i_prim_glob]*(1.0/(Beta*Beta) - 1.0);
                    }

                    Float_t mass2_use = mass_array[8]; // pion mass as default value
                    if(fabs(nSigmaProton) < 3.0) mass2_use = mass_array[14]; // proton
                    if(fabs(nSigmaKaon)   < 3.0) mass2_use = mass_array[11]; // kaon
                    if(fabs(nSigmaPion)   < 3.0) mass2_use = mass_array[8];  // pion
                    mass2_use *= mass2_use; // squared mass

                    Float_t E = TMath::Sqrt(p_x*p_x + p_y*p_y + p_z*p_z + mass2_use);
                    TLV_track[i_prim_glob].SetPxPyPzE(p_x,p_y,p_z,E);
                }


                JetTrackParticle = JetTrackEvent.createParticle();
                JetTrackParticle ->set_dca_to_prim(dcaA);
                JetTrackParticle ->set_Particle_m2(mass2[1]);
                JetTrackParticle ->set_Particle_nSigmaPi(nSigmaPion);
                JetTrackParticle ->set_Particle_nSigmaK(nSigmaKaon);
                JetTrackParticle ->set_Particle_nSigmaP(nSigmaProton);
                JetTrackParticle ->set_Particle_qp(MomentumA[1]*PolarityA);
                JetTrackParticle ->set_Particle_hits_fit(nHitsFitA);
                JetTrackParticle ->set_TLV_Particle_prim(TLV_track[0]);
                JetTrackParticle ->set_TLV_Particle_glob(TLV_track[1]);

                //if(p_t > 10.0)
                //{
                //    cout << "pt = " << p_t << ", phi = " << TLV_track.Phi() << ", eta = " << TLV_track.Eta()
                //        << ", dca = " << dcaA << ", nHitsFitA = " << nHitsFitA
                //        << ", nHitsPossA = " << nHitsPossA << endl;
                //}
            }
        }

        Tree_JetTrackEvent  ->Fill();

    }
}
//----------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------
// Event plane function to calculate Q-vectors for different harmonics + full TPC EP + eta-sub (gap) EP with different eta gap values + BBC EP
Int_t Event_Analysis(Int_t PID_counter_Array[][N_max_PIDs], Int_t PID_Array[][N_max_PIDs][N_max_tracks], StPicoAlexEvent* picoDst_A,
                     Int_t ParticleA, Int_t Ana_Num,Int_t run_events,Int_t event_number)
{
    //cout << "New event" << endl;
    // Event vertex information
    event_A_ana            = picoDst_A;
    StThreeVectorF pVertex = event_A_ana   ->primaryVertex();
    Float_t EventVertexX   = pVertex.x();
    Float_t EventVertexY   = pVertex.y();
    Float_t EventVertexZ   = pVertex.z();
    Int_t   refMult        = event_A_ana->refMult();
    Int_t   RunId          = event_A_ana->runId();
    Float_t ZDCx           = event_A_ana->ZDCx(); // ZDC coincidence rate -> luminosity
    Float_t BBCx           = event_A_ana->BBCx(); // BBC coincidence rate -> luminosity
    Float_t vzVpd          = event_A_ana->vzVpd();
    Int_t   eventId        = event_A_ana->eventId();

    StThreeVectorF vector_prim, vector_prim_new, vectoratsA, dirA;
    StPhysicalHelixD helixA;

    Float_t mean_pt[2]    = {0.0,0.0}; // [+,-] charge
    Int_t   N_charge[2]   = {0,0};     // [+,-] charge
    Float_t mean_ptB[2]   = {0.0,0.0}; // [+,-] charge
    Int_t   N_chargeB[2]  = {0,0};     // [+,-] charge
    Float_t charge_asymm  = 0.0;

    if(
       PID_counter_Array[Ana_Num][ParticleA] > 0
       && (EventVertexX*EventVertexX + EventVertexY*EventVertexY) < Event_radius_cut*Event_radius_cut
       && fabs(EventVertexZ) < vertex_z_cut
       && event_A_ana->isMinBias()
      )
    {
        vector_prim.set(EventVertexX,EventVertexY,EventVertexZ);
        //----------------------------------------------------------------------------------------------------
        // Loop over all particles
        for(Int_t i = 0; i < PID_counter_Array[Ana_Num][ParticleA]; i++) //
        {
            //cout << "i = " << i << endl;
            // Get the tracks and calculate the direction and base vectors
            Int_t trackA_num = PID_Array[Ana_Num][ParticleA][i];

            //StPicoAlexTrack trackA = *event->track( trackA_num );
            StPicoAlexTrack trackA = *picoDst_A->track( trackA_num );
            StThreeVectorF vectorA = trackA.origin();

            // Requires: momentum, origin, signed Magnetic Field
            //           and Charge of particle (+/- 1)
            helixA = StPhysicalHelixD(trackA.gMom(),trackA.origin(),event_A_ana->bField()*MAGFIELDFACTOR, trackA.charge());
            //helixA.setParameters(trackA.getcurvature(),trackA.getdipAngle(),trackA.getphase(),vectorA,trackA.geth());

            Float_t dcaA        = trackA.dca();
            Float_t nHitsPossA  = trackA.nHitsMax();
            Float_t nHitsFitA   = trackA.nHitsFit();
            Float_t MomentumA   = trackA.gMom().mag();
            Float_t PolarityA   = trackA.charge(); // > 0 -> positive charged particle, < 0 -> negative charged particle
            //Float_t BetaA       = trackA.btofBeta();
            Float_t nSigmaPion  = trackA.nSigmaPion();

            Int_t charge = 0;
            if(PolarityA < 0) charge = 1;

            if(
               nHitsFitA     > nHitsFitA_EP_cut  // 14
               && nHitsPossA > nHitsPossA_EP_cut
               && MomentumA  > MomentumA_EP_low_cut
               && MomentumA  < MomentumA_EP_high_cut
               && ((Float_t)nHitsFitA)/((Float_t)nHitsPossA)    > nHits_ratio_EP_cut
               && dcaA < 2.0
              )
            {
                Float_t pathA = -999.0;
                Float_t dcaAB = -999.0;

                fHelixAtoPointdca(vector_prim,helixA,pathA,dcaAB);

                vectoratsA  = helixA.cat(pathA);
                Float_t eta = vectoratsA.pseudoRapidity();
                dirA = helixA.cat(pathA); // direction vector
                Double_t phiA = dirA.phi(); // phiA has values from -pi..pi

                //if(fabs(eta) > 2.0) cout << "eta = " << eta << endl;

                Float_t p_x = MomentumA*dirA.x();
                Float_t p_y = MomentumA*dirA.y();
                Float_t p_z = MomentumA*dirA.z();
                Float_t p_t = sqrt(p_x*p_x + p_y*p_y);

                if(eta > -1.0 && eta < 1.0)
                {
                    N_charge[charge]++;
                    mean_pt[charge] += p_t;
                    if(p_t > 0.1 && p_t < 0.5 && nSigmaPion > -2.0 && nSigmaPion < 2.0)
                    {
                        N_chargeB[charge]++;
                        mean_ptB[charge] += p_t;
                    }
                }

            }
        }
        // End of particle loop


        if((N_charge[0] + N_charge[1]) > 0)
        {
            charge_asymm = ((Float_t)(N_charge[0] - N_charge[1]))/((Float_t)(N_charge[0] + N_charge[1]));
        }

        for(Int_t k = 0; k < 2; k++)
        {
            if(N_charge[k] > 0)
            {
                mean_pt[k] /= (Float_t)N_charge[k];
                p_mean_pt_vs_charge_asymm[k] ->Fill(charge_asymm,mean_pt[k]);
            }
            if(N_chargeB[k] > 0)
            {
                mean_ptB[k] /= (Float_t)N_chargeB[k];
            }
        }


        EventAnalysis_NTDataArray[0]    =(Float_t)EventVertexX;
        EventAnalysis_NTDataArray[1]    =(Float_t)EventVertexY;
        EventAnalysis_NTDataArray[2]    =(Float_t)EventVertexZ;
        EventAnalysis_NTDataArray[3]    =(Float_t)refMult;
        EventAnalysis_NTDataArray[4]    =(Float_t)N_charge[0];
        EventAnalysis_NTDataArray[5]    =(Float_t)N_charge[1];
        EventAnalysis_NTDataArray[6]    =(Float_t)mean_pt[0];
        EventAnalysis_NTDataArray[7]    =(Float_t)mean_pt[1];
        EventAnalysis_NTDataArray[8]    =(Float_t)N_chargeB[0];
        EventAnalysis_NTDataArray[9]    =(Float_t)N_chargeB[1];
        EventAnalysis_NTDataArray[10]   =(Float_t)mean_ptB[0];
        EventAnalysis_NTDataArray[11]   =(Float_t)mean_ptB[1];
        EventAnalysis_NTDataArray[12]   =(Float_t)charge_asymm;

        EventAnalysis_NT->Fill(EventAnalysis_NTDataArray);
        //----------------------------------------------------------------------------------------------------

        return 1;
    }
    else return 0;
}
//----------------------------------------------------------------------------------------------------------------------------------------------




#endif /* ANALYSIS_SNURF_FUNC_H */
