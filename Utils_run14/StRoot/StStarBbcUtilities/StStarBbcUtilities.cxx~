// - Implement phi angle of BBC (codes from Yadav)
// - Read, get calibration constants (recentering)

#include <fstream>

#include "TMath.h"
#include "TRandom.h"
#include "StStarBbcUtilities.h"

//using std::cout ;
using std::ifstream ;
//using std::endl ;

ClassImp(StStarBbcUtilities)

// Default constructor
StStarBbcUtilities::StStarBbcUtilities()
{

}

//StStarBbcUtilities* StStarBbcUtilities::fInstance = 0;

//______________________________________________________________________________
//StStarBbcUtilities::StStarBbcUtilities()
//  : TQVectorUtilities("Bbc", 2, kNHar, kNHar) // 0:East, 1:West
//{
//  gRandom->SetSeed(0);
//  SetTitle(0, "BBC East");
//  SetTitle(1, "BBC West");
//}

//______________________________________________________________________________
StStarBbcUtilities::~StStarBbcUtilities()
{}

//______________________________________________________________________________
//StStarBbcUtilities* StStarBbcUtilities::Instance()
//{
//  if( fInstance ) return fInstance ;

//  fInstance = new StStarBbcUtilities() ;

//  return fInstance ;
//}

//______________________________________________________________________________
UShort_t StStarBbcUtilities::GetAdcCut(const TString energy) const
{
  if ( energy.Contains("200") ) return 10 ;
  else if ( energy.Contains("11.5") ) return 40 ;
  else return 10 ;
}


//______________________________________________________________________________
Float_t StStarBbcUtilities::GetXY(const Int_t tileId, const Int_t idXY) const
{
    Float_t pos = 0.0;

    if((tileId >= 16 || tileId < 0) || !(idXY == 0 || idXY == 1))
    {
        return -1;
    }

    Float_t tile_Xpos[16] = {0.0,7.23,7.23,0.0,-7.23,-7.23,-7.23,0.0,14.46,14.46,14.46,7.23,0.0,-14.46,-14.46,-14.46};
    Float_t tile_Ypos[16] = {8.34848,4.17424,-4.17424,-8.34848,-4.17424,4.17424,12.5227,16.697,8.34848,0.0,-8.34848,-12.5227,-16.697,-8.34848,0.0,8.34848};

    if(tileId == 6 || tileId == 11) // tiles are connected to one PMT
    {
        if(gRandom->Rndm() > 0.5)
        {
            tile_Xpos[tileId] = -tile_Xpos[tileId];
        }
    }

    if(idXY == 0)
    {
        pos = tile_Xpos[tileId];
    }
    if(idXY == 1)
    {
        pos = tile_Ypos[tileId];
    }

    return pos;
}


//______________________________________________________________________________
void StStarBbcUtilities::GetRandomXY(const Int_t tileId, Float_t &x_pos_random, Float_t &y_pos_random) const
{
    Float_t x_pos = GetXY(tileId,0);
    Float_t y_pos = GetXY(tileId,1);

    const Float_t outer_radius = 9.64/2.0;
    const Float_t slope = TMath::Tan(TMath::Pi()/3.0);
    const Float_t t_val = slope*outer_radius;
    const Float_t slope_array[4] = {slope,-slope,slope,-slope};
    const Float_t t_array[4]     = {-t_val,t_val,t_val,-t_val};

    Int_t pos_out_of_hexagon = 0;

    while(pos_out_of_hexagon != 1)
    {
        x_pos_random = (gRandom->Rndm()-0.5)*outer_radius;
        y_pos_random = (gRandom->Rndm()-0.5)*outer_radius*TMath::Cos(TMath::Pi()/6.0);

        if(
           (y_pos_random > slope_array[0]*x_pos_random + t_array[0]) &&
           (y_pos_random < slope_array[1]*x_pos_random + t_array[1]) &&
           (y_pos_random < slope_array[2]*x_pos_random + t_array[2]) &&
           (y_pos_random > slope_array[3]*x_pos_random + t_array[3])
          )
        {
            pos_out_of_hexagon = 1;
        }
    }

    x_pos_random += x_pos;
    y_pos_random += y_pos;
}


//______________________________________________________________________________
// From Yadav's code
Float_t StStarBbcUtilities::GetPhi(const Int_t eastWest, const Int_t tileId) const
{
  //Float_t GetPhiInBBC(int eastWest, int bbcN) { //tileId=0 to 23
  const Float_t Pi = TMath::Pi() ;
  const Float_t phi_div=Pi/6.0; // 30 degree
  Float_t bbc_phi=phi_div;
  switch(tileId) {
  case 0: bbc_phi=3.0*phi_div;
          break;
  case 1: bbc_phi=phi_div;
          break;
  case 2: bbc_phi=-1.0*phi_div;
          break;
  case 3: bbc_phi=-3.0*phi_div;
          break;
  case 4: bbc_phi=-5.0*phi_div;
          break;
  case 5: bbc_phi=5.0*phi_div;
          break;
  case 6: bbc_phi= (gRandom->Rndm()>0.5) ? 2.0*phi_div:4.0*phi_div;
          break;
  case 7: bbc_phi=3.0*phi_div;
          break;
  case 8: bbc_phi=phi_div;
          break;
  case 9: bbc_phi=0.0;
          break;
  case 10: bbc_phi=-phi_div;
           break;
  case 11: bbc_phi=(gRandom->Rndm()>0.5) ? -2.0*phi_div:-4.0*phi_div;
           break;
  case 12: bbc_phi=-3.0*phi_div;
           break;
  case 13: bbc_phi=-5.0*phi_div;
           break;
  case 14: bbc_phi=Pi;
           break;
  case 15: bbc_phi=5.0*phi_div;
           break;
  case 16: bbc_phi=3.0*phi_div;
           break;
  case 17: bbc_phi=0.0;
           break;
  case 18: bbc_phi=-3.0*phi_div;
           break;
  case 19: bbc_phi= Pi;
           break;
  case 20: bbc_phi=3.0*phi_div;
           break;
  case 21: bbc_phi=0.0;
           break;
  case 22: bbc_phi=-3.0*phi_div;
           break;
  case 23: bbc_phi= Pi;
           break;
  }

  if(eastWest==0) {
    if(bbc_phi > -0.001) {bbc_phi = Pi-bbc_phi;}
    else {bbc_phi = -Pi-bbc_phi;}
  }
  if (bbc_phi < 0.0) { bbc_phi += 2.0*Pi;}

  return bbc_phi;
}

//______________________________________________________________________________
UInt_t StStarBbcUtilities::EastWest(const UInt_t ipmt_all) const
{
  return (ipmt_all<24) ? 0 : 1 ;
}

//______________________________________________________________________________
UInt_t StStarBbcUtilities::PmtId(const UInt_t ipmt_all) const
{
  return (ipmt_all%24) ;
}

//______________________________________________________________________________
UInt_t StStarBbcUtilities::InnerOuter(const UInt_t ipmt_each) const
{
  return (ipmt_each<16) ? 0 : 1 ;
}


