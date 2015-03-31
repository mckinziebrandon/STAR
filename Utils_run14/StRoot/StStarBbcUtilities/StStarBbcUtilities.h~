
#ifndef __StStarBbcUtilities__h_
#define __StStarBbcUtilities__h_

//#include "TQVectorUtilities.h"

//______________________________________________________________________________
//class StStarBbcUtilities : public TQVectorUtilities
class StStarBbcUtilities
{
public:

    StStarBbcUtilities();
    //static StStarBbcUtilities* Instance(); // return singleton pointer
    virtual ~StStarBbcUtilities();

    /// Get low ADC cut off for event planes
    UShort_t GetAdcCut(const TString energy) const ;

    /// Get phi position of BBC (code from Yadav)
    Float_t GetPhi(const Int_t eastWest, const Int_t tileId) const ;

    /// Get XY position of BBC
    Float_t GetXY(const Int_t tileId, const Int_t idXY) const ;

    /// Get random XY position of BBC
    void GetRandomXY(const Int_t tileId, Float_t &x_pos_random, Float_t &y_pos_random) const ;

    /// East or west
    UInt_t EastWest(const UInt_t ipmt_all) const ; /// East:0-23, West:24-47

    /// PMT id for each arm
    UInt_t PmtId(const UInt_t ipmt_all) const ; /// PMT id for each arm: 0-23

    /// Innter of outer
    UInt_t InnerOuter(const UInt_t ipmt_each) const ; /// Innter:0-15, Outer:16-23

    /// Read event plane resolution
    //    void ReadEventPlaneResolution(const UInt_t dayOrRun,
    //        const Char_t* targetDir="./table") ; // File name should be like "resolution_{det}_%d.txt"

private:
    //static StStarBbcUtilities* fInstance ; // singleton

    ClassDef(StStarBbcUtilities, 0)
};

#endif

