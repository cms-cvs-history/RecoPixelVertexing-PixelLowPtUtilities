#ifndef _LowPtThirdHitRZPrediction_h_
#define _LowPtThirdHitRZPrediction_h_

/** predicts a range in r-z for the position of third hit.
 *  the predicted reange is defined by stright line extrapolation/interpolation
 *  from hit pair and the margin defined by hit errors and multiple scattering
 */

#include "RecoTracker/TkMSParametrization/interface/PixelRecoRange.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoLineRZ.h"
#include "RecoTracker/TkTrackingRegions/interface/TkTrackingRegionsMargin.h"

#include "DataFormats/GeometryVector/interface/GlobalTag.h"
#include "DataFormats/GeometryVector/interface/Vector2DBase.h"
typedef Vector2DBase<float,GlobalTag> Global2DVector;

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

class DetLayer;
class OrderedHitPair;

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <vector>
#include <utility>
using namespace std;

class LowPtThirdHitRZPrediction {
public:
  typedef PixelRecoRange<float> Range;
  typedef TkTrackingRegionsMargin<float> Margin;

  LowPtThirdHitRZPrediction
    (float originRBound, float ptMin, GlobalPoint inner, GlobalPoint outer, const edm::EventSetup& es);
  ~LowPtThirdHitRZPrediction();

  void getRanges(const DetLayer *layer, float phi[],float rz[]);
  void getRanges(float rORz           , float phi[],float rz[]);

  bool isCompatibleWithMultipleScattering
    (GlobalPoint g3, vector<GlobalVector>& localDirs, const edm::EventSetup& es);

private:
  void initLayer(const DetLayer *layer);

  void printOut(char *text);

  void invertCircle(Global2DVector& c,float& r);
  void invertPoint (Global2DVector& p);

  pair<float,float> findMinimalCircles (float r);
  pair<float,float> findTouchingCircles(float r);

  pair<float,float> findArcIntersection
    (pair<float,float> a, pair<float,float> b, bool& keep);

  void fitParabola
    (const float x[3], const float y[3], float par[3]);
  void findRectangle
    (const float x[3], const float y[3], const float par[3],
     float phi[2],float z[2]);

  float areaParallelogram(const Global2DVector& a , const Global2DVector& b);
  float angleRatio       (const Global2DVector& p3, const Global2DVector& c);

  void spinCloser(float phi[3]);

  void calculateRangesBarrel (float  r3, float phi[2],float z[2], bool keep);
  void calculateRangesForward(float  z3, float phi[2],float r[2], bool keep);
  void calculateRanges       (float rz3, float phi[2],float rz[2]);

  bool theBarrel, theForward;
  Range theDetRange;
  Margin theTolerance;
  PixelRecoLineRZ theLine;

  const DetLayer * theLayer;

  // Data
  float Bz, r0,rm;
  GlobalPoint    g1,g2;
  Global2DVector p1,p2,dif;
  pair<float,float> arc_0m;

  bool keep;
};
#endif
