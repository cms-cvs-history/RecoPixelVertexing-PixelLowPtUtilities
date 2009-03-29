#include "RecoPixelVertexing/PixelLowPtUtilities/interface/TripletGenerator.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ThirdHitPrediction.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/TripletFilter.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/HitInfo.h"

#include "RecoTracker/TkHitPairs/interface/LayerHitMap.h"
#include "RecoTracker/TkHitPairs/interface/LayerHitMapLoop.h"

#include "RecoTracker/TkMSParametrization/interface/PixelRecoPointRZ.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#undef Debug

using namespace std;
using namespace ctfseeding;

/*****************************************************************************/
void TripletGenerator::init(const HitPairGenerator & pairs,
      const vector<SeedingLayer> & layers,
      LayerCacheType* layerCache)
{
  thePairGenerator = pairs.clone();
  theLayers        = layers;
  theLayerCache    = layerCache;

  checkMultipleScattering = ps.getParameter<bool>("checkMultipleScattering");
  nSigMultipleScattering  = ps.getParameter<double>("nSigMultipleScattering");
  checkClusterShape       = ps.getParameter<bool>("checkClusterShape"); 
  rzTolerance             = ps.getParameter<double>("rzTolerance");
  maxAngleRatio           = ps.getParameter<double>("maxAngleRatio");
  builderName             = ps.getParameter<string>("TTRHBuilder");
}

/*****************************************************************************/
void TripletGenerator::getTracker
  (const edm::EventSetup& es)
{
  if(theTracker == 0)
  {
    // Get tracker geometry
    edm::ESHandle<TrackerGeometry> tracker;
    es.get<TrackerDigiGeometryRecord>().get(tracker);

    theTracker = tracker.product();
  }

  if(theFilter == 0)
  {
    theFilter = new TripletFilter(es); 
  }
}

/*****************************************************************************/
GlobalPoint TripletGenerator::getGlobalPosition
  (const TrackingRecHit* recHit)
{
  DetId detId = recHit->geographicalId();

  return
    theTracker->idToDet(detId)->toGlobal(recHit->localPosition());
}

/*****************************************************************************/
void TripletGenerator::hitTriplets(
    const TrackingRegion& region,
    OrderedHitTriplets & result,
    const edm::Event & ev,
    const edm::EventSetup& es) 
{
  // Generate pairs
  OrderedHitPairs pairs; pairs.reserve(30000);
  thePairGenerator->hitPairs(region,pairs,ev,es);

  if (pairs.size() == 0) return;

  int size = theLayers.size(); 

  // Set aliases
  const LayerHitMap **thirdHitMap = new const LayerHitMap* [size];
  for(int il=0; il<size; il++)
    thirdHitMap[il] = &(*theLayerCache)(&theLayers[il], region, ev, es);

  // Get tracker
  getTracker(es);

  // Look at all generated pairs
  for(OrderedHitPairs::const_iterator ip = pairs.begin();
                                      ip!= pairs.end(); ip++)
  {
    // Fill rechits and points
    vector<const TrackingRecHit*> recHits(3);
    vector<GlobalPoint> points(3);

    recHits[0] = (*ip).inner()->hit();
    recHits[1] = (*ip).outer()->hit();

#ifdef Debug
    cerr << " RecHits " + HitInfo::getInfo(*recHits[0]) +
                          HitInfo::getInfo(*recHits[1]) << endl;
#endif

    for(int i=0; i<2; i++)
      points[i] = getGlobalPosition(recHits[i]);

    // Initialize helix prediction
    ThirdHitPrediction
      thePrediction(region.originRBound(), region.ptMin(),
                    points[0],points[1], es,
                    nSigMultipleScattering,maxAngleRatio,builderName);

    // Look at all layers
    for(int il=0; il<size; il++)
    {
      const SeedingLayer & layerwithhits = theLayers[il];
      const DetLayer * layer = layerwithhits.detLayer();

#ifdef Debug
      cerr << "  check layer " << layer->subDetector()
                        << " " << layer->location() << endl;
#endif

      // Get ranges for the third hit
      float phi[2],rz[2];
      thePrediction.getRanges(layer, phi,rz);

      PixelRecoRange<float> phiRange(phi[0]              , phi[1]             );
      PixelRecoRange<float>  rzRange( rz[0] - rzTolerance, rz[1] + rzTolerance);

      // Get third hit candidates from cache
      LayerHitMapLoop thirdHits = thirdHitMap[il]->loop(phiRange, rzRange);
      const SeedingHit * th;
      while( (th = thirdHits.getHit()) )
      {
        // Fill rechit and point
        recHits[2] = *th;
        points[2]  = getGlobalPosition(recHits[2]);

#ifdef Debug
        cerr << "  third hit " + HitInfo::getInfo(*recHits[2]) << endl;
#endif

        // Check if third hit is compatible with multiple scattering
        vector<GlobalVector> globalDirs;
        if(thePrediction.isCompatibleWithMultipleScattering
             (points[2], recHits, globalDirs, es) == false)
        {
#ifdef Debug
          cerr << "  not compatible: multiple scattering" << endl;
#endif
          if(checkMultipleScattering) continue;
        }

        // Convert to localDirs
/*
        vector<LocalVector> localDirs;
        vector<GlobalVector>::const_iterator globalDir = globalDirs.begin();
        for(vector<const TrackingRecHit *>::const_iterator
                                            recHit  = recHits.begin();
                                            recHit != recHits.end(); recHit++)
        {
          localDirs.push_back(theTracker->idToDet(
                             (*recHit)->geographicalId())->toLocal(*globalDir));
          globalDir++;
        }
*/

        // Check if the cluster shapes are compatible with thrusts
        if(checkClusterShape)
        {
          if(! theFilter->checkTrack(recHits,globalDirs))
          {
#ifdef Debug
            cerr << "  not compatible: cluster shape" << endl;
#endif
            continue;
          }
        }

        // All checks passed, put triplet back
        result.push_back(OrderedHitTriplet((*ip).inner(),(*ip).outer(),*th));
      }
    }
  } 
  delete [] thirdHitMap;

  return;
}


