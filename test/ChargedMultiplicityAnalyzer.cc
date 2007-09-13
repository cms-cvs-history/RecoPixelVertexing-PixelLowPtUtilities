#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterData.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShape.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"

#include <utility>
#include <vector>
#include <fstream>

using namespace std;
using namespace reco;

/*****************************************************************************/
class ChargedMultiplicityAnalyzer : public edm::EDAnalyzer
{
 public:
   explicit ChargedMultiplicityAnalyzer(const edm::ParameterSet& pset);
   ~ChargedMultiplicityAnalyzer();
   virtual void beginJob(const edm::EventSetup& es);
   virtual void analyze(const edm::Event& ev, const edm::EventSetup& es);
   virtual void endJob() { }

 private:
   bool isAtEdge(const RectangularPixelTopology * topology,
                 const SiPixelRecHit & recHit);
   int getNumberOfPixelBarrelHits(const TrackingParticle& simTrack);
   void checkSimTracks(const VertexCollection * vertices,
                       edm::Handle<TrackingParticleCollection>& simCollection);

   TFile* file;
   TNtuple* trackSim;
   TNtuple* multi;

   const TrackerGeometry* theTracker;
   const SiPixelRecHitCollection* theHits;
};

/*****************************************************************************/
void ChargedMultiplicityAnalyzer::beginJob(const edm::EventSetup& es)
{
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker); theTracker =
tracker.product(); }

/*****************************************************************************/
ChargedMultiplicityAnalyzer::ChargedMultiplicityAnalyzer(const edm::ParameterSet& pset)
{ 
  file = new TFile("multi.root","RECREATE");
  file->cd();
  multi    = new TNtuple("multi","multi", "eta:eloss:type");
  trackSim = new TNtuple("trackSim","trackSim", "nvtx:q:id:eta:pt:rho:npb");
}

/*****************************************************************************/
ChargedMultiplicityAnalyzer::~ChargedMultiplicityAnalyzer()
{ 
  file->cd();
  multi->Write();
  trackSim->Write();
  file->Close();
}

/*****************************************************************************/
int ChargedMultiplicityAnalyzer::getNumberOfPixelBarrelHits
  (const TrackingParticle& simTrack)
{
  // How many pixel barrel hits?
  int numberOfPixelBarrelHits = 0;

  for(std::vector<PSimHit>::const_iterator simHit = simTrack.pSimHit_begin();
                                           simHit!= simTrack.pSimHit_end();
                                           simHit++)
  {
    DetId id = DetId(simHit->detUnitId());

    if(theTracker->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelBarrel)
      numberOfPixelBarrelHits++;
  }

  return numberOfPixelBarrelHits;
}


/*****************************************************************************/
void ChargedMultiplicityAnalyzer::checkSimTracks
  (const VertexCollection * vertices,
   edm::Handle<TrackingParticleCollection>& simCollection)
{
  for(TrackingParticleCollection::size_type i=0;
               i < simCollection.product()->size(); ++i)
  {
    const TrackingParticleRef simTrack(simCollection, i);

    vector<float> result;

    // sim
    result.push_back(vertices->size());                // nvtx
    result.push_back(simTrack->charge());              // q
    result.push_back(simTrack->pdgId());               // id
    result.push_back(simTrack->eta());                 // eta
    result.push_back(simTrack->pt());                  // pt
    result.push_back(simTrack->vertex().Rho());        // rho
    result.push_back(getNumberOfPixelBarrelHits(*simTrack)); // npb

    // fill
    trackSim->Fill(&result[0]);
  }
}

/*****************************************************************************/
bool ChargedMultiplicityAnalyzer::isAtEdge
  (const RectangularPixelTopology* topology, const SiPixelRecHit & recHit)
{
  return topology->isItEdgePixelInX(recHit.cluster()->minPixelCol()) ||
         topology->isItEdgePixelInX(recHit.cluster()->maxPixelCol()) ||
         topology->isItEdgePixelInY(recHit.cluster()->minPixelRow()) ||
         topology->isItEdgePixelInY(recHit.cluster()->maxPixelRow());
}

/*****************************************************************************/
void ChargedMultiplicityAnalyzer::analyze
  (const edm::Event& ev, const edm::EventSetup& es)
{
  // Get vertices
  edm::Handle<VertexCollection> vertexCollection;
  ev.getByType(vertexCollection);
  const VertexCollection * vertices = vertexCollection.product();

  // Get simulated
  edm::Handle<TrackingParticleCollection> simCollection;
  ev.getByType(simCollection);
  checkSimTracks(vertices, simCollection);

  if(vertices->size() == 1)
  {
    // Get the first and only vertex
    VertexCollection::const_iterator vertex = vertices->begin();

    // Get simulated
    edm::Handle<TrackingParticleCollection> simCollection;
    ev.getByType(simCollection);
  
    // Get pixel hits
    edm::Handle<SiPixelRecHitCollection> pixelHits;
    ev.getByType(pixelHits);
    theHits = pixelHits.product();

    // Get associator
    TrackerHitAssociator theHitAssociator(ev);
  
    // Take all units
    for(SiPixelRecHitCollection::id_iterator id = theHits->id_begin();
                                             id!= theHits->id_end(); id++)
    {
      SiPixelRecHitCollection::range range = theHits->get(*id);
      const PixelGeomDetUnit* pixelDet =
        dynamic_cast<const PixelGeomDetUnit*> (theTracker->idToDet(*id));

      const RectangularPixelTopology* topology =
        dynamic_cast<const RectangularPixelTopology*>
          (&(pixelDet->specificTopology()));
  
      // Take all hits
      if((*id).subdetId() == PixelSubdetector::PixelBarrel)
      { 
        PXBDetId pid(*id);
  
        // Look at the first barrel only
        if(pid.layer() == 1)
        for(SiPixelRecHitCollection::const_iterator recHit = range.first;
                             recHit!= range.second; recHit++)
        {
          GlobalPoint gpos = pixelDet->toGlobal(recHit->localPosition());
  
          // eta and eloss
          double theta = atan2(double(gpos.perp()),
                               double(gpos.z() - vertex->position().z()));
          float eta   = -log(tan(theta/2));
          float eloss = recHit->cluster()->charge();

          // associate 
          vector<PSimHit> simHits = theHitAssociator.associateHit(*recHit);
          const PSimHit * bestSimHit = 0;

          // background = 1, primary = 2, primary looper = 3
          int type = 1;
          for(vector<PSimHit>::const_iterator simHit = simHits.begin();
                      simHit!= simHits.end(); simHit++)
            if(simHit->processType() == 2)
            {
              bestSimHit = &(*simHit);
              type = 2;
            }

          // Is it a looping hit?
          if(type == 2)
          for(TrackingParticleCollection::size_type i=0;
               i < simCollection.product()->size(); ++i)
          {
            const TrackingParticleRef simTrack(simCollection, i);

            const PSimHit * firstSimHit = 0;
            float pmax = 0;
            for(std::vector<PSimHit>::const_iterator
                          simHit = simTrack->pSimHit_begin();
                          simHit!= simTrack->pSimHit_end(); simHit++) 
              if(simHit->pabs() > pmax)
              { pmax = simHit->pabs(); firstSimHit = &(*simHit); }

            if(firstSimHit != 0)
              if(bestSimHit->trackId() == firstSimHit->trackId())
              {
                if(!(bestSimHit->pabs() == firstSimHit->pabs()))
                  type = 3;
  
                break;
              }
          }

          vector<float> result;

          result.push_back(eta);
          result.push_back(eloss);
          result.push_back(type * (isAtEdge(topology,*recHit) ? -1 : 1));

          multi->Fill(&result[0]);
        }
      }
    }
  }
}

DEFINE_FWK_MODULE(ChargedMultiplicityAnalyzer);

