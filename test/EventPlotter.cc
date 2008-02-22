#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"


#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/VZero/interface/VZero.h"
#include "DataFormats/VZero/interface/VZeroFwd.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

/*
class TransientTrackFromFTSFactory {
 public:

    reco::TransientTrack build (const FreeTrajectoryState & fts) const;
    reco::TransientTrack build (const FreeTrajectoryState & fts,
        const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
};
*/

#include "RecoVertex/KalmanVertexFit/interface/SingleTrackVertexConstraint.h"

#include <fstream>

#include <vector>
using namespace std;

class EventPlotter : public edm::EDAnalyzer
{
  public:
//    EventPlotter(const edm::EventSetup& es, vector<string> trackCollectionLabels);
    explicit EventPlotter(const edm::ParameterSet& pset);
    ~EventPlotter();
    virtual void beginJob(const edm::EventSetup& es);
    virtual void analyze(const edm::Event& ev, const edm::EventSetup& es);
    virtual void endJob();

    void printEvent(const edm::Event& ev);

  private:
    void printHelix(const GlobalPoint& p1, const GlobalPoint& p2, const
GlobalVector& n2, ofstream& outFile, int charge);
    void getGlobal(const PSimHit& simHit, GlobalPoint& p, GlobalVector& v, bool& isPixel);
    void printSimTracks(const TrackingParticleCollection* simTracks);
    void printRecTracks(const reco::TrackCollection* recTracks,
                        const vector<Trajectory>*    recTrajes,
                        const reco::VertexCollection* vertices);
    pair<float,float> refitWithVertex (const reco::Track & recTrack,
                                       const reco::VertexCollection* vertices);

    void printPixelRecHit
      (const SiPixelRecHit * recHit, ofstream& pixelDetUnits,
                                     ofstream& pixelHits);
    void printPixelRecHits(const edm::Event& ev);

    void printStripRecHit
      (const SiStripRecHit2D * recHit, ofstream& stripDetUnits,
                                       ofstream& stripHits);
    void printStripRecHits(const edm::Event& ev);

    FreeTrajectoryState getTrajectory(const reco::Track& track);
    void printVZeros(const reco::VZeroCollection* vZeros);
    void printVertices(const reco::VertexCollection* vertices);

    const TrackerGeometry* theTracker;
    const MagneticField* theMagField;
    const TransientTrackBuilder * theTTBuilder;

    vector<string> trackCollectionLabels;

    int sign;

//   vector<string> trackCollectionLabels;
   bool zipFiles;
};

#include "FWCore/Framework/interface/ESHandle.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"


#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include <fstream>
#include <string>

using namespace std;
using namespace reco;
using namespace edm;

/*****************************************************************************/
struct sortByPabs
{
  bool operator() (const PSimHit& a, const PSimHit& b) const
  {
    return (a.pabs() > b.pabs());
  }
};

/*****************************************************************************/
//EventPlotter::EventPlotter(const edm::EventSetup& es, vector<string> trackCollectionLabels_)
EventPlotter::EventPlotter(const edm::ParameterSet& pset)
{
  trackCollectionLabels = pset.getParameter<vector<string> >("trackCollection");
  zipFiles  = pset.getParameter<bool>("zipFiles");
}

/*****************************************************************************/
EventPlotter::~EventPlotter()
{
}

/*****************************************************************************/
void EventPlotter::beginJob(const edm::EventSetup& es)
{
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  theTracker = tracker.product();

  // Get magnetic field
  edm::ESHandle<MagneticField> magField;
  es.get<IdealMagneticFieldRecord>().get(magField);
  theMagField = magField.product();

  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> builder;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  theTTBuilder = builder.product();
}

/*****************************************************************************/
void EventPlotter::endJob()
{
}

/*****************************************************************************/
void EventPlotter::printHelix
  (const GlobalPoint& p1, const GlobalPoint& p2,
   const GlobalVector& v2, ofstream& outFile, int charge)
{
  GlobalVector dp = p2 - p1;
  GlobalVector n2(-v2.y(),v2.x(),0.);
  n2 = n2.unit();

  double r = -0.5 * (dp.x()*dp.x() + dp.y()*dp.y()) /
                    (dp.x()*n2.x() + dp.y()*n2.y());
  GlobalPoint c = p2 + r * n2;

  double dphi = sqrt(2 * 0.1 / fabs(r)); // allowed deflection: 0.1 cm

  double phi = acos(( (p1-c).x()*(p2-c).x() +
                     (p1-c).y()*(p2-c).y() )/(r*r));

  if(dp.x()*v2.x() + dp.y()*v2.y() < 0) phi = 2*M_PI - phi;

  int nstep = (int)(phi/dphi) + 1;

  if(nstep > 1)
  {
    dphi = phi / nstep;
    double dz = (p2 - p1).z() / nstep;

  
    GlobalPoint P0 = p1;
    GlobalPoint P1;

    charge = ((p1 - c).x() * (p2 - c).y() - (p1 - c).y() * (p2 - c).x() > 0 ?
              -1 : 1);
    if(dp.x()*v2.x() + dp.y()*v2.y() < 0) charge = -charge;
  
    for(int i = 0; i < nstep; i++)
    {
      double a = -charge * (i+1)*dphi;
      double z = p1.z() + (i+1)*dz;

      double x = c.x() + cos(a)*(p1 - c).x() - sin(a)*(p1 - c).y();
      double y = c.y() + sin(a)*(p1 - c).x() + cos(a)*(p1 - c).y();
  
      P1 = GlobalPoint(x,y,z);
  
      outFile << ", Line[{{"<<P0.x()<<","<<P0.y()<<","<<P0.z()/2<<"}, {"
                            <<P1.x()<<","<<P1.y()<<","<<P1.z()/2<<"}}]" << endl; 
      P0 = P1;
    }
  }
  else
  {
    GlobalPoint P0 = p1;
    GlobalPoint P1 = p2;

    outFile << ", Line[{{"<<P0.x()<<","<<P0.y()<<","<<P0.z()/2<<"}, {"
                          <<P1.x()<<","<<P1.y()<<","<<P1.z()/2<<"}}]" << endl; 
  }
}

/*****************************************************************************/
void EventPlotter::getGlobal
  (const PSimHit& simHit, GlobalPoint& p, GlobalVector& v, bool& isPixel)
{
  DetId id = DetId(simHit.detUnitId());

  p = theTracker->idToDetUnit(id)->toGlobal(simHit.localPosition());
  v = theTracker->idToDetUnit(id)->toGlobal(simHit.localDirection());

  isPixel =
   (theTracker->idToDet(id)->subDetector() == GeomDetEnumerators::PixelBarrel ||
    theTracker->idToDet(id)->subDetector() == GeomDetEnumerators::PixelEndcap);
}

/*****************************************************************************/
void EventPlotter::printSimTracks
  (const TrackingParticleCollection* simTracks)
{
  ofstream pixelFile("event/pixelSimHits.m");
  ofstream stripFile("event/stripSimHits.m");
  ofstream trajeFile("event/simTracks.m");

  for(TrackingParticleCollection::const_iterator simTrack = simTracks->begin();
                                                 simTrack!= simTracks->end();
                                                 simTrack++)
  {
    vector<PSimHit> simHits = simTrack->trackPSimHit();

    // reorder with help of momentum
    sort(simHits.begin(), simHits.end(), sortByPabs());

    for(vector<PSimHit>::const_iterator simHit = simHits.begin();
                                        simHit!= simHits.end(); simHit++)
    {
      GlobalPoint p1; GlobalVector v1; bool isPixel;
      getGlobal(*simHit, p1,v1, isPixel);

      if(isPixel)
        pixelFile << ", Point[{" << p1.x() << "," << p1.y()
                  << "," << p1.z()/2 << "}]"
//         << " Text[StyleForm[\"s\", FontFamily->\"Helvetica\", URL -> \"position" << p1 << "\"],"
//         << "{" << p1.x() << "," << p1.y()
//                  << "," << p1.z()/2 << "}, {0,1}]"
                  << endl;
      else
        stripFile << ", Point[{" << p1.x() << "," << p1.y()
                  << "," << p1.z()/2 << "}]"
//          << " Text[StyleForm[\"s\", FontFamily->\"Helvetica\", URL -> \"position" << p1 << "\"],"
//         << "{" << p1.x() << "," << p1.y()
//                  << "," << p1.z()/2 << "}, {0,1}]"
                  << endl;

      if(simHit == simHits.begin()) // vertex to first point
      {
        GlobalPoint p0(simTrack->vertex().x(),
                       simTrack->vertex().y(),
                       simTrack->vertex().z());
        printHelix(p0,p1,v1, trajeFile, simTrack->charge());
      }

      if(simHit+1 != simHits.end()) // if not last
      {
        GlobalPoint p2; GlobalVector v2;
        getGlobal(*(simHit+1), p2,v2, isPixel);
        printHelix(p1,p2,v2, trajeFile, simTrack->charge());
      }
    }
  }

  pixelFile.close();
  stripFile.close();
  trajeFile.close();
}

/*****************************************************************************/
pair<float,float> EventPlotter::refitWithVertex
  (const reco::Track & recTrack,
   const reco::VertexCollection* vertices)
{
  TransientTrack theTransientTrack = theTTBuilder->build(recTrack);

  // If there are vertices found
  if(vertices->size() > 0)
  {
    float dzmin = -1.;
    const reco::Vertex * closestVertex = 0;

    // Look for the closest vertex in z
    for(reco::VertexCollection::const_iterator
        vertex = vertices->begin(); vertex!= vertices->end(); vertex++)
    {
      float dz = fabs(recTrack.vertex().z() - vertex->position().z());
      if(vertex == vertices->begin() || dz < dzmin)
      { dzmin = dz ; closestVertex = &(*vertex); }
    }


    // Get vertex position and error matrix
    GlobalPoint vertexPosition(closestVertex->position().x(),
                               closestVertex->position().y(),
                               closestVertex->position().z());

    float beamSize = 15e-4; // 15 um
    GlobalError vertexError(beamSize*beamSize, 0,
                            beamSize*beamSize, 0,
                            0,closestVertex->covariance(2,2));

    // Refit track with vertex constraint
    SingleTrackVertexConstraint stvc;
    pair<TransientTrack, float> result =
      stvc.constrain(theTransientTrack, vertexPosition, vertexError);

    return pair<float,float>(result.first.impactPointTSCP().pt(),
                             result.second);
  }
  else
    return pair<float,float>(recTrack.pt(), -9999);
}

/*****************************************************************************/
void EventPlotter::printRecTracks
  (const reco::TrackCollection* recTracks,
   const vector<Trajectory>*    recTrajes,
   const reco::VertexCollection* vertices)
{
  // Hits
  ofstream pixelFile; pixelFile.open("event/pixelRecHits.m",ios_base::app);
  ofstream stripFile; stripFile.open("event/stripRecHits.m",ios_base::app);

  int i = 0;
  for(reco::TrackCollection::const_iterator recTrack = recTracks->begin();
                                            recTrack!= recTracks->end();
                                            recTrack++, i++)
  {
    int j = 0;

/*
    float ptv  = refitWithVertex(*recTrack,vertices).first;
    float chi2 = refitWithVertex(*recTrack,vertices).second;
    if(fabs(recTrack->eta()) < 1.5 && fabs(ptv-1) < 0.001 && log10(chi2) < 1.5)
*/

    for(trackingRecHit_iterator recHit = recTrack->recHitsBegin();
                                recHit!= recTrack->recHitsEnd();
                                recHit++, j++)
    if((*recHit)->isValid())
    {
      DetId id = (*recHit)->geographicalId();
      LocalPoint lpos = (*recHit)->localPosition();
      GlobalPoint p = theTracker->idToDet(id)->toGlobal(lpos);

      if(theTracker->idToDet(id)->subDetector() ==
           GeomDetEnumerators::PixelBarrel ||
         theTracker->idToDet(id)->subDetector() ==
           GeomDetEnumerators::PixelEndcap) 
      { 
        const SiPixelRecHit* pixelRecHit =
          dynamic_cast<const SiPixelRecHit *>(&(**recHit));
        SiPixelRecHit::ClusterRef const& cluster = pixelRecHit->cluster();
        vector<SiPixelCluster::Pixel> pixels = cluster->pixels();

if(0)
{
  DetId id = pixelRecHit->geographicalId();
         
  // DetUnit
  double x = theTracker->idToDet(id)->surface().bounds().width() /2;
  double y = theTracker->idToDet(id)->surface().bounds().length()/2;
  double z = 0.;

  GlobalPoint p00 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x,-y,z));
  GlobalPoint p01 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x, y,z));
  GlobalPoint p10 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x,-y,z));
  GlobalPoint p11 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x, y,z));
  
  pixelFile
    << ", Line[{{"<< p00.x()<<","<<p00.y()<<","<<p00.z()/2<<"}, {"
                  << p01.x()<<","<<p01.y()<<","<<p01.z()/2<<"}}]"
    << ", Line[{{"<< p01.x()<<","<<p01.y()<<","<<p01.z()/2<<"}, {"
                  << p11.x()<<","<<p11.y()<<","<<p11.z()/2<<"}}]"
    << ", Line[{{"<< p11.x()<<","<<p11.y()<<","<<p11.z()/2<<"}, {"
                  << p10.x()<<","<<p10.y()<<","<<p10.z()/2<<"}}]"
    << ", Line[{{"<< p10.x()<<","<<p10.y()<<","<<p10.z()/2<<"}, {"
                  << p00.x()<<","<<p00.y()<<","<<p00.z()/2<<"}}]" << endl;
}

        ostringstream o; o << (sign > 0 ? "" : "-") << i; 
        string pixelInfo = " Text[StyleForm[\""  + o.str() +
                           "\", FontFamily->\"Helvetica\", URL -> \"Track " + o.str();


        if(theTracker->idToDet(id)->subDetector() ==
             GeomDetEnumerators::PixelBarrel)
        {
          // 0 + (layer-1)<<1 + (ladder-1)%2 : 0-5
          PXBDetId pid(id);
          ostringstream o;
          o << " (" << pid.layer()  << "|" << pid.ladder()
            <<  "|" << pid.module() << ")";
          pixelInfo += o.str();
        }
        else
        {
          // 6 + (disk-1)<<1 + (panel-1)%2
          PXFDetId pid(id);
          ostringstream o;
          o << " (" << pid.side()   << "|" << pid.disk()
            <<  "|" << pid.blade()  << "|" << pid.panel()
            <<  "|" << pid.module() << ")";
          pixelInfo += o.str();
        }
 
        for(vector<SiPixelCluster::Pixel>::const_iterator
          pixel = pixels.begin(); pixel!= pixels.end(); pixel++)
        {
          ostringstream o; o << "[" << int(pixel->x)
                             << " " << int(pixel->y)
                             << " " << int(pixel->adc/135) << "]"; 
          pixelInfo += " " + o.str();
        }
        pixelInfo += "\"],";

        pixelFile
         << ", Point[{" << p.x() << "," << p.y() << "," << p.z()/2 << "}],"
         << pixelInfo
//         << " Text[StyleForm[\"" << i << "\", URL -> \"Track " << i << "\"],"
         << " {" << p.x() << "," << p.y() << "," << p.z()/2 << "},"
         << " {" << 0 << "," << -1 << "}]" << endl;
      }
      else
      {
        stripFile
         << ", Point[{" << p.x() << "," << p.y() << "," << p.z()/2 << "}],"
         << " Text[StyleForm[\"" << i << "\", FontFamily->\"Helvetica\", URL -> \"Track " << i << "\"],"
         << " {" << p.x() << "," << p.y() << "," << p.z()/2 << "},"
         << " {" << 3-i%6 << "," << -1 << "}]" << endl;
      }
    }
  }
  pixelFile.close();
  stripFile.close();

  // Trajectory
  ofstream trajeFile; trajeFile.open("event/recTracks.m",ios_base::app);

  reco::TrackCollection::const_iterator recTrack = recTracks->begin();

  for(vector<Trajectory>::const_iterator it = recTrajes->begin();
                                         it!= recTrajes->end(); it++)
  {
    vector<TrajectoryMeasurement> meas = it->measurements();

/*
    float ptv  = refitWithVertex(*recTrack,vertices).first;
    float chi2 = refitWithVertex(*recTrack,vertices).second;
    if(fabs(ptv-1) < 0.01 && log10(chi2) > 2.0) // 1.5
*/

    for(vector<TrajectoryMeasurement>::reverse_iterator im = meas.rbegin();
                                                        im!= meas.rend(); im++)
    {
      if(im == meas.rbegin())
      {
        GlobalPoint  p2 = (*(im  )).updatedState().globalPosition();
        GlobalVector v2 = (*(im  )).updatedState().globalDirection();
        GlobalPoint  p1 = GlobalPoint(recTrack->vertex().x(),
                                      recTrack->vertex().y(),
                                      recTrack->vertex().z());

        printHelix(p1,p2,v2, trajeFile, recTrack->charge());
      }

      if(im+1 != meas.rend())
      {
        GlobalPoint  p1 = (*(im  )).updatedState().globalPosition();
        GlobalPoint  p2 = (*(im+1)).updatedState().globalPosition();
        GlobalVector v2 = (*(im+1)).updatedState().globalDirection();
  
        printHelix(p1,p2,v2, trajeFile, recTrack->charge());
      }
    }
   
    recTrack++;
  }
  trajeFile.close();
}

/*****************************************************************************/
FreeTrajectoryState EventPlotter::getTrajectory(const reco::Track& track)
{
  GlobalPoint position(track.vertex().x(),
                       track.vertex().y(),
                       track.vertex().z());

  GlobalVector momentum(track.momentum().x(),
                        track.momentum().y(),
                        track.momentum().z());

  GlobalTrajectoryParameters gtp(position,momentum,
                                 track.charge(),theMagField);

  FreeTrajectoryState fts(gtp);

  return fts;
}

/*****************************************************************************/
void EventPlotter::printVZeros
  (const reco::VZeroCollection* vZeros)
{
  ClosestApproachInRPhi theApproach;

  ofstream outFile;
  outFile.open("event/vZeros.m");

  for(reco::VZeroCollection::const_iterator vZero = vZeros->begin();
                                            vZero!= vZeros->end();
                                            vZero++)
  {
    // Calculate closest approach to beam-line
    GlobalPoint crossing(vZero->vertex().position().x(),
                         vZero->vertex().position().y(),
                         vZero->vertex().position().z());

    // Calculate momentum at dca
    FreeTrajectoryState posTraj = getTrajectory(*(vZero->positiveDaughter()));
    FreeTrajectoryState negTraj = getTrajectory(*(vZero->negativeDaughter()));
    theApproach.calculate(posTraj,negTraj);
    pair<GlobalTrajectoryParameters, GlobalTrajectoryParameters>
      gtp = theApproach.trajectoryParameters();
    pair<GlobalVector,GlobalVector>
      momenta(gtp.first.momentum(), gtp.second.momentum());

    GlobalVector momentum = momenta.first + momenta.second;

    GlobalVector r_(crossing.x(),crossing.y(),0);
    GlobalVector p_(momentum.x(),momentum.y(),0);

    GlobalVector r (crossing.x(),crossing.y(),crossing.z());
    GlobalVector p (momentum.x(),momentum.y(),momentum.z());
    GlobalVector b  = r  - (r_*p_)*p  / p_.mag2();

    outFile << ", Line[{{" << vZero->vertex().position().x()
                    << "," << vZero->vertex().position().y()
                    << "," << vZero->vertex().position().z()/2
                 << "}, {" << b.x()
                    << "," << b.y()
                    << "," << b.z()/2 << "}}]" << endl;
  }

  outFile.close();
}

/*****************************************************************************/
void EventPlotter::printVertices
  (const reco::VertexCollection* vertices)
{
  ofstream outFile("event/vertices.m");

  for(reco::VertexCollection::const_iterator vertex = vertices->begin();
                                             vertex!= vertices->end();
                                             vertex++)
    outFile << ", Point[{" << vertex->position().x()
                    << "," << vertex->position().y()
                    << "," << vertex->position().z()/2 << "}]" << endl;

  outFile.close();
}

/*****************************************************************************/
void EventPlotter::printStripRecHit
 (const SiStripRecHit2D * recHit, ofstream& stripDetUnits, ofstream& stripHits)
{
  DetId id = recHit->geographicalId();

  // DetUnit
  double x = theTracker->idToDet(id)->surface().bounds().width() /2;
  double y = theTracker->idToDet(id)->surface().bounds().length()/2;
  double z = 0.;

  GlobalPoint p00 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x,-y,z));
  GlobalPoint p01 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x, y,z));
  GlobalPoint p10 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x,-y,z));
  GlobalPoint p11 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x, y,z));

  stripDetUnits
    << ", Line[{{"<< p00.x()<<","<<p00.y()<<","<<p00.z()/2<<"}, {"
                  << p01.x()<<","<<p01.y()<<","<<p01.z()/2<<"}}]"
    << ", Line[{{"<< p01.x()<<","<<p01.y()<<","<<p01.z()/2<<"}, {"
                  << p11.x()<<","<<p11.y()<<","<<p11.z()/2<<"}}]"
    << ", Line[{{"<< p11.x()<<","<<p11.y()<<","<<p11.z()/2<<"}, {"
                  << p10.x()<<","<<p10.y()<<","<<p10.z()/2<<"}}]"
    << ", Line[{{"<< p10.x()<<","<<p10.y()<<","<<p10.z()/2<<"}, {"
                  << p00.x()<<","<<p00.y()<<","<<p00.z()/2<<"}}]" << endl;
 
  // RecHit
  LocalPoint lpos; GlobalPoint p;

  lpos = LocalPoint(recHit->localPosition().x(),
                 y, recHit->localPosition().z());
  p = theTracker->idToDet(id)->toGlobal(lpos);
  stripHits << ", Line[{{"<<p.x()<<","<<p.y()<<","<<p.z()/2<<"}, {";

  lpos = LocalPoint(recHit->localPosition().x(),
                -y, recHit->localPosition().z());
  p = theTracker->idToDet(id)->toGlobal(lpos);
  stripHits << ""<<p.x()<<","<<p.y()<<","<<p.z()/2<<"}}]" << endl;
}

/*****************************************************************************/
void EventPlotter::printPixelRecHit
 (const SiPixelRecHit * recHit, ofstream& pixelDetUnits, ofstream& pixelHits)
{
  DetId id = recHit->geographicalId();

  // DetUnit
  double x = theTracker->idToDet(id)->surface().bounds().width() /2;
  double y = theTracker->idToDet(id)->surface().bounds().length()/2;
  double z = 0.;

  GlobalPoint p00 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x,-y,z));
  GlobalPoint p01 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x, y,z));
  GlobalPoint p10 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x,-y,z));
  GlobalPoint p11 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x, y,z));

  pixelDetUnits
    << ", Line[{{"<< p00.x()<<","<<p00.y()<<","<<p00.z()/2<<"}, {"
                  << p01.x()<<","<<p01.y()<<","<<p01.z()/2<<"}}]"
    << ", Line[{{"<< p01.x()<<","<<p01.y()<<","<<p01.z()/2<<"}, {"
                  << p11.x()<<","<<p11.y()<<","<<p11.z()/2<<"}}]"
    << ", Line[{{"<< p11.x()<<","<<p11.y()<<","<<p11.z()/2<<"}, {"
                  << p10.x()<<","<<p10.y()<<","<<p10.z()/2<<"}}]"
    << ", Line[{{"<< p10.x()<<","<<p10.y()<<","<<p10.z()/2<<"}, {"
                  << p00.x()<<","<<p00.y()<<","<<p00.z()/2<<"}}]" << endl;

  // RecHit
  LocalPoint lpos; GlobalPoint p;
  
  lpos = LocalPoint(recHit->localPosition().x(),
                    recHit->localPosition().y(),
                    recHit->localPosition().z());
  p = theTracker->idToDet(id)->toGlobal(lpos);
  pixelHits << ", Point[{"<<p.x()<<","<<p.y()<<","<<p.z()/2<<"}]" << endl;
}

/*****************************************************************************/
void EventPlotter::printPixelRecHits(const edm::Event& ev)
{
  ofstream pixelDetUnits, pixelHitsSingle;
  pixelHitsSingle.open("event/pixelHitsSingle.m");
  pixelDetUnits.open  ("event/pixelDetUnits.m");

  // Get pixel hit collections
  vector<edm::Handle<SiPixelRecHitCollection> > pixelColls;
  ev.getManyByType(pixelColls);

  for(vector<edm::Handle<SiPixelRecHitCollection> >::const_iterator
      pixelColl = pixelColls.begin();
      pixelColl!= pixelColls.end(); pixelColl++)
  {
    const SiPixelRecHitCollection* thePixelHits = (*pixelColl).product();

    for(SiPixelRecHitCollection::id_iterator
          id = thePixelHits->id_begin(); id!= thePixelHits->id_end(); id++)
    {
      SiPixelRecHitCollection::range range = thePixelHits->get(*id);

      // Take all hits
      for(SiPixelRecHitCollection::const_iterator
            recHit = range.first; recHit!= range.second; recHit++)
      if(recHit->isValid())
        printPixelRecHit(&(*recHit), pixelDetUnits, pixelHitsSingle);
    }
  }
}

/*****************************************************************************/
void EventPlotter::printStripRecHits(const edm::Event& ev)
{
  ofstream stripDetUnits, stripHitsSingle, stripHitsMatched;
  stripDetUnits.open   ("event/stripDetUnits.m");
  stripHitsSingle.open( "event/stripHitsSingle.m");
  stripHitsMatched.open("event/stripHitsMatched.m");

  {
  // Get strip hit collections
  vector<edm::Handle<SiStripRecHit2DCollection> > stripColls;
  ev.getManyByType(stripColls);

  for(vector<edm::Handle<SiStripRecHit2DCollection> >::const_iterator
      stripColl = stripColls.begin();
      stripColl!= stripColls.end(); stripColl++)
  {
    const SiStripRecHit2DCollection* theStripHits = (*stripColl).product();

    for(SiStripRecHit2DCollection::id_iterator
          id = theStripHits->id_begin(); id!= theStripHits->id_end(); id++)
    {
      SiStripRecHit2DCollection::range range = theStripHits->get(*id);

      // Take all hits
      for(SiStripRecHit2DCollection::const_iterator
            recHit = range.first; recHit!= range.second; recHit++)
      if(recHit->isValid())
        printStripRecHit(&(*recHit), stripDetUnits, stripHitsSingle);
    }
  }
  }

  // Get matched strip hit collections
  {
  vector<edm::Handle<SiStripMatchedRecHit2DCollection> > stripColls;
  ev.getManyByType(stripColls);

  for(vector<edm::Handle<SiStripMatchedRecHit2DCollection> >::const_iterator
      stripColl = stripColls.begin();
      stripColl!= stripColls.end(); stripColl++)
  {
    const SiStripMatchedRecHit2DCollection* theStripHits = (*stripColl).product();

    for(SiStripMatchedRecHit2DCollection::id_iterator
          id = theStripHits->id_begin(); id!= theStripHits->id_end(); id++)
    {
      SiStripMatchedRecHit2DCollection::range range = theStripHits->get(*id);

      // Take all hits
      for(SiStripMatchedRecHit2DCollection::const_iterator
            recHit = range.first; recHit!= range.second; recHit++)
      {
        if(recHit->monoHit()->isValid())
          printStripRecHit((recHit->monoHit()), stripDetUnits,stripHitsSingle);
        if(recHit->stereoHit()->isValid())
          printStripRecHit((recHit->stereoHit()), stripDetUnits,stripHitsSingle);

        DetId id = recHit->geographicalId();
        LocalPoint lpos = recHit->localPosition();
        GlobalPoint p = theTracker->idToDet(id)->toGlobal(lpos);

        stripHitsMatched
         << ", Point[{" << p.x() << "," << p.y() << "," << p.z()/2 << "}]"
         << endl;
      }
    }
  }
  }

  stripDetUnits.close();
  stripHitsSingle.close();
  stripHitsMatched.close();
}

/*****************************************************************************/
//void EventPlotter::printEvent(const edm::Event& ev)
void EventPlotter::analyze(const edm::Event& ev, const edm::EventSetup& es)
{
  // Get simulated
  edm::Handle<TrackingParticleCollection> simCollection;
//  ev.getByType(simCollection);
  ev.getByLabel("trackingtruthprod",simCollection);

  // Get reconstructed
//  edm::Handle<reco::TrackCollection>  recCollection;
//  ev.getByLabel(trackCollectionLabel, recCollection);

//  edm::Handle<vector<Trajectory> >    trajCollection;
//  ev.getByLabel(trackCollectionLabel, trajCollection);

  // Get vzeros
  edm::Handle<reco::VZeroCollection> vZeroCollection;
  ev.getByType(vZeroCollection);
  
  // Get vertices
  edm::Handle<reco::VertexCollection> vertexCollection;
  ev.getByLabel("pixelVertices",vertexCollection);

  // Get products
  const TrackingParticleCollection* simTracks = simCollection.product();
//  const reco::TrackCollection*      recTracks = recCollection.product();
  const reco::VZeroCollection*      vZeros    = vZeroCollection.product();
//  const vector<Trajectory>*         recTrajes = trajCollection.product();
  const reco::VertexCollection*     vertices  = vertexCollection.product();

  ofstream labelFile("event/label.m");
  labelFile << "run:"    << ev.id().run()
            << " event:" << ev.id().event()
//            << "   rec=" << recTracks->size()
            << " v0="    << vZeros->size()
            << " vtx="   << vertices->size();
//            << endl;
//  labelFile.close();

  printSimTracks(simTracks);

  { ofstream file("event/pixelRecHits.m"); file.close(); }
  { ofstream file("event/stripRecHits.m"); file.close(); }
  { ofstream file("event/recTracks.m");    file.close(); }

  sign = 1;
  for(vector<string>::iterator label = trackCollectionLabels.begin();
                               label!= trackCollectionLabels.end(); label++)
  {
    // Get reconstructed
    edm::Handle<reco::TrackCollection>  recCollection;
    ev.getByLabel(*label, recCollection);

    edm::Handle<vector<Trajectory> >    trajCollection;
    ev.getByLabel(*label, trajCollection);

    const reco::TrackCollection*      recTracks = recCollection.product();
    const vector<Trajectory>*         recTrajes = trajCollection.product();

    labelFile  << " rec=" << recTracks->size();

    cerr << "[EventPlotter] " << *label << " = " << recTracks->size() << endl;

    printRecTracks(recTracks, recTrajes, vertices);

    sign = -sign;
  }

  labelFile << endl;
  labelFile.close();

  printVZeros   (vZeros);
  printVertices (vertices);

  printPixelRecHits(ev);
  printStripRecHits(ev);

  if(zipFiles)
  {
    system("zip -q -m event.zip *.m");
    cerr << "[EventPlotter] event plotted (event.zip)" << endl;
  }
  else
    cerr << "[EventPlotter] event plotted (*.m)" << endl;

  while(getchar() == 0);
} 

DEFINE_FWK_MODULE(EventPlotter);
