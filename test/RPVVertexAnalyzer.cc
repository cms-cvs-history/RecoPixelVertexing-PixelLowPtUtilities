#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"

#include <fstream>
using namespace std;

/*****************************************************************************/
class RPVVertexAnalyzer : public edm::EDAnalyzer
{
 public:
   explicit RPVVertexAnalyzer(const edm::ParameterSet& pset);
   ~RPVVertexAnalyzer();
   virtual void beginJob(const edm::EventSetup& es);
   virtual void analyze(const edm::Event& ev, const edm::EventSetup& es);
   virtual void endJob();

 private:
   int NEvents;
   int NTrkMin;
   double ZSeparation;

   TFile * resultFile;
   TNtuple * ntuple;
};

/*****************************************************************************/
RPVVertexAnalyzer::RPVVertexAnalyzer(const edm::ParameterSet& pset)
{
  NEvents     = pset.getParameter<int>("NEvents");
  NTrkMin     = pset.getParameter<int>("NTrkMin");
  ZSeparation = pset.getParameter<double>("ZSeparation"); 
}

/*****************************************************************************/
RPVVertexAnalyzer::~RPVVertexAnalyzer()
{
}

/*****************************************************************************/
void RPVVertexAnalyzer::beginJob(const edm::EventSetup& es)
{
  // Root
  char fileName[256];
  sprintf(fileName,"out/vertices_%d_%d_%.2f.root",NEvents,NTrkMin,ZSeparation);

  resultFile = new TFile(fileName,"recreate");
  resultFile->cd();

  ntuple = new TNtuple("vertices","vertices", "nrec:nvtx");
}

/*****************************************************************************/
void RPVVertexAnalyzer::endJob()
{
  resultFile->cd();
  ntuple->Write();
  resultFile->Close();
}

/*****************************************************************************/
void RPVVertexAnalyzer::analyze
  (const edm::Event& ev, const edm::EventSetup& es)
{
  // Get reconstructed
  edm::Handle<reco::TrackCollection> recCollection;
  ev.getByLabel("pixel3ProtoTracks",recCollection);
  const reco::TrackCollection* recTracks = recCollection.product();

  // Get vertices
  edm::Handle<reco::VertexCollection> vertexCollection;
  ev.getByType(vertexCollection);
  const reco::VertexCollection* vertices = vertexCollection.product();

  vector<float> result;

  result.push_back(recTracks->size());
  result.push_back(vertices->size());

  ntuple->Fill(&result[0]);
}

DEFINE_FWK_MODULE(RPVVertexAnalyzer);
