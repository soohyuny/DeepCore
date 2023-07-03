// -*- C++ -*-
//
// Package:    RecoTracker/DeepCoreTraining
// Class:      DeepCoreNtuplizer
//
/**\class DeepCoreNtuplizer DeepCoreNtuplizer.cc RecoTracker/DeepCoreTraining/plugins/DeepCoreNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Valerio Bertacchi
//         Created:  Fri, 09 Jul 2021 10:43:34 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"



#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"

#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

#include "TrackingTools/GeomPropagators/interface/StraightLinePlaneCrossing.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"



#include <boost/range.hpp>
#include <boost/foreach.hpp>
#include "boost/multi_array.hpp"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimG4Core/Notification/interface/G4SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

// #include "SimDataFormats/Vertex/interface/SimVertex.h"


#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"



#include "TTree.h"


//
// class declaration
//

class DeepCoreNtuplizer : public edm::stream::EDProducer<> {
public:
  explicit DeepCoreNtuplizer(const edm::ParameterSet&);
  ~DeepCoreNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  struct TrackAndState
  {
    TrackAndState(const reco::Track *aTrack, TrajectoryStateOnSurface aState) :
      track(aTrack), state(aState) {}
    const reco::Track*      track;
    TrajectoryStateOnSurface state;
  };
  
  template<typename Cluster>
  struct ClusterWithTracks
  {
    ClusterWithTracks(const Cluster &c) : cluster(&c) {}
    const Cluster* cluster;
    std::vector<TrackAndState> tracks;
  };
  
  typedef ClusterWithTracks<SiPixelCluster> SiPixelClusterWithTracks;

  typedef boost::sub_range<std::vector<SiPixelClusterWithTracks> > SiPixelClustersWithTracks;

  TFile* DeepCoreNtuplizer_out;
  TTree* DeepCoreNtuplizerTree;
  static const int jetDimX =30;//200
  static const int jetDimY =30;//200
  static const int Nlayer =7;//4;
  static const int Ntrack = 100;
  static const int Npar = 5; //added 1/pt
  static const int Nover = 3;
  double clusterMeas[jetDimX][jetDimY][Nlayer];
  double trackPar[jetDimX][jetDimY][Nover][Npar+1]; //NOFLAG
  double trackProb[jetDimX][jetDimY][Nover];
  double clusterSplit[Ntrack][Nlayer][jetDimX][jetDimY];
  double track_pt[Ntrack] = {0.0};
  double track_pz[Ntrack] = {0.0};
  double track_eta[Ntrack] = {0.0};
  double track_phi[Ntrack] = {0.0};
  double jet_pt;
  double jet_p;
  double jet_eta;
  double jet_phi;
  double jet_Ntrack = 0;
  int NPixLay[Nlayer] = {0};
  double absEtaMin = 0.;
  double absEtaMax= 1000.;


  int eventID;


  double pitchX = 0.01;
  double pitchY = 0.015;
//  static const int distThr = 2;//4
  const int distThr = 2;//4



private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------
  
  std::string propagatorName_;
//  edm::ESHandle<MagneticField>          magfield_;
//  edm::ESHandle<GlobalTrackingGeometry> geometry_.;
//  edm::ESHandle<Propagator>             propagator_;

  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>          magfield_;
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geometryToken_;
  edm::ESGetToken<Propagator,TrackingComponentsRecord>             propagator_;
  
  edm::ESGetToken<PixelClusterParameterEstimator, TkPixelCPERecord> parEstHandle;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoHandle; 

  edm::EDGetTokenT<std::vector<reco::Vertex> > vertices_;
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelClusters_;
  std::vector<SiPixelClusterWithTracks> allSiPixelClusters;
  std::map<uint32_t, SiPixelClustersWithTracks> siPixelDetsWithClusters;
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > inputPixelClusters;
  edm::EDGetTokenT<edm::View<reco::Candidate> > cores_;
  edm::EDGetTokenT<std::vector<SimTrack> > simtracksToken;
  edm::EDGetTokenT<std::vector<PSimHit> > PSimHitToken;
  edm::Handle<std::vector<PSimHit> > simhits;
  edm::EDGetTokenT<std::vector<PSimHit> > PSimHitECToken;
  edm::Handle<std::vector<PSimHit> > simhitsEC;

  std::map<int, double [Nlayer][jetDimX][jetDimY]> trackMap;
    
  double ptMin_;
  double pMin_;
  double deltaR_;
  double chargeFracMin_;
  double centralMIPCharge_;
  std::string pixelCPE_;
  
  bool barrelTrain_;
  bool endcapTrain_;
  bool fullTrain_;



  std::pair<bool, Basic3DVector<float>> findIntersection(const GlobalVector & , const reco::Candidate::Point & ,const GeomDet*);

  void fillPixelMatrix(const SiPixelCluster &, int, auto, const GeomDet*);

  std::pair<int,int> local2Pixel(double, double, const GeomDet*);
  
  LocalPoint pixel2Local(int, int, const GeomDet*);

  // Pixel size fix:
 // std::pair<bool,bool> pixel2Size(int, int, const GeomDet*);
  
  int pixelFlipper(const GeomDet*);

  void fillTrackInfo(const reco::Candidate&, const auto &, auto, auto, auto, const GeomDet*, const PixelClusterParameterEstimator*, std::vector<PSimHit>, const TrackerTopology* const);

  const GeomDet* DetectorSelector(int ,const reco::Candidate& jet, GlobalVector,  const reco::Vertex& jetVertex, const TrackerTopology* const, const PixelClusterParameterEstimator*, const auto &);

  std::vector<GlobalVector> splittedClusterDirections(const reco::Candidate&, const TrackerTopology* const, auto pp, const reco::Vertex& jetVertex, int );


};

DeepCoreNtuplizer::DeepCoreNtuplizer(const edm::ParameterSet& iConfig) :
	  magfield_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
//	  geometry_(esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>()),
	  geometryToken_(esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>()),
	  propagator_(esConsumes<Propagator,TrackingComponentsRecord>()),

	  parEstHandle(esConsumes<PixelClusterParameterEstimator, TkPixelCPERecord>()),
	  tTopoHandle(esConsumes<TrackerTopology, TrackerTopologyRcd>()),

      vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      pixelClusters_(consumes<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("pixelClusters"))),
      cores_(consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("cores"))),
      simtracksToken(consumes<std::vector<SimTrack> >(iConfig.getParameter<edm::InputTag>("simTracks"))),
      PSimHitToken(consumes<std::vector<PSimHit> >(iConfig.getParameter<edm::InputTag>("simHit"))),
      PSimHitECToken(consumes<std::vector<PSimHit> >(iConfig.getParameter<edm::InputTag>("simHitEC"))),
      ptMin_(iConfig.getParameter<double>("ptMin")),
      pMin_(iConfig.getParameter<double>("pMin")),
      deltaR_(iConfig.getParameter<double>("deltaR")),
      chargeFracMin_(iConfig.getParameter<double>("chargeFractionMin")),
      centralMIPCharge_(iConfig.getParameter<double>("centralMIPCharge")),
      pixelCPE_(iConfig.getParameter<std::string>("pixelCPE")),
      barrelTrain_(iConfig.getParameter<bool>("barrelTrain")),
      endcapTrain_(iConfig.getParameter<bool>("endcapTrain")),
      fullTrain_(iConfig.getParameter<bool>("fullTrain"))
{
  
  //  usesResource("TFileService");
   edm::Service<TFileService> fileService;

   DeepCoreNtuplizerTree= fileService->make<TTree>("DeepCoreNtuplizerTree","DeepCoreNtuplizerTree");
   DeepCoreNtuplizerTree->Branch("cluster_measured",clusterMeas,"cluster_measured[30][30][7]/D");
   DeepCoreNtuplizerTree->Branch("trackPar", trackPar, "trackPar[30][30][3][6]/D"); //NOFLAG
   DeepCoreNtuplizerTree->Branch("trackProb", trackProb, "trackProb[30][30][3]/D");
   DeepCoreNtuplizerTree->Branch("jet_eta",&jet_eta);
   DeepCoreNtuplizerTree->Branch("jet_pt",&jet_pt);
   DeepCoreNtuplizerTree->Branch("jet_p",&jet_p);
   DeepCoreNtuplizerTree->Branch("jet_phi",&jet_phi);
   DeepCoreNtuplizerTree->Branch("jet_Ntrack",&jet_Ntrack);
   DeepCoreNtuplizerTree->Branch("track_pt",track_pt,"track_pt[100]/D");
   DeepCoreNtuplizerTree->Branch("track_pz",track_pz,"track_pz[100]/D");
   DeepCoreNtuplizerTree->Branch("track_eta",track_eta,"track_eta[100]/D");
   DeepCoreNtuplizerTree->Branch("track_phi",track_phi,"track_phi[100]/D");
   DeepCoreNtuplizerTree->Branch("event_ID",&eventID);

   //debug branches
   DeepCoreNtuplizerTree->Branch("NPixLay",NPixLay, "NPixLay[7]/I");





   /// dichiarare cosa produce  produces<asd
     for(int i=0; i<Nlayer; i++){ //NOFLAG
       for(int j=0; j<jetDimX; j++){
         for(int k=0; k<jetDimY; k++){
           if(j<jetDimX && k<jetDimY && i< Nlayer) clusterMeas[j][k][i] = 0.0;
           for(int m=0; m<Nover; m++){
             if(j<jetDimX && k<jetDimY && i< Npar+1 && m<Nover) trackPar[j][k][m][i] =0.0;
             if(j<jetDimX && k<jetDimY && m<Nover) trackProb[j][k][m] =0.0;
           }
         }
       }
     }
  
  if(barrelTrain_) {
    absEtaMax = 1.4;
  }
  if(endcapTrain_) {
    absEtaMin=1.4;
  }

  
}

DeepCoreNtuplizer::~DeepCoreNtuplizer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
#define foreach BOOST_FOREACH

// ------------ method called to produce the data  ------------
void DeepCoreNtuplizer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // for(jet){ //soglia pt 500 Gev
  //   save(jet->pt,eta)
  //   for(layer){
  //     find(inter = intersezione jet-layer)
  //     for(cluster on layer){
  //       if(cluser in finestra 30x30 intorno a inter){
  //         save(cluster->sig in posizione x,y,layer)
  //         (trackcluster[n],tracks[n])=split(cluster)
  //         for(trackcluster){ //cioè per ogni traccia che ha originato cluser
  //           save(trackcluster->sig in posizione x,y,layer,track)//con ampiezza tot cluster/Ntrack?
  //           save(trackcluster->track pt,eta,phi) //geom det, specific surface
  //         }
  //       }
  //     }
  //
  //   }
  // }

//  iSetup.get<IdealMagneticFieldRecord>().get( magfield_ );
//  iSetup.get<GlobalTrackingGeometryRecord>().get(geometry_.);
//  iSetup.get<TrackingComponentsRecord>().get( "AnalyticalPropagator", propagator_ ); //AnalyticalPropagator???

  iEvent.getByToken(pixelClusters_, inputPixelClusters);
  allSiPixelClusters.clear(); siPixelDetsWithClusters.clear();
  allSiPixelClusters.reserve(inputPixelClusters->dataSize()); // this is important, otherwise push_back invalidates the iterators

  edm::Handle<std::vector<SimTrack> > simtracks;
  iEvent.getByToken(simtracksToken, simtracks);

  Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(vertices_, vertices);

  iEvent.getByToken(PSimHitToken, simhits);
  iEvent.getByToken(PSimHitECToken, simhitsEC);


  Handle<edm::View<reco::Candidate> > cores;
  iEvent.getByToken(cores_, cores);


//  edm::ESHandle<PixelClusterParameterEstimator> parEstHandle;
//  const PixelClusterParameterEstimator* parEst;
//  iSetup.get<TkPixelCPERecord>().get(pixelCPE_, parEstHandle);
//  parEst = parEstHandle.product();
  const PixelClusterParameterEstimator* parEst = &iSetup.getData(parEstHandle);  //pixelCPE_???
  
//  edm::ESHandle<TrackerTopology> tTopoHandle;
//  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
//  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerTopology* tTopo = &iSetup.getData(tTopoHandle);

  const GlobalTrackingGeometry* geometry_ = &iSetup.getData(geometryToken_);


  auto output = std::make_unique<edmNew::DetSetVector<SiPixelCluster>>();



  edm::LogInfo("EventMonitor|seeding") << "number of jets in the event=" << cores->size();

  for (unsigned int ji = 0; ji < cores->size(); ji++) { //loop jet
    if ((*cores)[ji].pt() > ptMin_  && (*cores)[ji].p() > pMin_ && std::abs((*cores)[ji].eta())>absEtaMin && std::abs((*cores)[ji].eta())<absEtaMax ) { 
      edm::LogInfo("EventMonitor|seeding") << "jet in acceptance" << ", pt= "<< (*cores)[ji].pt() << ", eta=" << (*cores)[ji].eta() << ", p=" <<(*cores)[ji].p(); 
      
      const reco::Candidate& jet = (*cores)[ji];
      const reco::Vertex& jetVertex = (*vertices)[0];
      GlobalVector jetDirection(jet.px(), jet.py(), jet.pz());


      std::vector<GlobalVector> splitClustDirSet = splittedClusterDirections(jet, tTopo, parEst, jetVertex, 1);
      edm::LogInfo("EventMonitor|seeding") << "number of cluster to split on lay1="<<splitClustDirSet.size() << "+jetDir" ;
      if(splitClustDirSet.size()==0) {
        splitClustDirSet = splittedClusterDirections(jet, tTopo, parEst, jetVertex, 2);
        edm::LogInfo("EventMonitor|seeding") << "missing clusters on lay1, used lay2 to find directions, number of cluster=" << splitClustDirSet.size() << "+jetDir" << std::endl;
      }
      splitClustDirSet.push_back(jetDirection);

      for(int cc=0; cc<(int)splitClustDirSet.size(); cc++){
      for(int lp=0; lp<Nlayer; lp++) NPixLay[lp] = 0;

      GlobalVector bigClustDir = splitClustDirSet.at(cc);

      std::set<int> trkIDset;
      LocalPoint jetInter(0,0,0);
      const auto & simtracksVector = simtracks.product();

      jet_pt = jet.pt();
      jet_p = jet.p();
      jet_eta = jet.eta();
      jet_phi = jet.phi();
      eventID= iEvent.id().event();

      auto jetVert = jetVertex; //trackInfo filling
      std::vector<std::map<int,SiPixelCluster>> clusterMapVector;
      std::vector<PSimHit> goodSimHit;

      edmNew::DetSetVector<SiPixelCluster>::const_iterator detIt = inputPixelClusters->begin();

      const GeomDet* globDet = DetectorSelector(2, jet, bigClustDir, jetVertex, tTopo,parEst, simtracksVector); //select detector mostly hitten by the jet

      //fill the good det vector
      std::vector<const GeomDet*> goodDets;
      for(int l=1; l<Nlayer+1; l++){
        const GeomDet* goodDet = DetectorSelector(l, jet, bigClustDir, jetVertex, tTopo,parEst, simtracksVector);
        goodDets.push_back(goodDet);

      }

      for (; detIt != inputPixelClusters->end(); detIt++) { //loop deset
        const edmNew::DetSet<SiPixelCluster>& detset = *detIt;
        const GeomDet* det = geometry_->idToDet(detset.id()); 
//        const GeomDet* det = geometry_.idToDet(detset.id()); 

        bool goodDetBool = false;
        for(int l=0; l<(int)goodDets.size(); l++){
          if(goodDets.at(l)==det) goodDetBool = true;
        }

        for (auto cluster = detset.begin(); cluster != detset.end(); cluster++) { //loop cluster

          const SiPixelCluster& aCluster = *cluster;
          det_id_type aClusterID= detset.id();

          // if(DetId(aClusterID).subdetId()!=1) continue; #select only  barrel clusters
          int lay = tTopo->layer(det->geographicalId());
          if(det->geographicalId().subdetId()==PixelSubdetector::PixelEndcap) {
            lay=lay+4; //endcap layer counting = 5,6,7
            }

          std::pair<bool, Basic3DVector<float>> interPair = findIntersection(bigClustDir,(reco::Candidate::Point)jetVertex.position(), det);

          if(interPair.first==false) continue;

          Basic3DVector<float> inter = interPair.second;

          auto localInter = det->specificSurface().toLocal((GlobalPoint)inter);

          LocalPoint cPos_local = parEst->localParametersV(aCluster,(*geometry_->idToDetUnit(detIt->id())))[0].first;
//          LocalPoint cPos_local = parEst->localParametersV(aCluster,(*geometry_.idToDetUnit(detIt->id())))[0].first;

           edm::LogInfo("EventMonitor|seeding").log([&](auto & logger){
              GlobalPoint pointVertex(jetVertex.position().x(), jetVertex.position().y(), jetVertex.position().z());
              GlobalPoint cPos = det->surface().toGlobal(parEst->localParametersV(aCluster,(*geometry_->idToDetUnit(detIt->id())))[0].first);
//              GlobalPoint cPos = det->surface().toGlobal(parEst->localParametersV(aCluster,(*geometry_.idToDetUnit(detIt->id())))[0].first);
              // GlobalVector intersectionDir = (GlobalPoint) inter - pointVertex;
              // auto deltaR_glob = Geom::deltaR(bigClustDir, intersectionDir);
              GlobalVector clusterDir = cPos - pointVertex;
              auto deltaR_clust = Geom::deltaR(bigClustDir, clusterDir);
              logger << "cluster with DeltaR(clust,window center direction)" << deltaR_clust;
           });

          if(std::abs(cPos_local.x()-localInter.x())/pitchX<=jetDimX/2 && std::abs(cPos_local.y()-localInter.y())/pitchY<=jetDimY/2){ //cPos_local.x,y use barycenter, other solutions can be developed
            NPixLay[lay-1]= NPixLay[lay-1]+1;
            if(goodDetBool) {
              fillPixelMatrix(aCluster,lay,localInter, det);
              }
          }//cluster in ROI
        } //cluster
      } //detset


      // ------------------------  good sim hit selection (hit inside window) --------------------------------------//
      std::vector<PSimHit> simhitsALL;
      std::vector<PSimHit>::const_iterator shItB = simhits->begin();
      for (; shItB != simhits->end(); shItB++) { //loop on simit to find correspondent det (barrel)
        simhitsALL.push_back((*shItB));
      }
      // std::vector<PSimHit>::const_iterator shItEC = simhitsEC->begin(); #needed if you want the target with simhitsEC 
      // for (; shItEC != simhitsEC->end(); shItEC++) { //loop on simit to find correspondent det (endCaps)
      //   simhitsALL.push_back((*shItEC));
      // }

      std::vector<PSimHit>::const_iterator shIt = simhitsALL.begin();
      // std::set<const GeomDet*> simhitsDetSet;
      for (; shIt != simhitsALL.end(); shIt++) { //loop deset
        // const edmNew::DetSet<PSimHit>& detset = *shIt;
        const GeomDet* det = geometry_->idToDet((*shIt).detUnitId());
//        const GeomDet* det = geometry_.idToDet((*shIt).detUnitId());
        if(det!=globDet) continue;
        std::pair<bool, Basic3DVector<float>> interPair = findIntersection(bigClustDir,(reco::Candidate::Point)jetVertex.position(), det);
        if(interPair.first==false) continue;
        Basic3DVector<float> inter = interPair.second;
        auto localInter = det->specificSurface().toLocal((GlobalPoint)inter);
        if(jetInter.x()==0 && jetInter.y()==0 && jetInter.z()==0) jetInter = localInter; //filling infoTracks

        if(std::abs(((*shIt).localPosition()).x()-localInter.x())/pitchX<=jetDimX/2 && std::abs(((*shIt).localPosition()).y()-localInter.y())/pitchY<=jetDimY/2){
          edm::LogInfo("EventMonitor|seeding") << "good sim hit on layer 2,  track ID" << (*shIt).trackId();
          goodSimHit.push_back((*shIt));

        }
      }

    fillTrackInfo(jet, simtracksVector, jetInter, bigClustDir, jetVert, globDet, parEst, goodSimHit,tTopo);



    clusterMapVector.clear();
    DeepCoreNtuplizerTree->Fill();
    edm::LogInfo("EventMonitor|seeding") << "N.pixel on layer 0=" <<NPixLay[0] << ", 1=" <<NPixLay[1] <<", 2=" <<NPixLay[1] << ", 3=" <<NPixLay[2]<<", 4=" <<NPixLay[4] <<", 5=" <<NPixLay[5] << ", 6=" <<NPixLay[6]<<std::endl;



    for(int i=0; i<Nlayer; i++){
      for(int j=0; j<jetDimX; j++){
        for(int k=0; k<jetDimY; k++){
          if(j<jetDimX && k<jetDimY && i< Nlayer) clusterMeas[j][k][i] = 0.0;
          for(int m=0; m<Nover; m++){
            if(trackPar[j][k][m][i]!=0 && j<jetDimX && k<jetDimY && i< Npar+1 && m<Nover) trackPar[j][k][m][i] =0.0; //NOFLAG
            if(trackProb[j][k][m]!=0 && j<jetDimX && k<jetDimY && m<Nover) trackProb[j][k][m] =0.0;
          }
          for(int h=0; h<Ntrack; h++){
            if(h<Ntrack) track_pt[h] = 0.0;
            if(h<Ntrack) track_pz[h] = 0.0;
            if(h<Ntrack) track_eta[h] = 0.0;
            if(h<Ntrack) track_phi[h] = 0.0;
          }
        }
      }
    }
    trkIDset.clear();
  } //bigcluster
  } //jet > pt
 } //jet
trackMap.clear();
}




  std::pair<bool, Basic3DVector<float>> DeepCoreNtuplizer::findIntersection(const GlobalVector & dir,const  reco::Candidate::Point & vertex, const GeomDet* det){
     StraightLinePlaneCrossing vertexPlane(Basic3DVector<float>(vertex.x(),vertex.y(),vertex.z()), Basic3DVector<float>(dir.x(),dir.y(),dir.z()));

     std::pair<bool, Basic3DVector<float>> pos = vertexPlane.position(det->specificSurface());

     return pos;
  }



  std::pair<int,int> DeepCoreNtuplizer::local2Pixel(double locX, double locY, const GeomDet* det){
    LocalPoint locXY(locX,locY);
    float pixX=(dynamic_cast<const PixelGeomDetUnit*>(det))->specificTopology().pixel(locXY).first;
    float pixY=(dynamic_cast<const PixelGeomDetUnit*>(det))->specificTopology().pixel(locXY).second;
    std::pair<int, int> out(pixX,pixY);
    return out;
  }

  LocalPoint DeepCoreNtuplizer::pixel2Local(int pixX, int pixY, const GeomDet* det){
    float locX=(dynamic_cast<const PixelGeomDetUnit*>(det))->specificTopology().localX(pixX);
    float locY=(dynamic_cast<const PixelGeomDetUnit*>(det))->specificTopology().localY(pixY);
    LocalPoint locXY(locX,locY);
    return locXY;
  }
 
  // pixel size fix
  /*
    std::pair<bool,bool> DeepCoreNtuplizer::pixel2Local(int pixX, int pixY, const GeomDet* det){
    bool pix_x_size=(dynamic_cast<const PixelGeomDetUnit*>(det))->specificTopology().localX(pixX);
    bool pix_y_size=(dynamic_cast<const PixelGeomDetUnit*>(det))->specificTopology().localY(pixY);
    std::pair<bool,bool> pix_xy_size(pix_x_size,pix_y_size);
    return pix_xy_size;
  }
  */
  int DeepCoreNtuplizer::pixelFlipper(const GeomDet* det){
    int out =1;
    LocalVector locZdir(0,0,1);
    GlobalVector globZdir  = det->specificSurface().toGlobal(locZdir);
    GlobalPoint globDetCenter = det->position();
    float direction = globZdir.x()*globDetCenter.x()+ globZdir.y()*globDetCenter.y()+ globZdir.z()*globDetCenter.z();
    //float direction = globZdir.dot(globDetCenter);
    if(direction<0) out =-1;
    // out=1;
    return out;
   }



  void DeepCoreNtuplizer::fillPixelMatrix(const SiPixelCluster & cluster, int layer, auto inter, const GeomDet* det){

    int flip = pixelFlipper(det); // 1=not flip, -1=flip

    for(int i=0; i<cluster.size();i++){
      SiPixelCluster::Pixel pix = cluster.pixel(i);
      std::pair<int,int> pixInter = local2Pixel(inter.x(),inter.y(),det);
      int nx = pix.x-pixInter.first;
      int ny = pix.y-pixInter.second;
      nx=flip*nx;

      if(abs(nx)<jetDimX/2 && abs(ny)<jetDimY/2){
        nx = nx+jetDimX/2;
        ny = ny+jetDimY/2;
        if(nx<jetDimX && ny<jetDimY && layer-1< Nlayer && layer-1>=0 && nx>=0 && ny>=0) clusterMeas[nx][ny][layer-1] += (pix.adc)/(float)(14000);
      }
    }

  }

  void DeepCoreNtuplizer::fillTrackInfo(const reco::Candidate& jet, const auto & stVector, auto jti, auto jetDir, auto jVert, const GeomDet* det, const PixelClusterParameterEstimator* cpe,  std::vector<PSimHit> goodSimHit, const TrackerTopology* const tTopo){

    bool oneHitInfo = false; //track par given only where prob=1

    struct trkInfoObj
    {
      int prob;
      double dist;
      double xpos;
      double ypos;
      double xangle;
      double yangle;
      int zero_flag; //useful for CNN training only: 0 if x,y,eta,phi==0
      double one_over_pt;
      trkInfoObj(int pp, double dd, double xx, double yy, double tx, double ty, int zf, double ptInv) : //, double jeta, double jpt) : // double xd, double yd, double n) :
        prob(pp),
        dist(dd),
        xpos(xx),
        ypos(yy),
        xangle(tx),
        yangle(ty),
        zero_flag(zf),
        one_over_pt(ptInv) {}

      bool operator < (const trkInfoObj& trk) const
      {
        return (dist < trk.dist);
      }
    };

   std::vector<SimTrack> goodSimTrk;

    int trk =0;
    for(uint j=0; j<stVector->size(); j++){ //matched tracks selection and vertex assigment
      for(std::vector<PSimHit>::const_iterator it=goodSimHit.begin(); it!=goodSimHit.end(); ++it) {
        SimTrack st = stVector->at(j);
        if(st.trackId()==(*it).trackId()) {
          if(st.momentum().Pt()<1 or st.momentum().Pt()>100000) continue; //higher cut is to avoid problematic tracks
          goodSimTrk.push_back(st);
        }
      }
    }
    jet_Ntrack=goodSimTrk.size();
    for(uint j=0; j<goodSimTrk.size(); j++){
      SimTrack st = goodSimTrk.at(j);
      if(j<100){
        if(trk<Ntrack) track_pt[j] = st.momentum().Pt();
        if(trk<Ntrack) track_pz[j] = st.momentum().Pz();
        if(trk<Ntrack) track_eta[j] = st.momentum().Eta();
        if(trk<Ntrack) track_phi[j] = st.momentum().Phi();
      }
    }

    for(int x=0; x<jetDimX; x++){
      for(int y=0; y<jetDimY; y++){
        std::vector<trkInfoObj> tracksInfo;
        for(uint j=0; j<goodSimTrk.size(); j++){
          SimTrack st = goodSimTrk.at(j);
          GlobalVector trkMom(st.trackerSurfaceMomentum().x(),st.trackerSurfaceMomentum().y(), st.trackerSurfaceMomentum().z());
          int flip = pixelFlipper(det); // 1=not flip, -1=flip

          PSimHit theSimHit;
          for(uint sh=0; sh<goodSimHit.size(); sh++){
            if(goodSimHit[sh].trackId()==st.trackId() && (int)goodSimHit[sh].detUnitId()==(int)det->geographicalId()){
              theSimHit = goodSimHit[sh];
            }
          }

          LocalPoint localTrkInter(flip*(theSimHit.localPosition()).x(),(theSimHit.localPosition()).y(),(theSimHit.localPosition()).z());
          LocalPoint Jinter(flip*jti.x(),jti.y(),jti.z());

          std::pair<int,int> pixJInter = local2Pixel(Jinter.x(),Jinter.y(),det);
          std::pair<int,int> pixTrkInter = local2Pixel(localTrkInter.x(),localTrkInter.y(),det); //ADDED PITCH

          int pixX = (pixTrkInter.first-pixJInter.first);
          int pixY = (pixTrkInter.second-pixJInter.second);

          pixX = pixX+jetDimX/2;
          pixY = pixY+jetDimY/2;

          double info[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

          if (x==pixX && y==pixY) {
            info[0]= 1;
          }
          else info[0] = 0;
          
          LocalPoint pix2loc = pixel2Local(x+pixJInter.first-jetDimX/2,y+pixJInter.second-jetDimY/2,det);
          double distX =  localTrkInter.x()-pix2loc.x()-0.5*pitchX;//ADDED PITCH
          double distY =  localTrkInter.y()-pix2loc.y()-0.5*pitchY;
         // debug
          std::cout << " dx = " << distX << " Track crossing ( " << localTrkInter.x() << ") - pixel location ( " << pix2loc.x() << " ) - 50 mu" << std::endl;
          std::cout << " dy = " << distY << " Track crossing ( " << localTrkInter.y() << ") - pixel location ( " << pix2loc.y() << " ) - 75 mu" << std::endl;

           if(distX<0.000001 && distX> -0.000001) {
             distX = 0.0;
           }
           if(distY<0.000001 && distY> - 0.000001){
             distY = 0.0;
           }
          info[1] = sqrt(distX*distX+distY*distY);

          if(fabs(distX)<distThr*pitchX && fabs(distY)<distThr*pitchY){

            info[2] = distX;
            info[3] = distY;

            info[4] = st.momentum().Eta()-jetDir.eta();
            info[5] = deltaPhi(st.momentum().Phi(),jetDir.phi());
            info[7] = 1/st.momentum().Pt();
          }
          else{
          if(info[0]== 1) {
             edm::LogError("EventMonitor|seeding") << "prob=1 outside of radious of " << distThr << "pixel from track intersection, probable bug";
          }
          info[2] = 99999999;
          info[3] = 99999999;
          }
          
          if(info[0]==0 && info[2]==0.0 && info[3]==0.0 && info[4]==0.0 && info[5]==0.0) {
            info[6] = 0;
          }
          else if(info[2] > 9999999 && info[3] > 9999999) {
            info[6] = 0;
          }
          else{
             info[6]=1;
           }

          tracksInfo.push_back(trkInfoObj(info[0],info[1],info[2],info[3],info[4],info[5],info[6],info[7]));

        }//good tracks

        std::sort(tracksInfo.begin(), tracksInfo.end());

        int trkLim=Nover;
        if (tracksInfo.size()<3) trkLim = tracksInfo.size();

        for(trk=0; trk<trkLim; trk++){
          if(x<jetDimX && y<jetDimY && trk<Nover && x>=0 && y>=0) {
            if(tracksInfo.at(trk).xpos!=99999999) trackPar[x][y][trk][0]=100*tracksInfo.at(trk).xpos;
            else trackPar[x][y][trk][0]=0.0;
            if(tracksInfo.at(trk).ypos!=99999999) trackPar[x][y][trk][1]=100*tracksInfo.at(trk).ypos;
            else trackPar[x][y][trk][1]=0.0;
            trackPar[x][y][trk][2]=100*tracksInfo.at(trk).xangle;
            trackPar[x][y][trk][3]=100*tracksInfo.at(trk).yangle;
            trackPar[x][y][trk][5]=tracksInfo.at(trk).zero_flag; //NOFLAG
            trackPar[x][y][trk][4]=tracksInfo.at(trk).one_over_pt;
            trackProb[x][y][trk]=tracksInfo.at(trk).prob;

            if(oneHitInfo){
              if(trackProb[x][y][trk]==0) {
                for(int pp=0; pp<Npar+1; pp++){ //NOFLAG
                  trackPar[x][y][trk][pp] = 0.0;
                }
              }
            }

            if(tracksInfo.at(trk).zero_flag!=0){
              edm::LogInfo("EventMonitor|seeding") << "Filling track info with zero-flag!=0, track=" << trk << ", prob="<<tracksInfo.at(trk).prob<< ", x bin="<< x << ", y bin=" << y<< ", "<<  "x position=" << tracksInfo.at(trk).xpos << " y position=" << tracksInfo.at(trk).ypos << std::endl;
            }
         }
        }//trk (filling)
      tracksInfo.clear();
    } // y
    }//x
    goodSimTrk.clear();
  }


  const GeomDet* DeepCoreNtuplizer::DetectorSelector(int llay, const reco::Candidate& jet, GlobalVector jetDir, const reco::Vertex& jetVertex, const GlobalTrackingGeometry* geometry_, const TrackerTopology* const tTopo, const PixelClusterParameterEstimator* parEst, const auto & simtracksVector){

    struct distCompare {
    bool operator()(std::pair<double,const GeomDet*> x, std::pair<double,const GeomDet*> y) const
    {return x.first < y.first;} //order in ascending distance
    };

    std::set<std::pair<double,const GeomDet*>, distCompare> distDetSet;//distance square from center of det

    LocalPoint jetInter(0,0,0);

    std::vector<PSimHit>::const_iterator shIt = simhits->begin();
    std::vector<PSimHit>::const_iterator shItEC = simhitsEC->begin();

    std::set<const GeomDet*> simhitsDetSet;

    for (; shIt != simhits->end(); shIt++) { //loop on simit to find correspondent det (barrel)
      const GeomDet* det = geometry_->idToDet((*shIt).detUnitId());
//      const GeomDet* det = geometry_.idToDet((*shIt).detUnitId());
      simhitsDetSet.insert(det);
    }

    for (; shItEC != simhitsEC->end(); shItEC++) { //loop on simit to find correspondent det (endCaps)
      const GeomDet* det = geometry_->idToDet((*shItEC).detUnitId());
//      const GeomDet* det = geometry_.idToDet((*shItEC).detUnitId());
      simhitsDetSet.insert(det);
    }

    for (std::set<const GeomDet*>::iterator detIt=simhitsDetSet.begin(); detIt != simhitsDetSet.end(); detIt++) { //loop deset

      const GeomDet* det = *detIt;

        int lay = tTopo->layer(det->geographicalId());
        if(det->geographicalId().subdetId()==PixelSubdetector::PixelEndcap){
           lay=lay+4; //endcap layer counting = 5,6,7
        }
        std::pair<bool, Basic3DVector<float>> interPair = findIntersection(jetDir,(reco::Candidate::Point)jetVertex.position(), det);
        if(interPair.first==false) continue;
        Basic3DVector<float> inter = interPair.second;
        auto localInter = det->specificSurface().toLocal((GlobalPoint)inter);

          if(lay==llay) {
            double dist2 = localInter.x()*localInter.x()+localInter.y()*localInter.y();
            std::pair<double,const GeomDet*> distDet(dist2, det);
            distDetSet.insert(distDet);
          }

    } //detset

    if(distDetSet.size()!=0) edm::LogInfo("EventMonitor|seeding") << "Intersecting detector, layer ="<< llay << "  det=" << distDetSet.begin()->second->gdetIndex() << ", distance =" << distDetSet.begin()->first << std::endl;
    if(distDetSet.size()==0) return (GeomDet*)0;
    else {
      return distDetSet.begin()->second;
    }
  }

std::vector<GlobalVector> DeepCoreNtuplizer::splittedClusterDirections(const reco::Candidate& jet, const GlobalTrackingGeometry* geometry_, const TrackerTopology* const tTopo, auto parEst, const reco::Vertex& jetVertex , int layer){
  std::vector<GlobalVector> clustDirs;

  edmNew::DetSetVector<SiPixelCluster>::const_iterator detIt_int = inputPixelClusters->begin();


  for (; detIt_int != inputPixelClusters->end(); detIt_int++) {

    const edmNew::DetSet<SiPixelCluster>& detset_int = *detIt_int;
    const GeomDet* det_int = geometry_->idToDet(detset_int.id());
//    const GeomDet* det_int = geometry_.idToDet(detset_int.id());
    int lay = tTopo->layer(det_int->geographicalId());
    if(det_int->geographicalId().subdetId()==PixelSubdetector::PixelEndcap) lay=lay+4; //endcap layer counting = 5,6,7
    if(lay != layer) continue; //NB: saved bigclusetr on all the layers

    for (auto cluster = detset_int.begin(); cluster != detset_int.end(); cluster++) {
      const SiPixelCluster& aCluster = *cluster;
      GlobalPoint cPos = det_int->surface().toGlobal(parEst->localParametersV(aCluster,(*geometry_->idToDetUnit(detIt_int->id())))[0].first);
//      GlobalPoint cPos = det_int->surface().toGlobal(parEst->localParametersV(aCluster,(*geometry_.idToDetUnit(detIt_int->id())))[0].first);
      GlobalPoint ppv(jetVertex.position().x(), jetVertex.position().y(), jetVertex.position().z());
      GlobalVector clusterDir = cPos - ppv;
      GlobalVector jetDir(jet.px(), jet.py(), jet.pz());
      if (Geom::deltaR(jetDir, clusterDir) < deltaR_) {
            // check if the cluster has to be splitted
            bool isEndCap =
                (std::abs(cPos.z()) > 30.f);  // FIXME: check detID instead!
            float jetZOverRho = jet.momentum().Z() / jet.momentum().Rho();
            if (isEndCap)
              jetZOverRho = jet.momentum().Rho() / jet.momentum().Z();
            float expSizeY =
                std::sqrt((1.3f*1.3f) + (1.9f*1.9f) * jetZOverRho*jetZOverRho);
            if (expSizeY < 1.f) expSizeY = 1.f;
            float expSizeX = 1.5f;
            if (isEndCap) {
              expSizeX = expSizeY;
              expSizeY = 1.5f;
            }  // in endcap col/rows are switched
            float expCharge =
                std::sqrt(1.08f + jetZOverRho * jetZOverRho) * centralMIPCharge_;
            if (aCluster.charge() > expCharge * chargeFracMin_ && (aCluster.sizeX() > expSizeX + 1 ||  aCluster.sizeY() > expSizeY + 1)) {
              edm::LogInfo("EventMonitor|seeding") << "found mergedcluster with deltaR" << Geom::deltaR(jetDir, clusterDir)<< ", on layer=" <<lay << std::endl;
              clustDirs.push_back(clusterDir);
            }
          }
        }
      }
      return clustDirs;

}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DeepCoreNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("pixelClusters", edm::InputTag("siPixelClustersPreSplitting"));
  desc.add<edm::InputTag>("cores", edm::InputTag("ak4CaloJets"));
  desc.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("simHit", edm::InputTag("TrackerHitsPixelBarrelLowTof"));
  desc.add<edm::InputTag>("simHitEC", edm::InputTag("TrackerHitsPixelEndcapLowTof"));
  desc.add<double>("ptMin", 500);
  desc.add<double>("pMin", 0);
  desc.add<double>("deltaR", 0.1);
  desc.add<double>("centralMIPCharge", 18000.0);
  desc.add<double>("chargeFractionMin", 2);
  desc.add<std::string>("pixelCPE", "PixelCPEGeneric");
  desc.add<bool>("barrelTrain", true);
  desc.add<bool>("endcapTrain", false);
  desc.add<bool>("fullTrain", false);
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepCoreNtuplizer);
