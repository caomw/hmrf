#include <common.h>
#include <utility.h>

twister_base_gen_type mygenerator(42u);

unsigned ComputeNumSubs(std::string srcpath);

int BuildGraph(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	       unsigned numSubs, 
	       ImageType3DChar::Pointer maskPtr);

int BuildDataMap(lemon::SmartGraph & theGraph, 
		 lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		 lemon::SmartGraph::NodeMap<vnl_vector<float> > & tsMap,
		 std::string srcpath,
		 ParStruct & par);

int BuildEdgeMap(lemon::SmartGraph & theGraph, 
		 lemon::SmartGraph::EdgeMap<double> & edgeMap,
		 lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		 lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
		 ParStruct & par);

int InitSubSamples(std::string srcpath,
		   lemon::SmartGraph & theGraph, 
		   lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		   lemon::SmartGraph::NodeMap<unsigned short> & initSubjectMap,
		   ParStruct & par);

int EstimateMu(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	       lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	       lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	       ParStruct & par);

int EstimateKappa(lemon::SmartGraph & theGraph, 
		  ParStruct & par);

double EstimatePriorPar(lemon::SmartGraph & theGraph, 
			lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
			lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
			ParStruct & par);

int Sampling(lemon::SmartGraph & theGraph, 
	     lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	     lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
	     lemon::SmartGraph::EdgeMap<double> & edgeMap,
	     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	     ParStruct & par);

void *SamplingThreads(void * threadArgs);

int CompBetaDrv(lemon::SmartGraph & theGraph, 
		lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		ParStruct & par,
		double & drv1, double & drv2);

int EstimateBeta(lemon::SmartGraph & theGraph, 
		lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		lemon::SmartGraph::EdgeMap<double> & edgeMap,
		 ParStruct & par);

int CompSampleEnergy(lemon::SmartGraph & theGraph, 
		     lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		     lemon::SmartGraph::EdgeMap<double> & edgeMap,
		     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
		     ParStruct & par);

int PlotBeta(lemon::SmartGraph & theGraph, 
	     lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
	     lemon::SmartGraph::EdgeMap<double> & edgeMap,
	     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	     ParStruct & par);

using namespace lemon;
int main(int argc, char* argv[])
{
     ParStruct par;
     unsigned seed = 0;
     unsigned int emIter = 1;

     std::string dataFile, iLabelFile, oLabelFile, sampleFile;
     std::string initsubpath;
     std::string fmriPath, initGrplabel, cumusampledir, rsampledir;

     bool estprior = false;
     bool initsame = false;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")
	  ("nthreads", po::value<unsigned>(&par.nthreads)->default_value(100),
	   "number of threads to run for sampling")
	  ("burnin,b", po::value<unsigned>(&par.burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("sweepperthread", po::value<unsigned>(&par.sweepPerThread)->default_value(1),
	   "number of scans for burn-in period. ")
	  ("inittemp", po::value<float>(&par.initTemp)->default_value(1),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<float>(&par.finalTemp)->default_value(1),
	   "Final temperature for annealing to stop")
	  ("numSamples,n", po::value<unsigned int>(&par.numSamples)->default_value(20),
	   "Number of Monte Carlo samples. ")
	  ("emiter", po::value<unsigned int>(&emIter)->default_value(30),
	   "Number of EM iteration. ")

	  ("alpha", po::value<float>(&par.alpha)->default_value(0.7),
	   "connection between group lebel and individual subject level.")

	  ("beta", po::value<float>(&par.beta)->default_value(0.5),
	   "pairwise interation term")

	  ("gamma", po::value<float>(&par.gamma)->default_value(1),
	   " weights of lilelihood term.")

	  ("alpha0", po::value<float>(&par.alpha0)->default_value(0.7),
	   "Mean of prior on alpha.")

	  ("sigma2", po::value<float>(&par.sigma2)->default_value(1),
	   "variance of prior on alpha")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")

	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	   "Number of clusters. Default is 6.")

	  ("estprior", po::bool_switch(&estprior), 
	   "whether to estimate alpha and beta, the prior parameter. Default no.")
	  ("initsame", po::bool_switch(&initsame), 
	   "whether to init subject same to group label. Default is no.")
	  ("weightbetadata", po::bool_switch(&par.weightbetadata), 
	   "whether to weight beta by voxel correlation. Default is no.")
	  ("weightbetadist", po::bool_switch(&par.weightbetadist), 
	   "whether to weight beta by voxel's spatial distance. Default is no.")
	   ("initgrouplabel,i", po::value<std::string>(&initGrplabel)->default_value("grouplabel.nii"), 
	    "Initial group level label map (Intensity value 1-K. Also used as mask file.")
	   ("initsubpath,t", po::value<std::string>(&initsubpath)->default_value("./subpath"), 
	    "Initial subject label maps path")
	  ("fmripath,f", po::value<std::string>(&fmriPath)->default_value("."), 
	   "noised image file")

	  ("cumusampledir", po::value<std::string>(&cumusampledir)->default_value("./cumusamples"), 
	   "Monte Carlo samples file dir.")
	  ("rsampledir", po::value<std::string>(&rsampledir)->default_value("./rsamples"), 
	   "running samples file dir. Also can be the ultimate output dir.")

	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: gropumrf [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));

     // read in initial grp label. Also used as mask. 1-based label converted to
     // 0-based label.
     ReaderType3DChar::Pointer grpReader = ReaderType3DChar::New();
     grpReader->SetFileName(initGrplabel);
     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();
     myfilter->InPlaceOff();
     myfilter->SetInput( grpReader->GetOutput());
     myfilter->SetConstant(-1);
     myfilter->Update();

     // output label init'd.
     ImageType3DChar::Pointer labelPtr = myfilter->GetOutput();

     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = grpReader->GetOutput();
     par.maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     par.cumusampledir = cumusampledir;
     par.rsampledir = rsampledir;
     par.numSubs = ComputeNumSubs(fmriPath);
     printf("number of subjects: %i.\n", par.numSubs);
     par.temperature = 1;

     // init parameters so BuildDataMap knows what to do.
     par.vmm.sub.resize(par.numSubs);

     // Construct graph.
     lemon::SmartGraph theGraph;
     lemon::SmartGraph::NodeMap<SuperCoordType> coordMap(theGraph); // node --> coordinates.
     BuildGraph(theGraph, coordMap, par.numSubs, maskPtr);

     // Build Data map.
     lemon::SmartGraph::NodeMap<vnl_vector<float>> tsMap(theGraph); // node --> time series.
     BuildDataMap(theGraph, coordMap, tsMap, fmriPath, par);

     // Define edge map.
     lemon::SmartGraph::EdgeMap<double> edgeMap(theGraph);
     BuildEdgeMap(theGraph, edgeMap, coordMap, tsMap, par);

     // At this point, we know time series length so do the remaining
     // initializaion of par.vmm.
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  par.vmm.sub[subIdx].comp.resize(par.numClusters);
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       par.vmm.sub[subIdx].comp[clsIdx].mu.set_size(par.tsLength);
	       par.vmm.sub[subIdx].comp[clsIdx].mu = 0;
	       par.vmm.sub[subIdx].comp[clsIdx].meanNorm = 0;
	       par.vmm.sub[subIdx].comp[clsIdx].kappa = 0;
	       par.vmm.sub[subIdx].comp[clsIdx].numPts = 0;
	       par.vmm.sub[subIdx].comp[clsIdx].prop = 0;
	  }
     }

     PrintPar(0, par);

     // define a subject initial label map used for initialization of subject
     // label map.
     lemon::SmartGraph::NodeMap<unsigned short> initSubjectMap(theGraph);     
     if (!initsame) {
	  InitSubSamples(initsubpath, theGraph, coordMap, initSubjectMap, par);
     }

     // define cumulative sample map. Each node is assigned a vector. The k'th
     // elements of the vector is the number of samples that falls into cluster
     // k.
     lemon::SmartGraph::NodeMap<std::vector<unsigned short> > cumuSampleMap(theGraph);
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
     	  cumuSampleMap[nodeIt].resize(par.numClusters, 0);
     }

     // init the cumuSampleMap just to estimate parameters. The k'th bit is set
     // to the number of samples M, since I suppose the M samples are init'd
     // same to the initial group or subject label map.
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  if (coordMap[nodeIt].subid < par.numSubs && (!initsame) ) {
	       if (initSubjectMap[nodeIt] < par.numClusters) {
		    cumuSampleMap[nodeIt].at(initSubjectMap[nodeIt]) = par.numSamples;
	       }
	       else {
		    cumuSampleMap[nodeIt].at(0) = par.numSamples;
	       }
	  }
	  else {
	       if (labelPtr->GetPixel(coordMap[nodeIt].idx) < par.numClusters) {
		    cumuSampleMap[nodeIt].at(labelPtr->GetPixel(coordMap[nodeIt].idx)) = par.numSamples;
	       }
	       else {
		    cumuSampleMap[nodeIt].at(0) = par.numSamples;
	       }
	  }
     }
     if (par.verbose >= 2) {
	  SaveCumuSamples(theGraph, coordMap, cumuSampleMap, par);
     }
     
     // define a running sample map and init it with the input group
     // labels. This is used for sampling only. And when the whoel scan is done,
     // this map is saved into the cumuSampleMap.
     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > rSampleMap(theGraph);
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  rSampleMap[nodeIt].resize(par.numClusters, 0);
	  if (coordMap[nodeIt].subid < par.numSubs && (!initsame) ) {
	       if (initSubjectMap[nodeIt] < par.numClusters) {
		    rSampleMap[nodeIt].set(initSubjectMap[nodeIt], true);
	       }
	       else {
		    rSampleMap[nodeIt].set(0, true);
	       }
	  }
	  else {
	       if (labelPtr->GetPixel(coordMap[nodeIt].idx) < par.numClusters) {
		    rSampleMap[nodeIt].set(labelPtr->GetPixel(coordMap[nodeIt].idx), true);
	       }
	       else {
		    rSampleMap[nodeIt].set(0, true);
	       }
	  }
     }
     if (par.verbose >= 2) {
	  SaveRunningSamples(theGraph, coordMap, rSampleMap, par);
     }

     // as initial step, esimate mu and kappa.
     EstimateMu(theGraph, coordMap, cumuSampleMap, tsMap, par);
     EstimateKappa(theGraph, par);
     PrintPar(2, par);

     // EM starts here.
     for (unsigned short emIterIdx = 0; emIterIdx < emIter; emIterIdx ++) {     
     	  printf("EM iteration %i begin:\n", emIterIdx + 1);
	  par.temperature = par.initTemp * pow( (par.finalTemp / par.initTemp), float(emIterIdx) / float(emIter) );
	  Sampling(theGraph, coordMap, cumuSampleMap, rSampleMap, edgeMap, tsMap, par);
	  if (par.verbose >= 2) {
	       SaveCumuSamples(theGraph, coordMap, cumuSampleMap, par);

	  }


     	  // estimate vMF parameters mu, kappa.
     	  printf("EM iteration %i, parameter estimation begin. \n", emIterIdx + 1);
	  // estimate prior parameter, now beta.
	  if(estprior) {
	       CompSampleEnergy(theGraph, coordMap, rSampleMap, edgeMap, tsMap, par); 
	       EstimateBeta(theGraph, coordMap, rSampleMap, edgeMap, par);
	       // PlotBeta(theGraph, coordMap, rSampleMap, edgeMap, tsMap, par);
	  }

	  EstimateMu(theGraph, coordMap, cumuSampleMap, tsMap, par);
	  EstimateKappa(theGraph, par);
	  PrintPar(2, par);
	  CompSampleEnergy(theGraph, coordMap, rSampleMap, edgeMap, tsMap, par); 
     } // emIterIdx

     SaveCumuSamples(theGraph, coordMap, cumuSampleMap, par);
     SaveRunningSamples(theGraph, coordMap, rSampleMap, par);

     return 0;
}     

unsigned ComputeNumSubs(std::string srcpath)
{
     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(srcpath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );

     // Get number of subjects.
     unsigned numSubs  = sortedEntries.size();
     return numSubs;
}


int BuildGraph(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	       unsigned numSubs, 
	       ImageType3DChar::Pointer maskPtr)
{
     //mask
     IteratorType3DCharIdx maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     ImageType3DChar::IndexType neiIdx;
     
     // map xyz to graph node id.
     ImageType4Int::IndexType nodeMapIdx;
     nodeMapIdx.Fill(0);
     ImageType3DChar::SizeType maskSize =  maskPtr->GetLargestPossibleRegion().GetSize();
     ImageType4Int::SizeType nodeMapSize;
     nodeMapSize[0] = maskSize[0];
     nodeMapSize[1] = maskSize[1];
     nodeMapSize[2] = maskSize[2];
     nodeMapSize[3] = numSubs + 1;

     ImageType4Int::RegionType nodeMapRegion;
     nodeMapRegion.SetSize(nodeMapSize);
     nodeMapRegion.SetIndex(nodeMapIdx);
     ImageType4Int::Pointer nodeMapPtr = ImageType4Int::New();
     nodeMapPtr->SetRegions(nodeMapRegion);
     nodeMapPtr->Allocate();
     nodeMapPtr->FillBuffer(0);

     
     // Add nodes.

     lemon::SmartGraph::Node curNode, neiNode;

     // fill in subjects nodes.
     for (unsigned short subIdx = 0; subIdx < numSubs; subIdx ++) {
	  for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    curNode = theGraph.addNode();
		    coordMap[curNode].idx = maskIt.GetIndex();
		    coordMap[curNode].subid = subIdx;
		    nodeMapIdx[0] = coordMap[curNode].idx[0];
		    nodeMapIdx[1] = coordMap[curNode].idx[1];
		    nodeMapIdx[2] = coordMap[curNode].idx[2];
		    nodeMapIdx[3] = subIdx;
		    nodeMapPtr->SetPixel(nodeMapIdx, theGraph.id(curNode) );
	       }
	  }
     }

     // Fill in group level nodes.
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
     	  if (maskIt.Get() > 0) {
     	       curNode = theGraph.addNode();
	       coordMap[curNode].idx = maskIt.GetIndex();
	       coordMap[curNode].subid = numSubs;
	       
	       nodeMapIdx[0] = coordMap[curNode].idx[0];
	       nodeMapIdx[1] = coordMap[curNode].idx[1];
	       nodeMapIdx[2] = coordMap[curNode].idx[2];
	       nodeMapIdx[3] = numSubs;
	       nodeMapPtr->SetPixel(nodeMapIdx, theGraph.id(curNode) );
	  }
     }

     printf("BuildGraph(): number of Nodes: %i\n", theGraph.maxNodeId()+1 );

     // Add edges. 

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);     
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType maskNeiIt(radius, maskPtr, maskPtr->GetLargestPossibleRegion() );
     maskNeiIt.OverrideBoundaryCondition(&constCondition);
     
     // xplus, xminus, yplus, yminus, zplus, zminus
     // std::array<unsigned int, 6 > neiIdxSet = {{14, 12, 16, 10, 22, 4}}; 
     std::array<unsigned int, 26> neiIdxSet = {{ 0,1,2,3,4,5,6,7,8,9,10,11,12,//no 13
						 14,15,16,17,18,19,20,21,22,23,24,25,26 }};
     // std::array<unsigned int, 18 > neiIdxSet = {{1,3,4,5,7,9,10,11,12,
     // 						14,15,16,17,19,21,22,23,25}};
     ImageType3DChar::IndexType maskIdx;
     int curNodeId = 0, neiNodeId = 0;
     // std::array<short, 6>::const_iterator neiIdxIt;
     auto neiIdxIt = neiIdxSet.begin();

     // edge within grp.
     for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	  if (maskNeiIt.GetCenterPixel() > 0) {
	       maskIdx = maskNeiIt.GetIndex();
	       nodeMapIdx[0] = maskIdx[0];
	       nodeMapIdx[1] = maskIdx[1];
	       nodeMapIdx[2] = maskIdx[2];
	       nodeMapIdx[3] = numSubs;
	       curNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
	       // curNode is group node.
	       curNode = theGraph.nodeFromId( curNodeId );

	       for (neiIdxIt=neiIdxSet.begin(); neiIdxIt < neiIdxSet.end(); neiIdxIt++) {
		    if (maskNeiIt.GetPixel(*neiIdxIt) > 0) {
			 neiIdx = maskNeiIt.GetIndex(*neiIdxIt);
			 nodeMapIdx[0] = neiIdx[0];
			 nodeMapIdx[1] = neiIdx[1];
			 nodeMapIdx[2] = neiIdx[2];
			 nodeMapIdx[3] = numSubs;
			 neiNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
			 // make sure each edge is added only once.
			 if (neiNodeId > curNodeId) {
			      neiNode = theGraph.nodeFromId( neiNodeId );
			      theGraph.addEdge(curNode, neiNode);
			 } // if
		    } // maskIt
	       } // neiIt
	  } // maskNeiIt > 0
     } // maskNeiIt
     printf("BuildGraph(): After adding edges within group, number of edges: %i\n", theGraph.maxEdgeId() + 1);

     // edge between grp and subjects.
     for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	  if (maskNeiIt.GetCenterPixel() > 0) {
	       maskIdx = maskNeiIt.GetIndex();
	       nodeMapIdx[0] = maskIdx[0];
	       nodeMapIdx[1] = maskIdx[1];
	       nodeMapIdx[2] = maskIdx[2];

	       // get current node in group.
	       nodeMapIdx[3] = numSubs;
	       curNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
	       curNode = theGraph.nodeFromId( curNodeId );
	       
	       // get subject node and connect them
	       for (unsigned short subIdx = 0; subIdx < numSubs; subIdx ++) {
		    nodeMapIdx[3] = subIdx;
		    neiNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
		    neiNode = theGraph.nodeFromId( neiNodeId );
		    theGraph.addEdge(curNode, neiNode);
	       } // subIdx
	  } // maskNeiIt > 0
     } // MaskNeiIt

     printf("BuildGraph(): After adding edges btw grp and subs, number of edges: %i\n", theGraph.maxEdgeId() + 1);

     // edge within each subject.
     for (unsigned short subIdx = 0; subIdx < numSubs; subIdx ++) {
	  nodeMapIdx[3] = subIdx;

	  for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	       if (maskNeiIt.GetCenterPixel() > 0) {
		    maskIdx = maskNeiIt.GetIndex();
		    nodeMapIdx[0] = maskIdx[0];
		    nodeMapIdx[1] = maskIdx[1];
		    nodeMapIdx[2] = maskIdx[2];
		    curNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
		    // curNode is in subject map.
		    curNode = theGraph.nodeFromId( curNodeId );

		    for (neiIdxIt=neiIdxSet.begin(); neiIdxIt < neiIdxSet.end(); neiIdxIt++) {
			 if (maskNeiIt.GetPixel(*neiIdxIt) > 0) {
			      neiIdx = maskNeiIt.GetIndex(*neiIdxIt);
			      nodeMapIdx[0] = neiIdx[0];
			      nodeMapIdx[1] = neiIdx[1];
			      nodeMapIdx[2] = neiIdx[2];
			      neiNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
			      // make sure each edge is added only once.
			      if (neiNodeId > curNodeId) {
				   neiNode = theGraph.nodeFromId( neiNodeId );
				   theGraph.addEdge(curNode, neiNode);
			      } // if
			 } // maskNeiIt
		    } // neiIdxIt
	       } // maskNeiIt >0 
	  }
     } // subIdx

     printf("BuildGraph(): number of edges: %i\n", theGraph.maxEdgeId() + 1);

     return 0;
}


int BuildEdgeMap(lemon::SmartGraph & theGraph, 
		 lemon::SmartGraph::EdgeMap<double> & edgeMap,
		 lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		 lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
		 ParStruct & par)
{
     unsigned weightIdx = 0;
     std::array<double, 3> distWeight = {{BETAWEIGHT0, BETAWEIGHT1, BETAWEIGHT2}};
     for (SmartGraph::EdgeIt edgeIt(theGraph); edgeIt != INVALID; ++ edgeIt) {
	  if (coordMap[theGraph.u(edgeIt)].subid < par.numSubs && 
	      coordMap[theGraph.v(edgeIt)].subid < par.numSubs) {
	       // within subjects.
	       if (par.weightbetadata) {
		    edgeMap[edgeIt] = 0.5 * par.beta * pow(inner_product(tsMap[theGraph.u(edgeIt)], tsMap[theGraph.v(edgeIt)]) + 1, 2);
	       }
	       else if (par.weightbetadist){
		    weightIdx = abs(coordMap[theGraph.u(edgeIt)].idx[0] - coordMap[theGraph.v(edgeIt)].idx[0]) 
			 + abs(coordMap[theGraph.u(edgeIt)].idx[1] - coordMap[theGraph.v(edgeIt)].idx[1]) 
			 + abs(coordMap[theGraph.u(edgeIt)].idx[2] - coordMap[theGraph.v(edgeIt)].idx[2])
			 -1;
		    edgeMap[edgeIt] = distWeight[weightIdx] * par.beta;
	       }
	       else {
	       edgeMap[edgeIt] = par.beta;
	       }
	  }
	  else if (coordMap[theGraph.u(edgeIt)].subid == par.numSubs && 
		   coordMap[theGraph.v(edgeIt)].subid == par.numSubs) {
	       // within group.
	       edgeMap[edgeIt] = par.beta;
	  }
	  else {
	       // must be between group and subjects.
	       edgeMap[edgeIt] = par.alpha;
	  }
     }
     return 0;

}
int BuildDataMap(lemon::SmartGraph & theGraph, 
		 lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		 lemon::SmartGraph::NodeMap<vnl_vector<float> > & tsMap,
		 std::string srcpath,
		 ParStruct & par)
{
     // some of the code below refers to boost filesystem tutorial tut04.cpp.
     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(srcpath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );

     // Get number of subjects.
     unsigned numSubs  = sortedEntries.size();

     // init file reader pointers.
     std::vector<ReaderType4DFloat::Pointer> readerVec(numSubs);
     std::vector<ImageType4DF::Pointer> srcPtrVec(numSubs);
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  readerVec[subIdx] = ReaderType4DFloat::New();
	  readerVec[subIdx]->SetFileName( sortedEntries[subIdx].string() );
	  readerVec[subIdx]->Update();
	  srcPtrVec[subIdx]= readerVec[subIdx]->GetOutput();

	  // also save file name into par.
	  boost::filesystem::path thispath = sortedEntries[subIdx];
	  if(thispath.extension().string().compare(".gz") == 0 ) {
	       // must be a .nii.gz file. Remove .gz suffix so the filename is
	       // abc.nii.
	       thispath = thispath.stem();
	  }
	  par.vmm.sub[subIdx].name.assign(thispath.filename().stem().string() );
     }

     ImageType4DF::SizeType srcSize = readerVec[0]->GetOutput()->GetLargestPossibleRegion().GetSize();
     par.tsLength = srcSize[3];

     printf("BuildDataMap(): Time series length: %i.\n", par.tsLength);

     SuperCoordType superCoord;
     ImageType4DF::IndexType srcIdx;
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
     	  superCoord = coordMap[nodeIt];
	  // Only when not in group map, do this step, because only subject map
	  // has data term, while group nodes do not.
	  if (superCoord.subid < numSubs) {
	       tsMap[nodeIt].set_size(par.tsLength);
	       srcIdx[0] = superCoord.idx[0];
	       srcIdx[1] = superCoord.idx[1];
	       srcIdx[2] = superCoord.idx[2];
	       for (srcIdx[3] = 0; srcIdx[3] < par.tsLength; srcIdx[3] ++) {
		    tsMap[nodeIt][srcIdx[3]] = srcPtrVec[superCoord.subid]->GetPixel(srcIdx);
	       }
	  }
     }

     return 0;
}


/* Read initial subject label map from files. */
int InitSubSamples(std::string srcpath,
		   lemon::SmartGraph & theGraph, 
		   lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		   lemon::SmartGraph::NodeMap<unsigned short> & initSubjectMap,
		   ParStruct & par)
{
     // some of the code below refers to boost filesystem tutorial tut04.cpp.
     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(srcpath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );
     if (sortedEntries.size() != par.numSubs) {
	  std::cout << "number of files in " << srcpath << "is not equal to number of files in fmri path.\n";
	  exit(1);
     }

     // init file reader pointers.
     std::vector<ReaderType3DChar::Pointer> readerVec(par.numSubs);
     std::vector<ImageType3DChar::Pointer> srcPtrVec(par.numSubs);
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  readerVec[subIdx] = ReaderType3DChar::New();
	  readerVec[subIdx]->SetFileName( sortedEntries[subIdx].string() );
	  readerVec[subIdx]->Update();
	  srcPtrVec[subIdx]= readerVec[subIdx]->GetOutput();
     }

     SuperCoordType superCoord;
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
     	  superCoord = coordMap[nodeIt];
	  // when not in group map, do this step, since this function only read
	  // in subject's initial label map.
	  if (superCoord.subid < par.numSubs) {
	       initSubjectMap[nodeIt] = srcPtrVec[superCoord.subid]->GetPixel(superCoord.idx) - 1;
	  }
	  else {
	       initSubjectMap[nodeIt] = -1;
	  }
     }
     return 0;
}


int EstimateMu(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	       lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	       lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	       ParStruct & par)
{
     // reset all mu and numPts to zero.
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       par.vmm.sub[subIdx].comp[clsIdx].mu = 0;
	       par.vmm.sub[subIdx].comp[clsIdx].numPts = 0;
	  }
     }

     SuperCoordType superCoord;
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  superCoord = coordMap[nodeIt];
	  // if this is a subject node.
	  if (superCoord.subid < par.numSubs) {
	       for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
		    par.vmm.sub[superCoord.subid].comp[clsIdx].mu += tsMap[nodeIt] * cumuSampleMap[nodeIt][clsIdx];
		    par.vmm.sub[superCoord.subid].comp[clsIdx].numPts += cumuSampleMap[nodeIt][clsIdx];
	       }
	  } // subid < par.numSubs
     }

     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {     
	       if(par.vmm.sub[subIdx].comp[clsIdx].numPts > 0) {
		    // printf("numPts: %d. ", par.vmm.sub[subIdx].comp[clsIdx].numPts);
		    // printVnlVector(par.vmm.sub[subIdx].comp[clsIdx].mu, 5);
		    par.vmm.sub[subIdx].comp[clsIdx].meanNorm = par.vmm.sub[subIdx].comp[clsIdx].mu.two_norm() / par.vmm.sub[subIdx].comp[clsIdx].numPts;
		    par.vmm.sub[subIdx].comp[clsIdx].mu.normalize();
		    // printVnlVector(par.vmm.sub[subIdx].comp[clsIdx].mu, 5);
	       }
	       else {
		    par.vmm.sub[subIdx].comp[clsIdx].meanNorm = 0;
	       }
	  }
     }
	       
     printf("EstimateMu(). Done.\n");
     return 0;
}

int EstimateKappa(lemon::SmartGraph & theGraph, 
		  ParStruct & par)
{
     float kappa = 0, kappa_new = 0;
     double Ap = 0;
     float RBar = 0;
     float Dim = par.tsLength;

     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       if (par.vmm.sub[subIdx].comp[clsIdx].numPts > 0){

		    RBar = par.vmm.sub[subIdx].comp[clsIdx].meanNorm;
		    kappa_new = RBar * (Dim - RBar * RBar) / (1 - RBar * RBar);
		    unsigned iter = 0;
		    do {
			 iter ++;
			 kappa = kappa_new;
			 Ap = exp(logBesselI(Dim/2, kappa) - logBesselI(Dim/2 - 1, kappa));
			 kappa_new = kappa - ( (Ap - RBar) / (1 - Ap * Ap - (Dim - 1) * Ap / kappa)  );
			 if (par.verbose >= 3) {
			      printf("    sub[%i] cls[%i] kappa: %3.1f -> %3.1f\n", subIdx + 1, clsIdx + 1, kappa, kappa_new);
			 }
		    }
		    while(vnl_math_abs(kappa_new - kappa) > 0.01 * kappa && iter < 5);

		    // the temperature parameter change the definition of the
		    // posterior distribution (and also the vmf likelihood
		    // function. This is equivalnet to divide the kappa paramter
		    // by T.  The kappa_new we estimate is actually kappa/T. So
		    // we need kappa_new * T to recover the original kappa
		    // paramter.
		    par.vmm.sub[subIdx].comp[clsIdx].kappa = kappa_new * par.temperature;
	       } // numPts > 0
	       else {
		    par.vmm.sub[subIdx].comp[clsIdx].kappa = 0;
	       }

	  } // clsIdx
     } // subIdx

     printf("EstimateKappa(). Done.\n");
     return 0;
}

double EstimatePriorPar(lemon::SmartGraph & theGraph, 
			lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
			lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
			ParStruct & par)
{

     return 0;
}

int Sampling(lemon::SmartGraph & theGraph, 
	     lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	     lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
	     lemon::SmartGraph::EdgeMap<double> & edgeMap,
	     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	     ParStruct & par)
{
     unsigned taskid = 0;
     pthread_t thread_ids[par.nthreads];
     ThreadArgs threadArgs[par.nthreads];

     // reset cumuSampleMap to zero.
     printf("reset cumuSampleMap begin...");
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  cumuSampleMap[nodeIt].assign(par.numClusters, 0);
     }
     printf("done.\n");

     // Compute the number of nodes that each thread computes.
     unsigned totalNumNodes = theGraph.maxNodeId()+1;     
     unsigned numNodesPerThread = ceil( double(totalNumNodes) / double(par.nthreads) );

     // compute normalization constant of von Mises
     // Fisher. vmfLogConst[i] is the log (c_d (kappa) ) for the i'th
     // clusters. See page 1350 of "Clustering on the Unit-Sphere
     // using von Mises Fisher..."
     std::vector< vnl_vector<double> >vmfLogConst(par.numSubs);
     float myD = par.tsLength;
     double const Pi = 4 * atan(1);
     unsigned subIdx = 0, clsIdx = 0;
     for (subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  vmfLogConst[subIdx].set_size(par.numClusters);
	  for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       if (par.vmm.sub[subIdx].comp[clsIdx].kappa > 1) {
		    vmfLogConst[subIdx][clsIdx] = (myD/2 - 1) * log (par.vmm.sub[subIdx].comp[clsIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, par.vmm.sub[subIdx].comp[clsIdx].kappa);
	       }
	       else {
		    vmfLogConst[subIdx][clsIdx] = 25.5; // the limit when kappa -> 0
	       }
	  } // for clsIdx
     } // for subIdx

     for (unsigned taskid = 0; taskid < par.nthreads; taskid ++) {
	  threadArgs[taskid].taskid = taskid;
	  threadArgs[taskid].theGraphPtr = & theGraph;
	  threadArgs[taskid].coordMapPtr = & coordMap;
	  threadArgs[taskid].edgeMapPtr = & edgeMap;
	  threadArgs[taskid].rSampleMapPtr = & rSampleMap;
	  threadArgs[taskid].tsMapPtr = & tsMap;
	  threadArgs[taskid].parPtr = & par;
	  threadArgs[taskid].numThreads = par.nthreads;
	  threadArgs[taskid].startNodeid = taskid * numNodesPerThread;
	  threadArgs[taskid].vmfLogConstPtr = & vmfLogConst;
	  
	  // make sure the last thread have correct node id to process.
	  if (taskid == par.nthreads - 1) {
	       threadArgs[taskid].endNodeid = totalNumNodes - 1;
	  }
	  else {
	       threadArgs[taskid].endNodeid = (taskid + 1) * numNodesPerThread -1;
	  }
     }

     // sampling starts here.
     printf("Sampling() starts scan...\n");
     for (unsigned scanIdx = 0; scanIdx < par.burnin + par.numSamples; scanIdx ++) {
	  printf("%i, ", scanIdx+1);
	  fflush(stdout);
	  for (taskid = 0; taskid < par.nthreads; taskid ++) {
	       pthread_create(&thread_ids[taskid], NULL, SamplingThreads, (void*) &threadArgs[taskid]);
	  }

	  for (taskid = 0; taskid < par.nthreads; taskid ++) {
	       pthread_join(thread_ids[taskid], NULL);
	  }

	  if (par.verbose >= 3) {
	       CompSampleEnergy(theGraph, coordMap, rSampleMap, edgeMap, tsMap, par); 
	  }
	  
	  // after burnin peroid, save it to correct place.
	  if (scanIdx >= par.burnin) {
	       for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
		    cumuSampleMap[nodeIt][rSampleMap[nodeIt].find_first()] ++;
	       }
	  }

	  if (par.verbose >=3) {
	       SaveRunningSamples(theGraph, coordMap, rSampleMap, par);
	  }

     } // scanIdx
	    
     printf("Sampling done.\n");
     return 0;
}


void *SamplingThreads(void * threadArgs)
{
     ThreadArgs * args = (ThreadArgs *) threadArgs;
     lemon::SmartGraph * theGraphPtr = args->theGraphPtr;
     lemon::SmartGraph::NodeMap<SuperCoordType> * coordMapPtr = args->coordMapPtr;
     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > * rSampleMapPtr = args->rSampleMapPtr;
     lemon::SmartGraph::EdgeMap<double> * edgeMapPtr = args->edgeMapPtr;
     lemon::SmartGraph::NodeMap<vnl_vector<float>> * tsMapPtr = args->tsMapPtr;
     std::vector< vnl_vector<double> > * vmfLogConstPtr = args->vmfLogConstPtr;
     ParStruct * parPtr = args->parPtr;

     // define random generator.
     boost::random::mt19937 gen(std::time(0));
     boost::random::uniform_int_distribution<> roll_die(0, (*parPtr).numClusters - 1);

     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(gen, uni_dist);

     // Sampling starts here.
     unsigned sweepIdx = 0;
     float oldPriorEngy = 0, newPriorEngy = 0, oldDataEngy = 0, newDataEngy = 0, denergy;
     boost::dynamic_bitset<> cand;
     cand.resize((*parPtr).numClusters);
     float p_acpt = 0;
     // cl is current label, and nl is new label. s is subject id.
     unsigned short cl = 0, nl = 0, s = 0;
     lemon::SmartGraph::Node curNode;


     for (sweepIdx = 0 ;sweepIdx < parPtr->sweepPerThread; sweepIdx ++) {
	  for (unsigned nodeId = args->startNodeid; nodeId <= args->endNodeid; nodeId++) {
	       curNode = (*theGraphPtr).nodeFromId(nodeId);
	       // compute old and new prior energy.
	       oldPriorEngy = 0;
	       newPriorEngy = 0;
	       cand.reset();
	       cand[roll_die(gen)] = 1;
	       for (SmartGraph::IncEdgeIt edgeIt(*theGraphPtr, curNode); edgeIt != INVALID; ++ edgeIt) {
	       	    oldPriorEngy += (*edgeMapPtr)[edgeIt] * (double)(
	       		 (*rSampleMapPtr)[(*theGraphPtr).baseNode(edgeIt)] != 
	       		 (*rSampleMapPtr)[(*theGraphPtr).runningNode(edgeIt)]);
	       	    newPriorEngy += (*edgeMapPtr)[edgeIt] * (double)(
	       		 cand !=
	       		 (*rSampleMapPtr)[(*theGraphPtr).runningNode(edgeIt)]);
	       }
	       denergy = newPriorEngy - oldPriorEngy;

	       // if current node is in subject, compute data term. Otherwise,
	       // do not compute data term.
	       if ( (*coordMapPtr)[curNode].subid < (*parPtr).numSubs) {
	       	    cl = (*rSampleMapPtr)[curNode].find_first(); // current label.
	       	    nl = cand.find_first(); // new label.
	       	    s = (*coordMapPtr)[curNode].subid;
	       	    oldDataEngy = - (*vmfLogConstPtr)[s][cl] - (*parPtr).vmm.sub[s].comp[cl].kappa * inner_product((*tsMapPtr)[curNode], (*parPtr).vmm.sub[s].comp[cl].mu);
	       	    newDataEngy = - (*vmfLogConstPtr)[s][nl] - (*parPtr).vmm.sub[s].comp[nl].kappa * inner_product((*tsMapPtr)[curNode], (*parPtr).vmm.sub[s].comp[nl].mu);
	       	    denergy = denergy + parPtr->gamma * (newDataEngy - oldDataEngy);
	       }
	       else {
	       }
	       denergy = denergy / (*parPtr).temperature;

	       if (denergy < 0) {
	       	    (*rSampleMapPtr)[curNode] = cand;
	       }
	       else if (denergy > 0) {
	       	    p_acpt = exp(-denergy);
	       	    if (uni() < p_acpt) {
	       		 (*rSampleMapPtr)[curNode] = cand;
	       	    }
	       } // else if
	       else {
	       	    // candidate = current label. No change.
	       	    (*rSampleMapPtr)[curNode] = cand;
	       }
	  } // curNode
     } // sweepIdx
     return 0;
} 


int CompSampleEnergy(lemon::SmartGraph & theGraph, 
		     lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		     lemon::SmartGraph::EdgeMap<double> & edgeMap,
		     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
		     ParStruct & par)
{
     double samplePriorEng = 0, sampleDataEng = 0;
     boost::dynamic_bitset<> curBit(par.numClusters), runningBit(par.numClusters);
     double eta_cur = 0, M0 = 0;
     unsigned runningLabel = 0;
     vnl_vector<double> eta(par.numClusters, 0);
     double eta_offset = 0;
     unsigned pCount = 0;
     double priorEngOld = 0;
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  curBit = rSampleMap[nodeIt];
	  eta = 0;
	  for (runningLabel = 0; runningLabel < par.numClusters; runningLabel ++) {
	       runningBit.reset(), runningBit[runningLabel] = 1;
	       for (SmartGraph::IncEdgeIt edgeIt(theGraph, nodeIt); edgeIt != INVALID; ++ edgeIt) {
		    eta[runningLabel] = eta[runningLabel] - edgeMap[edgeIt] * (double)(runningBit != rSampleMap[theGraph.runningNode(edgeIt)]) ;
	       } // incEdgeIt

	       //The temperature parameter that change the definition of
	       //the posterior distribution.
	       eta[runningLabel] = eta[runningLabel] / par.temperature;

	       if (runningBit == curBit) {
		    eta_cur = eta[runningLabel];
	       }
	  } // runningBit.
	  eta_offset = eta.max_value();

	  M0 = 0;	  
	  double M0_old = 0;

	  for (runningLabel = 0; runningLabel < par.numClusters; runningLabel ++) {	  
	       M0 += exp (eta[runningLabel] - eta_offset);
	       M0_old += exp(eta[runningLabel]);

	  }

	  samplePriorEng += (eta_cur - eta_offset - log(M0));

	  // if(pCount < 30 && exp(eta_offset)*M0 != M0_old) {
	  if (pCount < 30 && theGraph.id(nodeIt)%1 == 0) {
	       // printf("node = %i, M0 (offset)=%E, M0_old=%E\n", theGraph.id(nodeIt), exp(eta_offset)*M0, M0_old);
	       printf("node = %i, M0 = %E, eta_cur = %f, eta_offset = %f, priorEngy = %E\n", theGraph.id(nodeIt), M0, eta_cur, eta_offset, samplePriorEng);
	       printVnlVector(eta, par.numClusters);
	       pCount ++;
	  }

	  priorEngOld += (eta_cur - log(M0_old) );
     } // nodeIt

     printf("samplePriorEng: %E. priorEngOld: %E\n", samplePriorEng, priorEngOld);

     // convert log likelihood to energy.
     samplePriorEng = - samplePriorEng;

     // for data energy.
     std::vector< vnl_vector<double> >vmfLogConst(par.numSubs);
     double myD = par.tsLength;
     double const Pi = 4 * atan(1);
     unsigned subIdx = 0, clsIdx = 0;
     for (subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  vmfLogConst[subIdx].set_size(par.numClusters);
	  for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       if (par.vmm.sub[subIdx].comp[clsIdx].kappa > 1) {
		    // since the data energy is also divided by T --
		    // temperature, it's equlivalent to divide the kappa.
		    vmfLogConst[subIdx][clsIdx] = 
			 (myD/2 - 1) * log (par.vmm.sub[subIdx].comp[clsIdx].kappa / par.temperature)
			 - myD/2 * log(2*Pi) 
			 - logBesselI(myD/2 - 1, par.vmm.sub[subIdx].comp[clsIdx].kappa / par.temperature);
	       }
	       else {
		    vmfLogConst[subIdx][clsIdx] = 25.5; // the limit when kappa -> 0
	       }
	  }
     } // for subIdx.

     // cl is current label, s is subject id.
     unsigned short s = 0;
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  if (coordMap[nodeIt].subid < par.numSubs) {
	       s = coordMap[nodeIt].subid;
	       clsIdx = rSampleMap[nodeIt].find_first();
	       sampleDataEng += vmfLogConst[s][clsIdx] + par.vmm.sub[s].comp[clsIdx].kappa * inner_product(tsMap[nodeIt], par.vmm.sub[s].comp[clsIdx].mu);


	       // check if NaN.
	       if (sampleDataEng != sampleDataEng) {
		    printf("NAN happend at sub: %d, [%d %d %d]\n", s, (int)coordMap[nodeIt].idx[0], (int)coordMap[nodeIt].idx[1], (int)coordMap[nodeIt].idx[2]);
		    exit(1);
	       }

	  } // subid < numSubs
     } // nodeIt

     //The temperature parameter that change the definition of
     //the posterior distribution.
     sampleDataEng = sampleDataEng / par.temperature;

     // convert log likelihood to energy.
     sampleDataEng = - sampleDataEng;

     printf("samplePriorEng: %E, sampleDataEng: %E, single sample energy (- log P(X,Y) ): %E.\n", samplePriorEng, sampleDataEng, samplePriorEng + sampleDataEng);
     return 0;
}


int CompBetaDrv(lemon::SmartGraph & theGraph, 
		lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		ParStruct & par,
		double & drv1, double & drv2)
{
     double a = 0, b = 0, M0 = 0, M1 = 0, M2 = 0, acur = 0, bcur = 0;
     drv1 = 0, drv2 = 0;
     double priorEngy = 0;
     boost::dynamic_bitset<> curBit(par.numClusters), runningBit(par.numClusters);
     std::array<double, 3> distWeight = {{BETAWEIGHT0, BETAWEIGHT1, BETAWEIGHT2}};
     unsigned runningLabel = 0, weightIdx = 0;

     // save a*alpha+b*beta for all x = {1...L}
     vnl_vector<double> eta(par.numClusters, 0);
     double eta_offset = 0;

     unsigned pCount = 0;
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  curBit = rSampleMap[nodeIt];
	  eta = 0;
	  for (runningLabel = 0; runningLabel < par.numClusters; runningLabel ++) {
	       runningBit.reset(), runningBit[runningLabel] = 1;
	       // here we can not use edgeMap since the derivative with alpha
	       // and beta need to compute a and b separately.
	       a = 0, b = 0;
	       for (SmartGraph::IncEdgeIt edgeIt(theGraph, nodeIt); edgeIt != INVALID; ++ edgeIt) {
		    if ( (coordMap[theGraph.u(edgeIt)].subid == par.numSubs && coordMap[theGraph.v(edgeIt)].subid == par.numSubs)
			 || (coordMap[theGraph.u(edgeIt)].subid < par.numSubs && coordMap[theGraph.v(edgeIt)].subid < par.numSubs) ){
			 // within group or within subject
			 if (par.weightbetadist){
			      weightIdx = abs(coordMap[theGraph.u(edgeIt)].idx[0] - coordMap[theGraph.v(edgeIt)].idx[0]) 
				   + abs(coordMap[theGraph.u(edgeIt)].idx[1] - coordMap[theGraph.v(edgeIt)].idx[1]) 
				   + abs(coordMap[theGraph.u(edgeIt)].idx[2] - coordMap[theGraph.v(edgeIt)].idx[2])
				   -1;

			      b = b - distWeight[weightIdx] * double(runningBit != rSampleMap[theGraph.runningNode(edgeIt)]) ;
			 }
			 else {
			      b = b - double(runningBit != rSampleMap[theGraph.runningNode(edgeIt)]) ;
			 } // weightbetadist
		    } // within grp or sub.
		    else {
			 // btw group and subjects.
			 a = a - double(runningBit != rSampleMap[theGraph.runningNode(edgeIt)]) ;

		    }
	       } // incEdgeIt

	       a = a / par.temperature;
	       b = b / par.temperature;

	       // since the T -- temperature, 
	       eta[runningLabel] = a*par.alpha + b*par.beta;

	       if (runningBit == curBit) {
		    acur = a;
		    bcur = b;
	       }
	  } // runningLabel

	  // substract the max value such that the new max would be zero. This
	  // is to prevent underflow. Tehre may be still underflow for some
	  // smaller values, but we just ignore them.
	  eta_offset = eta.max_value();
	  M0 = 0, M1 = 0, M2 = 0;
	  for (runningLabel = 0; runningLabel < par.numClusters; runningLabel ++) {
	       M0 += exp(eta[runningLabel] - eta_offset);
	       M1 += b * exp(eta[runningLabel] - eta_offset);
	       M2 += b * b * exp(eta[runningLabel] - eta_offset);
	  } // runningLabel

	  drv1 += (bcur - M1/M0);
	  drv2 -= ( M2/M0 - pow(M1/M0, 2) );


	  // need compensate the min of exponential a*alpha+b*beta.
	  priorEngy += (acur * par.alpha + bcur * par.beta - eta_offset - log (M0) );

	  if (pCount < 30 && theGraph.id(nodeIt)%1 == 0) {
	       printf("node: %i, M0 = %E, a = %f, b = %f, eta_cur = %f, eta_offset = %f, priorEngy = %E\n",  theGraph.id(nodeIt), M0, a, b, acur * par.alpha + bcur * par.beta, eta_offset, priorEngy );
	       printVnlVector(eta, par.numClusters);
	       pCount ++;
	  }

     } // nodeIt
	  
     // convert log likelihood to energy.
     drv1 = - drv1, drv2 = - drv2;
     priorEngy = - priorEngy;

     return priorEngy;
}

int EstimateBeta(lemon::SmartGraph & theGraph, 
		lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		lemon::SmartGraph::EdgeMap<double> & edgeMap,
		ParStruct & par)
{
     double beta_o = 0, priorEngy = 0;;
     unsigned numIter = 0;
     double drv1 = 0, drv2 = 0;
     do {
	  numIter ++;
	  beta_o = par.beta;
	  priorEngy = CompBetaDrv(theGraph, coordMap, rSampleMap, par, drv1, drv2);
	  par.beta = par.beta - (drv1/drv2);
	  if (par.verbose >= 0) {
	       printf("drv1 = %E, drv2 = %E, beta_o = %1.3f, par.beta = %1.3f, priorEngy = %E\n", drv1, drv2, beta_o, par.beta, priorEngy);
	  }
     } while (fabs(beta_o - par.beta) > 1e-5 && numIter <= 0);
     return 0;
}


int PlotBeta(lemon::SmartGraph & theGraph, 
	     lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
	     lemon::SmartGraph::EdgeMap<double> & edgeMap,
	     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	     ParStruct & par)
{
 
     double beta_old = par.beta;
     for (par.beta = 0.01; par.beta < 1; par.beta = par.beta + 0.05) {
	  BuildEdgeMap(theGraph, edgeMap, coordMap, tsMap, par);
	  printf("beta = %2.2f, ", par.beta);
	  CompSampleEnergy(theGraph, coordMap, rSampleMap, edgeMap, tsMap, par); 
     }
     par.beta = beta_old;
     return 0;
}
