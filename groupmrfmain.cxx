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

int InitSubSamples(std::string srcpath,
		   lemon::SmartGraph & theGraph, 
		   lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		   lemon::SmartGraph::NodeMap<unsigned short> & initSubjectMap);

int EstimateMu(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	       lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	       lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	       ParStruct & par);

int EstiamteKappa(lemon::SmartGraph & theGraph, 
		  ParStruct & par);

double EstimatePriorPar(lemon::SmartGraph & theGraph, 
			lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
			lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
			ParStruct & par);

int Sampling(lemon::SmartGraph & theGraph, 
	     lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	     lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	     ParStruct & par);

void *SamplingThreads(void * threadArgs);
void *PrintHello(void *threadarg);

using namespace lemon;
int main(int argc, char* argv[])
{
     ParStruct par;
     unsigned seed = 0;
     unsigned int emIter = 1;

     std::string dataFile, iLabelFile, oLabelFile, sampleFile;
     std::string fmriPath, initGrplabel, outGrpLabel, outGrpProb, outSubBase, 
	  grpSampleName, samplePrefix, initsubpath,  groupprob, subbaseprob;

     bool estprior = false;
     bool initsame = false;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")
	  ("burnin,b", po::value<unsigned>(&par.burnin)->default_value(20),
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
	   "whether to estimate alpha and beta, the prior parameter.")
	  ("initsame", po::bool_switch(&initsame), 
	   "whether to init subject same to group label. choose yes to init same, otherwise init with subject label map.")

	   ("initgrouplabel,i", po::value<std::string>(&initGrplabel)->default_value("grouplabel.nii"), 
	    "Initial group level label map (Intensity value 1-K. Also used as mask file.")
	   ("initsubpath,t", po::value<std::string>(&initsubpath)->default_value("subpath"), 
	    "Initial subject label maps path")
	  ("fmripath,f", po::value<std::string>(&fmriPath)->default_value("."), 
	   "noised image file")

	  ("sampleprefix", po::value<std::string>(&samplePrefix)->default_value("subject"), 
	   "Monte Carlo samples file name prefix")
	   ("subbasename", po::value<std::string>(&outSubBase)->default_value("subject"), 
	    "Output individual subject's base name (with path).")

	  ("grouplabel,g", po::value<std::string>(&outGrpLabel)->default_value("outgrplabel.nii"), 
	   "output group label file. nii format.")
	  ("grpsamplename", po::value<std::string>(&grpSampleName)->default_value("grpSample.nii.gz"), 
	   "Monte Carlo group label samples file name.")

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

     // get number of subjects.
     par.numSubs = ComputeNumSubs(fmriPath);
     printf("number of subjects: %i.\n", par.numSubs);

     // Construct graph.
     lemon::SmartGraph theGraph;
     lemon::SmartGraph::NodeMap<SuperCoordType> coordMap(theGraph); // node --> coordinates.

     BuildGraph(theGraph, coordMap, par.numSubs, maskPtr);

     // Build Data map.
     lemon::SmartGraph::NodeMap<vnl_vector<float>> tsMap(theGraph); // node --> time series.
     BuildDataMap(theGraph, coordMap, tsMap, fmriPath, par);

     // define a subject initial label map used for initialization of subject
     // label map.
     lemon::SmartGraph::NodeMap<unsigned short> initSubjectMap(theGraph);     
     InitSubSamples(initsubpath, theGraph, coordMap, initSubjectMap);

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
	       cumuSampleMap[nodeIt].at(labelPtr->GetPixel(coordMap[nodeIt].idx)) = par.numSamples;
	  }
	  else {
	       // only when not initsame for subjects nodes, we init it to
	       // subjects initial label map.
	       cumuSampleMap[nodeIt].at(initSubjectMap[nodeIt]) = par.numSamples;
	  }
     }
     
     // define a running sample map and init it with the input group
     // labels. This is used for sampling only. And when the whoel scan is done,
     // this map is saved into the cumuSampleMap.
     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > rSampleMap(theGraph);
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  rSampleMap[nodeIt].resize(par.numClusters, 0);
	  if (coordMap[nodeIt].subid < par.numSubs && (!initsame) ) {
	       rSampleMap[nodeIt].set(labelPtr->GetPixel(coordMap[nodeIt].idx), true);
	  }
	  else {
	       // only when not initsame for subjects nodes, we init it to
	       // subjects initial label map.
	       rSampleMap[nodeIt].set(initSubjectMap[nodeIt], true);
	  }
     }

     // init parameters.
     par.vmm.sub.resize(par.numSubs);
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  par.vmm.sub[subIdx].comp.resize(par.numClusters);
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       par.vmm.sub[subIdx].comp[clsIdx].mu.set_size(par.tsLength);
	       par.vmm.sub[subIdx].comp[clsIdx].meanNorm = 0;
	       par.vmm.sub[subIdx].comp[clsIdx].kappa = 0;
	       par.vmm.sub[subIdx].comp[clsIdx].numPts = 0;
	       par.vmm.sub[subIdx].comp[clsIdx].prop = 0;
	  }
     }
     
     Sampling(theGraph, coordMap, cumuSampleMap, tsMap, par);

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
     std::array<unsigned int, 6 > neiIdxSet = {14, 12, 16, 10, 22, 4}; 
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
			      nodeMapIdx[3] = numSubs;
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
		   lemon::SmartGraph::NodeMap<unsigned short> & initSubjectMap)
{
     // some of the code below refers to boost filesystem tutorial tut04.cpp.
     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(srcpath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );
     unsigned numSubs  = sortedEntries.size();

     printf("InitSubSamples(): numSubs = %i\n", numSubs);

     // init file reader pointers.
     std::vector<ReaderType3DChar::Pointer> readerVec(numSubs);
     std::vector<ImageType3DChar::Pointer> srcPtrVec(numSubs);
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
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
	  if (superCoord.subid < numSubs) {
	       initSubjectMap[nodeIt] = srcPtrVec[superCoord.subid]->GetPixel(superCoord.idx) - 1;
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

     // reset all mu to zero.
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       par.vmm.sub[subIdx].comp[clsIdx].mu = 0;
	  }
     }

     SuperCoordType superCoord;
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  superCoord = coordMap[nodeIt];
	  for (unsigned clsIdx = 0; clsIdx < par.numSubs; clsIdx ++) {
	       par.vmm.sub[superCoord.subid].comp[clsIdx].mu += tsMap[nodeIt] * cumuSampleMap[nodeIt][clsIdx];
	       par.vmm.sub[superCoord.subid].comp[clsIdx].numPts += cumuSampleMap[nodeIt][clsIdx];
	  }
     }

     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {     
	       if(par.vmm.sub[superCoord.subid].comp[clsIdx].numPts > 0) {
		    par.vmm.sub[superCoord.subid].comp[clsIdx].meanNorm = par.vmm.sub[superCoord.subid].comp[clsIdx].mu.two_norm() / par.vmm.sub[superCoord.subid].comp[clsIdx].numPts;
		    par.vmm.sub[superCoord.subid].comp[clsIdx].mu.normalize();
	       }
	       else {
		    par.vmm.sub[superCoord.subid].comp[clsIdx].meanNorm = 0;
	       }
	  }
     }
	       
     return 0;
}

int EstiamteKappa(lemon::SmartGraph & theGraph, 
		  ParStruct & par)
{
     float kappa = 0, kappa_new = 0;
     double Ap = 0;
     float RBar = 0;
     float Dim = par.tsLength;

     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {

	       RBar = par.vmm.sub[subIdx].comp[clsIdx].meanNorm;
	       kappa_new = RBar * (Dim - RBar * RBar) / (1 - RBar * RBar);
	       if (kappa_new  == 0 ) {
		    par.vmm.sub[subIdx].comp[clsIdx].kappa = kappa_new;
		    continue;
	       }
	       
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
	       par.vmm.sub[subIdx].comp[clsIdx].kappa = kappa_new;
	  } // clsIdx
     } // subIdx

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
	     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	     ParStruct & par)
{
     unsigned numThreads = 1;
     unsigned taskid = 0;
     pthread_t thread_ids[numThreads];
     ThreadArgs threadArgs[numThreads];

     for (unsigned taskid = 0; taskid < numThreads; taskid ++) {
	  threadArgs[taskid].taskid = taskid;
	  threadArgs[taskid].theGraph = & theGraph;
	  threadArgs[taskid].coordMap = & coordMap;
	  threadArgs[taskid].cumuSampleMap = & cumuSampleMap;
	  threadArgs[taskid].tsMap = & tsMap;
	  threadArgs[taskid].par = & par;
	  threadArgs[taskid].numThreads = numThreads;
     }

     // for (taskid = 0; taskid < numThreads; taskid ++) {
     // 	  printf("threadArgs[%i].taskid = %i\n", taskid, threadArgs[taskid].taskid);
     // 	  printf("threadArgs[%i].numThreads = %i\n", taskid, threadArgs[taskid].numThreads);
     // }

     // for (extaskid = 0; taskid < numThreads; taskid ++) {
     // 	  pthread_create(&thread_ids[taskid], NULL, SamplingThreads, (void*) &threadArgs[taskid]);
     // }


// ------------ test ----------------------
     unsigned NUM_THREADS = 10;
     struct thread_data thread_data_array[NUM_THREADS];
     pthread_t threads[NUM_THREADS];
     int sum = 50;
     
     for (taskid = 0; taskid < NUM_THREADS; taskid ++) {
	  thread_data_array[taskid].thread_id = taskid;
	  thread_data_array[taskid].sum = sum;
	  thread_data_array[taskid].message = "this is hello world";

	  printf("main, taskid = %i, sum = %i\n", taskid, sum);

	  pthread_create(&threads[taskid], NULL, PrintHello, 
			      (void *) &thread_data_array[taskid]);
     }

     for (taskid = 0; taskid < NUM_THREADS; taskid ++) {
	  pthread_join(threads[taskid], NULL);
     }
     return 0;
}

void *PrintHello(void *threadarg)
{
   struct thread_data *my_data;
   my_data = (struct thread_data *) threadarg;
   int taskid = my_data->thread_id;
   int sum = my_data->sum;
   char *hello_msg = my_data->message;
   
   printf("Thread, taskid = %i, sum = %i\n", taskid, sum);
}

void *SamplingThreads(void * threadArgs)
{
     ThreadArgs * args = (ThreadArgs *) threadArgs;
     unsigned taskid = args->taskid;
     lemon::SmartGraph * theGraph = args->theGraph;
     lemon::SmartGraph::NodeMap<SuperCoordType> * coordMap = args->coordMap;
     lemon::SmartGraph::NodeMap<std::vector<unsigned short> > * cumuSampleMap = args->cumuSampleMap;
     lemon::SmartGraph::NodeMap<vnl_vector<float>> * tsMap = args->tsMap;
     ParStruct * par = args->par;
     unsigned numThreads = args->numThreads;

     printf("Thread id %ld, task id %i. NumThreads = %i\n", pthread_self(), taskid, numThreads);

     printf("Thread id %ld, task id %i. NumThreads = %i\n", pthread_self(), args->taskid, args->numThreads);
}





