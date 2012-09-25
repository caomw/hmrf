double logBesselI(float nu, double x);

int SaveGrpCumuSamples(lemon::SmartGraph & theGraph, 
		       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		       lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
		       ImageType3DChar::Pointer maskPtr,
		       ParStruct & par,
		       std::string outFile);

int SaveSubCumuSamples(lemon::SmartGraph & theGraph, 
		       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		       lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
		       ImageType3DChar::Pointer maskPtr,
		       ParStruct & par,
		       std::string basename);

int PrintPar(unsigned prlevel, ParStruct & par);

int save3dcharInc(ImageType3DChar::Pointer ip, std::string fname);

int printVnlVector(vnl_vector<float> vec, unsigned numElements);
