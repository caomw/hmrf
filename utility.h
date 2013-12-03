double logBesselI(float nu, double x);

int SaveCumuSamples(lemon::SmartGraph & theGraph, 
		       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		       lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
		    ParStruct & par);


int SaveRunningSamples(lemon::SmartGraph & theGraph, 
		       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		       lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		       ParStruct & par);

int PrintPar(unsigned prlevel, ParStruct & par);

int save3dcharInc(ImageType3DChar::Pointer ip, std::string fname);
int save3dchar(ImageType3DChar::Pointer ip, std::string fname);
int printVnlVector(vnl_vector<float> vec, unsigned numElements);
int printVnlVector(vnl_vector<double> vec, unsigned numElements);

