#include <common.h>

int PrintMat(std::vector< std::vector< double > > mat);

int PrintPathMat(PathMat pathMat);

int LoadData(std::string datapath,
	     ImageTypeMat::Pointer dataPtr,
	     ImageType3DChar::Pointer maskPtr);

double ComputeICC(const std::vector< std::vector< double > > & mat);

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     std::string inpath, iccufile, iccfile, maskfile;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Compute intraclass correlation coefficient by assuming a one-way ANOVA model.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("inpath,i", po::value<std::string>(&inpath)->default_value("."), 
	   "Input data path. This path contains multiple session folder, which contains subject files.")
	  ("mask,m", po::value<std::string>(&maskfile)->default_value("mask.nii.gz"), 
	   "mask file with in-mask value > 0.")
	  ("icc,c", po::value<std::string>(&iccfile)->default_value("iccmap.nii.gz"), 
	   "Output icc map.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: compareclusterings [options]\n";
	       std::cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // std::vector< ImageType3DFloat::Pointer > D;

     ImageTypeMat::Pointer dataPtr = ImageTypeMat::New();

     // mask file
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskfile);
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();
     maskPtr->Update();
     ImageType3DChar::SizeType maskSize =  maskPtr->GetLargestPossibleRegion().GetSize();     
     IteratorType3DCharIdx maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // load data.
     LoadData(inpath, dataPtr, maskPtr);
     IteratorTypeMat dataIt(dataPtr, dataPtr->GetLargestPossibleRegion() );

     
     // output icc map file.
     ImageType3DFloat::IndexType iccIndex;
     iccIndex.Fill( 0 );
     ImageType3DFloat::RegionType iccRegion;
     iccRegion.SetSize(maskSize);
     iccRegion.SetIndex(iccIndex);
     ImageType3DFloat::Pointer iccPtr = ImageType3DFloat::New();
     iccPtr->SetRegions(iccRegion);
     iccPtr->Allocate();
     iccPtr->FillBuffer( 0 );
     
     IteratorType3DFloat iccIt(iccPtr, iccPtr->GetLargestPossibleRegion() );

     // compute ICC_u and ICC_c
     double icc;
     for (iccIt.GoToBegin(), maskIt.GoToBegin(), dataIt.GoToBegin(); !maskIt.IsAtEnd(); ++ iccIt, ++ maskIt, ++ dataIt) {
	  if(maskIt.Get() > 0) {
	       icc = ComputeICC( dataIt.Get() );
	       iccIt.Set( icc );
	  }
     }

     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
     writer->SetInput( iccPtr);
     writer->SetFileName(iccfile);

     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     }      
     return 0;
}

int LoadData(std::string datapath,
	     ImageTypeMat::Pointer outdataPtr,
	     ImageType3DChar::Pointer maskPtr)
{
     // create itk image for the data.
     ImageType3DChar::SizeType maskSize =  maskPtr->GetLargestPossibleRegion().GetSize();
     ImageTypeMat::IndexType dataIdx;
     dataIdx.Fill( 0 );
     ImageTypeMat::RegionType dataRegion;
     dataRegion.SetSize(maskSize);
     dataRegion.SetIndex(dataIdx);
     outdataPtr->SetRegions( dataRegion );
     outdataPtr->Allocate();
     IteratorTypeMat outdataIt(outdataPtr, outdataPtr->GetLargestPossibleRegion() );
     IteratorType3DCharIdx maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // prepare for sorting directory entries.
     PathVec  sessionEntries;
     copy(boost::filesystem::directory_iterator(datapath), boost::filesystem::directory_iterator(), std::back_inserter(sessionEntries) );
     sort(sessionEntries.begin(), sessionEntries.end() );

     unsigned numSessions = sessionEntries.size();
     unsigned numSubs = 0;

     PathMat sessionSubPath;
     PathVec subEntries;
     for (PathVec::const_iterator sessionDirIt(sessionEntries.begin()), sessionDirIt_end(sessionEntries.end()); sessionDirIt != sessionDirIt_end; ++ sessionDirIt)
     {
	  // obtain sub path for current session path.
	  subEntries.clear();
	  copy(boost::filesystem::directory_iterator(*sessionDirIt), boost::filesystem::directory_iterator(), std::back_inserter(subEntries) );
	  sort(subEntries.begin(), subEntries.end());
	  sessionSubPath.push_back(subEntries);
     }

     // PrintPathMat(sessionSubPath);
     numSubs = sessionSubPath[0].size();

     // allocate memeory for thisMat.
     std::vector< std::vector< double > >  thisMat;
     thisMat.resize(numSubs);
     for (std::vector< std::vector< double> >::iterator it = thisMat.begin(); it != thisMat.end(); ++ it) {
	  (*it).resize(numSessions, 0);
     }

     // init outdata to all zero at each voxel.
     for (outdataIt.GoToBegin(); !outdataIt.IsAtEnd(); ++ outdataIt) {
	  outdataIt.Set( thisMat );
     } // dataIt

     ReaderType3DFloat::Pointer indataReader = ReaderType3DFloat::New();
     ImageType3DFloat::Pointer indataPtr = indataReader->GetOutput();

     
     // read data from file in each session, each subject.
     for (unsigned sessionIdx = 0; sessionIdx < numSessions; sessionIdx ++) {
	  for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	       std::cout << "Reading " << sessionSubPath[sessionIdx][subIdx] << '\n';
	       indataReader->SetFileName( sessionSubPath[sessionIdx][subIdx].string() );
	       indataPtr->Update();
	       IteratorType3DFloat indataIt(indataPtr, indataPtr->GetLargestPossibleRegion() );
	       for (indataIt.GoToBegin(), outdataIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ indataIt, ++ outdataIt, ++ maskIt) {
		    // retrieve the matrix, update one value and save it back.
		    if (maskIt.Get() > 0) {
			 thisMat = outdataIt.Get();
			 thisMat[subIdx][sessionIdx] = indataIt.Get();
			 outdataIt.Set( thisMat );
			 // std::cout << '\n';
		    } // maskIt
	       } // indataIt
	  } // subIdx
	  std::cout << "\n";
     } // sessionIdx

     return 0;
}

double ComputeICC(const std::vector< std::vector< double > > & mat)
{
     double icc;
     unsigned numSubs = mat.size();
     unsigned numSessions = mat[0].size();
     std::vector<double> subMean(numSubs, 0), sessionMean(numSessions, 0);
     double groundMean = 0;
     double ssp = 0;
     double colVarHat = 0, errorVarHat = 0;
     
     // compute subMean (Y_idot)
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  for (unsigned sessionIdx = 0; sessionIdx < numSessions; sessionIdx ++) {
	       subMean[subIdx] += mat[subIdx][sessionIdx];
	  }
	  subMean[subIdx] = subMean[subIdx] / numSessions;
     }

     // compute sessionMean (Y_dotj)
     for (unsigned sessionIdx = 0; sessionIdx < numSessions; sessionIdx ++) {
	  for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	       sessionMean[sessionIdx] += mat[subIdx][sessionIdx];
	  }
	  sessionMean[sessionIdx] = sessionMean[sessionIdx] / numSubs;
     }
     
     // compute groundMean (Y_dotdot)
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  for (unsigned sessionIdx = 0; sessionIdx < numSessions; sessionIdx ++) {
	       groundMean += mat[subIdx][sessionIdx];
	  }
     }
     groundMean = groundMean / (numSubs * numSessions);

     // SS. Only used for debugging.
     double ss = 0;
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  for (unsigned sessionIdx = 0; sessionIdx < numSessions; sessionIdx ++) {
	       ss += pow((mat[subIdx][sessionIdx] - groundMean), 2);
	  }
     }

     // compute SS_p
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  ssp += numSessions * pow((subMean[subIdx] - groundMean), 2);
     }

     // SS_w. Only used for debugging.
     double ssw = 0;
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  for (unsigned sessionIdx = 0; sessionIdx < numSessions; sessionIdx ++) {
	       ssw += pow((mat[subIdx][sessionIdx] - subMean[subIdx]), 2);
	  }
     }

     // compute sigma_p estimates.
     double MS_p = ssp / (numSubs - 1);
     double MS_w = ssw / ( numSubs * (numSessions-1) );

     // variance estimation.
     colVarHat = (MS_p - MS_w) / numSessions;
     errorVarHat = MS_w;

     if (MS_p == 0 && MS_w == 0) {
	  icc = 0;
     }
     else {
	  icc = (MS_p - MS_w) / ( MS_p + (numSessions-1) * MS_w );
     }

     printf("Y_dd=%.1f, SS=%.1f, SS_p=%.1f, SS_w=%.1f, MS_p=%.1f, MS_w=%.1f, sigma2_p=%.1f, sigma2_w=%.1f, ICC=%1.2f\n", 
     	    groundMean, ss, ssp, ssw, MS_p, MS_w, colVarHat, errorVarHat, icc);
     return MS_w;
}

int PrintMat(std::vector< std::vector< double > >  mat)
{
     unsigned numSubs = mat.size();
     unsigned numSessions = mat[0].size();

     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  for (unsigned sessionIdx = 0; sessionIdx < numSessions; sessionIdx ++) {
	       printf("%1.1f ", mat[subIdx][sessionIdx]);
	  }
	  printf("\n");
     }
     
     return 0;
}

int PrintPathMat(PathMat pathMat)
{
     for (PathMat::iterator it = pathMat.begin(); it != pathMat.end(); ++it) {

	  std::cout << "this session has sub: " << (*it).size() << '\n';
	  for (PathVec::iterator subIt = (*it).begin(); subIt != (*it).end(); ++ subIt) {
	       std::cout << *subIt << '\n';
	  }
	  std::cout << '\n';
     }
     
     return 0;
}
	     
