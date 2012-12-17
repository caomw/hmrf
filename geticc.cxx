#include <common.h>

struct ICC{
     double rho_u;
     double rho_c;
};

int PrintMat(std::vector< std::vector< double > > mat);

int PrintPathMat(PathMat pathMat);

int LoadData(std::string datapath,
	     ImageTypeMat::Pointer dataPtr,
	     ImageType3DChar::Pointer maskPtr);

ICC ComputeICC(const std::vector< std::vector< double > > & mat);


int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     std::string inpath, iccufile, icccfile, maskfile;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Compute intraclass correlation coefficient by using Zuo's test-retest method.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("inpath,i", po::value<std::string>(&inpath)->default_value("."), 
	   "Input data path. This path contains multiple session folder, which contains subject files.")
	  ("mask,m", po::value<std::string>(&maskfile)->default_value("mask.nii.gz"), 
	   "mask file with in-mask value > 0.")
	  ("iccu,u", po::value<std::string>(&iccufile)->default_value("iccumap.nii.gz"), 
	   "Output icc map that does penalize system error")
	  ("iccc,c", po::value<std::string>(&icccfile)->default_value("icccmap.nii.gz"), 
	   "Output icc map that does NOT penalize system error.");

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
     ImageType3DFloat::Pointer icccPtr = ImageType3DFloat::New();
     ImageType3DFloat::Pointer iccuPtr = ImageType3DFloat::New();
     icccPtr->SetRegions(iccRegion), iccuPtr->SetRegions(iccRegion);
     icccPtr->Allocate(), iccuPtr->Allocate();
     icccPtr->FillBuffer( 0 ), iccuPtr->FillBuffer( 0 );
     
     IteratorType3DFloat icccIt(icccPtr, icccPtr->GetLargestPossibleRegion() );
     IteratorType3DFloat iccuIt(iccuPtr, iccuPtr->GetLargestPossibleRegion() );

     // compute ICC_u and ICC_c
     ImageType3DChar::IndexType maskIdx;
     maskIdx[0] = 44;
     maskIdx[1] = 33;
     maskIdx[2] = 30;
     ICC icc;
     for (icccIt.GoToBegin(), iccuIt.GoToBegin(), maskIt.GoToBegin(), dataIt.GoToBegin(); !maskIt.IsAtEnd(); ++ icccIt, ++ iccuIt, ++ maskIt, ++ dataIt) {
	  if(maskIt.Get() > 0) {
	       icc = ComputeICC( dataIt.Get() );
	       icccIt.Set( icc.rho_c );
	       iccuIt.Set( icc.rho_u );

	       if (maskIt.GetIndex() == maskIdx) {
		    PrintMat( dataIt.Get() );
	       }
	  }
     }

     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
     writer->SetInput( icccPtr);
     writer->SetFileName(icccfile);

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

     writer->SetInput( iccuPtr);
     writer->SetFileName(iccufile);

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
	       std::cout << "Reading file " << sessionSubPath[sessionIdx][subIdx] << '\n';
	       indataReader->SetFileName( sessionSubPath[sessionIdx][subIdx].string() );
	       // dataReader->Update();
	       indataPtr->Update();
	       IteratorType3DFloat indataIt(indataPtr, indataPtr->GetLargestPossibleRegion() );
	       for (indataIt.GoToBegin(), outdataIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ indataIt, ++ outdataIt, ++ maskIt) {
		    // retrieve the matrix, update one value and save it back.
		    if (maskIt.Get() > 0) {
			 thisMat = outdataIt.Get();
			 thisMat[subIdx][sessionIdx] = indataIt.Get();
			 outdataIt.Set( thisMat );
			 // PrintMat(outdataIt.Get());
			 // std::cout << '\n';
		    } // maskIt
	       } // indataIt
	  } // subIdx
	  std::cout << "\n";
     } // sessionIdx

     return 0;
}

ICC ComputeICC(const std::vector< std::vector< double > > & mat)
{
     ICC icc;
     unsigned numSubs = mat.size();
     unsigned numSessions = mat[0].size();
     std::vector<double> subMean(numSubs, 0), sessionMean(numSessions, 0);
     double groundMean = 0;
     double ssp = 0, sst = 0, sse = 0;

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
     
     // compute SS_t
     for (unsigned sessionIdx = 0; sessionIdx < numSessions; sessionIdx ++) {
	  sst += numSubs * pow((sessionMean[sessionIdx] - groundMean), 2);
     }

     // compute SS_e
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  for (unsigned sessionIdx = 0; sessionIdx < numSessions; sessionIdx ++) {
	       sse += pow((mat[subIdx][sessionIdx] - subMean[subIdx] - sessionMean[sessionIdx] + groundMean), 2);
	  }
     }

     // compute sigma_p estimates.
     double MS_p = ssp / (numSubs - 1);
     double MS_t = sst / (numSessions - 1);
     double MS_e = sse / ( (numSubs-1)*(numSessions-1) );

     icc.rho_u = (MS_p - MS_e) / ( MS_p + (numSessions-1)*MS_e + ((double)(numSessions) / (double)(numSubs))*(MS_t - MS_e) );
     icc.rho_c = (MS_p - MS_e) / ( MS_p + (numSessions-1)*MS_e );
     
     printf("Y_dd=%.1f, SS=%.1f, SS_p=%.1f, SS_w=%.1f, SS_t=%.1f, SS_e=%.1f, MS_p=%.1f, MS_t=%.1f, MS_e=%.1f, ICC_u=%.1f, ICC_c=%.1f\n", 
     	    groundMean, ss, ssp, ssw, sst, sse, MS_p, MS_t, MS_e, icc.rho_u, icc.rho_c);
     
     return icc;
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
	     
