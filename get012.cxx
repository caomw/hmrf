#include <common.h>

int PrintMat(std::vector< std::vector< double > > mat);


int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     std::string s1file, s2file, s3file, maskfile, outfile;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", " Compute a 0-1-2 map. take value 1 if all 3 labels are equal, takes 1 if two are euqal, and take 2 if all diferent.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("s1", po::value<std::string>(&s1file)->default_value("s1file.nii.gz"), 
	   "segmentation label map for session 1.")
	  ("s2", po::value<std::string>(&s2file)->default_value("s2file.nii.gz"), 
	   "segmentation label map for session 2.")
	  ("s3", po::value<std::string>(&s3file)->default_value("s3file.nii.gz"), 
	   "segmentation label map for session 3.")
	  ("mask,m", po::value<std::string>(&maskfile)->default_value("mask.nii.gz"), 
	   "mask file with in-mask value > 0.")
	  ("out,o", po::value<std::string>(&outfile)->default_value("consistency.nii.gz"), 
	   "Output consistency map.");

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
     
     // output icc map file.
     ImageType3DFloat::IndexType conIndex;
     conIndex.Fill( 0 );
     ImageType3DFloat::RegionType conRegion;
     conRegion.SetSize(maskSize);
     conRegion.SetIndex(conIndex);
     ImageType3DFloat::Pointer conPtr = ImageType3DFloat::New();
     conPtr->SetRegions(conRegion);
     conPtr->Allocate();
     conPtr->FillBuffer( 0 );

     ReaderType3DFloat::Pointer s1Reader = ReaderType3DFloat::New();
     s1Reader->SetFileName( s1file);
     s1Reader->Update();
     ImageType3DFloat::Pointer s1Ptr = s1Reader->GetOutput();

     ReaderType3DFloat::Pointer s2Reader = ReaderType3DFloat::New();
     s2Reader->SetFileName( s2file);
     s2Reader->Update();
     ImageType3DFloat::Pointer s2Ptr = s2Reader->GetOutput();

     ReaderType3DFloat::Pointer s3Reader = ReaderType3DFloat::New();
     s3Reader->SetFileName( s3file);
     s3Reader->Update();
     ImageType3DFloat::Pointer s3Ptr = s3Reader->GetOutput();
     
     IteratorType3DFloat s1It(s1Ptr, s1Ptr->GetLargestPossibleRegion() );
     IteratorType3DFloat s2It(s2Ptr, s2Ptr->GetLargestPossibleRegion() );
     IteratorType3DFloat s3It(s3Ptr, s3Ptr->GetLargestPossibleRegion() );
     IteratorType3DFloat conIt(conPtr, conPtr->GetLargestPossibleRegion() );

     for (s1It.GoToBegin(), s2It.GoToBegin(), s3It.GoToBegin(), conIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ s1It, ++s2It, ++s3It, ++conIt, ++maskIt) {
	  if (maskIt.Get() > 0) {
	       if ( (s1It.Get() == s2It.Get()) && (s1It.Get() == s3It.Get()) ) {
		    conIt.Set(0);
	       }
	       else if ((s1It.Get()==s2It.Get()) || (s1It.Get() == s3It.Get()) || (s2It.Get() == s3It.Get()) ) {
		    conIt.Set(1);
	       }
	       else {
		    conIt.Set(2);
	       }
	  }
     }

     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
     writer->SetInput( conPtr);
     writer->SetFileName(outfile);

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

     printf("get012(): file %s saved.\n", outfile.c_str() );
     return 0;
}
