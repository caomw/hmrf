#include <common.h>
#include <utility.h>

int main(int argc, char* argv[])
{
     std::string infile, maskfile, outfile;
     unsigned verbose = 0;
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Apply mask to image.")
	   ("inputimage,i", po::value<std::string>(&infile)->default_value("input.nii.gz"), 
	    "Input image to be masked.")
	   ("maskimae,m", po::value<std::string>(&maskfile)->default_value("mask.nii.gz"), 
	    "mask image file.")
	   ("outputimage,o", po::value<std::string>(&outfile)->default_value("outputimage.nii.gz"), 
	    "output image file.")
	  ("verbose,v", po::value<unsigned>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: applymask [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    


     ReaderType3DShort::Pointer inReader = ReaderType3DShort::New();
     inReader->SetFileName(infile);
     inReader->Update();
     ImageType3DShort::Pointer inPtr = inReader->GetOutput();

     ReaderType3DShort::Pointer maskReader = ReaderType3DShort::New();
     maskReader->SetFileName(maskfile);
     maskReader->Update();
     ImageType3DShort::Pointer maskPtr = maskReader->GetOutput();
     
     IteratorType3DShort inIt(inPtr, inPtr->GetLargestPossibleRegion());
     ImageType3DShort::PointType pixelPoint;
     ImageType3DShort::IndexType inIdx, maskIdx;
     for (inIt.GoToBegin(); !inIt.IsAtEnd(); ++ inIt) {
	  if (inIt.Get() > 0) {
	       inIdx = inIt.GetIndex();
	       inPtr->TransformIndexToPhysicalPoint(inIdx, pixelPoint);
	       maskPtr->TransformPhysicalPointToIndex(pixelPoint, maskIdx);
	       
	       if (maskPtr->GetPixel(maskIdx) <= 0) {
		    inIt.Set( 0 );
	       }
	  }    

     }

     WriterType3DShort::Pointer writer = WriterType3DShort::New();
     writer->SetInput( inPtr );
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
     std::cout << "applymask(): File  " << outfile << " saved.\n";

     return 0;
}
