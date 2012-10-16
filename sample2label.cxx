#include <common.h>
#include <utility.h>
using namespace lemon;
int main(int argc, char* argv[])
{
     std::string inFile, outFile;
     
     
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "convert cumulative sample (4D) file to label file. at each voxel, label with largest posterior probability (i.e. most number of samples) are selected.")
	  ("in,i", po::value<std::string>(&inFile)->default_value("infile.nii"),
	   "input file.")
	  ("out,o", po::value<std::string>(&outFile)->default_value("outfile.nii"),
	   "output file.");

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

     ReaderType4DS::Pointer reader = ReaderType4DS::New();
     reader->SetFileName(inFile);
     reader->Update();
     ImageType4DS::Pointer inPtr = reader->GetOutput();
     ImageType4DS::SizeType inSize = inPtr->GetLargestPossibleRegion().GetSize();


     ImageType3DShort::IndexType outImageIdx;
     outImageIdx.Fill(0);
     ImageType3DShort::SizeType outImageSize;
     outImageSize[0] = inSize[0];
     outImageSize[1] = inSize[1];
     outImageSize[2] = inSize[2];

     ImageType3DShort::RegionType outImageRegion;
     outImageRegion.SetSize(outImageSize);
     outImageRegion.SetIndex(outImageIdx);
     ImageType3DShort::Pointer outImagePtr = ImageType3DShort::New();
     outImagePtr->SetRegions(outImageRegion);
     outImagePtr->Allocate();
     outImagePtr->FillBuffer(0);


     IteratorType3DShort outIt(outImagePtr, outImagePtr->GetLargestPossibleRegion() ) ;
     ImageType3DShort::IndexType outIdx;
     ImageType4DS::IndexType inIdx;
     unsigned bestCls, mostSamples;
     for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++ outIt) {
	  outIdx = outIt.GetIndex();
	  inIdx[0] = outIdx[0];
	  inIdx[1] = outIdx[1];
	  inIdx[2] = outIdx[2];
	  mostSamples = 0;
	  for (inIdx[3] = 0; inIdx[3] < inSize[3]; inIdx[3] ++) {
	       if (inPtr->GetPixel(inIdx) > mostSamples) {
		    mostSamples = inPtr->GetPixel(inIdx);
		    bestCls = inIdx[3];
	       }
	  }
	  if (mostSamples > 0) {
	       outIt.Set(bestCls + 1);
	  }
	  else {
	       // must be outside of mask.
	       outIt.Set(0);
	  }
     }

     WriterType3DShort::Pointer writer = WriterType3DShort::New();
	  
     writer->SetInput(outImagePtr);
     writer->SetFileName(outFile);
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

     std::cout << "sample2label(): File " << outFile << " saved.\n";

     return 0;
}
