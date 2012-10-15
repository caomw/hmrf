#include <common.h>
#include <utility.h>
using namespace lemon;
int main(int argc, char* argv[])
{
     std::string inFile, outFile;
     
     
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "convert indicator image file to label file.")
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


     ImageType3DChar::IndexType outImageIdx;
     outImageIdx.Fill(0);
     ImageType3DChar::SizeType outImageSize;
     outImageSize[0] = inSize[0];
     outImageSize[1] = inSize[1];
     outImageSize[2] = inSize[2];

     ImageType3DChar::RegionType outImageRegion;
     outImageRegion.SetSize(outImageSize);
     outImageRegion.SetIndex(outImageIdx);
     ImageType3DChar::Pointer outImagePtr = ImageType3DChar::New();
     outImagePtr->SetRegions(outImageRegion);
     outImagePtr->Allocate();
     outImagePtr->FillBuffer(0);


     IteratorType4DS inIt(inPtr, inPtr->GetLargestPossibleRegion() ) ;
     ImageType4DS::IndexType inIdx = inIt.GetIndex();
     for (inIt.GoToBegin(); !inIt.IsAtEnd(); ++ inIt) {
	  if (inIt.Get()  >= 1) {
	       inIdx = inIt.GetIndex();
	       outImageIdx[0] = inIdx[0];
	       outImageIdx[1] = inIdx[1];
	       outImageIdx[2] = inIdx[2];
	       outImagePtr->SetPixel(outImageIdx, inIdx[3] + 1);
	  } // inIt > 0
     }

     WriterType3DChar::Pointer writer = WriterType3DChar::New();
	  
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

     std::cout << "label2ind(): File " << outFile << " saved.\n";

     return 0;
}
