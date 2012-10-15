#include <common.h>
#include <utility.h>
using namespace lemon;
int main(int argc, char* argv[])
{
     std::string inFile, outFile;
     
     
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")
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

     ReaderType3DChar::Pointer reader = ReaderType3DChar::New();
     reader->SetFileName(inFile);
     reader->Update();
     ImageType3DChar::Pointer inPtr = reader->GetOutput();
     ImageType3DChar::SizeType inSize = inPtr->GetLargestPossibleRegion().GetSize();



     // Compute the maximum intensity value of input image.
     typedef itk::MinimumMaximumImageCalculator <ImageType3DChar>  ImageCalculatorFilterType;
     ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
     imageCalculatorFilter->SetImage(inPtr);
     imageCalculatorFilter->Compute();
     unsigned numClusters = imageCalculatorFilter->GetMaximum();
     printf("num of clusters: %d\n", numClusters);

     ImageType4DS::IndexType outImageIdx;
     outImageIdx.Fill(0);
     ImageType4DS::SizeType outImageSize;
     outImageSize[0] = inSize[0];
     outImageSize[1] = inSize[1];
     outImageSize[2] = inSize[2];
     outImageSize[3] = numClusters;

     ImageType4DS::RegionType outImageRegion;
     outImageRegion.SetSize(outImageSize);
     outImageRegion.SetIndex(outImageIdx);
     ImageType4DS::Pointer outImagePtr = ImageType4DS::New();
     outImagePtr->SetRegions(outImageRegion);
     outImagePtr->Allocate();
     outImagePtr->FillBuffer(0);


     IteratorType3DChar inIt(inPtr, inPtr->GetLargestPossibleRegion() ) ;
     ImageType3DChar::IndexType inIdx = inIt.GetIndex();
     for (inIt.GoToBegin(); !inIt.IsAtEnd(); ++ inIt) {
	  if (inIt.Get() > 0) {
	       inIdx = inIt.GetIndex();
	       outImageIdx[0] = inIdx[0];
	       outImageIdx[1] = inIdx[1];
	       outImageIdx[2] = inIdx[2];
	       outImageIdx[3] = inIt.Get() - 1;
	       outImagePtr->SetPixel(outImageIdx, 1);
	  } // inIt > 0
     }

     WriterType4DS::Pointer writer = WriterType4DS::New();
	  
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
