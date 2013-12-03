#include <common.h>
#include <utility.h>

twister_base_gen_type mygenerator(42u);

using namespace lemon;
int main(int argc, char* argv[])
{
     unsigned seed = 0;
     unsigned tslength = 0;
     unsigned verbose = 0;
     std::string infile, outfile;

     // program options.
     po::options_description mydesc("Options can be used at commandline");
     mydesc.add_options()
	  ("help,h", "resample time points of each time course with replacement.")
	   ("input,i", po::value<std::string>(&infile)->default_value("fmri.nii.gz"), 
	    "Input fmri file for resampling.")
	   ("output,o", po::value<std::string>(&outfile)->default_value("outsample.nii.gz"), 
	    "output of resampled fmri file.")
	  ("tslength,k", po::value<unsigned>(&tslength)->default_value(300),
	   "Length of output time series. Should be same with length of input time series.")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("verbose,v", po::value<unsigned>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: fmriresample [options]\n";
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


     // read input fmri file.
     
     ReaderType4DFloat::Pointer inReader = ReaderType4DFloat::New();
     inReader->SetFileName(infile);
     inReader->Update();
     ImageType4DF::Pointer inPtr = inReader->GetOutput();
     ImageType4DF::SizeType inSize = inPtr->GetLargestPossibleRegion().GetSize();
     ImageType4DF::IndexType inIdx;

     // output fmri
     ImageType4DF::RegionType outRegion;
     ImageType4DF::SizeType outSize = inSize;
     
     if (outSize[3] >= tslength) {
	  outSize[3] = tslength;
     }
     else {
	  printf("input time series length even larger than input fmri length. \nI just use the input time series length.\n");
     }
     outRegion.SetSize(outSize);
     ImageType4DF::IndexType outIdx;
     outIdx.Fill( 0 );
     outRegion.SetSize( outSize );
     outRegion.SetIndex(outIdx);
     ImageType4DF::Pointer outPtr = ImageType4DF::New();
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(0);


     // define random generator.
     boost::random::mt19937 gen(seed);
     boost::random::uniform_int_distribution<> roll_die(0, inSize[3]-1);

     // sample[i] = j means output's i'th point would be equal to input's j'th
     // time point.
     std::vector<unsigned> sample(outSize[3]);

     for (unsigned i = 0; i < sample.size(); i ++) {
	  sample[i] = roll_die(gen);
	  if (verbose >= 1)   printf("sample[%d] = %i\n", i, sample[i]);
     }

     IteratorType4DFloat outIt(outPtr, outPtr->GetLargestPossibleRegion() );

     for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++ outIt) {
	  outIdx = outIt.GetIndex();
	  inIdx[0] = outIdx[0];
	  inIdx[1] = outIdx[1];
	  inIdx[2] = outIdx[2];
	  inIdx[3] = sample[outIdx[3]];

	  outIt.Set ( inPtr->GetPixel(inIdx) );
     }

     WriterType4DFloat::Pointer writer = WriterType4DFloat::New();
     writer->SetInput( outPtr );
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
     std::cout << "fmriresample(): File  " << outfile << " saved.\n";

     return 0;
}

