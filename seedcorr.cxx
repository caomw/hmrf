#include <common.h>
#include <utility.h>
using namespace lemon;
int main(int argc, char* argv[])
{
     std::string inFile, outFile, maskFile;
     unsigned seedx = 0, seedy = 0, seedz = 0, radius = 0;
     
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")
	  ("in,i", po::value<std::string>(&inFile)->default_value("infile.nii"),
	   "fMRI file for correlation. Must be zero-mean, and standard deviation 1.")
	  ("mask,m", po::value<std::string>(&maskFile)->default_value("maskfile.nii"),
	   "mask file. Only used for seed region computation. Not for correlation computation.")
	  ("xvoxel,x", po::value<unsigned>(&seedx)->default_value(0),
	   "The x (voxel) coordinate of seed voxel.")
	  ("yvoxel,y", po::value<unsigned>(&seedy)->default_value(0),
	   "The y (voxel) coordinate of seed voxel.")
	  ("zvoxel,z", po::value<unsigned>(&seedz)->default_value(0),
	   "The z (voxel) coordinate of seed voxel.")
	  ("radius,r", po::value<unsigned>(&radius)->default_value(0),
	   "The radius of the seed region. The voxels within r voxel with the center are averaged")

	  ("out,o", po::value<std::string>(&outFile)->default_value("outfile.nii.gz"),
	   "output file. Each voxel is the correlation with seed region.");
     

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: seedcorr [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     ImageType4DF::Pointer fmriPtr;
     ReaderType4DFloat::Pointer fmriReader = ReaderType4DFloat::New();
     fmriReader->SetFileName(inFile);
     fmriReader->Update();
     fmriPtr = fmriReader->GetOutput();
     ImageType4DF::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();
     ImageType4DF::IndexType fmriIdx;

     ImageType3DChar::Pointer maskPtr;
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskFile);
     maskReader->Update();
     maskPtr = maskReader->GetOutput();
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     ImageType3DChar::IndexType maskIdx;
     maskIdx.Fill(0);

     ImageType3DFloat::Pointer corrPtr = ImageType3DFloat::New();
     ImageType3DFloat::RegionType corrRegion;
     corrRegion.SetIndex(maskIdx);
     corrRegion.SetSize(maskSize);
     corrPtr->SetRegions(corrRegion);
     corrPtr->Allocate();
     corrPtr->FillBuffer( 0 );
     IteratorType3DFloat corrIt(corrPtr, corrPtr->GetLargestPossibleRegion() );
     ImageType3DFloat::IndexType corrIdx;

     // get average mean time series for seed region.
     // IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() ) ;
     vnl_vector<double> seedAvg(fmriSize[3], 0);
     unsigned numPts = 0;
     for (maskIdx[0] = seedx - radius; maskIdx[0] <= seedx + radius; maskIdx[0] ++) {
	  for (maskIdx[1] = seedy - radius; maskIdx[1] <= seedy + radius; maskIdx[1] ++) {
	       for (maskIdx[2] = seedz - radius; maskIdx[2] <= seedz + radius; maskIdx[2] ++) {
		    if ( (pow(maskIdx[0] - seedx, 2) + pow(maskIdx[1] - seedy, 2) + pow(maskIdx[2] - seedz, 2) <= radius * radius ) && maskPtr->GetPixel(maskIdx) > 0) {
			 fmriIdx[0] = maskIdx[0];
			 fmriIdx[1] = maskIdx[1];
			 fmriIdx[2] = maskIdx[2];
			 numPts ++;
			 for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
			      seedAvg[fmriIdx[3]] += fmriPtr->GetPixel(fmriIdx);
			 }
		    }
	       } // z
	  } // y
     } // x

     // get average time series.
     seedAvg = seedAvg / numPts;
     
     // make sure standard deviation is one.
     seedAvg = seedAvg /  sqrt( seedAvg.squared_magnitude() / (fmriSize[3] - 1) );

     // Now compute correlation.
     double corr = 0;
     vnl_vector<double> thists(fmriSize[3], 0);
     for (corrIt.GoToBegin(); !corrIt.IsAtEnd(); ++ corrIt) {
	  corrIdx = corrIt.GetIndex();
	  fmriIdx[0] = corrIdx[0];
	  fmriIdx[1] = corrIdx[1];
	  fmriIdx[2] = corrIdx[2];
	  for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
	       thists[fmriIdx[3]] = fmriPtr->GetPixel(fmriIdx);
	  }
	  
	  thists = thists - thists.mean();
	  thists = thists / sqrt( thists.squared_magnitude() / (thists.size() - 1) );
	  corr = inner_product(thists, seedAvg) / (thists.size() - 1);
	  corrIt.Set( corr );
     }

     // copy header
     corrPtr->SetOrigin( maskPtr->GetOrigin() );
     corrPtr->SetSpacing(maskPtr->GetSpacing() );
     corrPtr->SetDirection(maskPtr->GetDirection() );

     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
	  
     writer->SetInput(corrPtr);
     writer->SetFileName(outFile);
     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  // return EXIT_FAILURE;
	  exit(1);
     } 

     std::cout << "seedcorr(): File " << outFile << " saved.\n";

     return 0;
}

	       




			 
     



     
