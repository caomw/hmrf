#include <common.h>
using namespace lemon;

int SaveGraphToImage(lemon::SmartGraph & theGraph, 
		     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
		     lemon::SmartGraph::NodeMap<std::vector<unsigned char> > & labelMap,
		     ImageType3DChar::Pointer maskPtr,
		     std::string outFile)
{
     unsigned M = labelMap[theGraph.nodeFromId( 0 )].size();

     ImageType3DChar::IndexType outImageIdx;
     outImageIdx.Fill(0);
     ImageType3DChar::SizeType outImageSize = maskPtr->GetLargestPossibleRegion().GetSize();;

     ImageType3DChar::RegionType outImageRegion;
     outImageRegion.SetSize(outImageSize);
     outImageRegion.SetIndex(outImageIdx);
     ImageType3DChar::Pointer outImagePtr = ImageType3DChar::New();
     outImagePtr->SetRegions(outImageRegion);
     outImagePtr->Allocate();
     outImagePtr->FillBuffer(0);
     
     // write labels (plus 1) at each node to image.
     for (lemon::SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=lemon::INVALID; ++ nodeIt) {
	  outImageIdx = coordMap[nodeIt];
	  outImagePtr->SetPixel(outImageIdx, labelMap[nodeIt][M-1]+1);
     }

     WriterType3DChar::Pointer writer = WriterType3DChar::New();
     writer->SetInput( outImagePtr );
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

     std::cout << "SaveGraphToImage: File " << outFile << " saved.\n";

     return 0;
}


int SaveSamplesToImage(lemon::SmartGraph & theGraph, 
		       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
		       lemon::SmartGraph::NodeMap<std::vector<unsigned char> > & labelMap,
		       ImageType3DChar::Pointer maskPtr,
		       std::string outFile)
{
     unsigned M = labelMap[theGraph.nodeFromId( 0 )].size();

     ImageType4DS::IndexType outImageIdx;
     outImageIdx.Fill(0);
     ImageType4DS::SizeType outImageSize;
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     outImageSize[0] = maskSize[0];
     outImageSize[1] = maskSize[1];
     outImageSize[2] = maskSize[2];
     outImageSize[3] = M;

     ImageType4DS::RegionType outImageRegion;
     outImageRegion.SetSize(outImageSize);
     outImageRegion.SetIndex(outImageIdx);
     ImageType4DS::Pointer outImagePtr = ImageType4DS::New();
     outImagePtr->SetRegions(outImageRegion);
     outImagePtr->Allocate();
     outImagePtr->FillBuffer(0);
     
     // write labels (plus 1) at each node to image.
     ImageType3DChar::IndexType imageIdx;
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  imageIdx = coordMap[nodeIt];
	  outImageIdx[0] = imageIdx[0];
	  outImageIdx[1] = imageIdx[1];
	  outImageIdx[2] = imageIdx[2];

	  for (unsigned mcIdx = 0; mcIdx < M; mcIdx ++) {
	       outImageIdx[3] = mcIdx;
	       outImagePtr->SetPixel(outImageIdx, labelMap[nodeIt][mcIdx] + 1);
	  }
     }

     WriterType4DS::Pointer writer = WriterType4DS::New();
     writer->SetInput( outImagePtr );
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

     std::cout << "SaveGraphToImage: File " << outFile << " saved.\n";

     return 0;
}

double logBesselI(float nu, double x)
{
     // From
     // http://people.kyb.tuebingen.mpg.de/suvrit/work/progs/movmf/. See
     // matlab code logbesseli.m, and
     // http://people.math.sfu.ca/~cbm/aands/page_378.htm (eq 9.7.7)
     double const Pi = 4 * atan(1);

     double frac = x/nu;
     double square = 1 + frac * frac;
     double root = sqrt(square);
     double eta = root + log(frac) - log(1+root);
     double logb = - log(sqrt(2 * Pi * nu)) + nu * eta - 0.25 * log(square);
     return (logb);
}
