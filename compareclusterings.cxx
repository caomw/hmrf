#include <common.h>

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     std::string refimage, inputimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Compare two clusterings.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("ref,r", po::value<std::string>(&refimage)->default_value("reference.nii.gz"), 
	   "reference image.")
	   ("input,i", po::value<std::string>(&inputimage)->default_value("image.nii.gz"), 
	    "test image.");
	       
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

     ReaderType3DShort::Pointer refReader = ReaderType3DShort::New();
     refReader->SetFileName(refimage);
     ImageType3DShort::Pointer refPtr = refReader->GetOutput();
     refPtr->Update();

     ReaderType3DShort::Pointer inputImageReader = ReaderType3DShort::New();
     inputImageReader->SetFileName(inputimage);
     ImageType3DShort::Pointer inputImagePtr = inputImageReader->GetOutput();
     inputImagePtr->Update();

     // Compute the maximum intensity value of input image.
     typedef itk::MinimumMaximumImageCalculator <ImageType3DShort>  ImageCalculatorFilterType;
     ImageCalculatorFilterType::Pointer inputImageCalculatorFilter = ImageCalculatorFilterType::New ();
     inputImageCalculatorFilter->SetImage(inputImagePtr);
     inputImageCalculatorFilter->Compute();
     unsigned numClustersInputImage = inputImageCalculatorFilter->GetMaximum();

     ImageCalculatorFilterType::Pointer refCalculatorFilter = ImageCalculatorFilterType::New ();
     refCalculatorFilter->SetImage(refPtr);
     refCalculatorFilter->Compute();
     unsigned numClustersRef = refCalculatorFilter->GetMaximum();

     // init contingency table.
     unsigned U = numClustersInputImage;
     unsigned V = numClustersRef;
     vnl_matrix <double> conTable(U, V, 0);

     IteratorType3DShort inputIt(inputImagePtr, inputImagePtr->GetLargestPossibleRegion() );
     IteratorType3DShort refIt(refPtr, refPtr->GetLargestPossibleRegion() );


     unsigned N = 0; 
     for (inputIt.GoToBegin(), refIt.GoToBegin(); !inputIt.IsAtEnd(); ++ inputIt, ++ refIt) {
	  if (inputIt.Get() > 0 && refIt.Get() > 0) {
	       conTable[inputIt.Get()-1][refIt.Get()-1] ++;
	       N ++;
	  }
     }

     // Ranjith Unnikrishnan, "measures of similarity"
     vnl_vector<double> udot(V,0);
     vnl_vector<double> dotv(U,0);
     for (unsigned u = 0; u < U; u++) {
	  for (unsigned v = 0; v < V; v ++) {
	       udot(u) += conTable(u,v);
	  }
     }

     for (unsigned v = 0; v < V; v++) {
	  for (unsigned u = 0; u < U; u ++) {
	       dotv(v) += conTable(u,v);
	  }
     }

     
     double ri = 1 - (0.5 * (udot.squared_magnitude() + dotv.squared_magnitude() ) - pow(conTable.array_two_norm(), 2) ) / (double(N) * (double(N)-1) / 2);

     // printf("udot.squared_magnitude(): %f\n", udot.squared_magnitude() );
     // printf("dotv.squared_magnitude(): %f\n", dotv.squared_magnitude() );
     // printf("conTable.array_two_norm(): %f\n", conTable.array_two_norm() );
     // printf("(N * (N-1) / 2: %f\n", (double(N) * (double(N)-1) / 2) );
     // printf("N = %d\n", N);
     if(verbose >=1 ) {
	  printf("Num of clusters, input image: %i. ref image: %i\n", numClustersInputImage, numClustersRef);
     }
     
     if(verbose >=1 ) {
	  printf("RI: %2.2f%%.   ", ri*100);
     }


     // following is from "Comparing Partitions" of Lawrence Hubert. The
     // Unnikrishnan's paper seems wrong for ajusted Rand index.
     double Sudot = 0, Sdotv = 0;
     for (unsigned u = 0; u < U; u++) {
	  Sudot += udot(u) * (udot(u) - 1) / 2;
     }

     for (unsigned v = 0; v < V; v ++) {
	  Sdotv += dotv(v) * (dotv(v) - 1) / 2;
     }

     double E = Sudot * Sdotv / ( N*(N-1)/2 );
     double Suv = 0;
     for (unsigned u = 0; u < U; u++) {
	  for (unsigned v = 0; v < V; v ++) {
	       Suv += conTable(u, v) * (conTable(u, v) - 1) / 2;
	  }
     }
     
     double ari = (Suv - E) / ( 0.5 * (Sudot + Sdotv) - E);
     // printf("Sudot = %f, Sdov = %f, E = %f, Suv = %f\n", Sudot, Sdotv, E, Suv);
     if (verbose >= 1) {
	  printf("ARI: %2.2f%%   ", ari*100);
     }

     // Nguyen Xuan Vinh, Information Theoretic Measurs for Clustering Comparison...
     double HU = 0;
     for (unsigned u = 0; u < U; u++) {
	  if (udot(u) > 0) {
	       HU -= (udot(u) / N) * log (udot(u) / N);
	  }
     }

     double HV = 0;
     for (unsigned v = 0; v < V; v++) {
	  if (dotv(v) > 0) {
	       HV -= (dotv(v) / N) * log (dotv(v) / N);
	  }
     }
     
     double HUV = 0;
     for (unsigned u = 0; u < U; u++) {
	  for (unsigned v = 0; v < V; v ++) {
	       if (conTable(u,v) > 0) {
		    HUV -= ( conTable(u,v) / N ) * log ( conTable(u,v) / N );
	       }
	  }
     }

     double HUgivenV = 0;
     for (unsigned u = 0; u < U; u++) {
	  for (unsigned v = 0; v < V; v ++) {
	       if (conTable(u,v) > 0) {
		    HUgivenV -= ( conTable(u,v) / N ) * log ( (conTable(u,v) / N) / (dotv(v) / N) );
	       }
	  }
     }

     double IUV = 0;
     for (unsigned u = 0; u < U; u++) {
	  for (unsigned v = 0; v < V; v ++) {
	       if (conTable(u,v) > 0) {
		    IUV += ( conTable(u,v) / N ) * log ( (conTable(u,v) / N) / (udot(u) * dotv(v) / (N*N) ) );
	       }
	  }
     }

     double mni_joint = IUV / HUV;
     double mni_max = IUV / (HU>HV?HU:HV);
     double mni_sum = 2 * IUV / (HU + HV);
     double mni_sqrt = IUV / sqrt(HU*HV);
     double mni_min = IUV / (HU<HV?HU:HV);
     
     // printf("HU = %f, HV = %f, HUV = %f, HUgivenV = %f, IUV = %f\n", HU, HV, HUV, HUgivenV, IUV);
     if (verbose >= 1) {
	  printf("mni_joint = %.2f%%, mni_max = %.2f%%, mni_sum = %.2f%%, mni_sqrt = %.2f%%, mni_min = %.2f%%\n", mni_joint*100, mni_max*100, mni_sum*100, mni_sqrt*100, mni_min*100);
     }

     if (verbose == 0) {
	  printf("%1.3f", ri);
     }

}
