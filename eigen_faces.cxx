/*
This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.

eigen_faces.cxx   ---  ECS 189H Homework 1 -------  January 7, 2007
This program uses ITK components to implement the EigenFaces algorithm for modeling the appearance of human faces in photographs.

Written by Owen Carmichael, January 2007

*/
// Include files for ITK content:

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImagePCAShapeModelEstimator.h"
#include "itkResampleImageFilter.h"
#include "itkImagePCADecompositionCalculator.h"
#include "itkAddImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkMinimumImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkIndent.h"
#include "itkImageDuplicator.h"
#include "itkCollectionOfPatchesHandLabeledIterator.h"
#include "itkAffineTransform.h"

// Convenience typedefs for ITK classes:
typedef unsigned char PixelType;
typedef double FloatPixelType;
typedef itk::Image<PixelType,2> InputImageType;
typedef itk::Image<FloatPixelType,2> FloatImageType;
typedef itk::CollectionOfPatchesHandLabeledIterator< InputImageType , FloatImageType > PatchesIteratorType;
typedef itk::ResampleImageFilter< FloatImageType, FloatImageType, FloatPixelType > PatchResamplerType;
typedef itk::ImagePCAShapeModelEstimator<FloatImageType,FloatImageType> PCAEstimatorType;
typedef itk::ImagePCADecompositionCalculator<FloatImageType,FloatImageType> PCAProjectionCalculatorType;
typedef itk::ImageFileWriter< FloatImageType > ImageWriterType;
typedef itk::Functor::IntensityLinearTransform< FloatImageType::PixelType, FloatImageType::PixelType > ImageRescalerFunctorType;
typedef itk::UnaryFunctorImageFilter< FloatImageType, FloatImageType,ImageRescalerFunctorType > ImageRescalerType;
typedef itk::AddImageFilter< FloatImageType, FloatImageType, FloatImageType > ImageAdderType;
typedef itk::SquaredDifferenceImageFilter< FloatImageType, FloatImageType, FloatImageType > SquaredDifferenceCalculatorType;
typedef itk::StatisticsImageFilter< FloatImageType > ImageStatisticsCalculatorType;
typedef itk::Functor::Minimum<float,float,float> FMinType;
typedef itk::Functor::Maximum<float,float,float> FMaxType;
typedef itk::Functor::Cast< float, int > F2ICast;
typedef itk::ImageDuplicator<FloatImageType> ImagePatchDuplicatorType;
typedef itk::RescaleIntensityImageFilter<FloatImageType,FloatImageType> IntensityNormalizerType;
typedef itk::ImageFileWriter< InputImageType > InputImageWriterType;
typedef itk::CastImageFilter< FloatImageType, InputImageType > Float2IntCasterType;

// Use these #defines to control whether or not the program prints out the
// original image patches, the resized image patches, the eigenfaces, and the
// reprojections to image files as it executes.  Comment them out to
// prevent one or the other type of output from being produced.  Note that
// printing the repojections can take up a lot of time and disk space...
#define EIGENFACES_OUTPUT_INPUT_PATCHES 1
#define EIGENFACES_OUTPUT_RESIZED_PATCHES 1
#define EIGENFACES_OUTPUT_EIGENFACES 1
//#define EIGENFACES_OUTPUT_REPROJ 1

#define IMAGE_EXT "png"

int main(int argc,char* argv[]) {
  try {  // This allows the program to catch exceptions thrown by ITK code

    char *image_path=argv[1]; // Directory of face images
    int eigenimage_size_x=32; // Face image patches are re-sized to this size
    int eigenimage_size_y=32; //   prior to modeling

    FMinType fmin;            // simple functions that compute the min and max
    FMaxType fmax;            //   of two floats
    F2ICast f2icast;          // casts a float to an int

    PatchResamplerType::Pointer patchResampler = PatchResamplerType::New();
    PCAEstimatorType::Pointer pcaEstimator = PCAEstimatorType::New();

    // create output folder
    system( "rm -rf out/" );
    system( "mkdir out" );


    // The PatchResampler takes an arbitrarily-sized image patch with a face in
    //   it, and re-sizes the patch so that it is eigenimage_size_x by
    //   eigenimage_size_y pixels.  Here we tell it what size the output image
    //   patches should be:

    FloatImageType::SizeType eigenimage_size;
    eigenimage_size[0]=eigenimage_size_x;  eigenimage_size[1]=eigenimage_size_y;
    patchResampler->SetSize( eigenimage_size );

    // This is the type of geometric transformation it should use to re-size
    //   the image patch:
    typedef itk::AffineTransform< FloatPixelType , 2 > TransformType;
    TransformType::Pointer transform = TransformType::New();

    // When an input image patch is re-sized into a new image,
    //   the input image pixels usually do not land exactly on the locations of
    //   the new pixels.  The LinearInterpolateImageFunction fills in the values
    //   of the new pixels based on where the old pixels landed after re-sizing.
    //   (We will talk a lot more about interpolation later on in the course...)

    typedef itk::LinearInterpolateImageFunction< FloatImageType, FloatPixelType > InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    patchResampler->SetInterpolator( (itk::InterpolateImageFunction< FloatImageType, FloatPixelType >*) interpolator.GetPointer() );

    // The PatchesIterator goes through the file called ground_truth in the
    //   image directory one line at a time.  Each line describes the location
    //   of an image patch that has a face in it.  Calling Get() causes the
    //   image to be read from file, and the face region to be extracted and
    //   returned.

    PatchesIteratorType::Pointer patchesIt = PatchesIteratorType::New();
    patchesIt->SetImagePath( image_path );  // Point the iterator to the directory with the images in it
    int num_patches=patchesIt->GetNumberOfPatches(); // number of face patches in the ground truth file

    // The patches are read in from file, re-sized, and fed as input to the
    // PCAEstimator.  Here we tell the PCAEstimator how many input image patches
    // to expect and how many principal components we expect out of it:

    pcaEstimator->SetNumberOfTrainingImages( num_patches );
    pcaEstimator->SetNumberOfPrincipalComponentsRequired( f2icast(fmin(num_patches,eigenimage_size_x*eigenimage_size_y)) );

    // These classes are for printing the principal components (eigenimages)
    // and re-sized image patches, which are floating point images, out to
    // image files.  The IntensityNormalizer shifts the values in the images to
    // lie between 0 and 255, the Float2IntCaster typecasts the float images
    // to integer type, and the InputImagerWriter writes them out to file.

    IntensityNormalizerType::Pointer rescaler = IntensityNormalizerType::New();
    Float2IntCasterType::Pointer f2icaster = Float2IntCasterType::New();
    InputImageWriterType::Pointer inputImageWriter = InputImageWriterType::New();

    int numberOfRemovedPatches = 0;
    // Loop through all of the patches described in the ground truth file:

    for(patchesIt->GoToBegin();!patchesIt->IsAtEnd();++(*patchesIt)) {
      // Calling Get() on the iterator is like dereferencing a pointer; it
      // gives you a pointer to the current image patch.

      if( NULL == patchesIt->Get() )
      {
          // update the number of removed patches
          numberOfRemovedPatches++;
          // upadte the total number of patches because the estimator will need to be informed of the change
          num_patches--;

          // skip this iteration because this patch was decided to be skipped
          continue;
      }

      FloatImageType::Pointer current_patch=patchesIt->Get();

#ifdef EIGENFACES_OUTPUT_INPUT_PATCHES
      // This code outputs the original image patches as fetched by the
      // PatchesIterator.
      rescaler->SetInput( current_patch );
      rescaler->SetOutputMinimum( 0 );    rescaler->SetOutputMaximum( 255 );
      f2icaster->SetInput ( rescaler->GetOutput() );
      inputImageWriter->SetInput( f2icaster->GetOutput() );
      char original_patch_filename[1024];
      sprintf(original_patch_filename,"out/original_patch%03d.%s",patchesIt->GetCurrentPatchNum(), IMAGE_EXT );
      inputImageWriter->SetFileName(original_patch_filename );
      inputImageWriter->Update();
#endif // #ifdef EIGENFACES_OUTPUT_INPUT_PATCHES

      // This is one of the clunky things about ITK: to get an image's
      //  dimensions, you have to call GetLargestPossibleRegion(), which gives
      //  you a Region covering the entire image, and then call GetSize() on
      //  that region.
      FloatImageType::SizeType patch_size=current_patch->GetLargestPossibleRegion().GetSize();

      // The dimensions of the image patch start out as x1-by-y1, and we want
      // to re-size it to be x2-by-y2, so we scale the dimensions by x2/x1 in
      // the x direction and y2/y1 in the y direction.

      TransformType::OutputVectorType scale;
      scale[0]=((float) patch_size[0]-1) / ((float) eigenimage_size_x);
      scale[1]=((float) patch_size[1]-1) / ((float) eigenimage_size_y);
      transform->SetIdentity();  transform->Scale( scale );

      // Set the current image patch as the input to the PatchResampler and
      // set its transform so it does the rescaling

      patchResampler->SetInput( current_patch );
      patchResampler->SetTransform( (itk::Transform< FloatPixelType , 2, 2 >*) transform.GetPointer() );

      // calling Update() causes the PatchResampler to perform the resizing.

      patchResampler->Update();

 #ifdef EIGENFACES_OUTPUT_RESIZED_PATCHES
      // This code writes the resized patch out to a file called patchX.tif,
      //  where this is the X'th patch in the ground_truth file.  First the
      //  intensities in the patch are rescaled to a 0-to-255 range by the
      //  IntensityRescaler, then cast to an integer data type, then written
      // out to file by the ImageWriter.

      rescaler->SetInput( patchResampler->GetOutput() );
      rescaler->SetOutputMinimum( 0 );    rescaler->SetOutputMaximum( 255 );
      f2icaster->SetInput ( rescaler->GetOutput() );
      inputImageWriter->SetInput( f2icaster->GetOutput() );
      char resized_patch_filename[1024];
      sprintf(resized_patch_filename,"out/resized_patch%03d.%s",patchesIt->GetCurrentPatchNum(), IMAGE_EXT );
      inputImageWriter->SetFileName( resized_patch_filename );
      inputImageWriter->Update();
#endif //  #ifdef EIGENFACES_OUTPUT_RESIZED_PATCHES

      // The PatchResampler has a static output buffer to which it writes the
      //  resampled patches, meaning that the 1st resampled patch will be
      //  overwritten by the 2nd one when Update() is called the 2nd time.
      //  Therefore, we use the ImagePatchDuplicator to make a copy of the
      //  resized image patch.
      ImagePatchDuplicatorType::Pointer dup = ImagePatchDuplicatorType::New();

      dup->SetInputImage( patchResampler->GetOutput() );
          dup->Update();

      //  It might seem odd to allocate a new
      //  ImagePatchDuplicator in the middle of the for loop: won't we lose
      //  track of our pointer to it after this iteration of the loop
      //  completes?  Usually, yes, but remember that we are doing memory
      //  management by reference counting.  Specifically, in the following
      //  line the PCAEstimator calls dup->GetOutput() to grab a pointer to the
      //  ImagePatchDuplicator output buffer; therefore the number of references
      //  to the ImagePatchDuplicator output buffer is incremented to 1, thus
      //  insuring it won't be deallocated.  When
      //  the PCAEstimator is deallocated, its pointer to the
      //  ImagePatchDuplicator output buffer goes away, and the number of
      //  references to the ImagePatchDuplicator output buffer is decremented
      //  to 0, triggering its deallocation.

//      pcaEstimator->SetInput( patchesIt->GetCurrentPatchNum() , dup->GetOutput() );
      // modify the patch number because some of the patches have been skipped
      pcaEstimator->SetInput( patchesIt->GetCurrentPatchNum() - numberOfRemovedPatches, dup->GetOutput() );

    }

    // Redo the estimator setting with updated values
    pcaEstimator->SetNumberOfTrainingImages( num_patches );
    pcaEstimator->SetNumberOfPrincipalComponentsRequired( f2icast(fmin(num_patches,eigenimage_size_x*eigenimage_size_y)) );

    std::cout << numberOfRemovedPatches << " patches skipped based on my face recognition algorithm" << std::endl;

    //  Inputs to the PCAEstimator were set up in the previous for loop, so we
    //  just need to call Update() to get it to compute the principal components
    std::cout << std::endl;
    std::cout << "doing PCA estimation..."; std::cout.flush();
    // do the PCA estimation
    pcaEstimator->Update();
    std::cout << "done." << std::endl;

    // Uncomment this to dump ITK's information for this object
    // std::cout << pcaEstimator << std::endl;

#ifdef EIGENFACES_OUTPUT_EIGENFACES
  // This code writes the principal components (eigenfaces) to file after
  //   rescaling them to a 0-to-255 intensity range and casting them to "int."
  //   The file name for the Xth principal component is eigen_outX.tif
  //   The 0th output is the average of all input images

  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  std::cout << "writing PCA basis vectors to file..."; std::cout.flush();
  for (int i=0;i<fmin(num_patches,eigenimage_size_x*eigenimage_size_y);i++) {
    char eigenface_filename[1024];
    rescaler->SetInput( pcaEstimator->GetOutput(i)  );
    rescaler->SetOutputMinimum( 0 );    rescaler->SetOutputMaximum( 255 );
    f2icaster->SetInput ( rescaler->GetOutput() );
    inputImageWriter->SetInput( f2icaster->GetOutput() );
    sprintf(eigenface_filename,"out/eigen_out%03d.%s",i, IMAGE_EXT );
    inputImageWriter->SetFileName( eigenface_filename );
    inputImageWriter->Update();
  }
  std::cout << "done." << std::endl;
#endif // #ifdef EIGENFACES_OUTPUT_EIGENFACES

  // This segment computes the reconstruction error for the input images using
  //   different numbers of eigenfaces.  The PCA coefficients for each input
  //   image are computed using the PCAProjectionCalculator.  Then the input
  //   image is reconstructed using 1, 2, 3... of its PCA coefficients.  To do
  //   so, each eigenface is multiplied by the corresponding PCA coefficient
  //   using an ImageRescaler, and the scaled images for the first k eigenfaces
  //   are added together using an ImageAdder.  The difference between the
  //   original input image and the reconstructed one is evaluated by using
  //   a SquaredDifferenceCalculator to compute the squared difference in pixel
  //   intensity at each pixel location, and summing up all the squared
  //   differences with an ImageStatisticsCalculator.

  std::cout << "computing reconstruction error for training examples..."; std::cout.flush();
  // Allocate the PCAProjectionCalculator and pass it the eigenfaces from the
  //  PCAEstimator:
  PCAProjectionCalculatorType::Pointer pcaProjector = PCAProjectionCalculatorType::New();
  pcaProjector->SetBasisFromModel( pcaEstimator );
  // This is the number of eigenfaces output by the PCAEstimator:
  int ncomponents=pcaEstimator->GetNumberOfPrincipalComponentsRequired();

  // We sum the reprojection errors for K eigenfaces over all images, and
  // store the summed reprojection error for each K in reproj_error:
  float *reproj_error=(float*) malloc(ncomponents*sizeof(float));
  memset(reproj_error,0,sizeof(float)*ncomponents);

  // Allocate the modules to rescale the eigenfaces by the PCA coefficients,
  // add together rescaled eigenfaces, take the pixel-wise squared difference
  // in intensities between images, and sum up the squared differences over the
  // entire image:
  ImageRescalerType::Pointer imageRescaler = ImageRescalerType::New();
  ImageAdderType::Pointer imageAdder = ImageAdderType::New();
  SquaredDifferenceCalculatorType::Pointer squaredDiffer = SquaredDifferenceCalculatorType::New();
  ImageStatisticsCalculatorType::Pointer statsCalculator = ImageStatisticsCalculatorType::New();
  ImageRescalerFunctorType rescaleFunctor;

  // Loop over all image patches:
  for (int i=0;i<num_patches;i++) {
    // Get the PCA coefficients for this patch and store them in proj:

    pcaProjector->SetImage( pcaEstimator->GetInput(i) );
    pcaProjector->Compute();

    // Uncomment this to dump ITK's information for this object
    // std::cout<<pcaProjector<<std::endl;

    PCAProjectionCalculatorType::BasisVectorType proj=pcaProjector->GetProjection();

    // Scale the first eigenimage by the first PCA coefficient. The IntensityLinearTransform functor has bounds checking, so the minimum and maximum must be set:
    rescaleFunctor.SetFactor(proj[0]);    rescaleFunctor.SetOffset(0);
    rescaleFunctor.SetMinimum(-1e+100);    rescaleFunctor.SetMaximum(1e+100);
    imageRescaler->SetFunctor( rescaleFunctor );
    imageRescaler->SetInput( pcaEstimator->GetOutput(1) );
    imageRescaler->Update();

    // Initializes the output buffer of the ImageAdder to just have the first scaled eigenface in it:
    ImagePatchDuplicatorType::Pointer duplicator = ImagePatchDuplicatorType::New();

    // pcaEstimator Output 0 is the average image (average of all input images).
    // pcaProjector gives the projection (coefficients of the eigenvectors) that produces
    //  an input image when added to the average image. Thus, we add the average image first,
    //  then begin adding eigenfaces with the found coefficients.
    imageAdder->SetInput(0,pcaEstimator->GetOutput(0));
    imageAdder->SetInput(1,imageRescaler->GetOutput() );

    // compute the squared difference between the first scaled eigenface and
    //  the input image:
    squaredDiffer->SetInput( 0, imageAdder->GetOutput() );
    squaredDiffer->SetInput( 1, pcaEstimator->GetInput(i) );
    squaredDiffer->Update();
    // sum the squared differences over all pixels:
    statsCalculator->SetInput( squaredDiffer->GetOutput() );
    statsCalculator->Update();
    reproj_error[0]+=statsCalculator->GetSum();

#ifdef EIGENFACES_OUTPUT_REPROJ
    //     write reprojections to image files:

    char reproj_filename[1024];

    // write reprojection estimate with zero eigenfaces (the average input face)
    rescaler->SetInput( imageAdder->GetInput(0) );
    rescaler->SetOutputMinimum( 0 );    rescaler->SetOutputMaximum( 255 );
    f2icaster->SetInput ( rescaler->GetOutput() );
    inputImageWriter->SetInput( f2icaster->GetOutput() );
    sprintf(reproj_filename,"out/reproj.%03d.%03d.%s",i,0, IMAGE_EXT );
    inputImageWriter->SetFileName( reproj_filename );
    inputImageWriter->Update();

    // write reprojection estimate with 1 eigenface
    rescaler->SetInput( imageAdder->GetOutput() );
    rescaler->SetOutputMinimum( 0 );    rescaler->SetOutputMaximum( 255 );
    f2icaster->SetInput ( rescaler->GetOutput() );
    inputImageWriter->SetInput( f2icaster->GetOutput() );

    sprintf(reproj_filename,"out/reproj.%03d.%03d.%s",i,1, IMAGE_EXT );
    inputImageWriter->SetFileName( reproj_filename );
    inputImageWriter->Update();
#endif

    duplicator->SetInputImage( imageAdder->GetOutput() );

    // Loop over the remaining eigenfaces:
    for (int j=1;j<ncomponents;j++) {
      // scale the eigenface by the PCA coefficient:
      rescaleFunctor.SetFactor(proj[j]);
      imageRescaler->SetFunctor( rescaleFunctor );
      imageRescaler->SetInput( pcaEstimator->GetOutput(j+1) );
      imageRescaler->Update();
      duplicator->Update();
      FloatImageType* acc=duplicator->GetOutput();

      // add it to the current reprojection estimate:
      imageAdder->SetInput(0,acc);

      // implicitly, this is still set from having reused the object so it is not necessary:
      // imageAdder->SetInput(1,imageRescaler->GetOutput() );

      imageAdder->Update();
      // compute the difference with the input image:
      squaredDiffer->Update();
      // sum the squared differences
      statsCalculator->Update();

#ifdef EIGENFACES_OUTPUT_REPROJ
      //     write current reprojection estimate to an image file:
      rescaler->SetInput( imageAdder->GetOutput() );
      rescaler->SetOutputMinimum( 0 );    rescaler->SetOutputMaximum( 255 );
      f2icaster->SetInput ( rescaler->GetOutput() );
      inputImageWriter->SetInput( f2icaster->GetOutput() );
      sprintf(reproj_filename,"out/reproj.%03d.%03d.%s",i,j+1, IMAGE_EXT );
      inputImageWriter->SetFileName( reproj_filename );
      inputImageWriter->Update();
#endif // #ifdef EIGENFACES_OUTPUT_REPROJ
      // add squared difference for this image to the sum over all images:
      reproj_error[j]+=statsCalculator->GetSum();
    }
  }
  // write the reprojection errors to file:
  std::cout << "done." << std::endl;
  char reproj_error_filename[1024];
  sprintf(reproj_error_filename,"%s%s",image_path,"/reprojection_error");
  std::cout << "writing reprojection errors to file: " << reproj_error_filename << std::endl; std::cout.flush();
  std::ofstream reproj_error_file(reproj_error_filename);
  for (int i=0;i<ncomponents-1;i++) {
    // scale by the number of inputs so that we can compare errors across datasets of different size
    reproj_error[i] /= ncomponents;
    reproj_error_file << i+1 << "\t" << reproj_error[i] << std::endl;
  }
  std::cout << "done." << std::endl;
  }
  catch( itk::ExceptionObject & exp )
    {
      std::cerr << "Exception caught ! " << std::endl;
      std::cerr << exp << std::endl;
    }
}
