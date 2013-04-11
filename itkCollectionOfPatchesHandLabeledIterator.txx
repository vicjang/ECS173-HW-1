/*
This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.

itkCollectionOfPatchesHandLabeledIterator.txx   ---  ECS 189H Homework 1
January 7, 2007

This file contains the implementation of itk::CollectionOfPatchesHandLabeledIterator.  The class iterates over a list of image regions that are extracted from images in a directory.  The image regions to be extracted are specified in a text file called ground_truth in the directory.  The first line in ground_truth is a single number: the number of image regions listed in the file.  Each following line in ground_truth specifies a different face region in a particular image.  The format for one line of ground_truth is:

<image filename> <left eye x> <left eye y> <right eye x> <right eye y> <nose x> <nose y> <left corner of mouth x> <left corner of mouth y> <center of mouth x> <center of mouth y> <right corner of mouth x> <right corner of mouth y>

For example, "left eye x" and "left eye y" are the X and Y coordinates of the left eye of the face in the image called "image filename."  This class converts each region description into a rectangular image region by taking the bounding box of the face features and adding some extra pixels to the periphery of the bounding box ("padding").  How much padding to add around the face features is controlled by the SetCushionPctX() and SetCushionPctY() methods.

The class instantiates in itk::ImageFileReader to read the images and an itk::RegionOfInterestImageFilter to extract the specified regions.

Calling SetImagePath() triggers the reading of ground_truth in the specified directory, but regions themselves are not read in from the specified image files until Get() is called.  Changing the position of the iterator by calling GoToBegin(), operator++(), etc just changes which image region is the current one but don't cause the image to be read in from file.  Get() calls the ImageFileReader to read in the image, specifies the rectangular image region based on the entry in ground_truth,  calls the RegionOfInterestImageFilter to extract that region from the image, and returns the region to the user.

Currently, the iterator is implemented as a simple array of structs, each of which stores information on a particular image region, and an integer index that points to an element in the array.  Methods that change the iterator position (operator++(), operator--(), operator+=(), operator-=(), GoToBegin(), GoToPosition()) just change the value of the array index in obvious ways.  Being at the end of the list of regions-- being PAST the final region in the list, not AT the final region in the list-- is represented by making the array index one greater than the maximal index of image regions in the array.  If that's where the iterator is, IsAtEnd() returns true.  GoToEnd() puts the iterator there.
Written by Owen Carmichael, January 2007

*/

#ifndef _itkCollectionOfPatchesHandLabeledIterator_txx
#define _itkCollectionOfPatchesHandLabeledIterator_txx

// Include files for the class declaration and an ITK header that handles exceptions:
#include "itkCollectionOfPatchesHandLabeledIterator.h"
#include "itkExceptionObject.h"

namespace itk
{

  // The default constructor sets internal variables to obvious defaults and 
  // allocates the ImageFileReader and RegionOfInterest filter.  The default
  // padding around the face features is 20% in both directions.

  template<class TInputImage,class TOutputImage>	
  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::CollectionOfPatchesHandLabeledIterator() {
    patches=NULL;
    num_patches=0;
    min_patchnum=0;
    max_patchnum=0;
    current_patchnum=0;
    patches_is_allocated=false;
    imageReader=ImageReaderType::New();
    patchExtractor=PatchExtractorType::New();
    sprintf(image_path,".");
    cushion_pct_x=.20;    cushion_pct_y=.20;
  }
  
  // This internal function reads the image region descriptions in the ground_truth file and stores the image regions in the patches[] data structure.  It allocates patches[] as needed.  It is called by SetImagePath().

  template<class TInputImage,class TOutputImage>	
  void CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::ReadGroundTruthFile() {
    
    std::ifstream ground_truth_file;
    ground_truth_file.open(ground_truth_filename, std::ios::in);

    // Throw an exception if the ground_truth file doesn't exist:

    if( !ground_truth_file)
      {
        itk::ExceptionObject exception(__FILE__, __LINE__);
        char msg[1024];
        sprintf(msg,"File %s cannot be read.",ground_truth_filename);
        exception.SetDescription(msg);
        char loc[1024];
        sprintf(loc,"CollectionOfPatchesHandLabeledIterator::ReadGroundTruthFile(char *).");
        exception.SetLocation(loc);
        throw exception;
      }
    
    // The first line of the ground truth file tells us how many patches are in it:

    ground_truth_file >> num_patches;

    // Throw an exception if the specified number of patches is invalid:

    if (num_patches<=0) {
      itk::ExceptionObject exception(__FILE__, __LINE__);
      char msg[1024];
      sprintf(msg,"The reported number of patches in %s is %d",ground_truth_filename,num_patches);
      exception.SetDescription(msg);
      char loc[1024];
      sprintf(loc,"CollectionOfPatchesHandLabeledIterator::ReadGroundTruthFile(char *).");
      exception.SetLocation(loc);
      throw exception;
    }

    // If we are re-reading the ground_truth file, deallocate the previous image region description storage:
    if (patches_is_allocated) { delete patches; }

    // Allocate storage for the image region descriptions:
    patches=new patch_info[num_patches];
    patches_is_allocated=true;
    
    // These are the array indices for the start and end of the iterator:
    min_patchnum=0;  max_patchnum=num_patches;

    std::cout << "reading in " << num_patches << " image patches from " << ground_truth_filename << " : ";

    // Read in the entire file of image region descriptions:
    int patch_counter=0;
    while(ground_truth_file.peek()!=std::ifstream::traits_type::eof()) {
      // print out our progress every 10 lines:
      if (patch_counter%10==0) { std::cout << patch_counter << " ";  std::cout.flush(); }
      // Throw an exception if the file has more regions in it than it said in the first line:
      if (patch_counter>num_patches) { 
        itk::ExceptionObject exception(__FILE__, __LINE__);
        char msg[1024];
        sprintf(msg,"File %s claimed to have %d entries but appears to have more than that",ground_truth_filename,num_patches);
        exception.SetDescription(msg);
        char loc[1024];
        sprintf(loc,"CollectionOfPatchesHandLabeledIterator::CollectionOfPatchesHandLabeledIterator(char *)");
        exception.SetLocation(loc);
        throw exception;
      }

      // Read in the file name and coordinates for each face feature in the image region.  I should probably make sure that each of these are of a valid data type (floats), but I don't...

      std::string image_filename;
      ground_truth_file >> image_filename;
      strcpy(patches[patch_counter].image_filename,image_filename.c_str());
      ground_truth_file >> patches[patch_counter].left_eye_x; 
      ground_truth_file >> patches[patch_counter].left_eye_y; 
      ground_truth_file >> patches[patch_counter].right_eye_x;  
      ground_truth_file >> patches[patch_counter].right_eye_y;
      ground_truth_file >> patches[patch_counter].nose_x;       
      ground_truth_file >> patches[patch_counter].nose_y;   
      ground_truth_file >> patches[patch_counter].left_corner_mouth_x;
      ground_truth_file >> patches[patch_counter].left_corner_mouth_y;
      ground_truth_file >> patches[patch_counter].center_mouth_x; 
      ground_truth_file >> patches[patch_counter].center_mouth_y; 
      ground_truth_file >> patches[patch_counter].right_corner_mouth_x; 
      ground_truth_file >> patches[patch_counter].right_corner_mouth_y;
      patch_counter++;
    }
    std::cout << num_patches << " done" << std::endl; std::cout.flush();

    // Throw an exception if the first line in the file told us there were num_patches patches, but there were fewer than that many listed in the file.
    if (patch_counter<num_patches) { 
      itk::ExceptionObject exception(__FILE__, __LINE__);
      char msg[1024];
      sprintf(msg,"CollectionOfPatchesHandLabeledIterator::CollectionOfPatchesHandLabeledIterator(char *):  File %s claimed to have %d entries but appears to only have %d",ground_truth_filename,num_patches,patch_counter-1);
      exception.SetDescription(msg);
      throw exception;
    }

  }

  //  The destructor just deallocates patches[] if it has been allocated.

  template<class TInputImage,class TOutputImage>	
  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::~CollectionOfPatchesHandLabeledIterator() {
    if (patches_is_allocated) { delete patches; }
  }

  // Get() opens the currently-pointed-to image, extracts the requested region, and returns it.

  template<class TInputImage,class TOutputImage> 
  SmartPointer<TOutputImage>  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::Get() {

    // Make sure the patches[] data structure is allocated and the iterator actually points to a patch.  Note that this is a little different than the iterator being at a valid position, because the end of the array-- one past the last element-- is a valid position but there are no patches there.

    AssertAllocated();  AssertPatchPosition(current_patchnum,std::string("Get()"));

    // Convenience reference to the currently-pointed-to element in patches[]:

    patch_info& pi=patches[current_patchnum];

    // Read in the image:
    imageReader=ImageReaderType::New();
    char path_and_filename[2048];
    sprintf(path_and_filename,"%s/%s",image_path,pi.image_filename);
    imageReader->SetFileName ( path_and_filename );
    imageReader->Update();

    // Get the image size:

    SizeType image_size = imageReader->GetOutput()->GetLargestPossibleRegion().GetSize();

    // Calculate the coordinates of the rectangular region we will extract: take the bounding box around all face features in this region:

    float min_x=fmin(fmin(fmin(pi.left_eye_x,pi.right_eye_x),fmin(pi.nose_x,pi.left_corner_mouth_x)),fmin(pi.center_mouth_x,pi.right_corner_mouth_x));
    float max_x=fmax(fmax(fmax(pi.left_eye_x,pi.right_eye_x),fmax(pi.nose_x,pi.left_corner_mouth_x)),fmax(pi.center_mouth_x,pi.right_corner_mouth_x));
    float min_y=fmin(fmin(fmin(pi.left_eye_y,pi.right_eye_y),fmin(pi.nose_y,pi.left_corner_mouth_y)),fmin(pi.center_mouth_y,pi.right_corner_mouth_y));
    float max_y=fmax(fmax(fmax(pi.left_eye_y,pi.right_eye_y),fmax(pi.nose_y,pi.left_corner_mouth_y)),fmax(pi.center_mouth_y,pi.right_corner_mouth_y));
  
    // Calculate how many pixels of buffer to add around the periphery of the bounding box: 

    float cushion_x=cushion_pct_x*(max_x-min_x);    float cushion_y=cushion_pct_y*(max_y-min_y);

    // Add the X% buffer to the bounding box.  The image region we will extract has (patch_start[0],patch_start[1]) as its upper-lefthand corner, and it stretches to the right by patch_size[0] pixels and downward by patch_size[1] pixels

    IndexType patch_start;
    patch_start[0] = f2icast(fmax(min_x-cushion_x,0));
    patch_start[1] = f2icast(fmax(min_y-cushion_y,0));
    SizeType patch_size;    
    patch_size[0] = f2icast(fmin(max_x-min_x + 2*cushion_x,image_size[0]-min_x));
    patch_size[1] = f2icast(fmin(max_y-min_y + 2*cushion_y,image_size[1]-min_y));  

    // Set the output of the ImageFileReader as the input to the RegionOfInterest filter: 

    patchExtractor->SetInput( imageReader->GetOutput() );

    // Tell it what region to extract:

    RegionType patch;
    patch.SetIndex( patch_start );    patch.SetSize( patch_size );
    patchExtractor->SetRegionOfInterest( patch );

    // Press "Go":

    patchExtractor->Update();

    // Get the output of the RegionOfInterest filter:

    PatchExtractorOutputType* patch_image = patchExtractor->GetOutput();

    // Re-set the image region's coordinate frame to have (0,0) as its origin
    // (This is a bit obscure-- don't worry about it...)

    typename PatchExtractorOutputType::PointType origin;
    origin[0]=0; origin[1]=0;
    patch_image->SetOrigin( origin );

    return patch_image;
  }

  // Set the current patch index to the first one in the array:
  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::GoToBegin() { 
    //    Make sure patches[] is allocated
    AssertAllocated(); 
    current_patchnum=min_patchnum;
  }

  // Set the current patch index to ONE PAST the last one in the array:

  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::GoToEnd() {
    //    Make sure patches[] is allocated
    AssertAllocated(); 
    current_patchnum=max_patchnum;
  }

  // True if the current patch index is one past the last one in the array:

  template<class TInputImage,class TOutputImage> 
  bool  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::IsAtEnd() {
    //    Make sure patches[] is allocated and the array index is at a valid position
    AssertAllocated();  AssertValidPosition(current_patchnum,std::string("IsAtEnd()")); 
    return current_patchnum==max_patchnum;
  }

  // True if the current patch index is at the first one in the array:

  template<class TInputImage,class TOutputImage> 
  bool  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::IsAtBegin() { 
    //    Make sure patches[] is allocated and the array index is at a valid position
    AssertAllocated();  AssertValidPosition(current_patchnum,std::string("IsAtBegin()"));
    return current_patchnum==min_patchnum;
  }

  // Increments the region index to point to the next one in the array

  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::operator++() {
    //    Make sure patches[] is allocated and the array index is at a valid position
    AssertAllocated();  AssertValidPosition(current_patchnum,std::string("operator++()")); 
    if (IsAtEnd()) {      
      itk::ExceptionObject exception(__FILE__, __LINE__);
      char msg[1024];
      sprintf(msg,"CollectionOfPatchesHandLabeledIterator::operator++():  Iterator is already at the end.");
      exception.SetDescription(msg);
      throw exception;
    }
    current_patchnum++;
  }
  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::operator--() {
    //    Make sure patches[] is allocated and the array index is at a valid position
    AssertAllocated();  AssertValidPosition(current_patchnum,std::string("operator--()")); 
    if (IsAtBegin()) {
      itk::ExceptionObject exception(__FILE__, __LINE__);
      char msg[1024];
      sprintf(msg,"CollectionOfPatchesHandLabeledIterator::operator--():  Iterator is already at the beginning.");
      exception.SetDescription(msg);
      throw exception;
    }
    current_patchnum--;
  }
  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::operator+=(int increment) {
    //    Make sure patches[] is allocated and the array index is at a valid position
    AssertAllocated();  AssertValidPosition(current_patchnum,std::string("operator+=()"));
    AssertValidPosition(current_patchnum+increment,std::string("operator+=()")); 
    current_patchnum+=increment;
  }

  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::operator-=(int increment) {
    //    Make sure patches[] is allocated and the array index is at a valid position
    AssertAllocated();  AssertValidPosition(current_patchnum,std::string("operator-=()")); 
    AssertValidPosition(current_patchnum-increment,std::string("operator-=()")); 
    current_patchnum-=increment;
  }

  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::SetPosition( int new_patchnum) {
    //    Make sure patches[] is allocated and the array index is at a valid position
    AssertAllocated();  AssertValidPosition(new_patchnum,std::string("SetPosition( int )"));
    current_patchnum=new_patchnum;
  }
  
  // Throw an exception if patches[] is not allocated:

  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::AssertAllocated() {
    if (!patches_is_allocated) {
      itk::ExceptionObject exception(__FILE__, __LINE__);
      char msg[1024];
      sprintf(msg,"Invalid iterator; no patches have been read in from the ground truth file.");
      exception.SetDescription(msg);
      char loc[1024];
      sprintf(loc,"CollectionOfPatchesHandLabeledIterator::Get()");
      exception.SetLocation(loc);
      throw exception;
    }
  }

  // Throw an exception if the iterator is not set to a valid position:

  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::AssertValidPosition(int pos,std::string fname) { 
    if (!(pos <= max_patchnum && pos>=min_patchnum)) {
      itk::ExceptionObject exception(__FILE__, __LINE__);
      char msg[1024];
      sprintf(msg,"Invalid iterator position");
      exception.SetDescription(msg);
      char loc[1024];
      sprintf(loc,"CollectionOfPatchesHandLabeledIterator::%s",fname.c_str());
      exception.SetLocation(loc);
      throw exception;
    }
  }

  // Throw an exception if the iterator is not set to point to an image region

  template<class TInputImage,class TOutputImage> 
  void  CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::AssertPatchPosition(int pos,std::string fname) { 
    if (!(pos < max_patchnum && pos>=min_patchnum)) {
      itk::ExceptionObject exception(__FILE__, __LINE__);
      char msg[1024];
      sprintf(msg,"Iterator not pointed at a patch");
      exception.SetDescription(msg);
      char loc[1024];
      sprintf(loc,"CollectionOfPatchesHandLabeledIterator::%s",fname.c_str());
      exception.SetLocation(loc);
      throw exception;
    }
  }

  // Specify the path to the images and ground_truth file.  This triggers the 
  // reading of the ground_truth file so that image region information is stored in memory.

  template<class TInputImage,class TOutputImage> 
  void CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::SetImagePath(char *new_image_path) {

    // Make sure the specified path is valid:

    if (new_image_path==NULL) {
      ExceptionObject  e(__FILE__, __LINE__); 
      e.SetLocation("CollectionOfPatchesHandLabeledIterator::SetImagePath(char *)");
      e.SetDescription("Image path is NULL"); 
      throw e; 
    }
    strcpy(image_path,new_image_path);
    sprintf(ground_truth_filename,"%s/%s",image_path,"ground_truth");

    // Read the file with image region descriptions in it:
    ReadGroundTruthFile();

    // The default constructor should have allocated the ImageReader and RegionOfInterest filter for us, but just in case...

    if (imageReader.IsNull()) { imageReader=ImageReaderType::New();  }
    if (patchExtractor.IsNull()) { patchExtractor=PatchExtractorType::New(); }
  }

  // Return how many image regions are in the ground_truth file

  template<class TInputImage,class TOutputImage> 
  int CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::GetNumberOfPatches() {
    //    Make sure patches[] is allocated
    AssertAllocated();
    return num_patches;
  }

  // Get the index of the current patch

  template<class TInputImage,class TOutputImage> 
  int CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::GetCurrentPatchNum() {
    //    Make sure patches[] is allocated and the array index is at a valid position
    AssertAllocated();  AssertValidPosition(current_patchnum,"GetCurrentPatchNum");
    return current_patchnum - min_patchnum;
  }

  // Get the padding percentage in the X direction

  template<class TInputImage,class TOutputImage> 
  float CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::GetCushionPctX() 
  {
    return cushion_pct_x;
  }

  // Get the padding percentage in the Y direction

  template<class TInputImage,class TOutputImage> 
  float CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::GetCushionPctY() 
  {
    return cushion_pct_y;
  }

  // Set the padding around the region bounding box in the X direction:

  template<class TInputImage,class TOutputImage> 
  void CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::SetCushionPctX(float cpx) 
  {
    cushion_pct_x=cpx;
  }

  // Set the padding around the region bounding box in the Y direction:

  template<class TInputImage,class TOutputImage> 
  void CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage>::SetCushionPctY(float cpy) 
  {
     cushion_pct_y=cpy;
  }
} 
 
#endif  // #ifndef _itkCollectionOfPatchesHandLabeledIterator_txx


