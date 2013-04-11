/*
This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.

itkCollectionOfPatchesHandLabeledIterator.h   ---  ECS 189H Homework 1
January 7, 2007

This file contains the class declaration for itk::CollectionOfPatchesHandLabeledIterator.  The class iterates over a list of image regions that are extracted from images in a directory.  The image regions to be extracted are specified in a text file called ground_truth in the directory.  The first line in ground_truth is a single number: the number of image regions listed in the file.  Each following line in ground_truth specifies a different face region in a particular image.  The format for one line of ground_truth is:

<image filename> <left eye x> <left eye y> <right eye x> <right eye y> <nose x> <nose y> <left corner of mouth x> <left corner of mouth y> <center of mouth x> <center of mouth y> <right corner of mouth x> <right corner of mouth y>

For example, "left eye x" and "left eye y" are the X and Y coordinates of the left eye of the face in the image called "image filename."  This class converts each region description into a rectangular image region by taking the bounding box of the face features and adding some extra pixels to the periphery of the bounding box ("padding").  How much padding to add around the face features is controlled by the SetCushionPctX() and SetCushionPctY() methods.  The default for both is .2 (i.e. the patches add a 20% border all around).  If adding the padding to the image patch makes it extend outside the bounds of the image, the padding is cropped off so that the image patch stops at the image boundary.

The class instantiates in itk::ImageFileReader to read the images and an itk::RegionOfInterestImageFilter to extract the specified regions.

Calling SetImagePath() triggers the reading of ground_truth in the specified directory, but regions themselves are not read in from the specified image files until Get() is called.  Changing the position of the iterator by calling GoToBegin(), operator++(), etc just changes which image region is the current one but don't cause the image to be read in from file.  Get() calls the ImageFileReader to read in the image, specifies the rectangular image region based on the entry in ground_truth,  calls the RegionOfInterestImageFilter to extract that region from the image, and returns the region to the user.

Implementation:  Currently, the iterator is implemented as a simple array of structs, each of which stores information on a particular image region, and an integer index that points to an element in the array.  Methods that change the iterator position (operator++(), operator--(), operator+=(), operator-=(), GoToBegin(), GoToPosition()) just change the value of the array index in obvious ways.  Being at the end of the list of regions-- being PAST the final region in the list, not AT the final region in the list-- is represented by making the array index one greater than the maximal index of image regions in the array.  If that's where the iterator is, IsAtEnd() returns true.  GoToEnd() puts the iterator there.
Written by Owen Carmichael, January 2007

*/

#ifndef __itkCollectionOfPatchesHandLabeledIterator_h
#define __itkCollectionOfPatchesHandLabeledIterator_h

// Include files for ITK content:

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkRegionOfInterestImageFilter.h" 
#include "itkCollectionOfPatchesIterator.h"
#include "itkMinimumImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkCastImageFilter.h"

namespace itk
{
  // The class is templated over the image type of the images to be read from 
  // file, as well as the image type for the image regions being returned to the
  // user.  Its parent class, CollectionOfPatchesIterator, is the base class for
  // all iterators that iterate over sets of image regions of interest

  template<typename TInputImage,typename TOutputImage> 
    class ITK_EXPORT CollectionOfPatchesHandLabeledIterator : public CollectionOfPatchesIterator<TInputImage,TOutputImage>
  {
  public:
    // These three are standard typedefs that are part of most ITK classes.
    // They are useful for bookkeeping, memory management, and debugging.

    typedef CollectionOfPatchesHandLabeledIterator<TInputImage,TOutputImage> Self;
    typedef CollectionOfPatchesIterator<TInputImage,TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;

    // convenience typedefs for the data types used by the class:

    typedef typename TInputImage::IndexType  IndexType;
    typedef typename TInputImage::SizeType    SizeType;
    typedef typename TInputImage::RegionType   RegionType;
    typedef ImageFileReader< TInputImage > ImageReaderType;
    typedef RegionOfInterestImageFilter< TInputImage, TOutputImage > PatchExtractorType;
    typedef Function::Minimum<float,float,float> FMinType;
    typedef Function::Maximum<float,float,float> FMaxType;
    typedef Functor::Cast< float, int > F2ICast;
    typedef typename PatchExtractorType::OutputImageType PatchExtractorOutputType;
    typedef typename ImageReaderType::Pointer ImageReaderPointerType;
    typedef typename PatchExtractorType::Pointer PatchExtractorPointerType;

      /** Method for creation through the object factory. */
    itkNewMacro(Self);
    
  /** Run-time type information (and related methods). */
  itkTypeMacro(CollectionOfPatchesHandLabeledIterator, CollectionOfPatchesIterator);

  void SetImagePath(char *); // This sets the path to the directory with the images and the ground_truth file in it, and reads in the image region descriptions from ground_truth
  void SetCushionPctX(float);  // Sets the amount of padding surrounding the specified image features that should be included in the image region.  The amount of padding is specified as a percentage of the size of the bounding box around the features.  Specifically, if a set of face image features is enclosed in a bounding box going from min_x to max_x in the x direction and min_y to max_y in the y direction, the image region extracted from the image will go from min_x-cushion_pct_x*(max_x-min_x) to max_x+cushion_pct_x*(max_x-min_x) in the x direction and from min_y-cushion_pct_y*(max_y-min_y) to max_y+cushion_pct_y*(max_y-min_y) in the y direction.
  void SetCushionPctY(float); // Sets the amount of padding surrounding the specified image features that should be included in the image region.  The amount of padding is specified as a percentage of the size of the bounding box around the features.  Specifically, if a set of face image features is enclosed in a bounding box going from min_x to max_x in the x direction and min_y to max_y in the y direction, the image region extracted from the image will go from min_x-cushion_pct_x*(max_x-min_x) to max_x+cushion_pct_x*(max_x-min_x) in the x direction and from min_y-cushion_pct_y*(max_y-min_y) to max_y+cushion_pct_y*(max_y-min_y) in the y direction.
  
  int GetNumberOfPatches();  // Returns the number of patches listed in the ground_truth file
  int GetCurrentPatchNum();  // Returns the index of the image patch currently pointed to by the iterator.  The first patch in the ground_truth file has index 0.

  float GetCushionPctX();    //  Returns the amount of padding added to an image patch in the X direction, specified as a percentage of patch size in the X direction
  float GetCushionPctY();     //  Returns the amount of padding added to an image patch in the Y direction, specified as a percentage of patch size in the Y direction
  
  virtual SmartPointer<TOutputImage> Get(); // This function reads the current image in from file, extracts the desired image region, and returns a pointer to it.
  virtual void GoToBegin();  // Move the iterator to the first image region in the ground_truth file.
  virtual void GoToEnd();  // Move to the end of the list of regions.
  virtual bool IsAtEnd();  // Tests whether the iterator is at the end of the region list
  virtual bool IsAtBegin();   // Tests whether the iterator is at the beginning of the region list
  virtual void operator++();  // Increment the iterator to point to the next image region in the list
  virtual void operator--();  // Decrement the iterator to point to the previous image region in the list
  virtual void operator+=(int); // Increment the iterator to point to the image region that is some number of regions past of the current one
  virtual void operator-=(int); // Decrement the iterator to point to the image region that is some number of regions before the current one
  virtual void SetPosition( int );  // Set the iterator to point to the k-th image region in the list
  protected:
  CollectionOfPatchesHandLabeledIterator();  // Default constructor: allocates the ImageFileReader and RegionOfInterestFilter
  ~CollectionOfPatchesHandLabeledIterator(); // Destructor:  Deallocates internal storage for the image region descriptions that had been read in from file
  void ReadGroundTruthFile();  // Reads in the image region descriptions from the ground truth file
  void AssertAllocated();      // Checks that internal storage for the image region descriptions has been allocated; called prior to accessing those descriptions
  void AssertValidPosition(int,std::string);  // Checks that the iterator is currently pointing to a valid position, i.e. either pointed at a patch or at the END.
  void AssertPatchPosition(int,std::string);  // Checks that the iterator is currently pointing to a valid image region.  This differs from AssertValidPosition because the END is a valid position, but there are no image regions there.

  // patch_info is the data structure for a single image region description read in from file.  It keeps track of the X and Y locations of the face features and the name of the image file that the face is in.
  struct patch_info {
  public:
    char image_filename[1024];
    float left_eye_x,right_eye_x,nose_x,left_corner_mouth_x,center_mouth_x,right_corner_mouth_x;
    float left_eye_y,right_eye_y,nose_y,left_corner_mouth_y,center_mouth_y,right_corner_mouth_y;
  };
  patch_info *patches;  // The internal storage for the image region descriptions read by ground_truth is an array of patch_info structs.  It is allocated by ReadGroundTruthFile()
  int num_patches;  // The number of image regions specified in ground_truth
  int min_patchnum; // The index of the first patch in patches[].  This is always 0 for now...
  int max_patchnum; // The index that is ONE PAST the index of the final patch in patches[].  This is always set to num_patches for now...
  int current_patchnum; // The index of the current patch in patches[] that the iterator is pointing to.  
  bool patches_is_allocated; // True if we have allocated internal storage for the image region descriptions in patches[]
  char ground_truth_filename[1024]; // Name of the file with the image region descriptions in it-- for now, this is always image_path followed by "ground_truth"
  char image_path[1024];  // The path to the directory containing the images and the ground_truth file containing image region descriptions
  float cushion_pct_x,cushion_pct_y; // What percentage of region width (cushion_pct_x) and height (cushion_pct_y) should be added around the bounding box of the face features as padding when extracting the image patch for a face
  ImageReaderPointerType imageReader; // A pointer to the ITK module that reads images in from file
  PatchExtractorPointerType patchExtractor; // A pointer to the ITK module that extracts a rectangular region of interest from an image
  F2ICast f2icast;  // A function that casts floating point numbers to integers.  It is used to convert the the face feature locations (specified as floats in the ground_truth file) to pixel locations in the images
  FMinType fmin;  // Computes the minimum of 2 floats.  Used for computing the bounding box around a set of image features
  FMaxType fmax;  // Computes the maximum of 2 floats.  Used for computing the bounding box around a set of image features
  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCollectionOfPatchesHandLabeledIterator.txx"
#endif

#endif // #ifndef __itkCollectionOfPatchesHandLabeledIterator_h
