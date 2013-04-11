/*
This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.
*/
#ifndef __itkCollectionOfPatchesIterator_h
#define __itkCollectionOfPatchesIterator_h
#include "itkLightObject.h"
#include "itkImage.h"

namespace itk
{
  template<typename TInputImage,typename TOutputImage> 
    class ITK_EXPORT CollectionOfPatchesIterator : public LightObject
  {
  public:
    typedef CollectionOfPatchesIterator Self;
    typedef typename TInputImage::IndexType  IndexType;
    typedef typename TInputImage::SizeType    SizeType;
    typedef typename TInputImage::RegionType   RegionType;

    //    CollectionOfPatchesIterator(const Self& it) { }
    Self &operator=(const Self& it) { };
    virtual SmartPointer<TOutputImage> Get()=0;
    virtual void GoToBegin()=0;
    virtual void GoToEnd()=0;
    virtual bool IsAtEnd()=0;
    virtual bool IsAtBegin()=0; 
    virtual void operator++()=0;
    virtual void operator--()=0;
    virtual void operator+=(int)=0;
    virtual void operator-=(int)=0;
    virtual void SetPosition( int )=0;
    virtual int GetNumberOfPatches()=0;
    /** Method for creation through the object factory. */
    //    itkNewMacro(Self);
    
    /** Run-time type information (and related methods). */
    itkTypeMacro(CollectionOfPatchesIterator, LightObject);
  protected:
    CollectionOfPatchesIterator()  {  }
    ~CollectionOfPatchesIterator() {};
  };
}

#endif // #ifndef __itkCollectionOfPatchesIterator_h
