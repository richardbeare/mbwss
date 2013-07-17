#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"

template <class RawImType>
void cropNeck(typename RawImType::Pointer raw, float distance, int &outtop, int &outbottom)
{
  // simple neck blanking - find the top of head, go distance down and
  // then start blanking
  typedef typename itk::Image<unsigned char, RawImType::ImageDimension> MaskImType;

  itk::Instance<itk::OtsuThresholdImageFilter <RawImType, MaskImType> > Thresh;
  itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter<MaskImType> > discarder;
  Thresh->SetInput(raw);
  Thresh->SetInsideValue(0);
  Thresh->SetOutsideValue(1);
  discarder->SetNumberOfObjects(1);
  discarder->SetInput(Thresh->GetOutput());
  discarder->SetForegroundValue(1);

  typename MaskImType::Pointer mask = discarder->GetOutput();
  mask->Update();
  mask->DisconnectPipeline();
    
  // find the top
  typedef typename itk::BinaryImageToShapeLabelMapFilter<MaskImType> LabellerType;
  typename LabellerType::Pointer labeller = LabellerType::New();

  labeller->SetInput(mask);
  labeller->SetFullyConnected(true);
  labeller->SetInputForegroundValue(1);

  typedef typename LabellerType::OutputImageType LabelMapType;
  typedef typename LabellerType::OutputImagePointer LabelMapPointerType;
  typedef typename LabelMapType::LabelObjectType LabelObjectType;

  writeImDbg<MaskImType>(mask, "cropmask");

  LabelMapPointerType labmap = labeller->GetOutput();
  labmap->Update();
  labmap->DisconnectPipeline();

  // can only be one object
  LabelObjectType * labelObject = labmap->GetLabelObject(1);
  typename MaskImType::RegionType bb = labelObject->GetBoundingBox();


  typename RawImType::SpacingType sp = raw->GetSpacing();
  int top = bb.GetIndex()[2] + bb.GetSize()[2] - 1;
  int bottom = top - distance/sp[2];

  // std::cout << "bb= " << bb << std::endl;
  // std::cout << "bottom= " << bottom << std::endl;

  if (bottom > 0)
    {
    // blank all slices below
    typename RawImType::RegionType blank = raw->GetLargestPossibleRegion();
    typename RawImType::RegionType::SizeType s = blank.GetSize();
    s[2] = bottom;
    blank.SetSize(s);
    fillRegion<RawImType>(raw, blank, 0);
    fillRegion<MaskImType>(mask, blank, 0);
    }
  else 
    {
    bottom = 0;
    }

  // trial clipping extreme values
  std::vector<float> q;
  q.push_back(0.75);
  std::vector<typename RawImType::PixelType> v = computeImQuantile<RawImType, MaskImType>(raw, mask, q);

  itk::Instance<itk::ThresholdImageFilter<RawImType> > resetBright;
  resetBright->SetInput(raw);
  resetBright->ThresholdAbove(v[0]);
  resetBright->SetOutsideValue(v[0]);

  typename RawImType::Pointer r = resetBright->GetOutput();
  r->Update();
  r->DisconnectPipeline();
  raw=r;
  
  writeImDbg<RawImType>(raw, "cropraw");

  outbottom=bottom;
  outtop=top;
}

// alternative version that estimates x/y centroid
// distance is the distance from top to where we start blanking neck
// slice is the thickness of the region at the top we will use to
// estimate x/y centroid
template <class RawImType>
void cropNeck2(typename RawImType::Pointer raw, float distance, float slice, int &outtop, int &outbottom, int &x, int &y)
{
  // simple neck blanking - find the top of head, go distance down and
  // then start blanking
  typedef typename itk::Image<unsigned char, RawImType::ImageDimension> MaskImType;

  itk::Instance<itk::OtsuThresholdImageFilter <RawImType, MaskImType> > Thresh;
  itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter<MaskImType> > discarder;
  Thresh->SetInput(raw);
  Thresh->SetInsideValue(0);
  Thresh->SetOutsideValue(1);
  discarder->SetNumberOfObjects(1);
  discarder->SetInput(Thresh->GetOutput());
  discarder->SetForegroundValue(1);

  typename MaskImType::Pointer mask = discarder->GetOutput();
  mask->Update();
  mask->DisconnectPipeline();
    
  // find the top
  typedef typename itk::BinaryImageToShapeLabelMapFilter<MaskImType> LabellerType;
  typename LabellerType::Pointer labeller = LabellerType::New();

  labeller->SetInput(mask);
  labeller->SetFullyConnected(true);
  labeller->SetInputForegroundValue(1);

  typedef typename LabellerType::OutputImageType LabelMapType;
  typedef typename LabellerType::OutputImagePointer LabelMapPointerType;
  typedef typename LabelMapType::LabelObjectType LabelObjectType;

  writeImDbg<MaskImType>(mask, "cropmask");

  LabelMapPointerType labmap = labeller->GetOutput();
  labmap->Update();
  labmap->DisconnectPipeline();

  // can only be one object
  LabelObjectType * labelObject = labmap->GetLabelObject(1);
  typename MaskImType::RegionType bb = labelObject->GetBoundingBox();


  typename RawImType::SpacingType sp = raw->GetSpacing();
  int top = bb.GetIndex()[2] + bb.GetSize()[2] - 1;
  int bottom = top - distance/sp[2];

  // std::cout << "bb= " << bb << std::endl;
  // std::cout << "bottom= " << bottom << std::endl;

  if (bottom > 0)
    {
    // blank all slices below
    typename RawImType::RegionType blank = raw->GetLargestPossibleRegion();
    typename RawImType::RegionType::SizeType s = blank.GetSize();
    s[2] = bottom;
    blank.SetSize(s);
    fillRegion<RawImType>(raw, blank, 0);
    fillRegion<MaskImType>(mask, blank, 0);
    }
  else 
    {
    bottom = 0;
    }

  outbottom=bottom;
  outtop=top;

  {
  // estimate the centroid using the top 15mm
  typename MaskImType::RegionType blank = mask->GetLargestPossibleRegion();
  typename MaskImType::RegionType::SizeType s = blank.GetSize();
  s[2] = top - slice/sp[2];
  blank.SetSize(s);

  fillRegion<MaskImType>(mask, blank, 0);
  writeImDbg<MaskImType>(mask, "topmask");
  labeller->SetInput(mask);
  labeller->SetFullyConnected(true);
  labeller->SetInputForegroundValue(1);
  LabelMapPointerType labmap2 = labeller->GetOutput();
  labmap2->Update();
  labmap2->DisconnectPipeline();
  LabelObjectType * labelObject2 = labmap2->GetLabelObject(1);
  typename MaskImType::PointType cent = labelObject2->GetCentroid();
  typename MaskImType::IndexType ind;
  mask->TransformPhysicalPointToIndex(cent, ind);
  x=ind[0];
  y=ind[1];
  }
}

