#ifndef _RJBUTILITIES_H
#define _RJBUTILITIES_H

#include "itkLabelImageToStatisticsLabelMapFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include <itkBinaryThresholdImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageToHistogramFilter.h>
#include <itkMaskedImageToHistogramFilter.h>

#if 1
/////////////////////////////////////////////////////////
template <class RawIm, class MaskIm, class RealType>
std::vector<RealType>
computeMaskQuantiles(typename RawIm::Pointer raw, typename MaskIm::Pointer mask, 
		     std::vector<RealType> quantiles)
{
  // robust rescaling based on median and IQR
  typedef typename itk::LabelImageToStatisticsLabelMapFilter <MaskIm, RawIm> LabStatsType;
  
  typename LabStatsType::Pointer statsfilt = LabStatsType::New();
  statsfilt->SetFeatureImage(raw);
  statsfilt->SetInput(mask);
  statsfilt->Update();

  // typename LabStatsType::RealType Mn = statsfilt->GetMinimum(1);
  // typename LabStatsType::RealType Mx = statsfilt->GetMaximum(1);

  //statsfilt->SetHistogramParameters(100, Mn, Mx);
  statsfilt->SetComputeHistogram(true);
  statsfilt->Modified();
  statsfilt->Update();

  // the output image is actually a label map
  typename LabStatsType::OutputImageType * labelMap = statsfilt->GetOutput();
  const typename LabStatsType::LabelObjectType * labelObject = labelMap->GetLabelObject( 1 );
  const typename LabStatsType::LabelObjectType::HistogramType * hist = labelObject->GetHistogram();

  std::vector<RealType> result;
  for (unsigned K=0;K<quantiles.size();K++)
    {
    result.push_back((RealType)hist->Quantile(0, quantiles[K]));
    }
  return(result);
}
#endif
template <class RawIm, class MaskIm>
void
computeMaskMeanVar(typename RawIm::Pointer raw, typename MaskIm::Pointer mask,
		   typename MaskIm::PixelType Lab,
		   double &Mean,
		   double &Var)
{
  // robust rescaling based on median and IQR
  typedef typename itk::LabelStatisticsImageFilter<RawIm, MaskIm> StatsType;
  
  typename StatsType::Pointer statsfilt = StatsType::New();
  statsfilt->SetInput(raw);
  statsfilt->SetLabelInput(mask);
  statsfilt->Update();

  Mean = (double)statsfilt->GetMean(Lab);
  Var = (double)statsfilt->GetVariance(Lab);

}
//////////////////////////////////////////////

template <class TImage>
std::vector<typename TImage::PixelType> computeImQuantile(typename TImage::Pointer in, std::vector<float> quants)
{

  typedef typename itk::Statistics::ImageToHistogramFilter<TImage> HistMakerType;


  typename HistMakerType::Pointer HistMaker = HistMakerType::New();
  HistMaker->SetInput(in);

  HistMaker->SetAutoMinimumMaximum(true);
  typename HistMakerType::HistogramSizeType hsize(in->GetNumberOfComponentsPerPixel());
  hsize.Fill(1024);
  HistMaker->SetHistogramSize(hsize);
  HistMaker->Update();

  std::vector< typename TImage::PixelType > result(quants.size(), 0);

  const typename HistMakerType::HistogramType * hP = HistMaker->GetOutput();

  int size = hP->GetSize(0);

  int first = 0;
  while( first < size && hP->GetFrequency(first, 0) == 0 )
    {
    first++;
    }
  if (first == size)
    {
    std::cerr << "No data in histogram";
    return(result);
    }
  for (unsigned i = 0; i < quants.size(); i++)
    {
    result[i]=hP->Quantile(0, quants[i]);
    }
  
  return(result);

}
//////////////////////////////////////////////

template <class TImage, class MaskImage>
std::vector<typename TImage::PixelType> computeImQuantile(typename TImage::Pointer in, typename MaskImage::Pointer mask, std::vector<float> quants)
{
  //std::cout << "start hist create" << std::endl;
  typedef typename itk::Statistics::MaskedImageToHistogramFilter<TImage, MaskImage> HistMakerType;

  mask->Update();
  typename HistMakerType::Pointer HistMaker = HistMakerType::New();
  HistMaker->SetInput(in);
  HistMaker->SetMaskImage(mask);
  HistMaker->SetAutoMinimumMaximum(true);
  HistMaker->SetMaskValue(1);

  typename HistMakerType::HistogramSizeType hsize(in->GetNumberOfComponentsPerPixel());
  hsize.Fill(256);
  HistMaker->SetHistogramSize(hsize);
  HistMaker->Update();
  // std::cout << "start hist process" << std::endl;
 std::vector< typename TImage::PixelType > result(quants.size(), 0);

  const typename HistMakerType::HistogramType * hP = HistMaker->GetOutput();

  int size = hP->GetSize(0);

  int first = 0;
  while( first < size && hP->GetFrequency(first, 0) == 0 )
    {
    first++;
    }
  if (first == size)
    {
    std::cerr << "No data in histogram";
    return(result);
    }
  for (unsigned i = 0; i < quants.size(); i++)
    {
    result[i]=hP->Quantile(0, quants[i]);
    }
  
   return(result);
}

/////////////////////////////////////////////////////////
template <class RawIm, class MaskIm>
typename MaskIm::Pointer doThresh(typename RawIm::Pointer raw, float threshVal, float scale = 1.0)
{
  // a convenience function - sacrifices streaming
  typedef typename itk::BinaryThresholdImageFilter<RawIm, MaskIm> ThreshType;
  typename ThreshType::Pointer wthresh = ThreshType::New();
  wthresh->SetInput(raw);
  // take into account spm's scaling
  wthresh->SetUpperThreshold((typename RawIm::PixelType)(threshVal * scale));
  wthresh->SetLowerThreshold(0);
  wthresh->SetInsideValue(0);
  wthresh->SetOutsideValue(1);
  
  typename MaskIm::Pointer result = wthresh->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}

/////////////////////////////////////////////////////////
template <class RawIm, class MaskIm>
typename MaskIm::Pointer doThresh2(typename RawIm::Pointer raw, float threshVal, float scale = 1.0)
{
  // a convenience function - sacrifices streaming
  typedef typename itk::BinaryThresholdImageFilter<RawIm, MaskIm> ThreshType;
  typename ThreshType::Pointer wthresh = ThreshType::New();
  wthresh->SetInput(raw);
  // take into account spm's scaling
  wthresh->SetUpperThreshold((typename RawIm::PixelType)(threshVal * scale));
  wthresh->SetLowerThreshold(0);
  wthresh->SetInsideValue(1);
  wthresh->SetOutsideValue(0);
  
  typename MaskIm::Pointer result = wthresh->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}

///////////////////////////////////////////////////
template  <class RawIm>
typename RawIm::Pointer resampleIm(typename RawIm::Pointer input, typename RawIm::SpacingType NewSpacing, int interp=1)
{
  const int dim = RawIm::ImageDimension;
  typedef typename RawIm::PixelType PixelType;

  typedef typename itk::ResampleImageFilter<RawIm, RawIm >  ResampleFilterType;
  typedef typename itk::IdentityTransform< double, dim >  TransformType;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  input->Update();

  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  resampler->SetTransform( transform );
  typedef typename itk::LinearInterpolateImageFunction<RawIm, double >  LInterpolatorType;
  typedef typename itk::NearestNeighborInterpolateImageFunction<RawIm, double >  NNInterpolatorType;

  typename ResampleFilterType::InterpolatorPointerType interpolator;
  switch (interp)
    {
    case 0:
      interpolator = NNInterpolatorType::New();
      break;
    case 1:
      interpolator = LInterpolatorType::New();
      break;
    default:
      std::cout << "Unsupported interpolator" << std::endl;
    }

  resampler->SetInterpolator( interpolator );
  resampler->SetDefaultPixelValue( 0 );

  const typename RawIm::SpacingType& inputSpacing = input->GetSpacing();
  //typename RawIm::SpacingType spacing;
  typename RawIm::SizeType   inputSize = input->GetLargestPossibleRegion().GetSize();
  typename RawIm::SizeType   size;
  typename RawIm::PointType  newOrigin = input->GetOrigin();
  typedef typename RawIm::SizeType::SizeValueType SizeValueType;

  // use the continuous index concept to generate a world coordinate
  // for the origin
  typedef typename itk::ContinuousIndex<double, RawIm::ImageDimension> ContIndType;
  
  ContIndType NewContInd;


  // keep eye on this one for a while. Not sure it is correct
  typename RawIm::PointType oldPoint, newPoint;
  input->TransformIndexToPhysicalPoint(input->GetLargestPossibleRegion().GetIndex(), oldPoint);
  newPoint = oldPoint - (inputSpacing - NewSpacing)/2;
  input->TransformPhysicalPointToContinuousIndex(newPoint, NewContInd);
  
  typename RawIm::IndexType idx = input->GetLargestPossibleRegion().GetIndex();
  for (int i = 0; i < dim; i++)
    {
    float factor = inputSpacing[i]/NewSpacing[i];
    size[i] = static_cast< SizeValueType >(round( inputSize[i] * factor));
   newOrigin[i] -= (idx[i] - NewContInd[i])*inputSpacing[i];
    idx[i] *= factor;
    }
  resampler->SetSize( size );
  resampler->SetOutputSpacing( NewSpacing );
  resampler->SetOutputOrigin( newOrigin);

  // need to be careful setting the index
  resampler->SetOutputStartIndex ( idx );
  resampler->SetOutputDirection(input->GetDirection());
  resampler->SetInput(input);
  typename RawIm::Pointer result = resampler->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}

///////////////////////////////////////////////////
template  <class RawIm>
typename RawIm::Pointer resampleIm(typename RawIm::Pointer input, typename RawIm::Pointer exampleIm, int interp=1)
{
  const int dim = RawIm::ImageDimension;
  typedef typename RawIm::PixelType PixelType;

  typedef typename itk::ResampleImageFilter<RawIm, RawIm >  ResampleFilterType;
  typedef typename itk::IdentityTransform< double, dim >  TransformType;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  input->Update();

  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  resampler->SetTransform( transform );
  typedef typename itk::LinearInterpolateImageFunction<RawIm, double >  LInterpolatorType;
  typedef typename itk::NearestNeighborInterpolateImageFunction<RawIm, double >  NNInterpolatorType;

  typename ResampleFilterType::InterpolatorPointerType interpolator;
  switch (interp)
    {
    case 0:
      interpolator = NNInterpolatorType::New();
      break;
    case 1:
      interpolator = LInterpolatorType::New();
      break;
    default:
      std::cout << "Unsupported interpolator" << std::endl;
    }

  resampler->UseReferenceImageOn();
  resampler->SetReferenceImage(exampleIm);
  resampler->SetInterpolator( interpolator );
  resampler->SetDefaultPixelValue( 0 );

  resampler->SetInput(input);
  typename RawIm::Pointer result = resampler->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}


#endif
