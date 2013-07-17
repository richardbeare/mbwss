// the simplest possible bias correction
#include <itkBinaryThresholdImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkPowImageFilter.h>
#include <itkFastApproximateRankImageFilter.h>
#include <itkRankImageFilter.h>
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkDivideImageFilter.h"
#include "itkLogicOpsFunctors.h"
#include "itkBoxMeanImageFilter.h"

#include "rjbutilities.h"

template <class ImType>
typename ImType::Pointer crappyBiasCorrection(typename ImType::Pointer raw, float sigma = 10.0, float thresh=0.025)
{
  typedef typename itk::Image<float, ImType::ImageDimension> FloatImType;
  typedef typename itk::Image<unsigned char, ImType::ImageDimension> MaskImType;
  typedef typename itk::SmoothingRecursiveGaussianImageFilter< ImType, ImType > SmootherType;
  typename SmootherType::Pointer Smoother = SmootherType::New();

  typename SmootherType::SigmaArrayType sigarray;
  sigarray.Fill(sigma);
  sigarray[2] = sigma/4;

  Smoother->SetInput(raw);
  Smoother->SetSigmaArray(sigarray);

  // Smoother->SetSigma(sigma);
  Smoother->SetNormalizeAcrossScale(true);


  // Figure out max for scaling
  itk::Instance<itk::StatisticsImageFilter<ImType> > Stats;
  Stats->SetInput(Smoother->GetOutput());
  Stats->Update();

  writeImDbg<ImType>(Smoother->GetOutput(), "biasfield");

  itk::Instance<itk::ShiftScaleImageFilter <ImType, FloatImType> > Scale;
  Scale->SetInput(Smoother->GetOutput());
  Scale->SetScale( 1.0 / ((float)Stats->GetMaximum()));

  // Gamma correction
  itk::Instance< itk::PowImageFilter<FloatImType> > Pow;
  Pow->SetInput(Scale->GetOutput());
  Pow->SetConstant2(0.6);
  writeImDbg<FloatImType>(Pow->GetOutput(), "scale");
  //writeImDbg<FloatImType>(Scale->GetOutput(), "scale");

  itk::Instance< itk::DivideImageFilter<ImType, FloatImType, ImType> > Div;

  Div->SetInput(raw);
  Div->SetInput2(Pow->GetOutput());

  writeImDbg<ImType>(Div->GetOutput(), "div");
//  itk::Instance< itk::CastImageFilter <FloatImType, ImType> > Caster;
//  Caster->SetInput(Div->GetOutput());

  // mask out the noise
  itk::Instance< itk::BinaryThresholdImageFilter <FloatImType, MaskImType> > Thresh;
  Thresh->SetInput(Scale->GetOutput());
  Thresh->SetUpperThreshold(thresh);
  Thresh->SetInsideValue(0);
  Thresh->SetOutsideValue(1);

  itk::Instance<itk::MaskImageFilter<ImType, MaskImType> > Masker;
  Masker->SetInput(Div->GetOutput());
  Masker->SetInput2(Thresh->GetOutput());



  typename ImType::Pointer res = Masker->GetOutput();
//  typename ImType::Pointer res = Div->GetOutput();
  res->Update();
  res->DisconnectPipeline();
  //writeImDbg<ImType>(res, "res");
  return(res);

}

template <class ImType>
typename ImType::PixelType robustBGEst(typename ImType::Pointer im)
{
  typedef typename itk::Image<unsigned char, ImType::ImageDimension> MaskImType;
  // robust bg threshold - meant to account for blanked regions
  itk::Instance<itk::StatisticsImageFilter<ImType> > Stats;
  Stats->SetInput(im);
  Stats->Update();

  typename ImType::PixelType minIm = Stats->GetMinimum();

  // create a mask of > min
  itk::Instance< itk::BinaryThresholdImageFilter <ImType, MaskImType> > Thresh;
  Thresh->SetInput(im);
  Thresh->SetUpperThreshold(minIm);
  Thresh->SetInsideValue(0);
  Thresh->SetOutsideValue(1);

  // now compute the 1% threshold
  std::vector<float> quants;
  quants.push_back(0.01);
  std::vector<typename ImType::PixelType> qvals = computeImQuantile<ImType, MaskImType>(im, Thresh->GetOutput(), quants);
  return(qvals[0]);
}

template <class ImType>
typename ImType::Pointer crappyBiasCorrectionB(typename ImType::Pointer raw, float radius = 20.0, float thresh=0.05)
{
  // need to do something about edges, otherwise the brightness ramps
  // up too much at the edge because the background is underestimated.
  // Basically need to mask to define our ROI, but that is a chicken &
  // egg problem. The Haselgrove paper finds the top edge explicitly,
  // so we'll approximate that as follows. First do some basic
  // rejection of really faint stuff: Find the min, create a mask of
  // > min, and find the 1% quant. Ignore anything <= that. Create the
  // smoothed version, create a mask where input >= 0.5*local smoothed
  // value and > bg thresh. Compute the masked smoothed version from
  // this mask - this is simple if we use mean filters, as we can
  // compute the reweighting by using a mean filter on the mask.
  typedef typename itk::Image<float, ImType::ImageDimension> FloatImType;
  typedef typename itk::Image<unsigned char, ImType::ImageDimension> MaskImType;

  typedef typename itk::BoxMeanImageFilter<ImType, FloatImType> SmoothRawType;
  typedef typename itk::BoxMeanImageFilter<MaskImType, FloatImType> SmoothMaskType;

  // first compute a robust bg threshold
  typename ImType::PixelType bgthresh = robustBGEst<ImType>(raw);
  typename SmoothRawType::Pointer smooth1 = SmoothRawType::New();
  typename SmoothRawType::RadiusType rad;
  typename ImType::SpacingType spacing = raw->GetSpacing();
  for (unsigned r=0;r<ImType::ImageDimension;r++)
    {
    rad[r]=radius/spacing[r];
    }
  //rad[2] /= 4.0;
  smooth1->SetInput(raw);
  smooth1->SetRadius(rad);

  // now create a mask of locally bright parts for bias correction
  itk::Instance<itk::DivideImageFilter<FloatImType, FloatImType, FloatImType> > HalfMean;
  HalfMean->SetInput(smooth1->GetOutput());
  HalfMean->SetConstant(2);

  typedef itk::BinaryFunctorImageFilter<ImType, FloatImType, MaskImType,
    itk::Functor::GreaterEqual<typename ImType::PixelType,
                               typename FloatImType::PixelType,
                               typename MaskImType::PixelType> > GEType;

  typename GEType::Pointer GE = GEType::New();
  GE->SetInput(raw);
  GE->SetInput2(HalfMean->GetOutput());

  // combine with the rough bg mask
  itk::Instance< itk::BinaryThresholdImageFilter <ImType, MaskImType> > BGThresh;
  BGThresh->SetInput(raw);
  BGThresh->SetUpperThreshold(bgthresh);
  BGThresh->SetInsideValue(0);
  BGThresh->SetOutsideValue(1);

  itk::Instance< itk::MaskImageFilter <MaskImType, MaskImType, MaskImType> > Combine;
  Combine->SetInput(GE->GetOutput());
  Combine->SetInput2(BGThresh->GetOutput());
  typename MaskImType::Pointer biasmask = Combine->GetOutput();
  biasmask->Update();
  biasmask->DisconnectPipeline();
  // Now use this to mask the input, and recomput the smoothing
  itk::Instance< itk::MaskImageFilter <ImType, MaskImType, ImType> > MaskRaw;
  MaskRaw->SetInput(raw);
  MaskRaw->SetInput2(biasmask);

  smooth1->SetInput(MaskRaw->GetOutput());
  typename SmoothMaskType::Pointer masksmooth = SmoothMaskType::New();
  masksmooth->SetInput(biasmask);
  masksmooth->SetRadius(rad);

  // reweight the output of raw smoother to only include mask voxels
  itk::Instance<itk::DivideImageFilter<FloatImType, FloatImType, FloatImType> > Reweight;
  Reweight->SetInput(smooth1->GetOutput());
  Reweight->SetInput2(masksmooth->GetOutput());

  // Scale the bias field before dividing through
  // Figure out max for scaling
  itk::Instance<itk::StatisticsImageFilter<FloatImType> > FieldStats;
  FieldStats->SetInput(Reweight->GetOutput());
  FieldStats->Update();

  typename FloatImType::PixelType Fmax = FieldStats->GetMaximum();

  itk::Instance<itk::DivideImageFilter<FloatImType, FloatImType, FloatImType> > ScaleField;
  ScaleField->SetInput(Reweight->GetOutput());
  ScaleField->SetConstant(Fmax);

  // Gamma correction - this makes things more stable when
  // inhomogeneity is very strong.
  itk::Instance< itk::PowImageFilter<FloatImType> > ScaleFieldGamma;
  ScaleFieldGamma->SetInput(ScaleField->GetOutput());
  ScaleFieldGamma->SetConstant2(0.6);

  itk::Instance<itk::DivideImageFilter<ImType, FloatImType, ImType> > Correct;

  Correct->SetInput(raw);
  Correct->SetInput2(ScaleFieldGamma->GetOutput());

  // threshold the scale field and apply as a mask
  itk::Instance< itk::BinaryThresholdImageFilter <FloatImType, MaskImType> > FieldThresh;
  FieldThresh->SetInput(ScaleField->GetOutput());
  FieldThresh->SetUpperThreshold(thresh);
  FieldThresh->SetInsideValue(0);
  FieldThresh->SetOutsideValue(1);

  itk::Instance< itk::MaskImageFilter <ImType, MaskImType, ImType> > MaskFinal;
  MaskFinal->SetInput(Correct->GetOutput());
  MaskFinal->SetInput2(FieldThresh->GetOutput());

  writeImDbg<MaskImType>(Combine->GetOutput(), "ge");
  writeImDbg<FloatImType>(ScaleField->GetOutput(), "biasfield");
  typename ImType::Pointer res =MaskFinal->GetOutput();
  res->Update();
  res->DisconnectPipeline();
  //writeImDbg<ImType>(res, "res");
  return(res);

}

// use a rank filter instead of a smoother
template <class ImType>
typename ImType::Pointer crappyBiasCorrection2(typename ImType::Pointer raw, float radius = 10.0, float thresh=0.025)
{
  typedef typename itk::Image<float, ImType::ImageDimension> FloatImType;
  typedef typename itk::Image<unsigned char, ImType::ImageDimension> MaskImType;
#if 1
  typedef typename itk::FastApproximateRankImageFilter< ImType, ImType > SmootherType;
  typename SmootherType::Pointer Smoother = SmootherType::New();

  typename SmootherType::RadiusType rad;
  rad.Fill(radius);
  rad[2] /= 4;

  Smoother->SetInput(raw);
  Smoother->SetRadius(rad);
  Smoother->SetRank(0.5);
#else
  typedef typename itk::RankImageFilter< ImType, ImType > SmootherType;
  typename SmootherType::Pointer Smoother = SmootherType::New();

  typename SmootherType::RadiusType rad;
  rad.Fill(radius);
  rad[2] /= 4;

  Smoother->SetInput(raw);
  Smoother->SetRadius(rad);
  Smoother->SetRank(0.5);

#endif
  // Smoother->SetSigma(sigma);
  //Smoother->SetNormalizeAcrossScale(true);


  // Figure out max for scaling
  itk::Instance<itk::StatisticsImageFilter<ImType> > Stats;
  Stats->SetInput(Smoother->GetOutput());
  Stats->Update();

  writeImDbg<ImType>(Smoother->GetOutput(), "biasfield");

  itk::Instance<itk::ShiftScaleImageFilter <ImType, FloatImType> > Scale;
  Scale->SetInput(Smoother->GetOutput());
  Scale->SetScale( 1.0 / ((float)Stats->GetMaximum()));

  // Gamma correction
  itk::Instance< itk::PowImageFilter<FloatImType> > Pow;
  Pow->SetInput(Scale->GetOutput());
  Pow->SetConstant2(0.6);
  writeImDbg<FloatImType>(Pow->GetOutput(), "scale");
  //writeImDbg<FloatImType>(Scale->GetOutput(), "scale");

  itk::Instance< itk::DivideImageFilter<ImType, FloatImType, ImType> > Div;

  Div->SetInput(raw);
  Div->SetInput2(Pow->GetOutput());

  writeImDbg<ImType>(Div->GetOutput(), "div");
//  itk::Instance< itk::CastImageFilter <FloatImType, ImType> > Caster;
//  Caster->SetInput(Div->GetOutput());

  // mask out the noise
  itk::Instance< itk::BinaryThresholdImageFilter <FloatImType, MaskImType> > Thresh;
  Thresh->SetInput(Scale->GetOutput());
  Thresh->SetUpperThreshold(thresh);
  Thresh->SetInsideValue(0);
  Thresh->SetOutsideValue(1);

  itk::Instance<itk::MaskImageFilter<ImType, MaskImType> > Masker;
  Masker->SetInput(Div->GetOutput());
  Masker->SetInput2(Thresh->GetOutput());



  typename ImType::Pointer res = Masker->GetOutput();
//  typename ImType::Pointer res = Div->GetOutput();
  res->Update();
  res->DisconnectPipeline();
  //writeImDbg<ImType>(res, "res");
  return(res);

}
