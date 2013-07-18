// A marker based watershed scalper for macaque scans.
// Requires that the scan is in RL PA IS or LR PA IS orientation,
// as per the human mni space. Angle within this space doesn't matter
// much.
#include <algorithm>
#include <itkBinaryThresholdImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkMorphologicalWatershedFromMarkersImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkBinaryShapeOpeningImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkMaskNegatedImageFilter.h>
#include <itkMaximumImageFilter.h>
#include "itkBinaryShapeKeepNObjectsImageFilter.h"
#include "itkBinaryDilateParaImageFilter.h"
#include "itkBinaryErodeParaImageFilter.h"
#include "itkBinaryCloseParaImageFilter.h"
#include "itkBinaryOpenParaImageFilter.h"
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkMaskedRankImageFilter.h>
#include "itkBinaryFillholeImageFilter.h"

// Aidan's trick
#include <itkSmartPointer.h>
namespace itk
{
    template <typename T>
    class Instance : public T::Pointer {
    public:
        Instance() : SmartPointer<T>( T::New() ) {}
    };
}

// predeclarations
template <class TImage>
void writeImDbg(typename TImage::Pointer Im, std::string filename);


#include "tclap/CmdLine.h"
#include "ioutils.h"
#include "rjbutilities.h"
#include "morphutils.h"
#include "biascorrection.h"
#include "cropneck.h"


namespace TCLAP {
template<>
struct ArgTraits< std::vector<float> > {
  typedef StringLike ValueCategory;
};

template<>
void SetString< std::vector<float> >(std::vector<float> &v,
                                      const std::string &s)
{
  std::istringstream iss(s);
  std::vector<float> p;
  float tmp;
  while (iss >> tmp)
    {
    p.push_back(tmp);
    }
  // only replace the default if something happened
  if (p.size())
    {
    v = p;
    }
}
}



typedef class CmdLineType
{
public:
  std::string InputIm, OutputImPrefix, Suffix;
  float ClosingSize, resamplesize;
  std::vector<float> smoothGradSigma;
  int PreOpen;
  bool refine;
  bool debug;
  bool biascorrect;
  int threads;
} CmdLineType;


bool wssdebug = false;
std::string wsssuffix, wssprefix;

// defines the region in which we carry out refinement
const float refine_band_radius = 10.0;
// parameter for the directional gradient computation
float refine_grad_sigma = 0.25;
// size of band around median used
float refine_quantile_range = 0.5;

float macRadLR = 35;
float macRadAP = 45;
float macRadSI = 33;
float macBMarkRad = 5;

// used to find the initial brain marker
// const float mkmarker_gs_close_size = 10.0;
// //const float mkmarker_gs_erode_size = 5.0;
// const float mkmarker_binary_erode_size = 1.0;

// parameters for filtering the first stage rough segmentation
float rough_open_radius = 3.0;
float rough_close_radius = 3.5;


void ParseCmdLine(int argc, char* argv[],
                  CmdLineType &CmdLineObj
                  )
{
  using namespace TCLAP;
  try
    {
    // Define the command line object.
    CmdLine cmd("MBWSS - marker based watershed scalper - macaque version", ' ', "0.9");

    ValueArg<std::string> inArg("i","input","input image",true,"result","string");
    cmd.add( inArg );

    // ValueArg<std::string> markArg("m","marker","marker image",false,"result","string");
    // cmd.add( markArg );

    ValueArg<std::string> outArg("o","output","output image prefix", true,"","string");
    cmd.add( outArg );

    ValueArg<std::string> sufArg("s","suffix","output suffix", false,".nii.gz","string");
    cmd.add( sufArg );

    ValueArg<float> resArg("r","resample","resample", false, 0.0,"float");
    cmd.add( resArg );

    SwitchArg refineArg("f", "refine", "refine segmentation - don't use for infants", false);
    cmd.add( refineArg );

    SwitchArg debugArg("d", "debug", "save debug images", false);
    cmd.add( debugArg );

    ValueArg<float> closeArg("c","closeradius","size of the closing in (mm)", false, 15,"float");
    cmd.add(closeArg);

    std::vector<float> DefSig;
    DefSig.push_back(1.0);
    ValueArg< std::vector<float> > smoothGradArg("", "smoothGradSigma", "smoothing sizes for gradient in refinement step (mm) - include multiple scales in quotes",
                                                 false, DefSig,
                                                 "vector of smoothing sizes");
    cmd.add(smoothGradArg);

    SwitchArg biasArg("", "biascorrect", "Preform rough bias correction", false);
    cmd.add( biasArg );

    ValueArg<int> preopeningArg("", "preopening", "radius of opening before doing anything else - useful if eyes get included", false, 0, "int");
    cmd.add(preopeningArg);

    ValueArg<int> threadArg("", "threads", "How many threads to use. Default is max cores", false, -1, "int");
    cmd.add(threadArg);

   // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputImPrefix = outArg.getValue();
    // CmdLineObj.MarkerIm = markArg.getValue();
    CmdLineObj.ClosingSize = closeArg.getValue();
    CmdLineObj.Suffix = sufArg.getValue();
    CmdLineObj.resamplesize = resArg.getValue();
    CmdLineObj.refine = refineArg.getValue();
    CmdLineObj.smoothGradSigma = smoothGradArg.getValue();
    CmdLineObj.debug = debugArg.getValue();
    CmdLineObj.biascorrect = biasArg.getValue();
    CmdLineObj.PreOpen = preopeningArg.getValue();
    CmdLineObj.threads = threadArg.getValue();
    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}
////////////////////////////////////////////////////////
template <class TImage>
void writeImDbg(typename TImage::Pointer Im, std::string filename)
{
  if (wssdebug)
    {
    writeIm<TImage>(Im, wssprefix + "_" + filename + wsssuffix);
    }
}


template <class MaskImType>
typename MaskImType::Pointer mkMaskBorder(typename MaskImType::Pointer mask, float thickness, bool inside, typename MaskImType::Pointer &interior, bool UseImageSpacing=true)
{
  // create a border ring - mask must be 1/0 value
  if (inside)
    {
    itk::Instance<itk::BinaryErodeParaImageFilter<MaskImType, MaskImType> > ParaErode;
    ParaErode->SetInput(mask);
    ParaErode->SetUseImageSpacing(UseImageSpacing);
    ParaErode->SetRadius(thickness);
    itk::Instance<itk::SubtractImageFilter<MaskImType, MaskImType> > sub;
    sub->SetInput(mask);
    sub->SetInput2(ParaErode->GetOutput());
    typename MaskImType::Pointer result = sub->GetOutput();
    result->Update();
    result->DisconnectPipeline();
    interior = ParaErode->GetOutput();
    interior->Update();
    interior->DisconnectPipeline();
    return(result);
    }
  else
    {
    itk::Instance<itk::BinaryDilateParaImageFilter<MaskImType, MaskImType> > ParaDilate;
    ParaDilate->SetInput(mask);
    ParaDilate->SetUseImageSpacing(UseImageSpacing);
    ParaDilate->SetRadius(thickness);
    itk::Instance<itk::SubtractImageFilter<MaskImType, MaskImType> > sub;
    sub->SetInput(ParaDilate->GetOutput());
    sub->SetInput2(mask);
    typename MaskImType::Pointer result = sub->GetOutput();
    result->Update();
    result->DisconnectPipeline();
    interior = ParaDilate->GetOutput();
    interior->Update();
    interior->DisconnectPipeline();
    return(result);
    }
}

template <class RawImType, class LabImType>
typename LabImType::Pointer bigDarkAreas(typename RawImType::Pointer t1, int y)
{
  // try to find the eyeballs to create new background markers. Other
  // dark areas may be found too.
  // big tophat filter
  typename RawImType::Pointer closed = doClosingMM<RawImType>(t1, 4);
  //typename RawImType::Pointer bth = doBlackTopHatMM<RawImType>(closed, 7);
  //writeImDbg<RawImType>(bth, "eyebth");
  writeImDbg<RawImType>(closed, "eyeclosed");

  // threshold -- need to be careful because we can get interactions
  // with bias correction phase

  std::vector<float> quants;
  quants.push_back(0.01);

  std::vector<typename RawImType::PixelType> qvals = computeImQuantile<RawImType>(closed, quants);
  itk::Instance< itk::BinaryThresholdImageFilter <RawImType, LabImType> > Thresh;

  Thresh->SetInput(closed);
  Thresh->SetUpperThreshold(qvals[0]);
  Thresh->SetInsideValue(0);
  Thresh->SetOutsideValue(1);

  typename LabImType::Pointer zmask = Thresh->GetOutput();
  zmask->Update();
  zmask->DisconnectPipeline();
  typename LabImType::Pointer zmaskinv = Thresh->GetOutput();
  Thresh->SetInsideValue(2);
  Thresh->SetOutsideValue(0);
  zmaskinv->Update();
  zmaskinv->DisconnectPipeline();

  // now get back to the eyes
  quants[0]=0.25;

  qvals = computeImQuantile<RawImType, LabImType>(closed, zmask, quants);
  Thresh->SetUpperThreshold(qvals[0]);
  Thresh->SetInsideValue(2);
  Thresh->SetOutsideValue(0);
  // merge with zmask
  itk::Instance< itk::MaximumImageFilter <LabImType, LabImType, LabImType> > Comb;
  Comb->SetInput(zmaskinv);
  Comb->SetInput2(Thresh->GetOutput());

  writeImDbg<LabImType>(Comb->GetOutput(), "combmask");

  typename LabImType::Pointer tmp = Comb->GetOutput();
  {
  tmp->Update();
  tmp->DisconnectPipeline();
  // delete the  back, then dilate a bit
  typename LabImType::RegionType blank = tmp->GetLargestPossibleRegion();
  typename LabImType::RegionType::SizeType s  = blank.GetSize();
  typename LabImType::SpacingType sp = tmp->GetSpacing();
  //int bot = top - macRadSI/sp[2];
  //if (bot > 0)
    {
    //s[2] = bot;
    s[1] = y + 0.5*macRadAP/sp[1];
    blank.SetSize(s);
    fillRegion<LabImType>(tmp, blank, 0);
    }
  // blank = tmp->GetLargestPossibleRegion();
  // s=blank.GetSize();
  // s[2]=top - macRadSI/sp[2];
  // s[1]=y;
  // blank.SetSize(s);
  // fillRegion<LabImType>(tmp, blank, 0);

  }
  itk::Instance< itk::BinaryShapeOpeningImageFilter<LabImType> > AreaOpening;
  AreaOpening->SetInput(tmp);
  AreaOpening->SetLambda(100);
  AreaOpening->SetAttribute("PhysicalSize");
  AreaOpening->SetForegroundValue(2);


  typename LabImType::Pointer result = doDilateMM<LabImType>(AreaOpening->GetOutput(), 3, 3, 0);
  result->Update();
  result->DisconnectPipeline();

  writeImDbg<LabImType>(result, "thresheyebth");
  return(result);
}

////////////////////////////////////////////////////////
template <class RawImType>
typename RawImType::Pointer scaleSpaceSmooth(typename RawImType::Pointer input, std::vector<float> scales)
{
  // smooth at different scales and take a max
  typename RawImType::Pointer res, tmp;

  itk::Instance< itk::SmoothingRecursiveGaussianImageFilter< RawImType, RawImType > > Smoother;
  itk::Instance< itk::MaximumImageFilter<RawImType, RawImType, RawImType> > Max;

  if (scales.size() == 0)
    {
    std::cerr << "At least one smoothing scale must be specified" << std::endl;
    return(0);
    }
  Smoother->SetInput(input);
  Smoother->SetNormalizeAcrossScale(true);
  Smoother->SetSigma(scales[0]);
  res = Smoother->GetOutput();
  res->Update();
  res->DisconnectPipeline();

  for (unsigned i=1; i<scales.size(); i++)
    {
    Smoother->SetSigma(scales[i]);
    Max->SetInput(res);
    Max->SetInput2(Smoother->GetOutput());
    tmp=Max->GetOutput();
    tmp->Update();
    tmp->DisconnectPipeline();
    res=tmp;
    }
  return(res);
}
////////////////////////////////////////////////////////

template <class RawImType, class LabImType>
typename LabImType::Pointer mkMarker(typename RawImType::Pointer t1)
{
  // rough estimate of centroid in x/y and top of head
  int bottom, top, x, y;
  // this may modify the input
  cropNeck2<RawImType>(t1, 80, 15, top, bottom, x, y);
  writeImDbg<RawImType>(t1, "cropraw");

  typename LabImType::Pointer Eyes = bigDarkAreas<RawImType, LabImType>(t1, y);

  // construct a marker image
  typename LabImType::Pointer marker = LabImType::New();
  marker->SetRegions(t1->GetLargestPossibleRegion());
  marker->Allocate();
  marker->CopyInformation(t1);
  marker->FillBuffer(2);
  // create a brain marker
  typename LabImType::IndexType bmcent;
  typename LabImType::SpacingType bmspacing = t1->GetSpacing();
  typename LabImType::PointType bmcentwrld;
  bmcent[0]=x;
  bmcent[1]=y;
  bmcent[2]=top - 0.5*macRadSI/bmspacing[2];
  t1->TransformIndexToPhysicalPoint(bmcent, bmcentwrld);

  // create the brain "hole"
  // need to leave at least the border in place


  typename LabImType::RegionType blank = marker->GetLargestPossibleRegion();
  typename LabImType::SizeType bsz = blank.GetSize();
  typename LabImType::IndexType bInd = blank.GetIndex();
  typename LabImType::SizeType origsz = blank.GetSize();

  bInd[0]=(int)std::max(1.0, x-macRadLR/bmspacing[0]);
  bInd[1]=(int)std::max(1.0, y-macRadAP/bmspacing[1]);
  bInd[2]=(int)std::max(1, bottom);

  bsz[0]=2*macRadLR/bmspacing[0];
  bsz[1]=2*macRadAP/bmspacing[1];
  bsz[2]=2*macRadSI/bmspacing[2];

  // make sure the high side has a marker
  for (unsigned k=0; k < LabImType::ImageDimension; k++)
    {
    bsz[k] = std::min(bsz[k], origsz[k] - bInd[k] - 1);
    }
  blank.SetSize(bsz);
  blank.SetIndex(bInd);
  fillRegion<LabImType>(marker, blank, 0);
  fillBoxMM<LabImType>(marker, bmcentwrld, 1, macBMarkRad);
  writeImDbg<LabImType>(marker, "bmark");

  itk::Instance< itk::MaximumImageFilter <LabImType, LabImType, LabImType> > Comb;
  Comb->SetInput(marker);
  Comb->SetInput2(Eyes);
  typename LabImType::Pointer result = Comb->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}

/////////////////////////////////////////////////////////

template <class ImType>
typename ImType::Pointer mkControlImage(typename ImType::Pointer im)
{
  typename ImType::PixelType ImMax;

  itk::Instance < itk::StatisticsImageFilter<ImType> > StatsFilt;
  StatsFilt->SetInput(im);
  StatsFilt->Update();
  ImMax = StatsFilt->GetMaximum();


  // Invert intensity
  itk::Instance< itk::ShiftScaleImageFilter<ImType, ImType> > Inverter;
  Inverter->SetInput(im);
  Inverter->SetShift(-1 * ImMax);
  Inverter->SetScale(-1);
  typename ImType::Pointer result = Inverter->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

template <class RawImType, class LabImType>
typename RawImType::Pointer prefilt(typename RawImType::Pointer raw,
                                    typename LabImType::Pointer mask)
{
  itk::Instance< itk::MaskImageFilter< RawImType, LabImType> > ApplyStage1;

  ApplyStage1->SetInput(raw);
  ApplyStage1->SetInput2(mask);

  // now we apply a thresholding to reset WM intensity

  std::vector<float> quantlevs;
  quantlevs.push_back(0.5);

  std::vector<typename RawImType::PixelType> rthresh = computeImQuantile<RawImType, LabImType>(ApplyStage1->GetOutput(), mask, quantlevs);

  itk::Instance<itk::ThresholdImageFilter<RawImType> > resetWM;
  resetWM->SetInput(ApplyStage1->GetOutput());
  resetWM->ThresholdAbove(rthresh[0]);
  resetWM->SetOutsideValue(rthresh[0]);

  typename RawImType::Pointer result = resetWM->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
template <class RawImType, class MaskImType>
typename RawImType::Pointer doOldGrad(typename RawImType::Pointer raw, typename MaskImType::Pointer roughseg, double sigma1)
{
  itk::Instance< itk::GradientMagnitudeRecursiveGaussianImageFilter< RawImType> > GradFilt;
  GradFilt->SetInput(raw);
  GradFilt->SetSigma(sigma1);
  typename RawImType::Pointer result = GradFilt->GetOutput();
  itk::Instance<itk::MaskImageFilter<RawImType, MaskImType> > Masker;
  if (roughseg)
    {

    Masker->SetInput(GradFilt->GetOutput());
    Masker->SetInput2(roughseg);

    result = Masker->GetOutput();
    }
  result->Update();
  result->DisconnectPipeline();
  return(result);
}

template <class LabImType, class RawImType>
typename LabImType::Pointer doRefine(typename LabImType::Pointer mask,
                                     typename RawImType::Pointer raw0,
                                     int label,
                                     std::vector<float> sigma)
{
  // differences to human refinement - signal dropout means that
  // inferior parts of the rough mask are wrong by a long way. This
  // implies that we need to be careful if we use the rough mask to
  // derive the foreground marker

  // Don't need to worry about bright background markers - the skull
  // is much thicker and there don't seem to be problems with dural layers
  typedef typename RawImType::Pointer PRawImType;
  typedef typename LabImType::Pointer PLabImType;

  typedef typename itk::Image<float, RawImType::ImageDimension> FloatImType;
  typedef typename FloatImType::Pointer PFloatImType;

  float radius = 5.0;
  // strategy this time is to find dark areas in an eroded version of
  // the input. We are also going to use the eroded input to produce
  // the control image and then dilate the result.

  // prefiltering - tophat filtering to produce markers later.
  // reset brightness of everything above median to remove internal
  // edges.


  itk::Instance< itk::MaskImageFilter< RawImType, LabImType> > ApplyStage1B;
  ApplyStage1B->SetInput(raw0);
  ApplyStage1B->SetInput2(mask);
  PRawImType rawerode = ApplyStage1B->GetOutput();

  // now we apply a thresholding to reset WM intensity
  std::vector<float> quantlevs;
  quantlevs.push_back(0.5);

  std::vector<typename RawImType::PixelType> rthresh = computeImQuantile<RawImType, LabImType>(raw0, mask, quantlevs);
  itk::Instance<itk::ThresholdImageFilter<RawImType> > resetWM;
  resetWM->SetInput(ApplyStage1B->GetOutput());
  resetWM->ThresholdAbove(rthresh[0]);
  resetWM->SetOutsideValue(rthresh[0]);

  // this is a masked and upper threshold image
  PRawImType raw = resetWM->GetOutput();

  writeImDbg<RawImType>(raw, "csfprefilt");

  PLabImType eroded, erodedthin;
  PLabImType ring = mkMaskBorder<LabImType>(mask, radius, true, eroded);
  PLabImType thinring = mkMaskBorder<LabImType>(mask, radius/3, true, erodedthin);
  // construct new markers
  itk::Instance< itk::BinaryThresholdImageFilter <LabImType, LabImType> > Invert;
  Invert->SetInput(mask);
  Invert->SetInsideValue(0);
  Invert->SetOutsideValue(2);
  Invert->SetUpperThreshold(label);
  Invert->SetLowerThreshold(label);

  // make sure interior marker is only in bright areas (>= median used
  // for thresholding
  itk::Instance<itk::BinaryThresholdImageFilter<RawImType, LabImType> > BThresh;
  BThresh->SetInput(raw);
  BThresh->SetLowerThreshold(rthresh[0]);
  BThresh->SetInsideValue(1);
  BThresh->SetOutsideValue(0);

  itk::Instance<itk::MaskImageFilter< LabImType, LabImType > > SeedMask;
  SeedMask->SetInput(BThresh->GetOutput());
  SeedMask->SetInput2(eroded);

  itk::Instance< itk::MaximumImageFilter <LabImType, LabImType, LabImType> > CombMarker;
  CombMarker->SetInput(Invert->GetOutput());
  CombMarker->SetInput2(SeedMask->GetOutput());
  writeImDbg<LabImType>(CombMarker->GetOutput(), "refcsfmark");

  // dispose this
  erodedthin=0;

  // construct the control image. This is a combination of the
  // original border intensity and a scale space gradient.
  // To compute a gradient without overdue influence from the rough
  // edges we attempt to interpolate.

  // I think this stage is redundant, as masking happens during prefiltering
  itk::Instance <itk::MaskImageFilter<RawImType, LabImType> > MaskRaw;
  MaskRaw->SetInput(raw);
  MaskRaw->SetInput2(mask);

  // use the fine ring trick again
  typename RawImType::Pointer edgeintensity = 0;
  typename RawImType::Pointer modraw = 0;
  typename LabImType::Pointer finering = 0;
  {
  // extract the intensities along the border of the original
  typename LabImType::Pointer erodedfine = 0;
  finering = mkMaskBorder<LabImType>(mask, 2, true, erodedfine, false);
  erodedfine = 0;

  itk::Instance <itk::MaskImageFilter<RawImType, LabImType> > MaskRaw2;
  MaskRaw2->SetInput(raw);
  MaskRaw2->SetInput2(finering);
  edgeintensity = MaskRaw2->GetOutput();
  edgeintensity->Update();
  edgeintensity->DisconnectPipeline();

  }
  {
  // extract the intensities along the border of the original
  // subtract the median from the ring
  // we recreate edgeintensity here
  std::vector<float> quants;
  quants.push_back(0.5);
  std::vector<typename RawImType::PixelType> qt=computeImQuantile<RawImType,LabImType>(raw, finering, quants);

  itk::Instance< itk::SubtractImageFilter <RawImType> > subtractC;
  subtractC->SetInput(raw);
  subtractC->SetConstant(qt[0]);

  itk::Instance <itk::MaskImageFilter<RawImType, LabImType> > MaskRaw2;
  MaskRaw2->SetInput(subtractC->GetOutput());
  //MaskRaw2->SetInput(raw);
  MaskRaw2->SetInput2(finering);
  edgeintensity = MaskRaw2->GetOutput();
  edgeintensity->Update();
  edgeintensity->DisconnectPipeline();
  }

  typename RawImType::Pointer grad = 0;
  {
  typename RawImType::Pointer grad0 = doGradientOuter<RawImType>(raw,1);
  // mask
  itk::Instance <itk::MaskImageFilter<RawImType, LabImType> > MaskGrad;
  MaskGrad->SetInput(grad0);
  MaskGrad->SetInput2(mask);
  typename RawImType::Pointer gradR = scaleSpaceSmooth<RawImType>(MaskGrad->GetOutput(), sigma);
  grad0 = 0;

  // include edge intenstity
  itk::Instance<itk::MaximumImageFilter<RawImType, RawImType, RawImType> > CombGrad;
  CombGrad->SetInput(gradR);
  CombGrad->SetInput2(edgeintensity);
  grad = CombGrad->GetOutput();
  grad->Update();
  grad->DisconnectPipeline();

  writeImDbg<RawImType>(gradR, "refcsfgrad0");
  writeImDbg<RawImType>(grad, "refcsfgrad");
  }
  itk::Instance<itk::MorphologicalWatershedFromMarkersImageFilter<RawImType, LabImType> > WatershedFiltA;
  WatershedFiltA->SetInput(grad);
  WatershedFiltA->SetMarkerImage(CombMarker->GetOutput());
  WatershedFiltA->SetMarkWatershedLine(true);
  writeImDbg<LabImType>(CombMarker->GetOutput(), "refcsfmark");

// select the brain label and watershed line
  itk::Instance<itk::BinaryThresholdImageFilter<LabImType, LabImType> > Select;
  Select->SetInput(WatershedFiltA->GetOutput());
  Select->SetUpperThreshold(label);
  Select->SetLowerThreshold(0);
  Select->SetInsideValue(label);
  Select->SetOutsideValue(0);
  writeImDbg<LabImType>(WatershedFiltA->GetOutput(), "refcsfwshed");
  typename LabImType::Pointer result = Select->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
/////////////////////////////////////////////////////////

template <class PixType, class LabPixType, int dim>
void doWatershed(const CmdLineType &CmdLineObj)
{
  // the basic algorithm is simple.
  // 1) Rough alignment of brain with the standard template (done
  // outside)
  // 2) Transform a marker image to the subject space
  //    Make some adjustments to the markers based on thresholding,
  // because the registration is never precise enough.

  // 3) Apply marker based watershed, using inverse of original image
  // as the control
  // 4) Apply a large closing to smooth the mask
  // 5) (Optional) Use the result as an external bound and do a second watershed

  typedef typename itk::Image<PixType, dim> RawImType;
  typedef typename itk::Image<LabPixType, dim> LabImType;
  typedef typename RawImType::Pointer RawImPointer;
  typedef typename LabImType::Pointer LabImPointer;

  RawImPointer raworig = readIm<RawImType>(CmdLineObj.InputIm);
  if (CmdLineObj.biascorrect)
    {
    std::cout << "Applying rough bias correction" << std::endl;
    //raworig = crappyBiasCorrection<RawImType>(raworig, 10.0, 0.05);
    raworig = crappyBiasCorrectionB<RawImType>(raworig, 10.0, 0.05);
    writeImDbg<RawImType>(raworig, "biascor");
    }
  LabImPointer newmarker = mkMarker<RawImType, LabImType>(raworig);
  if (CmdLineObj.biascorrect)
    {
    // bias correction sets some parts to 0 (maybe). Include these as
    // a part of the background marker
    itk::Instance<itk::BinaryThresholdImageFilter<RawImType, LabImType> > SelectBiasZero;
    SelectBiasZero->SetInput(raworig);
    SelectBiasZero->SetUpperThreshold(0);
    SelectBiasZero->SetLowerThreshold(0);
    SelectBiasZero->SetInsideValue(2);
    SelectBiasZero->SetOutsideValue(0);
    itk::Instance< itk::MaximumImageFilter <LabImType, LabImType, LabImType> > CombMarker;

    CombMarker->SetInput(newmarker);
    CombMarker->SetInput2(SelectBiasZero->GetOutput());

    LabImPointer tmpmarker = CombMarker->GetOutput();
    tmpmarker->Update();
    tmpmarker->DisconnectPipeline();
    newmarker=tmpmarker;
    }

  RawImPointer raw = raworig;
  LabImPointer newmarker2 = newmarker;

  writeIm<LabImType>(newmarker, CmdLineObj.OutputImPrefix + "_marker" + CmdLineObj.Suffix);
  LabImPointer segPhase1;
  {
  if (CmdLineObj.PreOpen > 0)
    {
    std::cout << "Applying prefilter opening" << std::endl;
    raw = doOpening<RawImType>(raw, CmdLineObj.PreOpen);
    writeImDbg<RawImType>(raw, "preopened");
    }
  RawImPointer control = mkControlImage<RawImType>(raw);
  writeIm<RawImType>(control, CmdLineObj.OutputImPrefix + "_control" + CmdLineObj.Suffix);
  itk::Instance<itk::MorphologicalWatershedFromMarkersImageFilter<RawImType, LabImType> > WatershedFiltA;
  WatershedFiltA->SetInput(control);
  WatershedFiltA->SetMarkerImage(newmarker2);
  WatershedFiltA->SetMarkWatershedLine(false);
//  WatershedFiltA->SetMarkWatershedLine(true);
  segPhase1 = WatershedFiltA->GetOutput();
  segPhase1->Update();
  segPhase1->DisconnectPipeline();
  // std::cout << "Finished WS" << std::endl;
  // std::cout << segPhase1;
  // reverse the resampling

  }

  // select the brain label
  itk::Instance<itk::BinaryThresholdImageFilter<LabImType, LabImType> > Select;
  Select->SetInput(segPhase1);
  Select->SetUpperThreshold(1);
  Select->SetLowerThreshold(1);
  Select->SetInsideValue(1);
  Select->SetOutsideValue(0);

  // not really any point saving this as it is rarely much use, but it
  // keeps the structure the same as the human version
  writeIm<LabImType>(Select->GetOutput(), CmdLineObj.OutputImPrefix + "_brainmaskrough" + CmdLineObj.Suffix);
  //writeIm<LabImType>(segPhase1, CmdLineObj.OutputImPrefix + "_seg" + CmdLineObj.Suffix);
  // itk::Instance< itk::MaskImageFilter <RawImType, LabImType, RawImType> > Masker;
  // Masker->SetInput(raworig);
  // Masker->SetInput2(ParaClose->GetOutput());

  // writeIm<RawImType>(Masker->GetOutput(), CmdLineObj.OutputImPrefix + "_brain" + CmdLineObj.Suffix);

  typename LabImType::Pointer refined = 0;
  if (CmdLineObj.refine)
    {
    // reload raw, now that the marker has been dealt with
    //raw = readIm<RawImType>(CmdLineObj.InputIm);
    refined = doRefine<LabImType, RawImType>(Select->GetOutput(), raw, 1, CmdLineObj.smoothGradSigma);

    // smooth the refined version
    itk::Instance<itk::BinaryOpenParaImageFilter<LabImType, LabImType> > ParaOpenRef;
    itk::Instance<itk::BinaryCloseParaImageFilter<LabImType, LabImType> > ParaCloseRef;

    ParaOpenRef->SetInput(refined);
    ParaOpenRef->SetUseImageSpacing(true);
    ParaOpenRef->SetRadius(rough_open_radius/2);

    ParaCloseRef->SetInput(ParaOpenRef->GetOutput());
    ParaCloseRef->SetUseImageSpacing(true);
    ParaCloseRef->SetRadius(rough_close_radius/2);
    writeIm<LabImType>(ParaCloseRef->GetOutput(), CmdLineObj.OutputImPrefix + "_refinedmask1" + CmdLineObj.Suffix);
    writeIm<LabImType>(refined, CmdLineObj.OutputImPrefix + "_refinedmask0" + CmdLineObj.Suffix);
    refined = ParaCloseRef->GetOutput();
    refined->Update();
    refined->DisconnectPipeline();

    }

}
////////////////////////////////////////////////////////
int main(int argc, char * argv[])
{

  int dim1 = 0;
  CmdLineType CmdLineObj;
  ParseCmdLine(argc, argv, CmdLineObj);

  if (CmdLineObj.threads > 0)
    {
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(CmdLineObj.threads);
    }

  wssdebug = CmdLineObj.debug;
  wssprefix = CmdLineObj.OutputImPrefix;
  wsssuffix = CmdLineObj.Suffix;

  itk::ImageIOBase::IOComponentType ComponentType;
  const itk::ImageIOBase::IOComponentType MarkerComponentType = itk::ImageIOBase::UCHAR;

  if (!readImageInfo(CmdLineObj.InputIm, &ComponentType, &dim1))
    {
    std::cerr << "Failed to open " << CmdLineObj.InputIm << std::endl;
    return(EXIT_FAILURE);
    }

// this will be a big, ugly switch statement to handle all the cases
// we want char, short and int markers, char, short, unsigned short,
// and float raw images
  switch (dim1)
    {
    case 3:
    {
#define WSDIM 3
#include "nasty_switch.h"
#undef WSDIM
    }
    break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      break;
    }

  return EXIT_SUCCESS;
}
