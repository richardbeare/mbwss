// a marker based watershed scalper for human MRI.
#include <itkNumericTraits.h>
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
#include <itkChangeLabelImageFilter.h>
#include "itkBinaryShapeKeepNObjectsImageFilter.h"
#include "itkBinaryDilateParaImageFilter.h"
#include "itkBinaryErodeParaImageFilter.h"
#include "itkBinaryCloseParaImageFilter.h"
#include "itkBinaryOpenParaImageFilter.h"
#include <itkMaskedImageToHistogramFilter.h>
#include <itkMaskedRankImageFilter.h>
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include "itkBinaryFillholeImageFilter.h"
#include "itkBoxMeanImageFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"

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

template <class TImage>
void writeImDbg(const typename TImage::Pointer Im, std::string filename);

#include "tclap/CmdLine.h"
#include "ioutils.h"
#include "rjbutilities.h"
#include "morphutils.h"

typedef class CmdLineType
{
public:
  std::string InputIm, MaskIm, OutputImPrefix, Suffix;
  float ClosingSize;
  std::vector<float> smoothGradSigma;
  int PreOpen;
  bool refine;
  bool biascorrect;
  bool openSkull;
  bool InternalBGMarkers;
  bool debug;
  int threads;
} CmdLineType;

bool wssdebug = false;
std::string wsssuffix, wssprefix;

// configuration parameters - move these to command line or equivalent
// later

// used to find the initial brain marker
const float mkmarker_gs_erode_size = 10.0;
//const float mkmarker_gs_erode_size = 5.0;
const float mkmarker_binary_erode_size = 5.0;

// background marker parameters
const float mkbgmarker_open_radius = 30.0;
//const float mkbgmarker_erode_radius = 10.0;

// used in adjustMarker to filter the background marker (which is
// initially the inverse of the foreground marker. This is the first
// stage filtering to get the outside of the skull.
const float conf_adjust_radius = 5.0;
// second stage filtering in adjustMarker where we try to include skull.
const float conf_adjust_radius2 = 10.0;

// defines the region in which we carry out refinement
const float refine_band_radius = 20.0;
// parameter for the directional gradient computation
//float refine_grad_sigma = 0.25;
// size of band around median used
//float refine_quantile_range = 0.5;
// parameters for filtering the first stage rough segmentation
float rough_open_radius = 5.0;
float rough_close_radius = 6.5;

// find the top and then blank anything below this distance down.
float head_crop_size = 180;

// Distance from our top of head position to the marker centre
float TopToMarkCent  = 50;
float markerRadius = 20.0;

// how much do we erode around the fg marker before we start looking
// for BG markers
float mkbgmarker_erode_radius = 10.0;
float bias_correct_radius = 30.0;
float bias_gamma = 0.8;

// how much to we erode the raw image before refining.
float refine_erode_size = 1.0;

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


void ParseCmdLine(int argc, char* argv[],
                  CmdLineType &CmdLineObj
                  )
{
  using namespace TCLAP;
  try
    {
    // Define the command line object.
    CmdLine cmd("MBWSS - marker based watershed scalper", ' ', "0.9");

    ValueArg<std::string> inArg("i","input","input image",true,"result","string");
    cmd.add( inArg );

    ValueArg<std::string> markArg("m","mask","mask image to apply to markers - usually to remove background markers from the brain stem",false,"","string");
    cmd.add( markArg );

    ValueArg<std::string> outArg("o","output","output image prefix", true,"","string");
    cmd.add( outArg );

    ValueArg<std::string> sufArg("s","suffix","output suffix", false,".nii.gz","string");
    cmd.add( sufArg );

    SwitchArg refineArg("f", "refine", "refine segmentation", false);
    cmd.add( refineArg );

    SwitchArg debugArg("d", "debug", "save debug images", false);
    cmd.add( debugArg );

    SwitchArg biasArg("", "biascorrect", "Preform rough bias correction", false);
    cmd.add( biasArg );

    SwitchArg openSkullArg("", "openSkull", "Assume skull is open (or faint) - turns of a marker dilation", false);
    cmd.add( openSkullArg );

    SwitchArg internalBGMarkArg("", "nointernalDarkMarkers", "Do not use any background marker except the largest one", false);
    cmd.add( internalBGMarkArg );

    std::vector<float> DefSig;
    DefSig.push_back(1.0);
    ValueArg< std::vector<float> > smoothGradArg("", "smoothGradSigma", "smoothing sizes for gradient in refinement step (mm) - include multiple scales in quotes",
                                                 false, DefSig,
                                                 "vector of smoothing sizes");

    // ValueArg<float> smoothGradArg("","smoothGradSigma","smoothing size for gradient in refinement step (mm)", false, 1.0,"float");
    cmd.add(smoothGradArg);

    // shouldn't need to change these
    ValueArg<float> rough_open_radiusArg("", "roughOpenRadius", "size of opening applied to rough mask (mm)", false, 5.0, "float");
    cmd.add(rough_open_radiusArg);

    ValueArg<float> rough_close_radiusArg("", "roughCloseRadius", "size of closing applied to rough mask (mm). Should be bigger than the opening", false, 6.5, "float");
    cmd.add(rough_close_radiusArg);

    //ValueArg<float> refine_grad_sigmaArg("", "refineGradSigma", "variance of smoothing filter in directional gradient computation (mm)", false, 0.25, "float");
    //cmd.add(refine_grad_sigmaArg);

    ValueArg<int> preopeningArg("", "presmoothing", "radius of smoothing (in voxels), applied prior to segmentation. The smoothing is a morphological opening. Useful for restricting leaks, which can happen with high res data, or when arteries are very bright.", false, 0, "int");
    cmd.add(preopeningArg);

    ValueArg<int> threadArg("", "threads", "How many threads to use. Default is max cores", false, -1, "int");
    cmd.add(threadArg);

    ValueArg<float> neckCropArg("", "neckCropDistance", "distance from top of head to a point definitely containing no brain (mm)", false, head_crop_size, "float");
    cmd.add(neckCropArg);

    ValueArg<float> TopToMarkArg("", "topToMarkerDistance", "distance from top of head to the marker centre (mm)", false, TopToMarkCent, "float");
    cmd.add(TopToMarkArg);

    ValueArg<float> markerRadiusArg("", "markerRadius", "size of the box used to mark the brain (mm)", false, markerRadius, "float");
    cmd.add(markerRadiusArg);

    ValueArg<float> bgmarkerErodeRadiusArg("", "bgmarkerErodeRadius", "distance around the brain marker that we ignore when looking for bg markers (mm)", false, mkbgmarker_erode_radius, "float");
    cmd.add(bgmarkerErodeRadiusArg);

    ValueArg<float> refineErodeArg("", "refineErodeRadius", "Size of the grayscale erosion applied to the image prior to refinement (mm)", false, refine_erode_size, "float");
    cmd.add(refineErodeArg);

    ValueArg<float> biasCorrectRadiusArg("", "biasCorrectKernelRadius", "size of the mean kernel filter used in the bias correction (mm)", false, bias_correct_radius, "float");
    cmd.add(biasCorrectRadiusArg);

    ValueArg<float> biasCorrectGammaArg("", "biasCorrectGamma", "A gamma correction applied to the bias field so noise in extremely dark regions doesn't get blown out of all proportion. Values of 0.8 seem to perform OK in images with severe dropout. Otherwise it can be turned off. Leaks of brain segmentation are indications that it might need to be played with", false, bias_gamma, "float");
    cmd.add(biasCorrectGammaArg);

   // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputImPrefix = outArg.getValue();
    CmdLineObj.MaskIm = markArg.getValue();
    CmdLineObj.Suffix = sufArg.getValue();
    CmdLineObj.refine = refineArg.getValue();
    CmdLineObj.biascorrect = biasArg.getValue();
    CmdLineObj.openSkull = openSkullArg.getValue();
    CmdLineObj.smoothGradSigma = smoothGradArg.getValue();
    CmdLineObj.debug = debugArg.getValue();
    CmdLineObj.PreOpen = preopeningArg.getValue();
    CmdLineObj.InternalBGMarkers = !internalBGMarkArg.getValue();
    rough_open_radius  = rough_open_radiusArg.getValue();
    rough_close_radius = rough_close_radiusArg.getValue();
    //refine_grad_sigma  = refine_grad_sigmaArg.getValue();
    CmdLineObj.threads = threadArg.getValue();
    head_crop_size = neckCropArg.getValue();
    TopToMarkCent = TopToMarkArg.getValue();
    markerRadius = markerRadiusArg.getValue();
    mkbgmarker_erode_radius = bgmarkerErodeRadiusArg.getValue();
    bias_correct_radius = biasCorrectRadiusArg.getValue();
    bias_gamma = biasCorrectGammaArg.getValue();
    refine_erode_size = refineErodeArg.getValue();
    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}
////////////////////////////////////////////////////////
template <class TImage>
void writeImDbg(const typename TImage::Pointer Im, std::string filename)
{
  if (wssdebug)
    {
    writeIm<TImage>(Im, wssprefix + "_" + filename + wsssuffix);
    }
}

#include "cropneck.h"

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
/////////////////////////////////
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

#if 1
template <class RawImType, class LabImType>
typename LabImType::Pointer mkMarker(typename RawImType::Pointer t1, int top, int x, int y, float brainMarkRad)
{
  // trialing just banging a blob down as the foreground marker
  typename LabImType::Pointer marker = LabImType::New();
  marker->SetRegions(t1->GetLargestPossibleRegion());
  marker->Allocate();
  marker->CopyInformation(t1);
  marker->FillBuffer(0);

  typename LabImType::IndexType bmcent;
  typename LabImType::SpacingType bmspacing = t1->GetSpacing();
  typename LabImType::PointType bmcentwrld;
  bmcent[0]=x;
  bmcent[1]=y;
  bmcent[2]=top - TopToMarkCent/bmspacing[2];
  t1->TransformIndexToPhysicalPoint(bmcent, bmcentwrld);

  fillBoxMM<LabImType>(marker, bmcentwrld, 1, brainMarkRad, 2*brainMarkRad, brainMarkRad);

  // now compute some statistics of the marker. We will keep create a
  // mask of similar brightness and keep the component connected to
  // the marker.

  std::vector<float> qvals;
  qvals.push_back(0.5);
  qvals.push_back(0.95);
  std::vector<typename RawImType::PixelType> threshrange  = computeImQuantile<RawImType, LabImType>(t1, marker, qvals);
  itk::Instance<itk::BinaryThresholdImageFilter<RawImType, LabImType> > FGSelect;
  FGSelect->SetInput(t1);
  FGSelect->SetUpperThreshold(threshrange[0]*1.25); // possible multiplier here
  FGSelect->SetLowerThreshold(threshrange[0]);

  FGSelect->SetInsideValue(1);
  FGSelect->SetOutsideValue(0);
  // apply an opening to the mask
  itk::Instance< itk::BinaryOpenParaImageFilter <LabImType, LabImType> > ParaOpen;
  ParaOpen->SetInput(FGSelect->GetOutput());
  ParaOpen->SetUseImageSpacing(true);
  ParaOpen->SetRadius(2);

  // combine with the box
  itk::Instance< itk::MaximumImageFilter <LabImType, LabImType, LabImType> > Combine;
  Combine->SetInput(ParaOpen->GetOutput());
  Combine->SetInput2(marker);
  typename LabImType::Pointer comb = Combine->GetOutput();
  comb->Update();
  comb->DisconnectPipeline();
  // blank a slice so we don't get into the danger zone
#if 0
  {
    typename LabImType::RegionType blank = comb->GetLargestPossibleRegion();
    typename LabImType::RegionType::SizeType s = blank.GetSize();
    typename LabImType::RegionType::IndexType index = blank.GetIndex();

    s[2] = 1;
    index[2]= top - 2*(TopToMarkCent/bmspacing[2]);
    blank.SetSize(s);
    blank.SetIndex(index);
    fillRegion<LabImType>(comb, blank, 0);

  }
#endif
  // label the marker and keep the one corresponding to the box
  typedef typename itk::BinaryImageToShapeLabelMapFilter<LabImType> LabellerType;
  typename LabellerType::Pointer Labeller = LabellerType::New();
  Labeller->SetInput(comb);
  Labeller->SetInputForegroundValue(1);
  typedef typename LabellerType::OutputImageType LabelMapType;
  typedef typename LabellerType::OutputImagePointer LabelMapPointerType;
  typedef typename LabelMapType::LabelObjectType LabelObjectType;
  LabelMapPointerType labmap = Labeller->GetOutput();
  labmap->Update();
  labmap->DisconnectPipeline();

  // create a new label map with only the relevant object in it.
  // figure out the label including the box
  typename LabImType::PixelType boxlab = labmap->GetPixel(bmcent);
  std::vector<unsigned long> labelsToRemove;

  for(unsigned int i = 0; i < labmap->GetNumberOfLabelObjects(); i++)
    {
    // Get the ith region
    LabelObjectType* labelObject = labmap->GetNthLabelObject(i);
    if((LabImType::PixelType)labelObject->GetLabel() != boxlab)
      {
      labelsToRemove.push_back(labelObject->GetLabel());
      }
    }

  // Remove all regions that were marked for removal.
  for(unsigned int i = 0; i < labelsToRemove.size(); ++i)
    {
    labmap->RemoveLabel(labelsToRemove[i]);
    }

  // render it
  typedef typename itk::LabelMapToBinaryImageFilter<LabelMapType, LabImType> RenderType;
  typename RenderType::Pointer render = RenderType::New();
  render->SetInput(labmap);
  render->SetBackgroundValue(0);
  render->SetForegroundValue(1);
  marker=render->GetOutput();
  marker->Update();
  marker->DisconnectPipeline();
  writeImDbg<LabImType>(marker, "biggest");
  return(marker);
}

#else
template <class RawImType, class LabImType>
typename LabImType::Pointer mkMarker(typename RawImType::Pointer t1)
{
  // this takes the raw t1, and creates a brain marker
  typename RawImType::Pointer eroded = doErodeMM<RawImType>(t1, mkmarker_gs_erode_size);
  writeImDbg<RawImType>(eroded, "eroded");

  // threshold at 1% level to create a mask, which we will threshold
  // again
  std::vector<float> qvals;
  qvals.push_back(0.01);
  typename RawImType::PixelType onePC = computeImQuantile<RawImType>(eroded, qvals)[0];

  // this is a mask to which we apply our threshold calculation
  itk::Instance<itk::BinaryThresholdImageFilter<RawImType, LabImType> > FGSelect;
  FGSelect->SetInput(eroded);
  FGSelect->SetUpperThreshold(onePC);
  FGSelect->SetInsideValue(0);
  FGSelect->SetOutsideValue(1);
  writeImDbg<LabImType>(FGSelect->GetOutput(), "fgselect");

  itk::Instance< itk::OtsuThresholdImageFilter <RawImType, LabImType, LabImType> >Otsu;

  Otsu->SetOutsideValue(1);
  Otsu->SetInsideValue(0);
  Otsu->SetInput(eroded);
  Otsu->SetMaskImage(FGSelect->GetOutput());
  Otsu->SetMaskValue(1);

  writeImDbg<LabImType>(Otsu->GetOutput(), "fgthresh");

  // a binary erosion to get rid of fine paths to the outside and make
  // the internal marker more conservative
  itk::Instance<itk::BinaryErodeParaImageFilter<LabImType, LabImType> > ParaErode;
  ParaErode->SetInput(Otsu->GetOutput());
  ParaErode->SetUseImageSpacing(true);
  ParaErode->SetRadius(mkmarker_binary_erode_size);

  writeImDbg<LabImType>(ParaErode->GetOutput(), "fgthresherode");

  itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter<LabImType> > SizeFilter;

  SizeFilter->SetInput(ParaErode->GetOutput());
  SizeFilter->SetBackgroundValue(0);
  SizeFilter->SetForegroundValue(1);
  SizeFilter->SetNumberOfObjects(1);
  SizeFilter->SetAttribute("PhysicalSize");

  writeImDbg<LabImType>(SizeFilter->GetOutput(), "biggest");

  typename LabImType::Pointer result = SizeFilter->GetOutput();
  result = SizeFilter->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
#endif
template <class LabImType>
typename LabImType::Pointer mkBGMarker(typename LabImType::Pointer bmarker)
{
  // Turns the foreground marker into a background marker
  itk::Instance<itk::BinaryThresholdImageFilter<LabImType, LabImType> > Inverter;
  Inverter->SetInput(bmarker);
  Inverter->SetUpperThreshold(1);
  Inverter->SetLowerThreshold(1);
  Inverter->SetInsideValue(0);
  Inverter->SetOutsideValue(1);

  // stick an erosion in here - matches the size of the erosion used
  // in making the foreground marker
  itk::Instance< itk::BinaryErodeParaImageFilter <LabImType, LabImType> > ParaErode;
  ParaErode->SetInput(Inverter->GetOutput());
  ParaErode->SetUseImageSpacing(true);
  ParaErode->SetRadius(mkbgmarker_erode_radius);

  // now do an opening to make sure the background marker is nice and smooth
  itk::Instance< itk::BinaryOpenParaImageFilter <LabImType, LabImType> > ParaOpen;
  ParaOpen->SetInput(ParaErode->GetOutput());
  ParaOpen->SetUseImageSpacing(true);
  ParaOpen->SetRadius(mkbgmarker_open_radius);

  itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter<LabImType> > SizeFilter;

  SizeFilter->SetInput(ParaOpen->GetOutput());
  SizeFilter->SetBackgroundValue(0);
  SizeFilter->SetForegroundValue(1);
  SizeFilter->SetNumberOfObjects(1);
  SizeFilter->SetAttribute("PhysicalSize");

  typename LabImType::Pointer result = SizeFilter->GetOutput();
  result = SizeFilter->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}

////////////////////////////////////////////////////////
template <class RawImType, class LabImType>
typename LabImType::Pointer adjustMarker(typename LabImType::Pointer marker,
                                         typename RawImType::Pointer t1,
                                         int bottom,
                                         bool openskull=false,
                                         bool InternalBGMarkers=true,
                                         int origbglabel=1,
                                         int bglabel=2)
{

  typedef typename itk::Image<unsigned char, LabImType::ImageDimension> MaskImType;

  typename RawImType::Pointer t1a = doOpeningMM<RawImType>(t1, 5, 5, 5);

  itk::Instance<itk::OtsuThresholdImageFilter<RawImType, MaskImType, LabImType> > otsu;
  otsu->SetInput( t1a );
  otsu->SetMaskImage(marker);
  otsu->SetNumberOfHistogramBins( 200 );
  otsu->SetMaskValue(origbglabel);
  otsu->SetInsideValue( 1 );
  otsu->SetOutsideValue( 0 );

  otsu->Update();

  itk::Instance<itk::BinaryThresholdImageFilter<RawImType, MaskImType> > Thresh;
  Thresh->SetInput(t1);
  Thresh->SetUpperThreshold(otsu->GetThreshold());
  Thresh->SetInsideValue( 1 );
  Thresh->SetOutsideValue( 0 );



  typename MaskImType::Pointer dark = Thresh->GetOutput();

  writeImDbg<LabImType>(marker, "adj_marker");

  // pull bglabel out of the marker
  itk::Instance<itk::BinaryThresholdImageFilter<LabImType, MaskImType> > Select;
  Select->SetInput(marker);
  Select->SetUpperThreshold(origbglabel);
  Select->SetLowerThreshold(origbglabel);
  Select->SetInsideValue(1);
  Select->SetOutsideValue(0);

  itk::Instance<itk::MaskImageFilter<MaskImType, MaskImType> > Masker;
  Masker->SetInput(Select->GetOutput());
  Masker->SetInput2(dark);
  typename MaskImType::Pointer m = Masker->GetOutput();
  m->Update();
  m->DisconnectPipeline();

  // remove the posterior part cropped bit from the marker
  {
  if (bottom > 0)
    {
    // blank all slices below
    typename MaskImType::RegionType blank = m->GetLargestPossibleRegion();
    typename RawImType::RegionType::SizeType s = blank.GetSize();
    s[2] = bottom;
    s[1] *= 0.5;
    blank.SetSize(s);
    fillRegion<MaskImType>(m, blank, 0);


    }
  }

  writeImDbg<MaskImType>(Masker->GetOutput(), "adj_masker");
  // pull bglabel out of the marker

  const float radius = conf_adjust_radius;
  //const float radius2 = conf_adjust_radius2;

  // small spherical erosion here
  itk::Instance< itk::BinaryErodeParaImageFilter <MaskImType, MaskImType> > ParaErode;
  //ParaErode->SetInput(Masker->GetOutput());
  ParaErode->SetInput(m);
  ParaErode->SetUseImageSpacing(true);
  ParaErode->SetRadius(radius);
  writeImDbg<MaskImType>(ParaErode->GetOutput(), "erode1");


  // we are going to process the largest and other parts separately -
  // largest is background, others are holes in tissue. Use different dilations
  // keep the largest part
  itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter<MaskImType> > SizeFilter;
  SizeFilter->SetInput(ParaErode->GetOutput());
  SizeFilter->SetBackgroundValue(0);
  SizeFilter->SetForegroundValue(1);
  SizeFilter->SetNumberOfObjects(1);
  SizeFilter->SetAttribute("PhysicalSize");

  writeImDbg<MaskImType>(SizeFilter->GetOutput(), "adj_largest");
  typename MaskImType::Pointer sizefiltout = SizeFilter->GetOutput();
  sizefiltout->Update();
  sizefiltout->DisconnectPipeline();

  if (bottom > 0)
    {
    // fill the anterior parts - these can get removed by the
    // filtering, but we need them there
    typename LabImType::RegionType blank = sizefiltout->GetLargestPossibleRegion();
    typename LabImType::RegionType::IndexType index = blank.GetIndex();
    typename LabImType::RegionType::SizeType s = blank.GetSize();
    index[1]=s[1]*0.5;
    s[1] -= (index[1] + 1);
    s[2] = bottom+1;
    blank.SetIndex(index);
    blank.SetSize(s);
    fillRegion<MaskImType>(sizefiltout, blank, 1);

    // need to fill left and right sides too. Leave the middle third
    // blank
    blank = sizefiltout->GetLargestPossibleRegion();
    index = blank.GetIndex();
    s = blank.GetSize();
    s[2] = bottom + 1;
    s[0] *= 0.33;
    blank.SetIndex(index);
    blank.SetSize(s);
    fillRegion<MaskImType>(sizefiltout, blank, 1);

    blank = sizefiltout->GetLargestPossibleRegion();
    index = blank.GetIndex();
    s = blank.GetSize();
    index[0]=s[0]*0.66;
    s[2] = bottom + 1;
    s[0] -= index[0]+1;
    blank.SetIndex(index);
    blank.SetSize(s);
    fillRegion<MaskImType>(sizefiltout, blank, 1);



    }

  // now dilate it by the erosion size + a tiny bit
  int skullmargin=1;
  int skullmargin2=0;
  if (openskull)
    {
    std::cout << "Openskull modification" << std::endl;
    skullmargin=-1;
    skullmargin2=-1;
    }
  itk::Instance< itk::BinaryDilateParaImageFilter<MaskImType, MaskImType> > ParaDilate;
  // ParaDilate->SetInput(SizeFilter->GetOutput());
  ParaDilate->SetInput(sizefiltout);
  ParaDilate->SetRadius(radius+skullmargin);
  ParaDilate->SetUseImageSpacing(true);
#if 1
  // dilate the other parts by a little less than the erosion size, combine with a max
  itk::Instance< itk::BinaryDilateParaImageFilter<MaskImType, MaskImType> > ParaDilate2;
  ParaDilate2->SetInput(ParaErode->GetOutput());
  ParaDilate2->SetRadius(radius + skullmargin2);
  ParaDilate2->SetUseImageSpacing(true);

  itk::Instance< itk::MaximumImageFilter <MaskImType, MaskImType, MaskImType> > CombDilate;
  CombDilate->SetInput(ParaDilate->GetOutput());
  CombDilate->SetInput2(ParaDilate2->GetOutput());
#else
  itk::Instance< itk::MaximumImageFilter <MaskImType, MaskImType, MaskImType> > CombDilate;
  CombDilate->SetInput(ParaDilate->GetOutput());
  CombDilate->SetInput2(ParaErode->GetOutput());
#endif
  writeImDbg<MaskImType>(CombDilate->GetOutput(), "pdcomb");
  writeImDbg<MaskImType>(ParaDilate->GetOutput(), "pd1");
  writeImDbg<MaskImType>(ParaDilate2->GetOutput(), "pd2");

  // at this point we have a marker sitting just inside the outer
  // scalp. We need to deal with patients with thick skulls, so work on
  // the difference between our original marker and this one.

  typename MaskImType::Pointer BGmarker = CombDilate->GetOutput();
  if (!InternalBGMarkers)
    {
    BGmarker = ParaDilate->GetOutput();
    }


  itk::Instance< itk::BinaryThresholdImageFilter <MaskImType, MaskImType> > BThresh;
  BThresh->SetInput(BGmarker);
  BThresh->SetUpperThreshold(0);
  BThresh->SetInsideValue(0);
  BThresh->SetOutsideValue(bglabel);

  // now put the background marker back into the original image - new
  // bg marker is a subset of the old one.
  // knock out the bg label and put in the new one
  itk::Instance< itk::ChangeLabelImageFilter<LabImType, LabImType> > BGRemove;
  BGRemove->SetInput(marker);
  BGRemove->SetChange(origbglabel,0);

  itk::Instance< itk::MaximumImageFilter <MaskImType, LabImType, LabImType> > MaxFilt;
  MaxFilt->SetInput(BThresh->GetOutput());
  MaxFilt->SetInput2(BGRemove->GetOutput());

  typename LabImType::Pointer result = MaxFilt->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  // place a marker around all edges, except the bottom - no good if the brain is clipped.
  {
    typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<LabImType> FaceCalculatorType;
    FaceCalculatorType faceCalculator;
    typename MaskImType::SizeType mradius;
    mradius.Fill(1);
    typename FaceCalculatorType::FaceListType faceList;
    faceList = faceCalculator(result, result->GetLargestPossibleRegion(), mradius);
    typename FaceCalculatorType::FaceListType::iterator fit;
    fit = faceList.begin();
    ++fit;
    for (;fit != faceList.end(); ++fit)
      {
      if ((fit->GetIndex()[2] == 0) && (fit->GetSize()[2] == 1))
        {
        // this is the bottom one - don't fill it all
        typename LabImType::RegionType fitcopy = result->GetLargestPossibleRegion();
        typename LabImType::RegionType::IndexType index = fit->GetIndex();
        typename LabImType::RegionType::SizeType s = fit->GetSize();
        index[1]=s[1]*0.5;
        s[1] -= (index[1] + 1);
        fitcopy.SetIndex(index);
        fitcopy.SetSize(s);
        fillRegion<LabImType>(result, fitcopy, bglabel);
        }
      fillRegion<LabImType>(result, *fit, bglabel);
      }
  }
  // put the bit we cropped from the bottom back into the marker
  // but leave the upper-most slice blank. This is just a speed thing
  {
  if (bottom > 1)
    {
    // fill all slices below
    typename LabImType::RegionType blank = result->GetLargestPossibleRegion();
    typename LabImType::RegionType::SizeType s = blank.GetSize();
    s[2] = bottom-1;
    blank.SetSize(s);
    fillRegion<LabImType>(result, blank, bglabel);

    }

  }
  return(result);
}

/////////////////////////////////
#if 0
template <class LabImType>
typename LabImType::Pointer ErodeKeepBigDilate(typename LabImType::Pointer mask, float radius, float ReduceDilate=0.0)
{
  // Now for an erode/keep biggest/dilate combination to get rid of the crap
  itk::Instance< itk::BinaryErodeParaImageFilter <LabImType, LabImType> > ParaErode2;
  ParaErode2->SetInput(mask);
  ParaErode2->SetUseImageSpacing(true);
  ParaErode2->SetRadius(radius);
  // keep the largest part
  itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter<LabImType> > SizeFilter;

  SizeFilter->SetInput(ParaErode2->GetOutput());
  SizeFilter->SetBackgroundValue(0);
  SizeFilter->SetForegroundValue(1);
  SizeFilter->SetNumberOfObjects(1);
  SizeFilter->SetAttribute("PhysicalSize");

  // now dilate it by the erosion size
  itk::Instance< itk::BinaryDilateParaImageFilter<LabImType, LabImType> > ParaDilate;
  ParaDilate->SetInput(SizeFilter->GetOutput());
  ParaDilate->SetRadius(radius - ReduceDilate);
  ParaDilate->SetUseImageSpacing(true);

  typename LabImType::Pointer result = ParaDilate->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  return(result);
}
#endif
////////////////////////////////////////////////////////
template <class ImType>
typename ImType::Pointer mkControlImage(typename ImType::Pointer im)
{
  typename ImType::PixelType ImMax;

  itk::Instance < itk::StatisticsImageFilter<ImType> > StatsFilt;
  StatsFilt->SetInput(im);
  StatsFilt->Update();
  ImMax = StatsFilt->GetMaximum();

  itk::Instance< itk::ShiftScaleImageFilter<ImType, ImType> > Inverter;
  Inverter->SetInput(im);
  Inverter->SetShift(-1 * ImMax);
  Inverter->SetScale(-1);
  typename ImType::Pointer result = Inverter->GetOutput();
  result->Update();
  result->DisconnectPipeline();
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
////////////////////////////////////////////////////////////////

template <class LabImType, class RawImType, class OutImType>
typename OutImType::Pointer maskedSmooth(typename LabImType::Pointer mask,
                                         typename RawImType::Pointer raw,
                                         float x, float y, float z)
{
#if 1
  typedef typename itk::Image<float, RawImType::ImageDimension> FloatImType;
  typedef typename itk::BoxMeanImageFilter<RawImType, FloatImType> SmoothRawType;
  typedef typename itk::BoxMeanImageFilter<LabImType, FloatImType> SmoothMaskType;
  typename RawImType::SpacingType spacing = raw->GetSpacing();

  typename SmoothRawType::RadiusType rad;
  rad[0]=std::max(1, (int)(x/spacing[0]));
  rad[1]=std::max(1, (int)(y/spacing[1]));
  rad[2]=std::max(1, (int)(z/spacing[2]));

  // first apply the mask (may already have been done
  itk::Instance<itk::MaskImageFilter<RawImType, LabImType> > masker;
  masker->SetInput(raw);
  masker->SetInput2(mask);

  typename SmoothRawType::Pointer smooth1 = SmoothRawType::New();
  smooth1->SetInput(masker->GetOutput());
  smooth1->SetRadius(rad);

  typename SmoothMaskType::Pointer masksmooth = SmoothMaskType::New();
  masksmooth->SetInput(mask);
  masksmooth->SetRadius(rad);

  itk::Instance<itk::DivideImageFilter<FloatImType, FloatImType, OutImType> > Reweight;
  Reweight->SetInput(smooth1->GetOutput());
  Reweight->SetInput2(masksmooth->GetOutput());

  // a second masker to get rid of NAs
  itk::Instance<itk::MaskImageFilter<OutImType, LabImType> > masker2;
  masker2->SetInput(Reweight->GetOutput());
  masker2->SetInput2(mask);
  typename OutImType::Pointer res = masker2->GetOutput();
  res->Update();
  res->DisconnectPipeline();
  return(res);

#else
  // median version
  typedef typename itk::MaskedRankImageFilter< RawImType, LabImType, OutImType > LocalMedType;
  typename LocalMedType::Pointer BGEst = LocalMedType::New();
  typename LocalMedType::FlatKernelType::RadiusType rad;
  typename RawImType::SpacingType spacing = raw->GetSpacing();
  rad[0]=std::max(1, (int)(x/spacing[0]));
  rad[1]=std::max(1, (int)(y/spacing[1]));
  rad[2]=std::max(1, (int)(z/spacing[2]));
  typename LocalMedType::FlatKernelType kernel = LocalMedType::FlatKernelType::Box(rad);

  BGEst->SetInput(raw);
  BGEst->SetMaskImage(mask);
  BGEst->SetMaskValue(1);
  BGEst->SetKernel(kernel);

  typename OutImType::Pointer res = BGEst->GetOutput();
  res->Update();
  res->DisconnectPipeline();
  return(res);
#endif

}

////////////////////////////////////////////////////////////////

// refinement process
template <class LabImType, class RawImType>
typename LabImType::Pointer doRefineCSF(typename LabImType::Pointer mask,
                                        typename RawImType::Pointer raw0,
                                        float radius,
                                        std::vector<float> sigma)
{
  typedef typename RawImType::Pointer PRawImType;
  typedef typename LabImType::Pointer PLabImType;

  typedef typename itk::Image<float, RawImType::ImageDimension> FloatImType;
  typedef typename FloatImType::Pointer PFloatImType;

  // strategy this time is to find dark areas in an eroded version of
  // the input. We are also going to use the eroded input to produce
  // the control image and then dilate the result.


  // compute the median.
  std::vector<float> quantlevs;
  quantlevs.push_back(0.5);

  std::vector<typename RawImType::PixelType> rthresh = computeImQuantile<RawImType, LabImType>(raw0, mask, quantlevs);


  // prefiltering
  // reset brightness of everything above median to remove internal
  // edges.
  // now we apply a thresholding to reset WM intensity
  itk::Instance<itk::ThresholdImageFilter<RawImType> > resetWM;
  resetWM->SetInput(raw0);
  resetWM->ThresholdAbove(rthresh[0]);
  resetWM->SetOutsideValue(rthresh[0]);

  // want a small erosion, but need to be careful about voxel size.
  float ex, ey, ez;
  ex = std::max(refine_erode_size, (float)raw0->GetSpacing()[0]);
  ey = std::max(refine_erode_size, (float)raw0->GetSpacing()[1]);
  ez = std::max(refine_erode_size, (float)raw0->GetSpacing()[2]);
  // eroded version to eliminate thin edges
  PRawImType rawerode0 = doErodeMM<RawImType>(resetWM->GetOutput(), ex, ey, ez);
  // mask using stage 1 result
  itk::Instance< itk::MaskImageFilter< RawImType, LabImType> > ApplyStage1B;
  ApplyStage1B->SetInput(rawerode0);
  ApplyStage1B->SetInput2(mask);

  PRawImType rawerode = ApplyStage1B->GetOutput();

  // this is an eroded, masked and upper threshold image
  PRawImType raw = rawerode;

  writeImDbg<RawImType>(raw, "csfprefilt");

  PLabImType eroded, erodedthin;
  PLabImType ring = mkMaskBorder<LabImType>(mask, radius, true, eroded);
  PLabImType thinring = mkMaskBorder<LabImType>(mask, radius/2, true, erodedthin);

  writeImDbg<LabImType>(thinring, "thinring");

  // construct new markers
  itk::Instance< itk::BinaryThresholdImageFilter <LabImType, LabImType> > Invert;
  Invert->SetInput(mask);
  Invert->SetInsideValue(0);
  Invert->SetOutsideValue(2);
  Invert->SetUpperThreshold(1);
  Invert->SetLowerThreshold(1);

  // make sure interior marker is only in bright areas (>= median used
  // for thresholding
  itk::Instance<itk::BinaryThresholdImageFilter<RawImType, LabImType> > BThresh;
  BThresh->SetInput(raw);
  BThresh->SetLowerThreshold(rthresh[0]);
  BThresh->SetInsideValue(1);
  BThresh->SetOutsideValue(0);

  // now compute an approximate white matter median.
  std::vector<typename RawImType::PixelType> rthreshWM = computeImQuantile<RawImType, LabImType>(raw0, BThresh->GetOutput(), quantlevs);

  if (wssdebug)
    {
    std::cout << "Refine bright marker threshold = " << (float)rthreshWM[0] * 1.1 << std::endl;
    }

#if 0
  // dilate a bit to push over WM/GM boundary - hard to do well.
// possibly a closing instead
  itk::Instance< itk::BinaryDilateParaImageFilter<LabImType, LabImType> > ParaDilateFG;
  ParaDilateFG->SetInput(BThresh->GetOutput());
  ParaDilateFG->SetRadius(sigma[sigma.size() - 1]*2*2);
  ParaDilateFG->SetUseImageSpacing(true);

  itk::Instance<itk::MaskImageFilter< LabImType, LabImType > > SeedMask;
  SeedMask->SetInput(ParaDilateFG->GetOutput());
  SeedMask->SetInput2(eroded);
#else
  itk::Instance<itk::MaskImageFilter< LabImType, LabImType > > SeedMask;
  SeedMask->SetInput(BThresh->GetOutput());
  SeedMask->SetInput2(eroded);

#endif


  PLabImType newbgmarkers, newbrightmarkers;
  {
  // local (large) mean filter of stage 1 mask to set local
  // thresholds - 15voxels
  typename FloatImType::Pointer BGEst = maskedSmooth<LabImType, RawImType, FloatImType>(mask, resetWM->GetOutput(), 30, 30, 30);

  itk::Instance<itk::DivideImageFilter<RawImType, FloatImType, FloatImType> > Divider;
  Divider->SetInput(resetWM->GetOutput());
  Divider->SetInput2(BGEst);

  // debug only
  itk::Instance<itk::MaskImageFilter<FloatImType, LabImType> > DivMask;
  DivMask->SetInput(Divider->GetOutput());
  DivMask->SetInput2(mask);
  writeImDbg<FloatImType>(DivMask->GetOutput(), "refdiv");

  // find any voxels darker than 0.6 * mean
  itk::Instance<itk::BinaryThresholdImageFilter<FloatImType, LabImType> > SelDark;
  SelDark->SetInput(Divider->GetOutput());
  // magic number - consider providing configuration option
  SelDark->SetUpperThreshold(0.6);
  SelDark->SetInsideValue(2);
  SelDark->SetOutsideValue(0);

  // locally bright things
  itk::Instance<itk::BinaryThresholdImageFilter<RawImType, LabImType> > SelBright;
  SelBright->SetInput(raw0);
  SelBright->SetUpperThreshold(rthreshWM[0]* 1.1);
  SelBright->SetInsideValue(0);
  SelBright->SetOutsideValue(2);

  // mask with the thinring
  itk::Instance<itk::MaskImageFilter<LabImType, LabImType> > FineRingMask;
  FineRingMask->SetInput(SelBright->GetOutput());
  FineRingMask->SetInput2(thinring);

  itk::Instance< itk::BinaryShapeOpeningImageFilter<LabImType> > BrightAreaOpening;
  BrightAreaOpening->SetInput(FineRingMask->GetOutput());
  BrightAreaOpening->SetLambda(5);
  BrightAreaOpening->SetAttribute("PhysicalSize");
  BrightAreaOpening->SetForegroundValue(2);


  // blank the bottom
  PLabImType toblank = BrightAreaOpening->GetOutput();
  toblank->Update();
  toblank->DisconnectPipeline();
  {
  // another magic number - anatomy dependent
  float cereb_adjust = 100.0; //mm
  float ringthickness = radius;
  // discard markers that are low - we tend to pick up cerebellum
  // textures with the current marker method. Use the eroded stage 1
  // mask as a guide. Find the bottom and go up a bit.
  typedef typename itk::BinaryImageToShapeLabelMapFilter<LabImType> LabellerType;
  typename LabellerType::Pointer Labeller = LabellerType::New();
  Labeller->SetInput(eroded);
  Labeller->SetInputForegroundValue(1);
  typedef typename LabellerType::OutputImageType LabelMapType;
  typedef typename LabellerType::OutputImagePointer LabelMapPointerType;
  typedef typename LabelMapType::LabelObjectType LabelObjectType;
  LabelMapPointerType labmap = Labeller->GetOutput();
  labmap->Update();
  labmap->DisconnectPipeline();
  LabelObjectType * labelObject = labmap->GetLabelObject(1);
  typename LabelObjectType::RegionType bbox = labelObject->GetBoundingBox();

  // lower bound is the interesting bit
  int newz = bbox.GetIndex()[2] + (cereb_adjust-ringthickness)/toblank->GetSpacing()[2];

  typename LabImType::RegionType blank = toblank->GetLargestPossibleRegion();
  typename LabImType::RegionType::SizeType s = blank.GetSize();
  s[2] = newz;
  blank.SetSize(s);
  fillRegion<LabImType>(toblank, blank, 0);
  writeImDbg<LabImType>(toblank, "brightblanked");
  }

  itk::Instance<itk::MaximumImageFilter<LabImType, LabImType, LabImType> > CombBGMark;
  CombBGMark->SetInput(SelDark->GetOutput());
  CombBGMark->SetInput2(toblank);
  newbrightmarkers=toblank;
  itk::Instance<itk::MaskImageFilter<LabImType, LabImType > > MaskRing;
  MaskRing->SetInput(CombBGMark->GetOutput());
  MaskRing->SetInput2(ring);

  // size filter
  itk::Instance< itk::BinaryShapeOpeningImageFilter<LabImType> > AreaOpening;
  AreaOpening->SetInput(MaskRing->GetOutput());
  AreaOpening->SetLambda(10);
  AreaOpening->SetAttribute("PhysicalSize");
  AreaOpening->SetForegroundValue(2);

  newbgmarkers=AreaOpening->GetOutput();
  newbgmarkers->Update();
  newbgmarkers->DisconnectPipeline();

  writeImDbg<LabImType>(newbgmarkers, "divring");
  }
  itk::Instance< itk::MaximumImageFilter <LabImType, LabImType, LabImType> > CombMarker0;
  CombMarker0->SetInput(newbgmarkers);
  CombMarker0->SetInput2(Invert->GetOutput());
  itk::Instance< itk::MaximumImageFilter <LabImType, LabImType, LabImType> > CombMarker;
  CombMarker->SetInput(CombMarker0->GetOutput());
  CombMarker->SetInput2(SeedMask->GetOutput());
  writeImDbg<LabImType>(CombMarker->GetOutput(), "refcsfmark");

  // dispose this
  erodedthin=0;

  // construct the control image. This is a combination of the
  // original border intensity and a scale space gradient.
  // We'll use a morphological gradient as the base so we can
  // eliminate influence of the stage 1 mask.

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
  //typename RawImType::Pointer gradR = doOldGrad<RawImType,
  //LabImType>(modraw, mask, 0.2);
#if 1
  typename RawImType::Pointer grad0 = doGradientOuter<RawImType>(raw,1);
  // mask
  itk::Instance <itk::MaskImageFilter<RawImType, LabImType> > MaskGrad;
  MaskGrad->SetInput(grad0);
  MaskGrad->SetInput2(mask);
  typename RawImType::Pointer gradR = scaleSpaceSmooth<RawImType>(MaskGrad->GetOutput(), sigma);
#else
  typename RawImType::Pointer grad0 = doGradientMasked<RawImType>(raw, 1);
  typename RawImType::Pointer gradR = scaleSpaceSmooth<RawImType>(grad0, sigma);
#endif
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

// select the brain label and watershed line
  itk::Instance<itk::BinaryThresholdImageFilter<LabImType, LabImType> > Select;
  Select->SetInput(WatershedFiltA->GetOutput());
  Select->SetUpperThreshold(1);
//  Select->SetLowerThreshold(0);
  Select->SetLowerThreshold(1);
  Select->SetInsideValue(1);
  Select->SetOutsideValue(0);
  writeImDbg<LabImType>(WatershedFiltA->GetOutput(), "refcsfwshed");

  // This needs to match the original erosion
  typename LabImType::Pointer result0 = doDilateMM<LabImType>(Select->GetOutput(), ex, ey, ez);

  // only want to remove the bright markers
  itk::Instance< itk::MaskNegatedImageFilter<LabImType, LabImType> > RemoveMarkers;
  RemoveMarkers->SetInput(result0);
  RemoveMarkers->SetInput2(newbrightmarkers);
  typename LabImType::Pointer result = RemoveMarkers->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  writeImDbg<LabImType>(newbrightmarkers, "newbm");

  return(result);
}
///////////////////////////////////////////////
#include "biascorrection.h"
///////////////////////////////////////////////
template <class PixType, class LabPixType, int dim>
void doWatershed(const CmdLineType &CmdLineObj)
{
  // Marker based watershed scalping algorithm
  // Phase 1 is marker generation.
  // Phase 2 is watershed on inverse brightness image
  // Phase 3 is generation of gradient based control surface.
  // Phase 4 is another watershed
  typedef typename itk::Image<PixType, dim> RawImType;
  typedef typename itk::Image<LabPixType, dim> LabImType;
  typedef typename RawImType::Pointer RawImPointer;
  typedef typename LabImType::Pointer LabImPointer;

  RawImPointer raworig = readIm<RawImType>(CmdLineObj.InputIm);
  RawImPointer preopen;
  if (CmdLineObj.biascorrect)
    {
    std::cout << "Applying rough bias correction" << std::endl;
    //raworig = crappyBiasCorrection<RawImType>(raworig, 30.0);
    raworig = crappyBiasCorrectionB<RawImType>(raworig, bias_correct_radius, bias_gamma);
    writeImDbg<RawImType>(raworig, "biascor");
    }
  // do a little opening - not always needed
  if (CmdLineObj.PreOpen > 0)
    {
    std::cout << "Applying smoothing prefilter" << std::endl;
    preopen=raworig;
    raworig = doOpening<RawImType>(raworig, CmdLineObj.PreOpen);
    // std::vector<float> smth;
    // typename RawImType::SpacingType sp = raworig->GetSpacing();
    // smth.push_back(CmdLineObj.PreOpen * std::min(sp[0], std::min(sp[1], sp[2])));
    // raworig = scaleSpaceSmooth<RawImType>(raworig, smth);
    writeImDbg<RawImType>(raworig, "rawopen");
    }
  int bottom, top, x, y;
  cropNeck2<RawImType>(raworig, head_crop_size, 35, top, bottom, x, y);
  //cropNeck<RawImType>(raworig, 180, top, bottom);
  LabImPointer newmarker, newbgmarker;
  {
  //LabImPointer marker = mkMarker<RawImType, LabImType>(raworig);
  LabImPointer marker = mkMarker<RawImType, LabImType>(raworig, top, x, y, markerRadius);
  LabImPointer InitBGmarker = mkBGMarker<LabImType>(marker);
  writeImDbg<LabImType>(InitBGmarker, "initbg");
  // here's where we fiddle with the marker image. The basic problem
  // is the distance between the edge of the background marker and the
  // brain. We want the background marker to lie along the scalp, but
  // it is pretty much impossible to achieve this reliably with
  // registration alone. Instead we'll start with an oversize
  // background marker and prune it back by thresholding.
  //newmarker = adjustMarker<RawImType,LabImType>(Reorient->GetOutput(), raworig);
  newbgmarker = adjustMarker<RawImType,LabImType>(InitBGmarker, raworig, bottom, CmdLineObj.openSkull, CmdLineObj.InternalBGMarkers, 1);

  itk::Instance< itk::MaximumImageFilter <LabImType, LabImType, LabImType> > CombMarker;
  CombMarker->SetInput(marker);
  CombMarker->SetInput2(newbgmarker);
  newmarker = CombMarker->GetOutput();
  newmarker->Update();
  newmarker->DisconnectPipeline();

  }

  RawImPointer raw = raworig;
  LabImPointer newmarker2 = newmarker;

  writeImDbg<LabImType>(newmarker2, "marker");

  LabImPointer segPhase1;
  {
  RawImPointer control = mkControlImage<RawImType>(raw);
  //RawImPointer gradcontrol = doOldGrad<RawImType, LabImType>(raw, 0, 0.5);
  //RawImPointer control = scaleSpaceSmooth<RawImType>(gradcontrol, CmdLineObj.smoothGradSigma);

  writeImDbg<RawImType>(control, "control");
  itk::Instance<itk::MorphologicalWatershedFromMarkersImageFilter<RawImType, LabImType> > WatershedFiltA;
  WatershedFiltA->SetInput(control);
  WatershedFiltA->SetMarkerImage(newmarker2);
  //WatershedFiltA->SetMarkWatershedLine(false);
  WatershedFiltA->SetMarkWatershedLine(true);
  segPhase1 = WatershedFiltA->GetOutput();
  segPhase1->Update();
  segPhase1->DisconnectPipeline();
  }

  // select the brain label
  itk::Instance<itk::BinaryThresholdImageFilter<LabImType, LabImType> > Select;
  Select->SetInput(segPhase1);
  Select->SetUpperThreshold(1);
  Select->SetLowerThreshold(0);
  Select->SetInsideValue(1);
  Select->SetOutsideValue(0);

  itk::Instance<itk::BinaryOpenParaImageFilter<LabImType, LabImType> > ParaOpen;
  ParaOpen->SetInput(Select->GetOutput());
  ParaOpen->SetUseImageSpacing(true);
  ParaOpen->SetRadius(rough_open_radius);

  itk::Instance<itk::BinaryCloseParaImageFilter<LabImType, LabImType> > ParaClose;
  ParaClose->SetInput(ParaOpen->GetOutput());
  ParaClose->SetUseImageSpacing(true);
  ParaClose->SetRadius(rough_close_radius);

  // keep the biggest component
  itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter<LabImType> > SizeFilter;
  SizeFilter->SetInput(ParaClose->GetOutput());
  SizeFilter->SetBackgroundValue(0);
  SizeFilter->SetForegroundValue(1);
  SizeFilter->SetNumberOfObjects(1);
  SizeFilter->SetAttribute("PhysicalSize");


  writeIm<LabImType>(SizeFilter->GetOutput(), CmdLineObj.OutputImPrefix + "_brainmask" + CmdLineObj.Suffix);
  writeIm<LabImType>(Select->GetOutput(), CmdLineObj.OutputImPrefix + "_brainmaskrough" + CmdLineObj.Suffix);
  writeImDbg<LabImType>(ParaOpen->GetOutput(), "paraopen");

  typename LabImType::Pointer refinedCSF = 0;
  if (CmdLineObj.refine)
    {
    // the ground truth for the SVE dataset appears to exclude
    // peripheral CSF. Can we apply a simple post processing step to
    // do the same
    //refinedCSF = doRefineCSF<LabImType, RawImType>(refined, raw,
    // 10);
    if (CmdLineObj.PreOpen > 0)
      {
      // reload raw for stage 2.
      raw=preopen;
      }
    typename LabImType::Pointer refinedCSF0 = doRefineCSF<LabImType, RawImType>(SizeFilter->GetOutput(), raw, 10, CmdLineObj.smoothGradSigma);

    // cleaning up of refined segmentation - smooth with a closing
    itk::Instance<itk::BinaryCloseParaImageFilter<LabImType, LabImType> > ParaCloseFinal;
    ParaCloseFinal->SetInput(refinedCSF0);
    ParaCloseFinal->SetUseImageSpacing(true);
    ParaCloseFinal->SetRadius(rough_close_radius);
    // fill holes
    itk::Instance<itk::BinaryFillholeImageFilter<LabImType> > FillHoles;
    FillHoles->SetInput(ParaCloseFinal->GetOutput());
    FillHoles->SetForegroundValue(1);
    refinedCSF = FillHoles->GetOutput();
    refinedCSF->Update();
    refinedCSF->DisconnectPipeline();
#if 0
    if (CmdLineObj.SVEdilate)
      {
      // SVE seems to put edges a touch outside the gradient maximum -
      // this is a hack to increase the mask size
      refinedCSF = doSVEHack<LabImType, RawImType>(refinedCSF, raw, 1);
      }
#endif
    writeIm<LabImType>(refinedCSF, CmdLineObj.OutputImPrefix + "_refinedmasksmoothed2" + CmdLineObj.Suffix);
    writeIm<LabImType>(refinedCSF0, CmdLineObj.OutputImPrefix + "_refinedmask2" + CmdLineObj.Suffix);

    }

  // reload the original, as we might have messed with it.
  raworig = readIm<RawImType>(CmdLineObj.InputIm);

  itk::Instance< itk::MaskImageFilter <RawImType, LabImType, RawImType> > Masker;
  Masker->SetInput(raworig);
  Masker->SetInput2(ParaClose->GetOutput());

  writeIm<RawImType>(Masker->GetOutput(), CmdLineObj.OutputImPrefix + "_brain" + CmdLineObj.Suffix);
  if (CmdLineObj.refine)
    {
    Masker->SetInput2(refinedCSF);
    writeIm<RawImType>(Masker->GetOutput(), CmdLineObj.OutputImPrefix + "_brainrefined" + CmdLineObj.Suffix);
    }
}
////////////////////////////////////////////////////////
int main(int argc, char * argv[])
{

//  itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);

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
  // if (!readImageInfo(CmdLineObj.MarkerIm, &MarkerComponentType, &dim2))
  //   {
  //   std::cerr << "Failed to open " << CmdLineObj.MarkerIm << std::endl;
  //   return(EXIT_FAILURE);
  //   }

  // if (dim1 != dim2)
  //   {
  //   std::cerr << "Image dimensions must match " << dim1 << " " << dim2 << std::endl;
  //   return(EXIT_FAILURE);

  //   }

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
