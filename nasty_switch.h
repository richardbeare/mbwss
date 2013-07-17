
switch (MarkerComponentType)
    {
    case (itk::ImageIOBase::UCHAR):
#define WSMARKTYPE unsigned char
#include "watershed_switch.h"
#undef WSMARKTYPE
      break;
    case (itk::ImageIOBase::SHORT):
#define WSMARKTYPE short
#include "watershed_switch.h"
#undef WSMARKTYPE
      break;
    case (itk::ImageIOBase::INT):
#define WSMARKTYPE int
#include "watershed_switch.h"
#undef WSMARKTYPE
      break;
    default:
      std::cerr << "Unsupported marker pixel type" << std::endl;
      return(EXIT_FAILURE);
    }

