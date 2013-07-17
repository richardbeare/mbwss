// #define the macros in the files including this

switch (ComponentType)
  {
  case (itk::ImageIOBase::UCHAR):
    doWatershed<unsigned char, WSMARKTYPE, WSDIM>(CmdLineObj);
    break;
  case (itk::ImageIOBase::USHORT):
    doWatershed<unsigned short, WSMARKTYPE, WSDIM>(CmdLineObj);
    break;
  case (itk::ImageIOBase::SHORT):
    doWatershed<short, WSMARKTYPE, WSDIM>(CmdLineObj);
    break;
  case (itk::ImageIOBase::INT):
    doWatershed<short, WSMARKTYPE, WSDIM>(CmdLineObj);
    break;
  case (itk::ImageIOBase::FLOAT):
    doWatershed<float, WSMARKTYPE, WSDIM>(CmdLineObj);
    break;
  default:
    std::cerr << "Unsupported pixel type" << std::endl;
    return(EXIT_FAILURE);
    break;
  }
