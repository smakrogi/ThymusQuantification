#include "IO_Utils_ITK.h"

#include "itkImageRegionIteratorWithIndex.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkArchetypeSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRegionOfInterestImageFilter.h"

#include "itkOrientImageFilter.h"
#include <string>
#include <vector>

#include "itksys/SystemTools.hxx"


template <class ImageType>
itk::SmartPointer<ImageType>
ReorientImage( itk::SmartPointer<ImageType> inputImage,
               itk::SpatialOrientation::ValidCoordinateOrientationFlags desiredOrientation )
{
  // now orient and cast to output orientation and type
  typedef itk::OrientImageFilter< ImageType, ImageType > OrientingFilterType;
  typename OrientingFilterType::Pointer orientingFilter = OrientingFilterType::New();  
  orientingFilter->SetInput( inputImage );
  // get input direction from image
  orientingFilter->UseImageDirectionOn();
  orientingFilter->SetDesiredCoordinateOrientation( desiredOrientation );  
  orientingFilter->Update();
  typename ImageType::Pointer output = orientingFilter->GetOutput();
  orientingFilter = 0;

  return output;
}


void ZFlipIfNecessary( OrientedInputImageType::Pointer refImage, 
		       OrientedMaskImageType::Pointer image )
{

  OrientedInputImageType::DirectionType refDirection = refImage->GetDirection();
  OrientedMaskImageType::DirectionType imgDirection = image->GetDirection();

  float flip = 0;
  for (int k = 0; k < 3; k++)
    {
    flip += refDirection[2][k]*imgDirection[2][k];
    }

  if (flip > 0)
    {
    return;
    }

  // do ZFlip

  OrientedMaskImageType::PointType origin = image->GetOrigin();
  OrientedMaskImageType::SpacingType spacing = image->GetSpacing();

  OrientedMaskImageType::RegionType region = image->GetLargestPossibleRegion();
  OrientedMaskImageType::IndexType start = region.GetIndex();
  OrientedMaskImageType::SizeType size = region.GetSize();

  for (int k = 0; k < 3; k++)
    {
    start[k] += (size[k]-1);
    }

  OrientedMaskImageType::PointType pt;
  image->TransformIndexToPhysicalPoint( start, pt );

  // flip 
  int nZSlices = size[2];
  if ( size[2] % 2 == 0 )
    {
    size[2] = size[2]/2;
    }
  else
    {
    size[2] = (size[2]-1)/2;
    }
  region.SetSize( size );

  itk::ImageRegionIteratorWithIndex<OrientedMaskImageType> it( image, region );
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    OrientedMaskImageType::IndexType idx = it.GetIndex();
    OrientedMaskImageType::IndexType idxMirror = idx;
    idxMirror[2] = nZSlices-idx[2]-1;

    OrientedMaskImageType::PixelType a = image->GetPixel( idx );
    image->SetPixel( idx, image->GetPixel(idxMirror) );
    image->SetPixel( idxMirror, a );
    }

  OrientedMaskImageType::DirectionType direction = image->GetDirection();
  if (direction[2][2] < 0)
    {
    origin[2] = origin[2]-spacing[2]*(nZSlices-1);
    }
  else
    {
    origin[2] = origin[2]+spacing[2]*(nZSlices-1);
    }

  image->SetOrigin( origin );
  for (int k = 0; k < 3; k++)
    {
    direction[2][k] = -direction[2][k];
    }
  image->SetDirection( direction );

  return;

  }


 int ReadArchetypeImage(const std::string &inputfile,
                       OrientedInputImageType::Pointer &img,
                       std::string &subjectname)
{
  typedef itk::ImageFileReader< OrientedInputImageType >            ReaderType;
  typedef itk::GDCMImageIO                                  ImageIOType;
  typedef itk::GDCMSeriesFileNames                          SeriesFileNames;
  typedef itk::ImageSeriesReader< OrientedInputImageType >          SeriesReaderType;
  typedef itk::ArchetypeSeriesFileNames                     ArchetypeFileNames;
  
  if (itksys::SystemTools::FileIsDirectory(inputfile.c_str() ))    
    {
    // DICOM series input
    std::vector<std::string> lines;
    itksys::SystemTools::Split(inputfile.c_str(), lines, '/');
    subjectname = lines[lines.size()-1];
    
    ImageIOType::Pointer gdcmIO = ImageIOType::New();
    SeriesFileNames::Pointer it = SeriesFileNames::New();

    it->SetInputDirectory( inputfile.c_str() );

    SeriesReaderType::Pointer sReader = SeriesReaderType::New();

    const SeriesReaderType::FileNamesContainer & filenames =
      it->GetInputFileNames();

    unsigned int numberOfFilename = filenames.size();

    sReader->SetFileNames( filenames );
    sReader->SetImageIO( gdcmIO );

    try
      {
      sReader->Update();
      img = sReader->GetOutput();
      }
    catch (itk::ExceptionObject &excp)
      {
      std::cerr << "Exception thrown while reading the image" <<
        std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }

    return 1;
    }
  else
    {
    // Input is a file    
    itk::GDCMImageIO::Pointer dicomio = itk::GDCMImageIO::New();
    if ( dicomio->CanReadFile(inputfile.c_str()) ) 
      {
      // Input is a single DICOM file      
      std::string dicompathname
        = itksys::SystemTools::GetFilenamePath(inputfile);
      std::vector<std::string> lines;
      itksys::SystemTools::Split(dicompathname.c_str(), lines, '/');
      subjectname = lines[lines.size()-1];

      ImageIOType::Pointer gdcmIO = ImageIOType::New();
      SeriesFileNames::Pointer it = SeriesFileNames::New();
      
      it->SetInputDirectory( dicompathname.c_str() );
      
      SeriesReaderType::Pointer sReader = SeriesReaderType::New();
      
      const SeriesReaderType::FileNamesContainer & filenames =
        it->GetInputFileNames();
      
      unsigned int numberOfFilename = filenames.size();
      
      sReader->SetFileNames( filenames );
      sReader->SetImageIO( gdcmIO );
      
      try
        {
        sReader->Update();
        img = sReader->GetOutput();
        }
      catch (itk::ExceptionObject &excp)
        {
        std::cerr << "Exception thrown while writing the image" <<
          std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
        }
      
      return 1;
      }
    else
      {
      // Input is a single file but not DICOM
      subjectname
        = itksys::SystemTools::GetFilenameWithoutExtension(inputfile);

      ReaderType::Pointer reader2 = ReaderType::New();
      reader2->SetFileName( inputfile.c_str() );
      reader2->UpdateOutputInformation();
      OrientedInputImageType::RegionType region
        = reader2->GetOutput()->GetLargestPossibleRegion();

      if (region.GetSize()[2] > 1)
        {
        // Input is a single file representing a volume
        try
          {
          reader2->Update();
          img = reader2->GetOutput();
          }
        catch (itk::ExceptionObject &ex)
          {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
          }

        return 1;
        }
      else
        {
        // Input is file representing a slice, use it as an archetype
        // to collect all the slices.
        ArchetypeFileNames::Pointer it = ArchetypeFileNames::New();
        it->SetArchetype( inputfile.c_str() );

        const SeriesReaderType::FileNamesContainer &filenames =
          it->GetFileNames();

        // for (int i=0; i < it->GetFileNames().size(); ++i)
        //   std::cout << it->GetFileNames()[i] << std::endl;
        
        SeriesReaderType::Pointer sReader = SeriesReaderType::New();
        sReader->SetFileNames( filenames );

        try
          {
          sReader->Update();
          img = sReader->GetOutput();
          }
        catch (itk::ExceptionObject &excp)
          {
          std::cerr << "Exception thrown while reading the image" <<
            std::endl;
          std::cerr << excp << std::endl;
          return EXIT_FAILURE;
          }

        return 1;

        }
      }
    }
  
  return 1;
}

int WriteArchetypeImage(const std::string &file,
                        OrientedOutputImageType *img,
                        bool byslice, bool compress)
{
  typedef itk::ImageFileWriter<OrientedOutputImageType> WriterType;

  if (!byslice)
    {
    // write out the image using a file writer (single file)
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( file );
    writer->SetInput( img );

    if (compress)
      {
      writer->UseCompressionOn();
      }
    
    try
      {
      writer->Write();
      }
    catch (itk::ExceptionObject &exc)
      {
      std::cout << "Exception caught while writing image: " << exc
                << std::endl;
      return EXIT_FAILURE;
      }
    catch (std::exception &exc)
      {
      std::cout << "STD exception caught while writing image: " 
                << std::endl;
      return EXIT_FAILURE;
      }
    catch (...)
      {
      std::cout << "Unknown exception caughy while writing image."
                << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    // write out the image as if by a series writer.  we don't
    // actually use a series writer because it will not properly write
    // out the 3D meta information to a 2D image (except for DICOM).
    itk::NumericSeriesFileNames::Pointer outputNames
      = itk::NumericSeriesFileNames::New();
    outputNames->SetSeriesFormat( file );
    // could do these from StartIndex to StartIndex + Size - 1
    outputNames->SetStartIndex( 0 );
    outputNames->SetEndIndex(img->GetBufferedRegion().GetSize()[2]-1);

    WriterType::Pointer writer= WriterType::New();
    
    typedef itk::RegionOfInterestImageFilter<OrientedOutputImageType, OrientedOutputImageType>
      ExtractType;
    ExtractType::Pointer extract = ExtractType::New();
    extract->SetInput( img );
    
    ExtractType::RegionType region = img->GetBufferedRegion();
    ExtractType::IndexType tindex, originalindex=region.GetIndex();
    ExtractType::SizeType size = region.GetSize();

    // region will be a slice at a time
    size[2] = 1;
    region.SetSize( size );
    
    // loop over each slice and write a single image
    tindex = originalindex;
    for (int i=0; i < outputNames->GetFileNames().size(); ++i)
      {
      tindex[2] = originalindex[2] + i;
      region.SetIndex( tindex );
      extract->SetRegionOfInterest( region );
      extract->UpdateLargestPossibleRegion();
    
      writer->SetInput( extract->GetOutput() );
      writer->SetFileName( outputNames->GetFileNames()[i] );
      if (compress)
        {
        writer->UseCompressionOn();
        }

      try
        {
        writer->Write();
        }
      catch (itk::ExceptionObject &exc)
        {
        std::cout << "Exception caught while writing image "
                  << i << " of " << outputNames->GetFileNames().size()
                  << ": " << exc
                  << std::endl;
        return EXIT_FAILURE;
        }
      catch (std::exception &exc)
        {
        std::cout << "STD exception caught while writing image "
                  << i << " of " << outputNames->GetFileNames().size()
                  << std::endl;
        return EXIT_FAILURE;
        }
      catch (...)
        {
        std::cout << "Unknown exception caughy while writing image "
                  << i << " of " << outputNames->GetFileNames().size()
                  << std::endl;
        return EXIT_FAILURE;
        }
      }
    }

  return EXIT_SUCCESS;

}
