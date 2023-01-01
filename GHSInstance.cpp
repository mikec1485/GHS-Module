// ----------------------------------------------------------------------------
//      ______  __  __  ______
//     / ____/ / / / / / ____/
//    / / __  / /_/ / / /___
//   / / / / / __  / /___  /
//  / /_/ / / / / / ____/ /    GeneralisedHyperbolicStretch
//  \____/ /_/ /_/ /_____/     https://www.ghsastro.co.uk
// ----------------------------------------------------------------------------
// Standard GHS Process Module
// ----------------------------------------------------------------------------
// GHSInstance.cpp
// ----------------------------------------------------------------------------
//
// Copyright (C) 2022,2023 Mike Cranfield
//
// This product is based on software from the PixInsight project, developed
// by Pleiades Astrophoto and its contributors (https://pixinsight.com/).
// ----------------------------------------------------------------------------
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, version 3 of the License.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------------

#include "GHSInstance.h"

#include <pcl/Histogram.h>

#include <pcl/AutoViewLock.h>
#include <pcl/StandardStatus.h>
#include <pcl/View.h>

namespace pcl
{

// ----------------------------------------------------------------------------

GHSInstance::GHSInstance( const MetaProcess* m )
   : ProcessImplementation( m )
, p_ST( TheGHSSTParameter->DefaultValueIndex() )
, p_SC( TheGHSSCParameter->DefaultValueIndex() )
, p_Inv( TheGHSInvParameter->DefaultValue() )
, p_D( TheGHSParameter->DefaultValue() )
, p_b( TheGHSbParameter->DefaultValue() )
, p_SP( TheGHSSPParameter->DefaultValue() )
, p_LP( TheGHSLPParameter->DefaultValue() )
, p_HP( TheGHSHPParameter->DefaultValue() )
, p_BP( TheGHSBPParameter->DefaultValue() )
, p_WP( TheGHSWPParameter->DefaultValue() )
, p_CB( TheGHSCBParameter->DefaultValue() )
, p_CT( TheGHSCTParameter->DefaultValueIndex() )
, p_RGBWS( TheGHSRGBWSParameter->DefaultValue() )
{
    
}

// ----------------------------------------------------------------------------

GHSInstance::GHSInstance( const GHSInstance& x )
   : ProcessImplementation( x )
{
   Assign( x );
}

// ----------------------------------------------------------------------------

void GHSInstance::Assign( const ProcessImplementation& p )
{
   const GHSInstance* x = dynamic_cast<const GHSInstance*>( &p );
   if ( x != nullptr )
   {
       SetStretchType(x->p_ST, false);
       SetStretchChannel(x->p_SC, false);
       SetInverse(x->p_Inv, false);
       SetStretchFactor(x->p_D, false);
       SetLocalIntensity(x->p_b, false);
       SetSymmetryPoint(x->p_SP, false);
       SetShadowProtect(x->p_LP, false);
       SetHighlightProtect(x->p_HP, false);
       SetBlackpoint(x->p_BP, false);
       SetWhitepoint(x->p_WP, false);
       SetColourBlend(x->p_CB, false);
       SetClipType(x->p_CT, false);
       SetRGBWS(x->p_RGBWS, false);
       
       UpdateFlags();
       UpdateCoeffs();
   }
}

// ----------------------------------------------------------------------------

UndoFlags GHSInstance::UndoMode( const View& ) const
{
   /*
    * The following flag assumes that your process modifies pixel sample values
    * *exclusively*. If you are going to change anything else, such as image
    * geometry, keywords, etc., or maybe nothing at all (in case you are
    * implementing an image observer process), see the UpdateFlag enumeration
    * in pcl/ImageWindow.h for complete information.
    */
   return UndoFlag::PixelData;
}

// ----------------------------------------------------------------------------

bool GHSInstance::CanExecuteOn( const View& view, String& whyNot ) const
{
   if ( view.Image().IsComplexSample() )
   {
      whyNot = "GHS cannot be executed on complex images.";
      return false;
   }
   return true;
}

// ----------------------------------------------------------------------------

bool GHSInstance::CanExecuteOn( const ImageVariant& image, String& whyNot ) const
{
   if ( image.IsComplexSample() )
   {
      whyNot = "GHS cannot be executed on complex images.";
      return false;
   }
   return true;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

class GHSEngine
{
public:
    
    template <class P> static
    void Apply( GenericImage<P>& image, const GHSInstance& G, bool preview = false)
    {
        if ( image.IsEmptySelection() )
            return;
        
        image.EnsureUnique();

        Rect r = image.SelectedRectangle();
        int h = r.Height();
        
        Array<size_type> L = pcl::Thread::OptimalThreadLoads( h );
        size_type N = image.NumberOfSelectedSamples();
        if ( image.Status().IsInitializationEnabled() )
            image.Status().Initialize( "GHS transformation", N );
        
        ThreadData<P> data( image, G, N );

        ReferenceArray<GHSThread<P> > threads;
        for ( int i = 0, n = 0; i < int( L.Length() ); n += int( L[i++] ) )
            threads.Add( new GHSThread<P>( data, n, n + int( L[i] ) ) );
        
        AbstractImage::RunThreads( threads, data );
        
        if (G.p_CT == GHSCT::CT_RescaleGlobal)
        {
            double overallMax = 1.0;
            for ( int i = 0; i < int( L.Length() ); ++i)
                overallMax = Max( overallMax, threads[i].MaxAdjustment);
            for ( int i = 0; i < int( L.Length() ); ++i)
            {
                threads[i].MaxAdjustment = overallMax;
                threads[i].RescaleGlobalFirstPass = false;
            }
            AbstractImage::RunThreads( threads, data );
        }
        
        threads.Destroy();
        
        image.Status() = data.status;
    }
    
    // ----------------------------------------------------------------------------

    static void Apply( Histogram& dstH, const Histogram& srcH, const GHSInstance& G )
    {
        
        if ( srcH.IsEmpty() )
          return;
        
        Histogram::histogram_type histData(0, srcH.Resolution());
        
        double last = 0.0;
        G.Transform( last );

       for ( int i = 0; i < srcH.Resolution(); ++i )
       {
           double next = double( i + 1 )/srcH.Resolution();
           G.Transform( next );
           double span = next - last;
           
           if (next == last)
           {
               double start = last * dstH.Resolution();
               int startInt = Floor(start);
               if (startInt < srcH.Resolution())
                   histData[startInt] += uint64(srcH.Count(i));
           }
           else
           {
               double quota = srcH.Count(i) / (span * dstH.Resolution());
               double start = last * dstH.Resolution();
               int startInt = Floor(start);
               double end = next * dstH.Resolution();
               int endInt = Ceil(end);
               
               for (int i = startInt; i < Min(srcH.Resolution(), endInt); i++)
                   histData[i] += uint64(quota);
               
               histData[startInt] -= uint64((start - double(startInt)) * quota);
               if ( !(endInt > srcH.Resolution()) )
                   histData[endInt - 1] -= uint64((double(endInt) - end) * quota);
           }
           
           last = next;
       }
        
        dstH.SetHistogramData( histData );
        histData.Clear();
    }

private:

   template <class P>
   struct ThreadData : public AbstractImage::ThreadData
   {
      ThreadData( GenericImage<P>& a_image, const GHSInstance& a_instance, size_type a_count )
         : AbstractImage::ThreadData( a_image, a_count )
         , image( a_image )
         , instance( a_instance )
      {
      }

      GenericImage<P>&      image;
      const GHSInstance&    instance;
   };

   template <class P>
   class GHSThread : public pcl::Thread
   {
   public:

      GHSThread( ThreadData<P>& d, int startRow, int endRow )
         : m_data( d )
         , m_firstRow( startRow )
         , m_endRow( endRow ) // N.B. m_firstRow, m_endRow are relative to the current image selection
      {
      }
       
      double MaxAdjustment = 1.0;
      bool RescaleGlobalFirstPass = true;

      void Run() override
      {
          INIT_THREAD_MONITOR()

          Rect r = m_data.image.SelectedRectangle();
          size_type w = r.Width();
          int numChannels = m_data.image.NumberOfNominalChannels();
          RGBColorSystem rgbws = m_data.image.RGBWorkingSpace();
          
          pcl_enum m_stretchChannel = m_data.instance.p_SC;
          pcl_enum m_clipType = m_data.instance.p_CT;
          
          double m_LCR = 0.0;
          double m_LCG = 0.0;
          double m_LCB = 0.0;
          if (m_data.instance.p_RGBWS)
          {
              F32Vector lumCoeffs = rgbws.LuminanceCoefficients();
              m_LCR = lumCoeffs[0];
              m_LCG = lumCoeffs[1];
              m_LCB = lumCoeffs[2];
          }
          double LCTot = m_LCR + m_LCG + m_LCB;
          m_LCR = (LCTot > 0) ? m_LCR /= LCTot : 1.0/3.0;
          m_LCG = (LCTot > 0) ? m_LCG /= LCTot : 1.0/3.0;
          m_LCB = (LCTot > 0) ? m_LCB /= LCTot : 1.0/3.0;
          
          double m_CB = m_data.instance.p_CB;
             
          if (numChannels == 1)
          {
              for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
              {
                  typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                  for ( size_type i = 0; i < w; ++i )
                  {
                      double f0; P::FromSample( f0, *(p0 + i) );
                      m_data.instance.Transform( f0 );
                      *(p0 + i) = P::ToSample( f0 );
                      UPDATE_THREAD_MONITOR( 65536 )
                  }
              }
          }
          else
          {
             switch (m_stretchChannel)
             {
                 case GHSSC::SC_Red:    // Red stretch
                     for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                     {
                         typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                         for ( size_type i = 0; i < w; ++i )
                         {
                             double f0; P::FromSample( f0, *(p0 + i) );
                             m_data.instance.Transform( f0 );
                             *(p0 + i) = P::ToSample( f0 );
                             UPDATE_THREAD_MONITOR( 65536 )
                         }
                     }
                     break;
                 case GHSSC::SC_Green:    // Green stretch
                     for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                     {
                         typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                         for ( size_type i = 0; i < w; ++i )
                         {
                             double f1; P::FromSample( f1, *(p1 + i) );
                             m_data.instance.Transform( f1 );
                             *(p1 + i) = P::ToSample( f1 );
                             UPDATE_THREAD_MONITOR( 65536 )
                         }
                     }
                     break;
                 case GHSSC::SC_Blue:    // Blue stretch
                     for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                     {
                         typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                         for ( size_type i = 0; i < w; ++i )
                         {
                             double f2; P::FromSample( f2, *(p2 + i) );
                             m_data.instance.Transform( f2 );
                             *(p2 + i) = P::ToSample( f2 );
                             UPDATE_THREAD_MONITOR( 65536 )
                         }
                     }
                     break;
                 case GHSSC::SC_RGB:    // RGB stretch
                     for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                     {
                         typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                         typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                         typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                         for ( size_type i = 0; i < w; ++i )
                         {
                             double f0; P::FromSample( f0, *(p0 + i) );
                             double f1; P::FromSample( f1, *(p1 + i) );
                             double f2; P::FromSample( f2, *(p2 + i) );
                             m_data.instance.Transform( f0 );
                             m_data.instance.Transform( f1 );
                             m_data.instance.Transform( f2 );
                             *(p0 + i) = P::ToSample( f0 );
                             *(p1 + i) = P::ToSample( f1 );
                             *(p2 + i) = P::ToSample( f2 );
                             UPDATE_THREAD_MONITOR( 65536 )
                         }
                     }
                     break;
                 case GHSSC::SC_Lightness:    // Lightness stretch
                     for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                     {
                         typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                         typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                         typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                         for ( size_type i = 0; i < w; ++i )
                         {
                             double f0; P::FromSample( f0, *(p0 + i) );
                             double f1; P::FromSample( f1, *(p1 + i) );
                             double f2; P::FromSample( f2, *(p2 + i) );
                             double cieL = 0.0, cieA = 0.0, cieB = 0.0;
                             rgbws.RGBToCIELab(cieL, cieA, cieB, f0, f1, f2);
                             m_data.instance.Transform( cieL );
                             rgbws.CIELabToRGB(f0, f1, f2, cieL, cieA, cieB);
                             *(p0 + i) = P::ToSample( f0 );
                             *(p1 + i) = P::ToSample( f1 );
                             *(p2 + i) = P::ToSample( f2 );
                             UPDATE_THREAD_MONITOR( 65536 )
                         }
                     }
                     break;
                 case GHSSC::SC_Saturation:    // Saturation stretch
                     for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                     {
                         typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                         typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                         typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                         for ( size_type i = 0; i < w; ++i )
                         {
                             double f0; P::FromSample( f0, *(p0 + i) );
                             double f1; P::FromSample( f1, *(p1 + i) );
                             double f2; P::FromSample( f2, *(p2 + i) );
                             double cieL = 0.0, hue = 0.0, sat = 0.0, val = 0.0;
                             rgbws.RGBToHSVL(hue, sat, val, cieL, f0, f1, f2);
                             m_data.instance.Transform( sat );
                             rgbws.HSVLToRGB(f0, f1, f2, hue, sat, val, cieL);
                             *(p0 + i) = P::ToSample( f0 );
                             *(p1 + i) = P::ToSample( f1 );
                             *(p2 + i) = P::ToSample( f2 );
                             UPDATE_THREAD_MONITOR( 65536 )
                         }
                     }
                     break;
                 case GHSSC::SC_Colour:    // Colour stretch
                     switch (m_clipType)
                     {
                         case GHSCT::CT_Clip:
                            if (m_CB == 1.0)
                            {
                                for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                                {
                                    typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                                    typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                                    typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                                    for ( size_type i = 0; i < w; ++i )
                                    {
                                        double f0; P::FromSample( f0, *(p0 + i) );
                                        double f1; P::FromSample( f1, *(p1 + i) );
                                        double f2; P::FromSample( f2, *(p2 + i) );
                                             
                                        double fbar0 = (m_LCR * f0 + m_LCG * f1 +  m_LCB * f2);
                                        double fbar1 = fbar0;
                                        m_data.instance.Transform( fbar1 );
                                        double stretchFactor = (fbar0 == 0) ? 0 : fbar1 / fbar0;
                                        f0 = stretchFactor * f0;
                                        f1 = stretchFactor * f1;
                                        f2 = stretchFactor * f2;
                                        double fmax = Max(Max(f0, f1), f2);
                                        if (fmax > 1.0)
                                        {
                                            f0 = pcl::Range(f0, 0.0, 1.0);
                                            f1 = pcl::Range(f1, 0.0, 1.0);
                                            f2 = pcl::Range(f2, 0.0, 1.0);
                                        }
                                        *(p0 + i) = P::ToSample( f0 );
                                        *(p1 + i) = P::ToSample( f1 );
                                        *(p2 + i) = P::ToSample( f2 );
                                        UPDATE_THREAD_MONITOR( 65536 )
                                    }
                                }
                            }
                            else
                            {
                                for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                                {
                                    typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                                    typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                                    typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                                    for ( size_type i = 0; i < w; ++i )
                                    {
                                        double f0; P::FromSample( f0, *(p0 + i) );
                                        double f1; P::FromSample( f1, *(p1 + i) );
                                        double f2; P::FromSample( f2, *(p2 + i) );
                                        
                                        double sf0 = f0;
                                        double sf1 = f1;
                                        double sf2 = f2;
                                             
                                        double fbar = (m_LCR * f0 + m_LCG * f1 +  m_LCB * f2);
                                        double sfbar = fbar;
                                        
                                        m_data.instance.Transform( sf0 );
                                        m_data.instance.Transform( sf1 );
                                        m_data.instance.Transform( sf2 );
                                        m_data.instance.Transform( sfbar );
                                        
                                        double stretchFactor0 = (f0 == 0) ? 0 : sf0 / f0;
                                        double stretchFactor1 = (f1 == 0) ? 0 : sf1 / f1;
                                        double stretchFactor2 = (f2 == 0) ? 0 : sf2 / f2;
                                        double stretchFactor = (fbar == 0) ? 0 : sfbar / fbar;
                                        f0 *= (m_CB * stretchFactor + (1 - m_CB) * stretchFactor0);
                                        f1 *= (m_CB * stretchFactor + (1 - m_CB) * stretchFactor1);
                                        f2 *= (m_CB * stretchFactor + (1 - m_CB) * stretchFactor2);
                                        
                                        double fmax = Max(Max(f0, f1), f2);
                                        if (fmax > 1.0)
                                        {
                                            f0 = pcl::Range(f0, 0.0, 1.0);
                                            f1 = pcl::Range(f1, 0.0, 1.0);
                                            f2 = pcl::Range(f2, 0.0, 1.0);
                                        }
                                        
                                        *(p0 + i) = P::ToSample( f0 );
                                        *(p1 + i) = P::ToSample( f1 );
                                        *(p2 + i) = P::ToSample( f2 );
                                        UPDATE_THREAD_MONITOR( 65536 )
                                    }
                                }
                            }
                            break;
                         case GHSCT::CT_Rescale:
                             if (m_CB == 1.0)
                             {
                                 for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                                 {
                                     typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                                     typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                                     typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                                     for ( size_type i = 0; i < w; ++i )
                                     {
                                         double f0; P::FromSample( f0, *(p0 + i) );
                                         double f1; P::FromSample( f1, *(p1 + i) );
                                         double f2; P::FromSample( f2, *(p2 + i) );
                                       
                                         double fbar0 = (m_LCR * f0 + m_LCG * f1 +  m_LCB * f2);
                                         double fbar1 = fbar0;
                                         m_data.instance.Transform( fbar1 );
                                         double stretchFactor = (fbar0 == 0) ? 0 : fbar1 / fbar0;
                                         f0 *= stretchFactor;
                                         f1 *= stretchFactor;
                                         f2 *= stretchFactor;
                                         double fmax = Max(Max(f0, f1), f2);
                                         if (fmax > 1.0)
                                         {
                                             f0 /= fmax;
                                             f1 /= fmax;
                                             f2 /= fmax;
                                         }
                                         *(p0 + i) = P::ToSample( f0 );
                                         *(p1 + i) = P::ToSample( f1 );
                                         *(p2 + i) = P::ToSample( f2 );
                                         UPDATE_THREAD_MONITOR( 65536 )
                                     }
                                 }
                             }
                             else
                             {
                                 for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                                 {
                                     typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                                     typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                                     typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                                     for ( size_type i = 0; i < w; ++i )
                                     {
                                         double f0; P::FromSample( f0, *(p0 + i) );
                                         double f1; P::FromSample( f1, *(p1 + i) );
                                         double f2; P::FromSample( f2, *(p2 + i) );
                                         
                                         double sf0 = f0;
                                         double sf1 = f1;
                                         double sf2 = f2;
                                              
                                         double fbar = (m_LCR * f0 + m_LCG * f1 +  m_LCB * f2);
                                         double sfbar = fbar;
                                         
                                         m_data.instance.Transform( sf0 );
                                         m_data.instance.Transform( sf1 );
                                         m_data.instance.Transform( sf2 );
                                         m_data.instance.Transform( sfbar );
                                         
                                         double stretchFactor0 = (f0 == 0) ? 0 : sf0 / f0;
                                         double stretchFactor1 = (f1 == 0) ? 0 : sf1 / f1;
                                         double stretchFactor2 = (f2 == 0) ? 0 : sf2 / f2;
                                         double stretchFactor = (fbar == 0) ? 0 : sfbar / fbar;
                                         f0 *= (m_CB * stretchFactor + (1 - m_CB) * stretchFactor0);
                                         f1 *= (m_CB * stretchFactor + (1 - m_CB) * stretchFactor1);
                                         f2 *= (m_CB * stretchFactor + (1 - m_CB) * stretchFactor2);
                                         
                                         double fmax = Max(Max(f0, f1), f2);
                                         if (fmax > 1.0)
                                         {
                                             f0 /= fmax;
                                             f1 /= fmax;
                                             f2 /= fmax;
                                         }
                                         
                                         *(p0 + i) = P::ToSample( f0 );
                                         *(p1 + i) = P::ToSample( f1 );
                                         *(p2 + i) = P::ToSample( f2 );
                                         UPDATE_THREAD_MONITOR( 65536 )
                                     }
                                 }
                             }
                            break;
                         case GHSCT::CT_RescaleGlobal:
                             if (m_CB == 1.0)
                             {
                                 for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                                 {
                                     typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                                     typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                                     typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                                     for ( size_type i = 0; i < w; ++i )
                                     {
                                         double f0; P::FromSample( f0, *(p0 + i) );
                                         double f1; P::FromSample( f1, *(p1 + i) );
                                         double f2; P::FromSample( f2, *(p2 + i) );
                                       
                                         double fbar0 = (m_LCR * f0 + m_LCG * f1 +  m_LCB * f2);
                                         double fbar1 = fbar0;
                                         m_data.instance.Transform( fbar1 );
                                         double stretchFactor = (fbar0 == 0) ? 0 : fbar1 / fbar0;
                                         f0 *= stretchFactor;
                                         f1 *= stretchFactor;
                                         f2 *= stretchFactor;
                                         
                                         if (RescaleGlobalFirstPass)
                                             MaxAdjustment = Max(Max(Max(f0, f1), f2), MaxAdjustment);
                                         else
                                         {
                                             f0 /= MaxAdjustment;
                                             f1 /= MaxAdjustment;
                                             f2 /= MaxAdjustment;
                                             *(p0 + i) = P::ToSample( f0 );
                                             *(p1 + i) = P::ToSample( f1 );
                                             *(p2 + i) = P::ToSample( f2 );
                                         }
                                         UPDATE_THREAD_MONITOR( 65536 )
                                     }
                                 }
                             }
                             else
                             {
                                 for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                                 {
                                     typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                                     typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                                     typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                                     for ( size_type i = 0; i < w; ++i )
                                     {
                                         double f0; P::FromSample( f0, *(p0 + i) );
                                         double f1; P::FromSample( f1, *(p1 + i) );
                                         double f2; P::FromSample( f2, *(p2 + i) );
                                         
                                         double sf0 = f0;
                                         double sf1 = f1;
                                         double sf2 = f2;
                                              
                                         double fbar = (m_LCR * f0 + m_LCG * f1 +  m_LCB * f2);
                                         double sfbar = fbar;
                                         
                                         m_data.instance.Transform( sf0 );
                                         m_data.instance.Transform( sf1 );
                                         m_data.instance.Transform( sf2 );
                                         m_data.instance.Transform( sfbar );
                                         
                                         double stretchFactor0 = (f0 == 0) ? 0 : sf0 / f0;
                                         double stretchFactor1 = (f1 == 0) ? 0 : sf1 / f1;
                                         double stretchFactor2 = (f2 == 0) ? 0 : sf2 / f2;
                                         double stretchFactor = (fbar == 0) ? 0 : sfbar / fbar;
                                         f0 *= (m_CB * stretchFactor + (1 - m_CB) * stretchFactor0);
                                         f1 *= (m_CB * stretchFactor + (1 - m_CB) * stretchFactor1);
                                         f2 *= (m_CB * stretchFactor + (1 - m_CB) * stretchFactor2);
                                         
                                         if (RescaleGlobalFirstPass)
                                             MaxAdjustment = Max(Max(Max(f0, f1), f2), MaxAdjustment);
                                         else
                                         {
                                             f0 /= MaxAdjustment;
                                             f1 /= MaxAdjustment;
                                             f2 /= MaxAdjustment;
                                             *(p0 + i) = P::ToSample( f0 );
                                             *(p1 + i) = P::ToSample( f1 );
                                             *(p2 + i) = P::ToSample( f2 );
                                         }
                                         UPDATE_THREAD_MONITOR( 65536 )
                                     }
                                 }
                             }
                            break;
                         case GHSCT::CT_RGBBlend:
                             for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                            {
                                typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                                typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                                typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                                for ( size_type i = 0; i < w; ++i )
                                {
                                    double f0; P::FromSample( f0, *(p0 + i) );
                                    double f1; P::FromSample( f1, *(p1 + i) );
                                    double f2; P::FromSample( f2, *(p2 + i) );
                                    double sf0 = f0;
                                    double sf1 = f1;
                                    double sf2 = f2;
                           
                                    double fbar = (m_LCR * f0 + m_LCG * f1 +  m_LCB * f2);
                                    double sfbar = fbar;
                                    m_data.instance.Transform( sfbar );
                                    double stretchFactor = (fbar == 0) ? 0 : sfbar / fbar;
                                    sf0 *= stretchFactor;
                                    sf1 *= stretchFactor;
                                    sf2 *= stretchFactor;
                                    double sfmax = Max(Max(sf0, sf1), sf2);
                                
                                    m_data.instance.Transform( f0 );
                                    m_data.instance.Transform( f1 );
                                    m_data.instance.Transform( f2 );
                                    double tfmax = Max(Max(f0, f1), f2);
                                
                                    double d = sfmax - tfmax;
                                    if (tfmax + m_CB * d > 1)
                                    {
                                        double k = (d != 0) ? Min(m_CB, (1.0 - tfmax)/d) : m_CB;
                                        f0 = (1 - k) * f0 + k * sf0;
                                        f1 = (1 - k) * f1 + k * sf1;
                                        f2 = (1 - k) * f2 + k * sf2;
                                    } else {
                                        f0 = (1 - m_CB) * f0 + m_CB * sf0;
                                        f1 = (1 - m_CB) * f1 + m_CB * sf1;
                                        f2 = (1 - m_CB) * f2 + m_CB * sf2;
                                    }
                                    *(p0 + i) = P::ToSample( f0 );
                                    *(p1 + i) = P::ToSample( f1 );
                                    *(p2 + i) = P::ToSample( f2 );
                                    UPDATE_THREAD_MONITOR( 65536 )
                                }
                            }
                            break;
                    }
                    break;
                 default:
                     for ( int y = r.y0+m_firstRow, y1 = r.y0+m_endRow; y < y1; ++y )
                     {
                         typename P::sample* p0 = m_data.image.PixelAddress( r.x0, y, 0 );
                         typename P::sample* p1 = m_data.image.PixelAddress( r.x0, y, 1 );
                         typename P::sample* p2 = m_data.image.PixelAddress( r.x0, y, 2 );
                         for ( size_type i = 0; i < w; ++i )
                         {
                             double f0; P::FromSample( f0, *(p0 + i) );
                             double f1; P::FromSample( f1, *(p1 + i) );
                             double f2; P::FromSample( f2, *(p2 + i) );
                             m_data.instance.Transform( f0 );
                             m_data.instance.Transform( f1 );
                             m_data.instance.Transform( f2 );
                             *(p0 + i) = P::ToSample( f0 );
                             *(p1 + i) = P::ToSample( f1 );
                             *(p2 + i) = P::ToSample( f2 );
                             UPDATE_THREAD_MONITOR( 65536 )
                         }
                     }
             }
          }
      };

   private:

      ThreadData<P>& m_data;
      int            m_firstRow;
      int            m_endRow;
   };
};

// ----------------------------------------------------------------------------

bool GHSInstance::ExecuteOn( View& view )
{
    AutoViewLock lock( view );

    ImageVariant image = view.Image();
    return ExecuteOn(image, "");
}

// ----------------------------------------------------------------------------

bool GHSInstance::ExecuteOn( ImageVariant& image, const IsoString& hints )
{
   if ( image.IsComplexSample() )
      return false;

   StandardStatus status;
   image.SetStatusCallback( &status );

   Console().EnableAbort();
    
    this->UpdateFlags();
    this->UpdateCoeffs();

   if ( image.IsFloatSample() )
      switch ( image.BitsPerSample() )
      {
      case 32: GHSEngine::Apply( static_cast<Image&>( *image ), *this ); break;
      case 64: GHSEngine::Apply( static_cast<DImage&>( *image ), *this ); break;
      }
   else
      switch ( image.BitsPerSample() )
      {
      case  8: GHSEngine::Apply( static_cast<UInt8Image&>( *image ), *this ); break;
      case 16: GHSEngine::Apply( static_cast<UInt16Image&>( *image ), *this ); break;
      case 32: GHSEngine::Apply( static_cast<UInt32Image&>( *image ), *this ); break;
      }

   return true;
}

// ----------------------------------------------------------------------------

void GHSInstance::TransformHistogram( Histogram& dstH, const Histogram& srcH )
{
    this->UpdateFlags();
    this->UpdateCoeffs();
    GHSEngine::Apply( dstH, srcH, *this );
}

// ----------------------------------------------------------------------------

void GHSInstance::Preview( UInt16Image& image )
{
   try
   {
       this->UpdateFlags();
       this->UpdateCoeffs();
       GHSEngine::Apply( image, *this, true);
       
   }
   catch ( ... )
   {
   }
}

// ----------------------------------------------------------------------------

void* GHSInstance::LockParameter( const MetaParameter* p, size_type /*tableRow*/ )
{
    if ( p == TheGHSSTParameter )
        return &p_ST;
    if ( p == TheGHSSCParameter )
        return &p_SC;
    if ( p == TheGHSInvParameter )
        return &p_Inv;
    if ( p == TheGHSParameter )
        return &p_D;
    if ( p == TheGHSbParameter )
        return &p_b;
    if ( p == TheGHSSPParameter )
        return &p_SP;
    if ( p == TheGHSLPParameter )
        return &p_LP;
    if ( p == TheGHSHPParameter )
        return &p_HP;
    if ( p == TheGHSBPParameter )
        return &p_BP;
    if ( p == TheGHSWPParameter )
       return &p_WP;
    if ( p == TheGHSCBParameter )
       return &p_CB;
    if ( p == TheGHSCTParameter )
       return &p_CT;
    if ( p == TheGHSRGBWSParameter )
        return &p_RGBWS;

   return nullptr;
}

// ----------------------------------------------------------------------------

} // pcl

// ----------------------------------------------------------------------------
// EOF GHSInstance.cpp
