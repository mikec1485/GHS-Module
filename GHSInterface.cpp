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
// GHSInterface.cpp
// ----------------------------------------------------------------------------
//
// Copyright (C) 2022, Mike Cranfield
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

#include "GHSInterface.h"
#include "GHSParameters.h"
#include "GHSProcess.h"

#include <pcl/AutoViewLock.h>
#include <pcl/StandardStatus.h>
#include <pcl/MuteStatus.h>
#include <pcl/RealTimePreview.h>
#include <pcl/Vector.h>
#include <pcl/ColorSpace.h>

#include <pcl/GlobalSettings.h>
#include <pcl/Graphics.h>
#include <pcl/HistogramTransformation.h>
#include <pcl/ImageWindow.h>
#include <pcl/Settings.h>

#define WHEEL_STEP_ANGLE   PixInsightSettings::GlobalInteger( "ImageWindow/WheelStepAngle" )


namespace pcl
{

// ----------------------------------------------------------------------------

GHSInterface* TheGHSInterface = nullptr;

// ----------------------------------------------------------------------------

static const int s_maxZoom = 999;

// ----------------------------------------------------------------------------

GHSInterface::GHSInterface()
   : m_instance( TheGHSProcess )
{
    TheGHSInterface = this;
}

// ----------------------------------------------------------------------------

GHSInterface::~GHSInterface()
{
   if ( GUI != nullptr )
      delete GUI, GUI = nullptr;
}

// ----------------------------------------------------------------------------

InterfaceFeatures GHSInterface::Features() const
{
   return InterfaceFeature::Default | InterfaceFeature::RealTimeButton;
}

// ----------------------------------------------------------------------------

IsoString GHSInterface::Id() const
{
   return "GeneralizedHyperbolicStretch";
}

// ----------------------------------------------------------------------------

MetaProcess* GHSInterface::Process() const
{
   return TheGHSProcess;
}

// ----------------------------------------------------------------------------

String GHSInterface::IconImageSVGFile() const
{
   return "@module_icons_dir/GHS.svg";
}

// ----------------------------------------------------------------------------

void GHSInterface::ApplyInstance() const
{
   m_instance.LaunchOnCurrentView();
}

// ----------------------------------------------------------------------------

void GHSInterface::RealTimePreviewUpdated(bool active)
{
   if ( GUI != nullptr )
      if ( active )
         RealTimePreview::SetOwner( *this ); // implicitly updates the real time preview
      else
         RealTimePreview::SetOwner( ProcessInterface::Null() );
}

// ----------------------------------------------------------------------------

void GHSInterface::ResetInstance()
{
    
    GHSInstance defaultInstance( TheGHSProcess );
    ImportProcess( defaultInstance );
    
    m_readoutSource = FromNone;
    m_linearFAParameter = BP;
    m_nonlinFAParameter = SP;
    m_channel = 3;
    
    CalculateOutputHistograms();
    CalculateClippingCounts();
    UpdateControls();
    UpdateRealTimePreview();
    
}

// ----------------------------------------------------------------------------

bool GHSInterface::Launch( const MetaProcess& P, const ProcessImplementation*, bool& dynamic, unsigned& /*flags*/ )
{
   if ( GUI == nullptr )
   {
       GUI = new GUIData( *this );
       SetWindowTitle( "GeneralizedHyperbolicStretch" );
       
       ImageWindow theActiveWindow( ImageWindow::ActiveWindow() );
       
       if ( theActiveWindow != ImageWindow::Null() )
       {
           m_currentView = theActiveWindow.CurrentView();
           SynchronizeWithCurrentView();
       }
       
       ReadoutOptions roo = ReadoutOptions::GetCurrentOptions();
       m_readoutSize = roo.ProbeSize();
       m_readoutStatistic = roo.Mode();
       
       UpdateControls();
   }

   dynamic = false;
   return &P == TheGHSProcess;
}

// ----------------------------------------------------------------------------

ProcessImplementation* GHSInterface::NewProcess() const
{
   return new GHSInstance( m_instance );
}

// ----------------------------------------------------------------------------

bool GHSInterface::ValidateProcess( const ProcessImplementation& p, String& whyNot ) const
{
   if ( dynamic_cast<const GHSInstance*>( &p ) != nullptr )
      return true;
   whyNot = "Not a GHS instance.";
   return false;
}

// ----------------------------------------------------------------------------

bool GHSInterface::RequiresInstanceValidation() const
{
   return true;
}

// ----------------------------------------------------------------------------

bool GHSInterface::ImportProcess( const ProcessImplementation& p )
{
    m_instance.Assign( p );
    
    CalculateOutputHistograms();
    CalculateClippingCounts();

    UpdateControls();
    UpdateHistograms();
    UpdateRealTimePreview();
    
    return true;
}

// ----------------------------------------------------------------------------

bool GHSInterface::RequiresRealTimePreviewUpdate( const UInt16Image&, const View&, const Rect&, int ) const
{
   return true;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

GHSInterface::RealTimeThread::RealTimeThread()
    : m_instance( TheGHSProcess )
{
}

// ----------------------------------------------------------------------------

void GHSInterface::RealTimeThread::Reset( const UInt16Image& image, const GHSInstance& instance )
{
    image.ResetSelections();
    m_image.Assign( image );
    m_instance.Assign( instance );
}

// ----------------------------------------------------------------------------

void GHSInterface::RealTimeThread::Run()
{
    //MuteStatus status;
    //m_image.SetStatusCallback( &status );
    //m_instance.UpdateFlags();   Not needed because Reset()  includes m_instance.Assign() which
    //m_instance.UpdateCoeffs();  updates flags and coefficients
    m_instance.Preview( m_image );
    //m_image.ResetSelections();
}

// ----------------------------------------------------------------------------

bool GHSInterface::GenerateRealTimePreview( UInt16Image& image, const View& view, const Rect& rect, int z, String& info) const
{
    m_realTimeThread = new RealTimeThread();

    for ( ;; )
    {
        m_realTimeThread->Reset( image, m_instance );
        m_realTimeThread->Start();

        while ( m_realTimeThread->IsActive() )
        {
            ProcessEvents();

            if ( !IsRealTimePreviewActive() )
            {
                m_realTimeThread->Abort();
                m_realTimeThread->Wait();
                return false;
            }
        }

        if ( !m_realTimeThread->IsAborted() )
        {
            image.Assign( m_realTimeThread->m_image );
            return true;
        }
    }
}

// ----------------------------------------------------------------------------


void GHSInterface::UpdateControls()
{
    UpdateZoomControls();
    UpdateGraphicsControls();
    UpdateReadoutControls();
    UpdateClippingCountControls();
    UpdateTransformationControls();
    UpdateHistograms();
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateRealTimePreview()
{
   if ( IsRealTimePreviewActive() )
   {
      if ( m_realTimeThread.IsValid() )
         m_realTimeThread->Abort();
      GUI->UpdateRealTimePreview_Timer.Start();
   }
}

// ----------------------------------------------------------------------------

bool GHSInterface::WantsImageNotifications() const
{
   return true;
}

// ----------------------------------------------------------------------------

void GHSInterface::ImageUpdated( const View& v )
{
   if ( GUI != nullptr )
      if ( v == m_currentView )
          SynchronizeWithCurrentView();
}

// ----------------------------------------------------------------------------

void GHSInterface::ImageFocused( const View& v )
{
   if ( GUI != nullptr )
   {
       if (m_currentView != v )
       {
           m_readoutSource = FromNone;
           m_currentView = v;
           SynchronizeWithCurrentView();
           EnsureLayoutUpdated();
           AdjustToContents();
       }
    }
}

// ----------------------------------------------------------------------------

void GHSInterface::ImageDeleted( const View& v )
{
   if ( GUI != nullptr )
   {
       if (m_currentView == v )
       {
           m_readoutSource = FromNone;
           m_currentView = View::Null();
           SynchronizeWithCurrentView();
           EnsureLayoutUpdated();
           AdjustToContents();
       }
    }
}

// ----------------------------------------------------------------------------

bool GHSInterface::WantsViewPropertyNotifications() const
{
   return true;
}

// ----------------------------------------------------------------------------

void GHSInterface::ViewPropertyUpdated( const View& v, const IsoString& property )
{
   if ( GUI != nullptr )
      if ( v == m_currentView )
         SynchronizeWithCurrentView();
}

// ----------------------------------------------------------------------------

void GHSInterface::ViewPropertyDeleted( const View& v, const IsoString& property )
{
   if ( GUI != nullptr )
      if ( v == m_currentView )
         SynchronizeWithCurrentView();
}

// ----------------------------------------------------------------------------

bool GHSInterface::WantsReadoutNotifications() const
{
   return true;
}

// ----------------------------------------------------------------------------

void GHSInterface::BeginReadout( const View& v )
{
   if ( GUI != nullptr )
      if ( IsVisible() )
          m_readoutActive = true;
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateReadout( const View& v, const DPoint& p, double R, double G, double B, double A )
{
    if ( !m_readoutActive )
        return;
    
    if (v != m_currentView)
        return;
    
    if ( m_sourceData.IsEmpty() )
       if ( !GetSourceHistograms() )
           return;
    
    ReadoutOptions roo = ReadoutOptions::GetCurrentOptions();
    
    m_readoutSource = FromPreview;
    m_readoutSize = roo.ProbeSize();
    
    m_readoutPoint = p;
    double halfSide = ((double(roo.ProbeSize()) - 1.0) / 2.0);
    m_readoutRect.Set(p.x - halfSide, p.y - halfSide, p.x + halfSide + 1, p.y + halfSide + 1);
        
    int firstChannel = 0;
    int lastChannel = 0;
    ImageVariant img = v.Image();
    if ( v.Image().IsColor() )
    {
        switch (m_instance.p_SC)
        {
            case GHSSC::SC_Green: firstChannel = lastChannel = 1; break;
            case GHSSC::SC_Blue: firstChannel = lastChannel = 2; break;
            case GHSSC::SC_RGB: lastChannel = 2; break;
            case GHSSC::SC_Lightness: img = m_lightnessImage; break;
            case GHSSC::SC_Saturation: img = m_saturationImage; break;
            case GHSSC::SC_Colour: img = m_luminanceImage; break;
        }
    }
    
    double readoutMean = img.Mean(m_readoutRect, firstChannel, lastChannel);
    
    if (m_readoutStatistic == ReadoutOptions::readout_mode::Mean)
    {
        m_readoutValue = readoutMean;
    }
    else
    {
        switch (m_readoutStatistic)
        {
            case ReadoutOptions::readout_mode::Median:
                m_readoutValue = img.Median(m_readoutRect, firstChannel, lastChannel);
                break;
            case ReadoutOptions::readout_mode::Maximum:
                m_readoutValue = img.MaximumSampleValue(m_readoutRect, firstChannel, lastChannel);
                break;
            case ReadoutOptions::readout_mode::Minimum:
                m_readoutValue = img.MinimumSampleValue(m_readoutRect, firstChannel, lastChannel);
                break;
            default:
                m_readoutValue = readoutMean;
        }
    }
    
    img.Free();
     
    UpdateReadoutControls();
    UpdateRealTimePreview();
    UpdateHistograms();
}

// ----------------------------------------------------------------------------

void GHSInterface::EndReadout( const View& v )
{
   if ( m_readoutActive )
   {
       m_readoutActive = false;
       UpdateHistograms();
   }
}

// -----------------------------------------------------------------------------

bool GHSInterface::WantsGlobalNotifications() const
{
    return true;
}

// -----------------------------------------------------------------------------

void GHSInterface::ReadoutOptionsUpdated()
{
    ReadoutOptions roo = ReadoutOptions::GetCurrentOptions();
    
    m_readoutSize = roo.ProbeSize();
    m_readoutStatistic = roo.Mode();
    
    if (( GUI == nullptr ) || ( !IsVisible() ))
        return;
    
    if ( m_readoutSource == FromPreview )
    {
        m_readoutActive = true;
        UpdateReadout(m_currentView, m_readoutPoint, 0.0, 0.0, 0.0, 0.0);
        m_readoutActive = false;
        UpdateReadoutControls();
    }
}

// ----------------------------------------------------------------------------

bool GHSInterface::WantsRealTimePreviewNotifications() const
{
   return true;
}

void GHSInterface::RealTimePreviewOwnerChanged( ProcessInterface& iface )
{
   if ( GUI != nullptr )
   {
       
   }
}

// ----------------------------------------------------------------------------

#define KEY_HT    SettingsKey()

void GHSInterface::SaveSettings() const
{
    IsoString key = KEY_HT;

    Settings::Write( key + "RejectSaturated", m_rejectSaturated );
}

void GHSInterface::LoadSettings()
{
    IsoString key = KEY_HT;
    
    bool rejectSaturated;
    Settings::Read( key + "RejectSaturated", rejectSaturated );
    SetRejectSaturated( rejectSaturated );
}

#undef KEY_HT

// ----------------------------------------------------------------------------

bool GHSInterface::GetSourceHistograms()
{
   m_sourceData.Clear();
   if ( m_currentView.IsNull() )
      return false;
    
    ImageVariant fromImage(m_currentView.Image());
    
    int N = 1;
    if (fromImage.IsColor())
        N = 3;
    for (int i = 0; i < N; ++i)
    {
        Histogram h( m_plotResolution );
        h.SelectChannel(i);
        h << fromImage;
        m_sourceData.Add( h );
    }
       
    if (fromImage.IsColor())
    {
        // Lightness channel
        fromImage.GetLightness(m_lightnessImage);
        Histogram hLight( m_plotResolution );
        hLight << m_lightnessImage;
    
        // Luma channel (Luminance with gamma = 1)
        ImageVariant fromImageCopy;
        fromImageCopy.CopyImage(fromImage);
        RGBColorSystem rgbws = fromImage.RGBWorkingSpace();
        FVector x = rgbws.ChromaticityXCoordinates();
        FVector y = rgbws.ChromaticityYCoordinates();
        FVector Y = rgbws.LuminanceCoefficients();
        FVector newx(x[0], x[1], x[2]);
        FVector newy(y[0], y[1], y[2]);
        FVector newY(Y[0], Y[1], Y[2]);
        RGBColorSystem newRgbws(1, false, newx, newy, newY);
        fromImageCopy.SetRGBWorkingSpace( newRgbws );
        fromImageCopy.GetLuminance(m_luminanceImage);
        Histogram hLum( m_plotResolution );
        hLum << m_luminanceImage;
        
        // Saturation channel
        fromImageCopy.SetColorSpace(ColorSpace::HSV);
        Histogram hSat( m_plotResolution );
        hSat.SelectChannel(1);
        hSat << fromImageCopy;
        fromImageCopy.GetLuminance(m_saturationImage); // to create monochrome image with correct dimensions
        m_saturationImage.AssignImage(fromImageCopy, 0, 1, 1);
        
        m_sourceData.Add(hLight);
        m_sourceData.Add(hSat);
        m_sourceData.Add(hLum);
        
        fromImageCopy.Free();
        fromImage.Free();
    }
    
    return true;
}

// ----------------------------------------------------------------------------

void GHSInterface::CalculateInputHistograms()
{
    m_sourceData.Clear();
    m_inputData.Clear();

    if ( !GetSourceHistograms() )
        return;
    
    int chCount = m_currentView.IsColor() ? 6 : 1;
    
    for ( int c = 0; c < chCount; ++c )
    {
        Histogram h( m_plotResolution );
        m_sourceData[c].Resample(h);
        m_inputData.Add( h );
    }
}

// ----------------------------------------------------------------------------

void GHSInterface::CalculateOutputHistograms()
{
    m_outputData.Clear();

    if ( m_sourceData.IsEmpty() )
        if ( !GetSourceHistograms() )
            return;

    if (!m_currentView.IsColor())
    {
        Histogram h( m_plotResolution );
        m_instance.TransformHistogram( h, m_sourceData[0] );
        m_outputData.Add( h );
        return;
    }
    
    if (m_instance.p_SC < GHSSC::SC_RGB)
    {
        for ( int c = 0; c < 3; ++c )
        {
            Histogram h( m_plotResolution );
            if (c == m_instance.p_SC)
                m_instance.TransformHistogram( h, m_sourceData[c] );
            else
                m_sourceData[c].Resample(h);
            m_outputData.Add( h );
        }
    }
    else if (m_instance.p_SC == GHSSC::SC_RGB)
    {
        for ( int c = 0; c < 3; ++c )
        {
            Histogram h( m_plotResolution );
            m_instance.TransformHistogram( h, m_sourceData[c] );
            m_outputData.Add( h );
        }
    }
    else
    {
        for ( int c = 0; c < 3; ++c )
        {
            Histogram h( m_plotResolution );
            m_sourceData[c].Resample(h);
            m_outputData.Add( h );
        }
    }
    
    
    if (m_instance.p_SC == GHSSC::SC_Lightness)
    {
        Histogram h( m_plotResolution );
        m_instance.TransformHistogram( h, m_sourceData[3] );
        m_outputData.Add( h );
    }
    else
    {
        Histogram h( m_plotResolution );
        m_sourceData[3].Resample(h);
        m_outputData.Add( h );
    }
    
    if (m_instance.p_SC == GHSSC::SC_Saturation)
    {
        Histogram h( m_plotResolution );
        m_instance.TransformHistogram( h, m_sourceData[4] );
        m_outputData.Add( h );
    }
    else
    {
        Histogram h( m_plotResolution );
        m_sourceData[4].Resample(h);
        m_outputData.Add( h );
    }
    
    if (m_instance.p_SC == GHSSC::SC_Colour)
    {
        Histogram h( m_plotResolution );
        m_instance.TransformHistogram( h, m_sourceData[5] );
        m_outputData.Add( h );
    }
    else
    {
        Histogram h( m_plotResolution );
        m_sourceData[5].Resample(h);
        m_outputData.Add( h );
    }
}

// ----------------------------------------------------------------------------

void GHSInterface::CalculateClippingCounts()
{
    m_shadowsCount = 0;
    m_highlightsCount = 0;

    if ( m_sourceData.IsEmpty() )
        if ( !GetSourceHistograms() )
            return;

    ImageVariant image = m_currentView.Image();

    if (!image.IsColor() || (m_instance.p_SC != GHSSC::SC_RGB))
    {
        int cc = 0;
        if (image.IsColor() && (m_instance.p_SC < GHSSC::SC_RGB))
            cc = m_instance.p_SC;
        if (image.IsColor() && (m_instance.p_SC > GHSSC::SC_RGB))
            cc = m_instance.p_SC - 1;
                
        int r = m_sourceData[cc].Resolution();
        int i0 = Max(0, RoundInt( m_instance.p_BP * (r - 1)));
        int i1 = Min(r - 1, RoundInt( m_instance.p_WP * (r - 1)));
        
        for ( int i = 0; i < i0; ++i )
           m_shadowsCount += m_sourceData[cc].Count( i );
        for ( int i = i1; ++i < r; )
           m_highlightsCount += m_sourceData[cc].Count( i );
    }
    else
    {
        int r = m_sourceData[0].Resolution();
        int i0 = Max(0, RoundInt( m_instance.p_BP * (r - 1)));
        int i1 = Min(r - 1, RoundInt( m_instance.p_WP * (r - 1)));
        
        for (int c = 0; c < 3; ++c)
        {
            int m_shadowsCountc = 0;
            int m_highlightsCountc = 0;
            for ( int i = 0; i < i0; ++i )
                m_shadowsCountc += m_sourceData[c].Count( i );
            for ( int i = i1; ++i < r; )
                m_highlightsCountc += m_sourceData[c].Count( i );
            m_shadowsCount = Max(m_shadowsCount, m_shadowsCountc);
            m_highlightsCount = Max(m_highlightsCount, m_highlightsCountc);
        }
    }
    
    image.Free();
}

// ----------------------------------------------------------------------------

void GHSInterface::SetShadowProtect( double x )
{
    m_instance.SetShadowProtect( x );

    GUI->SP_NumericControl.SetValue( m_instance.p_SP );
    GUI->LP_NumericControl.SetValue( m_instance.p_LP );
    GUI->HP_NumericControl.SetValue( m_instance.p_HP );

    CalculateOutputHistograms();
    CalculateClippingCounts();

    UpdateClippingCountControls();
    UpdateHistograms();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

void GHSInterface::SetHighlightProtect( double x )
{
   m_instance.SetHighlightProtect( x );

    GUI->SP_NumericControl.SetValue( m_instance.p_SP );
    GUI->LP_NumericControl.SetValue( m_instance.p_LP );
    GUI->HP_NumericControl.SetValue( m_instance.p_HP );

    CalculateOutputHistograms();
    CalculateClippingCounts();

    UpdateClippingCountControls();
    UpdateHistograms();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

void GHSInterface::SetSymmetryPoint( double x )
{
    m_instance.SetSymmetryPoint( x );

    GUI->SP_NumericControl.SetValue( m_instance.p_SP );
    GUI->LP_NumericControl.SetValue( m_instance.p_LP );
    GUI->HP_NumericControl.SetValue( m_instance.p_HP );

    CalculateOutputHistograms();
    CalculateClippingCounts();

    UpdateClippingCountControls();
    UpdateHistograms();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

void GHSInterface::SetRejectSaturated( bool reject )
{
    m_rejectSaturated = reject;

    UpdateGraphicsControls();
    UpdateHistograms();
}

// ----------------------------------------------------------------------------

void GHSInterface::SetInputZoom( int hz, int vz, const Point* p )
{
    m_settingUp = true;

    m_inputZoomX = hz;
    m_inputZoomY = vz;

    bool hsb = m_inputZoomX > 1;
    bool vsb = m_inputZoomY > 1;

    GUI->InputHistogram_ScrollBox.ShowScrollBars( hsb, vsb );

    int visibleWidth = GUI->InputHistogram_ScrollBox.Viewport().Width();
    int visibleHeight = GUI->InputHistogram_ScrollBox.Viewport().Height();

    int sliderControlSize = RoundInt( LogicalPixelsToPhysical( m_sliderControlSize ) );
    int contentsWidth = visibleWidth * m_inputZoomX;
    int contentsHeight = (visibleHeight - sliderControlSize) * m_inputZoomY + sliderControlSize;

    if ( hsb )
    {
        int m = contentsWidth - visibleWidth;
        GUI->InputHistogram_ScrollBox.SetHorizontalScrollRange( 0, m );

        if ( p != nullptr )
            GUI->InputHistogram_ScrollBox.SetHorizontalScrollPosition(
                  Range( p->x*m_inputZoomX - (visibleWidth >> 1), 0, m ) );
    }
    else
        GUI->InputHistogram_ScrollBox.SetHorizontalScrollRange( 0, 0 );

    GUI->InputHistogram_ScrollBox.SetPageWidth( visibleWidth );

    if ( vsb )
    {
        int m = contentsHeight - visibleHeight;
        GUI->InputHistogram_ScrollBox.SetVerticalScrollRange( 0, m );

        if ( p != nullptr )
            GUI->InputHistogram_ScrollBox.SetVerticalScrollPosition(
                  Range( p->y*m_inputZoomY - (visibleHeight >> 1), 0, m ) );
    }
    else
        GUI->InputHistogram_ScrollBox.SetVerticalScrollRange( 0, 0 );

    GUI->InputHistogram_ScrollBox.SetPageHeight( visibleHeight );

    UpdateZoomControls();
    UpdateHistograms();

    m_settingUp = false;
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateZoomControls()
{
    GUI->HorizontalZoom_SpinBox.SetValue( m_inputZoomX );
    //GUI->VerticalZoom_SpinBox.SetValue( m_inputZoomY );
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateGraphicsControls()
{
    GUI->RejectSaturated_ToolButton.SetChecked( m_rejectSaturated );
    GUI->ShowCurve_ToolButton.SetChecked( m_showMTF );
    GUI->ShowGrid_ToolButton.SetChecked( m_showGrid );
    GUI->LogHistogram_ToolButton.SetChecked( m_logGraph );
    GUI->PreStretchHist_ToolButton.SetChecked( m_showHist );
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateReadoutControls()
{
    String readoutMode;
    switch (m_readoutStatistic)
    {
        case ReadoutOptions::readout_mode::Mean: readoutMode = "Mean value"; break;
        case ReadoutOptions::readout_mode::Median: readoutMode = "Median value"; break;
        case ReadoutOptions::readout_mode::Maximum: readoutMode = "Max value"; break;
        case ReadoutOptions::readout_mode::Minimum: readoutMode = "Min value"; break;
        default: readoutMode = "Undefined value";
    }
    
    
    String readoutChannel;
    if (m_currentView.Image().IsColor())
    {
        switch (m_instance.p_SC)
        {
            case GHSSC::SC_Red: readoutChannel = " in the Red channel"; break;
            case GHSSC::SC_Green: readoutChannel = " in the Green channel"; break;
            case GHSSC::SC_Blue: readoutChannel = " in the Blue channel"; break;
            case GHSSC::SC_RGB: readoutChannel = " in the RGB channels"; break;
            case GHSSC::SC_Lightness: readoutChannel = " in the Lightness channel"; break;
            case GHSSC::SC_Saturation: readoutChannel = " in the Saturation channel"; break;
            case GHSSC::SC_Colour: readoutChannel = " in the weighted average RGB channel"; break;
            default: readoutChannel = " in the Undefined channel";
        }
    }
    else
    {
        readoutChannel = "";
    }
    
    String soSize = String(m_readoutSize);
        
    Rect r = m_readoutRect;
    String sX0 = String(int(r.x0));
    String sY0 = String(int(r.y0));
    String sX1 = String(int(r.x1 - 1));
    String sY1 = String(int(r.y1 - 1));
    String readoutArea = "[" + soSize + "x" + soSize + "]";
    readoutArea += " | [(" + sX0 + ", " + sY0 + ") - (" + sX1 + ", " + sY1 + ")]";
    
    String readoutInfo = readoutMode + readoutChannel;
    
    switch (m_readoutSource)
    {
        case FromPreview:
            GUI->ReadoutSource_Label.SetText("Image");
            GUI->ReadoutSend_PushButton.Enable(true);
            GUI->ReadoutDesc_Label.SetText(readoutInfo);
            GUI->ReadoutArea_Label.SetText(readoutArea);
            break;
        case FromHistogram:
            GUI->ReadoutSource_Label.SetText("Histogram");
            GUI->ReadoutSend_PushButton.Enable(true);
            GUI->ReadoutDesc_Label.SetText("Selected histogram level");
            GUI->ReadoutArea_Label.SetText("Not applicable");
            break;
        default:
            GUI->ReadoutSource_Label.SetText("None");
            GUI->ReadoutSend_PushButton.Enable(false);
            GUI->ReadoutDesc_Label.SetText("No readout data");
            GUI->ReadoutArea_Label.SetText("None");
            m_readoutValue = 0.0;
            break;
    }
    GUI->Readout_Label.SetText(String().Format("%.6f", m_readoutValue));
    
    if (m_instance.p_ST == GHSST::ST_Linear)
    {
        GUI->ReadoutSend_PushButton.SetText("Send to BP");
        GUI->ReadoutSend_PushButton.SetToolTip( "Use readout value to set the value of BP." );
    }
    else
    {
        GUI->ReadoutSend_PushButton.SetText("Send to SP");
        GUI->ReadoutSend_PushButton.SetToolTip( "Use readout value to set the value of SP." );
    }
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateTransformationControls()
{
    GUI->ST_ComboBox.SetCurrentItem(m_instance.p_ST);
    GUI->SC_ComboBox.SetCurrentItem(m_instance.p_SC);
    
    if ((m_instance.p_SC == GHSSC::SC_Colour) && (m_instance.p_ST != GHSST::ST_Linear))
        GUI->ColourOptions_GroupBox.Enable(true);
    else
        GUI->ColourOptions_GroupBox.Enable(false);
    
    if (m_instance.IsInvertible() || (!m_currentView.IsNull() && !m_currentView.IsColor()))
    {
        GUI->Inv_CheckBox.SetChecked(m_instance.p_Inv);
        GUI->Inv_CheckBox.Enable(true);
    }
    else
    {
        GUI->Inv_CheckBox.SetChecked(false);
        GUI->Inv_CheckBox.Enable(false);
    }
    
    GUI->D_NumericControl.SetValue( m_instance.p_D );
    GUI->b_NumericControl.SetValue( m_instance.p_b );
    GUI->SP_NumericControl.SetValue( m_instance.p_SP );
    GUI->LP_NumericControl.SetValue( m_instance.p_LP );
    GUI->HP_NumericControl.SetValue( m_instance.p_HP );
    GUI->BP_NumericControl.SetValue( m_instance.p_BP );
    GUI->WP_NumericControl.SetValue( m_instance.p_WP );
    GUI->CT_ComboBox.SetCurrentItem( m_instance.p_CT );
    GUI->CB_NumericControl.SetValue( m_instance.p_CB );
    GUI->RGBWS_CheckBox.SetChecked( m_instance.p_RGBWS );

    if ( m_instance.p_ST != m_lastST )
    {
        if (m_instance.p_ST != GHSST::ST_GeneralisedHyperbolic)
        {
            GUI->b_NumericControl.SetVisible(false);
            GUI->bReset_ToolButton.SetVisible(false);
        }
        else
        {
            GUI->b_NumericControl.SetVisible(true);
            GUI->bReset_ToolButton.SetVisible(true);
        }
        
        if (m_instance.p_ST != GHSST::ST_Linear)
        {
            GUI->GHSParameters_Control.SetVisible(true);
            GUI->LinearParameters_Control.SetVisible(false);
        }
        else
        {
            GUI->GHSParameters_Control.SetVisible(false);
            GUI->LinearParameters_Control.SetVisible(true);
        }
        m_lastST = m_instance.p_ST;
        EnsureLayoutUpdated();
        AdjustToContents();
    }
    
    GUI->LinearFineAdjust_ComboBox.SetCurrentItem(m_linearFAParameter - 3);
    GUI->FineAdjust_ComboBox.SetCurrentItem(m_nonlinFAParameter);
    GUI->LinearFineAdjust_CheckBox.SetChecked(m_useFinestAdjustment);
    GUI->FineAdjust_CheckBox.SetChecked(m_useFinestAdjustment);
    
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateClippingCountControls()
{
   if ( m_currentView.IsNull() )
   {
       GUI->LCP_NumericControl.SetValue(0.0);
       GUI->LCP_NumericControl.Enable(false);
       GUI->HCP_NumericControl.SetValue(0.0);
       GUI->HCP_NumericControl.Enable(false);
   }
   else
   {
       count_type N  = m_currentView.Image()->NumberOfPixels();
       if (N == 0)
       {
           GUI->LCP_NumericControl.SetValue(0.0);
           GUI->LCP_NumericControl.Enable(false);
           GUI->HCP_NumericControl.SetValue(0.0);
           GUI->HCP_NumericControl.Enable(false);
       }
       else
       {
           count_type n0 = m_shadowsCount;
           count_type n1 = m_highlightsCount;
           GUI->LCP_NumericControl.SetValue(double(n0)/N );
           GUI->LCP_NumericControl.Enable(true);
           GUI->HCP_NumericControl.SetValue(double(n1)/N );
           GUI->HCP_NumericControl.Enable(true);
       }
       
   }
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateHistograms()
{
    UpdateHistogramSliders();
    UpdateHistogramInfo();

    m_inputDirty = true;
    GUI->InputHistogramPlot_Control.Update();
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateHistogramSliders()
{
   m_slidersDirty = true;
   GUI->HistogramSliders_Control.Update();
}

// ----------------------------------------------------------------------------

void GHSInterface::UpdateHistogramInfo()
{
    if ( m_cursorStatus == NoCursor )
    {
        GUI->InputHistogramPlot_Control.SetToolTip( "" );
        return;
    }

    String s;

    int w = GUI->InputHistogramPlot_Control.Width()*m_inputZoomX;
    int dx = (w == 0) ? 1 : Max( 1, RoundInt( double( m_plotResolution )/w ) );

    double x = m_histogramPos.x*(m_plotResolution - 1);
    int i = int( x );
    int j = Min( i+dx, m_plotResolution );
    
    m_instance.UpdateFlags();
    m_instance.UpdateCoeffs();
    double tx = Range(m_histogramPos.x, 0.0, 1.0);
    m_instance.Transform(tx);
    
    if ( !m_currentView.IsNull() )
    {
        if ( m_inputData.IsEmpty() )
            CalculateInputHistograms();
        
        if ( m_outputData.IsEmpty() )
            CalculateOutputHistograms();
    }
    
    // ### N.B. Detect rare cases where a UpdateHistogramInfo is called before we have
    //        regenerated histogram data.
    bool skipHistCalcs = false;
    if (m_currentView.Image().IsColor())
    {
        if ( (int ( m_inputData.Length() ) != 6) || (int ( m_outputData.Length() ) != 6) )
            skipHistCalcs = true;
    }
    else
    {
        if ( (int ( m_inputData.Length() ) != 1) || (int ( m_outputData.Length() ) != 1) )
            skipHistCalcs = true;
    }

    if ( !m_inputData.IsEmpty() && !m_outputData.IsEmpty() && !skipHistCalcs )
    {
        count_type inTotal = 0;
        count_type inCount = 0;
        count_type outTotal = 0;
        count_type outCount = 0;

        if ( !m_currentView.IsColor() )
        {
            int c0 = 0;
            for ( int k = i; k < j; ++k )
            {
                inCount += m_inputData[c0].Count( k );
                outCount += m_outputData[c0].Count( k );
            }
            inTotal += m_inputData[c0].Count();
            outTotal += m_outputData[c0].Count();
        }
        else if ( m_instance.p_SC < GHSSC::SC_RGB )
        {
            int c0 = int(m_instance.p_SC);
            for ( int k = i; k < j; ++k )
            {
                inCount += m_inputData[c0].Count( k );
                outCount += m_outputData[c0].Count( k );
            }
            inTotal += m_inputData[c0].Count();
            outTotal += m_outputData[c0].Count();
        }
        else if ( m_instance.p_SC == GHSSC::SC_RGB )
        {
            for ( int c = 0; c < 3; ++c )
            {
                for ( int k = i; k < j; ++k )
                {
                    inCount += m_inputData[c].Count( k );
                    outCount += m_outputData[c].Count( k );
                }
                inTotal += m_inputData[c].Count();
                outTotal += m_outputData[c].Count();
            }
        }
        else if ( m_instance.p_SC > GHSSC::SC_RGB )
        {
            int c0 = int(m_instance.p_SC - 1);
            for ( int k = i; k < j; ++k )
            {
                inCount += m_inputData[c0].Count( k );
                outCount += m_outputData[c0].Count( k );
            }
            inTotal += m_inputData[c0].Count();
            outTotal += m_outputData[c0].Count();
        }
        
        if (inTotal == 0)
        {
            inTotal = 1;
            inCount = 0;
        }
        
        if (outTotal == 0)
        {
            outTotal = 1;
            outCount = 0;
        }

        s.AppendFormat( "Cursor Position\n  x= %.4f \n  y= %.4f \nTransform(x)= %.4f \nLevel= (%d-%d) \nCurrent count= %llu | %.4f%%\nTransform count= %llu | %.4f%% ", m_histogramPos.x, m_histogramPos.y, tx, i, j-1, inCount, 100.0*inCount/inTotal, outCount, 100.0*outCount/outTotal );
    }

    if ( s.IsEmpty() )
        s.AppendFormat( "Cursor Position\n  x= %.4f \n  y= %.4f \nTransform(x)= %.4f" , m_histogramPos.x, m_histogramPos.y, tx );

    GUI->InputHistogramPlot_Control.SetToolTip( s );
}

// ----------------------------------------------------------------------------

void GHSInterface::SynchronizeWithCurrentView()
{
    if ( !m_currentView.IsNull() )
    {
        if ( m_currentView.IsColor() )
            GUI->ColourOptions_Control.Show();
        else
            GUI->ColourOptions_Control.Hide();
    }
    
    CalculateInputHistograms();
    CalculateOutputHistograms();
    CalculateClippingCounts();

    UpdateHistograms();
    UpdateClippingCountControls();
    UpdateTransformationControls();
    UpdateReadoutControls();
}

// ----------------------------------------------------------------------------

void GHSInterface::PlotGrid(
   Graphics& g, const Rect& r, int width, int height, int hZoom, int vZoom )
{
   int n = 8 * Max(1, Min( hZoom, vZoom ));
   double dx = double( width - 1 )/n;
   double dy = double( height - 1 )/n;
    
    if ((dx == 0) || (dy == 0))
        return;

   int ix0 = int( r.x0/dx );
   int ix1 = int( r.x1/dx );

   int iy0 = int( r.y0/dy );
   int iy1 = int( r.y1/dy );

   int w = r.Width();
   int h = r.Height();

   Pen p0( m_gridColour0, DisplayPixelRatio(), PenStyle::Solid );
   Pen p1( m_gridColour1, DisplayPixelRatio(), PenStyle::Dot );

   for ( int i = ix0; i <= ix1; ++i )
   {
      int x = RoundInt( dx*i ) - r.x0;
      if ( x >= w )
         break;

      g.SetPen( (i & 1) ? p1 : p0 );
      g.DrawLine( x, 0, x, h );
   }

   for ( int i = iy0; i <= iy1; ++i )
   {
      int y = RoundInt( dy*i ) - r.y0;
      if ( y >= h )
         break;

      g.SetPen( (i & 1) ? p1 : p0 );
      g.DrawLine( 0, y, w, y );
   }
}

// ----------------------------------------------------------------------------

void GHSInterface::PlotMonoHistogramValues(
   Graphics& g, const Rect& r, DVector plotValues, double peak, int width, int height, int hZoom, graph_style graphStyle, RGBA colour = RGBAColor(0x80, 0x80, 0x80))
{
    if (plotValues.Length() == 0)
        return;
    
    g.SetPen( colour, DisplayPixelRatio() );
    
    //double peakToUse = Max(peak, plotValues.MaxComponent());
    double peakToUse = peak;
    
    if (m_logGraph)
        peakToUse = Ln(1.0 + peakToUse);
    peakToUse = Max(1.0, peakToUse);
    
    Array<Point> points;
    for (int i = 0; i < plotValues.Length(); ++i)
    {
        int y;
        if (m_logGraph)
        {
            y = Max( -1, RoundInt( (height - 1)*(1 - Ln(1.0 + double( *plotValues.At(i) ) )/peakToUse) ) - r.y0 ); // -1 to allow clipping at top
        }
        else
        {
            y = Max( -1, RoundInt( (height - 1)*(1 - double( *plotValues.At(i) )/peakToUse) ) - r.y0 ); // -1 to allow clipping at top
        }
        points.Add( Point( i, y ) );
    }
    
    
    switch ( graphStyle )
    {
        default: // ?!
        case LineStyle:
            for ( Array<Point>::const_iterator i0 = points.Begin(), i = i0; ++i < points.End(); i0 = i )
                g.DrawLine( *i0, *i );
            break;

        case AreaStyle:
            points.Add( Point( (*points.ReverseBegin()).x, height-1 ) );
            points.Add( Point( (*points.Begin()).x, height-1 ) );
            g.SetBrush( g.Pen().Color() );
            g.DrawPolygon( points );
            break;

        case BarStyle:
            for ( Array<Point>::const_iterator i = points.Begin(); i != points.End(); ++i )
                g.DrawLine( *i, Point( i->x, height-1 ) );
            break;

        case DotStyle:
            for ( Array<Point>::const_iterator i = points.Begin(); i != points.End(); ++i )
                g.DrawPoint( *i );
            break;
    }
    
}

// ----------------------------------------------------------------------------

void GHSInterface::PlotColourHistogramValues(
   Graphics& g, const Rect& r, DVector plotRedValues, DVector plotGreenValues, DVector plotBlueValues, double peak, int width, int height, int hZoom, graph_style graphStyle, int frontChannel = -1)
{
    int valuesLength = Min(Min(plotRedValues.Length(), plotGreenValues.Length()), plotBlueValues.Length());
    double peakToUse = peak;
    if (valuesLength == 0)
        return;
    
    if (graphStyle == AreaStyle)
    {
        PlotMonoHistogramValues(g, r, plotRedValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0xff, 0x80, 0x80));
        PlotMonoHistogramValues(g, r, plotGreenValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0x80, 0xff, 0x80));
        PlotMonoHistogramValues(g, r, plotBlueValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0x80, 0x80, 0xff));
        
        DVector plotCyanValues(0.0, valuesLength);
        DVector plotMagentaValues(0.0, valuesLength);
        DVector plotYellowValues(0.0, valuesLength);
        DVector plotWhiteValues(0.0, valuesLength);
        for (int i = 0; i < valuesLength; ++i)
        {
            *plotCyanValues.At(i) = Min(*plotGreenValues.At(i), *plotBlueValues.At(i));
            *plotMagentaValues.At(i) = Min(*plotBlueValues.At(i), *plotRedValues.At(i));
            *plotYellowValues.At(i) = Min(*plotRedValues.At(i), *plotGreenValues.At(i));
            *plotWhiteValues.At(i) = Min(*plotCyanValues.At(i), *plotRedValues.At(i));
        }
        PlotMonoHistogramValues(g, r, plotCyanValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0x80, 0xff, 0xff));
        PlotMonoHistogramValues(g, r, plotMagentaValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0xff, 0x80, 0xff));
        PlotMonoHistogramValues(g, r, plotYellowValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0xff, 0xff, 0x80));
        PlotMonoHistogramValues(g, r, plotWhiteValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0x80, 0x80, 0x80));
        switch (frontChannel)
        {
            case 0:
                PlotMonoHistogramValues(g, r, plotRedValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0xff, 0x80, 0x80));
                break;
            case 1:
                PlotMonoHistogramValues(g, r, plotGreenValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0x80, 0xff, 0x80));
                break;
            case 2:
                PlotMonoHistogramValues(g, r, plotBlueValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0x80, 0x80, 0xff));
                break;
            default:
                break;
        }
    }
    else
    {
        PlotMonoHistogramValues(g, r, plotRedValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0xff, 0xc0, 0xc0));
        PlotMonoHistogramValues(g, r, plotGreenValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0xc0, 0xff, 0xc0));
        PlotMonoHistogramValues(g, r, plotBlueValues, peakToUse, width, height, hZoom, graphStyle, RGBAColor(0xc0, 0xc0, 0xff));
    }
}

// ----------------------------------------------------------------------------

DVector GHSInterface::GetHistogramPlotValues( Graphics& g, const Rect& r, const Histogram& H, int width, int hZoom, double& peak)
{
    int vLen = r.x1 - r.x0;
    DVector plotValues( 0.0, vLen );
    
    int histRes = H.Resolution();
    if (histRes == 0)
        return plotValues;
    double dx = double(hZoom * width)/double(histRes);
    if (dx == 0)
        return plotValues;
    
    for (int x = r.x0; x < r.x1; ++x)
    {
        double z0 = double(x) / dx;
        int i0 = int(Floor( z0 ));
        
        double z1 = double(x + 1) / dx;
        int i1 = int(Floor( z1 ));
        
        for (int j = i0; j < i1; ++j)
            *plotValues.At(x - r.x0) += double(H[j]);
        
        *plotValues.At(x - r.x0) -= double(H[i0]) * (z0 - i0);
        *plotValues.At(x - r.x0) += double(H[i1]) * (z1 - i1);
    }
    
    int histPeakLevel = 0;
    int x0 = 0;
    int x1 = 0;
    DVector maxAreaValues( 0.0, 20 );
    
    if ( m_rejectSaturated )
    {
        peak = *MaxItem(plotValues.Begin()+1, plotValues.End()-1);
        
        int margin = int(Ceil(1.0/dx));
        UI64Vector trimmedHistData(H.HistogramData().Begin() + margin, histRes - (2 * margin));
        histPeakLevel = trimmedHistData.IndexOfLargestComponent();
        int xm = int(Floor( histPeakLevel * dx ));
        x0 = Max(1, xm - 10);
        x1 = Min(xm + 10, hZoom * width - 2);
        
        for (int x = x0; x < x1; ++x)
        {
            double z0 = double(x) / dx;
            int i0 = int(Floor( z0 ));
            
            double z1 = double(x + 1) / dx;
            int i1 = int(Floor( z1 ));
            
            for (int j = i0; j < i1; ++j)
                *maxAreaValues.At(x - x0) += double(H[j]);
            
            *maxAreaValues.At(x - x0) -= double(H[i0]) * (z0 - i0);
            *maxAreaValues.At(x - x0) += double(H[i1]) * (z1 - i1);
        }
    }
    else
    {
        peak = plotValues.MaxComponent();
        
        histPeakLevel = H.PeakLevel();
        int xm = int(Floor( histPeakLevel * dx ));
        x0 = Max(0, xm - 10);
        x1 = Min(xm + 10, hZoom * width - 1);
        
        for (int x = x0; x < x1; ++x)
        {
            double z0 = double(x) / dx;
            int i0 = int(Floor( z0 ));
            
            double z1 = double(x + 1) / dx;
            int i1 = int(Floor( z1 ));
            
            for (int j = i0; j < i1; ++j)
                *maxAreaValues.At(x - x0) += double(H[j]);
            
            *maxAreaValues.At(x - x0) -= double(H[i0]) * (z0 - i0);
            *maxAreaValues.At(x - x0) += double(H[i1]) * (z1 - i1);
        }
    }
    
    peak = Max(peak, maxAreaValues.MaxComponent());
    return plotValues;
}

// ----------------------------------------------------------------------------

void GHSInterface::PlotScale( Graphics& g, const Rect& r, int width )
{
   int w = r.Width();
   int h = r.Height();

   bool c0 = m_channel == 0 || m_channel >= 3;
   bool c1 = m_channel == 1 || m_channel >= 3;
   bool c2 = m_channel == 2 || m_channel >= 3;
    
    if (!(w > 1))
        return;

   float v0 = float( r.x0 )/(width - 1);
   float v1 = float( r.x1 - 1 )/(width - 1);

   GradientBrush::stop_list stops;
   stops.Add( GradientBrush::Stop( 0.0, RGBAColor( c0 ? v0 : 0.0F, c1 ? v0 : 0.0F, c2 ? v0 : 0.0F ) ) );
   stops.Add( GradientBrush::Stop( 1.0, RGBAColor( c0 ? v1 : 0.0F, c1 ? v1 : 0.0F, c2 ? v1 : 0.0F ) ) );

   g.FillRect( 0, 0, w, h, LinearGradientBrush( 0, 0, w, 0, stops ) );
}

// ----------------------------------------------------------------------------

RGBA GHSInterface::HandlerColor( double v ) const
{
   // Ensure visibility of handlers on R, G, B and gray backgrounds.
   if ( (m_channel == 0 || m_channel == 1 || m_channel >= 3) && v > 0.5 )
      return 0xFF000000;
   return 0xFFFFFFFF;
}

// ----------------------------------------------------------------------------

void GHSInterface::PlotHandler( Graphics& g, double v, int x0, int width )
{
   int h = RoundInt( LogicalPixelsToPhysical( m_sliderControlSize ) );
   int h2 = (h >> 1) + 1;

   int x = RoundInt( v*(width - 1) ) - x0;

   GenericVector<Point> notch( 4 );
   notch[0] = Point( x,      h-h2 );
   notch[1] = Point( x-h2+1, h-1  );
   notch[2] = Point( x+h2-1, h-1  );
   notch[3] = Point( x,      h-h2 );

   g.SetPen( HandlerColor( v ), DisplayPixelRatio() );
   g.DrawLine( x, 0, x, h-h2-1 );
   g.DrawPolyline( notch );
}

// ----------------------------------------------------------------------------

void GHSInterface::PlotGHSCurve( Graphics& g, const Rect& r, int width, int height )
{
    int xc0 = 0;    //RoundInt( m_instance.ShadowsClipping( m_channel )*(width - 1) );
    if ( xc0 >= r.x1 )
        return;
    int xc1 = (width - 1);    //RoundInt( m_instance.HighlightsClipping( m_channel )*(width - 1) );
    if ( xc1 < r.x0 )
        return;
    g.SetPen( m_transformColour, 2 * DisplayPixelRatio() );
    
    
    if ( xc1 - xc0 < 2 )
    {
        g.DrawLine( xc0-r.x0, 0, xc1-r.x0, r.Height()-1 );
        return;
    }
    double dx = 1.0/(xc1 - xc0);

    Array<Point> points;

    int px0 = Max( r.x0-1, xc0 );
    int px1 = Min( r.x1, xc1 );
    
    m_instance.UpdateFlags();
    m_instance.UpdateCoeffs();

    for ( int xi = px0, x0 = -1, y0 = -1; xi < px1; ++xi )
    {
        int x = xi - r.x0;
        double yValue = ( xi - xc0 )*dx;
        m_instance.Transform( yValue );

        int y = RoundInt( (height - 1)*(1 - yValue) ) - r.y0;
        if ( x != x0 || y != y0 )
            points.Add( Point( x0 = x, y0 = y ) );
    }

    double yValue = ( px1 - xc0)*dx;
    m_instance.Transform( yValue );
    points.Add( Point( px1 - r.x0, RoundInt( (height - 1)*(1 - yValue) ) - r.y0 ) );

    g.DrawPolyline( points );
     
}

// ----------------------------------------------------------------------------

void GHSInterface::PlotReadouts( Graphics& g,
                        const Bitmap& bmp, const Rect& r, const DVector& readouts, int width, int height )
{
   int w = bmp.Width();
   int h = bmp.Height();
   float d = DisplayPixelRatio();
    
    int x = RoundInt( m_readoutValue*(width - 1) ) - r.x0;
    if ( x >= 0 && x < w )
    {
       g.SetPen( m_readoutColour, 2 * d );
       g.DrawLine( x, 0, x, h );
    }
}

// ----------------------------------------------------------------------------

void GHSInterface::PlotIdentityTransformation( Graphics& g,
                        const Bitmap& bmp, const Rect& r, int width, int height )
{
    float d = DisplayPixelRatio();

    int x0 = - r.x0;
    int y0 = height - 1 - r.y0;
    int x1 = width - 1 - r.x0;
    int y1 = - r.y0;
    g.SetPen( m_identityColour, d );
    g.DrawLine( x0, y0, x1, y1 );
      
}

// ----------------------------------------------------------------------------

void GHSInterface::PlotCursor( Graphics& g, const Rect& r )
{
    int w = r.Width();
    int h = r.Height();
    int x = m_cursorPos.x - r.x0;
    int y = m_cursorPos.y - r.y0;

    g.SetPen( m_cursorColour, DisplayPixelRatio() );
    if ( x >= 0 && x < w )
        g.DrawLine( x, 0, x, h );
    if ( y >= 0 && y < h )
        g.DrawLine( 0, y, w, y );
}

// ----------------------------------------------------------------------------

void GHSInterface::RegenerateInputViewport()
{
    Rect r0 = GUI->InputHistogramPlot_Control.BoundsRect();
    int w0 = r0.Width();
    int h0 = r0.Height();

    m_inputDirty = false;

    if ( m_inputBitmap.IsNull() )
        m_inputBitmap = Bitmap( w0, h0, BitmapFormat::RGB32 );
    m_inputBitmap.Fill( m_backgroundColour );

    Rect r( r0 + GUI->InputHistogram_ScrollBox.ScrollPosition() );
    int w = w0*m_inputZoomX;
    int h = h0*m_inputZoomY;

    if ( m_showGrid )
    {
        Graphics g( m_inputBitmap );
        g.EnableAntialiasing();
        PlotGrid( g, r, w, h, m_inputZoomX, m_inputZoomY );
    }

    if ( m_currentView.IsNull() )
        return;

    if ( m_inputData.IsEmpty() )
       CalculateInputHistograms();
    
    if ( m_outputData.IsEmpty() )
       CalculateOutputHistograms();

    histogram_list inH = m_inputData;
    histogram_list outH = m_outputData;
    
    ImageVariant image = m_currentView.Image();
    
    // ### N.B. Detect rare cases where a paint event is sent before we have
    //        regenerated histogram data; e.g. Paint() before ImageFocused().
    
    if (image.IsColor())
    {
        if ( (int ( inH.Length() ) != 6) || (int ( outH.Length() ) != 6) )
            return;
    }
    else
    {
        if ( (int ( inH.Length() ) != 1) || (int ( outH.Length() ) != 1) )
            return;
    }
     

   if ( !inH.IsEmpty() && !outH.IsEmpty() )
   {
      Bitmap bmp( w0, h0, BitmapFormat::RGB32 );
      {
          bmp.Fill( m_backgroundColour );
          Graphics g( bmp );
          g.EnableAntialiasing();
          if (m_showGrid)
              PlotGrid( g, r, w, h, m_inputZoomX, m_inputZoomY );
          g.SetCompositionOperator( CompositionOp::Source );
          
          if ( !image.IsColor() )
          {
              double peak = 0.0;
              DVector valuesToPlot = GetHistogramPlotValues(g, r, outH[0], w0, m_inputZoomX, peak);
              PlotMonoHistogramValues(g, r, valuesToPlot, peak, w, h, m_inputZoomX, m_outputGraphStyle, RGBAColor(0x80, 0x80, 0x80));
              
              if (m_showHist)
              {
                  valuesToPlot = GetHistogramPlotValues(g, r, inH[0], w0, m_inputZoomX, peak);
                  PlotMonoHistogramValues(g, r, valuesToPlot, peak, w, h, m_inputZoomX, m_inputGraphStyle, RGBAColor(0xC0, 0xC0, 0xC0));
              }
          }
          else
          {
              if ( m_instance.p_SC <= GHSSC::SC_RGB )
              {
                  int frontColour = -1;
                  if (m_instance.p_SC < GHSSC::SC_RGB)
                      frontColour = m_instance.p_SC;
                  
                  double peak = 0.0;
                  DVector redValuesToPlot = GetHistogramPlotValues(g, r, outH[0], w0, m_inputZoomX, peak);
                  double maxPeak = peak;
                  DVector greenValuesToPlot = GetHistogramPlotValues(g, r, outH[1], w0, m_inputZoomX, peak);
                  maxPeak = Max(maxPeak, peak);
                  DVector blueValuesToPlot = GetHistogramPlotValues(g, r, outH[2], w0, m_inputZoomX, peak);
                  maxPeak = Max(maxPeak, peak);
                  PlotColourHistogramValues(g, r, redValuesToPlot, greenValuesToPlot, blueValuesToPlot, maxPeak, w, h, m_inputZoomX, m_outputGraphStyle, frontColour);
                  
                  if (m_showHist)
                  {
                      peak = 0.0;
                      redValuesToPlot = GetHistogramPlotValues(g, r, inH[0], w0, m_inputZoomX, peak);
                      maxPeak = peak;
                      greenValuesToPlot = GetHistogramPlotValues(g, r, inH[1], w0, m_inputZoomX, peak);
                      maxPeak = Max(maxPeak, peak);
                      blueValuesToPlot = GetHistogramPlotValues(g, r, inH[2], w0, m_inputZoomX, peak);
                      maxPeak = Max(maxPeak, peak);
                      PlotColourHistogramValues(g, r, redValuesToPlot, greenValuesToPlot, blueValuesToPlot, maxPeak, w, h, m_inputZoomX, m_inputGraphStyle, frontColour);
                  }
              }
              else
              {
                  int channel = int(m_instance.p_SC) - 1;
                  
                  double peak = 0.0;
                  DVector valuesToPlot = GetHistogramPlotValues(g, r, outH[channel], w0, m_inputZoomX, peak);
                  PlotMonoHistogramValues(g, r, valuesToPlot, peak, w, h, m_inputZoomX, m_outputGraphStyle, RGBAColor(0x80, 0x80, 0x80));
                  
                  if (m_showHist)
                  {
                      valuesToPlot = GetHistogramPlotValues(g, r, inH[channel], w0, m_inputZoomX, peak);
                      PlotMonoHistogramValues(g, r, valuesToPlot, peak, w, h, m_inputZoomX, m_inputGraphStyle, RGBAColor(0xC0, 0xC0, 0xC0));
                  }
              }
          }
      }

      Graphics g( m_inputBitmap );
      g.SetCompositionOperator( CompositionOp::Source );
      g.DrawBitmap( 0, 0, bmp );
   }
    
    image.Free();
}

// ----------------------------------------------------------------------------

void GHSInterface::RegenerateSlidersViewport()
{
    Rect r0 = GUI->HistogramSliders_Control.BoundsRect();
    int w0 = r0.Width();
    int h0 = r0.Height();

    if ( m_slidersBitmap.IsNull() )
        m_slidersBitmap = Bitmap( w0, h0, BitmapFormat::RGB32 );

    m_slidersDirty = false;

    Graphics g( m_slidersBitmap );
    Rect r( r0 );
    r += GUI->InputHistogram_ScrollBox.ScrollPosition();
    PlotScale( g, r, w0*m_inputZoomX );
}

// ----------------------------------------------------------------------------

void GHSInterface::__Histogram_Paint( Control& sender, const pcl::Rect& updateRect )
{
   if ( GUI == nullptr )
      return;

   if ( sender == GUI->InputHistogramPlot_Control )
   {
      if ( m_inputDirty )
         RegenerateInputViewport();

      if ( m_showMTF || (m_readoutSource != FromNone) || m_cursorStatus == InputCursor )
      {
          Bitmap bmp = m_inputBitmap.Subimage( updateRect );
          {
              Graphics g( bmp );
              g.EnableAntialiasing();
              g.SetCompositionOperator( CompositionOp::Source );

             
              Rect r0 = sender.ClientRect();
              int w = r0.Width()*m_inputZoomX;
              int h = r0.Height()*m_inputZoomY;
                
              Rect r( updateRect + GUI->InputHistogram_ScrollBox.ScrollPosition() );
              
              if ( m_showMTF )
              {
                  PlotIdentityTransformation( g, bmp, r, w, h );
                  PlotGHSCurve( g, r, w, h );
              }
             
              if ( m_cursorStatus == InputCursor )
                  PlotCursor( g, updateRect );
               
              if (m_readoutSource != FromNone)
                  PlotReadouts( g, bmp, r, m_inputReadouts, w, h );
         }

         Graphics g( sender );
         g.DrawBitmap( updateRect.LeftTop(), bmp );
      }
      else
      {
         Graphics g( sender );
         g.DrawBitmapRect( updateRect.LeftTop(), m_inputBitmap, updateRect );
      }
   }
}

// ----------------------------------------------------------------------------

void GHSInterface::__Sliders_Paint( Control& sender, const pcl::Rect& updateRect )
{
   if ( GUI == nullptr )
      return;

   if ( m_slidersDirty )
      RegenerateSlidersViewport();
    
    double c0, c1;
    double m = 0.0;
    
    if (m_instance.p_ST == GHSST::ST_Linear)
    {
        c0 = m_instance.p_BP;
        c1 = m_instance.p_WP;
    }
    else
    {
        c0 = m_instance.p_LP;
        m = m_instance.p_SP;
        c1 = m_instance.p_HP;
    }

   

    int w = sender.Width()*m_inputZoomX;
    int x0 = GUI->InputHistogram_ScrollBox.HorizontalScrollPosition();

    Graphics g( sender );
    g.EnableAntialiasing();
    g.DrawBitmapRect( updateRect.LeftTop(), m_slidersBitmap, updateRect );
    PlotHandler( g, c0, x0, w );
    PlotHandler( g, c1, x0, w );
    if (m_instance.p_ST != GHSST::ST_Linear)
        PlotHandler( g, m, x0, w );
}

// ----------------------------------------------------------------------------

void GHSInterface::__Histogram_Resize( Control& sender,
                                             int/*newWidth*/, int/*newHeight*/, int/*oldWidth*/, int/*oldHeight*/ )
{
   if ( GUI == nullptr )
      return;

   if ( sender == GUI->InputHistogramPlot_Control )
   {
      m_inputBitmap = Bitmap::Null();
      m_inputDirty = true;

      if ( !m_settingUp )
         SetInputZoom( m_inputZoomX, m_inputZoomY );
   }
   else if ( sender == GUI->HistogramSliders_Control )
   {
      m_slidersBitmap = Bitmap::Null();
      m_slidersDirty = true;
   }
}

// ----------------------------------------------------------------------------

void GHSInterface::__Histogram_ScrollPosUpdated( ScrollBox& sender, int pos )
{
   if ( sender == GUI->InputHistogram_ScrollBox )
       UpdateHistograms();
}

// ----------------------------------------------------------------------------

void GHSInterface::__Histogram_Enter( Control& sender )
{
   if ( sender == GUI->InputHistogramPlot_Control )
      m_cursorStatus = InputCursor;
   m_cursorPos = -1;
}

// ----------------------------------------------------------------------------

void GHSInterface::__Histogram_Leave( Control& sender )
{
   m_cursorStatus = NoCursor;
   UpdateHistogramInfo();
   sender.Update();
}

// ----------------------------------------------------------------------------

void GHSInterface::__Histogram_MousePress( Control& sender,
                                             const pcl::Point& pos, int button, unsigned buttons, unsigned modifiers )
{
    m_readoutValue = m_histogramPos.x;
    m_readoutSource = FromHistogram;
    
    Rect r = sender.ClientRect();
    int w = r.Width();
    int h = r.Height();
    
    Point p = pos + GUI->InputHistogram_ScrollBox.ScrollPosition();
    if ((w > 1) && (h > 1))
    {
        m_histogramPos.x = Range( double( p.x )/(w - 1), 0.0, 1.0 );
        m_histogramPos.y = Range( 1 - double( p.y )/(h - 1), 0.0, 1.0 );

        UpdateHistogramInfo();
    }
    
    UpdateHistograms();
    UpdateReadoutControls();
}

// ----------------------------------------------------------------------------

void GHSInterface::__Histogram_MouseMove( Control& sender,
                                             const pcl::Point& pos, unsigned buttons, unsigned modifiers )
{
    Rect r = sender.ClientRect();
    int w = r.Width();
    int h = r.Height();

    for ( int i = 0; i < 2; ++i )
    {
        double f = DisplayPixelRatio();
        int ui1 = RoundInt( f );
        
        sender.Update( m_cursorPos.x-ui1-ui1, 0, m_cursorPos.x+ui1+ui1, h );
        sender.Update( 0, m_cursorPos.y-ui1-ui1, w, m_cursorPos.y+ui1+ui1 );
        
         if ( i == 0 )
            m_cursorPos = pos;
      }

      w *= m_inputZoomX;

      Point p = pos + GUI->InputHistogram_ScrollBox.ScrollPosition();

       if ((w > 1) && (h > 1))
       {
           m_histogramPos.x = Range( double( p.x )/(w - 1), 0.0, 1.0 );
           m_histogramPos.y = Range( 1 - double( p.y )/(h - 1), 0.0, 1.0 );

           UpdateHistogramInfo();
       }
}

// ----------------------------------------------------------------------------

GHSInterface::slider_id GHSInterface::FindHandler( double v ) const
{
    if (m_instance.p_ST == GHSST::ST_Linear)
    {
        double c0 = m_instance.p_BP;
        double c1 = m_instance.p_HP;

        double dc0 = Abs( v - c0 );
        double dc1 = Abs( v - c1 );

        if ( dc0 <= dc1 )
            return C0Slider;
        
        return C1Slider;
    }
    else
    {
        double c0 = m_instance.p_LP;
        double c1 = m_instance.p_HP;
        double cm = m_instance.p_SP;

        double dcm = Abs( v - cm );
        double dc0 = Abs( v - c0 );
        double dc1 = Abs( v - c1 );

        if ( dcm <= dc0 )
        {
            if ( dcm <= dc1 )
                return MSlider;
        }

        if ( dc0 <= dcm )
        {
            if ( dc0 <= dc1 )
                return C0Slider;
        }

        return C1Slider;
    }
    
}

// ----------------------------------------------------------------------------

double GHSInterface::SliderToHistogram( int x ) const
{
   return double( x + GUI->InputHistogram_ScrollBox.HorizontalScrollPosition() ) /
               (GUI->HistogramSliders_Control.Width()*m_inputZoomX - 1);
}

// ----------------------------------------------------------------------------

void GHSInterface::__Sliders_MousePress( Control& sender,
                                             const pcl::Point& pos, int button, unsigned buttons, unsigned modifiers )
{
   if ( button != MouseButton::Left )
      return;

   m_sliderBeingDragged = FindHandler( SliderToHistogram( pos.x ) );

   __Sliders_MouseMove( sender, pos, buttons, modifiers );
}

// ----------------------------------------------------------------------------

void GHSInterface::__Sliders_MouseRelease(
   Control& sender, const pcl::Point& pos, int button, unsigned buttons, unsigned modifiers )
{
   __Sliders_MouseMove( sender, pos, buttons, modifiers );
   m_sliderBeingDragged = NoSlider;
}

// ----------------------------------------------------------------------------

void GHSInterface::__Sliders_MouseMove( Control& sender,
                                             const pcl::Point& pos, unsigned buttons, unsigned modifiers )
{
   if ( m_sliderBeingDragged != NoSlider )
   {
      double v = SliderToHistogram( pos.x );

      if (m_instance.p_ST == GHSST::ST_Linear)
      {
          if ( m_sliderBeingDragged == C0Slider )
          {
             v = Round( v, TheGHSBPParameter->Precision() );
              m_instance.SetBlackpoint( v );
          }
          else if ( m_sliderBeingDragged == C1Slider )
          {
             v = Round( v, TheGHSWPParameter->Precision() );
              m_instance.SetWhitepoint( v );
          }
      }
      else
      {
          if ( m_sliderBeingDragged == MSlider )
          {
              v = Round( v, TheGHSSPParameter->Precision() );
              m_instance.SetSymmetryPoint( v );
          }
          else if ( m_sliderBeingDragged == C0Slider )
          {
             v = Round( v, TheGHSLPParameter->Precision() );
              m_instance.SetShadowProtect( v );
          }
          else if ( m_sliderBeingDragged == C1Slider )
          {
             v = Round( v, TheGHSHPParameter->Precision() );
              m_instance.SetHighlightProtect( v );
          }
      }

       CalculateOutputHistograms();
       CalculateClippingCounts();
       UpdateHistograms();
       UpdateTransformationControls();
       UpdateClippingCountControls();
       UpdateRealTimePreview();
   }
}

// ----------------------------------------------------------------------------

void GHSInterface::__Reset_ButtonClick( Button& sender, bool /*checked*/ )
{
    if ( sender == GUI->DReset_ToolButton )
        m_instance.SetStretchFactor( 0.0 );
    if ( sender == GUI->bReset_ToolButton )
        m_instance.SetLocalIntensity( 0.0 );
    if ( sender == GUI->SPReset_ToolButton )
        m_instance.SetSymmetryPoint( 0.0 );
    if ( sender == GUI->LPReset_ToolButton )
        m_instance.SetShadowProtect( 0.0 );
    if ( sender == GUI->HPReset_ToolButton )
        m_instance.SetHighlightProtect( 1.0 );
    if ( sender == GUI->CBReset_ToolButton )
        m_instance.SetColourBlend( 1.0 );
    if ( sender == GUI->BPReset_ToolButton )
        m_instance.SetBlackpoint( 0.0 );
    if ( sender == GUI->WPReset_ToolButton )
        m_instance.SetWhitepoint( 1.0 );
    if ( sender == GUI->LCPReset_ToolButton )
        SetBlackpointFromLCP( 0.0 );
    if ( sender == GUI->HCPReset_ToolButton )
        SetWhitepointFromHCP( 0.0 );
    
    CalculateOutputHistograms();
    CalculateClippingCounts();
    UpdateControls();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

void GHSInterface::__AutoZero_ButtonClick( Button& sender, bool /*checked*/ )
{
    
   if ( m_sourceData.IsEmpty() )
      if ( !GetSourceHistograms() )
         return;

   ImageVariant image = m_currentView.Image();
    
    if (!image.IsColor() || (m_instance.p_SC != GHSSC::SC_RGB))
    {
        int cc = 0;
        if (image.IsColor() && (m_instance.p_SC < GHSSC::SC_RGB))
            cc = m_instance.p_SC;
                
        int r1 = m_sourceData[cc].Resolution() - 1;
        if (r1 == 0)
            return;
        
        if ( sender == GUI->LCPReset_ToolButton )
        {
            int c0;
            for ( c0 = 0; c0 < r1; ++c0 )
                if ( m_sourceData[cc].Count( c0 ) != 0 )
                    break;
            m_instance.SetBlackpoint(double(c0)/r1);
        }
        
        else if ( sender == GUI->HCPReset_ToolButton )
        {
            int c1;
            for ( c1 = 0; c1 < r1; ++c1 )
                if ( m_sourceData[cc].Count( c1 ) != 0 )
                    break;
            m_instance.SetWhitepoint(double(c1)/r1);
        }
    }
    else
    {
        int r1 = m_sourceData[0].Resolution() - 1;
        
        if ( sender == GUI->LCPReset_ToolButton )
        {
            int c0;
            for ( c0 = 0; c0 < r1; ++c0 )
                if (( m_sourceData[0].Count( c0 ) != 0 ) || ( m_sourceData[1].Count( c0 ) != 0 ) || ( m_sourceData[2].Count( c0 ) != 0 ))
                    break;
            m_instance.SetBlackpoint(double(c0)/r1);
        }
        
        else if ( sender == GUI->HCPReset_ToolButton )
        {
            int c1;
            for ( c1 = 0; c1 < r1; ++c1 )
                if (( m_sourceData[0].Count( c1 ) != 0 ) || ( m_sourceData[1].Count( c1 ) != 0 ) || ( m_sourceData[2].Count( c1 ) != 0 ))
                    break;
            m_instance.SetWhitepoint(double(c1)/r1);
        }
    }

    CalculateOutputHistograms();
    CalculateClippingCounts();

    UpdateTransformationControls();
    UpdateHistograms();
    UpdateClippingCountControls();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

double GHSInterface::CalculateLowClipPoint( double clipRatio )
{
    if ( m_sourceData.IsEmpty() )
       if ( !GetSourceHistograms() )
          return 0.0;
    
    ImageVariant image = m_currentView.Image();
    count_type N = image->NumberOfPixels();
    count_type N0 = count_type( clipRatio*N );
    
    int c = 0;
    if (image.IsColor() && (m_instance.p_SC < GHSSC::SC_RGB))
        c = m_instance.p_SC;
    if (image.IsColor() && (m_instance.p_SC > GHSSC::SC_RGB))
        c = m_instance.p_SC - 1;
            
    int r1 = m_sourceData[c].Resolution() - 1;
    if (r1 == 0)
        return 0.0;
    
    int c0;

    
     
    if (!image.IsColor() || (m_instance.p_SC != GHSSC::SC_RGB))
    {
        count_type n = 0;
        for ( c0 = 0; c0 < r1; ++c0 )
            if ( (n += m_sourceData[c].Count( c0 )) > N0 )
               break;
            
        return (double(c0)/r1);
    }
    else
    {
        count_type n0 = 0;
        count_type n1 = 0;
        count_type n2 = 0;
        for ( c0 = 0; c0 < r1; ++c0 )
        {
            n0 += m_sourceData[0].Count( c0 );
            n1 += m_sourceData[1].Count( c0 );
            n2 += m_sourceData[2].Count( c0 );
            if ((n0 > N0) || (n1 > N0) || (n2 > N0))
                break;
        }
        return (double(c0)/r1);
             
     }
}

// ----------------------------------------------------------------------------

void GHSInterface::SetBlackpointFromLCP( double clipRatio )
{
    m_instance.SetBlackpoint( CalculateLowClipPoint( clipRatio ) );
}

// ----------------------------------------------------------------------------

double GHSInterface::CalculateHighClipPoint( double clipRatio )
{
    if ( m_sourceData.IsEmpty() )
       if ( !GetSourceHistograms() )
          return 0.0;
    
    ImageVariant image = m_currentView.Image();
    count_type N = image->NumberOfPixels();
    count_type N1 = count_type( clipRatio*N );
    
    int c = 0;
    if (image.IsColor() && (m_instance.p_SC < GHSSC::SC_RGB))
        c = m_instance.p_SC;
    if (image.IsColor() && (m_instance.p_SC > GHSSC::SC_RGB))
        c = m_instance.p_SC - 1;
            
    int r1 = m_sourceData[c].Resolution() - 1;
    if (r1 == 0)
        return 1.0;
    
    int c1;
     
    if (!image.IsColor() || (m_instance.p_SC != GHSSC::SC_RGB))
    {
        count_type n = 0;
        for ( c1 = r1; c1 > 0; --c1 )
           if ( (n += m_sourceData[c].Count( c1 )) > N1 )
              break;
            
        return (double(c1)/r1);
    }
    else
    {
        int r1 = m_sourceData[0].Resolution() - 1;
         
        count_type n0 = 0;
        count_type n1 = 0;
        count_type n2 = 0;
        for ( c1 = r1; c1 > 0; --c1 )
        {
            n0 += m_sourceData[0].Count( c1 );
            n1 += m_sourceData[1].Count( c1 );
            n2 += m_sourceData[2].Count( c1 );
            if ((n0 > N1) || (n1 > N1) || (n2 > N1))
                break;
        }
        return (double(c1)/r1);
     }
}

// ----------------------------------------------------------------------------

void GHSInterface::SetWhitepointFromHCP( double clipRatio )
{
    m_instance.SetWhitepoint( CalculateHighClipPoint( clipRatio ) );
}

// ----------------------------------------------------------------------------

void GHSInterface::__Zoom_ButtonClick( Button& sender, bool /*checked*/ )
{
    if ( sender == GUI->Zoom11_ToolButton )
    {
        SetInputZoom( 1, 1 );
    }
    if ( sender == GUI->ZoomIn_ToolButton )
    {
        SetInputZoom( Range(250*(1 + Floor(m_inputZoomX/250)), 1, 999), 1 );
    }
    if ( sender == GUI->ZoomOut_ToolButton )
    {
        SetInputZoom( Range(250*Floor(m_inputZoomX/251), 1, 999), 1 );
    }
}

// ----------------------------------------------------------------------------

void GHSInterface::__Readout_ButtonClick( Button& sender, bool /*checked*/ )
{
    if ( sender == GUI->ReadoutSend_PushButton )
        if (m_readoutSource != FromNone)
            if (m_instance.p_ST == GHSST::ST_Linear)
            {
                m_instance.SetBlackpoint( m_readoutValue );
            }
            else
            {
                m_instance.SetSymmetryPoint( m_readoutValue );
            }
    if ( sender == GUI->ReadoutClear_PushButton )
        m_readoutSource = FromNone;
    
    CalculateOutputHistograms();
    UpdateControls();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

void GHSInterface::__Zoom_ValueUpdated( SpinBox& sender, int value )
{
   if ( sender == GUI->HorizontalZoom_SpinBox )
   {
       if (m_inputZoomX > 0)
       {
           int visibleWidth = GUI->InputHistogram_ScrollBox.Viewport().Width();
           int hsPosition = GUI->InputHistogram_ScrollBox.HorizontalScrollPosition();
           Point p( (hsPosition + (visibleWidth>>1))/m_inputZoomX, 0);
           
           SetInputZoom( value, m_inputZoomY, &p );
       }
       else
       {
           SetInputZoom( value, m_inputZoomY );
       }
   }
}

// ----------------------------------------------------------------------------

void GHSInterface::__RejectSaturated_ButtonClick( Button& /*sender*/, bool checked )
{
    SetRejectSaturated( checked );
}

// ----------------------------------------------------------------------------

void GHSInterface::__ShowCurve_ButtonClick( Button& /*sender*/, bool checked )
{
    m_showMTF = checked;
    UpdateHistograms();
}

// ----------------------------------------------------------------------------

void GHSInterface::__ShowGrid_ButtonClick( Button& /*sender*/, bool checked )
{
   m_showGrid = checked;
   UpdateHistograms();
}

// ----------------------------------------------------------------------------

void GHSInterface::__LogHistogram_ButtonClick( Button& /*sender*/, bool checked )
{
   m_logGraph = checked;
   UpdateHistograms();
}

// ----------------------------------------------------------------------------

void GHSInterface::__PreStretchHist_ButtonClick( Button& /*sender*/, bool checked )
{
    m_showHist = checked;
    UpdateHistograms();
}

// ----------------------------------------------------------------------------
/*
void GHSInterface::__Web_ButtonClick( Button& sender, bool checked )
{
    GUI->GHSWebView.LoadContent("https://www.ghsastro.co.uk/");
    GUI->GHSWebView.SetVisible(true);
}

void GHSInterface::__Web_Close( Control& sender, bool& allowClose )
{
    GUI->GHSWebView.SetPlainText(String("loading..."));
}
*/
// ----------------------------------------------------------------------------


void GHSInterface::__RealValueUpdated( NumericEdit& sender, double value )
{
    if ( sender == GUI->D_NumericControl )
       m_instance.SetStretchFactor( value );
    if ( sender == GUI->b_NumericControl )
        m_instance.SetLocalIntensity( value );
    if ( sender == GUI->SP_NumericControl )
        m_instance.SetSymmetryPoint( value );
    if ( sender == GUI->LP_NumericControl )
        m_instance.SetShadowProtect( value );
    if ( sender == GUI->HP_NumericControl )
        m_instance.SetHighlightProtect( value );
    if ( sender == GUI->CB_NumericControl )
        m_instance.SetColourBlend( value );
    if ( sender == GUI->BP_NumericControl )
        m_instance.SetBlackpoint( value );
    if ( sender == GUI->WP_NumericControl )
        m_instance.SetWhitepoint( value );
    if ( sender == GUI->LCP_NumericControl )
        SetBlackpointFromLCP( value );
    if ( sender == GUI->HCP_NumericControl )
        SetWhitepointFromHCP( value );
    
    CalculateOutputHistograms();
    CalculateClippingCounts();
    UpdateControls();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

void GHSInterface::__ItemClicked( Button& sender, bool checked )
{
    if ( sender == GUI->Inv_CheckBox )
        m_instance.p_Inv = checked;
    else if ( sender == GUI->RGBWS_CheckBox )
        m_instance.p_RGBWS = checked;
    else if ( sender == GUI->FineAdjust_CheckBox )
        m_useFinestAdjustment = checked;
    else if ( sender == GUI->LinearFineAdjust_CheckBox )
        m_useFinestAdjustment = checked;
    
    CalculateOutputHistograms();
    CalculateClippingCounts();
    UpdateControls();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

void GHSInterface::__ItemSelected( ComboBox& sender, int itemIndex )
{
    if ( sender == GUI->ST_ComboBox )
       m_instance.p_ST = itemIndex;
    else if ( sender == GUI->SC_ComboBox )
    {
        m_instance.p_SC = itemIndex;
        m_channel = itemIndex;
        CalculateInputHistograms();
        
        if ( m_readoutSource == FromPreview )
        {
            m_readoutActive = true;
            UpdateReadout(m_currentView, m_readoutPoint, 0.0, 0.0, 0.0, 0.0);
            m_readoutActive = false;
        }
    }
    else if ( sender == GUI->CT_ComboBox )
       m_instance.p_CT = itemIndex;
    else if ( sender == GUI->FineAdjust_ComboBox )
        m_nonlinFAParameter = fa_parameter( itemIndex );
    else if ( sender == GUI->LinearFineAdjust_ComboBox )
        m_linearFAParameter = fa_parameter( 3 + itemIndex );
    
    CalculateOutputHistograms();
    CalculateClippingCounts();
    UpdateControls();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------

void GHSInterface::__UpdateRealTimePreview_Timer( Timer& sender )
{
   if ( m_realTimeThread.IsValid() )
      if ( m_realTimeThread->IsActive() )
         return;
   if ( IsRealTimePreviewActive() )
      RealTimePreview::Update();
}

// ----------------------------------------------------------------------------

void GHSInterface::__SliderValueUpdated( Slider& sender, int value )
{
    if ( sender == GUI->FineAdjust_Slider )
    {
        int d = m_useFinestAdjustment ? 1000000 : 10000;
        switch (m_nonlinFAParameter)
        {
            case SP: m_instance.SetSymmetryPoint( m_FAPBeforeAdjustment + double( value - 500 )/d ); break;
            case LP: m_instance.SetShadowProtect( m_FAPBeforeAdjustment + double( value - 500 )/d ); break;
            case HP: m_instance.SetHighlightProtect( m_FAPBeforeAdjustment + double( value - 500 )/d ); break;
            default:;
        }
    }
    else if ( sender == GUI->LinearFineAdjust_Slider )
    {
        int d = m_useFinestAdjustment ? 1000000 : 10000;
        switch (m_linearFAParameter)
        {
            case BP: m_instance.SetBlackpoint( m_FAPBeforeAdjustment + double( value - 500 )/d ); break;
            case LCP: SetBlackpointFromLCP( m_FAPBeforeAdjustment + double( value - 500 )/d ); break;
            case WP: m_instance.SetWhitepoint( m_FAPBeforeAdjustment + double( value - 500 )/d ); break;
            case HCP: SetWhitepointFromHCP( m_FAPBeforeAdjustment + double( value - 500 )/d ); break;
            default:;
        }
    }
    CalculateOutputHistograms();
    CalculateClippingCounts();
    UpdateControls();
    UpdateRealTimePreview();
}

// ----------------------------------------------------------------------------

void GHSInterface::__FineAdjustSliderEnter( Control& sender )
{
   __FineAdjustSliderGetFocus( sender );
}

// ----------------------------------------------------------------------------

void GHSInterface::__FineAdjustSliderLeave( Control& sender )
{
   __FineAdjustSliderLoseFocus( sender );
}

// ----------------------------------------------------------------------------

void GHSInterface::__FineAdjustSliderMousePress( Control& sender, const pcl::Point& pos, int button, unsigned buttons, unsigned modifiers )
{
   __FineAdjustSliderGetFocus( sender );
}

// ----------------------------------------------------------------------------

void GHSInterface::__FineAdjustSliderMouseRelease( Control& sender, const pcl::Point& pos, int button, unsigned buttons, unsigned modifiers )
{
   __FineAdjustSliderLoseFocus( sender );
}

// ----------------------------------------------------------------------------

void GHSInterface::__FineAdjustSliderLoseFocus( Control& sender )
{
    if ( sender == GUI->FineAdjust_Slider )
        GUI->FineAdjust_Slider.SetValue( 500 );
    if ( sender == GUI->LinearFineAdjust_Slider )
        GUI->LinearFineAdjust_Slider.SetValue( 500 );
}

// ----------------------------------------------------------------------------

void GHSInterface::__FineAdjustSliderGetFocus( Control& sender )
{
    if ( sender == GUI->FineAdjust_Slider )
    {
        if ( m_nonlinFAParameter == SP )
            m_FAPBeforeAdjustment = m_instance.p_SP;
        else if ( m_nonlinFAParameter == LP )
            m_FAPBeforeAdjustment = m_instance.p_LP;
        else if ( m_nonlinFAParameter == HP )
            m_FAPBeforeAdjustment = m_instance.p_HP;
    }
    else if ( sender == GUI->LinearFineAdjust_Slider )
    {
        if ( m_linearFAParameter == BP )
            m_FAPBeforeAdjustment = m_instance.p_BP;
        else if ( m_linearFAParameter == LCP )
            m_FAPBeforeAdjustment = GUI->LCP_NumericControl.Value();
        else if ( m_linearFAParameter == WP )
            m_FAPBeforeAdjustment = m_instance.p_WP;
        else if ( m_linearFAParameter == HCP )
            m_FAPBeforeAdjustment = GUI->HCP_NumericControl.Value();
    }
}

// ----------------------------------------------------------------------------

GHSInterface::GUIData::GUIData( GHSInterface& w )
{
    pcl::Font font = w.Font();
    int labelWidth1 = font.Width( String( "Protect highlights (HP):  " ) ); // the longest label text
    int editWidth1 = font.Width( String( '0', 10 ) );
    int ri16 = w.LogicalPixelsToResource( 16 );
    
    // ------------------------------------------------------------------
    // Histogram controls

    HorizontalZoom_SpinBox.SetRange( 1, s_maxZoom );
    HorizontalZoom_SpinBox.SetToolTip( "Histogram horizontal zoom" );
    HorizontalZoom_SpinBox.OnValueUpdated( (SpinBox::value_event_handler)&GHSInterface::__Zoom_ValueUpdated, w );

    ZoomIn_ToolButton.SetIcon( w.ScaledResource( ":/toolbar/view-zoom-in.png" ) );
    ZoomIn_ToolButton.SetScaledFixedSize( 20, 20 );
    ZoomIn_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    ZoomIn_ToolButton.SetToolTip( "Big zoom in" );
    ZoomIn_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Zoom_ButtonClick, w );
    
    ZoomOut_ToolButton.SetIcon( w.ScaledResource( ":/toolbar/view-zoom-out.png" ) );
    ZoomOut_ToolButton.SetScaledFixedSize( 20, 20 );
    ZoomOut_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    ZoomOut_ToolButton.SetToolTip( "Big zoom out" );
    ZoomOut_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Zoom_ButtonClick, w );
    
    Zoom11_ToolButton.SetIcon( w.ScaledResource( ":/toolbar/view-zoom-1-1.png" ) );
    Zoom11_ToolButton.SetScaledFixedSize( 20, 20 );
    Zoom11_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    Zoom11_ToolButton.SetToolTip( "Zoom 1:1" );
    Zoom11_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Zoom_ButtonClick, w );
    
    //
    
    RejectSaturated_ToolButton.SetIcon( Bitmap::FromSVGFile( "@module_icons_dir/reject_saturated.svg", ri16, ri16 ) );
    RejectSaturated_ToolButton.SetScaledFixedSize( 20, 20 );
    RejectSaturated_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    RejectSaturated_ToolButton.SetToolTip( "Reject saturated pixels for histogram representations" );
    RejectSaturated_ToolButton.SetCheckable();
    RejectSaturated_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__RejectSaturated_ButtonClick, w );

    ShowCurve_ToolButton.SetIcon( Bitmap::FromSVGFile( "@module_icons_dir/show_curve.svg", ri16, ri16 ) );
    ShowCurve_ToolButton.SetScaledFixedSize( 20, 20 );
    ShowCurve_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    ShowCurve_ToolButton.SetToolTip( "Show GHS curve" );
    ShowCurve_ToolButton.SetCheckable();
    ShowCurve_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__ShowCurve_ButtonClick, w );

    ShowGrid_ToolButton.SetIcon( Bitmap::FromSVGFile( "@module_icons_dir/show_grid.svg", ri16, ri16 ) );
    ShowGrid_ToolButton.SetScaledFixedSize( 20, 20 );
    ShowGrid_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    ShowGrid_ToolButton.SetToolTip( "Show grids" );
    ShowGrid_ToolButton.SetCheckable();
    ShowGrid_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__ShowGrid_ButtonClick, w );
    
    LogHistogram_ToolButton.SetIcon( Bitmap::FromSVGFile( "@module_icons_dir/log_histogram.svg", ri16, ri16 ) );
    LogHistogram_ToolButton.SetScaledFixedSize( 20, 20 );
    LogHistogram_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    LogHistogram_ToolButton.SetToolTip( "Plot histogram on a log scale. Useful for seeing smaller histogram counts relatingto stars, for example." );
    LogHistogram_ToolButton.SetCheckable();
    LogHistogram_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__LogHistogram_ButtonClick, w );
    
    PreStretchHist_ToolButton.SetIcon( Bitmap::FromSVGFile( "@module_icons_dir/show_hist.svg", ri16, ri16 ) );
    PreStretchHist_ToolButton.SetScaledFixedSize( 20, 20 );
    PreStretchHist_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    PreStretchHist_ToolButton.SetToolTip( "Show pre-stretch histogram." );
    PreStretchHist_ToolButton.SetCheckable();
    PreStretchHist_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__PreStretchHist_ButtonClick, w );

    Graphics_Sizer.SetSpacing( 4 );
    Graphics_Sizer.Add( RejectSaturated_ToolButton );
    Graphics_Sizer.Add( LogHistogram_ToolButton );
    Graphics_Sizer.Add( ShowCurve_ToolButton );
    Graphics_Sizer.Add( ShowGrid_ToolButton );
    Graphics_Sizer.Add( PreStretchHist_ToolButton );
    Graphics_Sizer.AddStretch();
    Graphics_Sizer.Add( HorizontalZoom_SpinBox );
    Graphics_Sizer.Add( ZoomIn_ToolButton );
    Graphics_Sizer.Add( ZoomOut_ToolButton );
    Graphics_Sizer.Add( Zoom11_ToolButton );
    
    // ------------------------------------------------------------------
    // Histogram graph

    InputHistogramPlot_Control.EnableMouseTracking();
    //InputHistogramPlot_Control.SetCursor( StdCursor::NoCursor );
    InputHistogramPlot_Control.OnPaint( (Control::paint_event_handler)&GHSInterface::__Histogram_Paint, w );
    InputHistogramPlot_Control.OnResize( (Control::resize_event_handler)&GHSInterface::__Histogram_Resize, w );
    InputHistogramPlot_Control.OnEnter( (Control::event_handler)&GHSInterface::__Histogram_Enter, w );
    InputHistogramPlot_Control.OnLeave( (Control::event_handler)&GHSInterface::__Histogram_Leave, w );
    InputHistogramPlot_Control.OnMouseMove( (Control::mouse_event_handler)&GHSInterface::__Histogram_MouseMove, w );
    InputHistogramPlot_Control.OnMousePress( (Control::mouse_button_event_handler)&GHSInterface::__Histogram_MousePress, w );

    HistogramSliders_Control.EnableMouseTracking();
    HistogramSliders_Control.SetCursor( StdCursor::UpArrow );
    HistogramSliders_Control.SetScaledFixedHeight( w.m_sliderControlSize );
    HistogramSliders_Control.OnPaint( (Control::paint_event_handler)&GHSInterface::__Sliders_Paint, w );
    HistogramSliders_Control.OnResize( (Control::resize_event_handler)&GHSInterface::__Histogram_Resize, w );
    HistogramSliders_Control.OnMouseMove( (Control::mouse_event_handler)&GHSInterface::__Sliders_MouseMove, w );
    HistogramSliders_Control.OnMousePress( (Control::mouse_button_event_handler)&GHSInterface::__Sliders_MousePress, w );
    HistogramSliders_Control.OnMouseRelease( (Control::mouse_button_event_handler)&GHSInterface::__Sliders_MouseRelease, w );

    InputHistogramViewport_Sizer.Add( InputHistogramPlot_Control, 100 );
    InputHistogramViewport_Sizer.Add( HistogramSliders_Control );

    InputHistogram_ScrollBox.DisableAutoScroll();
    InputHistogram_ScrollBox.SetScaledMinSize( w.m_minHistogramWidth, w.m_minHistogramHeight+w.m_sliderControlSize );
    InputHistogram_ScrollBox.OnHorizontalScrollPosUpdated( (ScrollBox::pos_event_handler)&GHSInterface::__Histogram_ScrollPosUpdated, w );
    InputHistogram_ScrollBox.OnVerticalScrollPosUpdated( (ScrollBox::pos_event_handler)&GHSInterface::__Histogram_ScrollPosUpdated, w );

    InputHistogram_ScrollBox.Viewport().SetSizer( InputHistogramViewport_Sizer );
    
    InputHistogram_Sizer.SetSpacing( 3 );
    InputHistogram_Sizer.Add( Graphics_Sizer );
    InputHistogram_Sizer.Add( InputHistogram_ScrollBox );
    InputHistogram_Control.SetSizer( InputHistogram_Sizer );
    
    InputHistogram_SectionBar.SetTitle( "Graph" );
    InputHistogram_SectionBar.SetSection( InputHistogram_Control );

    // ------------------------------------------------------------------
    // Colour mode controls
    
    SC_Label.SetText( "Mode:" );
    SC_Label.SetMinWidth( labelWidth1 );
    SC_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    
    SC_ComboBox.AddItem( "Red" );
    SC_ComboBox.AddItem( "Green" );
    SC_ComboBox.AddItem( "Blue" );
    SC_ComboBox.AddItem( "RGB" );
    SC_ComboBox.AddItem( "Lightness" );
    SC_ComboBox.AddItem( "Saturation" );
    SC_ComboBox.AddItem( "Colour" );
    
    SC_ComboBox.SetToolTip( "<p>This parameter specifies how the transformation is to be applied. It can be applied to each of the R, G and B channels individually or all three.  It can also be applied to the CIE Lightness channel or the Saturation.  The final Colour option will stretch all three RGB channels by the same factor calculated by reference to the weighted average of the RGB values.  This reflects the approach used in the standard arcsinh process.</p>" );
    SC_ComboBox.OnItemSelected( (ComboBox::item_event_handler)&GHSInterface::__ItemSelected, w );
    
    SC_Sizer.SetSpacing( 4 );
    SC_Sizer.Add( SC_Label );
    SC_Sizer.Add( SC_ComboBox );
    SC_Sizer.AddStretch();

    //
    
    CT_Label.SetText( "Clip type:" );
    CT_Label.SetMinWidth( labelWidth1 );
    CT_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    
    CT_ComboBox.AddItem( "Clip" );
    CT_ComboBox.AddItem( "Rescale" );
    CT_ComboBox.AddItem( "RGBBlend" );
    CT_ComboBox.AddItem( "RescaleGlobal" );
    
    CT_ComboBox.SetToolTip( "<p>This parameter specifies what to do when undertaking a colour stretch if any channels relating to a specific pixel stretch to a value greater than 1.  Clip: will simply clip all channels greater than 1 for that pixel back to 1. Rescale: will scale all three channels for that pixel down by the ratio 1/max(R,G,B).  RGBBlend will blend in sufficient of an RGB stretch into that pixel to ensure all channels are not greater than 1.</p>" );
    CT_ComboBox.OnItemSelected( (ComboBox::item_event_handler)&GHSInterface::__ItemSelected, w );
    
    CT_Sizer.SetSpacing( 4 );
    CT_Sizer.Add( CT_Label );
    CT_Sizer.Add( CT_ComboBox );
    CT_Sizer.AddStretch();
    
    //
    
    CB_NumericControl.label.SetText( "Colour blend:" );
    CB_NumericControl.label.SetFixedWidth( labelWidth1 );
    CB_NumericControl.slider.SetScaledMinWidth( 250 );
    CB_NumericControl.slider.SetRange( 1, 1000 );
    CB_NumericControl.SetReal();
    CB_NumericControl.SetRange( TheGHSCBParameter->MinimumValue(), TheGHSCBParameter->MaximumValue() );
    CB_NumericControl.SetPrecision( TheGHSCBParameter->Precision() );
    CB_NumericControl.SetToolTip( "<p>When using the colour stretch method this parameter will blend a proportion of a standard RGB stretch.  This will mute the colour saturation for a strong initial data stretch.  A value of 1 will not add any RGB stretch, a value of 0 is equivalent to a normal RGB stretch.</p>" );
    CB_NumericControl.edit.SetFixedWidth( editWidth1 );
    CB_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );
    
    //
    
    CBReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    CBReset_ToolButton.SetScaledFixedSize( 20, 20 );
    CBReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    CBReset_ToolButton.SetToolTip( "Reset CB parameter" );
    CBReset_ToolButton.SetCheckable(false);
    CBReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    //
    
    CB_Sizer.SetSpacing( 4 );
    CB_Sizer.Add( CB_NumericControl );
    CB_Sizer.Add( CBReset_ToolButton );
    
    //
    
    RGBWS_CheckBox.SetText( "Use RGB working space" );
    RGBWS_CheckBox.SetToolTip( "<p>Check here to use coefficients from the RGB working space of the active image to derive the RGB weighted average for a colour stretch.  If unchecked then equal weighted coefficients will be used.</p>" );
    RGBWS_CheckBox.OnClick( (pcl::Button::click_event_handler)&GHSInterface::__ItemClicked, w );

    RGBWS_Sizer.AddUnscaledSpacing( labelWidth1 + w.LogicalPixelsToPhysical( 4 ) );
    RGBWS_Sizer.Add( RGBWS_CheckBox );
    RGBWS_Sizer.AddStretch();
       
    ColourOptionsGB_Sizer.SetSpacing( 3 );
    ColourOptionsGB_Sizer.Add( CT_Sizer );
    ColourOptionsGB_Sizer.Add( CB_Sizer );
    ColourOptionsGB_Sizer.Add( RGBWS_Sizer );
    
    ColourOptions_GroupBox.SetTitle("Colour mode options");
    ColourOptions_GroupBox.SetSizer( ColourOptionsGB_Sizer );
    
    ColourOptions_Sizer.SetSpacing( 3 );
    ColourOptions_Sizer.Add( SC_Sizer );
    ColourOptions_Sizer.Add( ColourOptions_GroupBox );
    ColourOptions_Control.SetSizer( ColourOptions_Sizer );
    
    ColourOptions_SectionBar.SetTitle( "Colour Options" );
    ColourOptions_SectionBar.SetSection( ColourOptions_Control );
    
    // ------------------------------------------------------------------
    // GHS parameter controls
    
    ST_Label.SetText( "Transformation type:" );
    ST_Label.SetMinWidth( labelWidth1 );
    ST_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    
    ST_ComboBox.AddItem( "Generalised Hyperbolic" );
    ST_ComboBox.AddItem( "Midtone Transfer" );
    ST_ComboBox.AddItem( "Arcsinh" );
    ST_ComboBox.AddItem( "Linear" );
    ST_ComboBox.SetToolTip( "<p>Specifies the transformation equations to use. In most cases this will be the Generalised Hyperbolic Stretch equations designed for this module. Other options include the midtone transfer function used in the Histogram Transformation process, and the arcsinh function.  To both these the module brings additional functionality. Alternatively a linear stretch can be selected.</p>" );
    ST_ComboBox.OnItemSelected( (ComboBox::item_event_handler)&GHSInterface::__ItemSelected, w );
    
    ST_Sizer.SetSpacing( 4 );
    ST_Sizer.Add( ST_Label );
    ST_Sizer.Add( ST_ComboBox );
    ST_Sizer.AddStretch();

    //
    
    Inv_CheckBox.SetText( "Invert" );
    Inv_CheckBox.SetToolTip( "<p>Check here to use the inverse form of the transformation equations. <b>Warning:</b> Some transformations can involve clipping and are not truly invertible (e.g. Linear strteches and Non-linear saturation or colour stretches. Inversion is disabled for these cases.</p>" );
    Inv_CheckBox.OnClick( (pcl::Button::click_event_handler)&GHSInterface::__ItemClicked, w );

    Inv_Sizer.AddUnscaledSpacing( labelWidth1 + w.LogicalPixelsToPhysical( 4 ) );
    Inv_Sizer.Add( Inv_CheckBox );
    Inv_Sizer.AddStretch();
    
   //
    
    D_NumericControl.label.SetText( "Stretch factor (ln(D+1)):" );
    D_NumericControl.label.SetFixedWidth( labelWidth1 );
    D_NumericControl.slider.SetScaledMinWidth( 250 );
    D_NumericControl.slider.SetRange( 1, 1000 );
    D_NumericControl.SetReal();
    D_NumericControl.SetRange( TheGHSParameter->MinimumValue(), TheGHSParameter->MaximumValue() );
    D_NumericControl.SetPrecision( TheGHSParameter->Precision() );
    D_NumericControl.SetToolTip( "<p>Controls the amount of stretch. D is a variable that independently controls the contrast added (the slope of the stretch transform) at SP, thus adjusting the amount of stretch applied to the rest of the image.  D does not change the 'form' of the stretch, simply the amount.  D should be used in tandem with b to control the distribution of contrast and brightness. When D is set to zero, the stretch transform will be the identity (y=x) or 'no stretch' transform.</p>" );
    D_NumericControl.edit.SetFixedWidth( editWidth1 );
    D_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );
    
    DReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    DReset_ToolButton.SetScaledFixedSize( 20, 20 );
    DReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    DReset_ToolButton.SetToolTip( "Reset D parameter" );
    DReset_ToolButton.SetCheckable(false);
    DReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    D_Sizer.SetSpacing( 4 );
    D_Sizer.Add( D_NumericControl );
    D_Sizer.Add( DReset_ToolButton );
    
    //
    
    b_NumericControl.label.SetText( "Local intensity (b):" );
    b_NumericControl.label.SetFixedWidth( labelWidth1 );
    b_NumericControl.slider.SetScaledMinWidth( 250 );
    b_NumericControl.slider.SetRange( 1, 1000 );
    b_NumericControl.SetReal();
    b_NumericControl.SetRange( TheGHSbParameter->MinimumValue(), TheGHSbParameter->MaximumValue() );
    b_NumericControl.SetPrecision( TheGHSbParameter->Precision() );
    b_NumericControl.SetToolTip( "<p>Controls how tightly focused the stretch is around the Symetry point by changing the form of the transform itself. For concentrated stretches (such as initial stretches on linear images) a large +ve b factor should be employed to focus a stretch within a histogram peak while de-focusing the stretch away from the histogram peak (such as bright stars). For adjustment of non-linear images, lower or -ve b (and/or lower D) parameters should be employed to distribute contrast and brightness more evenly.  Large positive values of 'b' can be thought of as a histogram widener, ie spreading the histogram wider about the focus point, SP.  By contrast, lower and -ve values of b tend to shift the histogram to a brighter (or dimmer) position without affecting its width too greatly. As a general rule, the level of b employed will decrease as a stretch sequence nears completion, although larger +ve b values (with small D) can still be employed for precise placement of additional contrast.</p>" );
    b_NumericControl.edit.SetFixedWidth( editWidth1 );
    b_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );
    
    bReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    bReset_ToolButton.SetScaledFixedSize( 20, 20 );
    bReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    bReset_ToolButton.SetToolTip( "Reset b parameter" );
    bReset_ToolButton.SetCheckable(false);
    bReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    b_Sizer.SetSpacing( 4 );
    b_Sizer.Add( b_NumericControl );
    b_Sizer.Add( bReset_ToolButton );

    //

    SP_NumericControl.label.SetText( "Symmetry point (SP):" );
    SP_NumericControl.label.SetFixedWidth( labelWidth1 );
    SP_NumericControl.slider.SetScaledMinWidth( 250 );
    SP_NumericControl.slider.SetRange( 1, 1000 );
    SP_NumericControl.SetReal();
    SP_NumericControl.SetRange( TheGHSSPParameter->MinimumValue(), TheGHSSPParameter->MaximumValue() );
    SP_NumericControl.SetPrecision( TheGHSSPParameter->Precision() );
    SP_NumericControl.SetToolTip( "<p>Sets the focus point around which the stretch is applied - contrast will be distributed symmetrically about SP.  While 'b' provides the degree of focus of the stretch, SP determines where that focus is applied.  SP should generally be placed within a histogram peak so that the stretch will widen and lower the peak by adding the most contrast in the stretch at that point.  Pixel values will move away from the SP location. This parameter must be greater than or equal to the Shadow protection parameter and less than or equal to the Highlight protection parameter.</p>" );
    SP_NumericControl.edit.SetFixedWidth( editWidth1 );
    SP_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );
    
    SPReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    SPReset_ToolButton.SetScaledFixedSize( 20, 20 );
    SPReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    SPReset_ToolButton.SetToolTip( "Reset SP parameter" );
    SPReset_ToolButton.SetCheckable(false);
    SPReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    SP_Sizer.SetSpacing( 4 );
    SP_Sizer.Add( SP_NumericControl );
    SP_Sizer.Add( SPReset_ToolButton );
    
    //

    LP_NumericControl.label.SetText( "Protect shadows (LP):" );
    LP_NumericControl.label.SetFixedWidth( labelWidth1 );
    LP_NumericControl.slider.SetScaledMinWidth( 250 );
    LP_NumericControl.slider.SetRange( 1, 1000 );
    LP_NumericControl.SetReal();
    LP_NumericControl.SetRange( TheGHSLPParameter->MinimumValue(), TheGHSLPParameter->MaximumValue() );
    LP_NumericControl.SetPrecision( TheGHSLPParameter->Precision() );
    LP_NumericControl.SetToolTip( "<p>Sets a value below which the stretch is modified to preserve contrast in the shadows/lowlights. This is done by performing a linear stretch of the data below the 'LP' level by reserving contrast from the rest of the image. Moving the LP level towards the current setting of SP changes both the scope (range) and the amount of this contrast reservation, the net effect is to push the overal stretch to higher brightness levels while keeping the contrast and definition in the background.  The amount of contrast reserved for the lowlights is such that the continuity of the stretch is preserved. This parameter must be greater than or equal to 0 and not greater than the Symmetry point.</p>" );
    LP_NumericControl.edit.SetFixedWidth( editWidth1 );
    LP_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );
    
    LPReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    LPReset_ToolButton.SetScaledFixedSize( 20, 20 );
    LPReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    LPReset_ToolButton.SetToolTip( "Reset LP parameter" );
    LPReset_ToolButton.SetCheckable(false);
    LPReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    LP_Sizer.SetSpacing( 4 );
    LP_Sizer.Add( LP_NumericControl );
    LP_Sizer.Add( LPReset_ToolButton );

    //

    HP_NumericControl.label.SetText( "Protect highlights (HP):" );
    HP_NumericControl.label.SetFixedWidth( labelWidth1 );
    HP_NumericControl.slider.SetScaledMinWidth( 250 );
    HP_NumericControl.slider.SetRange( 1, 1000 );
    HP_NumericControl.SetReal();
    HP_NumericControl.SetRange( TheGHSHPParameter->MinimumValue(), TheGHSHPParameter->MaximumValue() );
    HP_NumericControl.SetPrecision( TheGHSHPParameter->Precision() );
    HP_NumericControl.SetToolTip( "<p>Sets a value above which the stretch is modified to preserve contrast in the highlights/stars. This is done by performing a linear stretch of the data above the 'HP' level by reserving contrast from the rest of the image. Moving the HP level towards the current setting of SP increases both the scope (range) and the amount of this contrast reservation, the net effect is to push the overal stretch to lower brightness levels while keeping the contrast and definition in the highlights.  The amount of contrast reserved for the highlights is such that the continuity of the stretch is preserved. This parameter must be less than or equal to 1 and not less than the Symmetry point.</p>" );
    HP_NumericControl.edit.SetFixedWidth( editWidth1 );
    HP_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );
    
    
    HPReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    HPReset_ToolButton.SetScaledFixedSize( 20, 20 );
    HPReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    HPReset_ToolButton.SetToolTip( "Reset HP parameter" );
    HPReset_ToolButton.SetCheckable(false);
    HPReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    HP_Sizer.SetSpacing( 4 );
    HP_Sizer.Add( HP_NumericControl );
    HP_Sizer.Add( HPReset_ToolButton );

    //

    FineAdjust_Label.SetText( "Adjust parameter: " );
    FineAdjust_Label.SetMinWidth( labelWidth1 );
    FineAdjust_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    
    FineAdjust_ComboBox.AddItem( "SP" );
    FineAdjust_ComboBox.AddItem( "LP" );
    FineAdjust_ComboBox.AddItem( "HP" );
    FineAdjust_ComboBox.SetToolTip( "<p>Specifies the parameter to be fine adjusted.</p>" );
    FineAdjust_ComboBox.OnItemSelected( (ComboBox::item_event_handler)&GHSInterface::__ItemSelected, w );

    FineAdjust_Slider.SetScaledMinWidth( 250 );
    FineAdjust_Slider.SetRange( 0, 1000 ); // Do not change without examining the chain of consequences i.e. UpdateSliderControls and __SliderValueUpdated
    FineAdjust_Slider.SetStepSize( 1 );
    FineAdjust_Slider.SetValue( 500 );
    FineAdjust_Slider.SetToolTip( "<p>Fine adjustment slider with re-centering. "
      "The slider re-centers after being used so the next adjustment can be up or down.</p>" );
    FineAdjust_Slider.OnValueUpdated( (Slider::value_event_handler)&GHSInterface::__SliderValueUpdated, w );
    FineAdjust_Slider.OnEnter( (Control::event_handler)&GHSInterface::__FineAdjustSliderEnter, w ); // This allows us to re-centre the slider
    FineAdjust_Slider.OnLeave( (Control::event_handler)&GHSInterface::__FineAdjustSliderLeave, w ); // This allows us to re-centre the slider
    
    //FineAdjust_Sizer.AddUnscaledSpacing( labelWidth1 + w.LogicalPixelsToPhysical( 4 ) );
    FineAdjust_Sizer.SetSpacing( 3 );
    FineAdjust_Sizer.Add( FineAdjust_Label );
    FineAdjust_Sizer.Add( FineAdjust_ComboBox );
    FineAdjust_Sizer.Add( FineAdjust_Slider );
    
    FineAdjust_CheckBox.SetText( "Use highest sensitivity" );
    FineAdjust_CheckBox.SetToolTip( "<p>Check to adjust the 6th decimal place, leave unchecked to adjust the 4th decimal place.</p>" );
    FineAdjust_CheckBox.OnClick( (pcl::Button::click_event_handler)&GHSInterface::__ItemClicked, w );

    FALevel_Sizer.AddUnscaledSpacing( labelWidth1 + w.LogicalPixelsToPhysical( 4 ) );
    FALevel_Sizer.Add( FineAdjust_CheckBox );
    FALevel_Sizer.AddStretch();
    
    FAGroupSizer.SetSpacing( 3 );
    FAGroupSizer.Add( FineAdjust_Sizer );
    FAGroupSizer.Add( FALevel_Sizer );
    
    FineAdjust_GroupBox.SetTitle("Fine adjustment");
    FineAdjust_GroupBox.SetSizer( FAGroupSizer );

    //

    GHSParameters_Sizer.SetSpacing(3);
    GHSParameters_Sizer.Add( Inv_Sizer );
    GHSParameters_Sizer.Add( D_Sizer );
    GHSParameters_Sizer.Add( b_Sizer );
    GHSParameters_Sizer.AddSpacing(3);
    GHSParameters_Sizer.Add( SP_Sizer );
    GHSParameters_Sizer.Add( LP_Sizer );
    GHSParameters_Sizer.Add( HP_Sizer );
    GHSParameters_Sizer.AddSpacing(3);
    GHSParameters_Sizer.Add( FineAdjust_GroupBox );
    GHSParameters_Control.SetSizer( GHSParameters_Sizer );

    // ------------------------------------------------------------------
    // Linear controls
    
    BP_NumericControl.label.SetText( "Blackpoint (BP):" );
    BP_NumericControl.label.SetFixedWidth( labelWidth1 );
    BP_NumericControl.slider.SetScaledMinWidth( 250 );
    BP_NumericControl.slider.SetRange( 1, 1000 );
    BP_NumericControl.SetReal();
    BP_NumericControl.SetRange( TheGHSBPParameter->MinimumValue(), TheGHSBPParameter->MaximumValue() );
    BP_NumericControl.SetPrecision( TheGHSBPParameter->Precision() );
    BP_NumericControl.SetToolTip( "<p>Sets the Blackpoint for a linear stretch of the image. Note that any pixel with values less than the Blackpoint input will be clipped and the data lost. Contrast gained by performing the linear stretch will be evenly distributed over the image, which will be dimmed.  Pixels with values less than the Blackpoint will appear black and have 0 value. Updating this parameter will automatically update the low clipping proportion. Setting this parameter to a value less than 0 will extend the dynamic range at the low end.</p>" );
    BP_NumericControl.edit.SetFixedWidth( editWidth1 );
    BP_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );
    
    BPReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    BPReset_ToolButton.SetScaledFixedSize( 20, 20 );
    BPReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    BPReset_ToolButton.SetToolTip( "Reset Blackpoint parameter" );
    BPReset_ToolButton.SetCheckable(false);
    BPReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    BP_Sizer.SetSpacing( 4 );
    BP_Sizer.Add( BP_NumericControl );
    BP_Sizer.Add( BPReset_ToolButton );
    
    //
    
    LCP_NumericControl.label.SetText( "Low clip (LCP):" );
    LCP_NumericControl.label.SetFixedWidth( labelWidth1 );
    LCP_NumericControl.slider.SetScaledMinWidth( 250 );
    LCP_NumericControl.slider.SetRange( 1, 1000 );
    LCP_NumericControl.SetReal();
    LCP_NumericControl.SetRange( 0.0, 1.0 );
    LCP_NumericControl.SetPrecision( TheGHSBPParameter->Precision() );
    LCP_NumericControl.SetToolTip( "<p>Sets the clipping level for linear stretch of the image. Updating this parameter will automatically update the Blackpoint parameter in such a way that LCP is the maximum fraction of pixels clipped in any channel and set to zero in the linear stretch.  Note that any pixel with values less than the Blackpoint input will be clipped and the data lost.  Contrast gained by performing the linear stretch will be evenly distributed over the image, which will be dimmed.  Pixels with values less than the Blackpoint will appear black and have 0 value.</p>" );
    LCP_NumericControl.edit.SetFixedWidth( editWidth1 );
    LCP_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );

    LCPReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    LCPReset_ToolButton.SetScaledFixedSize( 20, 20 );
    LCPReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    LCPReset_ToolButton.SetToolTip( "Set Blackpoint parameter to greatest value that results in no clipping" );
    LCPReset_ToolButton.SetCheckable(false);
    LCPReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    LCP_Sizer.SetSpacing( 4 );
    LCP_Sizer.Add( LCP_NumericControl );
    LCP_Sizer.Add( LCPReset_ToolButton );
    
    WP_NumericControl.label.SetText( "Whitepoint (WP):" );
    WP_NumericControl.label.SetFixedWidth( labelWidth1 );
    WP_NumericControl.slider.SetScaledMinWidth( 250 );
    WP_NumericControl.slider.SetRange( 1, 1000 );
    WP_NumericControl.SetReal();
    WP_NumericControl.SetRange( TheGHSWPParameter->MinimumValue(), TheGHSWPParameter->MaximumValue() );
    WP_NumericControl.SetPrecision( TheGHSWPParameter->Precision() );
    WP_NumericControl.SetToolTip( "<p>Sets the Whitepoint for a linear stretch of the image. Note that any pixel with value greater than the Whitepoint input will be clipped and the data lost. Contrast gained by performing the linear stretch will be evenly distributed over the image, which will be brightened.  Pixels with values greater than the Whitepoint will appear white and have a value of 1.0. Updating this parameter will automatically update the high clipping proportion. Setting this parameter to a value greater than 1 will extend the dynamic range at the high end.</p>" );
    WP_NumericControl.edit.SetFixedWidth( editWidth1 );
    WP_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );
    
    WPReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    WPReset_ToolButton.SetScaledFixedSize( 20, 20 );
    WPReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    WPReset_ToolButton.SetToolTip( "Reset Whitepoint parameter" );
    WPReset_ToolButton.SetCheckable(false);
    WPReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    WP_Sizer.SetSpacing( 4 );
    WP_Sizer.Add( WP_NumericControl );
    WP_Sizer.Add( WPReset_ToolButton );
    
    HCP_NumericControl.label.SetText( "High clip (HCP):" );
    HCP_NumericControl.label.SetFixedWidth( labelWidth1 );
    HCP_NumericControl.slider.SetScaledMinWidth( 250 );
    HCP_NumericControl.slider.SetRange( 1, 1000 );
    HCP_NumericControl.SetReal();
    HCP_NumericControl.SetRange( 0.0, 1.0 );
    HCP_NumericControl.SetPrecision( TheGHSWPParameter->Precision() );
    HCP_NumericControl.SetToolTip( "<p>Sets the upper clipping level for linear stretch of the image. Updating this parameter will automatically update the Whitepoint parameter in such a way that HCP is the maximum fraction of pixels clipped in any channel and set to 1.0 in the linear stretch. Note that any pixel with value greater than the white-point input will be clipped and the data lost. Contrast gained by performing the linear stretch will be evenly distributed over the image, which will be brightened.  Pixels with values greater than the Whitepoint will appear white and have a value of 1.0.</p>" );
    HCP_NumericControl.edit.SetFixedWidth( editWidth1 );
    HCP_NumericControl.OnValueUpdated( (NumericEdit::value_event_handler)&GHSInterface::__RealValueUpdated, w );

    HCPReset_ToolButton.SetIcon( w.ScaledResource( ":/icons/clear.png" ) );
    HCPReset_ToolButton.SetScaledFixedSize( 20, 20 );
    HCPReset_ToolButton.SetFocusStyle( FocusStyle::NoFocus );
    HCPReset_ToolButton.SetToolTip( "Set Whitepoint parameter to lowest value that results in no clipping" );
    HCPReset_ToolButton.SetCheckable(false);
    HCPReset_ToolButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Reset_ButtonClick, w );
    
    HCP_Sizer.SetSpacing( 4 );
    HCP_Sizer.Add( HCP_NumericControl );
    HCP_Sizer.Add( HCPReset_ToolButton );

    LinearFineAdjust_Label.SetText( "Adjust parameter: " );
    LinearFineAdjust_Label.SetMinWidth( labelWidth1 );
    LinearFineAdjust_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    
    LinearFineAdjust_ComboBox.AddItem( "BP" );
    LinearFineAdjust_ComboBox.AddItem( "LCP" );
    LinearFineAdjust_ComboBox.AddItem( "WP" );
    LinearFineAdjust_ComboBox.AddItem( "HCP" );
    LinearFineAdjust_ComboBox.SetToolTip( "<p>Specifies the parameter to be fine adjusted.</p>" );
    LinearFineAdjust_ComboBox.OnItemSelected( (ComboBox::item_event_handler)&GHSInterface::__ItemSelected, w );
    
    LinearFineAdjust_Slider.SetScaledMinWidth( 250 );
    LinearFineAdjust_Slider.SetRange( 0, 1000 ); // Do not change without examining the chain of consequences i.e. UpdateSliderControls and __SliderValueUpdated
    LinearFineAdjust_Slider.SetStepSize( 1 );
    LinearFineAdjust_Slider.SetValue( 500 );
    LinearFineAdjust_Slider.SetToolTip( "<p>Symmetry point fine adjustment slider with re-centering. "
      "The slider re-centers after being used so the next adjustment can be up or down.</p>" );
    LinearFineAdjust_Slider.OnValueUpdated( (Slider::value_event_handler)&GHSInterface::__SliderValueUpdated, w );
    LinearFineAdjust_Slider.OnEnter( (Control::event_handler)&GHSInterface::__FineAdjustSliderEnter, w ); // This allows us to re-centre the slider
    LinearFineAdjust_Slider.OnLeave( (Control::event_handler)&GHSInterface::__FineAdjustSliderLeave, w ); // This allows us to re-centre the slider
    
    //FineAdjust_Sizer.AddUnscaledSpacing( labelWidth1 + w.LogicalPixelsToPhysical( 4 ) );
    LinearFineAdjust_Sizer.SetSpacing( 3 );
    LinearFineAdjust_Sizer.Add( LinearFineAdjust_Label );
    LinearFineAdjust_Sizer.Add( LinearFineAdjust_ComboBox );
    LinearFineAdjust_Sizer.Add( LinearFineAdjust_Slider );
    
    LinearFineAdjust_CheckBox.SetText( "Use highest sensitivity" );
    LinearFineAdjust_CheckBox.SetToolTip( "<p>Check to adjust the 6th decimal place, leave unchecked to adjust the 4th decimal place.</p>" );
    LinearFineAdjust_CheckBox.OnClick( (pcl::Button::click_event_handler)&GHSInterface::__ItemClicked, w );

    LinearFALevel_Sizer.AddUnscaledSpacing( labelWidth1 + w.LogicalPixelsToPhysical( 4 ) );
    LinearFALevel_Sizer.Add( LinearFineAdjust_CheckBox );
    LinearFALevel_Sizer.AddStretch();
    
    LinearFAGroupSizer.SetSpacing( 3 );
    LinearFAGroupSizer.Add( LinearFineAdjust_Sizer );
    LinearFAGroupSizer.Add( LinearFALevel_Sizer );
    
    LinearFineAdjust_GroupBox.SetTitle("Fine adjustment");
    LinearFineAdjust_GroupBox.SetSizer( LinearFAGroupSizer );

    LinearParameters_Sizer.SetSpacing( 3 );
    LinearParameters_Sizer.Add( BP_Sizer );
    LinearParameters_Sizer.Add( LCP_Sizer );
    LinearParameters_Sizer.AddSpacing( 3 );
    LinearParameters_Sizer.Add( WP_Sizer );
    LinearParameters_Sizer.Add( HCP_Sizer );
    LinearParameters_Sizer.AddSpacing( 3 );
    LinearParameters_Sizer.Add( LinearFineAdjust_GroupBox );
    LinearParameters_Control.SetSizer( LinearParameters_Sizer );
    
    GHSParameters_Control.SetVisible(true);
    LinearParameters_Control.SetVisible(false);

    Parameters_Sizer.Add( ST_Sizer );
    Parameters_Sizer.AddSpacing( 6 );
    Parameters_Sizer.Add( GHSParameters_Control );
    Parameters_Sizer.Add( LinearParameters_Control );
    Parameters_Control.SetSizer( Parameters_Sizer );
    
    GHSParameters_SectionBar.SetTitle( "Transformation" );
    GHSParameters_SectionBar.SetSection( Parameters_Control );
    
    // ------------------------------------------------------------------
    // Readout controls
    
    ReadoutTitle_Label.SetText( "Value:" );
    ReadoutTitle_Label.SetMinWidth( labelWidth1 );
    ReadoutTitle_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    
    Readout_Label.SetTextAlignment( TextAlign::Left|TextAlign::VertCenter );
    Readout_Label.SetStyle( FrameStyle::Sunken );
    Readout_Label.SetLineWidth( 1 );
    
    ReadoutValue_Sizer.SetSpacing( 4 );
    ReadoutValue_Sizer.Add( ReadoutTitle_Label );
    ReadoutValue_Sizer.Add( Readout_Label );
    ReadoutValue_Sizer.AddStretch();

    //
    
    ReadoutSourceTitle_Label.SetText( "Source:" );
    ReadoutSourceTitle_Label.SetMinWidth( labelWidth1 );
    ReadoutSourceTitle_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    
    ReadoutSource_Label.SetTextAlignment( TextAlign::Left|TextAlign::VertCenter );
    ReadoutSource_Label.SetStyle( FrameStyle::Sunken );
    ReadoutSource_Label.SetLineWidth( 1 );
    
    ReadoutSource_Sizer.SetSpacing( 4 );
    ReadoutSource_Sizer.Add( ReadoutSourceTitle_Label );
    ReadoutSource_Sizer.Add( ReadoutSource_Label );
    ReadoutSource_Sizer.AddStretch();
    
    //
    
    ReadoutDescTitle_Label.SetText( "Description:" );
    ReadoutDescTitle_Label.SetMinWidth( labelWidth1 );
    ReadoutDescTitle_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    
    ReadoutDesc_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    ReadoutDesc_Label.SetStyle( FrameStyle::Sunken );
    ReadoutDesc_Label.SetLineWidth( 1 );
    ReadoutDesc_Label.SetToolTip("<p>Set the readout statistic by adjusting the Readout/Calculation Mode in the PixInsight Readout Options. The channel(s) reflect the currently selected colour mode in the Non-linear Parameters panel.</p>");
    
    ReadoutSend_PushButton.SetText( "Send to SP" );
    ReadoutSend_PushButton.SetFocusStyle( FocusStyle::NoFocus );
    ReadoutSend_PushButton.SetToolTip( "Use readout value to set the value of SP." );
    ReadoutSend_PushButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Readout_ButtonClick, w );
    
    ReadoutDesc_Sizer.SetSpacing( 4 );
    ReadoutDesc_Sizer.Add( ReadoutDescTitle_Label );
    ReadoutDesc_Sizer.Add( ReadoutDesc_Label );
    ReadoutDesc_Sizer.AddStretch();
    ReadoutDesc_Sizer.Add( ReadoutSend_PushButton );
    
    //
    
    ReadoutAreaTitle_Label.SetText( "Area: " );
    ReadoutAreaTitle_Label.SetMinWidth( labelWidth1 );
    ReadoutAreaTitle_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    
    ReadoutArea_Label.SetTextAlignment( TextAlign::Right|TextAlign::VertCenter );
    ReadoutArea_Label.SetStyle( FrameStyle::Sunken );
    ReadoutArea_Label.SetLineWidth( 1 );
    ReadoutArea_Label.SetToolTip("<p>Set the readout area size by adjusting the probe size in the PixInsight Readout Options. You may also wish to adjust the preview size and preview zoom in the PixInsight Readout Options to set the preview callout to your own preference.</p>");
    
    ReadoutClear_PushButton.SetText( "Clear" );
    ReadoutClear_PushButton.SetFocusStyle( FocusStyle::NoFocus );
    ReadoutClear_PushButton.SetToolTip( "Clear readout." );
    ReadoutClear_PushButton.OnClick( (ToolButton::click_event_handler)&GHSInterface::__Readout_ButtonClick, w );
    
    ReadoutArea_Sizer.SetSpacing( 4 );
    ReadoutArea_Sizer.Add( ReadoutAreaTitle_Label );
    //ReadoutSize_Sizer.Add( ReadoutSize_SpinBox );
    ReadoutArea_Sizer.Add( ReadoutArea_Label );
    ReadoutArea_Sizer.AddStretch();
    ReadoutArea_Sizer.Add( ReadoutClear_PushButton );
    
    //
    
    Readout_Sizer.SetSpacing( 3 );
    Readout_Sizer.Add( ReadoutValue_Sizer );
    Readout_Sizer.Add( ReadoutSource_Sizer );
    Readout_Sizer.Add( ReadoutDesc_Sizer );
    Readout_Sizer.Add( ReadoutArea_Sizer );
    
    Readout_Control.SetSizer( Readout_Sizer );
    
    Readout_SectionBar.SetTitle( "Readout Data" );
    Readout_SectionBar.SetSection( Readout_Control );
    
    // ------------------------------------------------------------------
    // Web view
    
    //GHSWebView.SetMinSize(600,400);
    //GHSWebView.SetVisible(false);
    //GHSWebView.SetPlainText(String("loading..."));
    //GHSWebView.OnClose( (Control::close_event_handler)&GHSInterface::__Web_Close, w );
    
    // ------------------------------------------------------------------
    // Global layout
    
    ColourOptions_Control.SetVisible( false );
    
    //
    
    Global_Sizer.SetMargin( 8 );
    Global_Sizer.SetSpacing( 3 );
    Global_Sizer.Add( InputHistogram_SectionBar );
    Global_Sizer.Add( InputHistogram_Control );
    Global_Sizer.Add( Readout_SectionBar );
    Global_Sizer.Add( Readout_Control );
    Global_Sizer.Add( ColourOptions_SectionBar );
    Global_Sizer.Add( ColourOptions_Control );
    Global_Sizer.Add( GHSParameters_SectionBar );
    Global_Sizer.Add( Parameters_Control );
    Global_Sizer.AddSpacing(6);
    

    w.SetSizer( Global_Sizer );

    w.EnsureLayoutUpdated();
    w.AdjustToContents();
    
    // ------------------------------------------------------------------
    // Real time preview timer
    
    UpdateRealTimePreview_Timer.SetSingleShot();
    UpdateRealTimePreview_Timer.SetInterval( 0.025 );
    UpdateRealTimePreview_Timer.OnTimer( (Timer::timer_event_handler)&GHSInterface::__UpdateRealTimePreview_Timer, w );
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------

} // pcl

// ----------------------------------------------------------------------------
// EOF GHSInterface.cpp
