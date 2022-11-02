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
// GHSInterface.h
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

#ifndef __GHSInterface_h
#define __GHSInterface_h

#include <pcl/CheckBox.h>
#include <pcl/ComboBox.h>
#include <pcl/Edit.h>
#include <pcl/Frame.h>
#include <pcl/GroupBox.h>
#include <pcl/Label.h>
#include <pcl/NumericControl.h>
#include <pcl/ProcessInterface.h>
#include <pcl/PushButton.h>
#include <pcl/Sizer.h>
#include <pcl/SpinBox.h>
#include <pcl/Thread.h>
#include <pcl/Timer.h>
#include <pcl/Console.h>
#include <pcl/SectionBar.h>
#include <pcl/ScrollBox.h>
#include <pcl/ToolButton.h>
#include <pcl/View.h>
#include <pcl/ViewList.h>
#include <pcl/WebView.h>
#include <pcl/Dialog.h>
#include <pcl/ReadoutOptions.h>

#include "GHSInstance.h"

namespace pcl
{

// ----------------------------------------------------------------------------

class PCL_CLASS Graphics;

class GHSInterface : public ProcessInterface
{
public:

   GHSInterface();
   virtual ~GHSInterface();
    
    enum graph_style   { LineStyle, AreaStyle, BarStyle, DotStyle, NumberOfGraphStyles };
    

   IsoString Id() const override;
   MetaProcess* Process() const override;
   String IconImageSVGFile() const override;
    
    InterfaceFeatures Features() const override;
    bool WantsRealTimePreviewNotifications() const override;
    void RealTimePreviewOwnerChanged(ProcessInterface& iface) override;
    void RealTimePreviewUpdated(bool active) override;
    
   void ApplyInstance() const override;
   void ResetInstance() override;
   bool Launch( const MetaProcess&, const ProcessImplementation*, bool& dynamic, unsigned& /*flags*/ ) override;
   ProcessImplementation* NewProcess() const override;
   bool ValidateProcess( const ProcessImplementation&, String& whyNot ) const override;
   bool RequiresInstanceValidation() const override;
   bool ImportProcess( const ProcessImplementation& ) override;
    
    bool RequiresRealTimePreviewUpdate(const UInt16Image&, const View&, const Rect&, int zoomLevel) const override;
    bool GenerateRealTimePreview(UInt16Image&, const View&, const Rect&, int zoomLevel, String& info) const override;
    
    bool WantsImageNotifications() const override;
    void ImageUpdated( const View& v ) override;
    void ImageFocused( const View& v ) override;
    void ImageDeleted( const View& v ) override;
    bool WantsViewPropertyNotifications() const override;
    void ViewPropertyUpdated( const View& v, const IsoString& property ) override;
    void ViewPropertyDeleted( const View& v, const IsoString& property ) override;
    bool WantsReadoutNotifications() const override;
    void BeginReadout( const View& v ) override;
    void UpdateReadout( const View& v, const DPoint& p, double R, double G, double B, double A ) override;
    void EndReadout( const View& v ) override;
    bool WantsGlobalNotifications() const override;
    void ReadoutOptionsUpdated() override;
    void SaveSettings() const override;
    void LoadSettings() override;

    int PlotResolution() const
    {
       return m_plotResolution;
    }

    bool RejectingSaturated() const
    {
       return m_rejectSaturated;
    }

private:

   /*
    * The instance being defined
    */
   GHSInstance m_instance;
    
    
    
    class RealTimeThread : public Thread
    {
    public:

        UInt16Image m_image;

        RealTimeThread();

        void Reset( const UInt16Image&, const GHSInstance& );

        void Run() override;

    private:

        GHSInstance m_instance;
    };

    /*
     * Child controls
     */
    struct GUIData
    {
        GUIData( GHSInterface& );
        VerticalSizer       Global_Sizer;
       
        HorizontalSizer     Graphics_Sizer;
            ToolButton          RejectSaturated_ToolButton;
            ToolButton          ShowCurve_ToolButton;
            ToolButton          ShowGrid_ToolButton;
            ToolButton          LogHistogram_ToolButton;
            SpinBox             HorizontalZoom_SpinBox;
            ToolButton          Zoom11_ToolButton;
            ToolButton          ZoomIn_ToolButton;
            ToolButton          ZoomOut_ToolButton;
            ToolButton          PreStretchHist_ToolButton;
            //ToolButton          WebView_ToolButton;


        SectionBar          InputHistogram_SectionBar;
        Control             InputHistogram_Control;
        VerticalSizer       InputHistogram_Sizer;
        ScrollBox           InputHistogram_ScrollBox;
        VerticalSizer       InputHistogramViewport_Sizer;
            Control             InputHistogramPlot_Control;
            Control             HistogramSliders_Control;
        
        
        SectionBar           GHSParameters_SectionBar;
        Control              GHSParameters_Control;
        VerticalSizer        GHSParameters_Sizer;
            HorizontalSizer         ST_Sizer;
                Label                   ST_Label;
                ComboBox                ST_ComboBox;
            HorizontalSizer         Inv_Sizer;
                CheckBox                Inv_CheckBox;
            HorizontalSizer         D_Sizer;
                NumericControl          D_NumericControl;
                ToolButton              DReset_ToolButton;
            HorizontalSizer         b_Sizer;
                NumericControl          b_NumericControl;
                ToolButton               bReset_ToolButton;
            GroupBox                FineAdjust_GroupBox;
            HorizontalSizer         FineAdjust_Sizer;
                Label                   FineAdjust_Label;
                ComboBox                FineAdjust_ComboBox;
                Slider                  FineAdjust_Slider;
            HorizontalSizer         FALevel_Sizer;
                CheckBox                FineAdjust_CheckBox;
            VerticalSizer           FAGroupSizer;
            HorizontalSizer         SP_Sizer;
                NumericControl          SP_NumericControl;
                ToolButton              SPReset_ToolButton;
            HorizontalSizer         LP_Sizer;
                NumericControl          LP_NumericControl;
                ToolButton              LPReset_ToolButton;
            HorizontalSizer         HP_Sizer;
                NumericControl          HP_NumericControl;
                ToolButton              HPReset_ToolButton;
        
        SectionBar           ColourOptions_SectionBar;
        Control              ColourOptions_Control;
        GroupBox             ColourOptions_GroupBox;
        VerticalSizer        ColourOptions_Sizer;
        VerticalSizer        ColourOptionsGB_Sizer;
            HorizontalSizer         SC_Sizer;
                Label                   SC_Label;
                ComboBox                SC_ComboBox;
            HorizontalSizer         CT_Sizer;
                Label                   CT_Label;
                ComboBox                CT_ComboBox;
            HorizontalSizer         CB_Sizer;
                NumericControl          CB_NumericControl;
                ToolButton              CBReset_ToolButton;
            HorizontalSizer         RGBWS_Sizer;
                CheckBox                RGBWS_CheckBox;
        
        SectionBar           LinearParameters_SectionBar;
        Control              LinearParameters_Control;
        VerticalSizer        LinearParameters_Sizer;
            HorizontalSizer         BP_Sizer;
                NumericControl          BP_NumericControl;
                ToolButton              BPReset_ToolButton;
            HorizontalSizer         LCP_Sizer;
                NumericControl          LCP_NumericControl;
                ToolButton              LCPReset_ToolButton;
            HorizontalSizer         WP_Sizer;
                NumericControl          WP_NumericControl;
                ToolButton              WPReset_ToolButton;
            HorizontalSizer         HCP_Sizer;
                NumericControl          HCP_NumericControl;
                ToolButton              HCPReset_ToolButton;
            
            GroupBox                LinearFineAdjust_GroupBox;
            HorizontalSizer         LinearFineAdjust_Sizer;
                Label                   LinearFineAdjust_Label;
                ComboBox                LinearFineAdjust_ComboBox;
                Slider                  LinearFineAdjust_Slider;
            HorizontalSizer         LinearFALevel_Sizer;
                CheckBox                LinearFineAdjust_CheckBox;
            VerticalSizer           LinearFAGroupSizer;
        
        VerticalSizer       Parameters_Sizer;
        Control             Parameters_Control;
        
        SectionBar           Readout_SectionBar;
        Control              Readout_Control;
        VerticalSizer        Readout_Sizer;
            HorizontalSizer         ReadoutSource_Sizer;
                Label                   ReadoutSourceTitle_Label;
                Label                   ReadoutSource_Label;
            HorizontalSizer         ReadoutDesc_Sizer;
                Label                   ReadoutDescTitle_Label;
                Label                   ReadoutDesc_Label;
            HorizontalSizer         ReadoutValue_Sizer;
                Label                   ReadoutTitle_Label;
                Label                   Readout_Label;
            HorizontalSizer         ReadoutArea_Sizer;
                Label                   ReadoutArea_Label;
                Label                   ReadoutAreaTitle_Label;
            PushButton              ReadoutClear_PushButton;
            PushButton              ReadoutSend_PushButton;
        
        //WebView             GHSWebView;
        
        Timer UpdateRealTimePreview_Timer;
    };

    GUIData* GUI = nullptr;
    mutable AutoPointer<RealTimeThread> m_realTimeThread;
    
    /*
     * Workbench
     */
    
    enum slider_id          { C0Slider, MSlider, C1Slider, NoSlider = -1 };
    enum cursor_status      { NoCursor, InputCursor, OutputCursor };
    enum readout_source     { FromPreview, FromHistogram, FromNone };
    enum fa_parameter       { SP, LP, HP, BP, LCP, WP, HCP };

    typedef Histogram::count_type       count_type;
    typedef Histogram::histogram_type   histogram_type;
    typedef Array<Histogram>            histogram_list;
    
    // Active view
    View            m_currentView;

    
    // Histogram data.
    histogram_list m_sourceData;        // source input histograms, 16-bit resolution
    histogram_list m_inputData;         // input histograms, rescaled to the plot resolution
    histogram_list m_outputData;        // output RGBA histograms (R,G,B include the combined RGB/K transformation)

    // Histogram plot resolution, or the number of discrete histogram levels to
    // be represented.
    int            m_plotResolution         = 65536;

    // Current histogram channel.
    // 0=R 1=G 2=B 3=RGB/K 4=Alpha
    int            m_channel                = 3;

    // Current graph style.
    graph_style    m_inputGraphStyle        = LineStyle;
    graph_style    m_outputGraphStyle       = AreaStyle;
    bool           m_logGraph               = false;
    bool           m_showHist               = true;

    // Histogram clippings.
    int m_shadowsCount                      = 0;
    int m_highlightsCount                   = 0;

    // Image readouts.
    bool           m_readoutActive          = false;
    DVector        m_inputReadouts          = DVector( 0.0, 5 ); // 0=R 1=G 2=B 3=notUsed 4=Alpha
    DVector        m_outputReadouts         = DVector( 0.0, 5 );
    DPoint                              m_readoutPoint;
    readout_source                      m_readoutSource         = FromNone;
    DRect                               m_readoutRect;
    ReadoutOptions::readout_mode        m_readoutStatistic      = ReadoutOptions::readout_mode::Mean;
    double                              m_readoutValue;
    double                              m_readoutSize           = 15.0;
    
    // Graph amplification factors.
    int            m_inputZoomX             = 1;
    int            m_inputZoomY             = 1;

    // Histogram representation options.
    bool           m_rejectSaturated        = true; // ignore the first and last histogram counts to compute peaks
    bool           m_showMTF                = true; // draw the midtones transfer function curve
    bool           m_showGrid               = true; // draw coordinate grids

    // Interactive states.
    slider_id      m_sliderBeingDragged     = NoSlider; // moving one of our little triangular things?

    // Graph cursor.
    cursor_status  m_cursorStatus           = NoCursor;
    Point          m_cursorPos              = -1; // cursor position in viewport crds.
    DPoint         m_histogramPos           = 0;  // cursor position in normalized crds.

    // Screen bitmap, input histogram viewport.
    Bitmap         m_inputBitmap            = Bitmap::Null();
    bool           m_inputDirty             = true;

    // Screen bitmap, slider area.
    Bitmap         m_slidersBitmap          = Bitmap::Null();
    bool           m_slidersDirty           = true;

    // Graph colours
    RGBA            m_gridColour0           = RGBAColor( 0x60, 0x60, 0x60 );
    RGBA            m_gridColour1           = RGBAColor( 0x20, 0x20, 0x20 );
    RGBA            m_backgroundColour      = RGBAColor( 0x20, 0x20, 0x20 );
    RGBA            m_cursorColour          = RGBAColor( 0xff, 0xff, 0xff );
    RGBA            m_identityColour        = RGBAColor( 0xff, 0xff, 0xff );
    RGBA            m_readoutColour         = RGBAColor( 0xff, 0xcc, 0x00 );
    RGBA            m_transformColour       = RGBAColor( 0xff, 0x20, 0x20 );

    // Minimum graph dimensions.
    int            m_minHistogramWidth      = 400;
    int            m_minHistogramHeight     = 200;
    int            m_sliderControlSize      = 12;

    // Flag true during viewport transitional states (e.g. resizing).
    bool           m_settingUp              = false;

    // Images
    ImageVariant    m_lightnessImage;
    ImageVariant    m_saturationImage;
    ImageVariant    m_luminanceImage;
    
    // Control update variables
    pcl_enum        m_STBeforeLinear;
    pcl_enum        m_lastST                = GHSST::ST_GeneralisedHyperbolic;
    bool            m_useFinestAdjustment   = true;
    
    double          m_FAPBeforeAdjustment;
    fa_parameter    m_linearFAParameter     = BP;
    fa_parameter    m_nonlinFAParameter     = SP;
    
    
    
    

    /*
     * Main calculation routines
     */
    bool GetSourceHistograms();
    void CalculateInputHistograms();
    void CalculateOutputHistograms();
    void DestroyOutputHistograms();
    void CalculateClippingCounts();
    double CalculateLowClipPoint(double);
    double CalculateHighClipPoint(double);

    /*
     * Setting parameters
     */
    void SetShadowProtect( double );
    void SetHighlightProtect( double );
    void SetSymmetryPoint( double );
    void SetRejectSaturated( bool );
    void SetInputZoom( int, int, const Point* = nullptr );
    void SetBlackpointFromLCP( double );
    void SetWhitepointFromHCP( double );

    /*
     * GUI Updates
     */
    void UpdateControls();
    void UpdateZoomControls();
    void UpdateGraphicsControls();
    void UpdateReadoutControls();
    void UpdateTransformationControls();
    void UpdateClippingCountControls();
    void UpdateHistograms();
    void UpdateHistogramSliders();
    void UpdateHistogramInfo();
    void UpdateRealTimePreview();
    void SynchronizeWithCurrentView();

    /*
     * Screen bitmap regeneration
     */
    void RegenerateInputViewport();
    void RegenerateSlidersViewport();

    /*
     * Histogram drawing primitives
     */
    DVector GetHistogramPlotValues( Graphics&, const Rect&, const Histogram&, int, int, double&);
    void PlotMonoHistogramValues(Graphics&, const Rect&, DVector, double, int, int, int, graph_style, RGBA);
    void PlotColourHistogramValues(Graphics&, const Rect&, DVector, DVector, DVector, double, int, int, int, graph_style, int);
    void PlotGrid( Graphics&, const Rect& viewport, int width, int height, int hZoom, int vZoom );
    void PlotScale( Graphics&, const Rect& viewport, int width );
    void PlotHandler( Graphics&, double value, int x0, int width );
    void PlotGHSCurve( Graphics& g, const Rect& viewport, int width, int height );
    void PlotIdentityTransformation( Graphics&, const Bitmap&, const Rect& viewport, int width, int height );
    void PlotReadouts( Graphics&, const Bitmap&, const Rect& viewport, const DVector&, int width, int height );
    void PlotInputSelection( Graphics&, const Bitmap&, const Rect& viewport, const DVector&, int width, int height );
    void PlotCursor( Graphics&, const Rect& viewport );

    /*
     * Miscellaneous drawing and GUI helpers
     */
    RGBA HandlerColor( double ) const;
    slider_id FindHandler( double ) const;
    double SliderToHistogram( int ) const;

    /*
    * Event handlers
    */
    void __RealValueUpdated( NumericEdit& sender, double value );
    void __ItemClicked( Button& sender, bool checked );
    void __ItemSelected( ComboBox& sender, int itemIndex );
    void __SliderValueUpdated( Slider& sender, int value );
    void __FineAdjustSliderEnter( Control& sender );
    void __FineAdjustSliderLeave( Control& sender );
    void __FineAdjustSliderGetFocus( Control& sender );
    void __FineAdjustSliderLoseFocus( Control& sender );
    void __FineAdjustSliderMousePress( Control& sender, const pcl::Point &pos, int button, unsigned buttons, unsigned modifiers );
    void __FineAdjustSliderMouseRelease( Control& sender, const pcl::Point &pos, int button, unsigned buttons, unsigned modifiers );
    
    void __Histogram_Paint( Control& sender, const pcl::Rect& updateRect );
    void __Sliders_Paint( Control& sender, const pcl::Rect& updateRect );
    void __Histogram_Resize( Control& sender, int newWidth, int newHeight, int oldWidth, int oldHeight );
    void __Histogram_ScrollPosUpdated( ScrollBox& sender, int pos );
    void __Histogram_Enter( Control& sender );
    void __Histogram_Leave( Control& sender );
    void __Histogram_MousePress( Control& sender, const pcl::Point& pos, int button, unsigned buttons, unsigned modifiers );
    void __Histogram_MouseMove( Control& sender, const pcl::Point& pos, unsigned buttons, unsigned modifiers );
    void __Sliders_MousePress( Control& sender, const pcl::Point& pos, int button, unsigned buttons, unsigned modifiers );
    void __Sliders_MouseRelease( Control& sender, const pcl::Point& pos, int button, unsigned buttons, unsigned modifiers );
    void __Sliders_MouseMove( Control& sender, const pcl::Point& pos, unsigned buttons, unsigned modifiers );
    void __Reset_ButtonClick( Button&, bool );
    void __AutoZero_ButtonClick( Button&, bool );
    void __Zoom_ButtonClick( Button&, bool );
    void __Readout_ButtonClick( Button&, bool );
    void __Zoom_ValueUpdated( SpinBox& sender, int value );
    void __RejectSaturated_ButtonClick( Button&, bool );
    void __ShowCurve_ButtonClick( Button&, bool );
    void __ShowGrid_ButtonClick( Button&, bool );
    void __LogHistogram_ButtonClick( Button&, bool );
    void __PreStretchHist_ButtonClick( Button&, bool );
    //void __Web_ButtonClick( Button&, bool );
    //void __Web_Close( Control&, bool& );
    void __UpdateRealTimePreview_Timer( Timer& );

    friend struct GUIData;
    friend class RealTimeThread;
};

// ----------------------------------------------------------------------------

PCL_BEGIN_LOCAL
extern GHSInterface* TheGHSInterface;
PCL_END_LOCAL

// ----------------------------------------------------------------------------

} // pcl

#endif   // __GHSInterface_h

// ----------------------------------------------------------------------------
// EOF GHSInterface.h
