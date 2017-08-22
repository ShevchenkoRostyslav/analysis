// Common style file for Hbb analysis
// --------------------------------------------------------------
// 
// Defines a set of functions to define the plotting style and
// provide commonly used style objects such as histogram title
// and legends. Can be used in compiled code and uncompiled scripts.
//
// Always call HbbStyle::set(PublicationStatus) in the beginning
// of your plotting script to adjust the gStyle options. Also,
// the behaviour of some other methods of this class depend on
// this, e.g. the histogram title displays "CMS preliminary" etc.
// depending on the PublicationStatus. Call HbbStyle::set before
// declaring any histograms (or call TH1::UseCurrentStyle()) to
// make sure the style settings are used.
// --------------------------------------------------------------

#ifndef HBB_STYLE_CC
#define HBB_STYLE_CC

#include "TColor.h"
#include "TError.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"


// Publication status: determines what is plotted in title
enum PublicationStatus { INTERNAL, INTERNAL_SIMULATION, PRELIMINARY, PUBLIC, SIMULATION, UNPUBLISHED, PRIVATE };
TString toTString(const PublicationStatus status) {
  TString str = "";
  if(      status == INTERNAL )            str = "internal";
  else if( status == INTERNAL_SIMULATION ) str = "simulation (internal)";
  else if( status == PRELIMINARY )         str = "preliminary";
  else if( status == PUBLIC      )         str = "public";
  else if( status == SIMULATION  )         str = "simulation (public)";
  else if( status == UNPUBLISHED )         str = "unpublished";
  else if( status == PRIVATE )			   str = "Private work";

  return str;
}

enum BkgTemplateType { B2B1C2bb, C1bb, Qbb, bbB2B1C2, bbC1Q };


class HbbStyle {
public:
  // Adjusts the gStyle settings and store the PublicationStatus
  static void set(const PublicationStatus status);
  static PublicationStatus status() { return publicationStatus_; }

  // Draws a title on the current pad
  //  <CMS label>,  #sqrt{s} = 8 TeV,  19.7fb^{-1}
  // where <CMS label> depends on the PublicationStatus
  //  INTERNAL    : no extra label (intended for AN-only plots with data)
  //  INTERNAL    : show "Simulation" label (intended for AN-only plots, no lumi, no "CMS")
  //  PRELIMINARY : show "CMS preliminary" label
  //  PUBLIC      : show "CMS" label
  //  SIMULATION  : show "CMS Simulation" label
  //  UNPUBLISHED : show "CMS (unpublished)" label (intended for additional material on TWiki)
  // Note that this method does not allow for easy memory
  // handling. For that, use standardTitle().
  static void drawStandardTitle() { standardTitle()->Draw("same"); }

  // Returns a TPaveText object that fits as a histogram title
  // with the current pad dimensions.
  // It has the same text as described in drawStandardTitle().
  // The idea of this method is that one has control over the
  // TPaveText object and can do proper memory handling.
  static TPaveText* standardTitle(const PublicationStatus status) {
    return title(header(status));
  }
  static TPaveText* standardTitle() {
    return standardTitle(publicationStatus_);
  }

  // Returns a TPaveText object that fits as a histogram title
  // with the current pad dimensions and displays the specified text txt.
  static TPaveText* customTitle(const TString& txt) { return title(txt); }

  // Returns a TLegend object that fits into the top-right corner
  // of the current pad. Its width, relative to the pad size (without
  // margins), can be specified. Its height is optimized for nEntries
  // entries.
  static TLegend* legend(const int nEntries, const double relWidth=0.5) {
    return legendTR(nEntries,relWidth);
  }
  static TLegend* legend(TString position, const int nEntries, const double relWidth=0.5) {
    position.ToLower();
    if( !( position.Contains("top") || position.Contains("bottom") ) )
      position += "top";
    if( !( position.Contains("left") || position.Contains("right") ) )
      position += "right";
    TLegend* leg = 0;
    if(        position.Contains("top")    && position.Contains("right") ) {
      leg = legendTR(nEntries,relWidth);
    } else if( position.Contains("top")    && position.Contains("left")  ) {
      leg = legendTL(nEntries,relWidth);
    } else if( position.Contains("bottom") && position.Contains("right") ) {
      leg = legendBR(nEntries,relWidth);
    } else if( position.Contains("bottom") && position.Contains("left")  ) {
      leg = legendBL(nEntries,relWidth);
    } else {
      leg = legendTR(nEntries,relWidth);
    }
    return leg;
  }
  // Same but explicitly state position on pad
  static TLegend* legendTL(const int nEntries, const double relWidth=0.5) {
    return legend(nEntries,relWidth,true,true);
  }
  static TLegend* legendTR(const int nEntries, const double relWidth=0.5) {
    return legend(nEntries,relWidth,false,true);
  }
  static TLegend* legendBL(const int nEntries, const double relWidth=0.5) {
    return legend(nEntries,relWidth,true,false);
  }
  static TLegend* legendBR(const int nEntries, const double relWidth=0.5) {
    return legend(nEntries,relWidth,false,false);
  }


  // Returns a TPaveText object that fits into the top-right corner
  // of the current pad and that can be used for additional labels.
  // Its width, relative to the pad size (without margins), can be
  // specified. Its height is optimized for nEntries entries.
  static TPaveText* label(const int nEntries, const double relWidth=0.5) {
    return labelTR(nEntries,relWidth);
  }

  static TPaveText* label(TString position, const int nEntries, const double relWidth=0.5) {
    position.ToLower();
    if( !( position.Contains("top") || position.Contains("bottom") ) )
      position += "top";
    if( !( position.Contains("left") || position.Contains("right") ) )
      position += "right";
    TPaveText* label = 0;
    if(        position.Contains("top")    && position.Contains("right") ) {
      label = labelTR(nEntries,relWidth);
    } else if( position.Contains("top")    && position.Contains("left")  ) {
      label = labelTL(nEntries,relWidth);
    } else if( position.Contains("bottom") && position.Contains("right") ) {
      label = labelBR(nEntries,relWidth);
    } else if( position.Contains("bottom") && position.Contains("left")  ) {
      label = labelBL(nEntries,relWidth);
    } else {
      label = labelTR(nEntries,relWidth);
    }

    return label;
  }

  // Same but explicitly state position on pad
  static TPaveText* labelTL(const int nEntries, const double relWidth=0.5) {
    return label(nEntries,relWidth,true,true);
  }
  static TPaveText* labelTR(const int nEntries, const double relWidth=0.5) {
    return label(nEntries,relWidth,false,true);
  }
  static TPaveText* labelBL(const int nEntries, const double relWidth=0.5) {
    return label(nEntries,relWidth,true,false);
  }
  static TPaveText* labelBR(const int nEntries, const double relWidth=0.5) {
    return label(nEntries,relWidth,false,false);
  }


  // Returns the nicely-formatted name of the background template.
  // Use in labels, legends, etc.
  static TString bkgTemplDisplayName(const BkgTemplateType type);

  // for axis titles
  static TString axisTitleMass() { return "M_{12} [GeV]"; }
  static TString axisTitleXTag() { return "X_{123}"; }


  static double lineHeight() { return lineHeight_; }


private:
  static PublicationStatus publicationStatus_;
  static double lineHeight_;
  static double margin_;
  static double textSize_;

  // creates a title
  static TPaveText* title(const TString& txt);

  // returns the standard-title (CMS label + sqrt{s} + L) depending
  // on the PublicationStatus 
  static TString header(const PublicationStatus status);

  // NDC coordinates for TPave, TLegend,...
  static void setXCoordinatesL(const double relWidth, double& x0, double& x1);
  static void setXCoordinatesR(const double relWidth, double& x0, double& x1);
  static void setYCoordinatesT(const int nEntries, double& y0, double& y1);
  static void setYCoordinatesB(const int nEntries, double& y0, double& y1);

  static TLegend* legend(const int nEntries, const double relWidth, const bool left, const bool top);
  static TPaveText* label(const int nEntries, const double relWidth, const bool leftt, const bool top);

};

PublicationStatus HbbStyle::publicationStatus_ = INTERNAL;
double HbbStyle::lineHeight_ = 0.042;
double HbbStyle::margin_ = 0.04;
//double HbbStyle::textSize_ = 0.035;
double HbbStyle::textSize_ = 0.05;

// --------------------------------------------------------------
inline void HbbStyle::setXCoordinatesL(const double relWidth, double& x0, double& x1) {
  x0 = gStyle->GetPadLeftMargin()+margin_;
  x1 = x0 + relWidth*(1.-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin()-2.*margin_);
}


// --------------------------------------------------------------
inline void HbbStyle::setXCoordinatesR(const double relWidth, double& x0, double& x1) {
  x0 = 1.-gStyle->GetPadRightMargin()-margin_-relWidth*(1.-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin()-2.*margin_);
  x1 = 1.-gStyle->GetPadRightMargin()-margin_;
}


// --------------------------------------------------------------
inline void HbbStyle::setYCoordinatesT(const int nEntries, double& y0, double& y1) {
  y1 = 1.-gStyle->GetPadTopMargin()-margin_;
  y0 = y1-nEntries*lineHeight_;
}


// --------------------------------------------------------------
inline void HbbStyle::setYCoordinatesB(const int nEntries, double& y0, double& y1) {
  y1 = gStyle->GetPadBottomMargin()+margin_;
  y0 = y1+nEntries*lineHeight_;
}


// --------------------------------------------------------------
inline TLegend* HbbStyle::legend(const int nEntries, const double relWidth, const bool left, const bool top) {
  double x0 = 0.;
  double x1 = 0.;
  double y0 = 0.;
  double y1 = 0.;
  if( left ) setXCoordinatesL(relWidth,x0,x1);
  else       setXCoordinatesR(relWidth,x0,x1);
  if( top  ) setYCoordinatesT(nEntries,y0,y1);
  else       setYCoordinatesB(nEntries,y0,y1);

  TLegend* leg = new TLegend(x0,y0,x1,y1);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(textSize_);

  return leg;
}


// --------------------------------------------------------------
inline TPaveText* HbbStyle::label(const int nEntries, const double relWidth, const bool left, const bool top) {
  double x0 = 0.;
  double x1 = 0.;
  double y0 = 0.;
  double y1 = 0.;
  if( left ) setXCoordinatesL(relWidth,x0,x1);
  else       setXCoordinatesR(relWidth,x0,x1);
  if( top  ) setYCoordinatesT(nEntries,y0,y1);
  else       setYCoordinatesB(nEntries,y0,y1);

  TPaveText* label = new TPaveText(x0,y0,x1,y1,"NDC");
  label->SetBorderSize(0);
  label->SetFillColor(0);
  label->SetFillStyle(0);
  label->SetTextFont(42);
  label->SetTextAlign(12);	// left adjusted and vertically centered
  label->SetTextSize(textSize_);
  label->SetMargin(0.);

  return label;
}


// --------------------------------------------------------------
inline TPaveText* HbbStyle::title(const TString& txt) {
  double x0 = gStyle->GetPadLeftMargin()*1.2;
  double x1 = 1.-gStyle->GetPadRightMargin();
  double y0 = 1.-gStyle->GetPadTopMargin()+0.005;
  double y1 = 1.;
  TPaveText* theTitle = new TPaveText(x0,y0,x1,y1,"NDC");
  theTitle->SetBorderSize(0);
  theTitle->SetFillColor(10);
  theTitle->SetFillStyle(1001);
  theTitle->SetTextFont(42);
  theTitle->SetTextAlign(12);	// left adjusted and vertically centered
  theTitle->SetTextSize(textSize_);
  theTitle->SetMargin(0.);
  theTitle->AddText(txt);

  return theTitle;
}


// --------------------------------------------------------------
inline TString HbbStyle::header(const PublicationStatus status) {
  TString txt = "35.7 fb^{-1} (13 TeV)";
  if( status == INTERNAL_SIMULATION ) {
    txt = "Simulation (8 TeV)";
  } else if( status == PRELIMINARY ) {
    txt = "CMS preliminary,  "+txt;
  } else if( status == PUBLIC ) {
    txt = "CMS,  "+txt;
  } else if( status == SIMULATION ) {
    txt = "CMS Simulation (8 TeV)";
  } else if( status == UNPUBLISHED ) {
    txt = "CMS (unpublished),  "+txt;
  } else if( status == PRIVATE ) {
	  txt = "Work in progress,  " + txt;
  }

  return txt;
}


// --------------------------------------------------------------
inline TString HbbStyle::bkgTemplDisplayName(const BkgTemplateType type) {
    TString name = "(B2+B1+C2,b)b";
    if(      type == C1bb     ) name = "(C1,b)b";
    else if( type == Qbb      ) name = "(Q,b)b";
    else if( type == bbB2B1C2 ) name = "bb(B2+B1+C2)";
    else if( type == bbC1Q    ) name = "bb(C1+Q)";
    
    return name;
}


// --------------------------------------------------------------
inline void HbbStyle::set(const PublicationStatus status) {
  // Store the PublicationStatus for later usage, e.g. in the title
  publicationStatus_ = status;

  // Suppress message when canvas has been saved
  gErrorIgnoreLevel = 1001;

  // Zero horizontal error bars
  gStyle->SetErrorX(0);

  //  For the canvas
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(800); //Height of canvas
  gStyle->SetCanvasDefW(800); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);

  //  For the frame
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(10);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth((Width_t) 1.);

  //  For the Pad
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  //  Margins
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.05);

  //  For the histo:
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineStyle(0);
//  gStyle->SetHistLineWidth(3);
  gStyle->SetMarkerSize(1.25);
  gStyle->SetMarkerStyle(20);
//  gStyle->SetEndErrorSize(4);
  gStyle->SetHatchesLineWidth(1);

  //  For the statistics box:
  gStyle->SetOptStat(0);

  //  For the axis
  gStyle->SetAxisColor(1,"XYZ");
  gStyle->SetTickLength(0.03,"XYZ");
  gStyle->SetNdivisions(510,"XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetStripDecimals(kFALSE);

    //  For the axis labels and titles
  gStyle->SetTitleColor(1,"XYZ");
  gStyle->SetLabelColor(1,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.4);



  //  For the legend
  gStyle->SetLegendBorderSize(0);
}

#endif

