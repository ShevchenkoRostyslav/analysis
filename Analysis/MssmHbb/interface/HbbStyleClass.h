/*
 * HbbStyleClass.h
 *
 *  Created on: 17 Aug 2017
 *      Author: shevchen
 */

#ifndef ANALYSIS_MSSMHBB_INTERFACE_HBBSTYLECLASS_H_
#define ANALYSIS_MSSMHBB_INTERFACE_HBBSTYLECLASS_H_

#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TBox.h"
#include "Analysis/Tools/interface/CMS_lumi.h"
#include "iostream"

class HbbStyle {
public:
  // Adjusts the gStyle settings and store the PublicationStatus
  static void set(const PublicationStatus status);
  static void setTDRstyle(const PublicationStatus status);
  static PublicationStatus status() { return publicationStatus_; }

  // Draws a title on the current pad
  //  <CMS label>,  #sqrt{s} = 13 TeV,  35.7fb^{-1}
  // where <CMS label> depends on the PublicationStatus
  //  INTERNAL    : no extra label (intended for AN-only plots with data)
  //  INTERNAL    : show "Simulation" label (intended for AN-only plots, no lumi, no "CMS")
  //  PRELIMINARY : show "CMS preliminary" label
  //  PUBLIC      : show "CMS" label
  //  SIMULATION  : show "CMS Simulation" label
  //  UNPUBLISHED : show "CMS (unpublished)" label (intended for additional material on TWiki)
  // Note that this method does not allow for easy memory
  // handling. For that, use standardTitle().
  static void drawStandardTitle();
  static void drawStandardTitle(const TString& CMS_string_position, const TString& CMS_string_alignment = "left", const int& iPeriod = 4);
//  static void drawStandardTitle() { standardTitle()->Draw("same"); }

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

  /*
   * Hack the TLegend header in case it's too long to fit into a single line
   */
  static void drawLegendSplittedHeader(TLegend * leg, const char* split_character);

  static void drawLegendOnTopOfThePad(TLegend *leg, TCanvas *can, const double& pad_size);


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


  // for axis titles
  static TString axisTitleMass() { return "M_{12} [GeV]"; }
  static TString axisTitleMAH() { return "m_{A/H} [GeV]"; }
  static TString axisTitleMA() { return "m_{A} [GeV]"; }
  
  //Luminosity string
  static TString lumiString() {return "35.7 fb^{-1} (13 TeV)"; }

  //Style for the limits plots
  static void setObservedLimitsStyle(TGraph *gr);
  static void setExpectedLimitsStyle(TGraph *gr);
  static void set1SigmaBandsStyle(TGraph *gr);
  static void set2SigmaBandsStyle(TGraph *gr);
  static void setFrameStyle(TH2 *frame);
  
  //Style for the Ratio plots
  static TPad* getRatioTopPad(const float& relative_size = 0.7);
  static TPad* getRatioBottomPad(const float& relative_size = 0.3);
  static void setRatioTopFrame(TH2D *h, const double& canva_h);
  static void setRatioTopFrame(TAxis *xaxis, TAxis *yaxis, const double& canva_h);
  static void setRatioBottomFrame(TH2D *h, const double& topPad_w, const double& canva_h);
  static void setRatioBottomFrame(TAxis *xaxis, TAxis *yaxis, const double& canva_h, const double& topPad_h);


  //Style for the TLegend
  static void setLegendStyle(TLegend *leg);
  static double lineHeight() { return lineHeight_; }

  //Style for the Histograms
  static void correctForTheBinWidth(TH1& h);
  static void applyStandardDataHistoStyle(TH1&h);
  static void applyStandardMCHistoStyle(TH1&h);


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

  // Change the TPad coordinates
  static void setCurrentPadCoordinates(const double& xlow, const double& ylow, const double& xup, const double& yup);
  // Move TLegend to a certain position
  static void moveTLegend(TLegend *pad, const double& xlow, const double& ylow, const double& xup, const double& yup);

  static TLegend* legend(const int nEntries, const double relWidth, const bool left, const bool top);
  static TPaveText* label(const int nEntries, const double relWidth, const bool leftt, const bool top);

};


#endif /* ANALYSIS_MSSMHBB_INTERFACE_HBBSTYLECLASS_H_ */
