/*
 * HbbStyleClass.h
 *
 *  Created on: 17 Aug 2017
 *      Author: shevchen
 */

#ifndef ANALYSIS_MSSMHBB_INTERFACE_HBBSTYLECLASS_H_
#define ANALYSIS_MSSMHBB_INTERFACE_HBBSTYLECLASS_H_

#include "TColor.h"
#include "TError.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"

// Publication status: determines what is plotted in title
enum PublicationStatus { INTERNAL, INTERNAL_SIMULATION, PRELIMINARY, PUBLIC, SIMULATION, UNPUBLISHED, PRIVATE };
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


#endif /* ANALYSIS_MSSMHBB_INTERFACE_HBBSTYLECLASS_H_ */
