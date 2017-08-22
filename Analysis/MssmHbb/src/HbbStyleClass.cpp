/*
 * HbbStyleClass.cpp
 *
 *  Created on: 17 Aug 2017
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/HbbStyleClass.h"

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

PublicationStatus HbbStyle::publicationStatus_ = INTERNAL;
double HbbStyle::lineHeight_ = 0.042;
double HbbStyle::margin_ = 0.04;
//double HbbStyle::textSize_ = 0.035;
double HbbStyle::textSize_ = 0.05;

// --------------------------------------------------------------
void HbbStyle::setXCoordinatesL(const double relWidth, double& x0, double& x1) {
  x0 = gStyle->GetPadLeftMargin()+margin_;
  x1 = x0 + relWidth*(1.-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin()-2.*margin_);
}


// --------------------------------------------------------------
void HbbStyle::setXCoordinatesR(const double relWidth, double& x0, double& x1) {
  x0 = 1.-gStyle->GetPadRightMargin()-margin_-relWidth*(1.-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin()-2.*margin_);
  x1 = 1.-gStyle->GetPadRightMargin()-margin_;
}


// --------------------------------------------------------------
void HbbStyle::setYCoordinatesT(const int nEntries, double& y0, double& y1) {
  y1 = 1.-gStyle->GetPadTopMargin()-margin_;
  y0 = y1-nEntries*lineHeight_;
}


// --------------------------------------------------------------
void HbbStyle::setYCoordinatesB(const int nEntries, double& y0, double& y1) {
  y1 = gStyle->GetPadBottomMargin()+margin_;
  y0 = y1+nEntries*lineHeight_;
}


// --------------------------------------------------------------
TLegend* HbbStyle::legend(const int nEntries, const double relWidth, const bool left, const bool top) {
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
TPaveText* HbbStyle::label(const int nEntries, const double relWidth, const bool left, const bool top) {
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
TPaveText* HbbStyle::title(const TString& txt) {
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
TString HbbStyle::header(const PublicationStatus status) {
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
TString HbbStyle::bkgTemplDisplayName(const BkgTemplateType type) {
    TString name = "(B2+B1+C2,b)b";
    if(      type == C1bb     ) name = "(C1,b)b";
    else if( type == Qbb      ) name = "(Q,b)b";
    else if( type == bbB2B1C2 ) name = "bb(B2+B1+C2)";
    else if( type == bbC1Q    ) name = "bb(C1+Q)";
    
    return name;
}


// --------------------------------------------------------------
void HbbStyle::set(const PublicationStatus status) {
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





