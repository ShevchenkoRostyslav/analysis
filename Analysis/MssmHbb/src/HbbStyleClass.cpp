/*
 * HbbStyleClass.cpp
 *
 *  Created on: 17 Aug 2017
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/HbbStyleClass.h"

TString toTString(const PublicationStatus status) {
  TString str = "";
  if(      status == INTERNAL )            			str = "internal";
  else if( status == INTERNAL_SIMULATION ) 			str = "simulation (internal)";
  else if( status == PRELIMINARY )         			str = "preliminary";
  else if( status == PRELIMINARY_SIMULATION )        str = "simulation preliminary";
  else if( status == PUBLIC      )         			str = "public";
  else if( status == SIMULATION  )         			str = "simulation (public)";
  else if( status == UNPUBLISHED )         			str = "unpublished";
  else if( status == WORKINPROGRESS )			   	str = "Work in progress";
  else if( status == PRIVATE )			   			str = "Private work";
  else if( status == SUPPLEMENTARY )				str = "Supplementary";

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
  TString txt = lumiString();
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
	  txt = "Private work  " + txt;
  } else if( status == WORKINPROGRESS){
	  txt = "Work in progress,  " + txt;
  }

  return txt;
}

void setPadCoordinates(const double& xlow, const double& ylow, const double& xup, const double& yup){
	gPad->SetPad(xlow, ylow, xup, yup);
}

void HbbStyle::moveTLegend(TLegend *pad, const double& xlow, const double& ylow, const double& xup, const double& yup){
	pad->SetX1(xlow);
	pad->SetX2(xup);
	pad->SetY1(ylow);
	pad->SetY2(yup);
}

void HbbStyle::drawLegendOnTopOfThePad(TLegend *leg, TCanvas *can, const double& pad_size){
	/*
	 * draw TLegend in a separate pad on top of the plot
	 */
	//Get the current pad and resize it
	auto *main_pad = static_cast<TPad*>(gPad);
	setPadCoordinates(main_pad->GetX1(),main_pad->GetY1(),main_pad->GetX2(),(1 - pad_size));
	main_pad->SetTopMargin(0.02);
	can->cd();

	//Make a Box for the upper part
	TPaveText *top_box = new TPaveText(main_pad->GetX1() + main_pad->GetLeftMargin(),1 - pad_size - 0.015, main_pad->GetX2() - main_pad->GetRightMargin(),1-gPad->GetTopMargin(),"NB");
	top_box->SetBorderSize(1);
	top_box->SetLineColor(kBlack);
	top_box->SetLineWidth(2);
	top_box->SetLineStyle(1);
	top_box->SetFillColor(kWhite);
	top_box->Draw();

	//Make TPad for the upper part
//	TPad *top_pad = new TPad("top_pad", "top_pad",main_pad->GetX1() + main_pad->GetLeftMargin(),1 - pad_size - 0.015, main_pad->GetX2() - main_pad->GetRightMargin(),1-gPad->GetTopMargin());
//	top_pad->Draw();
//	top_pad->cd();


	//Adjust legend style
//	moveTLegend(leg,leg->GetX1(),0.05,leg->GetX2(),0.9);
//	leg->SetTextSize(gStyle->GetTextSize());

	if(leg->GetListOfPrimitives()->GetSize() > 5) {
		leg->SetNColumns(3);
		moveTLegend(leg,main_pad->GetLeftMargin()*1.02,leg->GetY1()*1.02,1-main_pad->GetRightMargin()*1.08,leg->GetY2());
	}
	else {
		leg->SetNColumns(2);
//		auto left = leg->GetX1(), right = leg->GetX2();
		auto diff = leg->GetX2() - leg->GetX1();
		moveTLegend(leg,main_pad->GetLeftMargin()*1.02,leg->GetY1()*1.02,main_pad->GetLeftMargin()*1.02 + diff,leg->GetY2());
	}
	leg->SetBorderSize(0);
	leg->SetLineWidth(0);
	leg->SetLineColorAlpha(kWhite, 1);
//	leg->SetLineColor(kBlack);
	leg->Draw();
//	can->cd();
}


void HbbStyle::drawLegendSplittedHeader(TLegend * leg, const char* split_character){
	/*
	 * Hardcoded hack of the TLegend header
	 */
	std::string header = std::string(leg->GetHeader());
	leg->SetHeader();
	std::string header_p1, header_p2;
	//Search for the splitter
	auto pos = header.find(split_character);
	if(pos == std::string::npos) {
		std::stringstream msg;
		msg<<"No character "<<split_character<<" in a TLegend header";
		throw std::logic_error(msg.str());
	}
	else{
		header_p1 = header.substr(0,pos);
		header_p2 = header.substr(pos+2);
		header = "#splitline{" + header_p1 + "}{" + header_p2 + "}";

		//Define the height of one entry in the legend
		int number_of_leg_entries = leg->GetListOfPrimitives()->GetSize(); 			// -1 To get read of header
		float height_of_entry = (leg->GetY2() - leg->GetY1())/number_of_leg_entries;
		height_of_entry = leg->GetEntrySeparation();

		//Small box on top of the legend to draw Header
		float eps = 0.005;
//		auto *header_box = new TPaveText(leg->GetX1()* (1 + eps),leg->GetY2() + 0.95*height_of_entry, leg->GetX2()* (1 - eps),(leg->GetY2() + height_of_entry)*(1-2*eps),"NDC");
//		auto *header_box = new TPaveText(leg->GetX1(),leg->GetY2()*0.95, leg->GetX2()* (1 - eps),(leg->GetY2() + height_of_entry)*(1-2*eps),"NDC");
		auto *header_box = new TPaveText(leg->GetX1(),leg->GetY2()+ height_of_entry*0.5, leg->GetX2(),(leg->GetY2() + height_of_entry*0.6),"NDC");
		header_box->AddText(header.c_str());
		header_box->SetTextAlign(kHAlignLeft+kVAlignTop);
//		header_box->SetTextAlign(kHAlignRight+kVAlignTop);
		header_box->SetFillColorAlpha(0,1);
		header_box->SetTextSize(0.035);
		header_box->SetLineWidth(0);
		header_box->SetBorderSize(0);
		header_box->Draw();

		//Empty box to draw borders
//		auto *additional_box = new TPaveText(leg->GetX1()*(1-eps),leg->GetY1()*(1-0.5*eps),leg->GetX2()*(1+0.5*eps),(leg->GetY2() +  height_of_entry)*(1+0.5*eps),"NDC");
		auto *additional_box = new TPaveText(leg->GetX1()*(1-eps),leg->GetY1()*(1-0.52*eps),leg->GetX2()*(1+0.5*eps),(leg->GetY2() +  height_of_entry)*0.93,"NDC");
		additional_box->SetFillColor(leg->GetFillColor());
		additional_box->SetTextSize(leg->GetTextSize());
		additional_box->SetLineWidth(leg->GetLineWidth());
		additional_box->SetBorderSize(leg->GetBorderSize());
		additional_box->Draw("same");

		//Draw legend Prepare legend box
		leg->SetFillColorAlpha(0,1);
		leg->SetLineWidth(0);
		leg->SetBorderSize(0);
		leg->Draw("same");

		header_box->Draw("same");
	}

}

void HbbStyle::drawStandardTitle(const TString& CMS_string_position, const TString& CMS_string_alignment, const int& iPeriod){
	/*
	 * Standard title with defined position
	 */
	TPad * pad = static_cast<TPad*>(gPad);
	CMS_lumi l;
	l.drawCMSlumi(pad, publicationStatus_,iPeriod,CMS_string_position,CMS_string_alignment);
}

void HbbStyle::drawStandardTitle(){
	/*
	 * Re-implementation of the drawStandardTitle
	 */
	TPad * pad = static_cast<TPad*>(gPad);
	CMS_lumi l;
	l.drawCMSlumi(pad, publicationStatus_);
}

void HbbStyle::setTDRstyle(const PublicationStatus status){
	/*
	 * Set TStyle according to CMS recommendations
	 */
	  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
	  publicationStatus_ = status;

	// For the canvas:
	  tdrStyle->SetCanvasBorderMode(0);
	  tdrStyle->SetCanvasColor(kWhite);
	  tdrStyle->SetCanvasDefH(600); //Height of canvas
	  tdrStyle->SetCanvasDefW(600); //Width of canvas
	  tdrStyle->SetCanvasDefX(0);   //POsition on screen
	  tdrStyle->SetCanvasDefY(0);

	// For the Pad:
	  tdrStyle->SetPadBorderMode(0);
	  // tdrStyle->SetPadBorderSize(Width_t size = 1);
	  tdrStyle->SetPadColor(kWhite);
	  tdrStyle->SetPadGridX(false);
	  tdrStyle->SetPadGridY(false);
	  tdrStyle->SetGridColor(0);
	  tdrStyle->SetGridStyle(3);
	  tdrStyle->SetGridWidth(1);

	// For the frame:
	  tdrStyle->SetFrameBorderMode(0);
	  tdrStyle->SetFrameBorderSize(1);
	  tdrStyle->SetFrameFillColor(0);
	  tdrStyle->SetFrameFillStyle(0);
	  tdrStyle->SetFrameLineColor(1);
	  tdrStyle->SetFrameLineStyle(1);
	  tdrStyle->SetFrameLineWidth(1);

	// For the histo:
	  // tdrStyle->SetHistFillColor(1);
	  // tdrStyle->SetHistFillStyle(0);
	  tdrStyle->SetHistLineColor(1);
	  tdrStyle->SetHistLineStyle(0);
	  tdrStyle->SetHistLineWidth(1);
	  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
	  // tdrStyle->SetNumberContours(Int_t number = 20);

	  tdrStyle->SetEndErrorSize(2);
	  // tdrStyle->SetErrorMarker(20);
	  //tdrStyle->SetErrorX(0.);

	  tdrStyle->SetMarkerStyle(20);

	//For the fit/function:
	  tdrStyle->SetOptFit(1);
	  tdrStyle->SetFitFormat("5.4g");
	  tdrStyle->SetFuncColor(2);
	  tdrStyle->SetFuncStyle(1);
	  tdrStyle->SetFuncWidth(1);

	//For the date:
	  tdrStyle->SetOptDate(0);
	  // tdrStyle->SetDateX(Float_t x = 0.01);
	  // tdrStyle->SetDateY(Float_t y = 0.01);

	// For the statistics box:
	  tdrStyle->SetOptFile(0);
	  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
	  tdrStyle->SetStatColor(kWhite);
	  tdrStyle->SetStatFont(42);
	  tdrStyle->SetStatFontSize(0.025);
	  tdrStyle->SetStatTextColor(1);
	  tdrStyle->SetStatFormat("6.4g");
	  tdrStyle->SetStatBorderSize(1);
	  tdrStyle->SetStatH(0.1);
	  tdrStyle->SetStatW(0.15);
	  // tdrStyle->SetStatStyle(Style_t style = 1001);
	  // tdrStyle->SetStatX(Float_t x = 0);
	  // tdrStyle->SetStatY(Float_t y = 0);

	// Margins:
	  tdrStyle->SetPadTopMargin(0.05);
	  tdrStyle->SetPadBottomMargin(0.13);
	  tdrStyle->SetPadLeftMargin(0.16);
	  tdrStyle->SetPadRightMargin(0.02);

	// For the Global title:

	  tdrStyle->SetOptTitle(0);
	  tdrStyle->SetTitleFont(42);
	  tdrStyle->SetTitleColor(1);
	  tdrStyle->SetTitleTextColor(1);
	  tdrStyle->SetTitleFillColor(10);
	  tdrStyle->SetTitleFontSize(0.05);
	  // tdrStyle->SetTitleH(0); // Set the height of the title box
	  // tdrStyle->SetTitleW(0); // Set the width of the title box
	  // tdrStyle->SetTitleX(0); // Set the position of the title box
	  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
	  // tdrStyle->SetTitleStyle(Style_t style = 1001);
	  // tdrStyle->SetTitleBorderSize(2);

	// For the axis titles:

	  tdrStyle->SetTitleColor(1, "XYZ");
	  tdrStyle->SetTitleFont(42, "XYZ");
	  tdrStyle->SetTitleSize(0.06, "XYZ");
	  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
	  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
	  tdrStyle->SetTitleXOffset(0.9);
	  tdrStyle->SetTitleYOffset(1.25);
	  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

	// For the axis labels:

	  tdrStyle->SetLabelColor(1, "XYZ");
	  tdrStyle->SetLabelFont(42, "XYZ");
	  tdrStyle->SetLabelOffset(0.007, "XYZ");
	  tdrStyle->SetLabelSize(0.05, "XYZ");

	// For the axis:

	  tdrStyle->SetAxisColor(1, "XYZ");
	  tdrStyle->SetStripDecimals(kTRUE);
	  tdrStyle->SetTickLength(0.03, "XYZ");
	  tdrStyle->SetNdivisions(510, "XYZ");
	  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	  tdrStyle->SetPadTickY(1);

	// Change for log plots:
	  tdrStyle->SetOptLogx(0);
	  tdrStyle->SetOptLogy(0);
	  tdrStyle->SetOptLogz(0);

	// Postscript options:
	  tdrStyle->SetPaperSize(20.,20.);
	  // tdrStyle->SetLineScalePS(Float_t scale = 3);
	  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
	  // tdrStyle->SetHeaderPS(const char* header);
	  // tdrStyle->SetTitlePS(const char* pstitle);

	  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
	  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
	  // tdrStyle->SetPaintTextFormat(const char* format = "g");
	  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
	  // tdrStyle->SetTimeOffset(Double_t toffset);
	  // tdrStyle->SetHistMinimumZero(kTRUE);

//	  tdrStyle->SetHatchesLineWidth(5);
//	  tdrStyle->SetHatchesSpacing(0.05);

	  tdrStyle->cd();
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

void HbbStyle::setObservedLimitsStyle(TGraph *obsG) {
  	/*
  	* Set a style for the observed limits graph
  	*/
  	obsG->SetLineWidth(3);
	obsG->SetLineColor(1);
	obsG->SetLineWidth(2);
	obsG->SetMarkerColor(1);
	obsG->SetMarkerStyle(20);
	obsG->SetMarkerSize(1.4);
}

void HbbStyle::setExpectedLimitsStyle(TGraph *expG){
	/*
	* Set a style for the expected Limits graph
	*/
	expG->SetLineWidth(3);
	expG->SetLineColor(2);
	expG->SetLineStyle(2);
}

void HbbStyle::set1SigmaBandsStyle(TGraph *innerBand){
	/*
	* Set a style for the 1 Sigma expected Limits band
	*/
	innerBand->SetFillColor(kGreen+1);
	innerBand->SetLineColor(kGreen+1);
}

void HbbStyle::set2SigmaBandsStyle(TGraph *outerBand){
	/*
	* Set a style for the 1 Sigma expected Limits band
	*/ 
	outerBand->SetFillColor(kOrange);
	outerBand->SetLineColor(kOrange);
	
}
void HbbStyle::setFrameStyle(TH2 *frame){
	/*
	* Set a style for the 2 Sigma expected limits band 
	*/
	frame->GetXaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetNdivisions(206);
	frame->GetYaxis()->SetTitleOffset(1.3);
	frame->GetXaxis()->SetTitleOffset(1.05);
	frame->GetYaxis()->SetTitleSize(0.048);
	frame->SetStats(0);
}

void HbbStyle::setLegendStyle(TLegend *leg){
	/*
	* Set common TLegend style
	*/
	leg->SetFillColor(0);
	leg->SetTextSize(0.035);
//	leg->SetLineWidth(2);
//	leg->SetBorderSize(1);
	leg->SetLineWidth(2);
	leg->SetBorderSize(0);
}

TPad* HbbStyle::getRatioTopPad(const float& relative_size){
	/*
	 * Function to produce and set up the Top pad for the ratio plot
	 */
	TPad* pad1 = new TPad("pad_top","",0,1 - relative_size,1,1);
	pad1->SetTopMargin(gStyle->GetPadTopMargin() * 1. / (relative_size));
	pad1->SetBottomMargin(0.03);
	return pad1;
}
TPad* HbbStyle::getRatioBottomPad(const float& relative_size){
	/*
	 * Function to produce and set up the Bottom pad for the ratio plot
	 */
	TPad* pad1 = new TPad("pad_top","",0,0,1,relative_size);
	pad1->SetTopMargin(0.04);
	pad1->SetBottomMargin(gStyle->GetPadBottomMargin() * 0.8 / relative_size);
	return pad1;
}

void HbbStyle::setRatioTopFrame(TH2D *frame1, const double& canva_h){
	/*
	 * Set Top frame for the raito plot
	 */
	setRatioTopFrame(frame1->GetXaxis(),frame1->GetYaxis(),canva_h);
	frame1->GetYaxis()->SetRangeUser(frame1->GetMinimum(), frame1->GetMaximum());
}

void HbbStyle::setRatioTopFrame(TAxis *xaxis, TAxis *yaxis, const double& canva_h){
	/*
	 * Set Top frame for the raito plot
	 */

	double pad_hPixel = gPad->YtoPixel(gPad->GetY1());

	yaxis->SetTickLength(gStyle->GetTickLength("Y") );
	yaxis->SetTitleOffset(gStyle->GetTitleOffset("Y"));

	xaxis->SetTickLength(gStyle->GetTickLength("X") * canva_h / pad_hPixel);
    xaxis->SetTitleSize(0);
    xaxis->SetLabelSize(0);
	xaxis->SetTitleOffset(999); //Effectively turn off x axis title on main plot
	xaxis->SetLabelOffset(999); //Effectively turn off x axis label on main plot
}

void HbbStyle::setRatioBottomFrame(TH2D *frame2, const double& canva_h, const double& topPad_h){
	/*
	 * Set Bottom frame for the ratio plot
	 */
	setRatioBottomFrame(frame2->GetXaxis(),frame2->GetYaxis(),canva_h,topPad_h);
    frame2->SetMinimum(-3.9999);
    frame2->SetMaximum(+3.9999);
}

void HbbStyle::setRatioBottomFrame(TAxis *xaxis, TAxis *yaxis, const double& canva_h, const double& topPad_h){
	/*
	 * Set Bottom frame for the ratio plot
	 * More general implementation than for the TH2D
	 */
	double pad_hPixel = gPad->YtoPixel(gPad->GetY1());

	//Y-axis options
	yaxis->SetTitleFont(43);			// This gives the sizes in pixels!!!!
	yaxis->SetLabelFont(43);			// This gives the sizes in pixels!!!!
	yaxis->SetTickLength(gStyle->GetTickLength("Y") );
	yaxis->SetLabelOffset(gStyle->GetLabelOffset("Y"));
	yaxis->SetTitleOffset(gStyle->GetTitleOffset("Y") * 1.15);
	yaxis->SetLabelSize(gStyle->GetLabelSize("Y") * topPad_h);
	yaxis->SetTitleSize(gStyle->GetTitleSize("Y") * topPad_h);
    yaxis->CenterTitle();
    yaxis->SetNdivisions(105);

    //X-axis options
	xaxis->SetTitleFont(43);			// This gives the sizes in pixels!!!!
	xaxis->SetLabelFont(43);			// This gives the sizes in pixels!!!!
    xaxis->SetTitle(HbbStyle::axisTitleMass());
//    xaxis->SetTickLength(gStyle->GetTickLength("X") * canva_h / topPad_h);
    xaxis->SetTickLength(gStyle->GetTickLength("X") * canva_h / pad_hPixel);
    xaxis->SetLabelSize(gStyle->GetLabelSize("X") * topPad_h);
    xaxis->SetTitleSize(gStyle->GetTitleSize("X") * topPad_h);
    	xaxis->SetTitleOffset(3.5);
    	xaxis->SetLabelOffset(gStyle->GetLabelOffset("X") * topPad_h / pad_hPixel );
}

void HbbStyle::correctForTheBinWidth(TH1& h){
	/*
	 * Correct histo content for the bin width
	 */
	for(int i = 1; i<=h.GetNbinsX();++i){
		h.SetBinContent(i, h.GetBinContent(i) / h.GetBinWidth(i));
		h.SetBinError(i, h.GetBinError(i) / h.GetBinWidth(i));
	}
}

void HbbStyle::applyStandardDataHistoStyle(TH1&h){
	/*
	 * Apply Standard histo style for the data
	 */
	h.SetMarkerSize(gStyle->GetMarkerSize());
	h.SetMarkerStyle(gStyle->GetMarkerStyle());
	h.SetMarkerColor(gStyle->GetMarkerColor());
	h.SetLineColor(gStyle->GetLineColor());
	h.SetLineStyle(gStyle->GetLineStyle());
	h.SetLineWidth(gStyle->GetLineWidth());
}
void HbbStyle::applyStandardMCHistoStyle(TH1&h){
	/*
	 * Apply Standard histo style for MC
	 */
	h.SetMarkerSize(gStyle->GetMarkerSize());
	h.SetMarkerStyle(gStyle->GetMarkerStyle());
	h.SetMarkerColor(gStyle->GetMarkerColor());
	h.SetLineColor(gStyle->GetLineColor());
	h.SetLineStyle(gStyle->GetLineStyle());
	h.SetLineWidth(gStyle->GetLineWidth());
}
