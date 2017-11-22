/*
 * CMS_lumi.cpp
 *
 *  Created on: 7 Oct 2017
 *      Author: shevchen
 */

#include "Analysis/Tools/interface/CMS_lumi.h"
#include <iostream>

TString CMS_lumi::getExtraText(const PublicationStatus status){
	/*
	 * Translate Publication title to the extra text
	 */
	TString txt;
	if( status == INTERNAL_SIMULATION ) {
		txt = "Simulation";
	} else if( status == PRELIMINARY ) {
		txt = "Preliminary";
	} else if ( status == PRELIMINARY_SIMULATION) {
		txt = "Simulation Preliminary";
	} else if( status == PUBLIC ) {
		txt = "";
	} else if( status == SIMULATION ) {
		txt = "Simulation";
	} else if( status == UNPUBLISHED ) {
		txt = "(unpublished)";
	} else if( status == PRIVATE ) {
		txt = "Private work";
	} else if ( status == WORKINPROGRESS ) {
		txt = "Work in progress";
	}

	return txt;
}

void CMS_lumi::drawCMSlumi(TPad *pad, const PublicationStatus status, const int& iPeriod, const TString& CMS_string_position, const TString& CMS_string_alignment ){
	/*
	 * Overwrite drawCMSlumi method with a human readable CMS_string_position option.
	 * CMS_string_position can be : 'left-top', 'right-top', 'center-top', 'out'
	 * CMS_string_alignment can be: 'left','right','centered'
	 */
	TString loweredPosition 		= CMS_string_position;
	TString loweredAlignment 	= CMS_string_alignment;
	loweredAlignment.ToLower();
	loweredPosition.ToLower();
	float iPosx = 0;
	//Check position option specified by user
	if( !loweredPosition.Contains("top") && !loweredPosition.Contains("out") ) {
		std::cout<<"WARNING: Precise position of the CMS label was not specified. Use default 'top' location."<<std::endl;
		loweredPosition += "-top";
	}
	//Check alignment options specified by user
	if ( !loweredAlignment.Contains("left") && !loweredAlignment.Contains("right") && !loweredAlignment.Contains("centered")){
		std::cout<<"WARNING: Alignment of the CMS label was not specified. Use default 'left' alignment."<<std::endl;
		loweredAlignment += "left";
	}

	//Get iPos according to the user input
	// Setup position
	if ( loweredPosition.Contains("left-top")) iPosx += 1;		// Default one
	else if (loweredPosition.Contains("centert-top")) iPosx += 2;
	else if (loweredPosition.Contains("right-top")) iPosx += 3;
	else if (loweredPosition.Contains("out")) iPosx += 0;
	else iPosx += 1;

	// Setup alignment
	if ( loweredAlignment.Contains("left")) iPosx += 10;			// Default one
	else if ( loweredAlignment.Contains("center")) iPosx += 20;
	else if ( loweredAlignment.Contains("right")) iPosx += 30;
	else iPosx += 10;

	if(loweredPosition.Contains("out")) iPosx =0;
	drawCMSlumi(pad,status,iPeriod,iPosx);

}

void CMS_lumi::drawCMSlumi( TPad* pad, const PublicationStatus status, int iPeriod, int iPosX ){
	/*
	 * Draw Standard CMS collaboration sign with additional
	 * text, e.g. Preliminary, Provate Work, Simulation...
	 */
	bool outOfFrame    = false;
	if( iPosX/10==0 ) outOfFrame = true;
	auto extratext = getExtraText(status);

	float t = pad->GetTopMargin();
	float r = pad->GetRightMargin();
	//  float e = 0.025;

	pad->cd();

	TString lumiText = getLumiText(iPeriod,outOfFrame);

	TLatex latex;

	//Setup Lumi latex object
	drawLumiLatex(latex,lumiText,t,r,outOfFrame);

	//Setup CMS latex object
	drawCMSLatex(latex,extratext,iPosX,pad,outOfFrame);
}

int CMS_lumi::getXAlign_(const int& iPosx){
	int alignX = 2;
	if( iPosx/10==0 ) alignX=1;
	else if( iPosx==0    ) alignX=1;
	else if( iPosx/10==1 ) alignX=1;
	else if( iPosx/10==2 ) alignX=2;
	else if( iPosx/10==3 ) alignX=3;
	else alignX = 2;

	return alignX;
}

int CMS_lumi::getYAlign_(const int& iPosx){
	int alignY = 3;
	if( iPosx==0    ) alignY=1;
	return alignY;
}

int CMS_lumi::getTotAlign_(const int& alignX, const int& alignY){
	return 10*alignX + alignY;
}

TString CMS_lumi::get7TeVLumiText_(){
	/*
	 * 7 TeV Lumi Sign
	 */
	TString lumiText;
	lumiText += lumi_7TeV_;
	lumiText += " (7 TeV)";
	return lumiText;
}

TString CMS_lumi::get8TeVLumiText_(){
	/*
	 * 8 TeV Lumi Sign
	 */
	TString lumiText;
	lumiText += lumi_8TeV_;
	lumiText += " (8 TeV)";
	return lumiText;
}

TString CMS_lumi::get13TeVLumiText_(){
	/*
	 * 13 TeV Lumi sign
	 */
	TString lumiText;
	lumiText += lumi_13TeV_;
	lumiText += " (13 TeV)";
	return lumiText;

}
TString CMS_lumi::get7p8TeVLumiText_(){
	/*
	 * 8 + 7 TeV Lumi sign
	 */
	TString lumiText;
    lumiText = lumi_8TeV_;
    lumiText += " (8 TeV)";
    lumiText += " + ";
    lumiText += lumi_7TeV_;
    lumiText += " (7 TeV)";
    return lumiText;
}

TString CMS_lumi::get7p8p13TeVLumiText_(){
	/*
	 * 13 + 8 + 7 TeV Lumi sign
	 */
	TString lumiText;
    lumiText += lumi_13TeV_;
    lumiText += " (13 TeV)";
    lumiText += " + ";
    lumiText += lumi_8TeV_;
    lumiText += " (8 TeV)";
    lumiText += " + ";
    lumiText += lumi_7TeV_;
    lumiText += " (7 TeV)";
    return lumiText;
}

TString CMS_lumi::getLumiText(const int& iPeriod, const bool& outOfFrame){
	TString lumiText;
	if( iPeriod==1 ) lumiText = get7TeVLumiText_();
	else if ( iPeriod==2 ) lumiText = get8TeVLumiText_();
	else if ( iPeriod==3 ) lumiText = get7p8TeVLumiText_();
	else if ( iPeriod==4 ) lumiText = get13TeVLumiText_();
	else if ( iPeriod==7 ){
		if( outOfFrame) lumiText += "#scale[0.85]{" + get7p8p13TeVLumiText_() + "}";
		else lumiText = get7p8p13TeVLumiText_();
	}
	else if ( iPeriod==0 ){
		lumiText += lumi_sqrtS_;
	}
	return lumiText;
}

void CMS_lumi::drawLumiLatex(TLatex & latex,TString & lumiText, const float& t, const float& r, const bool& outOfFrame){
	/*
	 * Return Latex object with lumitext
	 */
	//Setup style
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);
	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(lumiTextSize_*t);

	//Draw Lumi latex
	latex.DrawLatex(1-r,1-t+lumiTextOffset_*t,lumiText);

}
float CMS_lumi::setCMSsignXpos_(const int& iPosX, const float& t, const float&  l, const float& r, const float& b, const float& W, const float& H, const bool& extra_text){
	/*
	 * position of the CMS sign
	 */
	float posX=0;
	if( iPosX%10<=1 ){
		posX =   l + relPosX_*(1-l-r);
	}
	else if( iPosX%10==2 ){
		posX =  l + 0.5*(1-l-r);
	}
	else if( iPosX%10==3 ){
		posX =  1-r - relPosX_*(1-l-r);
	}

	//Special sets for logo
	if( drawLogo_ ) posX =   l + 0.045*(1-l-r)*W/H;
	//Special sets if add more text
	if( extra_text && iPosX == 0 ) {
//		relPosX_ = 0.08;
//		posX =   l +  relPosX_*(1-l-r);
		posX =   l +  (relPosX_-0.025)*(1-l-r);
	}

	return posX;
}

float CMS_lumi::setCMSsignYpos_(const int& iPosX, const float& t, const float& b, const bool& extra_text){
	/*
	 * position of the CMS sign
	 */
	float posY=0;
	posY = 1-t - 0.045*(1-t-b);

	//Special sets for logo
	if( drawLogo_ ) posY = 1-t - 0.045*(1-t-b);
	//Special sets if add more text
	if( extra_text && iPosX == 0 ) posY =   1-t+lumiTextOffset_*t;

	return posY;
}

void CMS_lumi::drawCMSLatex(TLatex& latex, TString& extraText, const int& iPosX, TPad * pad, const bool& outOfFrame){
	/*
	 * Setup CMS sign latex preferences
	 */
	bool drawExtraText = true;
	TString cmsText = "CMS";
	if(extraText == "") drawExtraText = false;
	float extraTextSize = extraOverCmsTextSize_*cmsTextSize_;

	float H = pad->GetWh();
	float W = pad->GetWw();
	float l = pad->GetLeftMargin();
	float t = pad->GetTopMargin();
	float r = pad->GetRightMargin();
	float b = pad->GetBottomMargin();

	int alignX = getXAlign_(iPosX);
	int alignY = getYAlign_(iPosX);
	int align_ = getTotAlign_(alignX, alignY);

	pad->cd();

	auto posX = setCMSsignXpos_(iPosX, t, l, r, b, W, H, drawExtraText);
	auto posY = setCMSsignYpos_(iPosX, t, b, drawExtraText);

	TLatex* CMS_latex = nullptr;
	float CMS_latex_width = 0;
	if( outOfFrame ){
		latex.SetTextFont(cmsTextFont_);
		latex.SetTextAlign(11);
		latex.SetTextSize(cmsTextSize_*t);
		CMS_latex = latex.DrawLatex(l,1-t+lumiTextOffset_*t,cmsText);
		CMS_latex_width = CMS_latex->GetXsize();
	}
	pad->cd();
	if( !outOfFrame ){
		if( drawLogo_ )
		{
			float xl_0 = posX;
			float yl_0 = posY - 0.15;
			float xl_1 = posX + 0.15*H/W;
			float yl_1 = posY;
	//	  TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
			TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
			pad_logo->Draw();
			pad_logo->cd();
	//	  CMS_logo->Draw("X");
			pad_logo->Modified();
			pad->cd();
		}
		else {
			latex.SetTextFont(cmsTextFont_);
			latex.SetTextSize(cmsTextSize_*t);
			latex.SetTextAlign(align_);
			CMS_latex = latex.DrawLatex(posX, posY, cmsText);
			CMS_latex_width = CMS_latex->GetXsize();
			if( drawExtraText ){
				latex.SetTextFont(extraTextFont_);
				latex.SetTextAlign(align_);
				latex.SetTextSize(extraTextSize*t);
				latex.DrawLatex(posX, posY- relExtraDY_*cmsTextSize_*t, extraText);
			}
		}
	}
	else if( drawExtraText ){
		latex.SetTextFont(extraTextFont_);
		latex.SetTextSize(extraTextSize*t);
		latex.SetTextAlign(align_);
		latex.DrawLatex(posX + CMS_latex_width, posY, extraText);
	    }
}

