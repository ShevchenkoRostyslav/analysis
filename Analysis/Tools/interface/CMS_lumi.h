/*
 * CMS_lumi.h
 *
 *  Created on: 7 Oct 2017
 *      Author: shevchen
 */

#ifndef ANALYSIS_TOOLS_INTERFACE_CMS_LUMI_H_
#define ANALYSIS_TOOLS_INTERFACE_CMS_LUMI_H_

#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"
#include <memory>
//#include "Analysis/MssmHbb/interface/HbbStyleClass.h"

// Publication status: determines what is plotted in title
enum PublicationStatus { INTERNAL, INTERNAL_SIMULATION, PRELIMINARY, PUBLIC, SIMULATION, UNPUBLISHED, PRIVATE, PRELIMINARY_SIMULATION, WORKINPROGRESS };

class CMS_lumi{
public:
	void drawCMSlumi(TPad *pad, const PublicationStatus status, int iPeriod, int iPosX=10 );
	void drawCMSlumi(TPad *pad, const PublicationStatus status, const int& iPeriod = 4 , const TString& CMS_string_position  = "top-left", const TString& CMS_string_alignment = "left" );
	TString getLumiText(const int& iPeriod, const bool& outOfFrame = false);
	TString getExtraText(const PublicationStatus status);
	void drawLumiLatex(TLatex & latex,TString & lumiText, const float& t, const float& r, const bool& outOfFrame = false);
	void drawCMSLatex(TLatex& latex, TString & extratext, const int& iPosX, TPad * pad, const bool& outOfFrame = false);

	float getCmsTextFont() const;
	void setCmsTextFont(const float& cmsTextFont = 61);
	float getCmsTextOffset() const;
	void setCmsTextOffset(const float& cmsTextOffset = 0.1);
	float getCmsTextSize() const;
	void setCmsTextSize(const float& cmsTextSize = 0.75);
	bool isDrawLogo() const;
	void setDrawLogo(const bool& drawLogo = false);
	float getExtraOverCmsTextSize() const;
	void setExtraOverCmsTextSize(const float& extraOverCmsTextSize = 0.76);
	float getExtraTextFont() const;
	void setExtraTextFont(const float& extraTextFont = 52);
	const TString& getLumi13TeV() const;
	void setLumi13TeV(const TString& lumi13TeV = "35.7 fb^{-1}");
	const TString& getLumi7TeV() const;
	void setLumi7TeV(const TString& lumi7TeV = "5.1 fb^{-1}");
	const TString& getLumi8TeV() const;
	void setLumi8TeV(const TString& lumi8TeV = "19.7 fb^{-1}");
	const TString& getLumiSqrtS() const;
	void setLumiSqrtS(const TString& lumiSqrtS = "");
	float getLumiTextOffset() const;
	void setLumiTextOffset(const float& lumiTextOffset = 0.2);
	float getLumiTextSize() const;
	void setLumiTextSize(const float& lumiTextSize = 0.6);
	float getRelExtraDy() const;
	void setRelExtraDy(const float& relExtraDy = 1.2);
	float getRelPosX() const;
	void setRelPosX(const float& relPosX = 0.045);
	float getRelPosY() const;
	void setRelPosY(const float& relPosY = 0.035);

protected:

	int getXAlign_(const int& iPosx);
	int getYAlign_(const int& iPosx);
	int getTotAlign_(const int& alignX, const int& alignY);

	TString get7TeVLumiText_();
	TString get8TeVLumiText_();
	TString get13TeVLumiText_();
	TString get7p8TeVLumiText_();
	TString get7p8p13TeVLumiText_();

	float setCMSsignXpos_(const int& iPosX, const float& t, const float&  l, const float& r, const float& b, const float& W, const float& H, const bool& extraText);
	float setCMSsignYpos_(const int& iPosX, const float& t, const float& b, const bool& extraText);

	float cmsTextFont_   = 61;  // default is helvetic-bold
	float extraTextFont_ = 52;  // default is helvetica-italics

	// text sizes and text offsets with respect to the top frame
	// in unit of the top margin size
	float lumiTextSize_     = 0.6;
	float lumiTextOffset_   = 0.2;
	float cmsTextSize_      = 0.75;
	float cmsTextOffset_    = 0.1;  // only used in outOfFrame version

	float relPosX_    = 0.045;
	float relPosY_    = 0.035;
	float relExtraDY_ = 1.2;

	// ratio of "CMS" and extra text size
	float extraOverCmsTextSize_  = 0.76;

	TString lumi_13TeV_ = "35.7 fb^{-1}";
	TString lumi_8TeV_  = "19.7 fb^{-1}";
	TString lumi_7TeV_  = "5.1 fb^{-1}";
	TString lumi_sqrtS_ = "";

	bool drawLogo_      = false;
};

inline float CMS_lumi::getCmsTextFont() const {
	return cmsTextFont_;
}

inline void CMS_lumi::setCmsTextFont(const float& cmsTextFont) {
	cmsTextFont_ = cmsTextFont;
}

inline float CMS_lumi::getCmsTextOffset() const {
	return cmsTextOffset_;
}

inline void CMS_lumi::setCmsTextOffset(const float& cmsTextOffset) {
	cmsTextOffset_ = cmsTextOffset;
}

inline float CMS_lumi::getCmsTextSize() const {
	return cmsTextSize_;
}

inline void CMS_lumi::setCmsTextSize(const float& cmsTextSize) {
	cmsTextSize_ = cmsTextSize;
}

inline bool CMS_lumi::isDrawLogo() const {
	return drawLogo_;
}

inline void CMS_lumi::setDrawLogo(const bool& drawLogo) {
	drawLogo_ = drawLogo;
}

inline float CMS_lumi::getExtraOverCmsTextSize() const {
	return extraOverCmsTextSize_;
}

inline void CMS_lumi::setExtraOverCmsTextSize(const float& extraOverCmsTextSize) {
	extraOverCmsTextSize_ = extraOverCmsTextSize;
}

inline float CMS_lumi::getExtraTextFont() const {
	return extraTextFont_;
}

inline void CMS_lumi::setExtraTextFont(const float& extraTextFont) {
	extraTextFont_ = extraTextFont;
}

inline const TString& CMS_lumi::getLumi13TeV() const {
	return lumi_13TeV_;
}

inline void CMS_lumi::setLumi13TeV(const TString& lumi13TeV) {
	lumi_13TeV_ = lumi13TeV;
}

inline const TString& CMS_lumi::getLumi7TeV() const {
	return lumi_7TeV_;
}

inline void CMS_lumi::setLumi7TeV(const TString& lumi7TeV) {
	lumi_7TeV_ = lumi7TeV;
}

inline const TString& CMS_lumi::getLumi8TeV() const {
	return lumi_8TeV_;
}

inline void CMS_lumi::setLumi8TeV(const TString& lumi8TeV) {
	lumi_8TeV_ = lumi8TeV;
}

inline const TString& CMS_lumi::getLumiSqrtS() const {
	return lumi_sqrtS_;
}

inline void CMS_lumi::setLumiSqrtS(const TString& lumiSqrtS) {
	lumi_sqrtS_ = lumiSqrtS;
}

inline float CMS_lumi::getLumiTextOffset() const {
	return lumiTextOffset_;
}

inline void CMS_lumi::setLumiTextOffset(const float& lumiTextOffset) {
	lumiTextOffset_ = lumiTextOffset;
}

inline float CMS_lumi::getLumiTextSize() const {
	return lumiTextSize_;
}

inline void CMS_lumi::setLumiTextSize(const float& lumiTextSize) {
	lumiTextSize_ = lumiTextSize;
}

inline float CMS_lumi::getRelExtraDy() const {
	return relExtraDY_;
}

inline void CMS_lumi::setRelExtraDy(const float& relExtraDy) {
	relExtraDY_ = relExtraDy;
}

inline float CMS_lumi::getRelPosX() const {
	return relPosX_;
}

inline void CMS_lumi::setRelPosX(const float& relPosX) {
	relPosX_ = relPosX;
}

inline float CMS_lumi::getRelPosY() const {
	return relPosY_;
}

inline void CMS_lumi::setRelPosY(const float& relPosY) {
	relPosY_ = relPosY;
}

#endif /* ANALYSIS_TOOLS_INTERFACE_CMS_LUMI_H_ */
