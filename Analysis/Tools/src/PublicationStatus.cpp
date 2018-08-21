/*
 * PublicationStatus.cpp
 *
 *  Created on: 11 Jan 2018
 *      Author: shevchen
 */

#include "Analysis/Tools/interface/PublicationStatus.h"

std::string PublicationStatusToString(const PublicationStatus status){
	/* Convert publicationStatus to string
 	*
 	*/
	if(status == INTERNAL) return "INTERNAL";
	else if (status == INTERNAL_SIMULATION) return "INTERNAL_SIMULATION";
	else if (status == PRELIMINARY) return "PRELIMINARY";
	else if (status == PUBLIC) return "PUBLIC";
	else if (status == SIMULATION) return "SIMULATION";
	else if (status == UNPUBLISHED) return "UNPUBLISHED";
	else if (status == PRIVATE) return "PRIVATE";
	else if (status == PRIVATE_SIMULATION) return "PRIVATE_SIMULATION";
	else if (status == PRELIMINARY_SIMULATION) return "PRELIMINARY_SIMULATION";
	else if (status == WORKINPROGRESS) return "WORKINPROGRESS";
	else if (status == SUPPLEMENTARY) return "SUPPLEMENTARY";
	else return "";
}
