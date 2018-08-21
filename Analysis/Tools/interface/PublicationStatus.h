/*
 * PublicationStatus.h
 *
 *  Created on: 11 Jan 2018
 *      Author: shevchen
 */

#ifndef ANALYSIS_TOOLS_INTERFACE_PUBLICATIONSTATUS_H_
#define ANALYSIS_TOOLS_INTERFACE_PUBLICATIONSTATUS_H_

#include <string>

// Publication status: determines what is plotted in title
enum PublicationStatus { INTERNAL, INTERNAL_SIMULATION, PRELIMINARY, PUBLIC, SIMULATION, UNPUBLISHED, PRIVATE, PRELIMINARY_SIMULATION, WORKINPROGRESS, SUPPLEMENTARY,PRIVATE_SIMULATION };
std::string PublicationStatusToString(const PublicationStatus status);



#endif /* ANALYSIS_TOOLS_INTERFACE_PUBLICATIONSTATUS_H_ */
