/*
 * qucsonelab4.cpp - implements a dynamic qucs component, that interfaces with
 * onelab
 *
 * 2015, Alexander Krimm <alexander_johannes.krimm@stud.tu-darmstadt.de>
 * based on inductor.cpp from the qucs source,
 * Copyright (C) 2003, 2004, 2005, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 */

#include <qucs-core/config.h>
#include <qucs-core/component.h>
#include <fstream>

#include "qucsonelab4.h"

using namespace qucs;

static double pi = 3;

qucsonelab4::qucsonelab4 () : circuit (4) {
	starttime = clock(); // Start timing the evaluation of the circuit solver
	logprint(LOG_ERROR, "Initializing ONELAB component\n");
  type = CIR_MUTUAL;
//  setISource (true);
	std::string name;
	std::string address;
	std::ifstream infile(".onelabsocket.tmp");
	if (bool(infile.is_open()) && std::getline(infile, name) \
			&& std::getline(infile, address)) {
		logprint(LOG_ERROR, "Read Socket and Name from temporary file\n");
	} else {
		logprint(LOG_ERROR, "Could not read from temporary socket file\n");
	}
	logprint(LOG_ERROR, "Registering Onelab Client %s on socket %s\n", name.c_str(), address.c_str());
	// Read Simulation settings from onelab
	onelab_client = new onelab::remoteNetworkClient(name, address);
	onelab_client->sendInfo("qucs-onelab client connected");
	std::vector<onelab::string> OLResult;
	onelab_client->get(OLResult, "Circuit/Command to launch GetDP");
	cmd = OLResult.size() != 0 ? OLResult[0].getValue() : "";
	onelab_client->get(OLResult, "Circuit/Parameter for Flux");
	FluxParam = OLResult.size() != 0 ? OLResult[0].getValue() : "";
	onelab_client->get(OLResult, "Circuit/Parameter for differential Inductance");
	LdParam = OLResult.size() != 0 ? OLResult[0].getValue() : "";
	onelab_client->get(OLResult, "Circuit/Parameter to export current to");
	Iparam = OLResult.size() != 0 ? OLResult[0].getValue() : "";
	// setup logging if enabled
	onelab_client->get(OLResult, "Circuit/verbose");
	if (OLResult.size() != 0 && OLResult[0].getValue() == "yes"){
		verbose = true;
		std::stringstream msg;
		msg << "verbose output enabled:";
		msg << "realtime; ";
		msg << "it_duration; ";
		msg << "timestep; ";
		msg << "simulationtime; ";
		msg << "Ld12C; ";
		msg << "Flux12; ";
		msg << "I";
		onelab_client->sendInfo(msg.str());
	}
	// Disable Preinitialization with Field from last Run
	onelab::number InitFromRes("Input/42InitSolutionFromPrevious", 0);
	onelab_client->set(InitFromRes);
	// Clear Currents from Onelab
	std::vector<onelab::number> Is;
	onelab_client->get(Is,Iparam);
	onelab::number I = Is.size() ? Is[0] : onelab::number(Iparam,0);
	std::vector<double> choices = std::vector<double>();
	I.setChoices(choices);
	onelab_client->set(I);
}

qucsonelab4::~qucsonelab4 () {
	delete onelab_client;
}

void qucsonelab4::initDC (void) {
  setVoltageSources (2);
  allocMatrixMNA ();
  voltageSource (VSRC_1, NODE_1, NODE_4);
  voltageSource (VSRC_2, NODE_2, NODE_3);
}

void qucsonelab4::calcDC (void) {
  clearY ();
}

void qucsonelab4::initTR (void) {
  initDC ();
  clearY ();
  setStates (8);
}

#define fState11 0       // flux state
#define vState11 1       // voltage state
#define fState22 2
#define vState22 3
#define fState12 4
#define vState12 5
#define fState21 6
#define vState21 7

void qucsonelab4::calcTR (nr_double_t t) {
	time_t iterationtime = clock();
	std::vector<onelab::number> Is1;
	onelab_client->get(Is1,Iparam+"1");
	onelab::number I1 = Is1.size() ? Is1[0] : onelab::number(Iparam+"1",0);
	std::vector<onelab::number> Is2;
	onelab_client->get(Is2,Iparam+"2");
	onelab::number I2 = Is2.size() ? Is2[0] : onelab::number(Iparam+"2",0);
	// save time and current values for each timestep to onelab
	if (lastiterationtime != t) {
		std::vector<onelab::number> Ts;
		onelab_client->get(Ts,Tparam);
		onelab::number T = Ts.size() ? Ts[0] : onelab::number(Tparam,0);
		std::vector<double> choices = T.getChoices();
		choices.push_back(lastiterationtime);
		T.setChoices(choices);
		T.setValue(lastiterationtime);
		onelab_client->set(T);
		choices = I1.getChoices();
		choices.push_back(I1.getValue());
		I1.setChoices(choices);
		timestep++;
		lastiterationtime = t;
		std::stringstream msg;
		msg << "calcTR: Timestep: " << timestep << " @ t=" << t;
		onelab_client->sendInfo(msg.str());
	}
	double i1 = real(getJ(VSRC_1)); 
	double i2 = real(getJ(VSRC_2)); 
	I1.setValue(i1);
	onelab_client->set(I1);
	// run field computation and retrieve Inductance from Onelab
	onelab_client->runSubClient("getdp", cmd);
	// get Flux Phi
	std::vector<onelab::number> L;
	onelab_client->get(L, FluxParam+"1");
	nr_double_t l1;
	if (L.size()) {
		l1 = L[0].getValue();
	} else {
		l1 = 0;
		onelab_client->sendInfo("Error retreiving L1");
	}
	onelab_client->get(L, FluxParam+"2");
	nr_double_t l2;
	if (L.size()) {
		l2 = L[0].getValue();
	} else {
		l2 = 0;
		onelab_client->sendInfo("Error retreiving L2");
	}
	onelab_client->get(L, FluxParam+"C");
	nr_double_t lc;
	if (L.size()) {
		lc = L[0].getValue();
	} else {
		lc = 0;
		onelab_client->sendInfo("Error retreiving the common inductance");
	}
	// get differential inductance Ld = dPhi/di = L + Ldi * i approx Phi/i
	std::vector<onelab::number> Ld;
	onelab_client->get(Ld, LdParam+"1");
	nr_double_t ld1;
	if (Ld.size()) {
		ld1 = Ld[0].getValue();
	} else {
		ld1 = 0;
		onelab_client->sendInfo("Error retreiving Ld1");
	}
	onelab_client->get(Ld, LdParam+"2");
	nr_double_t ld2;
	if (Ld.size()) {
		ld2 = Ld[0].getValue();
	} else {
		ld2 = 0;
		onelab_client->sendInfo("Error retreiving Ld2");
	}
	onelab_client->get(Ld, LdParam+"C");
	nr_double_t ldc;
	if (Ld.size()) {
		ldc = Ld[0].getValue();
	} else {
		ldc = 0;
		onelab_client->sendInfo("Error retreiving Ldc");
	}
	if (verbose){
		std::stringstream msg;
		msg << "calcTR: ";
		msg	<< (float)(clock() - starttime)/CLOCKS_PER_SEC << "; ";  // Time since simulation start
		msg << (float)(clock()-iterationtime)/CLOCKS_PER_SEC << "; ";// Time for this iteration
		msg << timestep << "; ";                                     // timestep
		msg << t << "; ";                                            // simulation time
		msg << "[" << ld1 << ", " << ld2 << ", " << ldc << "]; ";    // Differential Inductance
		msg << "[" << l1 << ", " << l2 <<  ", " << lc <<  "]; ";     // Inductance 
		msg << "[" << i1 << ", " << i2 << "]\n";                     // Current
		onelab_client->sendInfo(msg.str());
	}
	// calculate matrices entries
  nr_double_t r11, r12, r21, r22, v11, v22, v12, v21;
  // self inductances
  setState  (fState11, l1*i1);
  integrate (fState11, ld1, r11, v11);
  setState  (fState22, l2*i2);
  integrate (fState22, ld2, r22, v22);

  // mutual inductances
  setState  (fState12, lc*i2);
  integrate (fState12, ldc, r12, v12);
  setState  (fState21, lc*i1);
  integrate (fState21, ldc, r21, v21);

  setD (VSRC_1, VSRC_1, -r11); setD (VSRC_1, VSRC_2, -r12);
  setD (VSRC_2, VSRC_1, -r21); setD (VSRC_2, VSRC_2, -r22);
  setE (VSRC_1, pol*(getState(vState11)-r11*i1 + getState(vState12)-r12*i1));
  setE (VSRC_2, pol*(getState(vState21)-r21*i2 + getState(vState22)-r22*i2));

	// Enable Preinitialization with Field from last Run
	onelab::number InitFromRes("Input/42InitSolutionFromPrevious", 1);
	onelab_client->set(InitFromRes);
}

// properties
PROP_REQ [] = {
	PROP_NO_PROP };
PROP_OPT [] = {
	{ "OLname", PROP_STR, { PROP_NO_VAL, "" }, PROP_NO_RANGE},
	{ "OLaddress", PROP_STR, { PROP_NO_VAL, "" }, PROP_NO_RANGE},
  PROP_NO_PROP };
struct define_t qucsonelab4::cirdef =
  { "OL4", 4, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_NONLINEAR, PROP_DEF };
