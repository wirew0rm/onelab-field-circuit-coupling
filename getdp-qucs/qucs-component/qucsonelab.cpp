/*
 * qucsonelab.cpp - implements a dynamic qucs component, that interfaces with
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

#include "qucsonelab.h"

using namespace qucs;

static double pi = 3;

qucsonelab::qucsonelab () : circuit (2) {
	starttime = clock(); // Start timing the evaluation of the circuit solver
	logprint(LOG_ERROR, "Initializing ONELAB component\n");
  type = CIR_ECVS;
  setISource (true);
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
		msg << "Ld; ";
		msg << "Flux; ";
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

qucsonelab::~qucsonelab () {
	delete onelab_client;
}

void qucsonelab::initDC (void) {
	logprint(LOG_ERROR, "initDC\n");
  setVoltageSources (1);
  allocMatrixMNA ();
  voltageSource (VSRC_1, NODE_1, NODE_2);
}

void qucsonelab::calcDC (void) {
	logprint(LOG_ERROR, "calcDC\n");
  clearY ();
}

void qucsonelab::initTR (void) {
	logprint(LOG_ERROR, "initTR\n");
  initDC ();
  clearY ();
  setStates (2);
}

#define fState 0       // flux state
#define vState 1       // voltage state
#define i_min 2e-10    // minimal excitation
#define i_delta 1e-5   // relative h for difference quotient

void qucsonelab::calcTR (nr_double_t t) {
	time_t iterationtime = clock();
	std::vector<onelab::number> Is;
	onelab_client->get(Is,Iparam);
	onelab::number I = Is.size() ? Is[0] : onelab::number(Iparam,0);
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
		choices = I.getChoices();
		choices.push_back(I.getValue());
		I.setChoices(choices);
		timestep++;
		lastiterationtime = t;
		std::stringstream msg;
		msg << "calcTR: Timestep: " << timestep << " @ t=" << t;
		onelab_client->sendInfo(msg.str());
	}
	double i = real(getJ(VSRC_1)); 
	double i_ol = i;
  /* apply initial condition if requested */
  if (getMode () == MODE_INIT && isPropertyGiven ("I")) {
    i = getPropertyDouble ("I");
		i_ol = i;
  }
	// Set excitation with a minimal excitation of 1e-10 to be able to calculate L
	// only if L = Lchord
	if (LdParam == "Lchord" || LdParam == "diffquot") {
		if (0 <= i && i < i_min) {
			i_ol = i_min;
		} else if (0 > i && i > -i_min) {
			i_ol = -i_min;
		}
	}
	I.setValue(i_ol);
	onelab_client->set(I);
	// run field computation and retrieve Inductance from Onelab
	onelab_client->runSubClient("getdp", cmd);
	// get Flux Phi
	std::vector<onelab::number> Flux;
	onelab_client->get(Flux, FluxParam);
	nr_double_t flux;
	if (Flux.size()) {
		flux= Flux[0].getValue();
	} else {
		flux= 0;
		onelab_client->sendInfo("Error retreiving Flux");
	}
	// get differential inductance Ld = dPhi/di = L + Ldi * i approx Phi/i
	nr_double_t ld = 0;
	if (LdParam == "diffquot") {
		onelab::number Idelta = onelab::number(Iparam, i_ol*(1+i_delta));
		onelab_client->set(Idelta);
		onelab_client->runSubClient("getdp", cmd);
		std::vector<onelab::number> FluxDelta;
		onelab_client->get(FluxDelta, FluxParam);
		nr_double_t fluxdelta;
		if (FluxDelta.size()) {
			fluxdelta = FluxDelta[0].getValue();
		} else {
			fluxdelta = 0;
			onelab_client->sendInfo("Error retreiving FluxDelta");
		}
		ld = (fluxdelta - flux)/i_delta/i_ol;
		// compensate for the minimal excitation
		flux = flux * (i / i_ol);
	} else if (LdParam == "Lchord") {
		ld = flux/i_ol;
		// compensate for the minimal excitation
		flux = flux * (i / i_ol);
	} else {
		std::vector<onelab::number> Ld;
		onelab_client->get(Ld, LdParam);
		if (Ld.size()) {
			ld = Ld[0].getValue();
		} else {
			ld = flux/I.getValue();
			onelab_client->sendInfo("Error retreiving Ld, fallback to Lchord");
		}
	}
	if (verbose){
		std::stringstream msg;
		msg << "calcTR: ";
		msg	<< (float)(clock() - starttime)/CLOCKS_PER_SEC << "; ";  // Time since simulation start
		msg << (float)(clock()-iterationtime)/CLOCKS_PER_SEC << "; ";// Time for this iteration
		msg << timestep << "; "; // timestep
		msg << t << "; ";        // simulation time
		msg << ld << "; ";       // Differential Inductance
		msg << flux << "; ";     // Flux
		msg << i << "\n";        // Current
		onelab_client->sendInfo(msg.str());
	}
	// calculate matrices entries
  nr_double_t r, v;
	setState (fState, flux);
	integrate(fState, ld, r, v);
	setD(VSRC_1, VSRC_1, -r);
	v = pol * (getState(vState) - r * i);
	setE(VSRC_1, v);

	// Enable Preinitialization with Field from last Run
	onelab::number InitFromRes("Input/42InitSolutionFromPrevious", 1);
	onelab_client->set(InitFromRes);
}

// properties
PROP_REQ [] = {
	PROP_NO_PROP };
PROP_OPT [] = {
  { "I", PROP_REAL, { 0, PROP_NO_STR }, PROP_NO_RANGE },
	{ "OLname", PROP_STR, { PROP_NO_VAL, "" }, PROP_NO_RANGE},
	{ "OLaddress", PROP_STR, { PROP_NO_VAL, "" }, PROP_NO_RANGE},
  PROP_NO_PROP };
struct define_t qucsonelab::cirdef =
  { "OL", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_NONLINEAR, PROP_DEF };
