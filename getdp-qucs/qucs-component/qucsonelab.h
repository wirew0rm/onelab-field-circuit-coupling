/*
 * inductor.h - inductor class definitions
 *
 * 2015, Alexander Krimm <alexander_johannes.krimm@stud.tu-darmstadt.de>
 * Copyright (C) 2003, 2004, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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
 * $Id$
 *
 */

#ifndef __QUCSONELAB_H__
#define __QUCSONELAB_H__
#include <gmsh/onelab.h>
#include <time.h>

class qucsonelab : public qucs::circuit
{
 public:
  CREATOR (qucsonelab);
	~qucsonelab ();
  void initDC (void);
  void calcDC (void);
  void initTR (void);
  void calcTR (nr_double_t);
 private:
	onelab::remoteNetworkClient *onelab_client;
	std::string cmd;
	std::string Iparam;
	std::string FluxParam;
	std::string LdParam;
	std::string Tparam = std::string("Circuit/55simulationtime");
	bool verbose;
	int timestep = 0;
	// timer to check
	clock_t starttime;
	// time of the last iteration to determine if we're in a new timestep
	nr_double_t lastiterationtime;
};

extern "C" {
	class va_proxy {
		public:
			va_proxy(){
				factorycreate["qucsonelab"] = qucsonelab::create;
				factorydef   ["qucsonelab"] = qucsonelab::definition;
			}
	};
	static va_proxy proxy_instance;
}

#endif /* __QUCSONELAB_H__ */
