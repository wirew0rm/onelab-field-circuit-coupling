# Coupling of Electromagnetic Fields with Electric Circuits Using Onelab

This repo contains the code for my 2015 BSc Thesis on strong field cirucit coupling.
The following tools where used:
- [QUCS](https://qucs.sourceforge.net): Quite universal circuit simulator. This repo contains code for a shared library which can be used to exchange data with a onelab server.
- [OCS](https://wiki.octave.org/Ocs_package): the octave circuit simulator.
  - [octave-onlab-interface](https://github.com/wirew0rm/octave-onelab-interface) is used for comunicating with getdp.
- [getdp](https://getdp.info): A finite element solver used for the field simulation part.
- [onelab](https://onelab.info): A framework for passing simulation data between different solvers/meshers/etc.

# License

All code in this repository is licensed as GPLv3+

    Copyright (C) 2015 Alexander Krimm <alex@wirew0rm.de>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

The thesis contained in the repository is licensed as CreativeCommons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0) 
for details see: https://creativecommons.org/licenses/by-sa/3.0/
