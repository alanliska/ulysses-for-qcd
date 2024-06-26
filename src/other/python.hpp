/*ULYSSES, a semi-empirical package
    Copyright (C) 2023- Filipe Menezes (filipe.menezes@helmholtz-munich.de)
                        Grzegorz Popowicz (grzegorz.popowicz@helmholtz-munich.de)

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software 
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA*/

#ifndef _Python_Package_
#define _Python_Package_
#include <vector>
#include <fstream>
#include <iostream>

//description:
//interface to some python functions

void SnapshotBox(size_t Nparticles, size_t ndim, double sigma, bool monochromatic, double timestamp, std::vector<std::vector<double> > * positions, std::vector<std::vector<double> > * velocities, std::string datafile, std::string folder) {
  //function to interface with python pylab to print snapshots of spheres in box
  //write to file
  std::ofstream myfile;
  std::cout << "t = " << timestamp << std::endl;
  myfile.open(folder+"/"+datafile);
  myfile << Nparticles << std::endl;
  myfile << sigma << std::endl;
  for (size_t idpart = 0; idpart < Nparticles; ++idpart) {
    myfile << idpart + 1 << " ; t = " << timestamp << ";  r = ";
	  for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
	    myfile << positions->at(idpart)[idcoord];
	    if (idcoord != ndim - 1) {myfile << ",";}
	  }
	  myfile << std::endl;
	  myfile << idpart + 1 << " ; t = " << timestamp << ";  v = ";
	  for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
	    myfile << velocities->at(idpart)[idcoord];
	    if (idcoord != ndim - 1) {myfile << ",";}
	  }
	  myfile << std::endl;
	  myfile << "-------" << std::endl;
  }
  myfile << "===========" << std::endl;
  myfile.close();
  //call python
  std::string command = "python other/event_snapshot.py "+folder+"/"+datafile+" "+folder;
  FILE * stream = popen(command.c_str(),"r");
  pclose(stream);
}
void SnapEvents(std::string datafile, std::string folder) {
  std::string command = "python other/event_snapshot.py "+folder+"/"+datafile+" "+folder;
  FILE * stream = popen(command.c_str(),"r");
  pclose(stream);
}

#endif //_Python_Package_
