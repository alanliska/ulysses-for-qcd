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

#ifndef _Global_Variables_
#define _Global_Variables_
#include <vector>

//description:
//global functions that determine some options
//not really in use

size_t size_print = 1;
std::vector<size_t> Prints(size_print,0);         //vector allocating the printing at some parts of the code. It works with >= condition
                                                  //0 = Symmetry functions: 1 - basic print with number of elements per symmetry operation
                                                  //                        2 - process to determine point group
                                                  //                        3 - detected axes and plane vectors
                                                  //                        4 - loop search information
                                                  //                        5 - comparison information with axes found

#endif //_Global_Variables_
