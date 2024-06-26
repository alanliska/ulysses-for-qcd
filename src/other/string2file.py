'''ULYSSES, a semi-empirical package
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
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA'''

import os, sys

#description:
#python function writing string to file

readfile = sys.argv[1]
output_dir = sys.argv[2]
string = sys.argv[3]
initiliaze = sys.argv[4]

if not os.path.isdir(output_dir): os.mkdir(output_dir)

mode = "a+"
if initiliaze == "true": mode = "w+"
fileread = open(output_dir+"/"+readfile,mode)
fileread.write(string)

fileread.close()