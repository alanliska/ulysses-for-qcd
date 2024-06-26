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

import scipy.special
import sys

#description:
#associated Legendre polynomials of first kind calculated via python

def legendrepmn(m,n,x):
    P_lmp, dP_lmp = scipy.special.lpmn(m,n,x)
    result = float(P_lmp[m][n])
    return result

m = int(sys.argv[1])
n = int(sys.argv[2])
x = float(sys.argv[3])
print legendrepmn(m,n,x)
