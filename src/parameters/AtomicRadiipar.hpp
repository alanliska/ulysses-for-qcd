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

#ifndef _Atomic_Radii_Parameters_
#define _Atomic_Radii_Parameters_

//description:
//atomic radii

double AtmRadii(int atmnumber) {
  //atomic radii of the elements in Angstroem according to
  //M. Mantina, R. Valero, C. J. Cramer, D. G. Truhlar, CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011), edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50; corrected Nov. 17, 2010 for the 92nd edition.
  double atmr = 0.0;
  switch (atmnumber) {
    case 1:
      atmr = 0.32;
      break;
    case 2:
      atmr = 0.37;
      break;
    case 3:
      atmr = 1.30;
      break;
    case 4:
      atmr = 0.99;
      break;
    case 5:
      atmr = 0.84;
      break;
    case 6:
      atmr = 0.75;
      break;
    case 7:
      atmr = 0.71;
      break;
    case 8:
      atmr = 0.64;
      break;
    case 9:
      atmr = 0.60;
      break;
    case 10:
      atmr = 0.62;
      break;
    case 11:
      atmr = 1.60;
      break;
    case 12:
      atmr = 1.40;
      break;
    case 13:
      atmr = 1.24;
      break;
    case 14:
      atmr = 1.14;
      break;
    case 15:
      atmr = 1.09;
      break;
    case 16:
      atmr = 1.04;
      break;
    case 17:
      atmr = 1.00;
      break;
    case 18:
      atmr = 1.01;
      break;
    case 19:
      atmr = 2.00;
      break;
    case 20:
      atmr = 1.74;
      break;
    case 21:
      atmr = 1.59;
      break;
    case 22:
      atmr = 1.48;
      break;
    case 23:
      atmr = 1.44;
      break;
    case 24:
      atmr = 1.30;
      break;
    case 25:
      atmr = 1.29;
      break;
    case 26:
      atmr = 1.24;
      break;
    case 27:
      atmr = 1.18;
      break;
    case 28:
      atmr = 1.17;
      break;
    case 29:
      atmr = 1.22;
      break;
    case 30:
      atmr = 1.20;
      break;
    case 31:
      atmr = 1.23;
      break;
    case 32:
      atmr = 1.20;
      break;
    case 33:
      atmr = 1.20;
      break;
    case 34:
      atmr = 1.18;
      break;
    case 35:
      atmr = 1.17;
      break;
    case 36:
      atmr = 1.16;
      break;
    case 37:
      atmr = 2.15;
      break;
    case 38:
      atmr = 1.90;
      break;
    case 39:
      atmr = 1.76;
      break;
    case 40:
      atmr = 1.64;
      break;
    case 41:
      atmr = 1.56;
      break;
    case 42:
      atmr = 1.46;
      break;
    case 43:
      atmr = 1.38;
      break;
    case 44:
      atmr = 1.36;
      break;
    case 45:
      atmr = 1.34;
      break;
    case 46:
      atmr = 1.30;
      break;
    case 47:
      atmr = 1.36;
      break;
    case 48:
      atmr = 1.40;
      break;
    case 49:
      atmr = 1.42;
      break;
    case 50:
      atmr = 1.40;
      break;
    case 51:
      atmr = 1.40;
      break;
    case 52:
      atmr = 1.37;
      break;
    case 53:
      atmr = 1.36;
      break;
    case 54:
      atmr = 1.36;
      break;
    case 55:
      atmr = 2.38;
      break;
    case 56:
      atmr = 2.06;
      break;
    case 57:
      atmr = 1.94;
      break;
    case 58:
      atmr = 1.84;
      break;
    case 59:
      atmr = 1.90;
      break;
    case 60:
      atmr = 1.88;
      break;
    case 61:
      atmr = 1.86;
      break;
    case 62:
      atmr = 1.85;
      break;
    case 63:
      atmr = 1.83;
      break;
    case 64:
      atmr = 1.82;
      break;
    case 65:
      atmr = 1.81;
      break;
    case 66:
      atmr = 1.80;
      break;
    case 67:
      atmr = 1.79;
      break;
    case 68:
      atmr = 1.77;
      break;
    case 69:
      atmr = 1.77;
      break;
    case 70:
      atmr = 1.78;
      break;
    case 71:
      atmr = 1.74;
      break;
    case 72:
      atmr = 1.64;
      break;
    case 73:
      atmr = 1.58;
      break;
    case 74:
      atmr = 1.50;
      break;
    case 75:
      atmr = 1.41;
      break;
    case 76:
      atmr = 1.36;
      break;
    case 77:
      atmr = 1.32;
      break;
    case 78:
      atmr = 1.30;
      break;
    case 79:
      atmr = 1.30;
      break;
    case 80:
      atmr = 1.32;
      break;
    case 81:
      atmr = 1.44;
      break;
    case 82:
      atmr = 1.45;
      break;
    case 83:
      atmr = 1.50;
      break;
    case 84:
      atmr = 1.42;
      break;
    case 85:
      atmr = 1.48;
      break;
    case 86:
      atmr = 1.46;
      break;
    case 87:
      atmr = 2.42;
      break;
    case 88:
      atmr = 2.11;
      break;
    case 89:
      atmr = 2.01;
      break;
    case 90:
      atmr = 1.90;
      break;
    case 91:
      atmr = 1.84;
      break;
    case 92:
      atmr = 1.83;
      break;
    case 93:
      atmr = 1.80;
      break;
    case 94:
      atmr = 1.80;
      break;
    case 95:
      atmr = 1.73;
      break;
    case 96:
      atmr = 1.68;
      break;
    case 97:
      atmr = 1.68;
      break;
    case 98:
      atmr = 1.68;
      break;
    case 99:
      atmr = 1.65;
      break;
    case 100:
      atmr = 1.67;
      break;
    case 101:
      atmr = 1.73;
      break;
    case 102:
      atmr = 1.76;
      break;
    case 103:
      atmr = 1.61;
      break;
    case 104:
      atmr = 1.57;
      break;
    case 105:
      atmr = 1.49;
      break;
    case 106:
      atmr = 1.43;
      break;
    case 107:
      atmr = 1.41;
      break;
    case 108:
      atmr = 1.34;
      break;
    case 109:
      atmr = 1.29;
      break;
    case 110:
      atmr = 1.28;
      break;
    case 111:
      atmr = 1.21;
      break;
    case 112:
      atmr = 1.22;
      break;
    case 113:
      atmr = 1.36;
      break;
    case 114:
      atmr = 1.43;
      break;
    case 115:
      atmr = 1.62;
      break;
    case 116:
      atmr = 1.75;
      break;
    case 117:
      atmr = 1.65;
      break;
    case 118:
      atmr = 1.57;
      break;
  }
  return atmr;
}
double vdWaalsRadii(size_t atmnumber) {
  //van der Waals radii of the elements in Angstroem
  double vdwrad = 0.0;
  switch (atmnumber) {
    case 1:
      vdwrad = 1.09155;
      break;
    case 2:
      vdwrad = 0.86735;
      break;
    case 3:
      vdwrad = 1.74780;
      break;
    case 4:
      vdwrad = 1.54910;
      break;
    case 5:
      vdwrad = 1.60800;
      break;
    case 6:
      vdwrad = 1.45515;
      break;
    case 7:
      vdwrad = 1.31125;
      break;
    case 8:
      vdwrad = 1.24085;
      break;
    case 9:
      vdwrad = 1.14980;
      break;
    case 10:
      vdwrad = 1.06870;
      break;
    case 11:
      vdwrad = 1.85410;
      break;
    case 12:
      vdwrad = 1.74195;
      break;
    case 13:
      vdwrad = 2.00530;
      break;
    case 14:
      vdwrad = 1.89585;
      break;
    case 15:
      vdwrad = 1.75085;
      break;
    case 16:
      vdwrad = 1.65535;
      break;
    case 17:
      vdwrad = 1.55230;
      break;
    case 18:
      vdwrad = 1.45740;
      break;
    case 19:
      vdwrad = 2.12055;
      break;
    case 20:
      vdwrad = 2.05175;
      break;
    case 21:
      vdwrad = 1.94515;
      break;
    case 22:
      vdwrad = 1.88210;
      break;
    case 23:
      vdwrad = 1.86055;
      break;
    case 24:
      vdwrad = 1.72070;
      break;
    case 25:
      vdwrad = 1.77310;
      break;
    case 26:
      vdwrad = 1.72105;
      break;
    case 27:
      vdwrad = 1.71635;
      break;
    case 28:
      vdwrad = 1.67310;
      break;
    case 29:
      vdwrad = 1.65040;
      break;
    case 30:
      vdwrad = 1.61545;
      break;
    case 31:
      vdwrad = 1.97895;
      break;
    case 32:
      vdwrad = 1.93095;
      break;
    case 33:
      vdwrad = 1.83125;
      break;
    case 34:
      vdwrad = 1.76340;
      break;
    case 35:
      vdwrad = 1.68310;
      break;
    case 36:
      vdwrad = 1.60480;
      break;
    case 37:
      vdwrad = 2.30880;
      break;
    case 38:
      vdwrad = 2.23820;
      break;
    case 39:
      vdwrad = 2.10980;
      break;
    case 40:
      vdwrad = 2.02985;
      break;
    case 41:
      vdwrad = 1.92980;
      break;
    case 42:
      vdwrad = 1.87715;
      break;
    case 43:
      vdwrad = 1.78450;
      break;
    case 44:
      vdwrad = 1.73115;
      break;
    case 45:
      vdwrad = 1.69875;
      break;
    case 46:
      vdwrad = 1.67625;
      break;
    case 47:
      vdwrad = 1.66540;
      break;
    case 48:
      vdwrad = 1.73100;
      break;
    case 49:
      vdwrad = 2.13115;
      break;
    case 50:
      vdwrad = 2.09370;
      break;
    case 51:
      vdwrad = 2.00750;
      break;
    case 52:
      vdwrad = 1.94505;
      break;
    case 53:
      vdwrad = 1.86900;
      break;
    case 54:
      vdwrad = 1.79445;
      break;
    case 55:
      vdwrad = 2.52835;
      break;
    case 56:
      vdwrad = 2.59070;
      break;
    case 57:
      vdwrad = 2.31305;
      break;
    case 58:
      vdwrad = 2.31005;
      break;
    case 59:
      vdwrad = 2.28510;
      break;
    case 60:
      vdwrad = 2.26355;
      break;
    case 61:
      vdwrad = 2.24480;
      break;
    case 62:
      vdwrad = 2.22575;
      break;
    case 63:
      vdwrad = 2.21170;
      break;
    case 64:
      vdwrad = 2.06215;
      break;
    case 65:
      vdwrad = 2.12135;
      break;
    case 66:
      vdwrad = 2.07705;
      break;
    case 67:
      vdwrad = 2.13970;
      break;
    case 68:
      vdwrad = 2.12250;
      break;
    case 69:
      vdwrad = 2.11040;
      break;
    case 70:
      vdwrad = 2.09930;
      break;
    case 71:
      vdwrad = 2.00650;
      break;
    case 72:
      vdwrad = 2.12250;
      break;
    case 73:
      vdwrad = 2.04900;
      break;
    case 74:
      vdwrad = 1.99275;
      break;
    case 75:
      vdwrad = 1.94775;
      break;
    case 76:
      vdwrad = 1.87450;
      break;
    case 77:
      vdwrad = 1.72280;
      break;
    case 78:
      vdwrad = 1.67625;
      break;
    case 79:
      vdwrad = 1.62820;
      break;
    case 80:
      vdwrad = 1.67995;
      break;
    case 81:
      vdwrad = 2.15635;
      break;
    case 82:
      vdwrad = 2.13820;
      break;
    case 83:
      vdwrad = 2.05875;
      break;
    case 84:
      vdwrad = 2.00270;
      break;
    case 85:
      vdwrad = 1.93220;
      break;
    case 86:
      vdwrad = 1.86080;
      break;
    case 87:
      vdwrad = 2.53980;
      break;
    case 88:
      vdwrad = 2.46470;
      break;
    case 89:
      vdwrad = 2.35215;
      break;
    case 90:
      vdwrad = 2.21260;
      break;
    case 91:
      vdwrad = 2.22970;
      break;
    case 92:
      vdwrad = 2.19785;
      break;
    case 93:
      vdwrad = 2.17695;
      break;
    case 94:
      vdwrad = 2.21705;
      break;
  }
  return vdwrad;
}
#endif //_Atomic_Radii_Parameters_
