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

#ifndef _Charge_Parameters_
#define _Charge_Parameters_
#include <vector>
#include <string>

//descrition:
//parameters to auxiliate the calculation of charges of atoms in molecules

double CNFACGFN0(size_t atmnr) {
  //function returning the CNFAK parameters from GFN0
  double cnfak = 0.0;
  switch (atmnr) {
    case 1:
      cnfak = 0.04916110;
      break;
    case 2:
      cnfak = 0.10937243;
      break;
    case 3:
      cnfak = -0.12349591;
      break;
    case 4:
      cnfak = -0.02665108;
      break;
    case 5:
      cnfak = -0.02631658;
      break;
    case 6:
      cnfak = 0.06005196;
      break;
    case 7:
      cnfak = 0.09279548;
      break;
    case 8:
      cnfak = 0.11689703;
      break;
    case 9:
      cnfak = 0.15704746;
      break;
    case 10:
      cnfak = 0.07987901;
      break;
    case 11:
      cnfak = -0.10002962;
      break;
    case 12:
      cnfak = -0.07712863;
      break;
    case 13:
      cnfak = -0.02170561;
      break;
    case 14:
      cnfak = -0.04964052;
      break;
    case 15:
      cnfak = 0.14250599;
      break;
    case 16:
      cnfak = 0.07126660;
      break;
    case 17:
      cnfak = 0.13682750;
      break;
    case 18:
      cnfak = 0.14877121;
      break;
    case 19:
      cnfak = -0.10219289;
      break;
    case 20:
      cnfak = -0.08979338;
      break;
    case 21:
      cnfak = -0.08273597;
      break;
    case 22:
      cnfak = -0.01754829;
      break;
    case 23:
      cnfak = -0.02765460;
      break;
    case 24:
      cnfak = -0.02558926;
      break;
    case 25:
      cnfak = -0.08010286;
      break;
    case 26:
      cnfak = -0.04163215;
      break;
    case 27:
      cnfak = -0.09369631;
      break;
    case 28:
      cnfak = -0.03774117;
      break;
    case 29:
      cnfak = -0.05759708;
      break;
    case 30:
      cnfak = 0.02431998;
      break;
    case 31:
      cnfak = -0.01056270;
      break;
    case 32:
      cnfak = -0.02692862;
      break;
    case 33:
      cnfak = 0.07657769;
      break;
    case 34:
      cnfak = 0.06561608;
      break;
    case 35:
      cnfak = 0.08006749;
      break;
    case 36:
      cnfak = 0.14139200;
      break;
    case 37:
      cnfak = -0.05351029;
      break;
    case 38:
      cnfak = -0.06701705;
      break;
    case 39:
      cnfak = -0.07377246;
      break;
    case 40:
      cnfak = -0.02927768;
      break;
    case 41:
      cnfak = -0.03867291;
      break;
    case 42:
      cnfak = -0.06929825;
      break;
    case 43:
      cnfak = -0.04485293;
      break;
    case 44:
      cnfak = -0.04800824;
      break;
    case 45:
      cnfak = -0.01484022;
      break;
    case 46:
      cnfak = 0.07917502;
      break;
    case 47:
      cnfak = 0.06619243;
      break;
    case 48:
      cnfak = 0.02434095;
      break;
    case 49:
      cnfak = -0.01505548;
      break;
    case 50:
      cnfak = -0.03030768;
      break;
    case 51:
      cnfak = 0.01418235;
      break;
    case 52:
      cnfak = 0.08953411;
      break;
    case 53:
      cnfak = 0.08967527;
      break;
    case 54:
      cnfak = 0.07277771;
      break;
    case 55:
      cnfak = -0.02129476;
      break;
    case 56:
      cnfak = -0.06188828;
      break;
    case 57:
      cnfak = -0.06568203;
      break;
    case 58:
      cnfak = -0.11000000;
      break;
    case 59:
      cnfak = -0.11000000;
      break;
    case 60:
      cnfak = -0.11000000;
      break;
    case 61:
      cnfak = -0.11000000;
      break;
    case 62:
      cnfak = -0.11000000;
      break;
    case 63:
      cnfak = -0.11000000;
      break;
    case 64:
      cnfak = -0.11000000;
      break;
    case 65:
      cnfak = -0.11000000;
      break;
    case 66:
      cnfak = -0.11000000;
      break;
    case 67:
      cnfak = -0.11000000;
      break;
    case 68:
      cnfak = -0.11000000;
      break;
    case 69:
      cnfak = -0.11000000;
      break;
    case 70:
      cnfak = -0.11000000;
      break;
    case 71:
      cnfak = -0.11000000;
      break;
    case 72:
      cnfak = -0.03585873;
      break;
    case 73:
      cnfak = -0.03132400;
      break;
    case 74:
      cnfak = -0.05902379;
      break;
    case 75:
      cnfak = -0.02827592;
      break;
    case 76:
      cnfak = -0.07606260;
      break;
    case 77:
      cnfak = -0.02123839;
      break;
    case 78:
      cnfak = 0.03814822;
      break;
    case 79:
      cnfak = 0.02146834;
      break;
    case 80:
      cnfak = 0.01580538;
      break;
    case 81:
      cnfak = -0.00894298;
      break;
    case 82:
      cnfak = -0.05864876;
      break;
    case 83:
      cnfak = -0.01817842;
      break;
    case 84:
      cnfak = 0.07721851;
      break;
    case 85:
      cnfak = 0.07936083;
      break;
    case 86:
      cnfak = 0.05849285;
      break;
  }
  return cnfak;
}
double XIGFN0(size_t atmnr) {
  //function returning the xi parameters from GFN0
  double xi = 0.0;
  switch (atmnr) {
    case 1:
      xi = 1.23695041;
      break;
    case 2:
      xi = 1.26590957;
      break;
    case 3:
      xi = 0.54341808;
      break;
    case 4:
      xi = 0.99666991;
      break;
    case 5:
      xi = 1.26691604;
      break;
    case 6:
      xi = 1.40028282;
      break;
    case 7:
      xi = 1.55819364;
      break;
    case 8:
      xi = 1.56866440;
      break;
    case 9:
      xi = 1.57540015;
      break;
    case 10:
      xi = 1.15056627;
      break;
    case 11:
      xi = 0.55936220;
      break;
    case 12:
      xi = 0.72373742;
      break;
    case 13:
      xi = 1.12910844;
      break;
    case 14:
      xi = 1.12306840;
      break;
    case 15:
      xi = 1.52672442;
      break;
    case 16:
      xi = 1.40768172;
      break;
    case 17:
      xi = 1.48154584;
      break;
    case 18:
      xi = 1.31062963;
      break;
    case 19:
      xi = 0.40374140;
      break;
    case 20:
      xi = 0.75442607;
      break;
    case 21:
      xi = 0.76482096;
      break;
    case 22:
      xi = 0.98457281;
      break;
    case 23:
      xi = 0.96702598;
      break;
    case 24:
      xi = 1.05266584;
      break;
    case 25:
      xi = 0.93274875;
      break;
    case 26:
      xi = 1.04025281;
      break;
    case 27:
      xi = 0.92738624;
      break;
    case 28:
      xi = 1.07419210;
      break;
    case 29:
      xi = 1.07900668;
      break;
    case 30:
      xi = 1.04712861;
      break;
    case 31:
      xi = 1.15018618;
      break;
    case 32:
      xi = 1.15388455;
      break;
    case 33:
      xi = 1.36313743;
      break;
    case 34:
      xi = 1.36485106;
      break;
    case 35:
      xi = 1.39801837;
      break;
    case 36:
      xi = 1.18695346;
      break;
    case 37:
      xi = 0.36273870;
      break;
    case 38:
      xi = 0.58797255;
      break;
    case 39:
      xi = 0.71961946;
      break;
    case 40:
      xi = 0.96158233;
      break;
    case 41:
      xi = 0.89585296;
      break;
    case 42:
      xi = 0.81360499;
      break;
    case 43:
      xi = 1.00794665;
      break;
    case 44:
      xi = 0.92613682;
      break;
    case 45:
      xi = 1.09152285;
      break;
    case 46:
      xi = 1.14907070;
      break;
    case 47:
      xi = 1.13508911;
      break;
    case 48:
      xi = 1.08853785;
      break;
    case 49:
      xi = 1.11005982;
      break;
    case 50:
      xi = 1.12452195;
      break;
    case 51:
      xi = 1.21642129;
      break;
    case 52:
      xi = 1.36507125;
      break;
    case 53:
      xi = 1.40340000;
      break;
    case 54:
      xi = 1.16653482;
      break;
    case 55:
      xi = 0.34125098;
      break;
    case 56:
      xi = 0.58884173;
      break;
    case 57:
      xi = 0.68441115;
      break;
    case 58:
      xi = 0.56999999;
      break;
    case 59:
      xi = 0.56999999;
      break;
    case 60:
      xi = 0.56999999;
      break;
    case 61:
      xi = 0.56999999;
      break;
    case 62:
      xi = 0.56999999;
      break;
    case 63:
      xi = 0.56999999;
      break;
    case 64:
      xi = 0.56999999;
      break;
    case 65:
      xi = 0.56999999;
      break;
    case 66:
      xi = 0.56999999;
      break;
    case 67:
      xi = 0.56999999;
      break;
    case 68:
      xi = 0.56999999;
      break;
    case 69:
      xi = 0.56999999;
      break;
    case 70:
      xi = 0.56999999;
      break;
    case 71:
      xi = 0.56999999;
      break;
    case 72:
      xi = 0.87936784;
      break;
    case 73:
      xi = 1.02761808;
      break;
    case 74:
      xi = 0.93297476;
      break;
    case 75:
      xi = 1.10172128;
      break;
    case 76:
      xi = 0.97350071;
      break;
    case 77:
      xi = 1.16695666;
      break;
    case 78:
      xi = 1.23997927;
      break;
    case 79:
      xi = 1.18464453;
      break;
    case 80:
      xi = 1.14191734;
      break;
    case 81:
      xi = 1.12334192;
      break;
    case 82:
      xi = 1.01485321;
      break;
    case 83:
      xi = 1.12950808;
      break;
    case 84:
      xi = 1.30804834;
      break;
    case 85:
      xi = 1.33689961;
      break;
    case 86:
      xi = 1.27465977;
      break;
  }
  return xi;
}
double ALPGFN0(size_t atmnr) {
  //function returning the alp parameters from EN charges in xTB
  double alp = 0.0;
  switch (atmnr) {
    case 1:
      alp = 0.55159092;
      break;
    case 2:
      alp = 0.66205886;
      break;
    case 3:
      alp = 0.90529132;
      break;
    case 4:
      alp = 1.51710827;
      break;
    case 5:
      alp = 2.86070364;
      break;
    case 6:
      alp = 1.88862966;
      break;
    case 7:
      alp = 1.32250290;
      break;
    case 8:
      alp = 1.23166285;
      break;
    case 9:
      alp = 1.77503721;
      break;
    case 10:
      alp = 1.11955204;
      break;
    case 11:
      alp = 1.28263182;
      break;
    case 12:
      alp = 1.22344336;
      break;
    case 13:
      alp = 1.70936266;
      break;
    case 14:
      alp = 1.54075036;
      break;
    case 15:
      alp = 1.38200579;
      break;
    case 16:
      alp = 2.18849322;
      break;
    case 17:
      alp = 1.36779065;
      break;
    case 18:
      alp = 1.27039703;
      break;
    case 19:
      alp = 1.64466502;
      break;
    case 20:
      alp = 1.58859404;
      break;
    case 21:
      alp = 1.65357953;
      break;
    case 22:
      alp = 1.50021521;
      break;
    case 23:
      alp = 1.30104175;
      break;
    case 24:
      alp = 1.46301827;
      break;
    case 25:
      alp = 1.32928147;
      break;
    case 26:
      alp = 1.02766713;
      break;
    case 27:
      alp = 1.02291377;
      break;
    case 28:
      alp = 0.94343886;
      break;
    case 29:
      alp = 1.14881311;
      break;
    case 30:
      alp = 1.47080755;
      break;
    case 31:
      alp = 1.76901636;
      break;
    case 32:
      alp = 1.98724061;
      break;
    case 33:
      alp = 2.41244711;
      break;
    case 34:
      alp = 2.26739524;
      break;
    case 35:
      alp = 2.95378999;
      break;
    case 36:
      alp = 1.20807752;
      break;
    case 37:
      alp = 1.65941046;
      break;
    case 38:
      alp = 1.62733880;
      break;
    case 39:
      alp = 1.61344972;
      break;
    case 40:
      alp = 1.63220728;
      break;
    case 41:
      alp = 1.60899928;
      break;
    case 42:
      alp = 1.43501286;
      break;
    case 43:
      alp = 1.54559205;
      break;
    case 44:
      alp = 1.32663678;
      break;
    case 45:
      alp = 1.37644152;
      break;
    case 46:
      alp = 1.36051851;
      break;
    case 47:
      alp = 1.23395526;
      break;
    case 48:
      alp = 1.65734544;
      break;
    case 49:
      alp = 1.53895240;
      break;
    case 50:
      alp = 1.97542736;
      break;
    case 51:
      alp = 1.97636542;
      break;
    case 52:
      alp = 2.05432381;
      break;
    case 53:
      alp = 3.80138135;
      break;
    case 54:
      alp = 1.43893803;
      break;
    case 55:
      alp = 1.75505957;
      break;
    case 56:
      alp = 1.59815118;
      break;
    case 57:
      alp = 1.76401732;
      break;
    case 58:
      alp = 1.63999999;
      break;
    case 59:
      alp = 1.63999999;
      break;
    case 60:
      alp = 1.63999999;
      break;
    case 61:
      alp = 1.63999999;
      break;
    case 62:
      alp = 1.63999999;
      break;
    case 63:
      alp = 1.63999999;
      break;
    case 64:
      alp = 1.63999999;
      break;
    case 65:
      alp = 1.63999999;
      break;
    case 66:
      alp = 1.63999999;
      break;
    case 67:
      alp = 1.63999999;
      break;
    case 68:
      alp = 1.63999999;
      break;
    case 69:
      alp = 1.63999999;
      break;
    case 70:
      alp = 1.63999999;
      break;
    case 71:
      alp = 1.63999999;
      break;
    case 72:
      alp = 1.47055223;
      break;
    case 73:
      alp = 1.81127084;
      break;
    case 74:
      alp = 1.40189963;
      break;
    case 75:
      alp = 1.54015481;
      break;
    case 76:
      alp = 1.33721475;
      break;
    case 77:
      alp = 1.57165422;
      break;
    case 78:
      alp = 1.04815857;
      break;
    case 79:
      alp = 1.78342098;
      break;
    case 80:
      alp = 2.79106396;
      break;
    case 81:
      alp = 1.78160840;
      break;
    case 82:
      alp = 2.47588882;
      break;
    case 83:
      alp = 2.37670734;
      break;
    case 84:
      alp = 1.76613217;
      break;
    case 85:
      alp = 2.66172302;
      break;
    case 86:
      alp = 2.82773085;
      break;
  }
  return alp;
}
double GAMMAGFN0(size_t atmnr) {
  //function returning the gamma parameters from GFN0
  double gamma = 0.0;
  switch (atmnr) {
    case 1:
      gamma = -0.35015861;
      break;
    case 2:
      gamma = 1.04121227;
      break;
    case 3:
      gamma = 0.09281243;
      break;
    case 4:
      gamma = 0.09412380;
      break;
    case 5:
      gamma = 0.26629137;
      break;
    case 6:
      gamma = 0.19408787;
      break;
    case 7:
      gamma = 0.05317918;
      break;
    case 8:
      gamma = 0.03151644;
      break;
    case 9:
      gamma = 0.32275132;
      break;
    case 10:
      gamma = 1.30996037;
      break;
    case 11:
      gamma = 0.24206510;
      break;
    case 12:
      gamma = 0.04147733;
      break;
    case 13:
      gamma = 0.11634126;
      break;
    case 14:
      gamma = 0.13155266;
      break;
    case 15:
      gamma = 0.15350650;
      break;
    case 16:
      gamma = 0.15250997;
      break;
    case 17:
      gamma = 0.17523529;
      break;
    case 18:
      gamma = 0.28774450;
      break;
    case 19:
      gamma = 0.42937314;
      break;
    case 20:
      gamma = 0.01896455;
      break;
    case 21:
      gamma = 0.07179178;
      break;
    case 22:
      gamma = -0.01121381;
      break;
    case 23:
      gamma = -0.03093370;
      break;
    case 24:
      gamma = 0.02716319;
      break;
    case 25:
      gamma = -0.01843812;
      break;
    case 26:
      gamma = -0.15270393;
      break;
    case 27:
      gamma = -0.09192645;
      break;
    case 28:
      gamma = -0.13418723;
      break;
    case 29:
      gamma = -0.09861139;
      break;
    case 30:
      gamma = 0.18338109;
      break;
    case 31:
      gamma = 0.08299615;
      break;
    case 32:
      gamma = 0.11370033;
      break;
    case 33:
      gamma = 0.19005278;
      break;
    case 34:
      gamma = 0.10980677;
      break;
    case 35:
      gamma = 0.12327841;
      break;
    case 36:
      gamma = 0.25345554;
      break;
    case 37:
      gamma = 0.58615231;
      break;
    case 38:
      gamma = 0.16093861;
      break;
    case 39:
      gamma = 0.04548530;
      break;
    case 40:
      gamma = -0.02478645;
      break;
    case 41:
      gamma = 0.01909943;
      break;
    case 42:
      gamma = 0.01402541;
      break;
    case 43:
      gamma = -0.03595279;
      break;
    case 44:
      gamma = 0.01137752;
      break;
    case 45:
      gamma = -0.03697213;
      break;
    case 46:
      gamma = 0.08009416;
      break;
    case 47:
      gamma = 0.02274892;
      break;
    case 48:
      gamma = 0.12801822;
      break;
    case 49:
      gamma = -0.02078702;
      break;
    case 50:
      gamma = 0.05284319;
      break;
    case 51:
      gamma = 0.07581190;
      break;
    case 52:
      gamma = 0.09663758;
      break;
    case 53:
      gamma = 0.09547417;
      break;
    case 54:
      gamma = 0.07803344;
      break;
    case 55:
      gamma = 0.64913257;
      break;
    case 56:
      gamma = 0.15348654;
      break;
    case 57:
      gamma = 0.05054344;
      break;
    case 58:
      gamma = 0.11000000;
      break;
    case 59:
      gamma = 0.11000000;
      break;
    case 60:
      gamma = 0.11000000;
      break;
    case 61:
      gamma = 0.11000000;
      break;
    case 62:
      gamma = 0.11000000;
      break;
    case 63:
      gamma = 0.11000000;
      break;
    case 64:
      gamma = 0.11000000;
      break;
    case 65:
      gamma = 0.11000000;
      break;
    case 66:
      gamma = 0.11000000;
      break;
    case 67:
      gamma = 0.11000000;
      break;
    case 68:
      gamma = 0.11000000;
      break;
    case 69:
      gamma = 0.11000000;
      break;
    case 70:
      gamma = 0.11000000;
      break;
    case 71:
      gamma = 0.11000000;
      break;
    case 72:
      gamma = -0.02786741;
      break;
    case 73:
      gamma = 0.01057858;
      break;
    case 74:
      gamma = -0.03892226;
      break;
    case 75:
      gamma = -0.04574364;
      break;
    case 76:
      gamma = -0.03874080;
      break;
    case 77:
      gamma = -0.03782372;
      break;
    case 78:
      gamma = -0.07046855;
      break;
    case 79:
      gamma = 0.09546597;
      break;
    case 80:
      gamma = 0.21953269;
      break;
    case 81:
      gamma = 0.02522348;
      break;
    case 82:
      gamma = 0.15263050;
      break;
    case 83:
      gamma = 0.08042611;
      break;
    case 84:
      gamma = 0.01878626;
      break;
    case 85:
      gamma = 0.08715453;
      break;
    case 86:
      gamma = 0.10500484;
      break;
  }
  return gamma;
}
double CM5alpha() {
  //alpha for CM5 in Angstroem^-1
  //A. V. Marenich, S. V. Jerome, C. J. Cramer, D. G. Truhlar, J. Chem, Theo. Comput, 8, 527, 2012
  return 2.474;
}
double CM5DZ(int atomnumb) {
  //DZ parameters required for CM5 (part of table 1); C coefficient already included in higher elements
  //A. V. Marenich, S. V. Jerome, C. J. Cramer, D. G. Truhlar, J. Chem, Theo. Comput, 8, 527, 2012
  double cm5dz = 0.0;
  if (atomnumb == 1) {cm5dz = 0.0056;}
  else if (atomnumb == 2) {cm5dz = -0.1543;}
  else if (atomnumb == 4) {cm5dz = 0.0333;}
  else if (atomnumb == 5) {cm5dz = -0.1030;}
  else if (atomnumb == 6) {cm5dz = -0.0446;}
  else if (atomnumb == 7) {cm5dz = -0.1072;}
  else if (atomnumb == 8) {cm5dz = -0.0802;}
  else if (atomnumb == 9) {cm5dz = -0.0629;}
  else if (atomnumb == 10) {cm5dz = -0.1088;}
  else if (atomnumb == 11) {cm5dz = 0.0184;}
  else if (atomnumb == 13) {cm5dz = -0.0726;}
  else if (atomnumb == 14) {cm5dz = -0.0790;}
  else if (atomnumb == 15) {cm5dz = -0.0756;}
  else if (atomnumb == 16) {cm5dz = -0.0565;}
  else if (atomnumb == 17) {cm5dz = -0.0444;}
  else if (atomnumb == 18) {cm5dz = -0.0767;}
  else if (atomnumb == 19) {cm5dz = 0.0130;}
  else if (atomnumb == 31) {cm5dz = -0.0512;}
  else if (atomnumb == 32) {cm5dz = -0.0557;}
  else if (atomnumb == 33) {cm5dz = -0.0533;}
  else if (atomnumb == 34) {cm5dz = -0.0399;}
  else if (atomnumb == 35) {cm5dz = -0.0313;}
  else if (atomnumb == 36) {cm5dz = -0.0541;}
  else if (atomnumb == 37) {cm5dz = 0.0092;}
  else if (atomnumb == 49) {cm5dz = -0.0361;}
  else if (atomnumb == 50) {cm5dz = -0.0393;}
  else if (atomnumb == 51) {cm5dz = -0.0376;}
  else if (atomnumb == 52) {cm5dz = -0.0281;}
  else if (atomnumb == 53) {cm5dz = -0.0220;}
  else if (atomnumb == 54) {cm5dz = -0.0381;}
  else if (atomnumb == 55) {cm5dz = 0.0065;}
  else if (atomnumb == 81) {cm5dz = -0.0255;}
  else if (atomnumb == 82) {cm5dz = -0.0277;}
  else if (atomnumb == 83) {cm5dz = -0.0265;}
  else if (atomnumb == 84) {cm5dz = -0.0198;}
  else if (atomnumb == 85) {cm5dz = -0.0155;}
  else if (atomnumb == 86) {cm5dz = -0.0269;}
  else if (atomnumb == 87) {cm5dz = 0.0046;}
  else if (atomnumb == 113) {cm5dz = -0.0179;}
  else if (atomnumb == 114) {cm5dz = -0.0195;}
  else if (atomnumb == 115) {cm5dz = -0.0187;}
  else if (atomnumb == 116) {cm5dz = -0.0140;}
  else if (atomnumb == 117) {cm5dz = -0.0110;}
  else if (atomnumb == 118) {cm5dz = -0.0189;}
  return cm5dz;
}
double TkkCM5(size_t atom1, size_t atom2) {
  //function returning the Tkk' parameter in CM5, given by equations 3, 4 and table 1
  //atom1 -> k         atom2 -> k'
  //A. V. Marenich, S. V. Jerome, C. J. Cramer, D. G. Truhlar, J. Chem, Theo. Comput, 8, 527, 2012
  size_t atomA = std::min(atom1,atom2);
  size_t atomB = std::max(atom1,atom2);
  double factor = 1.0 - 2.0*double(atom2 == atomA);
  double tkkcm5 = 0.0;
  if ((atomA == 1)&&(atomB == 6)) {tkkcm5 = factor*0.0502;}           //equation 3
  else if ((atomA == 1)&&(atomB == 7)) {tkkcm5 = factor*0.1747;}      //equation 3
  else if ((atomA == 1)&&(atomB == 8)) {tkkcm5 = factor*0.1671;}      //equation 3
  else if ((atomA == 6)&&(atomB == 7)) {tkkcm5 = factor*0.0556;}      //equation 3
  else if ((atomA == 6)&&(atomB == 8)) {tkkcm5 = factor*0.0234;}      //equation 3
  else if ((atomA == 7)&&(atomB == 8)) {tkkcm5 = -factor*0.0346;}     //equation 3
  else {tkkcm5 = CM5DZ(atom1) - CM5DZ(atom2);}                        //equation 4
  return tkkcm5;
}

#endif //_Charge_Parameters_
