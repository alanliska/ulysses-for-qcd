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

#ifndef _Isotope_Weight_Package_
#define _Isotope_Weight_Package_
#include <vector>

//descrition:
//base class containing weights of atom isotopes

class IsotopeWeights {
  std::vector<std::vector<double> > Weights;
public:
  IsotopeWeights() {
    //weights in kg/mol
    std::vector<double> temp(3,0.0);
    temp[0] = 1.007825032e-3;
    temp[1] = 2.014101778e-3;
    temp[2] = 3.016049278e-3;
    Weights.push_back(temp);          //H
    temp.resize(2);
    temp[0] = 3.016029320e-3;
    temp[1] = 4.002603254e-3;
    Weights.push_back(temp);          //He
    temp[0] = 6.015122887e-3;
    temp[1] = 7.016003437e-3;
    Weights.push_back(temp);          //Li 
    temp.resize(1);
    temp[0] = 9.012183066e-3;
    Weights.push_back(temp);          //Be 
    temp.resize(2);
    temp[0] = 10.01293695e-3;
    temp[1] = 11.00930536e-3;
    Weights.push_back(temp);          //B  
    temp.resize(3);
    temp[0] = 12.00000000e-3;
    temp[1] = 13.00335484e-3;
    temp[2] = 14.00324199e-3;
    Weights.push_back(temp);          //C  
    temp.resize(2);
    temp[0] = 14.00307400e-3;
    temp[1] = 15.00010890e-3;
    Weights.push_back(temp);          //N  
    temp.resize(3);
    temp[0] = 15.99491462e-3;
    temp[1] = 16.99913176e-3;
    temp[2] = 17.99915961e-3;
    Weights.push_back(temp);          //O  
    temp.resize(1);
    temp[0] = 18.99840316e-3;
    Weights.push_back(temp);          //F  
    temp.resize(3);
    temp[0] = 19.99244018e-3;
    temp[1] = 20.99384669e-3;
    temp[2] = 21.99138511e-3;
    Weights.push_back(temp);          //Ne 
    temp.resize(1);
    temp[0] = 22.98976928e-3;
    Weights.push_back(temp);          //Na 
    temp.resize(3);
    temp[0] = 23.98504170e-3;
    temp[1] = 24.98583698e-3;
    temp[2] = 25.98259297e-3;
    Weights.push_back(temp);          //Mg 
    temp.resize(1);
    temp[0] = 26.98153853e-3;
    Weights.push_back(temp);          //Al
    temp.resize(3);
    temp[0] = 27.97692653e-3;
    temp[1] = 28.97649466e-3;
    temp[2] = 29.97377014e-3;
    Weights.push_back(temp);          //Si 
    temp.resize(1);
    temp[0] = 30.97376200e-3;
    Weights.push_back(temp);          //P  
    temp.resize(4);
    temp[0] = 31.97207117e-3;
    temp[1] = 32.97145891e-3;
    temp[2] = 33.96786700e-3;
    temp[3] = 35.96708071e-3;
    Weights.push_back(temp);          //S 
    temp.resize(2);
    temp[0] = 34.96885268e-3;
    temp[1] = 36.96590260e-3;
    Weights.push_back(temp);          //Cl
    temp.resize(3);
    temp[0] = 35.96754511e-3;
    temp[1] = 37.96273211e-3;
    temp[2] = 39.96238312e-3;
    Weights.push_back(temp);          //Ar
    temp[0] = 38.96370649e-3;
    temp[1] = 39.96399817e-3;
    temp[2] = 40.96182526e-3;
    Weights.push_back(temp);          //K 
    temp.resize(6);
    temp[0] = 39.96259086e-3;
    temp[1] = 41.95861783e-3;
    temp[2] = 42.95876644e-3;
    temp[3] = 43.95548156e-3;
    temp[4] = 45.95368902e-3;
    temp[5] = 47.95252276e-3;
    Weights.push_back(temp);          //Ca
    temp.resize(1);
    temp[0] = 44.95590829e-3;
    Weights.push_back(temp);          //Sc
    temp.resize(5);
    temp[0] = 45.95262772e-3;
    temp[1] = 46.95175879e-3;
    temp[2] = 47.94794198e-3;
    temp[3] = 48.94786568e-3;
    temp[4] = 49.94478689e-3;
    Weights.push_back(temp);          //Ti
    temp.resize(2);
    temp[0] = 49.94715602e-3;
    temp[1] = 50.94395705e-3;
    Weights.push_back(temp);          //V 
    temp.resize(4);
    temp[0] = 49.94604184e-3;
    temp[1] = 51.94050624e-3;
    temp[2] = 52.94064816e-3;
    temp[3] = 53.93887917e-3;
    Weights.push_back(temp);          //Cr
    temp.resize(1);
    temp[0] = 54.93804391e-3;
    Weights.push_back(temp);          //Mn
    temp.resize(4);
    temp[0] = 53.93960900e-3;
    temp[1] = 55.93493633e-3;
    temp[2] = 56.93539284e-3;
    temp[3] = 57.93327444e-3;
    Weights.push_back(temp);          //Fe
    temp.resize(1);
    temp[0] = 58.93319430e-3;
    Weights.push_back(temp);          //Co
    temp.resize(5);
    temp[0] = 57.93534242e-3;
    temp[1] = 59.93078589e-3;
    temp[2] = 60.93105558e-3;
    temp[3] = 61.92834538e-3;
    temp[4] = 63.92796683e-3;
    Weights.push_back(temp);          //Ni
    emp.resize(2);
    temp[0] = 62.92959773e-3;
    temp[1] = 64.92778971e-3;
    Weights.push_back(temp);          //Cu
    emp.resize(5);
    temp[0] = 63.92914202e-3;
    temp[1] = 65.92603382e-3;
    temp[2] = 66.92712776e-3;
    temp[3] = 67.92484456e-3;
    temp[4] = 69.92531922e-3;
    Weights.push_back(temp);          //Zn
    temp.resize(2);
    temp[0] = 68.92557351e-3;
    temp[1] = 70.92470259e-3;
    Weights.push_back(temp);          //Ga
    temp.resize(5);
    temp[0] = 69.92424876e-3;
    temp[1] = 71.92207583e-3;
    temp[2] = 72.92345896e-3;
    temp[3] = 73.92117776e-3;
    temp[4] = 75.92140273e-3;
    Weights.push_back(temp);          //Ge
    temp.resize(1);
    temp[0] = 74.92159458e-3;
    Weights.push_back(temp);          //As
    temp.resize(6);
    temp[0] = 73.92247593e-3;
    temp[1] = 75.91921370e-3;
    temp[2] = 76.91991415e-3;
    temp[3] = 77.91730928e-3;
    temp[4] = 79.91652181e-3;
    temp[5] = 81.91669952e-3;
    Weights.push_back(temp);          //Se
    temp.resize(2);
    temp[0] = 78.91833761e-3;
    temp[1] = 80.91628971e-3;
    Weights.push_back(temp);          //Br
    temp.resize(6);
    temp[0] = 77.92036495e-3;
    temp[1] = 79.91637809e-3;
    temp[2] = 81.91348274e-3;
    temp[3] = 82.91412716e-3;
    temp[4] = 83.91149773e-3;
    temp[5] = 85.91061063e-3;
    Weights.push_back(temp);          //Kr
    temp.resize(2);
    temp[0] = 84.91178974e-3;
    temp[1] = 86.90918053e-3;
    Weights.push_back(temp);          //Rb
    temp.resize(4);
    temp[0] = 83.91341911e-3;
    temp[1] = 85.90926061e-3;
    temp[2] = 86.90887751e-3;
    temp[3] = 87.90561251e-3;
    Weights.push_back(temp);          //Sr
    temp.resize(1);
    temp[0] = 88.90584032e-3;
    Weights.push_back(temp);          //Y 
    temp.resize(5);
    temp[0] = 89.90469772e-3;
    temp[1] = 90.90563962e-3;
    temp[2] = 91.90503472e-3;
    temp[3] = 93.90631082e-3;
    temp[4] = 95.90827142e-3;
    Weights.push_back(temp);          //Zr
    temp.resize(1);
    temp[0] = 92.90637302e-3;
    Weights.push_back(temp);          //Nb
    temp.resize(7);
    temp[0] = 91.90680797e-3;
    temp[1] = 93.90508490e-3;
    temp[2] = 94.90583877e-3;
    temp[3] = 95.90467612e-3;
    temp[4] = 96.90601812e-3;
    temp[5] = 97.90540482e-3;
    temp[6] = 99.90747181e-3;
    Weights.push_back(temp);          //Mo
    temp.resize(3);
    temp[0] = 96.90636674e-3;
    temp[1] = 97.90721244e-3;
    temp[2] = 98.90625081e-3;
    Weights.push_back(temp);          //Tc
    temp.resize(7);
    temp[0] = 95.90759025e-3;
    temp[1] = 97.90528687e-3;
    temp[2] = 98.90593411e-3;
    temp[3] = 99.90421431e-3;
    temp[4] = 100.9057691e-3;
    temp[5] = 101.9043441e-3;
    temp[6] = 103.9054275e-3;
    Weights.push_back(temp);          //Ru
    temp.resize(1);
    temp[0] = 102.9054980e-3;
    Weights.push_back(temp);          //Rh
    temp.resize(6);
    temp[0] = 101.9056022e-3;
    temp[1] = 103.9040305e-3;
    temp[2] = 104.9050796e-3;
    temp[3] = 105.9034804e-3;
    temp[4] = 107.9038916e-3;
    temp[5] = 109.9051722e-3;
    Weights.push_back(temp);          //Pd
    temp.resize(2);
    temp[0] = 106.9050916e-3;
    temp[1] = 108.9047553e-3;
    Weights.push_back(temp);          //Ag
    temp.resize(8);
    temp[0] = 105.9064599e-3;
    temp[1] = 107.9041834e-3;
    temp[2] = 109.9030066e-3;
    temp[3] = 110.9041829e-3;
    temp[4] = 111.9027629e-3;
    temp[5] = 112.9044081e-3;
    temp[6] = 113.9033651e-3;
    temp[7] = 115.9047632e-3;
    Weights.push_back(temp);          //Cd
    temp.resize(2);
    temp[0] = 112.9040618e-3;
    temp[1] = 114.9038788e-3;
    Weights.push_back(temp);          //In
    temp.resize(10);
    temp[0] = 111.9048239e-3;
    temp[1] = 113.9027827e-3;
    temp[2] = 114.9033447e-3;
    temp[3] = 115.9017428e-3;
    temp[4] = 116.9029540e-3;
    temp[5] = 117.9016066e-3;
    temp[6] = 118.9033112e-3;
    temp[7] = 119.9022016e-3;
    temp[8] = 121.9034438e-3;
    temp[9] = 123.9052766e-3;
    Weights.push_back(temp);          //Sn
    temp.resize(2);
    temp[0] = 120.9038120e-3;
    temp[1] = 122.9042132e-3;
    Weights.push_back(temp);          //Sb
    temp.resize(8);
    temp[0] = 119.9040593e-3;
    temp[1] = 121.9030435e-3;
    temp[2] = 122.9042698e-3;
    temp[3] = 123.9028171e-3;
    temp[4] = 124.9044299e-3;
    temp[5] = 125.9033109e-3;
    temp[6] = 127.9044613e-3;
    temp[7] = 129.9062227e-3;
    Weights.push_back(temp);          //Te
    temp.resize(1);
    temp[0] = 126.9044719e-3;
    Weights.push_back(temp);          //I 
    temp.resize(9);
    temp[0] = 123.9058920e-3;
    temp[1] = 125.9042983e-3;
    temp[2] = 127.9035310e-3;
    temp[3] = 128.9047809e-3;
    temp[4] = 129.9035093e-3;
    temp[5] = 130.9050841e-3;
    temp[6] = 131.9041551e-3;
    temp[7] = 133.9053947e-3;
    temp[8] = 135.9072145e-3;
    Weights.push_back(temp);          //Xe
    temp.resize(1);
    temp[0] = 132.9054520e-3;
    Weights.push_back(temp);          //Cs
    temp.resize(7);
    temp[0] = 129.9063207e-3;
    temp[1] = 131.9050611e-3;
    temp[2] = 133.9045082e-3;
    temp[3] = 134.9056884e-3;
    temp[4] = 135.9045757e-3;
    temp[5] = 136.9058271e-3;
    temp[6] = 137.9052470e-3;
    Weights.push_back(temp);          //Ba
    temp.resize(2);
    temp[0] = 137.9071149e-3;
    temp[1] = 138.9063563e-3;
    Weights.push_back(temp);          //La
    temp.resize(4);
    temp[0] = 135.9071292e-3;
    temp[1] = 137.9059911e-3;
    temp[2] = 139.9054431e-3;
    temp[3] = 141.9092504e-3;
    Weights.push_back(temp);          //Ce
    temp.resize(1);
    temp[0] = 140.9076576e-3;
    Weights.push_back(temp);          //Pr
    temp.resize(7);
    temp[0] = 141.9077290e-3;
    temp[1] = 142.9098200e-3;
    temp[2] = 143.9100930e-3;
    temp[3] = 144.9125793e-3;
    temp[4] = 145.9131226e-3;
    temp[5] = 147.9168993e-3;
    temp[6] = 149.9209022e-3;
    Weights.push_back(temp);          //Nd
    temp.resize(2);
    temp[0] = 144.9127559e-3;
    temp[1] = 146.9151450e-3;
    Weights.push_back(temp);          //Pm
    temp.resize(7);
    temp[0] = 143.9120065e-3;
    temp[1] = 146.9149044e-3;
    temp[2] = 147.9148292e-3;
    temp[3] = 148.9171921e-3;
    temp[4] = 149.9172829e-3;
    temp[5] = 151.9197397e-3;
    temp[6] = 153.9222169e-3;
    Weights.push_back(temp);          //Sm
    temp.resize(2);
    temp[0] = 150.9198578e-3;
    temp[1] = 152.9212380e-3;
    Weights.push_back(temp);          //Eu
    temp.resize(7);
    temp[0] = 151.9197995e-3;
    temp[1] = 153.9208741e-3;
    temp[2] = 154.9226305e-3;
    temp[3] = 155.9221312e-3;
    temp[4] = 156.9239686e-3;
    temp[5] = 157.9241123e-3;
    temp[6] = 159.9270624e-3;
    Weights.push_back(temp);          //Gd   
    temp.resize(1);
    temp[0] = 158.9253547e-3;
    Weights.push_back(temp);          //Tb
    temp.resize(7);
    temp[0] = 155.9242847e-3;
    temp[1] = 157.9244159e-3;
    temp[2] = 159.9252046e-3;
    temp[3] = 160.9269405e-3;
    temp[4] = 161.9268056e-3;
    temp[5] = 162.9287383e-3;
    temp[6] = 163.9291819e-3;
    Weights.push_back(temp);          //Dy
    temp.resize(1);
    temp[0] = 164.9303288e-3;
    Weights.push_back(temp);          //Ho
    temp.resize(6);
    temp[0] = 161.9287884e-3;
    temp[1] = 163.9292088e-3;
    temp[2] = 165.9302995e-3;
    temp[3] = 166.9320546e-3;
    temp[4] = 167.9323767e-3;
    temp[5] = 169.9354702e-3;
    Weights.push_back(temp);          //Er
    temp.resize(1);
    temp[0] = 168.9342179e-3;
    Weights.push_back(temp);          //Tm
    temp.resize(7);
    temp[0] = 167.9338896e-3;
    temp[1] = 169.9347664e-3;
    temp[2] = 170.9363302e-3;
    temp[3] = 171.9363859e-3;
    temp[4] = 172.9382151e-3;
    temp[5] = 173.9388664e-3;
    temp[6] = 175.9425764e-3;
    Weights.push_back(temp);          //Yb
    temp.resize(2);
    temp[0] = 174.9407752e-3;
    temp[1] = 175.9426897e-3;
    Weights.push_back(temp);          //Lu
    temp.resize(6);
    temp[0] = 173.9400461e-3;
    temp[1] = 175.9414076e-3;
    temp[2] = 176.9432277e-3;
    temp[3] = 177.9437058e-3;
    temp[4] = 178.9458232e-3;
    temp[5] = 179.9465570e-3;
    Weights.push_back(temp);          //Hf
    temp.resize(2);
    temp[0] = 179.9474648e-3;
    temp[1] = 180.9479958e-3;
    Weights.push_back(temp);          //Ta
    temp.resize(5);
    temp[0] = 179.9467108e-3;
    temp[1] = 181.9482039e-3;
    temp[2] = 182.9502228e-3;
    temp[3] = 183.9509309e-3;
    temp[4] = 185.9543628e-3;
    Weights.push_back(temp);          //W 
    temp.resize(2);
    temp[0] = 184.9529545e-3;
    temp[1] = 186.9557501e-3;
    Weights.push_back(temp);          //Re
    temp.resize(7);
    temp[0] = 183.9524885e-3;
    temp[1] = 185.9538350e-3;
    temp[2] = 186.9557474e-3;
    temp[3] = 187.9558352e-3;
    temp[4] = 188.9581442e-3;
    temp[5] = 189.9584437e-3;
    temp[6] = 191.9614770e-3;
    Weights.push_back(temp);          //Os
    temp.resize(2);
    temp[0] = 190.9605893e-3;
    temp[1] = 192.9629216e-3;
    Weights.push_back(temp);          //Ir
    temp.resize(6);
    temp[0] = 189.9599298e-3;
    temp[1] = 191.9610387e-3;
    temp[2] = 193.9626809e-3;
    temp[3] = 194.9647917e-3;
    temp[4] = 195.9649521e-3;
    temp[5] = 197.9678949e-3;
    Weights.push_back(temp);          //Pt
    temp.resize(1);
    temp[0] = 196.9665688e-3;
    Weights.push_back(temp);          //Au
    temp.resize(7);
    temp[0] = 195.9658326e-3;
    temp[1] = 197.9667686e-3;
    temp[2] = 198.9682806e-3;
    temp[3] = 199.9683266e-3;
    temp[4] = 200.9703028e-3;
    temp[5] = 201.9706434e-3;
    temp[6] = 203.9734940e-3;
    Weights.push_back(temp);          //Hg
    temp.resize(2);
    temp[0] = 202.9723446e-3;
    temp[1] = 204.9744278e-3;
    Weights.push_back(temp);          //Tl
    temp.resize(4);
    temp[0] = 203.9730440e-3;
    temp[1] = 205.9744657e-3;
    temp[2] = 206.9758973e-3;
    temp[3] = 207.9766525e-3;
    Weights.push_back(temp);          //Pb
    temp.resize(1);
    temp[0] = 208.9803991e-3;
    Weights.push_back(temp);          //Bi
    temp.resize(2);
    temp[0] = 208.9824308e-3;
    temp[1] = 209.9828741e-3;
    Weights.push_back(temp);          //Po
    temp[0] = 209.9871480e-3;
    temp[1] = 210.9874966e-3;
    Weights.push_back(temp);          //At
    temp.resize(3);
    temp[0] = 210.9906012e-3;
    temp[1] = 220.0113941e-3;
    temp[2] = 222.0175782e-3;
    Weights.push_back(temp);          //Rn
    temp.resize(1);
    temp[0] = 223.0197360e-3;
    Weights.push_back(temp);          //Fr
    temp.resize(4);
    temp[0] = 223.0185023e-3;
    temp[1] = 224.0202120e-3;
    temp[2] = 226.0254103e-3;
    temp[3] = 228.0310707e-3;
    Weights.push_back(temp);          //Ra
    temp.resize(1);
    temp[0] = 227.0277523e-3;
    Weights.push_back(temp);          //Ac
    temp.resize(2);
    temp[0] = 230.0331341e-3;
    temp[1] = 232.0380558e-3;
    Weights.push_back(temp);          //Th
    temp.resize(1);
    temp[0] = 231.0358842e-3;
    Weights.push_back(temp);          //Pa
    temp.resize(5);
    temp[0] = 233.0396355e-3;
    temp[1] = 234.0409523e-3;
    temp[2] = 235.0439301e-3;
    temp[3] = 236.0455682e-3;
    temp[4] = 238.0507884e-3;
    Weights.push_back(temp);          //U 
    temp.resize(2);
    temp[0] = 236.0465705e-3;
    temp[1] = 237.0481736e-3;
    Weights.push_back(temp);          //Np
    temp.resize(6);
    temp[0] = 238.0495601e-3;
    temp[1] = 239.0521636e-3;
    temp[2] = 240.0538138e-3;
    temp[3] = 241.0568517e-3;
    temp[4] = 242.0587428e-3;
    temp[5] = 244.0642054e-3;
    Weights.push_back(temp);          //Pu
    temp.resize(2);
    temp[0] = 241.0568293e-3;
    temp[1] = 243.0613813e-3;
    Weights.push_back(temp);          //Am
    temp.resize(6);
    temp[0] = 243.0613893e-3;
    temp[1] = 244.0627528e-3;
    temp[2] = 245.0654915e-3;
    temp[3] = 246.0672238e-3;
    temp[4] = 247.0703541e-3;
    temp[5] = 248.0723500e-3;
    Weights.push_back(temp);          //Cm
    temp.resize(2);
    temp[0] = 247.0703074e-3;
    temp[1] = 249.0749877e-3;
    Weights.push_back(temp);          //Bk
    temp.resize(4);
    temp[0] = 249.0748539e-3;
    temp[1] = 250.0764062e-3;
    temp[2] = 251.0795886e-3;
    temp[3] = 252.0816273e-3;
    Weights.push_back(temp);          //Cf
    temp.resize(3);
    temp[0] = 252.0829805e-3;
    Weights.push_back(temp);          //Es
    temp[0] = 257.0951062e-3;
    Weights.push_back(temp);          //Fm
    temp.resize(2);
    temp[0] = 258.0984316e-3;
    temp[1] = 260.1036534e-3;
    Weights.push_back(temp);          //Md
    temp.resize(3);
    temp[0] = 19.99244018e-3;
    Weights.push_back(temp);          //No
    temp.resize(1);
    temp[0] = 259.1010311e-3
    Weights.push_back(temp);          //Lr
    temp[0] = 262.1096122e-3;
    Weights.push_back(temp);          //Rf
    temp[0] = 267.1217962e-3;
    Weights.push_back(temp);          //Db
    temp[0] = 268.1256757e-3;
    Weights.push_back(temp);          //Sg
    temp[0] = 271.1339363e-3;
    Weights.push_back(temp);          //Bh
    temp[0] = 272.1382658e-3;
    Weights.push_back(temp);          //Hs
    temp[0] = 270.1342927e-3;
    Weights.push_back(temp);          //Mt
    temp[0] = 276.1515959e-3;
    Weights.push_back(temp);          //Ds
    temp[0] = 281.1645159e-3;
    Weights.push_back(temp);          //Rg
    temp[0] = 280.1651461e-3;
    Weights.push_back(temp);          //Cn
    temp[0] = 285.1771260e-3;
    Weights.push_back(temp);          //Nh
    temp[0] = 284.1787362e-3;
    Weights.push_back(temp);          //Fl
    temp[0] = 289.1904260e-3;
    Weights.push_back(temp);          //Mc
    temp[0] = 288.1927462e-3;
    Weights.push_back(temp);          //Lv
    temp[0] = 293.2044960e-3;
    Weights.push_back(temp);          //Ts
    temp[0] = 292.2074675e-3;
    Weights.push_back(temp);          //Og
    temp[0] = 294.2139271e-3;
  }
  double Weight(size_t atomicnumber) {return Weights[atomicnumber-1];}
};

#endif //_Isotope_Weight_Package_
