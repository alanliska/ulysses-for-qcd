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

#ifndef _Molecule_Analyser_
#define _Molecule_Analyser_
#include <iostream>
#include <sstream>
#include <string>
#include "Utilities.hpp"
#include "../math/MatrixPackage.hpp"
#include "../atoms/AtomPackage.hpp"

//description:
//functions that check the molecular structure for consistency

void PrintGroupStat(const std::vector<size_t> & funcgroup, std::string group) {
  //function that prints functional group information
  size_t vectorsz = funcgroup.size();
  std::cout << group << ": " << vectorsz << "\n";
  std::cout << "           ";
  for (size_t idfg = 0; idfg < vectorsz; ++idfg) {
    std::cout << funcgroup[idfg] << " ";
  }
  std::cout << std::endl;
}
void GetGroupEvaluation(const std::vector<std::string> & funcgroups) {
  //function that prints functional group information
  size_t Natoms = funcgroups.size();
  //count groups
  std::vector<size_t> helium;
  std::vector<size_t> lithium;
  std::vector<size_t> liplus;
  std::vector<size_t> berylium;
  std::vector<size_t> boron;
  std::vector<size_t> bminus;
  std::vector<size_t> csp3;
  std::vector<size_t> csp2;
  std::vector<size_t> alkynes;
  std::vector<size_t> nitriles;
  std::vector<size_t> cyano;
  std::vector<size_t> aldehydes;
  std::vector<size_t> ketones;
  std::vector<size_t> imines;
  std::vector<size_t> enamines;
  std::vector<size_t> enols;
  std::vector<size_t> carboxyls;
  std::vector<size_t> carboxylates;
  std::vector<size_t> amides;
  std::vector<size_t> hydroxyamides;
  std::vector<size_t> alkoxides;
  std::vector<size_t> waters;
  std::vector<size_t> hydroxides;
  std::vector<size_t> carbonates;
  std::vector<size_t> carbamates;
  std::vector<size_t> alcohols;
  std::vector<size_t> ethers;
  std::vector<size_t> peroxides;
  std::vector<size_t> amines;
  std::vector<size_t> nitroso;
  std::vector<size_t> nitro;
  std::vector<size_t> nitrates;
  std::vector<size_t> guanidiniums;
  std::vector<size_t> nh3;
  std::vector<size_t> nplus;
  std::vector<size_t> oplus;
  std::vector<size_t> haloalkanes;
  std::vector<size_t> neon;
  std::vector<size_t> sodium;
  std::vector<size_t> naplus;
  std::vector<size_t> magnesium;
  std::vector<size_t> mgplus;
  std::vector<size_t> aluminium;
  std::vector<size_t> silicon;
  std::vector<size_t> phosphates;
  std::vector<size_t> phosphonates;
  std::vector<size_t> phosphinates;
  std::vector<size_t> phosphoryls;
  std::vector<size_t> aminophosphates;
  std::vector<size_t> aminophosphinates;
  std::vector<size_t> aminophosphoryls;
  std::vector<size_t> thiophosphates;
  std::vector<size_t> aminothiophosphates;
  std::vector<size_t> sulfonates;
  std::vector<size_t> sulfonateesthers;
  std::vector<size_t> h2s;
  std::vector<size_t> thiols;
  std::vector<size_t> splus;
  std::vector<size_t> sulfides;
  std::vector<size_t> disulfides;
  std::vector<size_t> sulfoxides;
  std::vector<size_t> argon;
  std::vector<size_t> potassium;
  std::vector<size_t> kplus;
  std::vector<size_t> calcium;
  std::vector<size_t> caplus;
  std::vector<size_t> scadium;
  std::vector<size_t> titanium;
  std::vector<size_t> vanadium;
  std::vector<size_t> chromium;
  std::vector<size_t> manganese;
  std::vector<size_t> iron;
  std::vector<size_t> cobalt;
  std::vector<size_t> nickel;
  std::vector<size_t> copper;
  std::vector<size_t> zink;
  std::vector<size_t> gallium;
  std::vector<size_t> germanium;
  std::vector<size_t> arsenic;
  std::vector<size_t> selenium;
  std::vector<size_t> kripton;
  std::vector<size_t> rubidium;
  std::vector<size_t> rbplus;
  std::vector<size_t> strontium;
  std::vector<size_t> srplus;
  std::vector<size_t> ruthenium;
  std::vector<size_t> rhodium;
  std::vector<size_t> palladium;
  std::vector<size_t> silver;
  std::vector<size_t> cadmium;
  std::vector<size_t> tin;
  std::vector<size_t> antimony;
  std::vector<size_t> tellurium;
  std::vector<size_t> osmium;
  std::vector<size_t> iridium;
  std::vector<size_t> platinum;
  std::vector<size_t> gold;
  std::vector<size_t> freddie;
  std::vector<size_t> tallium;
  std::vector<size_t> lead;
  std::vector<size_t> bismuth;
  for (size_t idfg = 0; idfg < Natoms; ++idfg) {
    if (funcgroups[idfg] == "He") {helium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Li") {lithium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Li+") {liplus.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Be") {berylium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "B") {boron.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "B-") {bminus.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-yne") {alkynes.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-nitrile") {nitriles.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-cyano") {cyano.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-sp2") {csp2.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-aldehyde") {aldehydes.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-ketone") {ketones.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "N-imine") {imines.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "N-enamine") {enamines.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "O-enol") {enols.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-carboxyl") {carboxyls.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-carboxylate") {carboxylates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-amide") {amides.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "N-hydroxyamide") {hydroxyamides.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-carbonate") {carbonates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "N-carbamate") {carbamates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-guanidinium") {guanidiniums.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-sp3") {csp3.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-alkoxyde") {alkoxides.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "C-alcohol") {alcohols.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "O-ether") {ethers.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "O-peroxide") {peroxides.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "N-amine") {amines.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "N-nitroso") {nitroso.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "N-nitro") {nitro.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "N-nitrate") {nitrates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "NH3") {nh3.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "N+") {nplus.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "OH2") {waters.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "hydroxide") {hydroxides.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "O+") {oplus.push_back(idfg + 1);}
    else if ((funcgroups[idfg] == "F-haloalkane")||(funcgroups[idfg] == "Cl-haloalkane")||(funcgroups[idfg] == "Br-haloalkane")||(funcgroups[idfg] == "I-haloalkane")) {haloalkanes.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ne") {neon.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Na") {sodium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Na+") {naplus.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Mg") {magnesium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Mg++") {mgplus.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Al") {aluminium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Si") {silicon.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "P-phosphate") {phosphates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "P-phosphonate") {phosphonates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "P-phosphinate") {phosphinates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "P-phosphoryl") {phosphoryls.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "P-aminophosphate") {aminophosphates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "P-aminophosphinate") {aminophosphinates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "P-aminophosphoryl") {aminophosphoryls.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "P-thiophosphate") {thiophosphates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "P-aminothiophosphate") {aminothiophosphates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "S-sulfonate") {sulfonates.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "S-sulfate") {sulfonateesthers.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "SH2") {h2s.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "S-thiol") {thiols.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "S+") {splus.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "sulfide") {sulfides.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "S-disulfide") {disulfides.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "S-sulfoxide") {sulfoxides.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ar") {argon.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "K") {potassium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "K+") {kplus.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ca") {calcium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ca++") {caplus.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Sc") {scadium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ti") {titanium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "V") {vanadium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Cr") {chromium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Mn") {manganese.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Fe") {iron.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Co") {cobalt.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ni") {nickel.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Cu") {copper.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Zn") {zink.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ga") {gallium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ge") {germanium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "As") {arsenic.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Se") {selenium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Kr") {kripton.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Rb") {rubidium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Rb+") {rbplus.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Sr") {strontium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Sr++") {srplus.push_back(idfg + 1);}
    //"Y"     still to implement
    //"Zr"     still to implement
    //"Nb"     still to implement
    //"Mo"     still to implement
    //"Tc"     still to implement
    else if (funcgroups[idfg] == "Ru") {ruthenium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Rh") {rhodium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Pd") {palladium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ag") {silver.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Cd") {cadmium.push_back(idfg + 1);}
    //"In"     still to implement
    else if (funcgroups[idfg] == "Sn") {tin.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Sb") {antimony.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Te") {tellurium.push_back(idfg + 1);}
    //"Xe"     still to implement
    //"Cs"     still to implement
    //"Cs+"     still to implement
    //"Ba"     still to implement
    //"Ba+"     still to implement
    //"La"     still to implement
    //"Ce"     still to implement
    //"Pr"     still to implement
    //"Nd"     still to implement
    //"Pm"     still to implement
    //"Sm"     still to implement
    //"Eu"     still to implement
    //"Gd"     still to implement
    //"Tb"     still to implement
    //"Dy"     still to implement
    //"Ho"     still to implement
    //"Er"     still to implement
    //"Tm"     still to implement
    //"Yb"     still to implement
    //"Lu"     still to implement
    //"Hf"     still to implement
    //"Ta"     still to implement
    //"W"     still to implement
    //"Re"     still to implement
    else if (funcgroups[idfg] == "Os") {osmium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Ir") {iridium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Pt") {platinum.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Au") {gold.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Hg") {freddie.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Tl") {tallium.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Pb") {lead.push_back(idfg + 1);}
    else if (funcgroups[idfg] == "Bi") {bismuth.push_back(idfg + 1);}
    //"Po"     still to implement
    //"At"     still to implement
    //"Rn"     still to implement
    //"Fr"     still to implement
    //"Ra"     still to implement
    //"Ac"     still to implement
  }
  //print stats
  if (helium.size() > 0) {PrintGroupStat(helium,"He");}
  if (lithium.size() > 0) {PrintGroupStat(lithium,"Li");}
  if (liplus.size() > 0) {PrintGroupStat(liplus,"Li+");}
  if (berylium.size() > 0) {PrintGroupStat(berylium,"Be");}
  if (boron.size() > 0) {PrintGroupStat(boron,"B");}
  if (bminus.size() > 0) {PrintGroupStat(bminus,"B-");}
  if (alkynes.size() > 0) {PrintGroupStat(alkynes,"alkynes");}
  if (nitriles.size() > 0) {PrintGroupStat(nitriles,"nitriles");}
  if (cyano.size() > 0) {PrintGroupStat(cyano,"cyano");}
  if (csp2.size() > 0) {PrintGroupStat(csp2,"sp2 carbons");}
  if (ketones.size() > 0) {PrintGroupStat(ketones,"ketones");}
  if (imines.size() > 0) {PrintGroupStat(imines,"imines");}
  if (enamines.size() > 0) {PrintGroupStat(enamines,"enamines");}
  if (aldehydes.size() > 0) {PrintGroupStat(aldehydes,"aldehydes");}
  if (enols.size() > 0) {PrintGroupStat(enols,"enols");}
  if (carboxyls.size() > 0) {PrintGroupStat(carboxyls,"carboxyls");}
  if (carboxylates.size() > 0) {PrintGroupStat(carboxylates,"carboxylates");}
  if (amides.size() > 0) {PrintGroupStat(amides,"amides");}
  if (hydroxyamides.size() > 0) {PrintGroupStat(hydroxyamides,"hydroxyamides");}
  if (carbonates.size() > 0) {PrintGroupStat(carbonates,"carbonates");}
  if (carbamates.size() > 0) {PrintGroupStat(carbamates,"carbamates");}
  if (guanidiniums.size() > 0) {PrintGroupStat(guanidiniums,"guanidiniums");}
  if (csp3.size() > 0) {PrintGroupStat(csp3,"sp3 carbons");}
  if (alkoxides.size() > 0) {PrintGroupStat(alkoxides,"alkoxides");}
  if (alcohols.size() > 0) {PrintGroupStat(alcohols,"alcohols");}
  if (ethers.size() > 0) {PrintGroupStat(ethers,"ethers");}
  if (peroxides.size() > 0) {PrintGroupStat(peroxides,"peroxides");}
  if (amines.size() > 0) {PrintGroupStat(amines,"amines");}
  if (nitroso.size() > 0) {PrintGroupStat(nitroso,"nitroso");}
  if (nitro.size() > 0) {PrintGroupStat(nitro,"nitro");}
  if (nitrates.size() > 0) {PrintGroupStat(nitrates,"nitrates");}
  if (nh3.size() > 0) {PrintGroupStat(nh3,"nh3");}
  if (nplus.size() > 0) {PrintGroupStat(nplus,"n+");}
  if (waters.size() > 0) {PrintGroupStat(waters,"waters");}
  if (hydroxides.size() > 0) {PrintGroupStat(hydroxides,"hydroxides");}
  if (oplus.size() > 0) {PrintGroupStat(oplus,"o+");}
  if (haloalkanes.size() > 0) {PrintGroupStat(haloalkanes,"haloalkanes");}
  if (neon.size() > 0) {PrintGroupStat(neon,"Ne");}
  if (sodium.size() > 0) {PrintGroupStat(sodium,"Na");}
  if (naplus.size() > 0) {PrintGroupStat(naplus,"Na+");}
  if (magnesium.size() > 0) {PrintGroupStat(magnesium,"Mg");}
  if (mgplus.size() > 0) {PrintGroupStat(mgplus,"M++");}
  if (aluminium.size() > 0) {PrintGroupStat(aluminium,"Al");}
  if (silicon.size() > 0) {PrintGroupStat(silicon,"Si");}
  if (phosphates.size() > 0) {PrintGroupStat(phosphates,"phosphates");}
  if (phosphonates.size() > 0) {PrintGroupStat(phosphonates,"phosphonates");}
  if (phosphinates.size() > 0) {PrintGroupStat(phosphinates,"phosphinates");}
  if (phosphoryls.size() > 0) {PrintGroupStat(phosphoryls,"phosphoryls");}
  if (aminophosphates.size() > 0) {PrintGroupStat(aminophosphates,"aminophosphates");}
  if (aminophosphinates.size() > 0) {PrintGroupStat(aminophosphinates,"aminophosphinates");}
  if (aminophosphoryls.size() > 0) {PrintGroupStat(aminophosphoryls,"aminophosphoryls");}
  if (thiophosphates.size() > 0) {PrintGroupStat(thiophosphates,"thiophosphates");}
  if (aminothiophosphates.size() > 0) {PrintGroupStat(aminothiophosphates,"aminothiophosphates");}
  if (h2s.size() > 0) {PrintGroupStat(h2s,"hydrogen sulfides");}
  if (thiols.size() > 0) {PrintGroupStat(thiols,"thiols");}
  if (splus.size() > 0) {PrintGroupStat(splus,"S+");}
  if (sulfides.size() > 0) {PrintGroupStat(sulfides,"sulfides");}
  if (disulfides.size() > 0) {PrintGroupStat(disulfides,"disulfides");}
  if (sulfoxides.size() > 0) {PrintGroupStat(sulfoxides,"sulfoxides");}
  if (sulfonates.size() > 0) {PrintGroupStat(sulfonates,"sulfonates");}
  if (sulfonateesthers.size() > 0) {PrintGroupStat(sulfonateesthers,"sulfate");}
  if (argon.size() > 0) {PrintGroupStat(argon,"Ar");}
  if (potassium.size() > 0) {PrintGroupStat(potassium,"K");}
  if (kplus.size() > 0) {PrintGroupStat(kplus,"K+");}
  if (calcium.size() > 0) {PrintGroupStat(calcium,"Ca");}
  if (caplus.size() > 0) {PrintGroupStat(caplus,"Ca++");}
  if (scadium.size() > 0) {PrintGroupStat(scadium,"Sc");}
  if (titanium.size() > 0) {PrintGroupStat(titanium,"Ti");}
  if (vanadium.size() > 0) {PrintGroupStat(vanadium,"V");}
  if (chromium.size() > 0) {PrintGroupStat(chromium,"Cr");}
  if (manganese.size() > 0) {PrintGroupStat(manganese,"Mn");}
  if (iron.size() > 0) {PrintGroupStat(iron,"Fe");}
  if (cobalt.size() > 0) {PrintGroupStat(cobalt,"Co");}
  if (nickel.size() > 0) {PrintGroupStat(nickel,"Ni");}
  if (copper.size() > 0) {PrintGroupStat(copper,"Cu");}
  if (zink.size() > 0) {PrintGroupStat(zink,"Zn");}
  if (gallium.size() > 0) {PrintGroupStat(gallium,"Ga");}
  if (germanium.size() > 0) {PrintGroupStat(germanium,"Ge");}
  if (arsenic.size() > 0) {PrintGroupStat(arsenic,"As");}
  if (selenium.size() > 0) {PrintGroupStat(selenium,"Se");}
  if (kripton.size() > 0) {PrintGroupStat(kripton,"Kr");}
  if (rubidium.size() > 0) {PrintGroupStat(rubidium,"Rb");}
  if (rbplus.size() > 0) {PrintGroupStat(rbplus,"Rb+");}
  if (strontium.size() > 0) {PrintGroupStat(strontium,"Sr");}
  if (srplus.size() > 0) {PrintGroupStat(srplus,"Sr++");}
  if (ruthenium.size() > 0) {PrintGroupStat(ruthenium,"Ru");}
  if (rhodium.size() > 0) {PrintGroupStat(rhodium,"Rh");}
  if (palladium.size() > 0) {PrintGroupStat(palladium,"Pd");}
  if (silver.size() > 0) {PrintGroupStat(silver,"Ag");}
  if (cadmium.size() > 0) {PrintGroupStat(cadmium,"Cd");}
  if (tin.size() > 0) {PrintGroupStat(tin,"Sn");}
  if (antimony.size() > 0) {PrintGroupStat(antimony,"Sb");}
  if (tellurium.size() > 0) {PrintGroupStat(tellurium,"Te");}
  if (osmium.size() > 0) {PrintGroupStat(osmium,"Os");}
  if (iridium.size() > 0) {PrintGroupStat(iridium,"Ir");}
  if (platinum.size() > 0) {PrintGroupStat(platinum,"Pt");}
  if (gold.size() > 0) {PrintGroupStat(gold,"Au");}
  if (freddie.size() > 0) {PrintGroupStat(freddie,"Hg");}
  if (tallium.size() > 0) {PrintGroupStat(tallium,"Tl");}
  if (lead.size() > 0) {PrintGroupStat(lead,"Pb");}
  if (bismuth.size() > 0) {PrintGroupStat(bismuth,"Bi");}
}
std::vector<std::string> CheckGroupsByNeighborhood(std::vector<size_t> & atoms, matrix<int> & neighbours, matrixE & geometry) {
  //function checking the functional groups for consistency
  //list of groups that might still be needed
  //groups contaning carbon
  //    aromatics
  //groups containing nitrogen
  //    amidine
  //    azo
  //    carbodiimides
  //    diazo
  //    guanidine
  //    heterocycles with N
  //    hydrazone
  //    nitrites
  //    nitrone
  //    oxime
  //    ureas
  //groups containing oxygen
  //    acid halide
  //    esther -> present but not differentiated from carboxyl acids
  //    anhydride -> present but not differentiated from carboxyl acids
  //    heterocycles with oxygen
  size_t maxnneigh = neighbours.cols();
  size_t Natoms = atoms.size();
  size_t nneighbours;
  size_t numberH;
  size_t numberC;
  size_t numberN;
  size_t numberO;
  size_t numberX;
  size_t numberP;
  size_t numberS;
  size_t numberother;
  size_t numberH_2;
  size_t numberC_2;
  size_t numberO_2;
  size_t numberother_2;
  size_t numberHC;
  size_t Oneighbours;
  size_t Nneighbours;
  size_t nneighboursneighbours;
  size_t nneighboursneighbours2;
  size_t nneighboursneighbours3;
  int idbtm;
  int idctm;
  int iddtm;
  int idetm;
  bool protonatedoxygens;
  bool atom1sp3;
  bool atom2sp3;
  bool atom3sp3;
  bool NObond;
  std::vector<std::string> atomtypes(Natoms,"");
  //first run, find groups based on neighbours
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    if (atomtypes[idAtm] + "blahblah" != "blahblah") {continue;}
    nneighbours = 0;
    for (size_t idn = 0; idn < maxnneigh; ++idn) {
      if (neighbours(idAtm + 1,idn + 1) == 0) {break;}
      ++nneighbours;
    }
    //though it makes the code very repetitive and extense, I believe that this approach on elements is the cleanest: all groups related to a single element are together
    if (atoms[idAtm] == 1) {atomtypes[idAtm] = "H";}                    //H
    else if (atoms[idAtm] == 2) {atomtypes[idAtm] = "He";}              //He
    else if (atoms[idAtm] == 3) {                                       //Li
      //assume neutral Li, except if it is bound to one atom more electronegative; more than 1 neighbour is weird
      if (nneighbours == 0) {atomtypes[idAtm] = "Li";}                        //no assumption on charge
      else if (nneighbours == 1) {
        idbtm = atoms[neighbours(idAtm + 1,1) - 1];
        //based on electronegativities
        if ((idbtm == 3)||(idbtm == 11)||(idbtm == 19)||(idbtm == 37)||(idbtm == 38)||(idbtm == 55)||(idbtm == 56)||(idbtm == 87)||(idbtm == 88)) {atomtypes[idAtm] = "Li";}
        else {atomtypes[idAtm] = "Li+";}
      }
      else {std::cout << "unknown functional group on (Li) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 4) {                                       //Be
      //no distinction
      atomtypes[idAtm] = "Be";
    }
    else if (atoms[idAtm] == 5) {                                       //B
      //if 2 neighbours, should be B with double bond to something; 3 neighbours is normal and expected; 4 neighbours should have a negative charge somewhere; for simplicity let B take it
      if (nneighbours == 2) {atomtypes[idAtm] = "B";}                         //there might be 
      else if (nneighbours == 3) {atomtypes[idAtm] = "B";}                    //don't want to detail
      else if (nneighbours == 4) {atomtypes[idAtm] = "B-";}                   //almost surely like this
      else {std::cout << "unknown functional group on (B) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 6) {                                       //C
      //check neighbours
      numberH = 0;
      numberC = 0;
      numberN = 0;
      numberO = 0;
      numberP = 0;
      numberS = 0;
      numberX = 0;
      numberother = 0;
      idctm = -1;
      iddtm = -1;
      idetm = -1;
      Oneighbours = 0;        //how many neighbours do eventual oxygen neighbours have?
      Nneighbours = 0;        //how many neighbours do eventual nitrogen neighbours have?
      NObond = false;
      for (size_t idn = 0; idn < maxnneigh; ++idn) {
        idbtm = neighbours(idAtm + 1,idn + 1) - 1;
        if (idbtm < 0) {break;}
        if (atoms[idbtm] == 1) {++numberH;}
        else if (atoms[idbtm] == 6) {++numberC;}
        else if (atoms[idbtm] == 7) {
          ++numberN;
          if (idctm < 0) {idctm = idbtm;}
          else {
            if (iddtm < 0) {iddtm = idbtm;}
            else {idetm = idbtm;}
          }
          //how many neighbours does amide nitrogen have?
          for (size_t idn2 = 0; idn2 < maxnneigh; ++idn2) {
            if (neighbours(idbtm + 1,idn2 + 1) == 0) {break;}
            if (neighbours(idbtm + 1,idn2 + 1) - 1 != idAtm) {
              ++Nneighbours;
              if (atoms[neighbours(idbtm + 1,idn2 + 1) - 1] == 8) {NObond = true;}
            }
          }
        }
        else if (atoms[idbtm] == 8) {
          ++numberO;
          if (idctm < 0) {idctm = idbtm;}
          else {
            if (iddtm < 0) {iddtm = idbtm;}
            else {idetm = idbtm;}
          }
          //is oxygen bound to something?
          for (size_t idn2 = 0; idn2 < maxnneigh; ++idn2) {
            if (neighbours(idbtm + 1,idn2 + 1) == 0) {break;}
            if (neighbours(idbtm + 1,idn2 + 1) - 1 != idAtm) {++Oneighbours;}
          }
        }
        else if (atoms[idbtm] == 15) {++numberP;}
        else if (atoms[idbtm] == 16) {++numberS;}
        else if ((atoms[idbtm] == 9)||(atoms[idbtm] == 17)||(atoms[idbtm] == 35)||(atoms[idbtm] == 53)) {++numberX;}
        else {++numberother;}
      }
      if (nneighbours == 2) {                                               //sp
        //look for alkyne, RCN and RNC
        idbtm = neighbours(idAtm + 1,1) - 1;
        idctm = neighbours(idAtm + 1,2) - 1;
        if (((atoms[idbtm] == 1)&&(atoms[idctm] == 6))||((atoms[idbtm] == 6)&&(atoms[idctm] == 1))||((atoms[idbtm] == 6)&&(atoms[idctm] == 6))) {
          atomtypes[idAtm] = "C-yne";
          if (atoms[idbtm] == 6) atomtypes[idbtm] = "C-yne";
          if (atoms[idctm] == 6) atomtypes[idctm] = "C-yne";
        }
        else if ((atoms[idbtm] == 7)||(atoms[idctm] == 7)) {
          //put nitrogen on iddtm
          if (atoms[idbtm] == 7) {
            iddtm = idbtm;
            idetm = idctm;
          }
          else {
            iddtm = idctm;
            idetm = idbtm;
          }
          nneighboursneighbours = 0;
          for (size_t idn = 0; idn < maxnneigh; ++idn) {
            if (neighbours(iddtm + 1,idn + 1) == 0) {break;}
            ++nneighboursneighbours;
          }
          if (nneighboursneighbours == 1) {        //check whether N has one neighbour (sp)
            if ((atoms[idetm] == 1)||(atoms[idetm] == 6)) {       //consider nitrile bound to carbon
              atomtypes[idAtm] = "C-nitrile";
              atomtypes[iddtm] = "N-nitrile";
            }
            else {                                                //consider cyano bound to other stuff
              atomtypes[idAtm] = "C-cyano";
              atomtypes[iddtm] = "N-cyano";
            }
          }
          else {std::cout << "unknown functional group on (CN) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
        else {std::cout << "unknown functional group on (Csp) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
      }
      else if (nneighbours == 3) {                                          //sp2
        //look for aldehydes, ketones, imines, enamines, enols, carboxylates, carboxyls, amides, carbamates, carbonates
        if ((numberH == 1)&&(numberC == 1)&&(numberO == 1)) {
          if (Oneighbours == 0) {          //aldehyde
            atomtypes[idAtm] = "C-aldehyde";
            atomtypes[idctm] = "O-aldehyde";
          }
          else if (Oneighbours == 1) {     //enol, assume there is a double bond somewhere else
            idbtm = neighbours(idctm + 1,1) - 1;
            if (idbtm == idAtm) {idbtm = neighbours(idctm + 1,2) - 1;}
            if ((atoms[idbtm] != 15)&&(atoms[idbtm] != 16)) {
              atomtypes[idAtm] = "C-enol";
              atomtypes[idctm] = "O-enol";
            }
          }
          else {std::cout << "unknown functional group on (Csp2-aldehyde) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
        else if ((numberC == 2)&&(numberO == 1)) {
          if (Oneighbours == 0) {          //ketone
            atomtypes[idAtm] = "C-ketone";
            atomtypes[idctm] = "O-ketone";
          }
          else if (Oneighbours == 1) {     //enol, assume there is a double bond somewhere else
            idbtm = neighbours(idctm + 1,1) - 1;
            if (idbtm == idAtm) {idbtm = neighbours(idctm + 1,2) - 1;}
            if ((atoms[idbtm] != 15)&&(atoms[idbtm] != 16)) {
              atomtypes[idAtm] = "C-enol";
              atomtypes[idctm] = "O-enol";
            }
          }
          else {std::cout << "unknown functional group on (Csp2-ketone) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
        else if ((numberC + numberH == 2)&&(numberN == 1)) {
          if (Nneighbours == 1) {          //imine
            atomtypes[idAtm] = "C-imine";
            atomtypes[idctm] = "N-imine";
          }
          else if (Nneighbours == 2) {     //enamine, assume there is a double bond somewhere else
            idbtm = neighbours(idctm + 1,1) - 1;
            if (idbtm == idAtm) {idbtm = neighbours(idctm + 1,2) - 1;}
            if ((atoms[idbtm] != 15)&&(atoms[idbtm] != 16)) {
              atomtypes[idAtm] = "C-enamine";
              atomtypes[idctm] = "N-enamine";
            }
          }
          else {std::cout << "unknown functional group on (Csp2-imine) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
        else if ((numberC + numberH == 1)&&(numberO == 2)) {
          //carboxylate, carboxyl?
          if (Oneighbours == 1) {
            idetm = neighbours(idctm + 1,1) - 1;
            if (idetm < 0) {idetm = neighbours(iddtm + 1,1) - 1;}
            if ((atoms[idetm] == 1)||(atoms[idetm] == 6)) {
              atomtypes[idAtm] = "C-carboxyl";
              atomtypes[idctm] = "O-carboxyl";
              atomtypes[iddtm] = "O-carboxyl";
            }
            else {std::cout << "unknown functional group on (Csp2-carboxyl) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
          }
          else if (Oneighbours == 0) {
            atomtypes[idAtm] = "C-carboxylate";
            atomtypes[idctm] = "O-carboxylate";
            atomtypes[iddtm] = "O-carboxylate";
          }
          else {std::cout << "unknown functional group on (Csp2-carboxyl[ate]) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
        else if ((numberC + numberH == 1)&&(numberN == 1)&&(numberO == 1)&&(Oneighbours == 0)) {
          if (Nneighbours == 2) {
            if (NObond) {
              atomtypes[idAtm] = "C-hydroxyamide";
              if (atoms[idctm] == 7) {
                atomtypes[idctm] = "N-hydroxyamide";
                atomtypes[iddtm] = "O-hydroxyamide";
              }
              else {
                atomtypes[idctm] = "O-hydroxyamide";
                atomtypes[iddtm] = "N-hydroxyamide";
              }
            }
            else {
              atomtypes[idAtm] = "C-amide";
              if (atoms[idctm] == 7) {
                atomtypes[idctm] = "N-amide";
                atomtypes[iddtm] = "O-amide";
              }
              else {
                atomtypes[idctm] = "O-amide";
                atomtypes[iddtm] = "N-amide";
              }
            }
          }
          else {std::cout << "unknown functional group on (Csp2-amide) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
        else if ((numberN == 1)&&(numberO == 2)&&(Nneighbours == 2)) {
          atomtypes[idAtm] = "C-carbamate";
          if (atoms[idctm] == 7) {
            atomtypes[idctm] = "N-carbamate";
            atomtypes[iddtm] = "O-carbamate";
            atomtypes[idetm] = "O-carbamate";
          }
          else if (atoms[iddtm] == 7) {
            atomtypes[idctm] = "O-carbamate";
            atomtypes[iddtm] = "N-carbamate";
            atomtypes[idetm] = "O-carbamate";
          }
          else if (atoms[idetm] == 7) {
            atomtypes[idctm] = "O-carbamate";
            atomtypes[iddtm] = "O-carbamate";
            atomtypes[idetm] = "N-carbamate";
          }
        }
        else if (numberO == 3) {
          atomtypes[idAtm] = "C-carbonate";
          atomtypes[idctm] = "O-carbonate";
          atomtypes[iddtm] = "O-carbonate";
          atomtypes[idetm] = "O-carbonate";
        }
        else if (numberN == 3) {
          if (Nneighbours == 6) {
            atomtypes[idAtm] = "C-guanidinium";
            atomtypes[idctm] = "N-guanidinium";
            atomtypes[iddtm] = "N-guanidinium";
            atomtypes[idetm] = "N-guanidinium";
          }
        }
        else if ((numberN == 2)&&(numberC == 1)) {
          //check for number of protons attached to nitrogens
          numberH = 0;
          for (size_t idn = 0; idn < maxnneigh; ++idn) {
            if (neighbours(idctm + 1,idn + 1) == 0) {break;}
            if (neighbours(idctm + 1,idn + 1) - 1 == idAtm) {continue;}
            if (atoms[neighbours(idctm + 1,idn + 1) - 1] == 1) {++numberH;}
          }
          for (size_t idn = 0; idn < maxnneigh; ++idn) {
            if (neighbours(iddtm + 1,idn + 1) == 0) {break;}
            if (neighbours(iddtm + 1,idn + 1) - 1 == idAtm) {continue;}
            if (atoms[neighbours(iddtm + 1,idn + 1) - 1] == 1) {++numberH;}
          }
          if (numberH == 4) {std::cout << "group with potential charge" << std::endl;}
        }
        else if (numberP > 0) {atomtypes[idAtm] = "C-sp2";}      //for the time being do so
        else if (numberS > 0) {atomtypes[idAtm] = "C-sp2";}      //for the time being do so
        else if (numberH + numberC == 3) {atomtypes[idAtm] = "C-sp2";}
        else {std::cout << "unknown functional group on (Csp2) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
      }
      else if (nneighbours == 4) {                                          //sp3
        //if here, then the carbon and potential neighbours not yet assigned a group
        if (numberH + numberC == 4) {atomtypes[idAtm] = "C-sp3";}
        //other groups checked from the side of the other atoms, too many possibilities
      }
      else {std::cout << "unknown functional group on (C) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 7) {                                       //N
      if (nneighbours == 1) {
        idbtm = neighbours(idAtm + 1,1) - 1;
        if (atoms[idbtm] == 6) {continue;}                                    //this is nitrile or cyano, check from there
        //currently not checking for diazo
        else {std::cout << "unknown functional group on (Nsp) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
      }
      else if (nneighbours == 2) {
        //imines, enamines on the carbon side; other groups not considered
        idbtm = neighbours(idAtm + 1,1) - 1;
        idctm = neighbours(idAtm + 1,2) - 1;
        if (atoms[idbtm] == 8) {iddtm = idbtm;}
        else if (atoms[idctm] == 8) {iddtm = idctm;}
        if (atoms[iddtm] == 8) {
          nneighboursneighbours = 0;   //count neighbours of oxygen
          for (size_t idn = 0; idn < maxnneigh; ++idn) {
            idetm = neighbours(iddtm + 1,idn + 1) - 1;
            if (idetm < 0) {break;}
            ++nneighboursneighbours;
          }
          if (nneighboursneighbours == 1) {
            atomtypes[idAtm] = "N-nitroso";
            atomtypes[iddtm] = "O-nitroso";
          }
        }
      }
      else if (nneighbours == 3) {
        //get to know the neighbours
        idbtm = neighbours(idAtm + 1,1) - 1;
        idctm = neighbours(idAtm + 1,2) - 1;
        iddtm = neighbours(idAtm + 1,3) - 1;
        numberC = (atoms[idbtm] == 6) + (atoms[idctm] == 6) + (atoms[iddtm] == 6);
        numberH = (atoms[idbtm] == 1) + (atoms[idctm] == 1) + (atoms[iddtm] == 1);
        numberO = (atoms[idbtm] == 8) + (atoms[idctm] == 8) + (atoms[iddtm] == 8);
        nneighboursneighbours = 0;
        nneighboursneighbours2 = 0;
        nneighboursneighbours3 = 0;
        for (size_t idn = 0; idn < maxnneigh; ++idn) {
          idetm = neighbours(idbtm + 1,idn + 1) - 1;
          if (idetm < 0) {break;}
          ++nneighboursneighbours;
        }
        for (size_t idn = 0; idn < maxnneigh; ++idn) {
          idetm = neighbours(idctm + 1,idn + 1) - 1;
          if (idetm < 0) {break;}
          ++nneighboursneighbours2;
        }
        for (size_t idn = 0; idn < maxnneigh; ++idn) {
          idetm = neighbours(iddtm + 1,idn + 1) - 1;
          if (idetm < 0) {break;}
          ++nneighboursneighbours3;
        }
        if (numberH == 3) {atomtypes[idAtm] = "NH3";}
        else if ((numberC == 3)||((numberC == 2)&&(numberH == 1))||((numberC == 1)&&(numberH == 2))) {
          //check for sp3 carbons and protons (stored as if sp3 for conveniency)
          atom1sp3 = ((atoms[idbtm] == 6)&&(nneighboursneighbours == 4));
          if (atoms[idbtm] == 1) {atom1sp3 = true;}
          atom2sp3 = ((atoms[idctm] == 6)&&(nneighboursneighbours2 == 4));
          if (atoms[idctm] == 1) {atom2sp3 = true;}
          atom3sp3 = ((atoms[iddtm] == 6)&&(nneighboursneighbours3 == 4));
          if (atoms[iddtm] == 1) {atom3sp3 = true;}
          if (atom1sp3 + atom2sp3 + atom3sp3 == 3) {          //3 sp3 carbons around nitrogen, assign amine
            atomtypes[idAtm] = "N-amine";
            if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-amine";}
            if (atoms[idctm] == 6) {atomtypes[idctm] = "C-amine";}
            if (atoms[iddtm] == 6) {atomtypes[iddtm] = "C-amine";}
          }
          else if (atom1sp3 + atom2sp3 + atom3sp3 == 2) {continue;}          //this is potential amide, look from carbon side, easier there
          else {std::cout << "unknown functional group on (Namine/amide) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
        else if (numberO == 2) {
          atomtypes[idAtm] = "N-nitro";
          if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-nitro";}
          if (atoms[idctm] == 8) {atomtypes[idctm] = "O-nitro";}
          if (atoms[iddtm] == 8) {atomtypes[iddtm] = "O-nitro";}
        }
        else if (numberO == 3) {
          atomtypes[idAtm] = "N-nitrate";
          atomtypes[idbtm] = "O-nitrate";
          atomtypes[idctm] = "O-nitrate";
          atomtypes[iddtm] = "O-nitrate";
        }
        else {std::cout << "unknown functional group on (N) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
      }
      else if (nneighbours == 4) {
        atomtypes[idAtm] = "N+";
        for (size_t idn = 0; idn < 4; ++idn) {
          idbtm = neighbours(idAtm + 1,idn + 1) - 1;
          if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-N+";}
        }
      }
    }
    else if (atoms[idAtm] == 8) {                                       //O
      if (nneighbours == 1) {
        idbtm = neighbours(idAtm + 1,1) - 1;
        if (atomtypes[idbtm] + "blahblah" == "blahblah") {
          if (atoms[idbtm] == 1) {atomtypes[idAtm] = "hydroxide";}
          else if (atoms[idbtm] == 6) {
            //check for neighbour's neighbours
            numberH = 0;
            numberC = 0;
            numberN = 0;
            numberO = 0;
            numberother = 0;
            nneighboursneighbours = 0;        //how many neighbours does carbon neighbour have
            nneighboursneighbours2 = 0;
            nneighboursneighbours3 = 0;
            iddtm = -1;
            for (size_t idn = 0; idn < maxnneigh; ++idn) {
              if (neighbours(idbtm + 1,idn + 1) == 0) {break;}
              ++nneighboursneighbours;
              if (neighbours(idbtm + 1,idn + 1) - 1 == idAtm) {continue;}
              idctm = neighbours(idbtm + 1,idn + 1) - 1;
              if (atoms[idctm] == 1) {++numberH;}
              else if (atoms[idctm] == 6) {++numberC;}
              else if (atoms[idctm] == 7) {
                ++numberN;
                iddtm = idctm;
                for (size_t idn2 = 0; idn2 < maxnneigh; ++idn2) {
                  if (neighbours(idctm + 1,idn2 + 1) == 0) {break;}
                  ++nneighboursneighbours3;
                }
              }
              else if (atoms[idctm] == 8) {
                ++numberO;
                iddtm = idctm;
                for (size_t idn2 = 0; idn2 < maxnneigh; ++idn2) {
                  if (neighbours(idctm + 1,idn2 + 1) == 0) {break;}
                  if (neighbours(idctm + 1,idn2 + 1) - 1 == idbtm) {continue;}
                  ++nneighboursneighbours2;
                }
              }
              else {++numberother;}
            }
            if ((nneighboursneighbours == 3)&&(numberH == 1)&&(numberC == 1)) {
              atomtypes[idbtm] = "C-aldehyde";
              atomtypes[idAtm] = "O-aldehyde";
            }
            else if ((nneighboursneighbours == 3)&&(numberC == 2)) {
              atomtypes[idbtm] = "C-ketone";
              atomtypes[idAtm] = "O-ketone";
            }
            else if ((nneighboursneighbours == 3)&&(numberO == 1)&&(numberC + numberH == 1)) {
              if (nneighboursneighbours2 == 1) {
                idetm = neighbours(iddtm + 1,1) - 1;
                if ((atoms[idetm] == 1)||(atoms[idetm] == 6)) {
                  atomtypes[idAtm] = "O-carboxyl";
                  atomtypes[idbtm] = "C-carboxyl";
                  atomtypes[iddtm] = "O-carboxyl";
                }
                else {std::cout << "unknown functional group on (O-carboxyl) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
              }
              else if (nneighboursneighbours2 == 0) {
                atomtypes[idAtm] = "O-carboxylate";
                atomtypes[idbtm] = "C-carboxylate";
                atomtypes[iddtm] = "O-carboxylate";
              }
              else {std::cout << "unknown functional group on (O-carboxyl[ate]) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
            }
            else if ((nneighboursneighbours == 3)&&(numberN == 1)&&(numberC + numberH == 1)) {
              if (nneighboursneighbours3 == 3) {
                atomtypes[idAtm] = "O-amide";
                atomtypes[idbtm] = "C-amide";
                atomtypes[iddtm] = "N-amide";
              }
              else {std::cout << "unknown functional group on (O-amide) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
            }
            else if (nneighboursneighbours == 4) {
              if (numberH + numberC == 3) {
                atomtypes[idbtm] = "C-alkoxyde";
                atomtypes[idAtm] = "O-alkoxyde";
              }
              else {std::cout << "unknown functional group on (O-alkoxide) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
            }
            else {std::cout << "unknown functional group on (O1C) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
          }
          else {std::cout << "unknown functional group on (O1) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
        //else {} already assigned
      }
      else if (nneighbours == 2) {
        //who are the neighbours?
        idbtm = neighbours(idAtm + 1,1) - 1;
        idctm = neighbours(idAtm + 1,2) - 1;
        if (atoms[idbtm] > atoms[idctm]) {        //idbtm should contain lowest atomic number
          idbtm = neighbours(idAtm + 1,2) - 1;
          idctm = neighbours(idAtm + 1,1) - 1;
        }
        if ((atomtypes[idbtm] + "blahblah" == "blahblah")&&(atomtypes[idctm] + "blahblah" == "blahblah")) {
          //if nothing previously assigned then lets do something about it
          //look for neighbour's neighbours
          numberH = 0;
          numberC = 0;
          numberO = 0;
          numberother = 0;
          nneighboursneighbours3 = 0;
          for (size_t idn = 0; idn < maxnneigh; ++idn) {
            iddtm = neighbours(idbtm + 1,idn + 1) - 1;
            if (iddtm < 0) {break;}
            ++nneighboursneighbours3;
            if (iddtm == idAtm) {continue;}
            if (atoms[iddtm] == 1) {++numberH;}
            else if (atoms[iddtm] == 6) {++numberC;}
            else if (atoms[iddtm] == 8) {++numberO;}
            else {++numberother;}
          }
          numberH_2 = 0;
          numberC_2 = 0;
          numberO_2 = 0;
          numberother_2 = 0;
          nneighboursneighbours2 = 0;
          for (size_t idn = 0; idn < maxnneigh; ++idn) {
            iddtm = neighbours(idctm + 1,idn + 1) - 1;
            if (iddtm < 0) {break;}
            ++nneighboursneighbours2;
            if (iddtm == idAtm) {continue;}
            if (atoms[iddtm] == 1) {++numberH_2;}
            else if (atoms[iddtm] == 6) {++numberC_2;}
            else if (atoms[iddtm] == 8) {++numberO_2;}
            else {++numberother_2;}
          }
          if (atoms[idctm] == 1) {atomtypes[idAtm] = "OH2";}
          else if ((atoms[idbtm] == 1)&&(atoms[idctm] == 6)) {
            if ((nneighboursneighbours2 == 4)&&(numberH_2 + numberC_2 == 3)) {
              atomtypes[idctm] = "C-alcohol";
              atomtypes[idAtm] = "O-alcohol";
            }
            else if ((nneighboursneighbours2 == 3)&&(numberH_2 + numberC_2 == 1)&&(numberO_2 == 1)) {continue;}     //potential carboxylic acid
            else {std::cout << "unknown functional group on (Oalcohol) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
          }
          else if ((atoms[idbtm] == 6)&&(atoms[idctm] == 6)) {
            if ((nneighboursneighbours2 == 4)&&(numberH_2 + numberC_2 == 3)&&(nneighboursneighbours3 == 4)&&(numberH + numberC == 3)) {
              atomtypes[idbtm] = "C-ether";
              atomtypes[idctm] = "C-ether";
              atomtypes[idAtm] = "O-ether";
            }
            else if ((nneighboursneighbours2 == 3)&&(numberH_2 + numberC_2 == 1)&&(numberO_2 == 1)&&(nneighboursneighbours3 == 4)&&(numberH + numberC == 3)) {continue;}     //potential esther
            else if ((nneighboursneighbours2 == 4)&&(numberH_2 + numberC_2 == 3)&&(nneighboursneighbours3 == 3)&&(numberH + numberC == 1)&&(numberO == 1)) {continue;}     //potential esther
            else {std::cout << "unknown functional group on (Oether) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
          }
          else if ((atoms[idbtm] == 6)&&(atoms[idctm] == 8)) {
            if ((nneighboursneighbours3 == 4)&&(numberH + numberC == 3)&&(nneighboursneighbours2 < 3)) {
              atomtypes[idbtm] = "C-peroxide";
              atomtypes[idctm] = "O-peroxide";
              atomtypes[idAtm] = "O-peroxide";
              //count neighbours of the second oxygen
              idbtm = -1;
              nneighboursneighbours2 = 0;
              for (size_t idn2 = 0; idn2 < maxnneigh; ++idn2) {
                if (neighbours(idctm + 1,idn2 + 1) == 0) {break;}
                if (neighbours(idctm + 1,idn2 + 1) - 1 != idAtm) {
                  ++nneighboursneighbours2;
                  idbtm = neighbours(idctm + 1,idn2 + 1) - 1;
                }
              }
              if (nneighboursneighbours2 == 1) {
                if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-peroxide";}
              }
            }
            else {std::cout << "unknown functional group on (Operoxide) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
          }
          else {std::cout << "unknown functional group on (O2) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
      }
      else if (nneighbours == 3) {atomtypes[idAtm] = "O+";}
      else if (nneighbours == 4) {std::cout << "unknown oxygen with four connections";}
    }
    else if (atoms[idAtm] == 9) {                                       //F
      if ((nneighbours == 1)&&(atoms[neighbours(idAtm + 1,1) - 1] == 6)) {
        atomtypes[idAtm] = "F-haloalkane";
        atomtypes[neighbours(idAtm + 1,1) - 1] = "C-haloalkane";
      }
      else {std::cout << "unknown functional group on (F) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 10) {atomtypes[idAtm] = "Ne";}             //Ne
    else if (atoms[idAtm] == 11) {                                      //Na
      if (nneighbours == 0) {atomtypes[idAtm] = "Na";}                        //no assumption on charge
      else if (nneighbours == 1) {
        idbtm = atoms[neighbours(idAtm + 1,1) - 1];
        //based on electronegativities
        if ((idbtm == 11)||(idbtm == 19)||(idbtm == 37)||(idbtm == 55)||(idbtm == 56)||(idbtm == 87)||(idbtm == 88)) {atomtypes[idAtm] = "Na";}
        else {atomtypes[idAtm] = "Na+";}
      }
      else {std::cout << "unknown functional group on (Na) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 12) {                                      //Mg
      if (nneighbours == 0) {atomtypes[idAtm] = "Mg";}                        //no assumption on charge
      else if (nneighbours == 1) {
        idbtm = atoms[neighbours(idAtm + 1,1) - 1];
        //based on electronegativities
        if ((idbtm == 3)||(idbtm == 11)||(idbtm == 12)||(idbtm == 19)||(idbtm == 20)||(idbtm == 37)||(idbtm == 38)||(idbtm == 55)||(idbtm == 56)||(idbtm == 87)||(idbtm == 88)) {atomtypes[idAtm] = "Mg";}
        else {atomtypes[idAtm] = "Mg++";}
      }
      else {std::cout << "unknown functional group on (Mg) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 13) {                                      //Al
      atomtypes[idAtm] = "Al";
    }
    else if (atoms[idAtm] == 14) {                                      //Si
      atomtypes[idAtm] = "Si";
    }
    else if (atoms[idAtm] == 15) {                                      //P
      if (nneighbours == 4) {
        numberC = 0;
        numberN = 0;
        numberO = 0;
        numberS = 0;
        protonatedoxygens = false;
        for (size_t idn = 1; idn < 5; ++idn) {
          idbtm = neighbours(idAtm + 1,idn) - 1;
          if (atoms[idbtm] == 6) {++numberC;}
          if (atoms[idbtm] == 7) {++numberN;}
          else if (atoms[idbtm] == 8) {
            ++numberO;
            //check for protonated phosphates
            if (!protonatedoxygens) {
              for (size_t idn2 = 0; idn2 < maxnneigh; ++idn2) {
                if (neighbours(idbtm + 1,idn2 + 1) == 0) {break;}
                if (atoms[neighbours(idbtm + 1,idn2 + 1) - 1] == 1) {
                  protonatedoxygens = true;
                  break;
                }
              }
            }
          }
          if (atoms[idbtm] == 16) {++numberS;}
        }
        if (numberO == 4) {
          atomtypes[idAtm] = "P-phosphate";
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            atomtypes[idbtm] = "O-phosphate";
          }
          if (protonatedoxygens) {std::cout << "potentially protonated phosphate: " << idAtm + 1 << "\n";}
        }
        else if ((numberO == 3)&&(numberC == 1)) {
          atomtypes[idAtm] = "P-phosphonate";
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-phosphonate";}
            else if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-phosphonate";}
          }
        }
        else if ((numberO == 2)&&(numberC == 2)) {
          atomtypes[idAtm] = "P-phosphinate";
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-phosphinate";}
            else if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-phosphinate";}
          }
        }
        else if ((numberO == 1)&&(numberC == 3)) {
          atomtypes[idAtm] = "P-phosphoryl";
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-phosphoryl";}
            else if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-phosphoryl";}
          }
        }
        else if ((numberN == 1)&&(numberO == 3)) {
          atomtypes[idAtm] = "P-aminophosphate";
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            if (atoms[idbtm] == 7) {atomtypes[idbtm] = "N-aminophosphate";}
            else if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-aminophosphate";}
          }
        }
        else if ((numberN == 1)&&(numberO == 2)&&(numberC == 1)) {
          atomtypes[idAtm] = "P-aminophosphinate";
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-aminophosphinate";}
            else if (atoms[idbtm] == 7) {atomtypes[idbtm] = "N-aminophosphinate";}
            else if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-aminophosphinate";}
          }
        }
        else if ((numberN == 1)&&(numberO == 1)&&(numberC == 2)) {
          atomtypes[idAtm] = "P-aminophosphoryl";
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-aminophosphoryl";}
            else if (atoms[idbtm] == 7) {atomtypes[idbtm] = "N-aminophosphoryl";}
            else if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-aminophosphoryl";}
          }
        }
        else if ((numberS == 1)&&(numberO == 3)) {
          atomtypes[idAtm] = "P-thiophosphate";
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-thiophosphate";}
            else if (atoms[idbtm] == 16) {atomtypes[idbtm] = "S-thiophosphate";}
          }
        }
        else if ((numberS == 1)&&(numberO == 2)&&(numberN == 1)) {
          atomtypes[idAtm] = "P-aminothiophosphate";
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            if (atoms[idbtm] == 7) {atomtypes[idbtm] = "N-aminothiophosphate";}
            else if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-aminothiophosphate";}
            else if (atoms[idbtm] == 16) {atomtypes[idbtm] = "S-aminothiophosphate";}
          }
        }
        else {std::cout << "unknown functional group on (P4) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
      }
      else {std::cout << "unknown functional group on (P) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 16) {                                      //S
      if (nneighbours == 2) {
        //look for thiols and sulfides
        idbtm = neighbours(idAtm + 1,1) - 1;
        idctm = neighbours(idAtm + 1,2) - 1;
        if (atoms[idbtm] > atoms[idctm]) {
          idbtm = neighbours(idAtm + 1,2) - 1;
          idctm = neighbours(idAtm + 1,1) - 1;
        }
        if (((atoms[idbtm] == 1)||(atoms[idbtm] == 6))&&(atoms[idctm] == 6)) {
          //look for neighbours of idbtm; can only be C or H, otherwise no real thiol of sulfide
          numberHC = 0;
          nneighboursneighbours = 0;
          for (size_t idn = 0; idn < maxnneigh; ++idn) {
            if (neighbours(idbtm + 1,idn + 1) == 0) {break;}
            if (neighbours(idbtm + 1,idn + 1) - 1 == idAtm) {continue;}
            if ((atoms[neighbours(idbtm + 1,idn + 1) - 1] == 1)||(atoms[neighbours(idbtm + 1,idn + 1) - 1] == 6)) {++numberHC;}
            ++nneighboursneighbours;
          }
          if ((numberHC == nneighboursneighbours)&&(nneighboursneighbours == 3)) {
            if (atoms[idctm] == 1) {
              atomtypes[idAtm] = "S-thiol";
              atomtypes[idbtm] = "C-thiol";
            }
            else {                                               //then idctm contains carbon
              //look for neighbours of idctm
              numberHC = 0;
              nneighboursneighbours = 0;
              for (size_t idn = 0; idn < maxnneigh; ++idn) {
                if (neighbours(idctm + 1,idn + 1) == 0) {break;}
                if (neighbours(idctm + 1,idn + 1) - 1 == idAtm) {continue;}
                if ((atoms[neighbours(idctm + 1,idn + 1) - 1] == 1)||(atoms[neighbours(idctm + 1,idn + 1) - 1] == 6)) {++numberHC;}
                ++nneighboursneighbours;
              }
              if ((numberHC == nneighboursneighbours)&&(nneighboursneighbours == 3)) {
                atomtypes[idAtm] = "S-sulfide";
                atomtypes[idbtm] = "C-sulfide";
                atomtypes[idctm] = "C-sulfide";
              }
              else {std::cout << "unknown functional group on (Ssulfide) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
            }
          }
          else {std::cout << "unknown functional group on (Sthiol/sulfide) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
        }
        else if (atoms[idctm] == 1) {atomtypes[idAtm] = "SH2";}
        else if (((atoms[idbtm] == 1)||(atoms[idbtm] == 6))&&(atoms[idctm] == 16)) {
          atomtypes[idAtm] = "S-disulfide";
          atomtypes[idctm] = "S-disulfide";
          if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-disulfide";}
          //count neighbours of the second sulfur
          idbtm = -1;
          nneighboursneighbours2 = 0;
          for (size_t idn2 = 0; idn2 < maxnneigh; ++idn2) {
            if (neighbours(idctm + 1,idn2 + 1) == 0) {break;}
            if (neighbours(idctm + 1,idn2 + 1) - 1 != idAtm) {
              ++nneighboursneighbours2;
              idbtm = neighbours(idctm + 1,idn2 + 1) - 1;
            }
          }
          if (nneighboursneighbours2 == 1) {
            if (atoms[idbtm] == 6) {atomtypes[idbtm] = "C-disulfide";}
          }
        }
        else {std::cout << "unknown functional group on (S2) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
      }
      else if (nneighbours == 3) {
        idbtm = neighbours(idAtm + 1,1) - 1;
        idctm = neighbours(idAtm + 1,2) - 1;
        iddtm = neighbours(idAtm + 1,3) - 1;
        numberC = (atoms[idbtm] == 6) + (atoms[idctm] == 6) + (atoms[iddtm] == 6);
        numberO = (atoms[idbtm] == 8) + (atoms[idctm] == 8) + (atoms[iddtm] == 8);
        if (numberO == 1) {
          Oneighbours = 0;
          idetm = idbtm;
          if (atoms[idctm] == 8) {idetm = idctm;}
          else if (atoms[iddtm] == 8) {idetm = iddtm;}
          for (size_t idn = 0; idn < maxnneigh; ++idn) {
            if (neighbours(idetm + 1,idn + 1) == 0) {break;}
            if (!(neighbours(idetm + 1,idn + 1) - 1 == idAtm)) {++Oneighbours;}
          }
          if ((Oneighbours == 0)&&(numberC == 2)) {
            atomtypes[idAtm] = "S-sulfoxide";
            for (size_t idn = 0; idn < maxnneigh; ++idn) {
              if (neighbours(idAtm + 1,idn + 1) == 0) {break;}
              if (atoms[neighbours(idAtm + 1,idn + 1) - 1] == 6) {atomtypes[idAtm] = "C-sp3";}
              else if (atoms[neighbours(idAtm + 1,idn + 1) - 1] == 8) {atomtypes[idAtm] = "O-sulfoxide";}
            }
          }
        }
        else if (numberC == 3) {
          atomtypes[idAtm] = "S+";
          atomtypes[idbtm] = "C-sp3";
          atomtypes[idctm] = "C-sp3";
          atomtypes[iddtm] = "C-sp3";
        }
      }
      else if (nneighbours == 4) {
        //look for sulfonates and sulfates
        numberO = 0;
        protonatedoxygens = false;
        nneighboursneighbours = 0;
        for (size_t idn = 1; idn < 5; ++idn) {
          idbtm = neighbours(idAtm + 1,idn) - 1;
          if (atoms[idbtm] == 8) {
            ++numberO;
            //check for protonated sulfonate and/or neighbours
            for (size_t idn2 = 0; idn2 < maxnneigh; ++idn2) {
              if (neighbours(idbtm + 1,idn2 + 1) == 0) {break;}
              if (neighbours(idbtm + 1,idn2 + 1) - 1 == idAtm) {continue;}
              if (atoms[neighbours(idbtm + 1,idn2 + 1) - 1] == 1) {protonatedoxygens = true;}
              ++nneighboursneighbours;
            }
          }
        }
        if (numberO == 3) {
          atomtypes[idAtm] = "S-sulfonate";
          //now the oxygens
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            if (atoms[idbtm] == 8) {atomtypes[idbtm] = "O-sulfonate";}
          }
          if (protonatedoxygens) {std::cout << "potentially protonated sulfonate: " << idAtm + 1 << "\n";}
        }
        else if (numberO == 4) {
          atomtypes[idAtm] = "S-sulfate";
          //now the oxygens
          for (size_t idn = 1; idn < 5; ++idn) {
            idbtm = neighbours(idAtm + 1,idn) - 1;
            atomtypes[idbtm] = "O-sulfate";
          }
          if (protonatedoxygens) {std::cout << "potential protonated sulfate: " << idAtm + 1 << "\n";}
        }
        else {std::cout << "unknown functional group on (S4) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
      }
      else {std::cout << "unknown functional group on (S) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 17) {                                      //Cl
      if ((nneighbours == 1)&&(atoms[neighbours(idAtm + 1,1) - 1] == 6)) {
        atomtypes[idAtm] = "Cl-haloalkane";
        atomtypes[neighbours(idAtm + 1,1) - 1] = "C-haloalkane";
      }
      else {std::cout << "unknown functional group on (Cl) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 18) {atomtypes[idAtm] = "Ar";}             //Ar
    else if (atoms[idAtm] == 19) {                                      //K
      if (nneighbours == 0) {atomtypes[idAtm] = "K";}                         //no assumption on charge
      else if (nneighbours == 1) {
        idbtm = atoms[neighbours(idAtm + 1,1) - 1];
        //based on electronegativities
        if ((idbtm == 19)||(idbtm == 37)||(idbtm == 55)||(idbtm == 87)) {atomtypes[idAtm] = "K";}
        else {atomtypes[idAtm] = "K+";}
      }
      else {std::cout << "unknown functional group on (K) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 20) {                                      //Ca
      if (nneighbours == 0) {atomtypes[idAtm] = "Ca";}                        //no assumption on charge
      else if (nneighbours == 1) {
        idbtm = atoms[neighbours(idAtm + 1,1) - 1];
        //based on electronegativities
        if ((idbtm == 3)||(idbtm == 11)||(idbtm == 19)||(idbtm == 20)||(idbtm == 37)||(idbtm == 38)||(idbtm == 55)||(idbtm == 56)||(idbtm == 87)||(idbtm == 88)) {atomtypes[idAtm] = "Ca";}
        else {atomtypes[idAtm] = "Ca++";}
      }
      else {std::cout << "unknown functional group on (Ca) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 21) {                                      //Sc
      atomtypes[idAtm] = "Sc";
    }
    else if (atoms[idAtm] == 22) {                                      //Ti
      atomtypes[idAtm] = "Ti";
    }
    else if (atoms[idAtm] == 23) {                                      //V
      atomtypes[idAtm] = "V";
    }
    else if (atoms[idAtm] == 24) {                                      //Cr
      atomtypes[idAtm] = "Cr";
    }
    else if (atoms[idAtm] == 25) {                                      //Mn
      atomtypes[idAtm] = "Mn";
    }
    else if (atoms[idAtm] == 26) {                                      //Fe
      atomtypes[idAtm] = "Fe";
    }
    else if (atoms[idAtm] == 27) {                                      //Co
      atomtypes[idAtm] = "Co";
    }
    else if (atoms[idAtm] == 28) {                                      //Ni
      atomtypes[idAtm] = "Ni";
    }
    else if (atoms[idAtm] == 29) {                                      //Cu
      atomtypes[idAtm] = "Cu";
    }
    else if (atoms[idAtm] == 30) {                                      //Zn
      atomtypes[idAtm] = "Zn";
    }
    else if (atoms[idAtm] == 31) {                                      //Ga
      atomtypes[idAtm] = "Ga";
    }
    else if (atoms[idAtm] == 32) {                                      //Ge
      atomtypes[idAtm] = "Ge";
    }
    else if (atoms[idAtm] == 33) {                                      //As
      atomtypes[idAtm] = "As";
    }
    else if (atoms[idAtm] == 34) {                                      //Se
      atomtypes[idAtm] = "Se";
    }
    else if (atoms[idAtm] == 35) {                                      //Br
      if ((nneighbours == 1)&&(atoms[neighbours(idAtm + 1,1) - 1] == 6)) {
        atomtypes[idAtm] = "Br-haloalkane";
        atomtypes[neighbours(idAtm + 1,1) - 1] = "C-haloalkane";
      }
      else {std::cout << "unknown functional group on (Br) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 36) {atomtypes[idAtm] = "Kr";}             //Kr
    else if (atoms[idAtm] == 37) {                                      //Rb
      if (nneighbours == 0) {atomtypes[idAtm] = "Rb";}                        //no assumption on charge
      else if (nneighbours == 1) {
        idbtm = atoms[neighbours(idAtm + 1,1) - 1];
        //based on electronegativities
        if ((idbtm == 37)||(idbtm == 55)||(idbtm == 87)) {atomtypes[idAtm] = "Rb";}
        else {atomtypes[idAtm] = "Rb+";}
      }
      else {std::cout << "unknown functional group on (Rb) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 38) {                                      //Sr
      if (nneighbours == 0) {atomtypes[idAtm] = "Sr";}                        //no assumption on charge
      else if (nneighbours == 1) {
        idbtm = atoms[neighbours(idAtm + 1,1) - 1];
        //based on electronegativities
        if ((idbtm == 11)||(idbtm == 19)||(idbtm == 37)||(idbtm == 38)||(idbtm == 55)||(idbtm == 56)||(idbtm == 87)||(idbtm == 88)) {atomtypes[idAtm] = "Sr";}
        else {atomtypes[idAtm] = "Sr++";}
      }
      else {std::cout << "unknown functional group on (Sr) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 39) {                                      //Y
      atomtypes[idAtm] = "Y";
    }
    else if (atoms[idAtm] == 40) {                                      //Zr
      atomtypes[idAtm] = "Zr";
    }
    else if (atoms[idAtm] == 41) {                                      //Nb
      atomtypes[idAtm] = "Nb";
    }
    else if (atoms[idAtm] == 42) {                                      //Mo
      atomtypes[idAtm] = "Mo";
    }
    else if (atoms[idAtm] == 43) {                                      //Tc
      atomtypes[idAtm] = "Tc";
    }
    else if (atoms[idAtm] == 44) {                                      //Ru
      atomtypes[idAtm] = "Ru";
    }
    else if (atoms[idAtm] == 45) {                                      //Rh
      atomtypes[idAtm] = "Rh";
    }
    else if (atoms[idAtm] == 46) {                                      //Pd
      atomtypes[idAtm] = "Pd";
    }
    else if (atoms[idAtm] == 47) {                                      //Ag
      atomtypes[idAtm] = "Ag";
    }
    else if (atoms[idAtm] == 48) {                                      //Cd
      atomtypes[idAtm] = "Cd";
    }
    else if (atoms[idAtm] == 49) {                                      //In
      atomtypes[idAtm] = "In";
    }
    else if (atoms[idAtm] == 50) {                                      //Sn
      atomtypes[idAtm] = "Sn";
    }
    else if (atoms[idAtm] == 51) {                                      //Sb
      atomtypes[idAtm] = "Sb";
    }
    else if (atoms[idAtm] == 52) {                                      //Te
      atomtypes[idAtm] = "Te";
    }
    else if (atoms[idAtm] == 53) {                                      //I
      if ((nneighbours == 1)&&(atoms[neighbours(idAtm + 1,1) - 1] == 6)) {
        atomtypes[idAtm] = "I-haloalkane";
        atomtypes[neighbours(idAtm + 1,1) - 1] = "C-haloalkane";
      }
      else {std::cout << "unknown functional group on (I) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 54) {atomtypes[idAtm] = "Xe";}             //Xe
    else if (atoms[idAtm] == 55) {                                      //Cs
      if (nneighbours == 0) {atomtypes[idAtm] = "Cs";}                        //no assumption on charge
      else if (nneighbours == 1) {
        idbtm = atoms[neighbours(idAtm + 1,1) - 1];
        //based on electronegativities
        if ((idbtm == 55)||(idbtm == 87)) {atomtypes[idAtm] = "Cs";}
        else {atomtypes[idAtm] = "Cs+";}
      }
      else {std::cout << "unknown functional group on (Cs) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 56) {                                      //Ba
      if (nneighbours == 0) {atomtypes[idAtm] = "Ba";}                        //no assumption on charge
      else if (nneighbours == 1) {
        idbtm = atoms[neighbours(idAtm + 1,1) - 1];
        //based on electronegativities
        if ((idbtm == 19)||(idbtm == 37)||(idbtm == 55)||(idbtm == 56)||(idbtm == 87)||(idbtm == 88)) {atomtypes[idAtm] = "Ba";}
        else {atomtypes[idAtm] = "Ba++";}
      }
      else {std::cout << "unknown functional group on (Ba) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
    }
    else if (atoms[idAtm] == 57) {                                      //La
      atomtypes[idAtm] = "La";
    }
    else if (atoms[idAtm] == 58) {                                      //Ce
      atomtypes[idAtm] = "Ce";
    }
    else if (atoms[idAtm] == 59) {                                      //Pr
      atomtypes[idAtm] = "Pr";
    }
    else if (atoms[idAtm] == 60) {                                      //Nd
      atomtypes[idAtm] = "Nd";
    }
    else if (atoms[idAtm] == 61) {                                      //Pm
      atomtypes[idAtm] = "Pm";
    }
    else if (atoms[idAtm] == 62) {                                      //Sm
      atomtypes[idAtm] = "Sm";
    }
    else if (atoms[idAtm] == 63) {                                      //Eu
      atomtypes[idAtm] = "Eu";
    }
    else if (atoms[idAtm] == 64) {                                      //Gd
      atomtypes[idAtm] = "Gd";
    }
    else if (atoms[idAtm] == 65) {                                      //Tb
      atomtypes[idAtm] = "Tb";
    }
    else if (atoms[idAtm] == 66) {                                      //Dy
      atomtypes[idAtm] = "Dy";
    }
    else if (atoms[idAtm] == 67) {                                      //Ho
      atomtypes[idAtm] = "Ho";
    }
    else if (atoms[idAtm] == 68) {                                      //Er
      atomtypes[idAtm] = "Er";
    }
    else if (atoms[idAtm] == 69) {                                      //Tm
      atomtypes[idAtm] = "Tm";
    }
    else if (atoms[idAtm] == 70) {                                      //Yb
      atomtypes[idAtm] = "Yb";
    }
    else if (atoms[idAtm] == 71) {                                      //Lu
      atomtypes[idAtm] = "Lu";
    }
    else if (atoms[idAtm] == 72) {                                      //Hf
      atomtypes[idAtm] = "Hf";
    }
    else if (atoms[idAtm] == 73) {                                      //Ta
      atomtypes[idAtm] = "Ta";
    }
    else if (atoms[idAtm] == 74) {                                      //W
      atomtypes[idAtm] = "W";
    }
    else if (atoms[idAtm] == 75) {                                      //Re
      atomtypes[idAtm] = "Re";
    }
    else if (atoms[idAtm] == 76) {                                      //Os
      atomtypes[idAtm] = "Os";
    }
    else if (atoms[idAtm] == 77) {                                      //Ir
      atomtypes[idAtm] = "Ir";
    }
    else if (atoms[idAtm] == 78) {                                      //Pt
      atomtypes[idAtm] = "Pt";
    }
    else if (atoms[idAtm] == 79) {                                      //Au
      atomtypes[idAtm] = "Au";
    }
    else if (atoms[idAtm] == 80) {                                      //Hg
      atomtypes[idAtm] = "Hg";
    }
    else if (atoms[idAtm] == 81) {                                      //Tl
      atomtypes[idAtm] = "Tl";
    }
    else if (atoms[idAtm] == 82) {                                      //Pb
      atomtypes[idAtm] = "Pb";
    }
    else if (atoms[idAtm] == 83) {                                      //Bi
      atomtypes[idAtm] = "Bi";
    }
    else if (atoms[idAtm] == 84) {                                      //Po
      atomtypes[idAtm] = "Po";
    }
    else if (atoms[idAtm] == 85) {                                      //At
      atomtypes[idAtm] = "At";
    }
    else if (atoms[idAtm] == 86) {atomtypes[idAtm] = "Rn";}             //Rn
    else if (atoms[idAtm] == 87) {                                      //Fr
      atomtypes[idAtm] = "Fr";
    }
    else if (atoms[idAtm] == 88) {                                      //Ra
      atomtypes[idAtm] = "Ra";
    }
    else if (atoms[idAtm] == 89) {                                      //Ac
      atomtypes[idAtm] = "Ac";
    }
    else {std::cout << "unknown functional group on (Atom) " << atoms[idAtm] << "(" << idAtm + 1 << ") \n";}
  }
  //second run to fix a few missing groups
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    if (atoms[idAtm] == 1) {
      atomtypes[idAtm] = "H";
      continue;
    }
    else if ((atomtypes[idAtm] == "")&&(atoms[idAtm] == 6)) {      //fix missing carbon assignments for now
      idbtm = neighbours(idAtm + 1,1) - 1;
      idctm = neighbours(idAtm + 1,2) - 1;
      if (idctm > 0) {
        if (atoms[idctm] < atoms[idbtm]) {swap(idctm,idbtm);}
        iddtm = neighbours(idAtm + 1,3) - 1;
        if (iddtm > 0) {
          if (atoms[iddtm] < atoms[idctm]) {swap(iddtm,idctm);}
          if (atoms[idctm] < atoms[idbtm]) {swap(idctm,idbtm);}
          idetm = neighbours(idAtm + 1,3) - 1;
          if (idetm > 0) {
            if (atoms[idetm] < atoms[iddtm]) {swap(idetm,iddtm);}
            if (atoms[iddtm] < atoms[idctm]) {swap(iddtm,idctm);}
            if (atoms[idctm] < atoms[idbtm]) {swap(idctm,idbtm);}
            nneighbours = 4;
          }
          else {nneighbours = 3;}
        }
        else {nneighbours = 2;}
      }
      else {nneighbours = 1;}
      if (nneighbours == 4) {
        if (atoms[idetm] == 7) {
          if (atomtypes[idetm] == "N-amide") {atomtypes[idAtm] = "C-sp3";}
          else {std::cout << "unknown functional group on (CN) " << atoms[idAtm] << "(" << idAtm + 1 << ")* \n";}
        }
        else if (atoms[idetm] == 8) {
          if (atomtypes[idetm] == "O-carboxyl") {atomtypes[idAtm] = "C-sp3";}
          else if (atomtypes[idetm] == "C-peroxide") {atomtypes[idAtm] = "C-sp3";}
          else if (atomtypes[idetm] == "O-phosphate") {atomtypes[idAtm] = "C-sp3";}
          else if (atomtypes[idetm] == "O-aminophosphate") {atomtypes[idAtm] = "C-sp3";}
          else if (atomtypes[idetm] == "O-thiophosphate") {atomtypes[idAtm] = "C-sp3";}
          else if (atomtypes[idetm] == "O-aminothiophosphate") {atomtypes[idAtm] = "C-sp3";}
          else if (atomtypes[idetm] == "O-sulfate") {atomtypes[idAtm] = "C-sp3";}
          else {std::cout << "unknown functional group on (CO) " << atoms[idAtm] << "(" << idAtm + 1 << ")* \n";}
        }
      }
      else if (nneighbours == 3) {}
      else {std::cout << "unknown functional group on (C) " << atoms[idAtm] << "(" << idAtm + 1 << ")* \n";}
    }
    else if (atomtypes[idAtm] == "C-enol") {
      //fix on enols -> two carbons
      idbtm = neighbours(idAtm + 1,1) - 1;
      idctm = neighbours(idAtm + 1,2) - 1;
      iddtm = neighbours(idAtm + 1,3) - 1;
      numberHC = (atomtypes[idbtm] == "C-sp2") + (atomtypes[idctm] == "C-sp2") + (atomtypes[iddtm] == "C-sp2");     //number of sp2 carbons around enol carbon
      numberO = (atomtypes[idbtm] == "O-enol") + (atomtypes[idctm] == "O-enol") + (atomtypes[iddtm] == "O-enol");   //number of sp2 carbons around enol carbon
      if ((numberHC == 1)&&(numberO == 1)) {
        if (atomtypes[idbtm] == "C-sp2") {atomtypes[idbtm] = "C-enol";}
        else if (atomtypes[idctm] == "C-sp2") {atomtypes[idctm] = "C-enol";}
        else if (atomtypes[iddtm] == "C-sp2") {atomtypes[iddtm] = "C-enol";}
      }
      else if ((numberHC == 2)&&(numberO == 1)) {
        nneighboursneighbours2 = 0;         //number of sp2 carbons
        nneighboursneighbours3 = 0;         //number of sp2 carbons
        if (atomtypes[iddtm] == "C-sp2") {        //put sp2 carbons on idbtm and idctm
          if (atomtypes[idctm] != "C-sp2") {idctm = iddtm;}
          else {idbtm = iddtm;}
        }
        for (size_t idn = 0; idn < 3; ++idn) {
          iddtm = neighbours(idbtm + 1,idn + 1) - 1;
          idetm = neighbours(idctm + 1,idn + 1) - 1;
          if (atomtypes[iddtm] == "C-sp2") {++nneighboursneighbours2;}
          if (atomtypes[idetm] == "C-sp2") {++nneighboursneighbours3;}
        }
        if ((nneighboursneighbours2 > 0)&&(nneighboursneighbours3 == 0)) {atomtypes[idctm] = "C-enol";}
        else if ((nneighboursneighbours2 == 0)&&(nneighboursneighbours3 > 0)) {atomtypes[idbtm] = "C-enol";}
        else {std::cout << "unknown functional group on (C-enol) " << atoms[idAtm] << "(" << idAtm + 1 << ")* \n";}
      }
      //else {}   no longer enol
    }
  }
  GetGroupEvaluation(atomtypes);
  return atomtypes;
}
int TotalChargeEstimator(const std::vector<std::string> & funcgroups, const std::vector<size_t> & atoms, matrix<int> & neighbours) {
  //function that estimates the total charge of a system based on the functional groups
  int tcharge = 0;
  int nneighbours;
  size_t Natoms = funcgroups.size();
  size_t maxneighbours = neighbours.cols();
  int atmA;
  int atmB;
  for (size_t idatm = 0; idatm < Natoms; ++idatm) {
    if (funcgroups[idatm] == "O-alkoxyde") {--tcharge;}
    else if (funcgroups[idatm] == "Li+") {++tcharge;}
    else if (funcgroups[idatm] == "B-") {--tcharge;}
    else if (funcgroups[idatm] == "N+") {++tcharge;}
    else if (funcgroups[idatm] == "C-guanidinium") {++tcharge;}
    else if (funcgroups[idatm] == "hydroxide") {--tcharge;}
    else if (funcgroups[idatm] == "O+") {++tcharge;}
    else if (funcgroups[idatm] == "C-carboxylate") {--tcharge;}
    else if ((funcgroups[idatm] == "P-phosphate")||(funcgroups[idatm] == "P-thiophosphate")) {           //further analyse
      nneighbours = 0;           //number of neighbours of oxygen atoms
      for (size_t idn = 0; idn < maxneighbours; ++idn) {
        atmA = neighbours(idatm + 1,idn + 1) - 1;
        if (atmA < 0) {break;}
        if (atoms[atmA] != 8) {continue;}
        //check for neighbour's neighbours
        for (size_t idn2 = 0; idn2 < maxneighbours; ++idn2) {
          atmB = neighbours(atmA + 1,idn2 + 1) - 1;
          if (atmB < 0) {break;}
          if (atmB == idatm) {continue;}
          ++nneighbours;
        }
      }
      if (nneighbours == 0) {tcharge -= 3;}
      else if (nneighbours == 1) {tcharge -= 2;}
      else if (nneighbours == 2) {tcharge -= 1;}
    }
    else if ((funcgroups[idatm] == "P-phosphonate")||(funcgroups[idatm] == "P-aminophosphate")||(funcgroups[idatm] == "P-aminothiophosphate")) {           //further analyse
      nneighbours = 0;           //number of neighbours of oxygen atoms
      for (size_t idn = 0; idn < maxneighbours; ++idn) {
        atmA = neighbours(idatm + 1,idn + 1) - 1;
        if (atmA < 0) {break;}
        //check for neighbour's neighbours
        if (atoms[atmA] == 8) {
          for (size_t idn2 = 0; idn2 < maxneighbours; ++idn2) {
            atmB = neighbours(atmA + 1,idn2 + 1) - 1;
            if (atmB < 0) {break;}
            if (atmB == idatm) {continue;}
            ++nneighbours;
          }
        }
      }
      if (nneighbours == 0) {tcharge -= 2;}
      else if (nneighbours == 1) {--tcharge;}
    }
    else if ((funcgroups[idatm] == "P-phosphinate")||(funcgroups[idatm] == "P-aminophosphinate")) {           //further analyse
      nneighbours = 0;           //number of neighbours of oxygen atoms
      for (size_t idn = 0; idn < maxneighbours; ++idn) {
        atmA = neighbours(idatm + 1,idn + 1) - 1;
        if (atmA < 0) {break;}
        //check for neighbour's neighbours
        if (atoms[atmA] == 8) {
          for (size_t idn2 = 0; idn2 < maxneighbours; ++idn2) {
            atmB = neighbours(atmA + 1,idn2 + 1) - 1;
            if (atmB < 0) {break;}
            if (atmB == idatm) {continue;}
            ++nneighbours;
          }
        }
      }
      if (nneighbours == 0) {--tcharge;}
    }
    else if (funcgroups[idatm] == "S-sulfonate") {           //further analyse
      nneighbours = 0;           //number of neighbours of oxygen atoms
      for (size_t idn = 0; idn < maxneighbours; ++idn) {
        atmA = neighbours(idatm + 1,idn + 1) - 1;
        if (atmA < 0) {break;}
        if (atoms[atmA] != 8) {continue;}
        //check for neighbour's neighbours
        for (size_t idn2 = 0; idn2 < maxneighbours; ++idn2) {
          atmB = neighbours(atmA + 1,idn2 + 1) - 1;
          if (atmB < 0) {break;}
          if (atmB == idatm) {continue;}
          ++nneighbours;
        }
      }
      if (nneighbours == 0) {--tcharge;}
    }
    else if (funcgroups[idatm] == "S-sulfate") {             //further analyse
      nneighbours = 0;           //number of neighbours of oxygen atoms
      for (size_t idn = 0; idn < maxneighbours; ++idn) {
        atmA = neighbours(idatm + 1,idn + 1) - 1;
        if (atmA < 0) {break;}
        //check for neighbour's neighbours
        for (size_t idn2 = 0; idn2 < maxneighbours; ++idn2) {
          atmB = neighbours(atmA + 1,idn2 + 1) - 1;
          if (atmB < 0) {break;}
          if (atmB == idatm) {continue;}
          ++nneighbours;
        }
      }
      if (nneighbours == 0) {tcharge -= 2;}
      else if (nneighbours == 1) {--tcharge;}
    }
    else if (funcgroups[idatm] == "S+") {++tcharge;}
    else if (funcgroups[idatm] == "Na+") {++tcharge;}
    else if (funcgroups[idatm] == "Mg++") {tcharge += 2;}
    else if (funcgroups[idatm] == "K+") {++tcharge;}
    else if (funcgroups[idatm] == "Ca++") {tcharge += 2;}
    else if (funcgroups[idatm] == "Rb+") {++tcharge;}
    else if (funcgroups[idatm] == "Sr++") {tcharge += 2;}
  }
  return tcharge;
}
void Ncoord2Connect(matrix<int> & Connectivity, matrixE & Ncoord, std::vector<int> & neigh, size_t maxcoord = 8, double threshbond = 0.3, double threshsave = 0.5) {
  //function converting a coordination number matrix into connectivity matrix
  int Natoms = Ncoord.rows();
  int connect;
  double maxnc;
  int neighmaxnc;
  Connectivity.resize(Natoms,maxcoord);
  neigh.resize(Natoms);
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    //zero data container
    for (size_t idneigh = 0; idneigh < maxcoord; ++idneigh) {
      Connectivity(idAtm + 1,idneigh + 1) = 0;
    }
    for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
      connect = (Ncoord(idAtm + 1,idBtm + 1) >= threshbond);
      if (connect > 0) {
        Connectivity(idAtm + 1,neigh[idAtm] + 1) = idBtm + 1;
        Connectivity(idBtm + 1,neigh[idBtm] + 1) = idAtm + 1;
      }
      neigh[idAtm] += connect;
      neigh[idBtm] += connect;
    }
  }
  //check the consistency
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    if (neigh[idAtm] == 0) {
      //no neighbours found, this must be related to overly elongated bonds
      maxnc = 0.0;
      neighmaxnc = -1;
      //find the maximum coordination number
      for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
        if (maxnc < Ncoord(idAtm + 1,idBtm + 1)) {
          maxnc = Ncoord(idAtm + 1,idBtm + 1);
          neighmaxnc = idBtm;
        }
      }
      //add the first neighbour
      Connectivity(idAtm + 1,neigh[idAtm] + 1) = neighmaxnc + 1;
      Connectivity(neighmaxnc + 1,neigh[neighmaxnc] + 1) = idAtm + 1;
      ++neigh[idAtm];
      ++neigh[neighmaxnc];
      //now find more potential neighbours
      for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
        if (idBtm == neighmaxnc) {continue;}
        if (Ncoord(idAtm + 1,idBtm + 1)/maxnc > threshsave) {
          Connectivity(idAtm + 1,neigh[idAtm] + 1) = idBtm + 1;
          Connectivity(idBtm + 1,neigh[idBtm] + 1) = idAtm + 1;
          ++neigh[idAtm];
          ++neigh[idBtm];
        }
      }
    }
  }
}

#endif //_Molecule_Analyser_
