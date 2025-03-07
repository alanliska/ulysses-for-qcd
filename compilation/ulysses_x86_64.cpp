#include <iostream>
#include "../src/MolecularDynamics.hpp"
#include "../src/GFN.hpp"
#include "../src/MNDOd.hpp"
#include "../src/math/SolverPackage.hpp"
#include "../src/Gas.hpp"
#include <stdlib.h>

// mndod_pddg: ../src/math/eigen/src/Core/Product.h:96: Eigen::Product<Lhs, Rhs, Option>::Product(const Lhs&, const Rhs&) [with _Lhs = Eigen::Matrix<double, -1, -1>; _Rhs = Eigen::Matrix<double, -1, -1>; int Option = 0; Eigen::Product<Lhs, Rhs, Option>::Lhs = Eigen::Matrix<double, -1, -1>; Eigen::Product<Lhs, Rhs, Option>::Rhs = Eigen::Matrix<double, -1, -1>]: Assertion `lhs.cols() == rhs.rows() && "invalid matrix product" && "if you wanted a coeff-wise or a dot product use the respective explicit functions"' failed.
// to avoid this error, it is necessary to keep the name of the metod consistent with the basis set used
// i.e. not to define the basis e.g. mndo-pddg and the method MNDOPDDG
// but basis has to be mndopddg for the method MNDOPDDG (see MNDO.hpp, MNDOd.hpp)

//basis can be also "slater", "pople", "burns", "gfn2l"
// in the next versions it would be necessary to separate the metod and basis set

int main(int argc, char** argv) {

  // banner - citation, program version
  std::cout << "*** ULYSSES ***" << "\n";
  std::cout << " " << "\n";
  std::cout << "Authors: Filipe Menezes and Grzegorz M. Popowicz" << "\n";
  std::cout << "" << "\n";
  std::cout << "Version for Android (x86_64, pie) powered by EIGEN, RAPIDJSON and XSUM libraries compiled by A. Liska & V. Ruzickova on March 7, 2025." << "\n";
  std::cout << "" << "\n";
  std::cout << "" << "\n";
  std::cout << "References:" << "\n";
  std::cout << "* ULYSSES: Menezes, F.; Popowicz, G. M. ULYSSES: An Efficient and Easy to Use Semiempirical Library for C++. J. Chem. Inf. Model. 2022, , . 10.1021/acs.jcim.2c00757" << "\n";
  std::cout << "* MNDO: a) Dewar, M. J. S.; Thiel, W. Ground states of molecules. 38. The MNDO method. Approximations and parameters. J. Am. Chem. Soc. 1977, 99, 4899. 10.1021/ja00457a004; b) Dewar, M. J. S.; Thiel, W. A semiempirical model for the two-center repulsion integrals in the NDDO approximation. Theor. Chim. Acta (Berl.) 1977, 46, 89. 10.1007/BF00548085" << "\n";
  std::cout << "* AM1: Dewar, M. J. S.; Zoebisch, E. G.; Healy, E. F.; Stewart, J. J. P. Development and use of quantum mechanical molecular models. 76. AM1: a new general purpose quantum mechanical molecular model. J. Am. Chem. Soc. 1985, 107, 3902. 10.1021/ja00299a024" << "\n";
  std::cout << "* PM3: a) Stewart, J. J. P. Optimization of Parameters for Semi-Empirical Methods I. Method. J. Comput. Chem. 1989, 10, 209. 10.1002/jcc.540100208; b) Stewart, J. J. P. Optimization of parameters for semiempirical methods. III Extension of PM3 to Be, Mg, Zn, Ga, Ge, As, Se, Cd, In, Sn, Sb, Te, Hg, Tl, Pb, and Bi. J. Comput. Chem. 1991, 12, 320. 10.1002/jcc.540120306" << "\n";
  std::cout << "* PM3-PDDG: a) Repasky, M. P.; Chandrasekhar, J.; Jorgensen, W. L. PDDG/PM3 and PDDG/MNDO: improved semiempirical methods. J. Comput. Chem. 2002, 23, 1601. 10.1002/jcc.10162; b) Tubert-Brohman, I.; Guimaraes, C. R. W.; Repasky, M. P.; Jorgensen, W. L. Extension of the PDDG/PM3 and PDDG/MNDO semiempirical molecular orbital methods to the halogens. J. Comput. Chem. 2004, 25, 138. 10.1002/jcc.10356; c) Tubert-Brohman, I.; Guimaraes, C. R. W.; Jorgensen, W. L. Extension of the PDDG/PM3 Semiempirical Molecular Orbital Method to Sulfur, Silicon, and Phosphorus. J. Chem. Theory Comput. 2005, 1, 817. 10.1021/ct0500287" << "\n";
  std::cout << "* MNDO-PDDG: a) Repasky, M. P.; Chandrasekhar, J.; Jorgensen, W. L. PDDG/PM3 and PDDG/MNDO: improved semiempirical methods. J. Comput. Chem. 2002, 23, 1601. 10.1002/jcc.10162; b) Tubert-Brohman, I.; Guimaraes, C. R. W.; Repasky, M. P.; Jorgensen, W. L. Extension of the PDDG/PM3 and PDDG/MNDO semiempirical molecular orbital methods to the halogens. J. Comput. Chem. 2004, 25, 138. 10.1002/jcc.10356" << "\n";
  std::cout << "* PM3-BP: Giese, T. J.; Sherer, E. C.; Cramer, C. J.; York, D. M. A Semiempirical Quantum Model for Hydrogen-Bonded Nucleic Acid Base Pairs. J. Chem. Theory Comput. 2005, 1, 1275. 10.1021/ct050102l" << "\n";
  std::cout << "* RM1: Rocha, G. B.; Freire, R. O.; Simas, A. M.; Stewart, J. J. P. RM1: a reparameterization of AM1 for H, C, N, O, P, S, F, Cl, Br, and I. J. Comput. Chem. 2006, 27, 1101. 10.1002/jcc.20425" << "\n";
  std::cout << "* MNDOd: a) Thiel, W.; Voityuk, A. A. Extension of MNDO to d Orbitals: Parameters and Results for the Second-Row Elements and for the Zinc Group. J. Phys. Chem. 1996, 100, 616-626. 10.1021/jp952148o; b) Thiel, W.; Voityuk, A. A. Extension of the MNDO formalism to d orbitals: Integral approximations and preliminary numerical results. Theor. Chim. Acta 1992, 81, 391. 10.1007/BF01134863" << "\n";
  std::cout << "* PM6: Stewart, J. J. P. Optimization of parameters for semiempirical methods V: Modification of NDDO approximations. J. Mol. Model. 2007, 13, 1173. 10.1007/s00894-007-0233-4" << "\n";
  std::cout << "* GFN2-xTB: Bannwarth, C.; Ehlert, S.; Grimme, S. GFN2-xTB - An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions. J. Chem. Theory Comput. 2019, 15, 1652-1671. 10.1021/acs.jctc.8b01176" << "\n";
  std::cout << "Non-Bonded Corrections:" << "\n";
  std::cout << "* D3: Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H. A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu. J. Chem. Phys. 2010, 132, 154104. 10.1063/1.3382344" << "\n";
  std::cout << "* D4: Caldeweyher, E.; Ehlert, S.; Hansen, A.; Neugebauer, H.; Spicher, S.; Bannwarth, C.; Grimme, S. A generally applicable atomic-charge dependent London dispersion correction. J. Chem. Phys. 2019, 150, 154122. 10.1063/1.5090222" << "\n";
  std::cout << "* H4: Rezac, J.; Hobza, P. Advanced Corrections of Hydrogen Bonding. J. Chem. Theory Comput. 2012, 8, 141-151. 10.1021/ct200751e" << "\n";
  std::cout << "* X: Rezac, J.; Hobza, P. A halogen-bonding correction for the semiempirical PM6 method. Chem. Phys. Lett. 2011, 506, 286-289. 10.1016/j.cplett.2011.03.009" << "\n";
  std::cout << "* H+: Kromann, J. C.; Christensen, A. S.; Steinmann, C.; Korth, M.; Jensen, J. H. A third-generation dispersion and third-generation hydrogen bonding corrected PM6 method: PM6-D3H+. Peer J. 2014, 2, e449. 10.7717/peerj.449" << "\n";
  std::cout << "" << "\n";

  std::ifstream is10 ("runtype.tmp");
  std::string runtype;
  is10 >> runtype;
  is10.close();
  
  
  
  
  
  if (runtype == "md") {
  
  
  
  
  
  //Test MD
  std::cout << "--- Running molecular dynamics simulation ---" << std::endl;

  //argument 1 to program is the geometry in xyz format
  //argument 2 to program is the charge
  //argument 3 to program is the temperature in Kelvin
  //argument 4 to program is the total simulation time (ps)
  //argument 5 to program is the time time (ps)

  std::cout << "running " << argv[1] << "\n";

  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  //char *s;
  //double Temperature = strtod(argv[3],&s);
  std::ifstream is13 ("temperature.tmp");
  double Temperature;
  is13 >> Temperature;
  is13.close();
  //char *q;
  //double tmax = strtod(argv[4],&q);
  std::ifstream is11 ("simulation_time.tmp");
  double tmax;
  is11 >> tmax;
  is11.close();
  //char *r;
  //double tstep = strtod(argv[5],&r);
  std::ifstream is12 ("time_step.tmp");
  double tstep;
  is12 >> tstep;
  is12.close();
  std::ifstream is4 ("basis_set.tmp");
  std::string basis_set;
  is4 >> basis_set;
  is4.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  
  // frequency (in ps) at which geometries are dumped
  std::ifstream is36 ("t_opt.tmp");
  double topt;
  is36 >> topt;
  is36.close();
  
  std::ifstream is34 ("md_do_equilibration.tmp");
  int do_equ;
  is34 >> do_equ;
  is34.close();
  
  std::ifstream is35 ("md_opt_geometry.tmp");
  int opt_geometry;
  is35 >> opt_geometry;
  is35.close();
  bool optgeometry = (opt_geometry > 0);
  
  //Molecule mol(argv[1],charge,1);
  Molecule mol(argv[1],charge,mult,point_group);
  //BSet basis(mol,"gfn2");
  BSet basis(mol,basis_set);
  GFN2 electron(basis,mol);
  
  std::ifstream is5 ("solvent.tmp");
  std::string solvent;
  is5 >> solvent;
  is5.close();
  // if (strcmp(solvent, "without_solvent") != 0) {
  if (solvent != "without_solvent") {
  electron.setSolvent(solvent);
  }
  
  std::ifstream is40 ("md_shake.tmp");
  int SHAKE;
  is40 >> SHAKE;
  is40.close();
  
  std::ifstream is48 ("md_print_freq.tmp");
  int printfreq;
  is48 >> printfreq;
  is48.close();
  // the folder must exist before, ULYSSES program does not create it
  std::ifstream is39 ("md_geometry_location.tmp");
  std::string geom_loc;
  is39 >> geom_loc;
  is39.close();
  
  ///
  std::cout << "Debug: 1" << std::endl;
  ///
  
  Dynamics MDobj(mol,tmax,do_equ,optgeometry,tstep,topt);
  MDobj.setGeometryFile(geom_loc);
  MDobj.ApplyConstraints(false,SHAKE);
  MDobj.setIntegration("LeapFrog");
  
  std::ifstream is37 ("md_trajectory_location.tmp");
  std::string traj_loc;
  is37 >> traj_loc;
  is37.close();
  
  std::ifstream is38 ("md_equilibration_location.tmp");
  std::string equi_loc;
  is38 >> equi_loc;
  is38.close();
  
  ///
  std::cout << "Debug: 2" << std::endl;
  ///
  
  //where to save files
  MDobj.setTrajectoryFile(traj_loc);
  MDobj.setEquilibrationFile(equi_loc);
  
  ///
  std::cout << "Debug: 3" << std::endl;
  ///
  
  //metadynamics?
  std::ifstream is41 ("md_mtd.tmp");
  int dometadyn;
  is41 >> dometadyn;
  is41.close();
  std::ifstream is42 ("md_numb_struct.tmp");
  int numbstruct;
  is42 >> numbstruct;
  is42.close();
  std::ifstream is43 ("md_kappa.tmp");
  double kappa;
  is43 >> kappa;
  is43.close();
  std::ifstream is44 ("md_alpha.tmp");
  double alpha;
  is44 >> alpha;
  is44.close();
  std::ifstream is45 ("md_mtdcollect.tmp");
  int mtdcollect;
  is45 >> mtdcollect;
  is45.close();
  // error - must have the .xyz extension
  // input the file name, not its content!
  // std::ifstream is46 ("md_bias_mol.tmp");
  // std::string bias_mol;
  // is46 >> bias_mol;
  // is46.close();
  std::ifstream is47 ("md_drift_thresh.tmp");
  double driftthresh;
  is47 >> driftthresh;
  is47.close();
  
  if (dometadyn > 0) {
    MDobj.setMetaDynamics(true,numbstruct,kappa,alpha,mtdcollect,0.03);
    std::vector<int> mtdrestr = getNumbersFromFile("md_restraints.tmp");
    for (size_t idatm = 0; idatm < mtdrestr.size(); ++idatm) {
      MDobj.addMetaAtom(mtdrestr[idatm]);
    }
    Molecule newmol("md_bias_mol.xyz",charge,mult,point_group);
    if ((newmol.Natoms() > 0)&&(kappa < 0.0)) {
      std::vector<matrixE> metaset;
      metaset.resize(numbstruct);
      matrixE geom = newmol.Geometry();
      for (size_t idmt = 0; idmt < numbstruct; ++idmt) {
        metaset[idmt] = geom;
      }
      MDobj.setMetaSet(metaset);
    }
  }
  // at this point, the MD folder must exist, otherwise an error occurs
  MDobj.runMD(electron,Temperature,driftthresh,printfreq);
  
  } else if (runtype == "nddo") {
  
  std::ifstream is4 ("basis_set.tmp");
  std::string basis_set;
  is4 >> basis_set;
  is4.close();
  
  if (basis_set == "mndo") {
  //program to run PM6-D3H4X calculations
  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "--- Running MNDO calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";

  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  //char *q;
  //double lshift = strtod(argv[3],&q);
  std::ifstream is8 ("level_shift.tmp");
  double lshift;
  is8 >> lshift;
  is8.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  //Molecule Mol1(argv[1],charge,1,"C1");
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //BSet basis(Mol1,"pm6");
  BSet basis(Mol1,basis_set);
  
  // this block causes the error:
  //...terminate called after throwing an instance of 'char const*'
  //std::string radical;
  //radical == "0";
  //if (mult == 1) {
  //radical == "0";
  //} else {
  //radical == "UHF";
  //}
  
  std::ifstream is9 ("correction.tmp");
  std::string correction;
  is9 >> correction;
  is9.close();
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //PM6 electron(basis,Mol1,"0","D3H4X");
  MNDO electron(basis,Mol1,radical,correction);
  
  electron.setEpsilonS(lshift);
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  //electron.Calculate(0);
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
  std::cout << "Geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);

  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  std::ifstream is19 ("charge_types.tmp");
  std::string charge_types;
  is19 >> charge_types;
  is19.close();
  
    //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful because most CMx charges only available for AM1 and PM3
  //std::vector<double> charges = electron.getCharges("Mulliken");
  std::cout << std::endl;
  if (charge_types == "Mulliken"){
  std::vector<double> charges_Mulliken = electron.getCharges("Mulliken");
  std::cout << "Mulliken charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Mulliken.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Mulliken[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Loewdin"){
  std::vector<double> charges_Loewdin = electron.getCharges("Loewdin");
  std::cout << "Loewdin charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Loewdin.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Loewdin[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Goedecker"){
  std::vector<double> charges_Goedecker = electron.getCharges("Goedecker");
  std::cout << "Goedecker charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Goedecker.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Goedecker[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM1"){
  std::vector<double> charges_CM1 = electron.getCharges("CM1");
  std::cout << "CM1 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM1.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM1[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM2"){
  std::vector<double> charges_CM2 = electron.getCharges("CM2");
  std::cout << "CM2 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM2.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM2[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM3"){
  std::vector<double> charges_CM3 = electron.getCharges("CM3");
  std::cout << "CM3 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM3.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM3[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM5"){
  std::vector<double> charges_CM5 = electron.getCharges("CM5");
  std::cout << "CM5 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM5.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM5[idAtm] << std::endl;
  }
  std::cout << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  std::ifstream is24 ("calc_density.tmp");
  int calc_density;
  is24 >> calc_density;
  is24.close();
  
  std::ifstream is25 ("density_file.tmp");
  std::string density_file;
  is25 >> density_file;
  is25.close();
  //get the density?
  if (calc_density > 0) {
    electron.ElectronicDensity(density_file);
  }
  } else if (basis_set == "am1") {
  //program to run PM6-D3H4X calculations
  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "--- Running AM1 calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";

  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  //char *q;
  //double lshift = strtod(argv[3],&q);
  std::ifstream is8 ("level_shift.tmp");
  double lshift;
  is8 >> lshift;
  is8.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  //Molecule Mol1(argv[1],charge,1,"C1");
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //BSet basis(Mol1,"pm6");
  BSet basis(Mol1,basis_set);
  
  // this block causes the error:
  //...terminate called after throwing an instance of 'char const*'
  //std::string radical;
  //radical == "0";
  //if (mult == 1) {
  //radical == "0";
  //} else {
  //radical == "UHF";
  //}
  
  std::ifstream is9 ("correction.tmp");
  std::string correction;
  is9 >> correction;
  is9.close();
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //PM6 electron(basis,Mol1,"0","D3H4X");
  AM1 electron(basis,Mol1,radical,correction);
  
  electron.setEpsilonS(lshift);
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  //electron.Calculate(0);
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
  std::cout << "Geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);

  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  std::ifstream is19 ("charge_types.tmp");
  std::string charge_types;
  is19 >> charge_types;
  is19.close();
  
    //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful because most CMx charges only available for AM1 and PM3
  //std::vector<double> charges = electron.getCharges("Mulliken");
  std::cout << std::endl;
  if (charge_types == "Mulliken"){
  std::vector<double> charges_Mulliken = electron.getCharges("Mulliken");
  std::cout << "Mulliken charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Mulliken.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Mulliken[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Loewdin"){
  std::vector<double> charges_Loewdin = electron.getCharges("Loewdin");
  std::cout << "Loewdin charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Loewdin.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Loewdin[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Goedecker"){
  std::vector<double> charges_Goedecker = electron.getCharges("Goedecker");
  std::cout << "Goedecker charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Goedecker.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Goedecker[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM1"){
  std::vector<double> charges_CM1 = electron.getCharges("CM1");
  std::cout << "CM1 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM1.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM1[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM2"){
  std::vector<double> charges_CM2 = electron.getCharges("CM2");
  std::cout << "CM2 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM2.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM2[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM3"){
  std::vector<double> charges_CM3 = electron.getCharges("CM3");
  std::cout << "CM3 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM3.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM3[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM5"){
  std::vector<double> charges_CM5 = electron.getCharges("CM5");
  std::cout << "CM5 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM5.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM5[idAtm] << std::endl;
  }
  std::cout << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  std::ifstream is24 ("calc_density.tmp");
  int calc_density;
  is24 >> calc_density;
  is24.close();
  
  std::ifstream is25 ("density_file.tmp");
  std::string density_file;
  is25 >> density_file;
  is25.close();
  //get the density?
  if (calc_density > 0) {
    electron.ElectronicDensity(density_file);
  }
  } else if (basis_set == "pm3") {
  //program to run PM6-D3H4X calculations
  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "--- Running PM3 calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";

  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  //char *q;
  //double lshift = strtod(argv[3],&q);
  std::ifstream is8 ("level_shift.tmp");
  double lshift;
  is8 >> lshift;
  is8.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  //Molecule Mol1(argv[1],charge,1,"C1");
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //BSet basis(Mol1,"pm6");
  BSet basis(Mol1,basis_set);
  
  // this block causes the error:
  //...terminate called after throwing an instance of 'char const*'
  //std::string radical;
  //radical == "0";
  //if (mult == 1) {
  //radical == "0";
  //} else {
  //radical == "UHF";
  //}
  
  std::ifstream is9 ("correction.tmp");
  std::string correction;
  is9 >> correction;
  is9.close();
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //PM6 electron(basis,Mol1,"0","D3H4X");
  PM3 electron(basis,Mol1,radical,correction);
  
  electron.setEpsilonS(lshift);
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  //electron.Calculate(0);
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
  std::cout << "Geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);

  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  std::ifstream is19 ("charge_types.tmp");
  std::string charge_types;
  is19 >> charge_types;
  is19.close();
  
    //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful because most CMx charges only available for AM1 and PM3
  //std::vector<double> charges = electron.getCharges("Mulliken");
  std::cout << std::endl;
  if (charge_types == "Mulliken"){
  std::vector<double> charges_Mulliken = electron.getCharges("Mulliken");
  std::cout << "Mulliken charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Mulliken.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Mulliken[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Loewdin"){
  std::vector<double> charges_Loewdin = electron.getCharges("Loewdin");
  std::cout << "Loewdin charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Loewdin.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Loewdin[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Goedecker"){
  std::vector<double> charges_Goedecker = electron.getCharges("Goedecker");
  std::cout << "Goedecker charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Goedecker.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Goedecker[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM1"){
  std::vector<double> charges_CM1 = electron.getCharges("CM1");
  std::cout << "CM1 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM1.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM1[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM2"){
  std::vector<double> charges_CM2 = electron.getCharges("CM2");
  std::cout << "CM2 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM2.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM2[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM3"){
  std::vector<double> charges_CM3 = electron.getCharges("CM3");
  std::cout << "CM3 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM3.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM3[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM5"){
  std::vector<double> charges_CM5 = electron.getCharges("CM5");
  std::cout << "CM5 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM5.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM5[idAtm] << std::endl;
  }
  std::cout << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  std::ifstream is24 ("calc_density.tmp");
  int calc_density;
  is24 >> calc_density;
  is24.close();
  
  std::ifstream is25 ("density_file.tmp");
  std::string density_file;
  is25 >> density_file;
  is25.close();
  //get the density?
  if (calc_density > 0) {
    electron.ElectronicDensity(density_file);
  }
  } else if (basis_set == "pm3pddg") {
  //program to run PM6-D3H4X calculations
  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "--- Running PM3-PDDG calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";

  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  //char *q;
  //double lshift = strtod(argv[3],&q);
  std::ifstream is8 ("level_shift.tmp");
  double lshift;
  is8 >> lshift;
  is8.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  //Molecule Mol1(argv[1],charge,1,"C1");
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //BSet basis(Mol1,"pm6");
  BSet basis(Mol1,basis_set);
  
  // this block causes the error:
  //...terminate called after throwing an instance of 'char const*'
  //std::string radical;
  //radical == "0";
  //if (mult == 1) {
  //radical == "0";
  //} else {
  //radical == "UHF";
  //}
  
  std::ifstream is9 ("correction.tmp");
  std::string correction;
  is9 >> correction;
  is9.close();
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //PM6 electron(basis,Mol1,"0","D3H4X");
  PM3PDDG electron(basis,Mol1,radical,correction);
  
  electron.setEpsilonS(lshift);
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  //electron.Calculate(0);
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
  std::cout << "Geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);

  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  std::ifstream is19 ("charge_types.tmp");
  std::string charge_types;
  is19 >> charge_types;
  is19.close();
  
    //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful because most CMx charges only available for AM1 and PM3
  //std::vector<double> charges = electron.getCharges("Mulliken");
  std::cout << std::endl;
  if (charge_types == "Mulliken"){
  std::vector<double> charges_Mulliken = electron.getCharges("Mulliken");
  std::cout << "Mulliken charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Mulliken.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Mulliken[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Loewdin"){
  std::vector<double> charges_Loewdin = electron.getCharges("Loewdin");
  std::cout << "Loewdin charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Loewdin.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Loewdin[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Goedecker"){
  std::vector<double> charges_Goedecker = electron.getCharges("Goedecker");
  std::cout << "Goedecker charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Goedecker.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Goedecker[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM1"){
  std::vector<double> charges_CM1 = electron.getCharges("CM1");
  std::cout << "CM1 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM1.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM1[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM2"){
  std::vector<double> charges_CM2 = electron.getCharges("CM2");
  std::cout << "CM2 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM2.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM2[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM3"){
  std::vector<double> charges_CM3 = electron.getCharges("CM3");
  std::cout << "CM3 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM3.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM3[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM5"){
  std::vector<double> charges_CM5 = electron.getCharges("CM5");
  std::cout << "CM5 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM5.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM5[idAtm] << std::endl;
  }
  std::cout << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  std::ifstream is24 ("calc_density.tmp");
  int calc_density;
  is24 >> calc_density;
  is24.close();
  
  std::ifstream is25 ("density_file.tmp");
  std::string density_file;
  is25 >> density_file;
  is25.close();
  //get the density?
  if (calc_density > 0) {
    electron.ElectronicDensity(density_file);
  }
  } else if (basis_set == "pm3bp") {
  //program to run PM6-D3H4X calculations
  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "--- Running PM3-BP calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";

  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  //char *q;
  //double lshift = strtod(argv[3],&q);
  std::ifstream is8 ("level_shift.tmp");
  double lshift;
  is8 >> lshift;
  is8.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  //Molecule Mol1(argv[1],charge,1,"C1");
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //BSet basis(Mol1,"pm6");
  BSet basis(Mol1,basis_set);
  
  // this block causes the error:
  //...terminate called after throwing an instance of 'char const*'
  //std::string radical;
  //radical == "0";
  //if (mult == 1) {
  //radical == "0";
  //} else {
  //radical == "UHF";
  //}
  
  std::ifstream is9 ("correction.tmp");
  std::string correction;
  is9 >> correction;
  is9.close();
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //PM6 electron(basis,Mol1,"0","D3H4X");
  PM3BP electron(basis,Mol1,radical,correction);
  
  electron.setEpsilonS(lshift);
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  //electron.Calculate(0);
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
  std::cout << "Geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);

  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  std::ifstream is19 ("charge_types.tmp");
  std::string charge_types;
  is19 >> charge_types;
  is19.close();
  
    //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful because most CMx charges only available for AM1 and PM3
  //std::vector<double> charges = electron.getCharges("Mulliken");
  std::cout << std::endl;
  if (charge_types == "Mulliken"){
  std::vector<double> charges_Mulliken = electron.getCharges("Mulliken");
  std::cout << "Mulliken charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Mulliken.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Mulliken[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Loewdin"){
  std::vector<double> charges_Loewdin = electron.getCharges("Loewdin");
  std::cout << "Loewdin charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Loewdin.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Loewdin[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Goedecker"){
  std::vector<double> charges_Goedecker = electron.getCharges("Goedecker");
  std::cout << "Goedecker charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Goedecker.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Goedecker[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM1"){
  std::vector<double> charges_CM1 = electron.getCharges("CM1");
  std::cout << "CM1 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM1.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM1[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM2"){
  std::vector<double> charges_CM2 = electron.getCharges("CM2");
  std::cout << "CM2 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM2.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM2[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM3"){
  std::vector<double> charges_CM3 = electron.getCharges("CM3");
  std::cout << "CM3 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM3.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM3[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM5"){
  std::vector<double> charges_CM5 = electron.getCharges("CM5");
  std::cout << "CM5 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM5.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM5[idAtm] << std::endl;
  }
  std::cout << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  std::ifstream is24 ("calc_density.tmp");
  int calc_density;
  is24 >> calc_density;
  is24.close();
  
  std::ifstream is25 ("density_file.tmp");
  std::string density_file;
  is25 >> density_file;
  is25.close();
  //get the density?
  if (calc_density > 0) {
    electron.ElectronicDensity(density_file);
  }
  } else if (basis_set == "rm1") {
  //program to run PM6-D3H4X calculations
  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "--- Running RM1 calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";
  
  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  //char *q;
  //double lshift = strtod(argv[3],&q);
  std::ifstream is8 ("level_shift.tmp");
  double lshift;
  is8 >> lshift;
  is8.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  //Molecule Mol1(argv[1],charge,1,"C1");
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //BSet basis(Mol1,"pm6");
  BSet basis(Mol1,basis_set);
  
  // this block causes the error:
  //...terminate called after throwing an instance of 'char const*'
  //std::string radical;
  //radical == "0";
  //if (mult == 1) {
  //radical == "0";
  //} else {
  //radical == "UHF";
  //}
  
  std::ifstream is9 ("correction.tmp");
  std::string correction;
  is9 >> correction;
  is9.close();
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //PM6 electron(basis,Mol1,"0","D3H4X");
  RM1 electron(basis,Mol1,radical,correction);
  
  electron.setEpsilonS(lshift);
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  //electron.Calculate(0);
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
  std::cout << "Geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);

  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  std::ifstream is19 ("charge_types.tmp");
  std::string charge_types;
  is19 >> charge_types;
  is19.close();
  
    //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful because most CMx charges only available for AM1 and PM3
  //std::vector<double> charges = electron.getCharges("Mulliken");
  std::cout << std::endl;
  if (charge_types == "Mulliken"){
  std::vector<double> charges_Mulliken = electron.getCharges("Mulliken");
  std::cout << "Mulliken charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Mulliken.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Mulliken[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Loewdin"){
  std::vector<double> charges_Loewdin = electron.getCharges("Loewdin");
  std::cout << "Loewdin charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Loewdin.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Loewdin[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Goedecker"){
  std::vector<double> charges_Goedecker = electron.getCharges("Goedecker");
  std::cout << "Goedecker charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Goedecker.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Goedecker[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM1"){
  std::vector<double> charges_CM1 = electron.getCharges("CM1");
  std::cout << "CM1 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM1.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM1[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM2"){
  std::vector<double> charges_CM2 = electron.getCharges("CM2");
  std::cout << "CM2 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM2.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM2[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM3"){
  std::vector<double> charges_CM3 = electron.getCharges("CM3");
  std::cout << "CM3 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM3.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM3[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM5"){
  std::vector<double> charges_CM5 = electron.getCharges("CM5");
  std::cout << "CM5 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM5.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM5[idAtm] << std::endl;
  }
  std::cout << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  std::ifstream is24 ("calc_density.tmp");
  int calc_density;
  is24 >> calc_density;
  is24.close();
  
  std::ifstream is25 ("density_file.tmp");
  std::string density_file;
  is25 >> density_file;
  is25.close();
  //get the density?
  if (calc_density > 0) {
    electron.ElectronicDensity(density_file);
  }
  } else if (basis_set == "mndod") {
  //program to run PM6-D3H4X calculations
  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "--- Running MNDOd calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";

  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  //char *q;
  //double lshift = strtod(argv[3],&q);
  std::ifstream is8 ("level_shift.tmp");
  double lshift;
  is8 >> lshift;
  is8.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  //Molecule Mol1(argv[1],charge,1,"C1");
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //BSet basis(Mol1,"pm6");
  BSet basis(Mol1,basis_set);
  
  // this block causes the error:
  //...terminate called after throwing an instance of 'char const*'
  //std::string radical;
  //radical == "0";
  //if (mult == 1) {
  //radical == "0";
  //} else {
  //radical == "UHF";
  //}
  
  std::ifstream is9 ("correction.tmp");
  std::string correction;
  is9 >> correction;
  is9.close();
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //PM6 electron(basis,Mol1,"0","D3H4X");
  MNDOd electron(basis,Mol1,radical,correction);
  
  electron.setEpsilonS(lshift);
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  //electron.Calculate(0);
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
  std::cout << "Geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);

  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  std::ifstream is19 ("charge_types.tmp");
  std::string charge_types;
  is19 >> charge_types;
  is19.close();
  
    //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful because most CMx charges only available for AM1 and PM3
  //std::vector<double> charges = electron.getCharges("Mulliken");
  std::cout << std::endl;
  if (charge_types == "Mulliken"){
  std::vector<double> charges_Mulliken = electron.getCharges("Mulliken");
  std::cout << "Mulliken charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Mulliken.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Mulliken[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Loewdin"){
  std::vector<double> charges_Loewdin = electron.getCharges("Loewdin");
  std::cout << "Loewdin charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Loewdin.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Loewdin[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Goedecker"){
  std::vector<double> charges_Goedecker = electron.getCharges("Goedecker");
  std::cout << "Goedecker charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Goedecker.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Goedecker[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM1"){
  std::vector<double> charges_CM1 = electron.getCharges("CM1");
  std::cout << "CM1 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM1.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM1[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM2"){
  std::vector<double> charges_CM2 = electron.getCharges("CM2");
  std::cout << "CM2 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM2.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM2[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM3"){
  std::vector<double> charges_CM3 = electron.getCharges("CM3");
  std::cout << "CM3 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM3.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM3[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM5"){
  std::vector<double> charges_CM5 = electron.getCharges("CM5");
  std::cout << "CM5 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM5.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM5[idAtm] << std::endl;
  }
  std::cout << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  std::ifstream is24 ("calc_density.tmp");
  int calc_density;
  is24 >> calc_density;
  is24.close();
  
  std::ifstream is25 ("density_file.tmp");
  std::string density_file;
  is25 >> density_file;
  is25.close();
  //get the density?
  if (calc_density > 0) {
    electron.ElectronicDensity(density_file);
  }
  } else if (basis_set == "mndopddg") {
  //program to run PM6-D3H4X calculations
  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "--- Running MNDO-PDDG calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";

  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  //char *q;
  //double lshift = strtod(argv[3],&q);
  std::ifstream is8 ("level_shift.tmp");
  double lshift;
  is8 >> lshift;
  is8.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  //Molecule Mol1(argv[1],charge,1,"C1");
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //BSet basis(Mol1,"pm6");
  BSet basis(Mol1,basis_set);
  
  // this block causes the error:
  //...terminate called after throwing an instance of 'char const*'
  //std::string radical;
  //radical == "0";
  //if (mult == 1) {
  //radical == "0";
  //} else {
  //radical == "UHF";
  //}
  
  std::ifstream is9 ("correction.tmp");
  std::string correction;
  is9 >> correction;
  is9.close();
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //PM6 electron(basis,Mol1,"0","D3H4X");
  MNDOd electron(basis,Mol1,radical,correction);
  
  electron.setEpsilonS(lshift);
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  //electron.Calculate(0);
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
  std::cout << "Geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);

  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  std::ifstream is19 ("charge_types.tmp");
  std::string charge_types;
  is19 >> charge_types;
  is19.close();
  
    //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful because most CMx charges only available for AM1 and PM3
  //std::vector<double> charges = electron.getCharges("Mulliken");
  std::cout << std::endl;
  if (charge_types == "Mulliken"){
  std::vector<double> charges_Mulliken = electron.getCharges("Mulliken");
  std::cout << "Mulliken charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Mulliken.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Mulliken[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Loewdin"){
  std::vector<double> charges_Loewdin = electron.getCharges("Loewdin");
  std::cout << "Loewdin charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Loewdin.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Loewdin[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Goedecker"){
  std::vector<double> charges_Goedecker = electron.getCharges("Goedecker");
  std::cout << "Goedecker charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Goedecker.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Goedecker[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM1"){
  std::vector<double> charges_CM1 = electron.getCharges("CM1");
  std::cout << "CM1 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM1.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM1[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM2"){
  std::vector<double> charges_CM2 = electron.getCharges("CM2");
  std::cout << "CM2 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM2.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM2[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM3"){
  std::vector<double> charges_CM3 = electron.getCharges("CM3");
  std::cout << "CM3 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM3.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM3[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM5"){
  std::vector<double> charges_CM5 = electron.getCharges("CM5");
  std::cout << "CM5 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM5.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM5[idAtm] << std::endl;
  }
  std::cout << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  std::ifstream is24 ("calc_density.tmp");
  int calc_density;
  is24 >> calc_density;
  is24.close();
  
  std::ifstream is25 ("density_file.tmp");
  std::string density_file;
  is25 >> density_file;
  is25.close();
  //get the density?
  if (calc_density > 0) {
    electron.ElectronicDensity(density_file);
  }
  } else if (basis_set == "pm6") {
  //program to run PM6-D3H4X calculations
  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "--- Running PM6 calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";

  //char *p;
  //int charge = strtol(argv[2],&p,10);
  std::ifstream is1 ("charge.tmp");
  int charge;
  is1 >> charge;
  is1.close();
  //char *q;
  //double lshift = strtod(argv[3],&q);
  std::ifstream is8 ("level_shift.tmp");
  double lshift;
  is8 >> lshift;
  is8.close();
  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  //Molecule Mol1(argv[1],charge,1,"C1");
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //BSet basis(Mol1,"pm6");
  BSet basis(Mol1,basis_set);
  
  // this block causes the error:
  //...terminate called after throwing an instance of 'char const*'
  //std::string radical;
  //radical == "0";
  //if (mult == 1) {
  //radical == "0";
  //} else {
  //radical == "UHF";
  //}
  
  std::ifstream is9 ("correction.tmp");
  std::string correction;
  is9 >> correction;
  is9.close();
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //PM6 electron(basis,Mol1,"0","D3H4X");
  PM6 electron(basis,Mol1,radical,correction);
  
  electron.setEpsilonS(lshift);
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  //electron.Calculate(0);
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
  std::cout << "Geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  //SolverOpt(electron,solve,4,0,5e-6,1e-3);
  SolverOpt(electron,solve,hessian_update,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);

  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  std::ifstream is19 ("charge_types.tmp");
  std::string charge_types;
  is19 >> charge_types;
  is19.close();
  
    //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful because most CMx charges only available for AM1 and PM3
  //std::vector<double> charges = electron.getCharges("Mulliken");
  std::cout << std::endl;
  if (charge_types == "Mulliken"){
  std::vector<double> charges_Mulliken = electron.getCharges("Mulliken");
  std::cout << "Mulliken charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Mulliken.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Mulliken[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Loewdin"){
  std::vector<double> charges_Loewdin = electron.getCharges("Loewdin");
  std::cout << "Loewdin charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Loewdin.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Loewdin[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "Goedecker"){
  std::vector<double> charges_Goedecker = electron.getCharges("Goedecker");
  std::cout << "Goedecker charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_Goedecker.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_Goedecker[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM1"){
  std::vector<double> charges_CM1 = electron.getCharges("CM1");
  std::cout << "CM1 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM1.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM1[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM2"){
  std::vector<double> charges_CM2 = electron.getCharges("CM2");
  std::cout << "CM2 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM2.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM2[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM3"){
  std::vector<double> charges_CM3 = electron.getCharges("CM3");
  std::cout << "CM3 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM3.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM3[idAtm] << std::endl;
  }
  std::cout << std::endl;
  } else if (charge_types == "CM5"){
  std::vector<double> charges_CM5 = electron.getCharges("CM5");
  std::cout << "CM5 charges:" << std::endl;
  for (size_t idAtm = 0; idAtm < charges_CM5.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges_CM5[idAtm] << std::endl;
  }
  std::cout << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  } 
  
  
  
  
  
  
  
  
  
  
  } else if (runtype == "xtb") {
  
  
  
  
  
  
  std::cout << "--- Running XTB-GFN2 calculation ---" << std::endl;
  std::cout << "loading " << argv[1] << "\n";
  
  std::ifstream is13 ("temperature.tmp");
  double temperature;
  is13 >> temperature;
  is13.close();
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  //argument 2 to program is the charge
  int charge;
  std::ifstream is1 ("charge.tmp");
  is1 >> charge;
  is1.close();
  //char *p;
  //int charge = strtol(charge_str,&p,10);
  //int charge = 0;
  //int charge = strtol("0",&p,10);
  std::cout << "charge = " << charge << std::endl;
  
  //argument 3 to program is the electronic temperature (default = 300 K)
  double Telec;
  std::ifstream is2 ("electronic_temperature.tmp");
  is2 >> Telec;
  is2.close();
  //char *q;
  //double Telec = strtod(Telec_str,&q);
  //double Telec = 300.0;
  //double Telec = strtod("300.0",&q);
  std::cout << "T electron = " << Telec << std::endl;

  // define the molecule symmetry;
  // DINFh, CINFv, C2h, C3h, C4h, C5h, C6h, C2v, C3v, C4v, C5v, C6v, C7v, C8v, C1, C2, C3, C4, C5, C6, C7, C8, Cs, Ci, D2h, D3h, D4h, D5h, D6h, D7h, D8h, D2d, D3d, D4d, D5d, D6d, D7d, D8d, D2, D3, D4, D5, D6, D7, D8, T, Th, Td, O, Oh, I, Ih, S2, S4, S6, S8, S10, S12
  std::ifstream is7 ("point_group.tmp");
  std::string point_group;
  is7 >> point_group;
  is7.close();
  
  //allocation of molecule with certain charge and multiplicity 2S+1 = 1
  // Molecule Mol1(argv[1],charge,1);
  std::ifstream is3 ("multiplicity.tmp");
  int mult;
  is3 >> mult;
  is3.close();
  Molecule Mol1(argv[1],charge,mult,point_group);
  
  //allocation of basis set object
  // BSet basis(Mol1,"gfn2");
  // basis_set = gfn2, pm6, am1, pm3, mndo, mndod, rm1, mndo-pddg, pm3-pddg, pm3-bp
  std::ifstream is4 ("basis_set.tmp");
  std::string basis_set;
  is4 >> basis_set;
  is4.close();
  BSet basis(Mol1,basis_set);
  
  //allocation of GFN2 object
  GFN2 electron(basis,Mol1);
  electron.setElectronTemp(Telec);
  
  //to use ALPB solvation model, uncomment one of the statements below
  //electron.setSolvent("water");
  //electron.setSolvent("acetone");
  //electron.setSolvent("acetonitrile");
  //electron.setSolvent("aniline");
  //electron.setSolvent("benzaldehyde");
  //electron.setSolvent("benzene");
  //electron.setSolvent("dichloromethane");
  //electron.setSolvent("chloroform");
  //electron.setSolvent("carbon disulfide");
  //electron.setSolvent("dioxane");
  //electron.setSolvent("dmf");
  //electron.setSolvent("dmso");
  //electron.setSolvent("ethanol");
  //electron.setSolvent("diethyl ether");
  //electron.setSolvent("ethyl acetate");
  //electron.setSolvent("furane");
  //electron.setSolvent("hexadecane");
  //electron.setSolvent("hexane");
  //electron.setSolvent("methanol");
  //electron.setSolvent("nitromethane");
  //electron.setSolvent("octanol");
  //electron.setSolvent("phenol");
  //electron.setSolvent("thf");
  //electron.setSolvent("toluene");
  //electron.setSolvent("water");
  //electron.setSolvent("octanol wet");
  std::ifstream is5 ("solvent.tmp");
  std::string solvent;
  is5 >> solvent;
  is5.close();
  // if (strcmp(solvent, "without_solvent") != 0) {
  if (solvent != "without_solvent") {
  electron.setSolvent(solvent);
  }
  
  //SCF, with no print (silent); 0 or 1
  std::ifstream is6 ("verbosity.tmp");
  int verbosity;
  is6 >> verbosity;
  is6.close();
  electron.Calculate(verbosity);
  
  std::ifstream is23 ("opt.tmp");
  int opt_option;
  is23 >> opt_option;
  is23.close();
  
  std::ifstream is123 ("ts.tmp");
  int ts_option;
  is123 >> ts_option;
  is123.close();
  
  if (opt_option == 1) {
  
  // BFGSd solve(hessian_update,bfgs_allocation);
  if (rfo == 0) {
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //BFGSd solve(4,6);
  BFGSd solve(hessian_update,bfgs_allocation);
   std::ifstream is18 ("energy_threshold.tmp");
  double energy_threshold;
  is18 >> energy_threshold;
  is18.close();
  
  std::ifstream is19 ("gradient_threshold.tmp");
  double gradient_threshold;
  is19 >> gradient_threshold;
  is19.close();
  
  std::cout << "Geometry optimization" << std::endl;
  
  //double energy_threshold = 5.0e-6;
  //double gradient_threshold = 1.0e-3;
  SolverOpt(electron,solve,4,0,energy_threshold,gradient_threshold);
  } else {
  //other solvers: 
  BakerRFO solve(solver,0,0);
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
   std::ifstream is18 ("energy_threshold.tmp");
  double energy_threshold;
  is18 >> energy_threshold;
  is18.close();
  
  std::ifstream is19 ("gradient_threshold.tmp");
  double gradient_threshold;
  is19 >> gradient_threshold;
  is19.close();
  
  std::cout << "Geometry optimization (Baker RFO solver)" << std::endl;
  
  //double energy_threshold = 5.0e-6;
  //double gradient_threshold = 1.0e-3;
  SolverOpt(electron,solve,4,0,energy_threshold,gradient_threshold);
  }
  
  std::cout << std::setprecision(10);
  
  //get the molecule with optimized geometry
  Molecule Mol2 = electron.Component();
  //write geometry to xyz
  std::string geomfilename = argv[1];
  size_t extpos = geomfilename.find(".xyz");
  size_t lengeom = geomfilename.size();
  geomfilename.replace(extpos,lengeom,"");
  geomfilename += "_opt";
  Mol2.WriteXYZ(geomfilename);
  //this writes the optimized geometry to the file adenine_opt.xyz
  
  }
  
  if (ts_option == 1) {
  
  std::ifstream is14 ("bfgs.tmp");
  int bfgs_allocation;
  is14 >> bfgs_allocation;
  is14.close();
  
  std::ifstream is15 ("hessian_update.tmp");
  int hessian_update;
  is15 >> hessian_update;
  is15.close();
  
  std::ifstream is16 ("solver.tmp");
  int solver;
  is16 >> solver;
  is16.close();
  
  std::ifstream is17 ("rfo.tmp");
  int rfo;
  is17 >> rfo;
  is17.close();
  
  std::cout << "Transition state geometry optimization" << std::endl;
  
  std::ifstream is20 ("energy_threshold.tmp");
  double energy_threshold;
  is20 >> energy_threshold;
  is20.close();
  
  std::ifstream is21 ("gradient_threshold.tmp");
  double gradient_threshold;
  is21 >> gradient_threshold;
  is21.close();
  
  OptimizeTS(electron,0,0,energy_threshold,gradient_threshold);
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  }  
  
  // properties
  
  // molecular orbitals
  
  size_t AOcount = electron.NAO();
  // std::cout << "Number of AOs: " << std::endl;
  // std::cout << AOcount << std::endl << std::endl;
  
  std::cout << std::endl;
  std::cout << "Energy levels: " << std::endl;
  for(int i=1; i<=AOcount; ++i)
  {
  std::vector<double> MOorb;
  double Eorb_i = electron.getOrbital(i, MOorb);
  std::cout << Eorb_i << std::endl;
  }
  
  std::cout << std::endl;
  
  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO: " << Ehomo << std::endl << std::endl;
  std::cout << "LUMO: " << Elumo << std::endl << std::endl;
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
      
  // std::ifstream is19 ("charge_types.tmp");
  // std::string charge_types;
  // is19 >> charge_types;
  // is19.close();
  
  //get the partial charges
  std::vector<double> AtmCharge = electron.getQAtoms();
  
  //get the atom list
  std::vector<size_t> atoms = Mol1.Atoms();
  size_t Natoms = atoms.size();
  
  //AO basis info
  std::vector<size_t> AOS = basis.AtomNAOs(atoms);
  
  //get the atomic polarizabilities
  std::vector<double> polarizabilities;
  electron.AtomicPolarizabilities(polarizabilities,AtmCharge);
  
  //print
  std::cout << std::setprecision(5) << "\n";
  std::cout << "atom       AOs          charge           pol\n";
  for (size_t idx = 0; idx < Natoms; ++idx) {
    std::cout << AtomNr2Symbol(atoms[idx]) << "          ";
    std::cout << AOS[idx] << "          ";
    if (AtmCharge[idx] > 0.0) {std::cout << " ";}
    std::cout << AtmCharge[idx] << "          " << polarizabilities[idx] << "\n";
  }
  std::cout << "\n";
  
  double polbity = 0.0;
  electron.TotalPolarizability(polbity,AtmCharge);
  std::cout << " Total Polarizability          " << polbity << "\n";
  
  std::ifstream is26 ("elec_rx.tmp");
  int elec_rx;
  is26 >> elec_rx;
  is26.close();
  
  std::ifstream is27 ("orb_rx.tmp");
  int orb_rx;
  is27 >> orb_rx;
  is27.close();
  
  std::ifstream is28 ("koopman_ip.tmp");
  int koopman_ip;
  is28 >> koopman_ip;
  is28.close();
  
  std::ifstream is29 ("ip.tmp");
  int ip;
  is29 >> ip;
  is29.close();
  
  std::ifstream is30 ("ea.tmp");
  int ea;
  is30 >> ea;
  is30.close();
  
  std::ifstream is31 ("electronegativity.tmp");
  int electronegativity;
  is31 >> electronegativity;
  is31.close();
  
  std::ifstream is32 ("hardness.tmp");
  int hardness;
  is32 >> hardness;
  is32.close();
  
  //additional properties
  matrixE RxData(1,1);
  if (elec_rx > 0) {
    electron.ReactivityIndices(RxData,false);
    std::cout << ">Electronic Reactivity indices" << std::endl;
    RxData.Print(4);
    std::cout << "<Electronic Reactivity indices" << std::endl;
  }

  if (orb_rx > 0) {
    electron.ReactivityIndices(RxData,true);
    std::cout << ">Orbital Reactivity indices" << std::endl;
    RxData.Print(4);
    std::cout << "<Orbital Reactivity indices" << std::endl;
  }

  if (koopman_ip > 0) {
  std::cout << "" << std::endl;
  std::cout << "Ionization Potential (Koopman): " << electron.IonizationPotential(true)*au2eV << "   eV" << std::endl;
  std::cout << "" << std::endl;
  }
  if (ip > 0) {
  std::cout << "empirical IP shift " << "4.8455" << "   eV" << std::endl;
  std::cout << "Ionization Potential (Definition): " << electron.IonizationPotential(false)*au2eV-4.8455 << "   eV" << std::endl;
  std::cout << "" << std::endl;
  }
  if (ea > 0) {
  std::cout << "empirical EA shift " << "4.8455" << "   eV" << std::endl;
  std::cout << "Electron Affinity (Definition): " << electron.ElectronAffinity()*au2eV-4.8455 << "   eV" << std::endl;
  std::cout << "" << std::endl;
  }
  
  if ((electronegativity > 0)||(hardness > 0)) {
    double chi;
    double eta;
    electron.HSABdata(chi,eta);
    std::cout << "Electronegativity: " << chi*au2eV << "   eV" << std::endl;
    std::cout << "Hardness: " << eta*au2eV << "   eV" << std::endl;
  }
  
  std::ifstream is22 ("thermo.tmp");
  int thermo_prop;
  is22 >> thermo_prop;
  is22.close();
  
  if (thermo_prop == 1){
  
  std::ifstream is18 ("radical.tmp");
  std::string radical;
  is18 >> radical;
  is18.close();
  
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  //double T = 298.15;
  std::ifstream is13 ("temperature.tmp");
  double T;
  is13 >> T;
  is13.close();
  std::ifstream is33 ("grimme_corr.tmp");
  int grimme_corr;
  is33 >> grimme_corr;
  is33.close();
  
  bool grimmecorrection;
  if (grimme_corr  == 1){
  	grimmecorrection = true;
  } else {
  	grimmecorrection = false;
  }
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,point_group,radical,grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)  U(J/mol)        P(Pa)           H(J/mol)        CV(J/[K.mol])   CP(J/[K.mol])   S(J/[K.mol])    A(J/mol)        G(J/mol)" << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.U() << "  " << IdealGas.P() << "  " << IdealGas.H() << "  " << IdealGas.CV() << "  " << IdealGas.CP() << "  " << IdealGas.S() << "  " << IdealGas.A() << "  " << IdealGas.G() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  }
  std::ifstream is24 ("calc_density.tmp");
  int calc_density;
  is24 >> calc_density;
  is24.close();
  
  std::ifstream is25 ("density_file.tmp");
  std::string density_file;
  is25 >> density_file;
  is25.close();
  //get the density?
  if (calc_density > 0) {
    electron.ElectronicDensity(density_file);
  }
  // end of gfn2
  
  }
  
  return 0;
}
