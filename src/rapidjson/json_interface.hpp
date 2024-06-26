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

#ifndef _json_interface_
#define _json_interface_
#include "document.h"
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

//description:
//interface to RapidJSON

//simple interfaces to rapidjson so that I can read json type files
//each interface is database specific
void DebyeJSON(std::string member, std::vector<std::string> * rstring, std::vector<double> * rdouble, int buffer = 256) {
  //interface function to read Debye parameters in /db/solids.json
  std::string jsonfile = "../db/debye.json";
  std::ifstream jfile(jsonfile, std::ios::in);
  if (!jfile.is_open()) {throw std::string("ERROR: json_interface.hpp: DebyeJSON(): Error, json file could not be open");}
  std::string readjson = "";
  char aux[buffer];
  while (!jfile.eof()) {
    jfile >> aux;
    readjson += aux;
  }
  jfile.close();
  rapidjson::Document doc;
  if (doc.Parse(readjson.c_str()).HasParseError()) {throw std::string("ERROR: json_interface.hpp: DebyeJSON(): Error, parse error.");}
  std::cout << "parsed file" << std::endl;
  assert(doc.IsObject());
  //check whether structure is in DB
  assert(doc.HasMember(member.c_str()));
  std::cout << member << " exists" << std::endl;
  rstring->push_back(doc[member.c_str()]["formula"].GetString());
  rstring->push_back(doc[member.c_str()]["date"].GetString());
  rstring->push_back(doc[member.c_str()]["source"].GetString());
  rdouble->push_back(doc[member.c_str()]["Phi0"].GetDouble());
  rdouble->push_back(doc[member.c_str()]["MW"].GetDouble());
  rdouble->push_back(doc[member.c_str()]["density"].GetDouble());
  rdouble->push_back(doc[member.c_str()]["speedsound"].GetDouble());
  rdouble->push_back(doc[member.c_str()]["thetaD"].GetDouble());
}
void PBlRRlHOEJSON(std::string database, std::vector<std::string> member, std::vector<std::string> * rstring, std::vector<double> * rdouble, std::vector<std::vector<double> > * rdvect, int & nconf, int buffer = 256) {
  //interface function to read Debye parameters in /db/X.json, where X is a molecule's name
  //rstring contains string data from database in a vector
  //rdouble contains double data from database in a vector
  //rdvect contains vector data from database in a vector
  //member = (file, method, basis set)       for method = exp then there is no basis set
  std::string jsonfile = "../db/" + database + "/" + member[0] + ".json";
  std::cout << "reading file " << jsonfile << std::endl;
  std::ifstream jfile(jsonfile, std::ios::in);
  if (!jfile.is_open()) {throw std::string("ERROR: json_interface.hpp: PBlRRlHOEJSON(): json file could not be open");}
  std::string readjson = "";
  char aux[buffer];
  while (!jfile.eof()) {
    jfile >> aux;
    readjson += aux;
  }
  jfile.close();
  rapidjson::Document doc;
  if (doc.Parse(readjson.c_str()).HasParseError()) {throw std::string("ERROR: json_interface.hpp: PBlRRlHOEJSON(): parsing error");}
  std::cout << "parsed file" << std::endl;
  assert(doc.IsObject());
  //check whether method is in DB
  assert(doc.HasMember(member[1].c_str()));
  std::cout << member[1] << " exists" << std::endl;
  if ((member[1] == "exp")||(member[1] == "EXP")) {
    nconf = 1;
    //for system with experimental data
    //strings
    rstring->push_back(doc["exp"]["source"].GetString());
    rstring->push_back(doc["exp"]["date"].GetString());
    rstring->push_back(doc["exp"]["geometry"].GetString());
    rstring->push_back(doc["exp"]["symmetry"].GetString());
    rstring->push_back(doc["exp"]["rotationalsymmetry"].GetString());
    //doubles
    rdouble->push_back(doc["exp"]["DfH"].GetDouble());
    //vectors
    const rapidjson::Value& vector1 = doc["exp"]["I"]; // Using a reference for consecutive access is handy and faster.
    std::vector<double> aux;
    assert(vector1.IsArray());
    for (rapidjson::SizeType i = 0; i < vector1.Size(); ++i) {
      aux.push_back(vector1[i].GetDouble());  
    }
    rdvect->push_back(aux);
    const rapidjson::Value& vector2 = doc["exp"]["vibr"];
    aux.clear();
    assert(vector2.IsArray());
    for (rapidjson::SizeType i = 0; i < vector2.Size(); ++i) {
      aux.push_back(vector2[i].GetDouble());  
    }
    rdvect->push_back(aux);
    const rapidjson::Value& vector3 = doc["exp"]["Eel"];
    aux.clear();
    assert(vector3.IsArray());
    for (rapidjson::SizeType i = 0; i < vector3.Size(); ++i) {
      aux.push_back(vector3[i].GetDouble());  
    }
    rdvect->push_back(aux);
    const rapidjson::Value& vector4 = doc["exp"]["gel"];
    aux.clear();
    assert(vector4.IsArray());
    for (rapidjson::SizeType i = 0; i < vector4.Size(); ++i) {
      aux.push_back(vector4[i].GetDouble());  
    }
    rdvect->push_back(aux);
  }
  else if ((member[1] == "qm9")||(member[1] == "QM9")) {
    nconf = 1;
    //for system with experimental data
    //strings
    rstring->push_back(doc["qm9"]["source"].GetString());
    rstring->push_back(doc["qm9"]["date"].GetString());
    rstring->push_back(doc["qm9"]["geometry"].GetString());
    rstring->push_back(doc["qm9"]["symmetry"].GetString());
    rstring->push_back(doc["qm9"]["rotationalsymmetry"].GetString());
    //doubles
    rdouble->push_back(doc["qm9"]["DfH"].GetDouble());
    //vectors
    const rapidjson::Value& vector1 = doc["qm9"]["I"]; // Using a reference for consecutive access is handy and faster.
    std::vector<double> aux;
    assert(vector1.IsArray());
    for (rapidjson::SizeType i = 0; i < vector1.Size(); ++i) {
      aux.push_back(vector1[i].GetDouble());  
    }
    rdvect->push_back(aux);
    const rapidjson::Value& vector2 = doc["qm9"]["vibr"];
    aux.clear();
    assert(vector2.IsArray());
    for (rapidjson::SizeType i = 0; i < vector2.Size(); ++i) {
      aux.push_back(vector2[i].GetDouble());  
    }
    rdvect->push_back(aux);
    const rapidjson::Value& vector3 = doc["qm9"]["Eel"];
    aux.clear();
    assert(vector3.IsArray());
    for (rapidjson::SizeType i = 0; i < vector3.Size(); ++i) {
      aux.push_back(vector3[i].GetDouble());  
    }
    rdvect->push_back(aux);
    const rapidjson::Value& vector4 = doc["qm9"]["gel"];
    aux.clear();
    assert(vector4.IsArray());
    for (rapidjson::SizeType i = 0; i < vector4.Size(); ++i) {
      aux.push_back(vector4[i].GetDouble());  
    }
    rdvect->push_back(aux);
  }
  else {
    //theoretical calculations
    //here there might be conformers
    int nconformer = doc[member[1].c_str()][member[2].c_str()]["conf"].GetInt();
    nconf = nconformer;
    rdouble->push_back(doc[member[1].c_str()][member[2].c_str()]["DfH"].GetDouble());
    //saving metadata
    rstring->push_back(doc[member[1].c_str()][member[2].c_str()]["source"].GetString());
    rstring->push_back(doc[member[1].c_str()][member[2].c_str()]["date"].GetString());
    std::string conf;
    std::ostringstream ostr;
    for (size_t idx = 0; idx < nconformer; ++idx) {
      ostr << idx+1;
      conf = "c"+ostr.str();
      rstring->push_back(doc[member[1].c_str()][member[2].c_str()][conf.c_str()]["geometry"].GetString());
      rstring->push_back(doc[member[1].c_str()][member[2].c_str()][conf.c_str()]["symmetry"].GetString());
      rstring->push_back(doc[member[1].c_str()][member[2].c_str()][conf.c_str()]["rotationalsymmetry"].GetString());
      const rapidjson::Value& vector1 = doc[member[1].c_str()][member[2].c_str()][conf.c_str()]["I"];
      std::vector<double> aux;
      assert(vector1.IsArray());
      for (rapidjson::SizeType i = 0; i < vector1.Size(); ++i) {
        aux.push_back(vector1[i].GetDouble());  
      }
      rdvect->push_back(aux);
      const rapidjson::Value& vector2 = doc[member[1].c_str()][member[2].c_str()][conf.c_str()]["vibr"];
      aux.clear();
      assert(vector2.IsArray());
      for (rapidjson::SizeType i = 0; i < vector2.Size(); ++i) {
        aux.push_back(vector2[i].GetDouble());  
      }
      rdvect->push_back(aux);
      const rapidjson::Value& vector3 = doc[member[1].c_str()][member[2].c_str()][conf.c_str()]["Eel"];
      aux.clear();
      assert(vector3.IsArray());
      for (rapidjson::SizeType i = 0; i < vector3.Size(); ++i) {
        aux.push_back(vector3[i].GetDouble());  
      }
      rdvect->push_back(aux);
      const rapidjson::Value& vector4 = doc[member[1].c_str()][member[2].c_str()][conf.c_str()]["gel"];
      aux.clear();
      assert(vector4.IsArray());
      for (rapidjson::SizeType i = 0; i < vector4.Size(); ++i) {
        aux.push_back(vector4[i].GetDouble());  
      }
      rdvect->push_back(aux);
      ostr.str(std::string());
    }
  }
}

#endif //_json_interface_
