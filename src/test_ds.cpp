#include "Datastruct.h"
#include <iostream>
#include <fstream>
#include <Rcpp.h>
#include "debug.h"

DataStruct DataFrameToDataStruct(Rcpp::DataFrame DF);
Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);

// [[Rcpp::export]]
Rcpp::DataFrame test_ds(Rcpp::DataFrame DF, Rcpp::IntegerVector In) {
  DataStruct DS = DataFrameToDataStruct(DF);
  DataStruct DS2(DS, In);
  return DataStructToDataFrame(DS2);
}

// [[Rcpp::export]]
Rcpp::DataFrame test_ds2() {
  // premier colonne : entiers
  std::vector<int> v;
  v.push_back(12);
  Column col(v);
  col.push_back(1);

  // deuxième colonne
  Column col2(datatype::DOUBLE);
  col2.push_back(3.14);
  col2.push_back(2.71);

  // troisième colonne
  Column col3(datatype::STRING);
  col3.push_back(std::string("choin"));
  col3.push_back(std::string("choin"));

  // data struct à 3 deux colonnes
  DataStruct DS;
  DS.push_back(col,  "a");
  DS.push_back(col2, "b");
  DS.push_back(col3, "c");

  DS.at(0).push_back(7);
  DS.at(1).push_back(7.);
  DS.at(2).push_back(std::string("poupou"));

  std::string s = "123.99";
  DS.at(0).push_back_convert(s);
  DS.at(1).push_back_convert(s);
  DS.at(2).push_back_convert(s);

  const char * a = "45e2";

  DS.at(0).push_back_convert(a);
  DS.at(1).push_back_convert(a);
  DS.at(2).push_back_convert(a);

  const char * b = " 74   12.74 lapin 022";
  char * z; 
  z = DS.at(0).push_back_token(b);
  z = DS.at(1).push_back_token(z);
  z = DS.at(2).push_back_token(z);

  return DataStructToDataFrame(DS);
}

