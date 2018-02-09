/*
 * NtupleConverter.h
 *
 *  Created on: 8 Feb 2018
 *      Author: ebols
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLECONVERTER_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLECONVERTER_H_


#include <vector>
#include <string>
#include "TString.h"
#include "TChain.h"
#include "LumiReweightingStandAlone.h"



bool isMC;

Float_t         Jet_DeepFlavourBDisc;   //[nJet]                                                                                                 
Float_t         Jet_DeepFlavourCvsLDisc;   //[nJet]                                                                                              
Float_t         Jet_DeepFlavourCvsBDisc;   //[nJet]                                                                                              
Float_t         Jet_DeepFlavourB;   //[nJet]                                                                                                     
Float_t         Jet_DeepFlavourBB;   //[nJet]                                                                                                    
Float_t         Jet_DeepFlavourLEPB;   //[nJet]                                                                                                  
Float_t         Jet_DeepFlavourC;   //[nJet]                                                                                                     
Float_t         Jet_DeepFlavourUDS;   //[nJet]                                                                                                   
Float_t         Jet_DeepFlavourG;   //[nJet]                                                                                                     
Float_t         Jet_DeepCSVBDisc;   //[nJet]                                                                                                     
Float_t         Jet_DeepCSVBDiscN;   //[nJet]                                                                                                    
Float_t         Jet_DeepCSVBDiscP;   //[nJet]                                                                                                    
Float_t         Jet_DeepCSVCvsLDisc;   //[nJet]                                                                                                  
Float_t         Jet_DeepCSVCvsLDiscN;   //[nJet]                                                                                                 
Float_t         Jet_DeepCSVCvsLDiscP;   //[nJet]                                                                                                 
Float_t         Jet_DeepCSVCvsBDisc;   //[nJet]                                                                                                  
Float_t         Jet_DeepCSVCvsBDiscN;   //[nJet]                                                                                                 
Float_t         Jet_DeepCSVCvsBDiscP;   //[nJet]                                                                                                 
Float_t         Jet_DeepCSVb;   //[nJet]                                                                                                         
Float_t         Jet_DeepCSVc;   //[nJet]                                                                                                         
Float_t         Jet_DeepCSVl;   //[nJet]                                                                                                         
Float_t         Jet_DeepCSVbb;   //[nJet]                                                                                                        
Float_t         Jet_DeepCSVcc;   //[nJet]                                                                                                        
Float_t         Jet_DeepCSVbN;   //[nJet]                                                                                                        
Float_t         Jet_DeepCSVcN;   //[nJet]                                                                                                        
Float_t         Jet_DeepCSVlN;   //[nJet]               
Float_t         Jet_DeepCSVbbN;   //[nJet]
Float_t         Jet_DeepCSVccN;   //[nJet]
Float_t         Jet_DeepCSVbP;   //[nJet]
Float_t         Jet_DeepCSVcP;   //[nJet]
Float_t         Jet_DeepCSVlP;   //[nJet]
Float_t         Jet_DeepCSVbbP;   //[nJet]
Float_t         Jet_DeepCSVccP;   //[nJet]
Float_t         Jet_pt;   //[nJet]

Int_t           isB;
Int_t           isC;
Int_t           isUDSG;


Int_t           nBitTrigger;
Int_t           BitTrigger[3];
Float_t         pthat;
Int_t           nPV;


TString gentype;
bool qcdtype;
int sqrtstev;

float WeightPU;
float weightXS;

double x_section[12]; 
double nmc_evt_vect[12];
float sum_xs;

reweight::LumiReWeighting LumiWeights;


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLECONVERTER_H_ */
