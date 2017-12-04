#ifndef EVENTBUFFER_H
#define EVENTBUFFER_H
//----------------------------------------------------------------------------
// File:        eventBuffer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Mon Dec  4 00:43:33 2017 by mkanalyzer.py
// Author:      Shakespeare's ghost
//----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <map>
#include <cassert>
#include "treestream.h"

struct eventBuffer
{
  //--------------------------------------------------------------------------
  // --- Declare variables
  //--------------------------------------------------------------------------
  std::vector<double>	Particle_barcode;
  std::vector<int>	Particle_d1;
  std::vector<int>	Particle_d2;
  std::vector<double>	Particle_energy;
  std::vector<double>	Particle_mass;
  std::vector<int>	Particle_pid;
  std::vector<double>	Particle_px;
  std::vector<double>	Particle_py;
  std::vector<double>	Particle_pz;
  std::vector<int>	Particle_status;

  double	Event_alphaQCD;
  double	Event_alphaQED;
  int	Event_barcodeBP1;
  int	Event_barcodeBP2;
  int	Event_barcodeSPV;
  int	Event_number;
  int	Event_numberMP;
  int	Event_numberP;
  int	Event_numberV;
  double	Event_scale;
  double	PDF_Q2;
  int	PDF_id1;
  int	PDF_id2;
  int	PDF_parton1;
  int	PDF_parton2;
  double	PDF_x1;
  double	PDF_x1f;
  double	PDF_x2;
  double	PDF_x2f;
  double	Xsection_error;
  double	Xsection_value;

  //--------------------------------------------------------------------------
  // --- Structs can be filled by calling fill(), or individual fill
  // --- methods, e.g., fillElectrons()
  // --- after the call to read(...)
  //----------- --------------------------------------------------------------
  struct Particle_s
  {
    double	barcode;
    int	d1;
    int	d2;
    double	energy;
    double	mass;
    int	pid;
    double	px;
    double	py;
    double	pz;
    int	status;

    std::ostream& operator<<(std::ostream& os)
    {
      char r[1024];
      os << "Particle" << std::endl;
      sprintf(r, "  %-32s: %f\n", "barcode", ( double)barcode); os << r;
      sprintf(r, "  %-32s: %f\n", "d1", ( double)d1); os << r;
      sprintf(r, "  %-32s: %f\n", "d2", ( double)d2); os << r;
      sprintf(r, "  %-32s: %f\n", "energy", ( double)energy); os << r;
      sprintf(r, "  %-32s: %f\n", "mass", ( double)mass); os << r;
      sprintf(r, "  %-32s: %f\n", "pid", ( double)pid); os << r;
      sprintf(r, "  %-32s: %f\n", "px", ( double)px); os << r;
      sprintf(r, "  %-32s: %f\n", "py", ( double)py); os << r;
      sprintf(r, "  %-32s: %f\n", "pz", ( double)pz); os << r;
      sprintf(r, "  %-32s: %f\n", "status", ( double)status); os << r;
      return os;
    }
  };


  void fillParticles()
  {
    Particle.resize(Particle_barcode.size());
    for(unsigned int i=0; i < Particle.size(); ++i)
      {
        Particle[i].barcode	= Particle_barcode[i];
        Particle[i].d1	= Particle_d1[i];
        Particle[i].d2	= Particle_d2[i];
        Particle[i].energy	= Particle_energy[i];
        Particle[i].mass	= Particle_mass[i];
        Particle[i].pid	= Particle_pid[i];
        Particle[i].px	= Particle_px[i];
        Particle[i].py	= Particle_py[i];
        Particle[i].pz	= Particle_pz[i];
        Particle[i].status	= Particle_status[i];
      }
  }


  std::vector<eventBuffer::Particle_s> Particle;

  void fillObjects()
  {
    fillParticles();
  }

  //--------------------------------------------------------------------------
  // Select objects for which the select function was called
  void saveObjects()
  {
    int n;

    n = 0;
    try
      {
         n = indexmap["Particle"].size();
      }
    catch (...)
      {}
    if ( n > 0 )
      {
        std::vector<int>& index = indexmap["Particle"];
        for(int i=0; i < n; ++i)
          {
            int j = index[i];
            Particle_barcode[i]	= Particle_barcode[j];
            Particle_d1[i]	= Particle_d1[j];
            Particle_d2[i]	= Particle_d2[j];
            Particle_energy[i]	= Particle_energy[j];
            Particle_mass[i]	= Particle_mass[j];
            Particle_pid[i]	= Particle_pid[j];
            Particle_px[i]	= Particle_px[j];
            Particle_py[i]	= Particle_py[j];
            Particle_pz[i]	= Particle_pz[j];
            Particle_status[i]	= Particle_status[j];
          }
      }
    Event_numberP = n;
  }

  //--------------------------------------------------------------------------
  // A read-only buffer 
  eventBuffer() : input(0), output(0) {}
  eventBuffer(itreestream& stream)
  : input(&stream),
    output(0)
  {
    if ( !input->good() ) 
      {
        std::cout << "eventBuffer - please check stream!" 
                  << std::endl;
	exit(0);
      }
    initBuffers();

    input->select("Event_alphaQCD", 	Event_alphaQCD);
    input->select("Event_alphaQED", 	Event_alphaQED);
    input->select("Event_barcodeBP1", 	Event_barcodeBP1);
    input->select("Event_barcodeBP2", 	Event_barcodeBP2);
    input->select("Event_barcodeSPV", 	Event_barcodeSPV);
    input->select("Event_number", 	Event_number);
    input->select("Event_numberMP", 	Event_numberMP);
    input->select("Event_numberP", 	Event_numberP);
    input->select("Event_numberV", 	Event_numberV);
    input->select("Event_scale", 	Event_scale);
    input->select("PDF_Q2", 	PDF_Q2);
    input->select("PDF_id1", 	PDF_id1);
    input->select("PDF_id2", 	PDF_id2);
    input->select("PDF_parton1", 	PDF_parton1);
    input->select("PDF_parton2", 	PDF_parton2);
    input->select("PDF_x1", 	PDF_x1);
    input->select("PDF_x1f", 	PDF_x1f);
    input->select("PDF_x2", 	PDF_x2);
    input->select("PDF_x2f", 	PDF_x2f);
    input->select("Particle_barcode", 	Particle_barcode);
    input->select("Particle_d1", 	Particle_d1);
    input->select("Particle_d2", 	Particle_d2);
    input->select("Particle_energy", 	Particle_energy);
    input->select("Particle_mass", 	Particle_mass);
    input->select("Particle_pid", 	Particle_pid);
    input->select("Particle_px", 	Particle_px);
    input->select("Particle_py", 	Particle_py);
    input->select("Particle_pz", 	Particle_pz);
    input->select("Particle_status", 	Particle_status);
    input->select("Xsection_error", 	Xsection_error);
    input->select("Xsection_value", 	Xsection_value);

  }

  // A write-only buffer
  eventBuffer(otreestream& stream)
  : input(0),
    output(&stream)
  {
    initBuffers();

    output->add("Event_numberP",	 Event_numberP);
  
    output->add("Event_alphaQCD", 	Event_alphaQCD);
    output->add("Event_alphaQED", 	Event_alphaQED);
    output->add("Event_barcodeBP1", 	Event_barcodeBP1);
    output->add("Event_barcodeBP2", 	Event_barcodeBP2);
    output->add("Event_barcodeSPV", 	Event_barcodeSPV);
    output->add("Event_number", 	Event_number);
    output->add("Event_numberMP", 	Event_numberMP);
    output->add("Event_numberV", 	Event_numberV);
    output->add("Event_scale", 	Event_scale);
    output->add("PDF_Q2", 	PDF_Q2);
    output->add("PDF_id1", 	PDF_id1);
    output->add("PDF_id2", 	PDF_id2);
    output->add("PDF_parton1", 	PDF_parton1);
    output->add("PDF_parton2", 	PDF_parton2);
    output->add("PDF_x1", 	PDF_x1);
    output->add("PDF_x1f", 	PDF_x1f);
    output->add("PDF_x2", 	PDF_x2);
    output->add("PDF_x2f", 	PDF_x2f);
    output->add("Particle_barcode[Event_numberP]", 	Particle_barcode);
    output->add("Particle_d1[Event_numberP]", 	Particle_d1);
    output->add("Particle_d2[Event_numberP]", 	Particle_d2);
    output->add("Particle_energy[Event_numberP]", 	Particle_energy);
    output->add("Particle_mass[Event_numberP]", 	Particle_mass);
    output->add("Particle_pid[Event_numberP]", 	Particle_pid);
    output->add("Particle_px[Event_numberP]", 	Particle_px);
    output->add("Particle_py[Event_numberP]", 	Particle_py);
    output->add("Particle_pz[Event_numberP]", 	Particle_pz);
    output->add("Particle_status[Event_numberP]", 	Particle_status);
    output->add("Xsection_error", 	Xsection_error);
    output->add("Xsection_value", 	Xsection_value);

  }

  void initBuffers()
  {
    Particle_barcode	= std::vector<double>(5937,0);
    Particle_d1	= std::vector<int>(5937,0);
    Particle_d2	= std::vector<int>(5937,0);
    Particle_energy	= std::vector<double>(5937,0);
    Particle_mass	= std::vector<double>(5937,0);
    Particle_pid	= std::vector<int>(5937,0);
    Particle_px	= std::vector<double>(5937,0);
    Particle_py	= std::vector<double>(5937,0);
    Particle_pz	= std::vector<double>(5937,0);
    Particle_status	= std::vector<int>(5937,0);
    Particle	= std::vector<eventBuffer::Particle_s>(5937);

  }
      
  void read(int entry)
  {
    if ( !input ) 
      { 
        std::cout << "** eventBuffer::read - first  call read-only constructor!"
                  << std::endl;
        assert(0);
      }
    input->read(entry);

    // clear indexmap
    for(std::map<std::string, std::vector<int> >::iterator
    item=indexmap.begin(); 
    item != indexmap.end();
    ++item)
    item->second.clear();
  }

  void select(std::string objname)
  {
    indexmap[objname] = std::vector<int>();
  }

  void select(std::string objname, int index)
  {
    try
     {
       indexmap[objname].push_back(index);
     }
    catch (...)
     {
       std::cout << "** eventBuffer::select - first call select(""" 
                 << objname << """)" 
                 << std::endl;
       assert(0);
    }
  }

 void ls()
 {
   if( input ) input->ls();
 }

 int size()
 {
   if( input ) 
     return input->size();
   else
     return 0;
 }

 void close()
 {
   if( input )   input->close();
   if( output ) output->close();
 }

 // --- indexmap keeps track of which objects have been flagged for selection
 std::map<std::string, std::vector<int> > indexmap;

 // to read events
 itreestream* input;

 // to write events
 otreestream* output;

}; 
#endif
