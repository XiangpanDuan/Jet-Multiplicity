#include "Particle.h"


Particle::Particle(){
  _name="parton"; _id=-273; _mass=0.0; _charge=0.0;
}

Particle::Particle(const std::string name){
  _name=name; _id=-273; _mass=0.0; _charge= 0.0;
  setParticle(_name);
}

void Particle::setParticle(const std::string name){
  //All values are taken from PDG 2024
  if      (name=="parton")   {_id=-273;   _mass=0.0;       _charge=0.0;}
  else if (name=="d")        {_id=1;      _mass=0.0;       _charge=-1.0/3.0;}
  else if (name=="u")        {_id=2;      _mass=0.0;       _charge= 2.0/3.0;}
  else if (name=="s")        {_id=3;      _mass=0.0;       _charge=-1.0/3.0;}
  else if (name=="c")        {_id=4;      _mass=0.0;       _charge= 2.0/3.0;}
  else if (name=="b")        {_id=5;      _mass=0.0;       _charge=-1.0/3.0;}
  else if (name=="t")        {_id=6;      _mass=173.0;     _charge= 2.0/3.0;}
  else if (name=="g")        {_id=21;     _mass=0.0;       _charge= 0.0;}
  else if (name=="photon")   {_id=22;     _mass=0.0;       _charge= 0.0;}
  else if (name=="Z0")       {_id=23;     _mass=91.188;    _charge= 0.0;}
  else if (name=="Wp")       {_id=24;     _mass=80.3692;   _charge= 1.0;}
  else if (name=="Higgs")    {_id=25;     _mass=125.2;     _charge= 0.0;}
}
