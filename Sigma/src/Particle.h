#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>


class Particle{

 private:
 
  std::string _name;  //particle name
  int    _id;         //particle id according to PDF particle codes
  double _mass;       //particle mass
  double _charge;     //electric charge in the unit of the elementary charge e
  
  
 public:

  // enum _type {d=1,u=2,s=3,c=4,b=5,t=6,g=21,photon=22,Z0=23,Wp=24,parton=0};

  Particle();
  Particle(const std::string name);
  ~Particle(){}
  void setParticle(const std::string name);

  inline std::string Name() const {return _name;}
  inline int    ID()        const {return _id;}
  inline double Mass()      const {return _mass;}
  inline double Charge()    const {return _charge;}
  
};

#endif  //PARTICLE_H