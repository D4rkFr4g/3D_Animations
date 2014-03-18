#pragma once
#ifndef RIGTFORM_H
#define RIGTFORM_H

#include <iostream>
#include <cassert>

#include "quat.h"
#include "matrix4.h"

class RigTForm {
  Cvec3 t_; // translation component
  Quat r_;  // rotation component represented as a quaternion

public:
  RigTForm() : t_(0) {
    assert(norm2(Quat(1,0,0,0) - r_) < CS175_EPS2);
  }

  RigTForm(const Cvec3& t, const Quat& r) 
  {
	  t_ = t;
	  r_ = r;
  }

  explicit RigTForm(const Cvec3& t) 
  {
    t_ = t;
	 r_ = Quat();
  }

  explicit RigTForm(const Quat& r) 
  {
	  t_ = Cvec3(0,0,0);
	  r_ = r;
  }

  Cvec3 getTranslation() const 
  {
    return t_;
  }

  Quat getRotation() const 
  {
    return r_;
  }

  RigTForm& setTranslation(const Cvec3& t) {
    t_ = t;
    return *this;
  }

  RigTForm& setRotation(const Quat& r) {
    r_ = r;
    return *this;
  }

  Cvec4 operator * (const Cvec4& a) const 
  {
	  return (r_ * a) + Cvec4(t_,0);
  }

  RigTForm operator * (const RigTForm& a) const 
  {
	  return RigTForm(t_ + Cvec3(r_ * Cvec4(a.t_,0)),r_ * a.r_);
  }

  static Matrix4 makeTRmatrix(const RigTForm& rbt, Matrix4 s)
	{
		Matrix4 t = Matrix4::makeTranslation(rbt.getTranslation());
		Matrix4 r = quatToMatrix(rbt.getRotation());
		
		return t * r * s;
	}
};

inline RigTForm inv(const RigTForm& tform) 
{
	Quat quatInverse = inv(tform.getRotation());
	return RigTForm(Cvec3(quatInverse * Cvec4(-tform.getTranslation(),1)),quatInverse);
}

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}

inline Matrix4 rigTFormToMatrix(const RigTForm& tform) 
{
	Matrix4 t = Matrix4::makeTranslation(tform.getTranslation());
	Matrix4 r =  quatToMatrix(tform.getRotation());
	Matrix4 m = t * r;

	return m;
}

#endif
