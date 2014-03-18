/****************************************************** 
* Copyright (c):   2013, All Rights Reserved. 
* Project:         CS 116A Homework #3 
* File:            RigTFormObj.cpp 
* Purpose:         Manage a object's transRtForm, scaleRigTForm, rotateRigTForm
* Start date:      10/28/13 
* Programmer:      Zane Melcho 
* 
****************************************************** 
*/
#include "cvec.h"
#include "matrix4.h"
#include "rigtform.h"

using namespace std;
/*
* RigidBody is a class that holds an RigTForm and a Matrix4 scale
*
*/
//Declarations
namespace
{
	class RigidBody
	{ 
	 public:
		RigidBody();
		RigidBody(RigTForm rtf, Matrix4 scale);
		RigTForm getRigidTForm();
		Matrix4 getScale();
		RigTForm getRigidTFormRespectToParent();

	 private:
		RigTForm rtf_;
		Matrix4 scale_;
		RigTForm parent_;
		Cvec3 color_;
		Geometry geom_;
	};

	//Implementations

	//Default Constructor
	RigidBody::RigidBody()
	{
		/* PURPOSE:		
			RETURNS:    
		*/
		rtf_ = RigTForm();
		scale_ = Matrix4();
		parent_ = RigTForm();
	}
	/*-----------------------------------------------*/
	
	RigidBody::RigidBody(RigTForm rtf, Matrix4 scale)
	{
		/* PURPOSE:		 
			RECEIVES:	 
							
			RETURNS:		
		*/

		rtf_ = rtf;
		scale_ = scale;
		parent_ = RigTForm();
	}

	RigidBody::RigidBody(RigTForm rtf, Matrix4 scale, RigTForm parent)
	{
		/* PURPOSE:		 
			RECEIVES:	 
							
			RETURNS:		
		*/

		rtf_ = rtf;
		scale_ = scale;
		parent_ = parent;
	}

	RigTForm RigidBody::getRigidTForm()
	{
		return rtf_;
	}
	
	Matrix4 RigidBody::getScale()
	{
		return scale_;
	}

	RigTForm RigidBody::getRigidTFormRespectToParent()
	{
		return parent_ * rtf_;
	}
}