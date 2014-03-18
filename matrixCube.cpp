/****************************************************** 
* Copyright (c):   2013, All Rights Reserved. 
* Project:         CS 116A Homework #2 
* File:            MatrixCube.cpp 
* Purpose:         Manage a cube object's transRbt, scaleRbt, rotateRbt
* Start date:      9/28/13 
* Programmer:      Zane Melcho 
* 
****************************************************** 
*/
#include "cvec.h"
#include "matrix4.h"

using namespace std;
/*
* MatrixCube is a class that holds an Object, translation, rotation, and scale frame
*
*/
//Declarations
namespace
{
	class MatrixCube
	{ 
	 public:
		MatrixCube();
		MatrixCube(Matrix4 translate, Matrix4 scale);
		MatrixCube::MatrixCube(Matrix4 translate, Matrix4 rotate, Matrix4 scale);
		Matrix4 getTransRbt();
		void setTransRbt(Matrix4 T);
		void resetTransRbt();
		Matrix4 getScaleRbt();
		void setScaleRbt(Matrix4 S);
		void resetScaleRbt();
		Matrix4 getRotateRbt();
		void setRotateRbt(Matrix4 R);
		void resetRotateRbt();
		Matrix4 getObjRbt();
		void update(Matrix4 objRbt);
		void updateRotateRbt(Matrix4 objRbt);
		Matrix4 updateWithRespectToA(Matrix4 A);
		void originalPosition();

	 private:
		Matrix4 myObjRbt;
		Matrix4 myTransRbt;
		Matrix4 myScaleRbt;
		Matrix4 myRotateRbt;
		Matrix4 myOrigObjRbt;
		Matrix4 myOrigScaleRbt;
		Matrix4 myOrigTransRbt;
		Matrix4 myOrigRotateRbt;
	};

	//Implementations

	//Default Constructor
	MatrixCube::MatrixCube()
	{
		/* PURPOSE:		Constructs a default MatrixCube 
			RETURNS:    Returns a default MatrixCube 
		*/
		myObjRbt = Matrix4();
		myTransRbt = Matrix4();
		myScaleRbt = Matrix4();
		myRotateRbt = Matrix4();
		myOrigScaleRbt = Matrix4();
		myOrigTransRbt = Matrix4();
		myOrigRotateRbt = Matrix4();
	}
	/*-----------------------------------------------*/
	MatrixCube::MatrixCube(Matrix4 translate, Matrix4 scale)
	{
		/* PURPOSE:		Constructs a MatrixCube with given transRbt and scaleRbt 
			RECEIVES:	Matrix4 translate - a transRbt Matrix that has the cube's initial translation 
							Matrix4 scale - a scaleRbt Matrix that has the cube's initial scale
			RETURNS:		Returns a MatrixCube with received transRbt and scaleRbt  
		*/

		MatrixCube();
		myTransRbt *= translate;
		myScaleRbt *= scale;

		//Save a copy of the original position
		myOrigObjRbt *= translate * scale;
		myOrigScaleRbt *= scale;
		myOrigTransRbt *= translate;
	}
	/*-----------------------------------------------*/
	MatrixCube::MatrixCube(Matrix4 translate, Matrix4 rotate, Matrix4 scale)
	{
		/* PURPOSE:		Constructs a MatrixCube with given transRbt and scaleRbt 
			RECEIVES:	Matrix4 translate - a transRbt Matrix that has the cube's initial translation 
							Matrix4 scale - a scaleRbt Matrix that has the cube's initial scale
			RETURNS:		Returns a MatrixCube with received transRbt and scaleRbt  
		*/

		MatrixCube();
		myTransRbt *= translate;
		myScaleRbt *= scale;
		myRotateRbt *= rotate;

		//Save a copy of the original position
		myOrigObjRbt *= translate * rotate * scale;
		myOrigScaleRbt *= scale;
		myOrigTransRbt *= translate;
		myOrigRotateRbt *= rotate;
	}
	/*-----------------------------------------------*/
	Matrix4 MatrixCube::getTransRbt()
	{
		/* PURPOSE:		Retrieves translationRbt 
			RETURNS:		the cube's translation Rbt
		*/

		return myTransRbt;
	}
	/*-----------------------------------------------*/
	void MatrixCube::setTransRbt(Matrix4 T)
	{
		/* PURPOSE:		Sets the cube's Translation Rbt with respect to the origin 
			RECEIVES:	Matrix4 T - The cube's new Translation Rbt 
			REMARKS:		If you wish to update the cube with respect to it's parent you will need to take (objRbt * T) as T  
		*/

		myTransRbt = Matrix4() * T;
	}
	/*-----------------------------------------------*/
	void MatrixCube::resetTransRbt()
	{
		/* PURPOSE:		Reset Translation Rbt to move cube back to starting point 
		*/

		myTransRbt = Matrix4() * myOrigTransRbt;
	}
	/*-----------------------------------------------*/
	Matrix4 MatrixCube::getScaleRbt()
	{
		/* PURPOSE:		Retrieves scaleRbt 
			RETURNS:		the cube's scale Rbt
		*/
		return myScaleRbt;
	}
	/*-----------------------------------------------*/
	void MatrixCube::setScaleRbt(Matrix4 S)
	{
		/* PURPOSE:		Sets the cube's Scale Rbt with respect to the origin 
			RECEIVES:	Matrix4 S - The cube's new Scale Rbt 
			REMARKS:		If you wish to update the cube with respect to it's parent you will need to take (objRbt * S) as S
		*/
		myScaleRbt = Matrix4() * S;
	}
	/*-----------------------------------------------*/
	void MatrixCube::resetScaleRbt()
	{
		/* PURPOSE:		Reset Scale Rbt to resize cube back to starting size
		*/

		myScaleRbt = Matrix4() * myOrigScaleRbt;
	}
	/*-----------------------------------------------*/
	Matrix4 MatrixCube::getRotateRbt()
	{
		/* PURPOSE:		Retrieves rotateRbt 
			RETURNS:		the cube's rotation Rbt
		*/
		return myRotateRbt;
	}
	/*-----------------------------------------------*/
	void MatrixCube::setRotateRbt(Matrix4 R)
	{
		/* PURPOSE:		Sets the cube's Rotate Rbt with respect to the origin 
			RECEIVES:	Matrix4 R - The cube's new Rotate Rbt 
			REMARKS:		If you wish to update the cube with respect to it's parent you will need to take (objRbt * R) as R  
		*/
		myRotateRbt = Matrix4() * R;
	}
	/*-----------------------------------------------*/
	void MatrixCube::resetRotateRbt()
	{
		/* PURPOSE:		Reset Rotate Rbt to rotate cube back to starting angle 
		*/
		myRotateRbt = Matrix4() * myOrigRotateRbt;
	}
	/*-----------------------------------------------*/
	Matrix4 MatrixCube::getObjRbt()
	{
		/* PURPOSE:		Retrieves cube's obj Rbt 
			RETURNS:		the cube's obj Rbt
		*/
		return myObjRbt * myTransRbt * myRotateRbt * myScaleRbt;
	}
	/*-----------------------------------------------*/
	void MatrixCube::update(Matrix4 objRbt)
	{
		/* PURPOSE:		Update cube's original position and size 
			RECEIVES:	Matrix4 objRbt - Matrix to adjust cube's objRbt by 
			REMARKS:		Used specifically to handle mouse based adjustments to cube 
		*/

		myObjRbt *= objRbt;
	}
	/*-----------------------------------------------*/
	void MatrixCube::updateRotateRbt(Matrix4 objRbt)
	{
		/* PURPOSE:		Update cube's original position and size 
			RECEIVES:	Matrix4 objRbt - Matrix to adjust cube's objRbt by 
			REMARKS:		Used specifically to handle mouse based adjustments to cube 
		*/

		myRotateRbt = objRbt;
		myOrigRotateRbt = objRbt;

	}
	/*-----------------------------------------------*/
	Matrix4 MatrixCube::updateWithRespectToA(Matrix4 A)
	{
		/* PURPOSE:		Update the cubes objRbt relative to Matrix A 
			RECEIVES:	Matrix4 A - Parent obj that update is relative to 
			RETURNS:		a Matrix4 that is the cube's new objRbt 
		*/
		return A * myObjRbt * myTransRbt * myRotateRbt * myScaleRbt;
	}
	/*-----------------------------------------------*/
	void MatrixCube::originalPosition()
	{
		/* PURPOSE:		Resets Cube's transRbt, scaleRbt, and RotateRbt back to original values 
		*/
		resetScaleRbt();
		resetTransRbt();
		resetRotateRbt();
	}
}