// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <vector>
namespace Mcgee
{
	struct vector
	{
		double x, y, z;
	};
}

typedef Mcgee::vector point;

class face
{
private:
	bool boundary; //is boundary?'
	std::vector<int> globalPointPointers_;
	std::vector<point>& globalPoints_;
	//Reorders points for right hand rule. 
	checkPointOrder()
	{

	};
public:
	//Constructor
	face(std::vector<point>& points, std::vector<int> globalPointPointers, boundary = false) :
		globalPoints_(points),
		globalPointPointers_(globalPointPointers)
	{
		checkPointOrder()
	};
};

class cell
{
private:
	std::vector<face> faces_; //Faces owned by cell
	uint64_t numFaces; // number of owned faces
public:
	//constructor
	cell(std::vector<face> faces) :
		faces_(faces)
	{};

	cell()
	{
	};

};
int main()
{	
	//Load file 
	std::vector<point> globalPoints;

	//getPoints(globalPoints);

	//Loop all faces to consturct
	for (int i = 0; i < AllFaces; i++)
	{
		std::vector<point> localPoints;
		std::vector<uint64_t> globalPointPointers;
		localPoints.resize(4);
		globalPointPointers.resize(4);
		globalPointPointers[0] = ? ;
		globalPointPointers[1] = ? ;
		globalPointPointers[2] = ? ;
		globalPointPointers[3] = ? ;

		face(globalPoints, globalPointPointers);
	}

	//create face from 4 points

	// would need to follow right hand rule on face
	// 
}

// 0 1 2
// p1 p1
// p2 p5
// p3  p6
// p4 p2

// 
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
