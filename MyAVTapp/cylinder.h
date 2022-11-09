#pragma once
#include "AABB.h"

class Rock;
class Rover;
class Cylinder {
public:

	float radius;
	float height;
	float pos[3] = { .0f, .0f, .0f };
	float speed = 0;
	float new_speed = 0;
	float direction[3] = { .0f, .0f, .0f };
	
	AABB aabb = AABB(0, 0, 0, 0);


	Cylinder(float x, float y);
	void move(float delta);
	void checkRockCollision(Rock rock);
	void checkRoverCollision(Rover rover);
	void checkCylinderCollision(Cylinder cylinder2);
	void setAABB();
};