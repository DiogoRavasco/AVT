#pragma once
#include <vector>
#include <AABB.h>

#define DEAD_TIME 5
#define SPEED_FACTOR 1.03

class Rover;
class Cylinder;
class Rock {
public:

	float radius;
	float pos[3] = { .0f, .0f, .0f };
	float speed = 0;
	float new_speed = 0;
	float direction[3] = { .0f, .0f, .0f };
	float spin = 0;
	bool dead = false;
	int deathTimer = 0;
	AABB aabb = AABB(0,0,0,0);


	Rock(float radius);
	void increaseSpeed();
	void reset();
	void move(float delta);
	void setNewPosition();
	void checkRockCollision(Rock rock);
	void checkRoverCollision(Rover rover);
	void checkCylinderCollision(Cylinder cylinder);
	void setAABB();
};