#pragma once
#define _USE_MATH_DEFINES


#define TURNING_ANGLE M_PI / 10
#define ANGL_ACCEL .1f

#define MAX_SPEED 4.f
#define ACCEL 0.05f

#define FWD 1
#define RVRS -1
#define RIGHT 1
#define LEFT -1

#define WHEELBASE .7f
#define RADS_TO_DEG 180 /M_PI

#define HEADLAMP
#include <vector>
#include <math.h>
#include "AABB.h"

class Rock;
class Cylinder;

class Rover {
	float angle = 0;
	float heading = 0;

	int fwd = 0;
	int right = 0;

public:
	float pos[3] = { .0f, .0f, .0f };
	AABB aabb = AABB(0, 0, 0, 0);
	float speed = 0;

	int lives = 5;
	bool died = false;

private:
	//void updateWheelPosition(float* ret, bool frontWheels, float* wheelDir, float delta);
	void updateWheelPosition(float* ret, bool frontWheels, float dir, float turnDir, float delta);
	void updatePosition(float* front, float* back);

public:
	Rover(float posX, float posZ, float facesX, float facesZ);

	void accel(int dir);
	void turn(int dir);
	void updateSpeedAndAngle(float delta);

	void stop();
	void straighten();
	void reset();
	void resetPos();
	void checkRockCollision(Rock rock);
	void checkCylinderCollision(Cylinder cylinder);
	void move(float delta);

	float getAngle();
	float getAngleDegs();
	float getTurningAngle();
	float getTurningAngleDegs();

	void setAABB();

};