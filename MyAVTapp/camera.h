#pragma once

#include "rover.h"

#define FIXED 0
#define ORTHO 1
#define FOLLOW 2

class Camera {
public:
	float camPos[3] = { 0, 0, 0 };
	float camTarget[3] = { 0, 0, 0 };
	float camUp[3] = { 0, 0, 0 };
	float camDirection[3] = { 0, 0, 0 };
	float camYaw = 0;
	float camPitch = 0;
	int type = 0; // 0:fixed, 1:ortho, 2:follow 

public:
	Camera(int type);
	Camera(int type, Rover rover);

	void updatePos(float x, float y, float z);
	void updateTarget(float x, float y, float z);
	void updateFollowPos(Rover rover);
	void updateFollowTarget(Rover rover, float a, float b);
};