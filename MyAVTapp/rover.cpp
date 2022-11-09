#define _USE_MATH_DEFINES

#include <algorithm>
#include <thread>
#include <chrono>
#include <math.h>
#include "rover.h"
#include "AVTmathLib.h"
#include "rock.h"
#include "cylinder.h"

Rover::Rover(float posX, float posZ, float facesX, float facesZ) {
	pos[0] = posX;
	pos[2] = posZ;
}

void Rover::accel(int dir) {
	fwd = dir;
}

void Rover::turn(int dir) {
	right = dir;
}

void Rover::updateSpeedAndAngle(float delta) {
	speed = speed + fwd * ACCEL * delta / 10;
	angle = angle + right * ANGL_ACCEL * delta / 100;

	speed = std::max(-MAX_SPEED, std::min(speed, MAX_SPEED));
	angle = std::max((float)-TURNING_ANGLE, std::min(angle, (float)TURNING_ANGLE));
}

void Rover::move(float delta) {
	updateSpeedAndAngle(delta);

	if (speed == 0) return;

	/* Movement is calculated through bicycle model
	*  View car as having one front wheel  and one back wheel, both in the middle
	*  Doesn't take wheel slippage into account
	*  Assumes wheels only move in direction they face
	*/

	float frontWheelPos[3] = {.0f, .0f, .0f};
	updateWheelPosition(frontWheelPos, true, heading, angle, delta);


	float backWheelPos[3] = {.0f, .0f, .0f};
	updateWheelPosition(backWheelPos, false, heading, .0f, delta);

	//updatePosition(frontWheelPos, backWheelPos);
	pos[0] = (frontWheelPos[0] + backWheelPos[0]) / 2;
	pos[2] = (frontWheelPos[2] + backWheelPos[2]) / 2;

	heading = atan2(frontWheelPos[2] - backWheelPos[2], frontWheelPos[0] - backWheelPos[0]);

}

void Rover::updateWheelPosition(float* ret, bool frontWheels, float dir, float turnDir,  float delta) {
	// Calculate the position of the passed wheels before the turn
	ret[0] = pos[0] + (frontWheels ? 1 : -1) * (WHEELBASE / 2) * cos(dir);
	ret[2] = pos[2] + (frontWheels ? 1 : -1) * (WHEELBASE / 2) * sin(dir);

	// Calculate the position of the passed wheels after the turn
	ret[0] = ret[0] + cos(dir + turnDir) * speed * (delta / 1000.f);
	ret[2] = ret[2] + sin(dir + turnDir) * speed * (delta / 1000.f);

}

void Rover::stop() {
	fwd = 0;
	right = 0;

	speed = 0;
	angle = 0;
}

void Rover::straighten() {
	right = 0;

	angle = 0;
}

void Rover::reset() {
	resetPos();
	lives = 5;
}

void Rover::resetPos() {
	stop();
	pos[0] = 0.f;
	pos[2] = 0.f;
	heading = 0;
}

void Rover::checkRockCollision(Rock rock) {
	rock.setAABB();
	setAABB();

	if (aabb.intersects(rock.aabb)) {
		printf("lost a life\n");
		printf("remaining: %d\n", lives);
		// The rover stops
		if (!--lives) {
			reset();
			died = true;
		} else {
			resetPos();
		}
	}
}

void Rover::checkCylinderCollision(Cylinder cylinder) {
	cylinder.setAABB();
	setAABB();

	if (aabb.intersects(cylinder.aabb)) {
		// The rover goes back to the beginning
		stop();

		pos[0] += (pos[0] - cylinder.pos[0]) * 0.01;
		pos[2] += (pos[2] - cylinder.pos[2]) * 0.01;
	}
}


float Rover::getAngle() {
	return heading;
}

float Rover::getAngleDegs() {
	return heading * RADS_TO_DEG;
}

float Rover::getTurningAngle() {
	return angle;
}

float Rover::getTurningAngleDegs() {
	return angle * RADS_TO_DEG;
}

void Rover::setAABB() {
	float c = cos(heading) * 0.5;
	float s = sin(heading) * 0.5;

	float cornerX1 = pos[0] + c - s;
	float cornerX2 = pos[0] + c + s;
	float cornerX3 = pos[0] - c - s;
	float cornerX4 = pos[0] - c + s;

	float cornerY1 = pos[2] + s - c;
	float cornerY2 = pos[2] + s + c;
	float cornerY3 = pos[2] - s - c;
	float cornerY4 = pos[2] - s + c;

	float xmin = std::min(std::min(cornerX1, cornerX2), std::min(cornerX3, cornerX4));
	float xmax = std::max(std::max(cornerX1, cornerX2), std::max(cornerX3, cornerX4));

	float ymin = std::min(std::min(cornerY1, cornerY2), std::min(cornerY3, cornerY4));
	float ymax = std::max(std::max(cornerY1, cornerY2), std::max(cornerY3, cornerY4));

	aabb = AABB(xmin, xmax, ymin, ymax);
}
