#include "camera.h"
#include "AVTmathLib.h"

#include <math.h>

Camera::Camera(int type) {

	this->type = type;
	
	if (type == FIXED) {
		camUp[0] = 1; camUp[1] = 0; camUp[2] = 0;
	} else if (type == ORTHO) {
		camUp[0] = 1; camUp[1] = 0; camUp[2] = 0;
	}

}

Camera::Camera(int type, Rover rover) {

	this->type = type;
	camUp[0] = 0; camUp[1] = 1; camUp[2] = 0;

}

void Camera::updatePos(float x, float y, float z) {
	camPos[0] = x;
	camPos[1] = y;
	camPos[2] = z;
}

void Camera::updateTarget(float x, float y, float z) {
	camTarget[0] = x;
	camTarget[1] = y;
	camTarget[2] = z;
}

void Camera::updateFollowPos(Rover rover) {

	float view_angle = 45 * M_PI / 180;

	float dir = rover.getAngle();

	camPos[0] = rover.pos[0] - 5 * sin(view_angle) * cos(dir);
	camPos[1] = 5.5 - 5 * cos(view_angle);
	camPos[2] = rover.pos[2] - 5 * sin(view_angle) * sin(dir);

	this->updateFollowTarget(rover, camYaw, camPitch);

}

void Camera::updateFollowTarget(Rover rover, float yaw, float pitch) {
	
	camYaw = yaw;
	camPitch = pitch;

	yaw *= 3.14f / 180;
	pitch *= 3.14f / 180;

	yaw += rover.getAngle();

	camDirection[0] = 5 * cos(yaw) * cos(pitch);
	camDirection[1] = 5 * sin(pitch);
	camDirection[2] = 5 * sin(yaw) * cos(pitch);
	normalize(camDirection);

	float camRight[3];
	float worldUp[3] = { 0, 1, 0 };
	crossProduct(camDirection, worldUp, camRight);
	normalize(camRight);

	crossProduct(camRight, camDirection, camUp);
	normalize(camUp);

	camTarget[0] = camDirection[0] + camPos[0];
	camTarget[1] = camDirection[1] + camPos[1];
	camTarget[2] = camDirection[2] + camPos[2];

}