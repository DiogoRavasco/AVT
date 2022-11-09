#include "cylinder.h"
#include "AVTmathLib.h"
#include "rock.h"
#include "rover.h"
#include <iostream>



Cylinder::Cylinder(float x, float y) {
	radius = 2;
	height = 5;
	pos[0] = x;
	pos[2] = y;
}


void Cylinder::move(float delta) {

	if (new_speed != 0) {
		speed = new_speed;
		new_speed = 0;
	}

	// Cylinder speed decrements with time --> Tipo atrito
	if (speed < 0.0001) {
		speed = 0;
		direction[0] = 0;
		direction[2] = 0;
	}
	else {
		speed = speed - speed * delta *  0.01;
	}

	pos[0] = pos[0] + direction[0] * speed * delta;
	pos[2] = pos[2] + direction[2] * speed * delta;

}


void Cylinder::checkRockCollision(Rock rock2) {
	setAABB();
	rock2.setAABB();

	if (aabb.intersects(rock2.aabb)) {
		// The cylinder moves in the direction opposite to the rock
		direction[0] = pos[0] - rock2.pos[0];
		direction[2] = pos[2] - rock2.pos[2];
		normalize(direction);

		// The speeds are averaged (The rock speed is fractioned tho) 
		new_speed = (speed + rock2.speed/5) / 2;

	}

}

void Cylinder::checkRoverCollision(Rover rover) {
	setAABB();
	rover.setAABB();

	if (aabb.intersects(rover.aabb)) {
		// The cylinder moves in the direction opposite to the rock
		direction[0] = pos[0] - rover.pos[0];
		direction[2] = pos[2] - rover.pos[2];
		normalize(direction);

		new_speed = speed/2 + rover.speed / 1000;
	}
}

void Cylinder::checkCylinderCollision(Cylinder cylinder2) {
	setAABB();
	cylinder2.setAABB();

	if (aabb.intersects(cylinder2.aabb)) {
		// The cylinder moves in the direction opposite to the rock
		direction[0] = pos[0] - cylinder2.pos[0];
		direction[2] = pos[2] - cylinder2.pos[2];
		normalize(direction);

		// The speeds are averaged
		new_speed = (cylinder2.speed + speed) * 2;

	}
}

void Cylinder::setAABB() {

	float x_min = pos[0] - radius;
	float x_max = pos[0] + radius;
	float y_min = pos[2] - radius;
	float y_max = pos[2] + radius;

	aabb = AABB(x_min, x_max, y_min, y_max);
}


