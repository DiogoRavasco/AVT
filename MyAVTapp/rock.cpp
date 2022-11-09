#include "rock.h"
#include "AVTmathLib.h"
#include "rover.h"
#include "cylinder.h"

Rock::Rock(float r) {
	setNewPosition();
	speed = ((rand() % 5) + 1) * 0.001;
	radius = r;
} 

void Rock::increaseSpeed() {
	new_speed = speed * SPEED_FACTOR;
}

void Rock::reset() {
	speed = ((rand() % 5) + 1) * 0.001;
	setNewPosition();
}

void Rock::move(float delta) {

	if (dead) {
		deathTimer += delta;
		if (deathTimer > DEAD_TIME * 1000) {
			dead = false;
			deathTimer = 0;
			setNewPosition();
		}
	}

	else {
		if (new_speed != 0) {
			speed = new_speed;
			new_speed = 0;
		}
		pos[0] = pos[0] + direction[0] * speed * delta;
		pos[2] = pos[2] + direction[2] * speed * delta;

		float spinRadians = (speed * delta) / radius;
		spin += (spinRadians * 180 / 3.14159265358);

		float dist = length(pos);

		if (dist > 30) {
			setNewPosition();
			dead = true;
		}

	}

	//std::cout << "rock position: " << pos[0] << " " << pos[2] << std::endl;

}

void Rock::setNewPosition() {

	if ((rand() % 2) == 0) {
		pos[0] = rand() % 81 + (-40);
		pos[2] = (rand() % 2) ? 40.f : -40.f;
	}
	else {
		pos[0] = (rand() % 2) ? 40.f : -40.f;
		pos[2] = rand() % 81 + (-40);
	}

	normalize(pos);

	constProduct(30.f, pos, pos);



	direction[0] = -pos[0] + rand() % 11 + -(5);
	direction[2] = -pos[2] + rand() % 11 + -(5);
	normalize(direction);

}

void Rock::checkRockCollision(Rock rock2) {
	setAABB();
	rock2.setAABB();

	if (aabb.intersects(rock2.aabb)) {
		// It goes in the opposite direction to the other rock
		direction[0] = pos[0] - rock2.pos[0];
		direction[2] = pos[2] - rock2.pos[2];
		normalize(direction);

		// The speeds are averaged
		new_speed = (speed + rock2.speed) / 2;


	}

}

void Rock::checkRoverCollision(Rover rover) {
	setAABB();
	rover.setAABB();

	if (aabb.intersects(rover.aabb)) {
		// It goes in the opposite direction to the rover rock
		direction[0] = pos[0] - rover.pos[0];
		direction[2] = pos[2] - rover.pos[2];
		normalize(direction);

	}
}

void Rock::checkCylinderCollision(Cylinder cylinder) {
	setAABB();
	cylinder.setAABB();

	if (aabb.intersects(cylinder.aabb)) {
		// It goes in the opposite direction of the cylinder
		direction[0] = pos[0] - cylinder.pos[0];
		direction[2] = pos[2] - cylinder.pos[2];
		normalize(direction);

		// The speed is the same
		new_speed = speed;

	}

}

void Rock::setAABB() {

	float x_min = pos[0] - radius;
	float x_max = pos[0] + radius;
	float y_min = pos[2] - radius;
	float y_max = pos[2] + radius;

	aabb = AABB(x_min, x_max, y_min, y_max);
}


