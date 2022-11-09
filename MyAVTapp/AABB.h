#pragma once

class AABB {
public:

	float xmin;
	float xmax;
	float ymin;
	float ymax;

	AABB(float xmin, float xmax, float ymin, float ymax);

	bool intersects(AABB aabb);
};