#include <iostream>
#include "AABB.h"

AABB::AABB(float xmin, float xmax, float ymin, float ymax) {
	AABB::xmin = xmin;
	AABB::xmax = xmax;
	AABB::ymin = ymin;
	AABB::ymax = ymax;
}

bool AABB::intersects(AABB aabb) {

	if ((xmax >= aabb.xmin && xmin <= aabb.xmin) || (aabb.xmax >= xmin && aabb.xmin <= xmin))
		if ((ymax >= aabb.ymin && ymin <= aabb.ymin) || (aabb.ymax >= ymin && aabb.ymin <= ymin))
			return true;

	return false;
}