#pragma once

void pushModel(GLint& loc, int objId);

void popModel(int objId);

void renderBillboards(GLint& loc, int objId);
void renderFireworks(GLint& loc, int objId);
void renderSkyBox(GLint& loc, int objId);