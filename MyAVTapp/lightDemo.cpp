//
// AVT: Phong Shading and Text rendered with FreeType library
// The text rendering was based on https://learnopengl.com/In-Practice/Text-Rendering
// This demo was built for learning purposes only.
// Some code could be severely optimised, but I tried to
// keep as simple and clear as possible.
//
// The code comes with no warranties, use it at your own risk.
// You may use it, or parts of it, wherever you want.
// 
// Author: João Madeiras Pereira
//

#include <math.h>
#include <iostream>
#include <sstream>
#include <string>

// include GLEW to access OpenGL 3.3 functions
#include <GL/glew.h>


// GLUT is the toolkit to interface with the OS
#include <GL/freeglut.h>

#include <IL/il.h>


// Use Very Simple Libs
#include "VSShaderlib.h"
#include "AVTmathLib.h"
#include "VertexAttrDef.h"
#include "geometry.h"
#include "rover.h"
#include "camera.h"
#include "lights.h"
#include "rock.h"
#include "cylinder.h"
#include "Texture_Loader.h"
#include "l3dBillboard.h"
#include "flare.h"

#include "avtFreeType.h"
#include "lightDemo.h"

#include "assimp/Importer.hpp"
#include "assimp/scene.h"
#include "meshFromAssimp.h"

using namespace std;

#define CAPTION "AVT Demo: Phong Shading and Text rendered with FreeType"

#define frand()			((float)rand()/RAND_MAX)
#define MAX_PARTICULAS  6000


inline double clamp(const double x, const double min, const double max) {
	return (x < min ? min : (x > max ? max : x));
}

inline int clampi(const int x, const int min, const int max) {
	return (x < min ? min : (x > max ? max : x));
}

int deltaMove = 0, deltaUp = 0, type = 0;
int fireworks = 0;

int counter = 0;

typedef struct {
	float	life;		// vida
	float	fade;		// fade
	float	r, g, b;    // color
	GLfloat x, y, z;    // posi‹o
	GLfloat vx, vy, vz; // velocidade 
	GLfloat ax, ay, az; // acelera‹o
} Particle;

Particle particula[MAX_PARTICULAS];
int dead_num_particles = 0;

int WindowHandle = 0;
int WinX = 1024, WinY = 768;
int quadSize = 80; // defines the size of the platform

unsigned int FrameCount = 0;
unsigned int oldTime = glutGet(GLUT_ELAPSED_TIME);
bool fog = false;
bool directionalOn = true;
bool pointsOn = true;
bool spotsOn = true;
bool rearView = false;
bool flareEffect = true;

GLuint FlareTextureArray[5];

//Flare effect
FLARE_DEF AVTflare;
float lightScreenPos[3];  //Position of the light in Window Coordinates


int rockCount = 0;

// 0: running, 1: paused, 2: reset, 3: game over
int paused = 0;

int bumpmap1 = 0;

int final_score = 0;
int score = 0;

// assimp
extern Assimp::Importer importer;
extern const aiScene* scene;
char model_dir[50];
extern float scaleFactor;

bool normalMapKey = TRUE;

//shaders
VSShaderLib shader;  //geometry
VSShaderLib shaderText;  //render bitmap text

//File with the font
const string font_name = "fonts/arial.ttf";

//Vector with meshes
vector<struct MyMesh> myMeshes;
//Vector with Rover assimp mesh
vector<struct MyMesh> roverMesh;

Rover rover = Rover(.0f, .0f, 1.f, .0f);
Light lights[9];

vector<Rock> rockList;
#define N_ROCKS 10
#define ROCK_RADIUS 1.f

vector<Cylinder> cylinderList;
#define N_CYLINDERS 4

/// The storage for matrices
extern float mMatrix[COUNT_MATRICES][16];
extern float mCompMatrix[COUNT_COMPUTED_MATRICES][16];

/// The normal matrix
extern float mNormal3x3[9];

GLint pvm_uniformId;
GLint vm_uniformId;
GLint normal_uniformId;
GLint model_uniformId;
GLint view_uniformId;
GLint directional_uniformId;
GLint point1_uniformId;
GLint point2_uniformId;
GLint point3_uniformId;
GLint point4_uniformId;
GLint point5_uniformId;
GLint point6_uniformId;
GLint left_uniformId;
GLint right_uniformId;
GLint lampDir_uniformId;
GLint tex_loc, tex_loc1, tex_loc2, tex_loc3, tex_loc4, tex_cube_loc, tex_normalMap_loc;
GLint texMode_uniformId;
GLint reflect_perFragment_uniformId;
GLint shadowMode_uniformId;
GLint normalMap_loc;
GLint specularMap_loc;
GLint diffMapCount_loc;
GLint bumpmap1_loc;

GLuint TextureArray[7];


// Cameras
Camera cams[3] = { Camera(FIXED), Camera(ORTHO), Camera(FOLLOW, rover) };
int activeCam = FOLLOW; // 0: fixed, 1: top, 2: follow

// Mouse Tracking Variables
int startX, startY, tracking = 0;

// Camera Spherical Coordinates
float alpha = .0f, beta = .0f;
float r = quadSize;

// Frame counting and FPS computation
long myTime, timebase = 0, frame = 0;
char s[32];
float lightPos[4] = { 10.0f, 6.0f, .0f, 0.0f };


void updateParticles()
{
	int i;
	float h;

	/* Método de Euler de integração de eq. diferenciais ordinárias
	h representa o step de tempo; dv/dt = a; dx/dt = v; e conhecem-se os valores iniciais de x e v */

	//h = 0.125f;
	h = 0.033;
	if (fireworks) {

		for (i = 0; i < MAX_PARTICULAS; i++)
		{
			particula[i].x += (h * particula[i].vx);
			particula[i].y += (h * particula[i].vy);
			particula[i].z += (h * particula[i].vz);
			particula[i].vx += (h * particula[i].ax);
			particula[i].vy += (h * particula[i].ay);
			particula[i].vz += (h * particula[i].az);
			particula[i].life -= particula[i].fade;
		}
	}
}


void iniParticles(void)
{
	GLfloat v, theta, phi;
	int i;

	for (i = 0; i < MAX_PARTICULAS; i++)
	{
		v = 0.8 * frand() + 0.2;
		phi = frand() * M_PI;
		theta = 6.0 * frand() * M_PI;

		if (i / (MAX_PARTICULAS / 4) == 0) {
			particula[i].x = 40.0f;
			particula[i].z = 0.0f;
		}
		else if (i / (MAX_PARTICULAS / 4) == 1) {
			particula[i].x = -40.0f;
			particula[i].z = 0.0f;
		}
		else if (i / (MAX_PARTICULAS / 4) == 2) {
			particula[i].x = 0.0f;
			particula[i].z = 40.0f;
		}
		else if (i / (MAX_PARTICULAS / 4) == 3) {
			particula[i].x = 0.0f;
			particula[i].z = -40.0f;
		}

		particula[i].y = 0.0f;
		particula[i].vx = v * cos(theta) * sin(phi);
		particula[i].vy = v * cos(phi);
		particula[i].vz = v * sin(theta) * sin(phi);
		particula[i].ax = 0.0f; /* simular um pouco de vento */
		particula[i].ay = 0.5f; /* simular a aceleração da gravidade */
		particula[i].az = 0.0f;

		/* tom amarelado que vai ser multiplicado pela textura que varia entre branco e preto */
		particula[i].r = 0.882f;
		particula[i].g = 0.552f;
		particula[i].b = 0.211f;

		particula[i].life = 1.0f;		/* vida inicial */
		particula[i].fade = 0.003f;	    /* step de decréscimo da vida para cada iteração */
	}
}

void timer(int value)
{
	std::ostringstream oss;
	oss << CAPTION << ": " << FrameCount << " FPS @ (" << WinX << "x" << WinY << ")";
	std::string s = oss.str();
	glutSetWindow(WindowHandle);
	glutSetWindowTitle(s.c_str());
	FrameCount = 0;
	if (++rockCount == 5) {
		rockCount = 0;
		for (int i = 0; i < N_ROCKS; i++) {
			rockList[i].increaseSpeed();
		}
	}
	if (!paused) score++;

	if (score % 10 == 0) {
		fireworks = 1;
		iniParticles();
	}

	glutTimerFunc(1000, timer, 0);
}

void pause() {
	paused = 1;
}

void unpause() {
	paused = 0;
	oldTime = glutGet(GLUT_ELAPSED_TIME);
}

void reset_game(int status) {
	paused = status; // 2: reset, 3: game over
	rover.reset();
	for (int i = 0; i < N_ROCKS; i++) rockList[i].reset();
	rover.died = false;
	final_score = score;
	score = 0;
}

void refresh(int value)
{

	if (rover.died) reset_game(3);

	if (paused) {
		glutTimerFunc(0, refresh, 0);
		return;
	}

	rover.move(glutGet(GLUT_ELAPSED_TIME) - oldTime);
	for (int i = 0; i < N_ROCKS; i++)
		rockList[i].move(glutGet(GLUT_ELAPSED_TIME) - oldTime);
	for (int i = 0; i < N_CYLINDERS; i++)
		cylinderList[i].move(glutGet(GLUT_ELAPSED_TIME) - oldTime);

	// Check collisions between rocks
	for (int i = 0; i < N_ROCKS; i++) {
		for (int j = 0; j < N_ROCKS; j++) {
			if (i == j) continue;
			rockList[i].checkRockCollision(rockList[j]);
		}
	}

	// Check rover-rock collision
	for (int i = 0; i < N_ROCKS; i++) {
		rover.checkRockCollision(rockList[i]);
		rockList[i].checkRoverCollision(rover);
	}

	// Check cylinder collisions
	for (int i = 0; i < N_CYLINDERS; i++) {
		cylinderList[i].checkRoverCollision(rover);
		rover.checkCylinderCollision(cylinderList[i]);
		for (int j = 0; j < N_ROCKS; j++) {
			cylinderList[i].checkRockCollision(rockList[j]);
			rockList[j].checkCylinderCollision(cylinderList[i]);
		}

		// cylinder-cylinder
		for (int k = 0; k < N_CYLINDERS; k++) {
			if (i == k) continue;
			cylinderList[i].checkCylinderCollision(cylinderList[k]);
		}
	}

	oldTime = glutGet(GLUT_ELAPSED_TIME);

	glutTimerFunc(0, refresh, 0);
}

void createLights() {
	// Create directional light
	Light directional = {
		{10.f, 5.f, 10.f, .0f}
	};

	lights[0] = directional;

	// Create point lights
	for (int i = 0; i < 6; i++) {
		lights[i + 1].src[0] =  ((i % 2 == 0) ? -1 : 1) * 10.f;
		lights[i + 1].src[1] = 2;
		lights[i + 1].src[2] = ((i % 2 == 0) ? i : (i - 1)) * 10 - 20;
		lights[i + 1].src[3] = 1.f; 
	}

	// Create spot lights
	Light leftLamp = {
		{.51f, .4f, .35f, 1.f}
	};
	
	Light rightLamp = {
		{.51f, .4f, -.35f, 1.f}
	};

	lights[7] = leftLamp;
	lights[8] = rightLamp;

}

// ------------------------------------------------------------
//
// Reshape Callback Function
//

void createFloorStencil(GLenum func) {
	GLint loc;

	glUseProgram(shader.getProgramIndex());

	glStencilFunc(GL_NEVER, 0x0, 0x1);
	glStencilOp(func, GL_KEEP, GL_KEEP);

	pushModel(loc, 0);
	rotate(MODEL, 90, -1, 0, 0);

	glUniform1i(texMode_uniformId, 2);
	popModel(0);
}

void createStencil(int w, int h, GLint ref) {
	/* create a diamond shaped stencil area */
	loadIdentity(PROJECTION);
	if (w <= h)
		ortho(-2.0, 2.0, -2.0 * (GLfloat)h / (GLfloat)w,
			2.0 * (GLfloat)h / (GLfloat)w, -10, 10);
	else
		ortho(-2.0 * (GLfloat)w / (GLfloat)h,
			2.0 * (GLfloat)w / (GLfloat)h, -2.0, 2.0, -10, 10);

	// load identity matrices for Model-View
	loadIdentity(VIEW);
	loadIdentity(MODEL);

	glUseProgram(shader.getProgramIndex());

	//não vai ser preciso enviar o material pois o cubo não é desenhado
	translate(MODEL, -1.f, .8f, -0.5f);
	scale(MODEL, 2.f, 1.f, 2.f);
	// send matrices to OGL
	computeDerivedMatrix(PROJ_VIEW_MODEL);
	//glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
	glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
	computeNormalMatrix3x3();
	glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

	glStencilFunc(GL_NEVER, 0x0, 0x1);
	glStencilOp(GL_REPLACE, GL_KEEP, GL_KEEP);

	glBindVertexArray(myMeshes[9].vao);
	glDrawElements(myMeshes[9].type, myMeshes[9].numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}

void changeSize(int w, int h) {

	float ratio;
	// Prevent a divide by zero, when window is too short
	if (h == 0)
		h = 1;
	// set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// set the projection matrix
	ratio = (1.0f * w) / h;
	WinX = w;
	WinY = h;
	loadIdentity(PROJECTION);
	perspective(53.13f, ratio, 0.1f, 1000.0f);
}


// ------------------------------------------------------------
//
// Render stufff
//
// Recursive render of the Assimp Scene Graph

void aiRecursive_render(const aiScene* sc, const aiNode* nd)
{
	GLint loc;

	// Get node transformation matrix
	aiMatrix4x4 m = nd->mTransformation;
	// OpenGL matrices are column major
	m.Transpose();

	// save model matrix and apply node transformation
	pushMatrix(MODEL);

	float aux[16];
	memcpy(aux, &m, sizeof(float) * 16);
	multMatrix(MODEL, aux);


	// draw all meshes assigned to this node
	for (unsigned int n = 0; n < nd->mNumMeshes; ++n) {


		// send the material
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.ambient");
		glUniform4fv(loc, 1, roverMesh[nd->mMeshes[n]].mat.ambient);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
		glUniform4fv(loc, 1, roverMesh[nd->mMeshes[n]].mat.diffuse);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.specular");
		glUniform4fv(loc, 1, roverMesh[nd->mMeshes[n]].mat.specular);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.emissive");
		glUniform4fv(loc, 1, roverMesh[nd->mMeshes[n]].mat.emissive);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.shininess");
		glUniform1f(loc, roverMesh[nd->mMeshes[n]].mat.shininess);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.texCount");
		glUniform1i(loc, roverMesh[nd->mMeshes[n]].mat.texCount);

		unsigned int  diffMapCount = 0;  //read 2 diffuse textures

		//devido ao fragment shader suporta 2 texturas difusas simultaneas, 1 especular e 1 normal map

		glUniform1i(normalMap_loc, false);   //GLSL normalMap variable initialized to 0
		glUniform1i(specularMap_loc, false);
		glUniform1ui(diffMapCount_loc, 0);

		if (roverMesh[nd->mMeshes[n]].mat.texCount != 0)
			for (unsigned int i = 0; i < roverMesh[nd->mMeshes[n]].mat.texCount; ++i) {
				if (roverMesh[nd->mMeshes[n]].texTypes[i] == DIFFUSE) {
					if (diffMapCount == 0) {
						diffMapCount++;
						loc = glGetUniformLocation(shader.getProgramIndex(), "texUnitDiff");
						glUniform1i(loc, roverMesh[nd->mMeshes[n]].texUnits[i]);
						glUniform1ui(diffMapCount_loc, diffMapCount);
					}
					else if (diffMapCount == 1) {
						diffMapCount++;
						loc = glGetUniformLocation(shader.getProgramIndex(), "texUnitDiff1");
						glUniform1i(loc, roverMesh[nd->mMeshes[n]].texUnits[i]);
						glUniform1ui(diffMapCount_loc, diffMapCount);
					}
					else printf("Only supports a Material with a maximum of 2 diffuse textures\n");
				}
				else if (myMeshes[nd->mMeshes[n]].texTypes[i] == SPECULAR) {
					loc = glGetUniformLocation(shader.getProgramIndex(), "texUnitSpec");
					glUniform1i(loc, roverMesh[nd->mMeshes[n]].texUnits[i]);
					glUniform1i(specularMap_loc, true);
				}
				else if (myMeshes[nd->mMeshes[n]].texTypes[i] == NORMALS) { //Normal map
					loc = glGetUniformLocation(shader.getProgramIndex(), "texUnitNormalMap");
					if (normalMapKey)
						glUniform1i(normalMap_loc, normalMapKey);
					glUniform1i(loc, roverMesh[nd->mMeshes[n]].texUnits[i]);

				}
				else printf("Texture Map not supported\n");
			}

		// send matrices to OGL
		computeDerivedMatrix(PROJ_VIEW_MODEL);
		glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
		glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
		computeNormalMatrix3x3();
		glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

		// bind VAO
		glBindVertexArray(roverMesh[nd->mMeshes[n]].vao);

		if (!shader.isProgramValid()) {
			printf("Program Not Valid!\n");
			exit(1);
		}
		// draw
		glDrawElements(roverMesh[nd->mMeshes[n]].type, roverMesh[nd->mMeshes[n]].numIndexes, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}

	// draw all children
	for (unsigned int n = 0; n < nd->mNumChildren; ++n) {
		aiRecursive_render(sc, nd->mChildren[n]);
	}
	popMatrix(MODEL);
}

void render_flare(FLARE_DEF* flare, int lx, int ly, int* m_viewport) {  //lx, ly represent the projected position of light on viewport

	int     dx, dy;          // Screen coordinates of "destination"
	int     px, py;          // Screen coordinates of flare element
	int		cx, cy;
	float    maxflaredist, flaredist, flaremaxsize, flarescale, scaleDistance;
	int     width, height, alpha;    // Piece parameters;
	int     i;
	float	diffuse[4];

	GLint loc;

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	int screenMaxCoordX = m_viewport[0] + m_viewport[2] - 1;
	int screenMaxCoordY = m_viewport[1] + m_viewport[3] - 1;

	//viewport center
	cx = m_viewport[0] + (int)(0.5f * (float)m_viewport[2]) - 1;
	cy = m_viewport[1] + (int)(0.5f * (float)m_viewport[3]) - 1;

	// Compute how far off-center the flare source is.
	maxflaredist = sqrt(cx * cx + cy * cy);
	flaredist = sqrt((lx - cx) * (lx - cx) + (ly - cy) * (ly - cy));
	scaleDistance = (maxflaredist - flaredist) / maxflaredist;
	flaremaxsize = (int)(m_viewport[2] * flare->fMaxSize);
	flarescale = (int)(m_viewport[2] * flare->fScale);

	// Destination is opposite side of centre from source
	dx = clampi(cx + (cx - lx), m_viewport[0], screenMaxCoordX);
	dy = clampi(cy + (cy - ly), m_viewport[1], screenMaxCoordY);

	// Render each element. To be used Texture Unit 0

	glUniform1i(texMode_uniformId, 8); // draw modulated textured particles 
	glUniform1i(tex_loc, 0);  //use TU 0

	for (i = 0; i < flare->nPieces; ++i)
	{
		// Position is interpolated along line between start and destination.
		px = (int)((1.0f - flare->element[i].fDistance) * lx + flare->element[i].fDistance * dx);
		py = (int)((1.0f - flare->element[i].fDistance) * ly + flare->element[i].fDistance * dy);
		px = clampi(px, m_viewport[0], screenMaxCoordX);
		py = clampi(py, m_viewport[1], screenMaxCoordY);

		// Piece size are 0 to 1; flare size is proportion of screen width; scale by flaredist/maxflaredist.
		width = (int)(scaleDistance * flarescale * flare->element[i].fSize);

		// Width gets clamped, to allows the off-axis flaresto keep a good size without letting the elements get big when centered.
		if (width > flaremaxsize)  width = flaremaxsize;

		height = (int)((float)m_viewport[3] / (float)m_viewport[2] * (float)width);
		memcpy(diffuse, flare->element[i].matDiffuse, 4 * sizeof(float));
		diffuse[3] *= scaleDistance;   //scale the alpha channel

		if (width > 1)
		{
			// send the material - diffuse color modulated with texture
			loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
			glUniform4fv(loc, 1, diffuse);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, FlareTextureArray[flare->element[i].textureId]);
			pushMatrix(MODEL);
			translate(MODEL, (float)(px - width * 0.0f), (float)(py - height * 0.0f), 0.0f);
			scale(MODEL, (float)width, (float)height, 1);
			computeDerivedMatrix(PROJ_VIEW_MODEL);
			glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
			glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
			computeNormalMatrix3x3();
			glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

			glBindVertexArray(myMeshes[11].vao);
			glDrawElements(myMeshes[11].type, myMeshes[11].numIndexes, GL_UNSIGNED_INT, 0);
			glBindVertexArray(0);
			popMatrix(MODEL);
		}
	}
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glDisable(GL_BLEND);
}

void loadLights() {

	loadIdentity(MODEL);

	for (int i = 0; i < 7; i++) {
		multMatrixPoint(VIEW, lights[i].src, lights[i].srcEye);
	}

	float res[4];
	float res2[4];

	float model[4];
	float model2[4];

	translate(MODEL, rover.pos[0], 0.f, rover.pos[2]);
	rotate(MODEL, -rover.getAngleDegs(), 0, 1, 0);

	multMatrixPoint(MODEL, lights[7].src, model);
	multMatrixPoint(MODEL, lights[8].src, model2);

	multMatrixPoint(VIEW, model, res);
	multMatrixPoint(VIEW, model2, res2);


	float headingDir[4] = { 1.f, .0f, .0f, .0f };


	float headingDirModel[4];
	float headingDirEye[4];


	loadIdentity(MODEL);

	rotate(MODEL, -rover.getAngleDegs(), 0, 1, 0);

	multMatrixPoint(MODEL, headingDir, headingDirModel);



	multMatrixPoint(VIEW, headingDirModel, headingDirEye);


	glUniform4fv(directional_uniformId, 1, lights[0].srcEye);
	glUniform4fv(point1_uniformId, 1, lights[1].srcEye);
	glUniform4fv(point2_uniformId, 1, lights[2].srcEye);
	glUniform4fv(point3_uniformId, 1, lights[3].srcEye);
	glUniform4fv(point4_uniformId, 1, lights[4].srcEye);
	glUniform4fv(point5_uniformId, 1, lights[5].srcEye);
	glUniform4fv(point6_uniformId, 1, lights[6].srcEye);
	glUniform4fv(left_uniformId, 1, res);
	glUniform4fv(right_uniformId, 1, res2);
	glUniform4fv(lampDir_uniformId, 1, headingDirEye);
}

void drawFloor(bool blend) {
	GLint loc;

	pushModel(loc, 0);
	if (!blend) glDisable(GL_BLEND);
	if (activeCam == FOLLOW && blend) {
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}	
	rotate(MODEL, 90, -1, 0, 0);
	glUniform1i(texMode_uniformId, 2);
	popModel(0);
}

void drawAllSceneElements(bool specialFx) {
	GLint loc;

	for (int objId = 1; objId < 10; objId++) {


		// Rover
		if (objId == 1) {

			// Rover Mesh (1178 meshes)
			pushMatrix(MODEL);
			// rover movement
			translate(MODEL, rover.pos[0], .0f, rover.pos[2]);
			rotate(MODEL, -rover.getAngleDegs(), 0, 1, 0);
			rotate(MODEL, 90, 0, 1, 0);
			float scalingFactor = 0.0065f;
			scale(MODEL, scalingFactor, scalingFactor, scalingFactor);

			glUniform1i(texMode_uniformId, 7);
			aiRecursive_render(scene, scene->mRootNode);
			popMatrix(MODEL);
		}

		// Rocks
		if (objId == 3) {
			for (int i = 0; i < N_ROCKS; i++) {
				pushModel(loc, objId);
				Rock rock = rockList[i];
				if (!rock.dead) {
					translate(MODEL, rock.pos[0], 1.f, rock.pos[2]);
					rotate(MODEL, rock.spin, rock.direction[2], 0, -rock.direction[0]);
				}
				else {
					translate(MODEL, 1000.f, -10.f, 1000.f);
				}
				glUniform1i(texMode_uniformId, 1);
				popModel(objId);

			}
			continue;
		}


		// Cylinders
		if (objId == 4) {
			for (int i = 0; i < N_CYLINDERS; i++) {
				pushModel(loc, objId);
				Cylinder cylinder = cylinderList[i];
				translate(MODEL, cylinder.pos[0], 2.5f, cylinder.pos[2]);
				glUniform1i(texMode_uniformId, 0);
				popModel(objId);
			}
		}

		if (specialFx) {
			// Billboards
			if (objId == 6) {
				renderBillboards(loc, objId);
			}


			// Render skybox
			if (objId == 7) {
				renderSkyBox(loc, objId);
			}

			// Fireworks
			if (objId == 8) {
				renderFireworks(loc, objId);
			}


        }

        // Environment Mapping Cube
        if (objId == 9) {
            pushModel(loc, objId);
            glEnable(GL_BLEND);
            translate(MODEL, -2.f, 1.5f, -0.5f);
            glUniform1i(reflect_perFragment_uniformId, 0);
            glUniform1i(texMode_uniformId, 7);
            popModel(objId);
        }

	}
}

void actualRendering(bool rearView) {
	GLint loc;
	// use our shader

	glUseProgram(shader.getProgramIndex());

	//send the light position in eye coordinates
	//glUniform4fv(lPos_uniformId, 1, lightPos); //efeito capacete do mineiro, ou seja lighPos foi definido em eye coord 
	
	
	

	loc = glGetUniformLocation(shader.getProgramIndex(), "fog");
	glUniform1f(loc, (activeCam == FOLLOW) ? fog : false);
	loc = glGetUniformLocation(shader.getProgramIndex(), "directionalOn");
	glUniform1f(loc, directionalOn);
	loc = glGetUniformLocation(shader.getProgramIndex(), "pointsOn");
	glUniform1f(loc, pointsOn);
	loc = glGetUniformLocation(shader.getProgramIndex(), "spotsOn");
	glUniform1f(loc, spotsOn);

	loadIdentity(MODEL);

	//Associar os Texture Units aos Objects Texture
	//stone.tga loaded in TU0; checker.tga loaded in TU1;  lightwood.tga loaded in TU2

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, TextureArray[0]);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, TextureArray[1]);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, TextureArray[2]);

	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D, TextureArray[3]);

	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_2D, TextureArray[4]);

	glActiveTexture(GL_TEXTURE5);
	glBindTexture(GL_TEXTURE_CUBE_MAP, TextureArray[5]);

	//Indicar aos tres samplers do GLSL quais os Texture Units a serem usados
	glUniform1i(tex_loc, 0);
	glUniform1i(tex_normalMap_loc, 1);
	glUniform1i(tex_loc2, 2);
	glUniform1i(tex_loc3, 3);
	glUniform1i(tex_loc4, 4);
	glUniform1i(tex_cube_loc, 5);

	float mat[16];
	GLfloat floor[4] = {0, 1, 0, 0};

	glEnable(GL_DEPTH_TEST);

	if (activeCam == FOLLOW) {
		if (rearView) {

			glStencilFunc(GL_EQUAL, 0x0, 0x1);
			glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
			loadLights();
			drawFloor(false);

		}
		else {
			glClear(GL_STENCIL_BUFFER_BIT);
			createFloorStencil(GL_INCR);
			glStencilFunc(GL_EQUAL, 0x2, 0x1);
			glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

			for (int i = 0; i < 7; i++) {
				lights[i].src[1] *= -1;
			}

			loadLights();
			loadIdentity(MODEL);
			pushMatrix(MODEL);
			scale(MODEL, 1.0f, -1.0f, 1.0f);
			glCullFace(GL_FRONT);
			drawAllSceneElements(true);
			glCullFace(GL_BACK);
			popMatrix(MODEL);


			for (int i = 0; i < 7; i++) {
				lights[i].src[1] *= -1;
			}

			glClear(GL_STENCIL_BUFFER_BIT);
			glStencilFunc(GL_EQUAL, 0x1, 0x1);

			drawFloor(true);

			float res[4];

			glUniform1i(shadowMode_uniformId, 1);
			constProduct(1, lights[0].src, res);

			shadow_matrix(mat, floor, res);
			glDisable(GL_DEPTH_TEST); //To force the shadow geometry to be rendered even if behind the floor

			//Dark the color stored in color buffer
			glEnable(GL_BLEND);
			glBlendFunc(GL_DST_COLOR, GL_ZERO);

			loadIdentity(MODEL);
			pushMatrix(MODEL);
			multMatrix(MODEL, mat);
			drawAllSceneElements(false);
			popMatrix(MODEL);

			glDisable(GL_BLEND);
			glEnable(GL_DEPTH_TEST);

			//render the geometry
			glUniform1i(shadowMode_uniformId, 0);

		}
	}
	else {
		drawFloor(false);
	}

	loadLights();
	loadIdentity(MODEL);
	drawAllSceneElements(true);

	
}

void renderScene(void) {
	// rover

	rearView = !rearView;
	createStencil(WinX, WinY, 0x0);

	if (activeCam == FOLLOW) {
		cams[activeCam].updateFollowPos(rover);
		glEnable(GL_STENCIL_TEST);
	}
	else {
		glDisable(GL_STENCIL_TEST);
	}

	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

	/* ------------- */


	FrameCount++;
	glClear(GL_DEPTH_BUFFER_BIT);
	// load identity matrices
	int vport[4];
	glGetIntegerv(GL_VIEWPORT, vport);
	loadIdentity(PROJECTION);

	loadIdentity(VIEW);
	loadIdentity(MODEL);

	if (rearView && activeCam == FOLLOW) {
		float rearDir[4] = { -1.f, .0f, .0f, .0f };
		float up[4] = { .0f, 1.f, .0f, .0f };

		float rearDirModel[4];
		float upModel[4];

		rotate(MODEL, -rover.getAngleDegs(), 0, 1, 0);
		rotate(MODEL, 45, 0, 0, 1);

		multMatrixPoint(MODEL, up, upModel);
		multMatrixPoint(MODEL, rearDir, rearDirModel);


		perspective(123.13f, (1.0f * WinX) / WinY, .1f, 1000.0f);

		lookAt(rover.pos[0], rover.pos[1] + 1.6f, rover.pos[2],
			rover.pos[0] + rearDirModel[0], rover.pos[1] + rearDirModel[1] + 1.6f, rover.pos[2] + rearDirModel[2],
			upModel[0], upModel[1], upModel[2]);


		glStencilFunc(GL_NOTEQUAL, 0x1, 0x1);

		actualRendering(true);
	}
	else {
		glClear(GL_COLOR_BUFFER_BIT);
	}

	if (activeCam == FIXED || activeCam == FOLLOW) perspective(53.13f, (1.0f * WinX) / WinY, (activeCam == FIXED ? 70.f : .1f), 1000.0f);
	if (activeCam == ORTHO) ortho(-(quadSize / 2) * vport[2] / vport[3], (quadSize / 2) * vport[2] / vport[3], -(quadSize / 2) * vport[3] / vport[3], (quadSize / 2) * vport[3] / vport[3], 70, 100);


	// set the camera position based on the active camera
	lookAt(cams[activeCam].camPos[0], cams[activeCam].camPos[1], cams[activeCam].camPos[2],
		cams[activeCam].camTarget[0], cams[activeCam].camTarget[1], cams[activeCam].camTarget[2],
		cams[activeCam].camUp[0], cams[activeCam].camUp[1], cams[activeCam].camUp[2]);



	if (activeCam == FOLLOW) {
		glStencilFunc(GL_NOTEQUAL, 0x0, 0x1);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
	}
	else {
		glStencilFunc(GL_ALWAYS, 0x1, 0x1);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
	}

	actualRendering(false);

	if (flareEffect && !rearView) {

		int flarePos[2];
		int m_viewport[4];
		glGetIntegerv(GL_VIEWPORT, m_viewport);

		pushMatrix(MODEL);
		loadIdentity(MODEL);
		computeDerivedMatrix(PROJ_VIEW_MODEL);  //pvm to be applied to lightPost. pvm is used in project function

		if (!project(lightPos, lightScreenPos, m_viewport))
			printf("Error in getting projected light in screen\n");  //Calculate the window Coordinates of the light position: the projected position of light on viewport
		flarePos[0] = clampi((int)lightScreenPos[0], m_viewport[0], m_viewport[0] + m_viewport[2] - 1);
		flarePos[1] = clampi((int)lightScreenPos[1], m_viewport[1], m_viewport[1] + m_viewport[3] - 1);

		popMatrix(MODEL);

		//viewer looking down at  negative z direction
		pushMatrix(PROJECTION);
		loadIdentity(PROJECTION);
		pushMatrix(VIEW);
		loadIdentity(VIEW);
		ortho(m_viewport[0], m_viewport[0] + m_viewport[2] - 1, m_viewport[1], m_viewport[1] + m_viewport[3] - 1, -1, 1);
		render_flare(&AVTflare, flarePos[0], flarePos[1], m_viewport);
		popMatrix(PROJECTION);
		popMatrix(VIEW);
		glEnable(GL_BLEND);

	}

	//glEnable(GL_BLEND);
	glClear(GL_STENCIL_BUFFER_BIT);

	//Render text (bitmap fonts) in screen coordinates. So use ortoghonal projection with viewport coordinates.
	glDisable(GL_DEPTH_TEST);
	//the glyph contains transparent background colors and non-transparent for the actual character pixels. So we use the 
	int m_viewport[4];
	glGetIntegerv(GL_VIEWPORT, m_viewport);

	//viewer at origin looking down at  negative z direction
	pushMatrix(MODEL);
	loadIdentity(MODEL);
	pushMatrix(PROJECTION);
	loadIdentity(PROJECTION);
	pushMatrix(VIEW);
	loadIdentity(VIEW);
	ortho(m_viewport[0], m_viewport[0] + m_viewport[2] - 1, m_viewport[1], m_viewport[1] + m_viewport[3] - 1, -1, 1);

	if (paused) {
		if (paused == 1) {
			RenderText(shaderText, "Game Paused", WinX / 2 - 150, WinY / 2, 1.0f, 0.3f, 0.7f, 0.9f);
		}
		else if (paused == 2) {
			RenderText(shaderText, "Press Pause [s] to Restart", WinX / 2 - 300, WinY / 2, 1.0f, 0.3f, 0.7f, 0.9f);
		}
		else if (paused == 3) {
			RenderText(shaderText, "Game Over", WinX / 2 - 260, WinY / 2, 2.0f, 1.0f, 0.2f, 0.2f);
			char display_score[20];
			snprintf(display_score, 20, "Final Score: %d", final_score);
			RenderText(shaderText, display_score, WinX / 2 - 160, WinY / 2 - 70, 1.0f, 1.0f, 1.0f, 1.0f);
			RenderText(shaderText, "Press Pause [s] to Restart", WinX / 2 - 210, WinY / 2 - 130, 0.7f, 0.6f, 0.6f, 0.6f);
		}
	}
	char display_lives[19];
	snprintf(display_lives, 19, "Remaining Lives: %d", rover.lives);
	RenderText(shaderText, display_lives, 25.0f, 25.0f, 0.6f, 0.5f, 0.8f, 0.2f);
	char display_current_score[10];
	snprintf(display_current_score, 10, "Score: %d", score);
	RenderText(shaderText, display_current_score, WinX - 150, 25.0f, 0.6f, 1.0f, 1.0f, 1.0f);
	//RenderText(shaderText, "AVT Light and Text Rendering Demo", 440.0f, 570.0f, 0.5f, 0.3, 0.7f, 0.9f);


	popMatrix(PROJECTION);
	popMatrix(VIEW);
	popMatrix(MODEL);
	glEnable(GL_DEPTH_TEST);

    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

	if (rearView) {
		glutSwapBuffers();
	}

	else {

	}


}


void renderFireworks(GLint& loc, int objId)
{
	if (fireworks) {

		float particle_color[4];

		updateParticles();

		glBindTexture(GL_TEXTURE_2D, TextureArray[4]); //particle.tga associated to TU4 

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glDepthMask(GL_FALSE);  //Depth Buffer Read Only

		glUniform1i(texMode_uniformId, 4); // draw modulated textured particles 

		for (int i = 0; i < MAX_PARTICULAS; i++)
		{
			if (particula[i].life > 0.0f) /* só desenha as que ainda estão vivas */
			{

				/* A vida da partícula representa o canal alpha da cor. Como o blend está activo a cor final é a soma da cor rgb do fragmento multiplicada pelo
				alpha com a cor do pixel destino */

				particle_color[0] = particula[i].r;
				particle_color[1] = particula[i].g;
				particle_color[2] = particula[i].b;
				particle_color[3] = particula[i].life;

				// send the material - diffuse color modulated with texture
				loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
				glUniform4fv(loc, 1, particle_color);

				pushMatrix(MODEL);
				translate(MODEL, particula[i].x, particula[i].y, particula[i].z);

				if (i / (MAX_PARTICULAS / 4) == 0) {
					rotate(MODEL, -90, 0, 1, 0);
				}
				else if (i / (MAX_PARTICULAS / 4) == 1) {
					rotate(MODEL, 90, 0, 1, 0);
				}
				else if (i / (MAX_PARTICULAS / 4) == 2) {
					rotate(MODEL, 180, 0, 1, 0);
				}
				else if (i / (MAX_PARTICULAS / 4) == 3) {
					;
				}


				// send matrices to OGL
				computeDerivedMatrix(PROJ_VIEW_MODEL);
				glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
				glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
				computeNormalMatrix3x3();
				glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

				glBindVertexArray(myMeshes[objId].vao);
				glDrawElements(myMeshes[objId].type, myMeshes[objId].numIndexes, GL_UNSIGNED_INT, 0);
				popMatrix(MODEL);
			}
			else dead_num_particles++;
		}

		glDepthMask(GL_TRUE); //make depth buffer again writeable

		if (dead_num_particles == MAX_PARTICULAS) {
			fireworks = 0;
			dead_num_particles = 0;
			printf("All particles dead\n");
		}

	}
}

void renderBillboards(GLint& loc, int objId)
{
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glUniform1i(texMode_uniformId, 3); // draw textured quads
	float pos[3];
	int type = 0;

	for (int i = -3; i < 3; i++)
		for (int j = -3; j < 3; j++) {

			pushMatrix(MODEL);
			translate(MODEL, 5 + i * 15.0, 20.0, 5 + j * 15.0);

			pos[0] = 5 + i * 15.0; pos[1] = 20.0; pos[2] = 5 + j * 15.0;

			if (type == 2)
				l3dBillboardSphericalBegin(cams[activeCam].camPos, pos);
			else if (type == 3)
				l3dBillboardCylindricalBegin(cams[activeCam].camPos, pos);

			//diffuse and ambient color are not used in the tree quads
			loc = glGetUniformLocation(shader.getProgramIndex(), "mat.specular");
			glUniform4fv(loc, 1, myMeshes[objId].mat.specular);
			loc = glGetUniformLocation(shader.getProgramIndex(), "mat.shininess");
			glUniform1f(loc, myMeshes[objId].mat.shininess);

			pushMatrix(MODEL);
			translate(MODEL, 0.0, 3.0, 0.0f);

			// send matrices to OGL
			if (type == 0 || type == 1) {     //Cheating matrix reset billboard techniques
				computeDerivedMatrix(VIEW_MODEL);

				//reset VIEW_MODEL
				if (type == 0) BillboardCheatSphericalBegin();
				else BillboardCheatCylindricalBegin();

				computeDerivedMatrix_PVM(); // calculate PROJ_VIEW_MODEL
			}
			else computeDerivedMatrix(PROJ_VIEW_MODEL);

			glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
			glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
			computeNormalMatrix3x3();
			glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);
			glBindVertexArray(myMeshes[objId].vao);
			glDrawElements(myMeshes[objId].type, myMeshes[objId].numIndexes, GL_UNSIGNED_INT, 0);
			popMatrix(MODEL);

			//	if (type==0 || type==1) // restore matrix VIEW_MODEL não é necessário pois a PVM é sempre calculada a pArtir da MODEL e da VIEW que não são ALTERADAS

			popMatrix(MODEL);
		}
}

void renderSkyBox(GLint& loc, int objId) {
	glUniform1i(texMode_uniformId, 5);

	//it won't write anything to the zbuffer; all subsequently drawn scenery to be in front of the sky box. 
	glDepthMask(GL_FALSE);
	glFrontFace(GL_CW); // set clockwise vertex order to mean the front

	pushMatrix(MODEL);
	pushMatrix(VIEW);  //se quiser anular a translação

	//  Fica mais realista se não anular a translação da câmara 
	// Cancel the translation movement of the camera - de acordo com o tutorial do Antons
	mMatrix[VIEW][12] = 0.0f;
	mMatrix[VIEW][13] = 0.0f;
	mMatrix[VIEW][14] = 0.0f;

	scale(MODEL, 100.0f, 100.0f, 100.0f);
	translate(MODEL, -0.5f, -0.5f, -0.5f);

	// send matrices to OGL
	glUniformMatrix4fv(model_uniformId, 1, GL_FALSE, mMatrix[MODEL]); //Transformação de modelação do cubo unitário para o "Big Cube"
	computeDerivedMatrix(PROJ_VIEW_MODEL);
	glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);

	glBindVertexArray(myMeshes[objId].vao);
	glDrawElements(myMeshes[objId].type, myMeshes[objId].numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	popMatrix(MODEL);
	popMatrix(VIEW);

	glFrontFace(GL_CCW); // restore counter clockwise vertex order to mean the front
	glDepthMask(GL_TRUE);
}

void popModel(int objId)
{
    glUniformMatrix4fv(view_uniformId, 1, GL_FALSE, mMatrix[VIEW]);
	computeDerivedMatrix(PROJ_VIEW_MODEL);
	glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
	glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
	computeNormalMatrix3x3();
	glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);


	// Render mesh
	glBindVertexArray(myMeshes[objId].vao);

	glDrawElements(myMeshes[objId].type, myMeshes[objId].numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	popMatrix(MODEL);
}

void pushModel(GLint& loc, int objId)
{
	//glDisable(GL_BLEND);

	loc = glGetUniformLocation(shader.getProgramIndex(), "mat.ambient");
	glUniform4fv(loc, 1, myMeshes[objId].mat.ambient);
	loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
	glUniform4fv(loc, 1, myMeshes[objId].mat.diffuse);
	loc = glGetUniformLocation(shader.getProgramIndex(), "mat.specular");
	glUniform4fv(loc, 1, myMeshes[objId].mat.specular);
	loc = glGetUniformLocation(shader.getProgramIndex(), "mat.shininess");
	glUniform1f(loc, myMeshes[objId].mat.shininess);

	pushMatrix(MODEL);
}

// ------------------------------------------------------------
//
// Events from the Keyboard
//

void processKeys(unsigned char key, int xx, int yy)
{
	switch (key) {

	case 27:
		glutLeaveMainLoop();
		break;

	case 'm': glEnable(GL_MULTISAMPLE); break;
	case ',': glDisable(GL_MULTISAMPLE); break;
		// Switch cameras
	case '1': activeCam = FIXED; break;
	case '2': activeCam = ORTHO; break;
	case '3': activeCam = FOLLOW; break;
		//Lights
	case 'c': pointsOn = !pointsOn; break;
	case 'n': directionalOn = !directionalOn; break;
	case 'h': spotsOn = !spotsOn; break;
		// Bumpmap 1 light
	case 'b': bumpmap1 = !bumpmap1;  glUniform1i(bumpmap1_loc, bumpmap1); break;
		// Pause
	case 's': paused ? unpause() : pause(); break;
		// Reset game
	case 'r': reset_game(2); break;
		// Rover movement
	case 'q': rover.accel(FWD); break;
	case 'a': rover.accel(RVRS); break;
	case 'p': rover.turn(RIGHT); break;
	case 'o': rover.turn(LEFT); break;
	case 'z': rover.stop(); break;
	case 'l': rover.straighten(); break;
	case 'f': fog = !fog; break;
	case 'x': flareEffect = !flareEffect; break;
	}
}


// ------------------------------------------------------------
//
// Mouse Events
//

void processMouseButtons(int button, int state, int xx, int yy)
{
	// start tracking the mouse
	if (state == GLUT_DOWN) {
		startX = xx;
		startY = yy;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
	}

	//stop tracking the mouse
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			alpha -= (xx - startX);
			beta += (yy - startY);
		}
		else if (tracking == 2) {
			r += (yy - startY) * 0.01f;
			if (r < 0.1f)
				r = 0.1f;
		}
		tracking = 0;
	}
}

// Track mouse motion while buttons are pressed

void processMouseMotion(int xx, int yy)
{

	int deltaX, deltaY;
	float alphaAux, betaAux;
	float rAux;

	deltaX = -xx + startX;
	deltaY = yy - startY;

	// left mouse button: move camera
	if (tracking == 1) {

		alphaAux = alpha + deltaX;
		betaAux = beta + deltaY;

		if (betaAux > 85.0f)
			betaAux = 85.0f;
		else if (betaAux < -85.0f)
			betaAux = -85.0f;
		rAux = r;
	}
	// right mouse button: zoom
	else if (tracking == 2) {

		alphaAux = alpha;
		betaAux = beta;
		rAux = r + (deltaY * 0.01f);
		if (rAux < 0.1f)
			rAux = 0.1f;
	}

	if (activeCam == FOLLOW) {
		cams[activeCam].updateFollowTarget(rover, alphaAux, betaAux);
	}

	//  uncomment this if not using an idle or refresh func
	//	glutPostRedisplay();
}


void mouseWheel(int wheel, int direction, int x, int y) {

	r += direction * 0.1f;
	if (r < 0.1f)
		r = 0.1f;

	if (r < 2) r = 2;

	cams[activeCam].camPos[1] = r;

	//  uncomment this if not using an idle or refresh func
	//	glutPostRedisplay();
}

// --------------------------------------------------------
//
// Shader Stuff
//


GLuint setupShaders() {

	// Shader for models
	shader.init();
	shader.loadShader(VSShaderLib::VERTEX_SHADER, "shaders/pointlight.vert");
	shader.loadShader(VSShaderLib::FRAGMENT_SHADER, "shaders/pointlight.frag");

	// set semantics for the shader variables
	glBindFragDataLocation(shader.getProgramIndex(), 0, "colorOut");
	glBindAttribLocation(shader.getProgramIndex(), VERTEX_COORD_ATTRIB, "position");
	glBindAttribLocation(shader.getProgramIndex(), NORMAL_ATTRIB, "normal");
	glBindAttribLocation(shader.getProgramIndex(), TEXTURE_COORD_ATTRIB, "texCoord");
	glBindAttribLocation(shader.getProgramIndex(), TANGENT_ATTRIB, "tangent");

	glLinkProgram(shader.getProgramIndex());
	printf("InfoLog for Model Rendering Shader\n%s\n\n", shaderText.getAllInfoLogs().c_str());

	if (!shader.isProgramValid()) {
		printf("GLSL Model Program Not Valid!\n");
		//exit(1);
	}

	texMode_uniformId = glGetUniformLocation(shader.getProgramIndex(), "texMode"); // different modes of texturing
	pvm_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_pvm");
	vm_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_viewModel");
	normal_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_normal");
    model_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_Model");
    view_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_View");
    reflect_perFragment_uniformId = glGetUniformLocation(shader.getProgramIndex(), "reflect_perFrag"); //reflection vector calculated in the frag shader

	
	directional_uniformId = glGetUniformLocation(shader.getProgramIndex(), "directional");
	point1_uniformId = glGetUniformLocation(shader.getProgramIndex(), "point1");
	point2_uniformId = glGetUniformLocation(shader.getProgramIndex(), "point2");
	point3_uniformId = glGetUniformLocation(shader.getProgramIndex(), "point3");
	point4_uniformId = glGetUniformLocation(shader.getProgramIndex(), "point4");
	point5_uniformId = glGetUniformLocation(shader.getProgramIndex(), "point5");
	point6_uniformId = glGetUniformLocation(shader.getProgramIndex(), "point6");
	left_uniformId = glGetUniformLocation(shader.getProgramIndex(), "left");
	right_uniformId = glGetUniformLocation(shader.getProgramIndex(), "right");
	lampDir_uniformId = glGetUniformLocation(shader.getProgramIndex(), "lampDirection");
	
	tex_loc = glGetUniformLocation(shader.getProgramIndex(), "texmap");
	tex_loc1 = glGetUniformLocation(shader.getProgramIndex(), "texmap1");
	tex_loc2 = glGetUniformLocation(shader.getProgramIndex(), "texmap2");
	tex_loc3 = glGetUniformLocation(shader.getProgramIndex(), "texmap3");
	tex_loc4 = glGetUniformLocation(shader.getProgramIndex(), "texmap4");
	tex_cube_loc = glGetUniformLocation(shader.getProgramIndex(), "cubeMap");
	tex_normalMap_loc = glGetUniformLocation(shader.getProgramIndex(), "texNormalMap");

	normalMap_loc = glGetUniformLocation(shader.getProgramIndex(), "normalMap");
	specularMap_loc = glGetUniformLocation(shader.getProgramIndex(), "specularMap");
	diffMapCount_loc = glGetUniformLocation(shader.getProgramIndex(), "diffMapCount");

	shadowMode_uniformId = glGetUniformLocation(shader.getProgramIndex(), "shadowMode");

	bumpmap1_loc = glGetUniformLocation(shader.getProgramIndex(), "bumpmap1");

	printf("InfoLog for Per Fragment Phong Lightning Shader\n%s\n\n", shader.getAllInfoLogs().c_str());

	// Shader for bitmap Text
	shaderText.init();
	shaderText.loadShader(VSShaderLib::VERTEX_SHADER, "shaders/text.vert");
	shaderText.loadShader(VSShaderLib::FRAGMENT_SHADER, "shaders/text.frag");

	glLinkProgram(shaderText.getProgramIndex());
	printf("InfoLog for Text Rendering Shader\n%s\n\n", shaderText.getAllInfoLogs().c_str());

	if (!shaderText.isProgramValid()) {
		printf("GLSL Text Program Not Valid!\n");
		exit(1);
	}

	return(shader.isProgramLinked() && shaderText.isProgramLinked());
}

// ------------------------------------------------------------
//
// Model loading and OpenGL setup
//

void init()
{
	MyMesh amesh;
	oldTime = glutGet(GLUT_ELAPSED_TIME);

	/* Initialization of DevIL */
	if (ilGetInteger(IL_VERSION_NUM) < IL_VERSION)
	{
		printf("wrong DevIL version \n");
		exit(0);
	}
	ilInit();

	/// Initialization of freetype library with font_name file
	freeType_init(font_name);

	// set camera positions
	cams[FIXED].updatePos(rover.pos[0], rover.pos[1] + r, rover.pos[2]);
	cams[FIXED].updateTarget(rover.pos[0], rover.pos[1], rover.pos[2]);
	cams[ORTHO].updatePos(rover.pos[0], rover.pos[1] + r, rover.pos[2]);
	cams[ORTHO].updateTarget(rover.pos[0], rover.pos[1], rover.pos[2]);
	cams[FOLLOW].updateFollowPos(rover);
	cams[FOLLOW].updateFollowTarget(rover, 0, 0);

	// Textures
	glGenTextures(6, TextureArray);
	Texture2D_Loader(TextureArray, "stone.tga", 0);
	Texture2D_Loader(TextureArray, "normal.tga", 1);
	Texture2D_Loader(TextureArray, "lightwood.tga", 2);
	Texture2D_Loader(TextureArray, "cloud.png", 3);
	Texture2D_Loader(TextureArray, "particle.tga", 4);

	//Sky Box Texture Object
	// Imagens tem q ser 2048x2048
	const char* filenames[] = { "front.tga", "back.tga", "up.tga", "down.tga", "right.tga", "left.tga" };
	TextureCubeMap_Loader(TextureArray, filenames, 5);

	//Flare elements textures
	glGenTextures(5, FlareTextureArray);
	Texture2D_Loader(FlareTextureArray, "crcl.tga", 0);
	Texture2D_Loader(FlareTextureArray, "flar.tga", 1);
	Texture2D_Loader(FlareTextureArray, "hxgn.tga", 2);
	Texture2D_Loader(FlareTextureArray, "ring.tga", 3);
	Texture2D_Loader(FlareTextureArray, "sun.tga", 4);

	float amb[] = { 0.2f, 0.15f, 0.1f, .5f };
	float diff[] = { 0.9f, 0.2f, 0.2f, .5f };
	float spec[] = { 0.9f, 0.2f, 0.2f, .5f };
	float emissive[] = { 0.0f, 0.0f, 0.0f, .5f };
	float shininess = 100.0f;
	int texcount = 0;

	// Create Flat Terrain
	amesh = createQuad(quadSize, quadSize);

	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
	amesh.mat.shininess = shininess;
	amesh.mat.texCount = texcount;
	myMeshes.push_back(amesh);

	//// create Rover
	float diff_rover[] = { 0.9f, 0.9f, 0.9f, 1.0f };
	float spec_rover[] = { 1.f, 1.f, 1.f, 1.0f };
	float emissive_rover[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	float shininess_rover = 10.0f;

	amesh = createCube();
	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff_rover, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec_rover, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive_rover, 4 * sizeof(float));
	amesh.mat.shininess = shininess_rover;
	amesh.mat.texCount = texcount;
	myMeshes.push_back(amesh);

	//// create wheels
	float diff_wheels[] = { 0.1f, 0.1f, 0.1f, 1.0f };
	float spec_wheels[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	float emissive_wheels[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	float shininess_wheels = 1.0f;
	amesh = createCylinder(0.25f, 0.15f, 20.f);
	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff_wheels, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec_wheels, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive_wheels, 4 * sizeof(float));
	amesh.mat.shininess = shininess_wheels;
	amesh.mat.texCount = texcount;
	myMeshes.push_back(amesh);

	createLights();

	// create rocks
	for (int i = 0; i < N_ROCKS; i++) rockList.push_back(Rock(ROCK_RADIUS));

	float diff_rock[] = { 0.1f, 0.1f, 0.1f, 1.0f };
	float spec_rock[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	float emissive_rock[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	float shininess_rock = 1.0f;
	amesh = createSphere(ROCK_RADIUS, 20.f);
	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff_rock, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec_rock, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive_rock, 4 * sizeof(float));
	amesh.mat.shininess = shininess_rock;
	amesh.mat.texCount = texcount;
	myMeshes.push_back(amesh);

	// CREATE CYLINDERS
	cylinderList.push_back(Cylinder(7, 7));
	cylinderList.push_back(Cylinder(7, -7));
	cylinderList.push_back(Cylinder(-7, -7));
	cylinderList.push_back(Cylinder(-7, 7));

	float diff_cylinder[] = { 0.4f, 0.2f, 0.5f, 1.f };
	float spec_cylinder[] = { 0.6f, 0.8f, 0.2f, 1.f };
	float emissive_cylinder[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	float shininess_cylinder = .2f;
	amesh = createCylinder(5.f, 2.f, 30.f);
	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff_cylinder, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec_cylinder, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive_cylinder, 4 * sizeof(float));
	amesh.mat.shininess = shininess_cylinder;
	amesh.mat.texCount = texcount;
	myMeshes.push_back(amesh);

	// Color Overlay when Paused
	amesh = createQuad(1, 1);

	float amb_pause[] = { 0, 0, 0, 0 };
	float diff_pause[] = { 0, 0, 0, 0 };
	float spec_pause[] = { 0, 0, 0, 0 };
	float emissive_pause[] = { 0, 0, 0, 0 };
	float shininess_pause = 0;

	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff_pause, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec_pause, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive_pause, 4 * sizeof(float));
	amesh.mat.shininess = shininess;
	amesh.mat.texCount = texcount;
	myMeshes.push_back(amesh);


	// Billboards, objId=6;
	float tree_spec[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	float tree_shininess = 10.0f;

	amesh = createQuad(6, 6);
	memcpy(amesh.mat.specular, tree_spec, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
	amesh.mat.shininess = tree_shininess;
	amesh.mat.texCount = texcount;
	myMeshes.push_back(amesh);


	// skybox cube, objId=7;
	amesh = createCube();
	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
	amesh.mat.shininess = shininess;
	amesh.mat.texCount = texcount;
	myMeshes.push_back(amesh);

	// Fireworks, objId = 8;
	amesh = createQuad(2, 2);
	amesh.mat.texCount = texcount;
    myMeshes.push_back(amesh);
	
	// Environment Mapping Cube, objId = 9;
    amesh = createCube();
    memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
    memcpy(amesh.mat.diffuse, diff, 4 * sizeof(float));
    memcpy(amesh.mat.specular, spec, 4 * sizeof(float));
    memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
    amesh.mat.shininess = shininess;
    amesh.mat.texCount = texcount;
    myMeshes.push_back(amesh);

	// create geometry and VAO of the cube
	amesh = createCube();
	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
	amesh.mat.shininess = shininess;
	amesh.mat.texCount = texcount;
	myMeshes.push_back(amesh);


	// Rover OBJ Mesh
	string roverMeshPath = "CuriosityRover/Rover3.obj";
	if (!Import3DFromFile(roverMeshPath))
		return;
	roverMesh = createMeshFromAssimp(scene);


	// create geometry and VAO of the quad for flare elements
	amesh = createQuad(1, 1);
	myMeshes.push_back(amesh);

	//Load flare from file
	loadFlareFile(&AVTflare, "flare.txt");


	// some GL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_MULTISAMPLE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBlendEquation(GL_FUNC_ADD);
    glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
	glClearStencil(1);
}

// ------------------------------------------------------------
//
// Main function
//


int main(int argc, char** argv) {

	//  GLUT initialization
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_MULTISAMPLE | GLUT_STENCIL);

	glutInitContextVersion(4, 3);
	glutInitContextProfile(GLUT_CORE_PROFILE);
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE | GLUT_DEBUG);

	glutInitWindowPosition(100, 100);
	glutInitWindowSize(WinX, WinY);
	WindowHandle = glutCreateWindow(CAPTION);


	//  Callback Registration
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);


	glutTimerFunc(0, timer, 0);
	glutIdleFunc(renderScene);  // Use it for maximum performance
	glutTimerFunc(0, refresh, 0);    //use it to to get 60 FPS whatever

//	Mouse and Keyboard Callbacks
	glutKeyboardFunc(processKeys);
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);
	glutMouseWheelFunc(mouseWheel);


	//	return from main loop
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

	//	Init GLEW
	glewExperimental = GL_TRUE;
	glewInit();

	printf("Vendor: %s\n", glGetString(GL_VENDOR));
	printf("Renderer: %s\n", glGetString(GL_RENDERER));
	printf("Version: %s\n", glGetString(GL_VERSION));
	printf("GLSL: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

	if (!setupShaders())
		return(1);

	init();

	//  GLUT main loop
	glutMainLoop();

	return(0);
}

unsigned int getTextureId(char* name) {
	int i;

	for (i = 0; i < NTEXTURES; ++i)
	{
		if (strncmp(name, flareTextureNames[i], strlen(name)) == 0)
			return i;
	}
	return -1;
}
void    loadFlareFile(FLARE_DEF* flare, char* filename)
{
	int     n = 0;
	FILE* f;
	char    buf[256];
	int fields;

	memset(flare, 0, sizeof(FLARE_DEF));

	f = fopen(filename, "r");
	if (f)
	{
		fgets(buf, sizeof(buf), f);
		sscanf(buf, "%f %f", &flare->fScale, &flare->fMaxSize);

		while (!feof(f))
		{
			char            name[8] = { '\0', };
			double          dDist = 0.0, dSize = 0.0;
			float			color[4];
			int				id;

			fgets(buf, sizeof(buf), f);
			fields = sscanf(buf, "%4s %lf %lf ( %f %f %f %f )", name, &dDist, &dSize, &color[3], &color[0], &color[1], &color[2]);
			if (fields == 7)
			{
				for (int i = 0; i < 4; ++i) color[i] = clamp(color[i] / 255.0f, 0.0f, 1.0f);
				id = getTextureId(name);
				if (id < 0) printf("Texture name not recognized\n");
				else
					flare->element[n].textureId = id;
				flare->element[n].fDistance = (float)dDist;
				flare->element[n].fSize = (float)dSize;
				memcpy(flare->element[n].matDiffuse, color, 4 * sizeof(float));
				++n;
			}
		}

		flare->nPieces = n;
		fclose(f);
	}
	else printf("Flare file opening error\n");
}


