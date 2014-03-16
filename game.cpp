#include <cstdlib>
#include <ctime>
#include <iostream>
#include <math.h>
#include <set>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <GL/freeglut.h>
#include <string.h>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <time.h>



#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "imageloader.h"
#include "vec3f.h"

#define ESC 27
#define DEG2RAD(deg) (deg * PI / 180)
#define RAD2DEG(deg) (deg * 180 / PI)
#define PI 3.141592653589

#define start 0
#define ingame 1
#define side 6
#define follow 5
#define heli 4
#define over 3
#define wheel 2
#define driver 1

using namespace std;

int game=start;
const int NUM_GUYS = 1;
const int NUM_FOSSILS = 36;
const float TERRAIN_WIDTH = 140.0f;

int ti=60;
int deltaMove = 0;
float deltaRotate = 0;
int ifrolled;
int TdeltaMove =0;
int view = wheel;
int pause=0;
time_t prev=time(0);
time_t curr;
char text[100000];
float randomFloat() { return (float)rand() / ((float)RAND_MAX + 1); }
void printtext(float x, float y, string String)
{
	//(x,y) is from the bottom left of the window
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_DEPTH_TEST);
	glDisable(GL_DEPTH_TEST);
	glRasterPos2f(x,y);
	for (size_t i=0; i<String.size(); i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, String[i]);
	}
	glPopAttrib();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

//Represents a terrain, by storing a set of heights and normals at 2D locations
class Terrain {
	private:
		int w; //Width
		int l; //Length
		float** hs; //Heights
		Vec3f** normals;
		bool computedNormals; //Whether normals is up-to-date
	public:
		Terrain(int w2, int l2) {
			w = w2;
			l = l2;

			hs = new float*[l];
			for(int i = 0; i < l; i++) {
				hs[i] = new float[w];
			}

			normals = new Vec3f*[l];
			for(int i = 0; i < l; i++) {
				normals[i] = new Vec3f[w];
			}

			computedNormals = false;
		}

		~Terrain() {
			for(int i = 0; i < l; i++) {
				delete[] hs[i];
			}
			delete[] hs;

			for(int i = 0; i < l; i++) {
				delete[] normals[i];
			}
			delete[] normals;
		}

		int width() {
			return w;
		}

		int length() {
			return l;
		}

		//Sets the height at (x, z) to y
		void setHeight(int x, int z, float y) {
			hs[z][x] = y;
			computedNormals = false;
		}

		//Returns the height at (x, z)
		float getHeight(int x, int z) {
			if(z>=0 && x>=0)
				return hs[z][x];
			else
				return 0;
		}

		//Computes the normals, if they haven't been computed yet
		void computeNormals() {
			if (computedNormals) {
				return;
			}

			//Compute the rough version of the normals
			Vec3f** normals2 = new Vec3f*[l];
			for(int i = 0; i < l; i++) {
				normals2[i] = new Vec3f[w];
			}

			for(int z = 0; z < l; z++) {
				for(int x = 0; x < w; x++) {
					Vec3f sum(0.0f, 0.0f, 0.0f);

					Vec3f out;
					if (z > 0) {
						out = Vec3f(0.0f, hs[z - 1][x] - hs[z][x], -1.0f);
					}
					Vec3f in;
					if (z < l - 1) {
						in = Vec3f(0.0f, hs[z + 1][x] - hs[z][x], 1.0f);
					}
					Vec3f left;
					if (x > 0) {
						left = Vec3f(-1.0f, hs[z][x - 1] - hs[z][x], 0.0f);
					}
					Vec3f right;
					if (x < w - 1) {
						right = Vec3f(1.0f, hs[z][x + 1] - hs[z][x], 0.0f);
					}

					if (x > 0 && z > 0) {
						sum += out.cross(left).normalize();
					}
					if (x > 0 && z < l - 1) {
						sum += left.cross(in).normalize();
					}
					if (x < w - 1 && z < l - 1) {
						sum += in.cross(right).normalize();
					}
					if (x < w - 1 && z > 0) {
						sum += right.cross(out).normalize();
					}

					normals2[z][x] = sum;
				}
			}

			//Smooth out the normals
			const float FALLOUT_RATIO = 0.5f;
			for(int z = 0; z < l; z++) {
				for(int x = 0; x < w; x++) {
					Vec3f sum = normals2[z][x];

					if (x > 0) {
						sum += normals2[z][x - 1] * FALLOUT_RATIO;
					}
					if (x < w - 1) {
						sum += normals2[z][x + 1] * FALLOUT_RATIO;
					}
					if (z > 0) {
						sum += normals2[z - 1][x] * FALLOUT_RATIO;
					}
					if (z < l - 1) {
						sum += normals2[z + 1][x] * FALLOUT_RATIO;
					}

					if (sum.magnitude() == 0) {
						sum = Vec3f(0.0f, 1.0f, 0.0f);
					}
					normals[z][x] = sum;
				}
			}

			for(int i = 0; i < l; i++) {
				delete[] normals2[i];
			}
			delete[] normals2;

			computedNormals = true;
		}

		//Returns the normal at (x, z)
		Vec3f getNormal(int x, int z) {
			if (!computedNormals) {
				computeNormals();
			}
			return normals[z][x];
		}
};

//Loads a terrain from a heightmap.  The heights of the terrain range from
//-height / 2 to height / 2.
Terrain* _terrain;
Terrain* loadTerrain(const char* filename, float height) {

	Image* image = loadBMP(filename);
	int  span=2;
	Terrain* t = new Terrain(image->width*span, image->height*span);
	for(int y = 0; y < image->height*span; y++) {
		for(int x = 0; x < image->width*span; x++) {
			unsigned char color =
				(unsigned char)image->pixels[3 * ((y%image->height) * image->width + x%image->width)];
			float h = height * ((color / 255.0f) - 0.5f);
			if(h<0)
				t->setHeight(x, y, 0);
			else
				t->setHeight(x, y, h);

		}
	}

	delete image;
	t->computeNormals();
	return t;
}

//Returns the approximate height of the terrain at the specified (x, z) position
float heightAt(Terrain* terrain, float x, float z) {
	//Make (x, z) lie within the bounds of the terrain
	if (x < 0) {
		x = 0;
	}
	else if (x > terrain->width() - 1) {
		x = terrain->width() - 1;
	}
	if (z < 0) {
		z = 0;
	}
	else if (z > terrain->length() - 1) {
		z = terrain->length() - 1;
	}

	//Compute the grid cell in which (x, z) lies and how close we are to the
	//left and outward edges
	int leftX = (int)x;
	if (leftX == terrain->width() - 1) {
		leftX--;
	}
	float fracX = x - leftX;

	int outZ = (int)z;
	if (outZ == terrain->width() - 1) {
		outZ--;
	}
	float fracZ = z - outZ;

	//Compute the four heights for the grid cell
	float h11 = terrain->getHeight(leftX, outZ);
	float h12 = terrain->getHeight(leftX, outZ + 1);
	float h21 = terrain->getHeight(leftX + 1, outZ);
	float h22 = terrain->getHeight(leftX + 1, outZ + 1);

	//Take a weighted average of the four heights
	return (1 - fracX) * ((1 - fracZ) * h11 + fracZ * h12) +
		fracX * ((1 - fracZ) * h21 + fracZ * h22);
}

class Fossil {
	Terrain* terrain;
	float terrainScale; //The scaling factor for the terrain
	float yscale;
	float angle;
	public:
	float x0;
	float z0;
	Fossil(Terrain* terrain1, float terrainScale1)
	{
		terrain = terrain1;
		terrainScale = terrainScale1;
		x0 = randomFloat() *
			(terrainScale * (terrain->width() - 1) );
		z0 = randomFloat() *
			(terrainScale * (terrain->length() - 1) );
		yscale = 0.5f;
		angle=0.0f;
	}
	void draw()
	{

		glColor3f(0.309804, 0.184314, 0.184314); // set drawing color to white
		//glColor3f(0.52, 0.37, 0.26); // set drawing color to white
		glPushMatrix();
		glTranslatef(x0, y()+0.5f , z0);
		//glRotatef(90.0,1.0,0.0,0.0);
		glRotatef(angle,0.0,1.0,0.0);
		// glTranslatef(0.0, 0.0, 0.75);
		glutSolidCylinder(0.5,0.2, 20, 20);
		glPopMatrix();

	}
	void update()
	{
		angle+=10;
		if(angle>360)
			angle=0;
	}
	float x() { return x0;}

	float z() { return z0;}

	//Returns the current height of the guy on the terrain
	float y() { return terrainScale*heightAt(terrain, x0 / terrainScale, z0 / terrainScale);}

};
//Returns a vector of numGuys new guys
vector<Fossil*> makeFossils(int numFossils,  Terrain* terrain) 
{
	vector<Fossil*> fossils;
	for(int i = 0; i < numFossils; i++)
		fossils.push_back(new Fossil( terrain, TERRAIN_WIDTH / (terrain->width() - 1)));
	return fossils;
}
float getmagnitude(float a,float b,float c)
{
	return sqrt((a*a)+(b*b)+(c*c));
}
float distance(float a1,float a2,float b1,float b2,float c1, float c2)
{
	return sqrt((a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2));
}
void drawBox(float len, int boxtype) {
	glPolygonMode(GL_FRONT_AND_BACK, boxtype);
	glBegin(GL_QUADS);
	glVertex2f(-len / 2, -len / 2);
	glVertex2f(len / 2, -len / 2);
	glVertex2f(len / 2, len / 2);
	glVertex2f(-len / 2, len / 2);
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}
vector<Fossil*> _fossils;

//Represents a guy
class Guy {
	private:
		Terrain* terrain;
		float terrainScale; //The scaling factor for the terrain
		float x0;
		float z0;
		float y0;
		float vel;
	public:
		int jump;
		int brake;
		float acc;
		float ret;
		float angle;
		float radius0;
		//float pitch;
		float yaw;
		float pitch;
		float prevpitch;
		float roll;
		float tyre_angle ; 
		float prevheight;
		float h1,h2;

		Guy( Terrain* terrain1, float terrainScale1) 
		{
			terrain = terrain1;
			terrainScale = terrainScale1;

			//Initialize certain fields to random values
			prevheight=-1;
			radius0 = 0.4f * randomFloat() + 0.25f;
			vel = 0.0f;
			angle = 0;
			roll = 0.0f;
			tyre_angle = 0.0f;
			yaw = 0.0f;
			pitch = 0.0f;
			prevpitch = 0.0f;
			acc = 0.0f;
			ret = 0.0f;
			jump=0;
			brake=0;
			x0=0.0f;
			z0=0.0f;
			y0=0.375;
		}

		void draw() 
		{
			glColor3f(1.0,1.0,0.0); 
			glPushMatrix();
			glTranslatef(x0,0,z0);
			glRotatef(yaw,0,1,0);
			glRotatef(-pitch,1,0,0);
			glRotatef(roll,0,0,1);
			glTranslatef(0,y0+0.05,0);//to keep it till ground
			glPushMatrix();
			glTranslatef(0.0,0.0,-0.125);
			glutSolidCube(0.25);
			glPopMatrix();
			glPushMatrix();
			glTranslatef(0.0,0.0,0.125);
			glutSolidCube(0.25);
			glTranslatef(0.0,0.0,0.125);

			//headlight

			//headlight
			glPushMatrix();
			glColor3f(1.0,1.0,1.0);
			glTranslatef(0.0,0.175,0.0);
			GLfloat light1_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
			GLfloat light1_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
			GLfloat light1_specular[] = { 1.0, 1.0, 1.0, 1.0 };
			GLfloat light1_position[] = { 0.0, 0.0, 0.0, 1.0 };
			GLfloat spot_direction[] = { 1.0, 0.0, 1.0 };

			glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
			glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
			glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
			glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
			glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.5);
			glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.2);
			glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.1);

			glLightf(GL_LIGHT1, GL_SPOT_CUTOFF, 90.0);
			glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, spot_direction);
			glLightf(GL_LIGHT1, GL_SPOT_EXPONENT, 2.0);

			glEnable(GL_LIGHT1);
			glutSolidCone(0.05,0.1,20,20);
			glPopMatrix();



			glPushMatrix();
			glColor3f(1.0,1.0,1.0); 
			glTranslatef(0.0,0.175,0.0);
			glutSolidCone(0.05,0.1,20,20);
			glPopMatrix();

			//handle
			glPushMatrix();
			glColor3f(1.0,0.0,0.0); 
			glLineWidth(3);
			glBegin(GL_LINES);
			glVertex3f(0.125,0.125,0.0);
			glVertex3f(0.375,0.375,0.0);
			glVertex3f(-0.125,0.125,0.0);
			glVertex3f(-0.375,0.375,0.0);
			glEnd();
			glPopMatrix();

			glPopMatrix();


			tyre_angle+=vel/0.06;

			//back wheel
			glPushMatrix();
			glColor3f(0.0,0.0,0.0); 
			glTranslatef(0.0,-0.125,-0.375);
			glRotatef(-90,1,0,0);
			glRotatef(90,0,1,0);
			glRotatef(RAD2DEG(tyre_angle),0,0,1);//wheel rotation
			glutWireTorus(0.06,0.12,10,10);
			glPopMatrix();

			//front wheel
			glPushMatrix();
			glColor3f(0.0,0.0,0.0);
			glTranslatef(0.0,-0.125,0.375);
			glRotatef(-90,1,0,0);
			glRotatef(90,0,1,0);
			glRotatef(RAD2DEG(tyre_angle),0,0,1);//wheel rotation
			glutWireTorus(0.06,0.12,10,10);
			glPopMatrix();

			glPopMatrix();
		}

		void step() {
			if(jump)
			{
				glutPostRedisplay(); // redisplay everything
				return;
			}
			else
			{
				if(prevpitch<0 && pitch==0)
				{
					ret=-0.005;
				}
				prevpitch=pitch;
				if(vel>0)
					vel=vel+acc+ret-(0.001*sinf(DEG2RAD(pitch)));
				else //if(vel==0)
				{
					ret=0;
					// 
					//if(pitch<0)
					// {
					vel=vel+acc-(0.001*sinf(DEG2RAD(pitch)));
					//  }
					//  else
					//vel=vel+acc;
				}
				if(deltaRotate==-1)
					yaw+=5;
				else if(deltaRotate==1)
					yaw-=5;
				if( -45< (roll+(5*ifrolled)) && roll+(5*ifrolled)<45)
					roll=roll+(5*ifrolled);
				if(ifrolled)
					if(vel!=0 && vel>0.05)
					{
						yaw+=-0.7*tanf(DEG2RAD(roll))/vel;

					}
			}
			for(int i=0;i<_fossils.size();i++)
			{
				if(_fossils[i]->x0<x0+0.4f && _fossils[i]->x0>x0-0.4f && _fossils[i]->z0<z0+0.4f && _fossils[i]->z0>z0-0.4f)
					_fossils.erase(_fossils.begin()+i);
			}

		}
		void pre_draw()
		{
			if(pause==1 or game!=ingame)
				return;
			if(brake)
			{
				vel=0;
				ret=0;
				acc=0;
			}
			if(vel>=0 || pitch!=0) 
			{
				z0+=vel*cosf(DEG2RAD(yaw));
				x0+=vel*sinf(DEG2RAD(yaw));
			}
			else if(vel<0)
			{
				vel=0;
				acc=0;
				ret=0;
			}
			if(z0<0)
				z0=0;
			if(x0<0)
				x0=0;
			if(z0>117)
				z0=117;
			if(x0>117)
				x0=117;

			h1=heightAt(_terrain,x0,z0);
			h2=h1;
			//cout<<prevheight<<" "<<h2<<" "<<jump<<"\n";
			if(jump==1)
				h1=prevheight-0.1;
			if(prevheight-(h2)>0.2)
				jump=1;
			else if(h1<0)
				jump=0;
			prevheight=h1;
			Vec3f tnormal = _terrain->getNormal(x0, z0);
			if(!jump)
			{
				float tryc,ab,ac,ad;//,ac,ab,ac;
				ab=sinf(DEG2RAD(yaw));
				ac=0;//sinf(DEG2RAD(pitch));
				ad=cosf(DEG2RAD(yaw));
				tryc=((tnormal[0]*ab+tnormal[2]*ad+tnormal[1]*ac)/(getmagnitude(tnormal[0],tnormal[1],tnormal[2])*getmagnitude(ab,ac,ad)));
				tryc=RAD2DEG(acos(tryc));
				pitch=-(90-tryc);
			}
			y0=h1+0.25;
		}

		float x() { return x0;}

		float z() { return z0;}

		float y() { return y0;}


};

//Returns a vector of numGuys new guys
vector<Guy*> makeGuys(int numGuys,  Terrain* terrain) 
{
	vector<Guy*> guys;
	for(int i = 0; i < numGuys; i++)
		guys.push_back(new Guy( terrain, TERRAIN_WIDTH / (terrain->width() - 1)));
	return guys;
}

//Draws the terrain
void drawTerrain(Terrain* terrain) 
{
	glDisable(GL_TEXTURE_2D);
	glColor3f(0.3f, 0.9f, 0.0f);
	for(int z = 0; z < terrain->length() - 1; z++) 
	{
		glBegin(GL_TRIANGLE_STRIP);
		for(int x = 0; x < terrain->width(); x++) 
		{
			Vec3f normal = terrain->getNormal(x, z);
			glNormal3f(normal[0], normal[1], normal[2]);
			glVertex3f(x, terrain->getHeight(x, z), z);
			normal = terrain->getNormal(x, z + 1);
			glNormal3f(normal[0], normal[1], normal[2]);
			glVertex3f(x, terrain->getHeight(x, z + 1), z + 1);
		}
		glEnd();
	}
}

vector<Guy*> _guys;
float _angle = 0;

void cleanup() 
{
	for(unsigned int i = 0; i < _guys.size(); i++)
		delete _guys[i];   
}

void handleKeypress(unsigned char key, int x, int y) 
{
	switch (key) 
	{
		case 27: //Escape key
			cleanup();
			exit(0);
	}
}

GLuint loadTexture(Image *image) 
{
	GLuint textureId;
	glGenTextures(1, &textureId);
	glBindTexture(GL_TEXTURE_2D, textureId);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, image->width, image->height, 0, GL_RGB, GL_UNSIGNED_BYTE, image->pixels);
	return textureId;
}

void initRendering() 
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glShadeModel(GL_SMOOTH);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
}

void handleResize(int w, int h) 
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (float)w / (float)h, 1.0, 200.0);
}

void Set_Camera()
{
	if(view==follow)		//Follow Cam
	{
		gluLookAt(
				_guys[0]->x()-(3.5*sin(DEG2RAD(_guys[0]->yaw))), _guys[0]->y()+1,_guys[0]->z()-(3.5*cos(DEG2RAD(_guys[0]->yaw))),
				_guys[0]->x(),      _guys[0]->y()+1,      _guys[0]->z(),
				0.0,    1.0,    0.0
			 );
	}
	else if(view==heli)		// Helicopter View
	{
		gluLookAt(
				_guys[0]->x()-3,     _guys[0]->y(),_guys[0]->z(),
				_guys[0]->x(), _guys[0]->y()+0.25, _guys[0]->z(),
				0.0,    1.0,    0.0
			 );
	}
	else if(view==wheel)		// Wheel View 
	{
		gluLookAt(
				_guys[0]->x(),      _guys[0]->y()+0.25,      _guys[0]->z(),
				_guys[0]->x()+2*sinf(DEG2RAD(_guys[0]->yaw)), _guys[0]->y()+0.25,_guys[0]->z()+2*cosf(DEG2RAD(_guys[0]->yaw)),
				0.0,    1.0,    0.0
			 );
	}
	else if(view==over)		// Overhead view
	{
		gluLookAt(
				_guys[0]->x(),      _guys[0]->y()+10,      _guys[0]->z(),
				_guys[0]->x()+sin(DEG2RAD(_guys[0]->yaw)), _guys[0]->y(),_guys[0]->z()+cos(DEG2RAD(_guys[0]->yaw)),
				0.0,    1.0,    0.0
			 );

	}
	else if(view==driver)		// Driver View
	{
		gluLookAt(
				_guys[0]->x()-(2*sin(DEG2RAD(_guys[0]->yaw))), _guys[0]->y()+1,_guys[0]->z()-(2*cos(DEG2RAD(_guys[0]->yaw))),
				_guys[0]->x(),      _guys[0]->y()+1,      _guys[0]->z(),
				0.0,    1.0,    0.0
			 );
	}
	else if(view==side)		// Side View
	{
		gluLookAt(
				//_guys[0]->x() , _guys[0]->y() ,_guys[0]->z()+4,
				_guys[0]->x() -5*sin(DEG2RAD(_guys[0]->yaw)), _guys[0]->y()+1.0f,_guys[0]->z()-5*cos(DEG2RAD(_guys[0]->yaw)),
				_guys[0]->x()  , _guys[0]->y()+1.0f,_guys[0]->z(),
				0.0,    1.0,    0.0
			 );
	}
	else if(view==0)
	{
		float scale = TERRAIN_WIDTH / (_terrain->width() - 1);
		glTranslatef(0.0f,0.0f,0.0f);
		glTranslatef(0, 0, -1.0f * scale * (_terrain->length() - 1));
		glRotatef(30, 1, 0, 0);
		glRotatef(90, 0, 1, 0);
		glTranslatef(-TERRAIN_WIDTH / 2, 0, -scale * (_terrain->length() - 1) / 2);
	}
}

void processNormalKeys(unsigned char key, int xx, int yy)
{
	if (key == ESC || key == 'q' || key == 'Q') {cleanup();exit(0);}
	if(key==13 && game==start) // spacebar 
	{
		game=ingame;
		//system("aplay explosion.wav &");
	}
	if(game!=ingame)	
		return;
	if(key==112) {pause=(pause+1)%2;}

	if(game!=ingame or pause==1)
		return;
	if(key == '0')
		view = 0;
	else if(key=='1')
		view=driver;
	else if(key=='2')
		view=wheel;
	else if(key=='3')
		view=over;
	else if(key=='4')
		view=heli;
	else if(key=='5')
		view=follow;
	else if(key=='6')
		view=side;
} 

void drawScene() 
{
	if (game==start)
	{
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glPushMatrix();


		// Draw Frame

		glColor3f(1.0, 1.0, 1.0);
		drawBox(4.0,GL_FILL);

		glColor3f(1.0,0.0,0.0);
		sprintf(text,"Start Game \n");
		printtext(.0f,.0f,string(text));

		glColor3f(0.309804, 0.184314, 0.184314);
		sprintf(text,"Press Enter to Start \n ");
		printtext(-0.1f,-0.07f,string(text));
		glPopMatrix();
		glutSwapBuffers();

	}
	else
	{
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		if(ti<=0)
			ti=0;
		glColor3f(1.0,0.0,0.0);
		sprintf(text,"Time Left: %d",ti);
		printtext(0.6,0.9,string(text));
		if(_fossils.size()==0)
		{
			glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			sprintf(text,"Game Over\nScore: %d",60-ti);
			printtext(.0f,.0f,string(text));
			game=over;
			glutSwapBuffers();
			return;
		}
		else if(ti<=0)
		{
			glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			int sz=_fossils.size();
			sprintf(text,"Game Over \n You Lost");
			printtext(.0f,.0f,string(text));
			game=over;
			glutSwapBuffers();
			return;
		}
		_guys[0]->pre_draw();	
		Set_Camera();
		float scale = TERRAIN_WIDTH / (_terrain->width() - 1);
		GLfloat ambientLight[] = {0.5f, 0.5f, 0.5f, 1.0f};
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);

		GLfloat lightColor[] = {0.5f, 0.5f, 0.5f, 1.0f};
		GLfloat lightPos[] = {-0.2f, 0.3f, -1, 0.0f};
		glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
		glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
		glColor3f(0.3f, 0.9f, 0.0f);
		for(int z = 0; z < _terrain->length() - 1; z++) {
			//Makes OpenGL draw a triangle at every three consecutive vertices
			glBegin(GL_TRIANGLE_STRIP);
			for(int x = 0; x < _terrain->width(); x++) {
				Vec3f normal = _terrain->getNormal(x, z);
				glNormal3f(normal[0], normal[1], normal[2]);
				glVertex3f(x, _terrain->getHeight(x, z), z);
				normal = _terrain->getNormal(x, z + 1);
				glNormal3f(normal[0], normal[1], normal[2]);
				glVertex3f(x, _terrain->getHeight(x, z + 1), z + 1);
			}
			glEnd();
		}

		//drawTerrain(_terrain);

		_guys[0]->draw();
		glColor3f(0.0,0.1,0.8);
		glBegin(GL_QUADS); glVertex3f(33.0f,0.2f,3.0f); glVertex3f(117.0f,0.2f,43.0f); glVertex3f(51.0f,0.2f,0.6f); glVertex3f(117.0f,0.2f,31.0f); 			glEnd();
		for (int i=0;i<_fossils.size();i++)
			_fossils[i]->draw();

		//glScalef(scale, scale, scale);

		glutSwapBuffers();
	}
}

void update() 
{
	if(game!=ingame or pause==1)
		return;
	_guys[0]->step();	
	for(int i=0;i<_fossils.size();i++)
	{
		_fossils[i]->update();	
	}
	time(&curr);
	if(difftime(curr,prev)>1)
	{	
		prev=curr;
		ti--;
	}
	glutPostRedisplay();
	//glutTimerFunc(25, update, 0);
}

void pressSpecialKey(int key, int xx, int yy)
{
	if(game!=ingame or pause==1)
		return;
	switch (key) 
	{
		case GLUT_KEY_UP : deltaMove = 1.0 ;_guys[0]->acc=0.005;break; 
		case GLUT_KEY_DOWN : deltaMove = -1.0;_guys[0]->brake=1; break;
		case GLUT_KEY_LEFT : ifrolled=-1;break; 
		case GLUT_KEY_RIGHT : ifrolled=1; break;
	}
} 

void releaseSpecialKey(int key, int x, int y) 
{
	if(game!=ingame or pause==1)
		return;
	switch (key) 
	{
		case GLUT_KEY_UP : 
			deltaMove = 0.0; 
			_guys[0]->ret=-0.005;
			_guys[0]->acc=0;
			break;
		case GLUT_KEY_DOWN : 
			deltaMove = 0.0; 
			_guys[0]->acc=0;
			_guys[0]->brake=0;
			break;
		case GLUT_KEY_LEFT :
			deltaRotate = 0.0; 
			_guys[0]->roll=0;
			ifrolled=0;
			break;
		case GLUT_KEY_RIGHT : 
			deltaRotate = 0.0; 
			_guys[0]->roll=0;
			ifrolled=0;
			break;
	}
} 

int main(int argc, char** argv) 
{
	srand((unsigned int)time(0)); //Seed the random number generator

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500, 500);

	glutCreateWindow("MotoGame");
	initRendering();

	_terrain = loadTerrain("heightmap.bmp", 30.0f); //Load the terrain
	_guys = makeGuys(NUM_GUYS, _terrain); //Create the guys
	_fossils = makeFossils(NUM_FOSSILS, _terrain); //Create the guys
	float scaledTerrainLength = TERRAIN_WIDTH / (_terrain->width() - 1) * (_terrain->length() - 1);

		glutSpecialUpFunc(releaseSpecialKey); // process special key release*/
	glutDisplayFunc(drawScene);
	glutReshapeFunc(handleResize);
	glutIdleFunc(update); // incremental update 
	glutIgnoreKeyRepeat(1); // ignore key repeat when holding key down
	glutKeyboardFunc(processNormalKeys); // process standard key clicks
	glutSpecialFunc(pressSpecialKey); // process special key pressed
	//glutKeyboardUpFunc(KeyUp);
	// Warning: Nonstandard function! Delete if desired.
	glutSpecialUpFunc(releaseSpecialKey); // process special key release

	glutMainLoop();
	return 0;
}
