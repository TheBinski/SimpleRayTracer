//Author: Caleb Melnychenko
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <stdio.h>

//--- BEGINNING OF STRUCTURES ---

typedef struct {
   float color[3];
   int reflective;
}Material;

typedef struct {
   //camera coords
   float cameraX;
   float cameraY;
   float cameraZ;
   //direction coords
   float dirX;
   float dirY;
   float dirZ;
}Ray;

typedef struct {
   float radius;
   float pos[3];
   Material mat;
}Sphere;

typedef struct {
   float v[3][3];
   Material mat;
}Triangle;

//--- END OF STRUCTURES ---

//--- BEGINNING OF GLOBALS ---

int width, height, channels, debug;
float camX,camY,camZ;

//material
Material refl = { .color = {0,0,0}, .reflective = 1 };
Material red = { .color = {1,0,0}, .reflective = 0 };
Material blue = { .color = {0,0,1}, .reflective = 0 };
Material white = { .color = {1,1,1}, .reflective = 0 };

//spheres
Sphere spheres[100];
int numSpheres;

//Triangles
Triangle triangles[100];
int numTriangles;

//light
float Light[3] = {3,5,-15};

//--- END OF GLOBALS ---

/* fills the image( img ) with the given rgb values */
void fillImage(unsigned char *img, int r, int g, int b)
{
   int totalValues = width * height * channels;
   for ( int i = 0; i < totalValues; i++ ) {
      if ( i % 3 == 0 ) { /* set Red value */
         img[i] = r;
      } else if ( i % 3 == 1 ) { /* set Green value */
         img[i] = g;
      } else { /* set Blue value */
         img[i] = b;
      }
   }
}

/* changes specific (x,y) pixel in image to rgb value */
void changePixel(unsigned char *img, int x, int y, int r, int g, int b)
{
   int position = 0; //keeps track of position

   position += y*(width)*channels; //adds each row that needs to go down
   position += x*channels; //adds rest of rgb values to skip over

   //change pixel values
   img[position] = r;
   img[position+1] = g;
   img[position+2] = b;
}

//returns the dot product of the given vectors
float dot(float x, float y, float z, float x2, float y2, float z2) {
   return (x*x2) + (y*y2) + (z*z2);
}

//returns the rays position and direction starting from Camera
//given the (X,Y) coordinates on the image that is on the Z-axis
//NOTE - Coords given are not actual coords in 3-dimensional plane
//ray = {startX, startY, startZ, directionX, directionY, directionZ}
Ray getRay(int x, int y) {

   //creates Ray that will be found and returned
   Ray RayTemp;
   
   //sets the camera's position in RayTemp
   RayTemp.cameraX = 0;
   RayTemp.cameraY = 0;
   RayTemp.cameraZ = 0;

   //converts pixel coords to the images X,Y,Z coords in 3D space
   //double ratio = 2/width; //ratio from image to pixel
   float newX = ((x-255.5) * 2.0) / 512.0;//converts into image X - has offset of 256 
   float newY = ((255.5-y) * 2.0) / 512.0;//converts into image Y - has offset of 256
   float newZ = -2.0; //default Z value

   //calculates the direction based on given pixel location
   //get magnitude
   float magnitude = sqrt(pow(0-newX,2)+pow(0-newY,2)+pow(0-newZ,2));

   //obtain unit vector (normalize)
   RayTemp.dirX = newX / magnitude;
   RayTemp.dirY = newY / magnitude;
   RayTemp.dirZ = newZ / magnitude;
 
   return RayTemp;
}

//returns a ray starting from the intersection found using r and t, 
//and going towards the light
Ray getLightRay( Ray r, float t ) {
   //creates Ray that will be found and returned
   Ray RayTemp;

   //sets the values for the direction it needs to go
   float newX = r.cameraX + r.dirX * t;
   float newY = r.cameraY + r.dirY * t;
   float newZ = r.cameraZ + r.dirZ * t;

   //sets the starting point
   RayTemp.cameraX = newX;
   RayTemp.cameraY = newY;
   RayTemp.cameraZ = newZ;

   //gets the mag
   float magnitude = sqrt(pow(Light[0]-newX,2)+pow(Light[1]-newY,2)+pow(Light[2]-newZ,2));

   //normalizes
   RayTemp.dirX = (Light[0]-newX) / magnitude;
   RayTemp.dirY = (Light[1]-newY) / magnitude;
   RayTemp.dirZ = (Light[2]-newZ) / magnitude;

   return RayTemp;
}

//returns the intersection point between the ray and the sphere
float sphere_intersect(Ray r, Sphere s) {
   float t[2];
   float dirVec[] = {r.dirX, r.dirY, r.dirZ};
   float camVec[] = {r.cameraX, r.cameraY, r.cameraZ};
   float sphereVec[] = {s.pos[0], s.pos[1], s.pos[2]};

   //(e-c)
   float eSubC[] = { camVec[0]-sphereVec[0], camVec[1]-sphereVec[1], camVec[2]-sphereVec[2] };
   //-d . (e-c)
   float t1 = dot(-dirVec[0], -dirVec[1], -dirVec[2], eSubC[0], eSubC[1], eSubC[2]);
   //d . (e-c)
   float t2 = dot(dirVec[0], dirVec[1], dirVec[2], eSubC[0], eSubC[1], eSubC[2]);
   //d . d
   float ddDot = dot(dirVec[0], dirVec[1], dirVec[2], dirVec[0], dirVec[1], dirVec[2]);
   //(e-c) . (e-c) - R^2
   float t3 = dot( eSubC[0], eSubC[1], eSubC[2], eSubC[0], eSubC[1], eSubC[2]) - (s.radius*s.radius);
   
   //positive value from sqrt
   t[0] = t1 + sqrt((t2*t2)-(ddDot*t3)) / ddDot;

   //negative value from sqrt
   t[1] = t1 - sqrt((t2*t2)-(ddDot*t3)) / ddDot;

   if ( t[0] > 0 && t[1] > 0 ) { //both positive
      if ( t[0] < t[1] ) {
         return t[0];
      } else {
         return t[1];
      }
   } else if ( t[0] <= 0 && t[1] > 0 ) { //t[1] +, t[0] -
      //printf("t0 = %f\n",t[0]);
      return t[1];
   } else if ( t[0] > 0 && t[1] <= 0 ) { //t[0] +, t[1] -
      //printf("t1 = %f\n",t[1]);
      return t[0];
   } else { //both negative
      return -1.0;
   }
}

//returns the intersection point between the ray and the triangle
float tri_intersect(Ray r, Triangle tria) {
   //float v1[3] = { tria.v[0][0], tria.v[0][1], tria.v[0][2] };
   //float v2[3] = { tria.v[1][0], tria.v[1][1], tria.v[1][2] };
   //float v3[3] = { tria.v[2][0], tria.v[2][1], tria.v[2][2] };

   float A = tria.v[0][0]-tria.v[1][0];
   float B = tria.v[0][1]-tria.v[1][1];
   float C = tria.v[0][2]-tria.v[1][2];

   float D = tria.v[0][0]-tria.v[2][0];
   float E = tria.v[0][1]-tria.v[2][1];
   float F = tria.v[0][2]-tria.v[2][2];

   float G = r.dirX;
   float H = r.dirY;
   float I = r.dirZ;

   float J = tria.v[0][0]-r.cameraX;
   float K = tria.v[0][1]-r.cameraY;
   float L = tria.v[0][2]-r.cameraZ;

   float M = A*(E*I-H*F) + B*(G*F-D*I) + C*(D*H-E*G);
   float t = -1.0; 
   t = (-1.0 *(F*(A*K-J*B)+ E*(J*C-A*L)+ D*(B*L-K*C))) / M;

   //check t
   if ( t < 0 ) {
      return -1.0; //no hit
   }
   float gamma = (I*(A*K-J*B)+ H*(J*C-A*L)+ G*(B*L-K*C)) / M;
   if ( gamma < 0 || gamma > 1 ) {
      return -1.0; //no hit
   }
   float beta = (J*(E*I-H*F)+ K*(G*F-D*I)+ L*(D*H-E*G)) / M;
   if (beta < 0 || beta > (1.0-gamma)) {
      return -1.0; //no hit
   }
   //hit is good
   //printf("Gamma = %f, Beta = %f \n", gamma, beta);
   return t; //hit
}

//applies the shadow from the object
int shadow(Ray r, float h) {
	Ray light = getLightRay(r,h);
	float t = -1;
	float t1;
	for ( int i = 0; i < numSpheres; i++ ) {
		t1 = sphere_intersect(light, spheres[i]);
		t = (t1 < 0.001 ? -1 : t1);
		//printf("t = %f, t1 = %f\n", t,t1);
		if ( t != -1 ) {
			return 1;
		}
	}
	for ( int i = 0; i < numTriangles; i++ ) {
		t1 = tri_intersect(light, triangles[i]);
		t = (t1 < 0.001 ? -1 : t1);
		if ( t != -1 ) {
			return 1;
		}
	}
	return 0; //false, no hits
}

//returns diffusion shading float value
float diffShade( Ray r, Ray lr, Sphere s, float t ) {
   //finds normal vector of sphere
   Ray norm;
   float intSect[3]; //intersection of sphere
   intSect[0] = r.cameraX + r.dirX * t;
   intSect[1] = r.cameraY + r.dirY * t;
   intSect[2] = r.cameraZ + r.dirZ * t;
   
   float magnitude = sqrt(pow(intSect[0]-s.pos[0],2)+pow(intSect[1]-s.pos[1],2)+pow(intSect[2]-s.pos[2],2));

   norm.dirX = (intSect[0]-s.pos[0]) / magnitude;
   norm.dirY = (intSect[1]-s.pos[1]) / magnitude;
   norm.dirZ = (intSect[2]-s.pos[2]) / magnitude;

   float d = dot( norm.dirX, norm.dirY, norm.dirZ, lr.dirX, lr.dirY, lr.dirZ );
   if ( shadow(r, t) == 1 ) {
	   return 0.2;
   }
   if ( d < 0.2 ) {
      return 0.2;
   } else {
      return d;
   }
}

float diffShadeT( Ray r, Ray lr, Triangle tria, float t ) {
  
   //get normal vector
   float ca[3] = {tria.v[2][0]-tria.v[0][0], tria.v[2][1]-tria.v[0][1], tria.v[2][2]-tria.v[0][2] };//vertices c - a
   float ba[3] = {tria.v[1][0]-tria.v[0][0], tria.v[1][1]-tria.v[0][1], tria.v[1][2]-tria.v[0][2] };//vertices b - a

   float normX = ba[1] * ca[2] - ba[2] * ca[1];   
   float normY = ba[2] * ca[0] - ba[0] * ca[2];
   float normZ = ba[0] * ca[1] - ba[1] * ca[0];

   float magnitude = sqrt(pow(normX,2) + pow(normY,2) + pow(normZ,2));

   normX = normX / magnitude;
   normY = normY / magnitude;
   normZ = normZ / magnitude;

   float d = dot(normX, normY, normZ, lr.dirX, lr.dirY, lr.dirZ); 
   if ( shadow(r, t) == 1 ) {
	   return 0.2;
   }
   if ( d < 0.2 ) {
      return 0.2;
   } else {
      return d;
   }
   
}

//applies the reflection for a reflective surface
void sphere_reflection(Ray r, Sphere s, float t, int x, int y, unsigned char *img) {
   //finds normal vector of sphere
   Ray norm;
   float intSect[3]; //intersection of sphere
   intSect[0] = r.cameraX + r.dirX * t;
   intSect[1] = r.cameraY + r.dirY * t;
   intSect[2] = r.cameraZ + r.dirZ * t;
   
   float magnitude = sqrt(pow(intSect[0]-s.pos[0],2)+pow(intSect[1]-s.pos[1],2)+pow(intSect[2]-s.pos[2],2));

   norm.dirX = (intSect[0]-s.pos[0]) / magnitude;
   norm.dirY = (intSect[1]-s.pos[1]) / magnitude;
   norm.dirZ = (intSect[2]-s.pos[2]) / magnitude;

   Ray reflect;
   reflect.cameraX = intSect[0];
   reflect.cameraY = intSect[1];
   reflect.cameraZ = intSect[2];
   float d[3] = { r.dirX + intSect[0], r.dirY + intSect[1], r.dirZ + intSect[2] };
   float magnitudeD = sqrt(pow(d[0],2)+pow(d[1],2)+pow(d[2],2));
   //normalize d
   d[0] = d[0] / magnitudeD;
   d[1] = d[1] / magnitudeD;
   d[2] = d[2] / magnitudeD;
   float dDotn = dot(d[0],d[1],d[2],norm.dirX,norm.dirY,norm.dirZ);
   reflect.dirX = d[0] + -2.0*dDotn*norm.dirX;
   reflect.dirY = d[1] + -2.0*dDotn*norm.dirY;
   reflect.dirZ = d[2] + -2.0*dDotn*norm.dirZ;
   
   int numRefl = 0; //keeps track of relections count to prevent infinite loops
   while (1) {
	   //printf("%d\n",numRefl);
	   numRefl++;
      if ( numRefl >= 10 ) {
         //change color to black
		 changePixel(img, x, y, 0, 0, 0);
		 break;
      }

      float h = -1; //tracks intersect value, no hit by default, replaced with g for now
	  float t1; //intersection value found
	  int numOfSphere; //tracks sphere # with closest hit
	  int obj = 0; //tracks current closest obj hit, 0 = none
	  //for each Sphere
 	  for ( int i = 0; i < numSpheres; i++ ) {
		 t1 = sphere_intersect(reflect,spheres[i]);
		 if ( t1 > 0 ) { //intersection exists
		   if ( h == -1 ) { //there is no current intersection
			  h = (t1 < 0.001 ? -1 : t1); //t1 is new intersection
			  obj = 1;
			  numOfSphere = i;
		   } else if ( t1 < h ) {
			  h = (t1 < 0.001 ? -1 : t1);
			  obj = 1;
			  numOfSphere = i;
		   }
		}
	 }
	 //for each Triangle
	 int numOfTri; //tracks triangle # with closest hit
	 for ( int i = 0; i < numTriangles; i++ ) {
		t1 = tri_intersect(reflect,triangles[i]);
		if ( t1 > 0 ) {
		   if ( h == -1 ) {
			  h = (t1 < 0.001 ? -1 : t1);
			  obj = 2;//signals it's a triangle that is hit
			  numOfTri = i;
		   } else if ( t1 < h ) {
			  h = (t1 < 0.001 ? -1 : t1);
			  obj = 2;
			  numOfTri = i;
		   }  
		}
	 }
	 if ( obj != 0 ) { //only execute if an object is hit
		if ( obj == 1 ) { //sphere is closest hit
		   if ( h > 0 ) {
			  if ( spheres[numOfSphere].mat.reflective == 0 ) { //do Diffuse Shading
				 Ray lray = getLightRay( reflect, h );
				 float diffuse = diffShade( reflect, lray, spheres[numOfSphere], h );
				 int r = diffuse * spheres[numOfSphere].mat.color[0] * 255;
				 int g = diffuse * spheres[numOfSphere].mat.color[1] * 255;
				 int b = diffuse * spheres[numOfSphere].mat.color[2] * 255;
				 changePixel(img, x, y, r, g, b);
				 break;
			  } else { //do Reflection
				//reset values for new reflection
				numRefl++;
				intSect[0] = reflect.cameraX + reflect.dirX * h;
				intSect[1] = reflect.cameraY + reflect.dirY * h;
				intSect[2] = reflect.cameraZ + reflect.dirZ * h;
				
				magnitude = sqrt(pow(intSect[0]-spheres[numOfSphere].pos[0],2)+pow(intSect[1]-spheres[numOfSphere].pos[1],2)+pow(intSect[2]-spheres[numOfSphere].pos[2],2));
				
				norm.dirX = (intSect[0]-spheres[numOfSphere].pos[0]) / magnitude;
				norm.dirY = (intSect[1]-spheres[numOfSphere].pos[1]) / magnitude;
				norm.dirZ = (intSect[2]-spheres[numOfSphere].pos[2]) / magnitude;
				
				d[0] = reflect.dirX + intSect[0];
				d[1] = reflect.dirY + intSect[1];
				d[2] = reflect.dirZ + intSect[2];
				
				magnitudeD = sqrt(pow(d[0],2)+pow(d[1],2)+pow(d[2],2));
				
				d[0] = d[0] / magnitudeD;
				d[1] = d[1] / magnitudeD;
				d[2] = d[2] / magnitudeD;
				
				dDotn = dot(d[0],d[1],d[2],norm.dirX,norm.dirY,norm.dirZ);
				reflect.cameraX = intSect[0];
			    reflect.cameraY = intSect[1];
			    reflect.cameraZ = intSect[2];
				reflect.dirX = d[0] + -2.0*dDotn*norm.dirX;
			    reflect.dirY = d[1] + -2.0*dDotn*norm.dirY;
			    reflect.dirZ = d[2] + -2.0*dDotn*norm.dirZ;
				
			  } 
		   } else {
			   changePixel(img, x, y, 0, 0, 0);
			   break;
		   }
		} else if ( obj == 2 ) { //Triangle is closest hit
		   if ( h > 0 ) {
			  if ( triangles[numOfTri].mat.reflective == 0 ) { //do Diffuse Shading
				 Ray lray = getLightRay( reflect, h );
				 float diffuse = diffShadeT(reflect,lray,triangles[numOfTri],h);
				 int r = diffuse * triangles[numOfTri].mat.color[0] * 255;
				 int g = diffuse * triangles[numOfTri].mat.color[1] * 255;
				 int b = diffuse * triangles[numOfTri].mat.color[2] * 255;
				 changePixel(img, x, y, r, g, b);
				 break;
			  } else { //do Reflection
				 //do nothing
				 numRefl++;
			  }
		   } 
		} else { //hit nothing
			changePixel(img, x, y, 0, 0, 0);
			break;
		}
	 } 
	 
   }

   //return 0.0; //default for now

}

int main() {
   //initializes image dimensions and channels
   width = 512;
   height = 512;
   channels = 3;

   //1 => Yes debug, 0 => Don't debug
   debug = 0;

   //initializes the camera's position coords
   camX = 0;
   camY = 0;
   camZ = 0;

   //creates image data
   char *filename = "reference.png";
   unsigned char *image = malloc(width*height*channels); 
   fillImage(image,0,0,0);//fills image with black by default 

   spheres[0] = (Sphere) { .pos = { 0,0,-16 }, .radius = 2, .mat = refl };
   spheres[1] = (Sphere) { .pos = { 3,-1,-14 }, .radius = 1, .mat = refl };
   spheres[2] = (Sphere) { .pos = { -3,-1,-14 }, .radius = 1, .mat = red };
   numSpheres = 3;

   //int numSpheres = 0;

   //Back Wall
   triangles[0] = (Triangle) { .v = { {-8,-2,-20}, {8,-2,-20}, {8,10,-20} }, .mat = blue };
   triangles[1] = (Triangle) { .v = { {-8,-2,-20}, {8,10,-20}, {-8,10,-20} }, .mat = blue };
   //Floor
   triangles[2] = (Triangle) { .v = { {-8,-2,-20}, {8,-2,-10}, {8,-2,-20} }, .mat = white };
   triangles[3] = (Triangle) { .v = { {-8,-2,-20}, {-8,-2,-10}, {8,-2,-10} }, .mat = white };
   //Right Red Triangle
   triangles[4] = (Triangle) { .v = { {8,-2,-20}, {8,-2,-10}, {8,10,-20} }, .mat = red };
   numTriangles = 5;

   /*printf("Triangle 1 = {%f,%f,%f} {%f,%f,%f} {%f,%f,%f}\n",
         triangles[0].v[0][0],triangles[0].v[0][1],triangles[0].v[0][2],
         triangles[0].v[1][0],triangles[0].v[1][1],triangles[0].v[1][2],
         triangles[0].v[2][0],triangles[0].v[2][1],triangles[0].v[2][2]);
   */
   
   Ray r; //Ray from Camera
   Ray lr; //Ray from Light

   for ( int y = 0; y < height; y++ ) {
      for ( int x = 0; x < width; x++ ) {
         float t = -1; //tracks intersect value, no hit by default
         float t1; //intersection value found
         int numOfSphere; //tracks sphere # with closest hit
         int obj = 0; //tracks current closest obj hit, 0 = none
         r = getRay(x,y);
         //for each Sphere
         for ( int i = 0; i < numSpheres; i++ ) {
            t1 = sphere_intersect(r,spheres[i]);
            if ( t1 > 0 ) { //intersection exists
               if ( t == -1 ) { //there is no current intersection
                  t = t1; //t1 is new intersection
                  obj = 1;
                  numOfSphere = i;
               } else if ( t1 < t ) {
                  t = t1;
                  obj = 1;
                  numOfSphere = i;
               }
            }
         }
         //for each Triangle
         int numOfTri; //tracks triangle # with closest hit
         for ( int i = 0; i < numTriangles; i++ ) {
            t1 = tri_intersect(r,triangles[i]);
            if ( t1 > 0 ) {
               if ( t == -1 ) {
                  t = t1;
                  obj = 2;//signals it's a triangle that is hit
                  numOfTri = i;
               } else if ( t1 < t ) {
                  t = t1;
                  obj = 2;
                  numOfTri = i;
               }  
            }
         }
         if ( obj != 0 ) { //only execute if an object is hit
            if ( obj == 1 ) { //sphere is closest hit
               if ( t > 0 ) {
                  if ( spheres[numOfSphere].mat.reflective == 0 ) { //do Diffuse Shading
                     lr = getLightRay( r, t );
                     float diffuse = diffShade( r, lr, spheres[numOfSphere], t );
					 if ( shadow(r, t) == 1 ) {
						 diffuse = 0.2;
					 }
                     int r = diffuse * spheres[numOfSphere].mat.color[0] * 255;
                     int g = diffuse * spheres[numOfSphere].mat.color[1] * 255;
                     int b = diffuse * spheres[numOfSphere].mat.color[2] * 255;
                     changePixel(image, x, y, r, g, b);
                  } else { //do Reflection
                     sphere_reflection(r, spheres[numOfSphere], t, x, y, image);
                  }
               }
            } else { //Triangle is closest hit
               if ( t > 0 ) {
                  if ( triangles[numOfTri].mat.reflective == 0 ) { //do Diffuse Shading
                     lr = getLightRay( r, t );
                     float diffuse = diffShadeT(r,lr,triangles[numOfTri],t);
					 if ( shadow(r, t) == 1 ) {
						 diffuse = 0.2;
					 }
                     int r = diffuse * triangles[numOfTri].mat.color[0] * 255;
                     int g = diffuse * triangles[numOfTri].mat.color[1] * 255;
                     int b = diffuse * triangles[numOfTri].mat.color[2] * 255;
                     changePixel(image, x, y, r, g, b);

                  } else { //do Reflection
                     //reflection();

                  }
               }
            }
         } 
      }
   }

   //writes out to gradient image 
   stbi_write_png(filename,width,height,3,image,width*3);

   //frees image
   free(image);
   return 0;
}
