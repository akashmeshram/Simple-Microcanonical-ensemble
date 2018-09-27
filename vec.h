#ifndef VEC_H
#define VEC_H
 
class createVector {
public:
   float x;
   float y;
   float z;
   createVector(float x = 0, float y = 0, float z = 0); 
   float getX() const;
   void setX(float x);
   float getY() const;
   void setY(float y);
   float getZ() const;
   void setZ(float z);
   void setXYZ(float x, float y, float z);
   double getMagnitude() const;
};
 
#endif