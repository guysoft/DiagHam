#include "AbstractZDensityProfile.h"
#include "InfiniteWellDensityProfile.h"

#include <iostream>
using std::cout;
using std::endl;

AbstractZDensityProfile::~AbstractZDensityProfile()
{
}



// a static class function to return an actual DensityProfile object of some type
AbstractZDensityProfile* AbstractZDensityProfile::CreateZDensityProfile (unsigned int type, double width)
{
  switch(type)
    {
    case AbstractZDensityProfile::InfiniteWell:
      return (AbstractZDensityProfile*) new InfiniteWellDensityProfile(width);
      break;
    default:
      cout << "This type of Density Profile is not defined, yet"<<endl;
      return 0;
      break;
    } 
}

char *AbstractZDensityProfile::DensityProfileName(unsigned int type)
{
  char * buffer = new char[1000];
  switch(type)
    {
    case AbstractZDensityProfile::InfiniteWell:
      sprintf(buffer,"Infinite Well Potential");
      break;
    default:
      sprintf(buffer,"Unknown");
      break;
    } 
  char *rst = new char[strlen(buffer)+1];
  strcpy(rst,buffer);
  delete [] buffer;
  return rst;

}
