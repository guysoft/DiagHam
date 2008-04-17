#include "GroupedPermutations.h"

class PermutationElement
{
public:

  SmallIntegerArray Value;
  unsigned long Multiplicity;
  
  PermutationElement();

  PermutationElement(SmallIntegerArray value, unsigned long multiplicity);
  
  PermutationElement(const PermutationElement &per);

  ~PermutationElement();

  // assignment operator
  //
  // per = permutation to assign
  // return value = reference on current permutation
  PermutationElement& operator = (const PermutationElement& array);

  // multiply on multiplicity
  PermutationElement& operator *= (const unsigned long mult);

  friend bool operator == (const PermutationElement& a1, const PermutationElement& a2);
  friend bool operator < (const PermutationElement& a1,const PermutationElement& a2);
  friend bool operator > (const PermutationElement& a1,const PermutationElement& a2);  
  
};

PermutationElement::PermutationElement()
{
  this->Multiplicity=1;
}

PermutationElement::PermutationElement(SmallIntegerArray value, unsigned long multiplicity)
{
  this->Value = value;
  this->Multiplicity = multiplicity;
}

PermutationElement::PermutationElement(const PermutationElement &per)
{
  this->Value = per.Value;
  this->Multiplicity = per.Multiplicity;
}

PermutationElement::~PermutationElement()
{
}


// assignment operator
//
// per = permutation to assign
// return value = reference on current permutation
PermutationElement& PermutationElement::operator = (const PermutationElement& per)
{
  this->Value = per.Value;
  this->Multiplicity = per.Multiplicity;
  return *this;
}

// multiply on multiplicity
PermutationElement& PermutationElement::operator *= (const unsigned long mult)
{
  this->Multiplicity*=mult;
  return *this;
}


bool operator == (const PermutationElement& a1, const PermutationElement& a2)
{
  return (a1.Value==a2.Value);
}

bool operator < (const PermutationElement& a1,const PermutationElement& a2)
{
  return (a1.Value<a2.Value);
}

bool operator > (const PermutationElement& a1,const PermutationElement& a2)
{
  return (a1.Value>a2.Value);
}


// default constructor
// nbrGroups = number of groups
// elementsPerGroup = number of elements per group
// orderedGroups = flag indicating whether order of groups matters
GroupedPermutations::GroupedPermutations(int nbrGroups, unsigned elementsPerGroup, bool orderedGroups)
{
  this->NbrGroups=nbrGroups;
  this->ElementsPerGroup=elementsPerGroup;
  this->NbrElements = nbrGroups*elementsPerGroup;
  this->OrderedGroups = orderedGroups;
  this->NbrBitsForElements = getHighestBit(NbrElements);
  this->NbrBitsPerGroup = getHighestBit(NbrGroups);  
  this->PermutationList = OrderedList<PermutationElement>(/* eliminateDuplicates */ true);
  this->MyArray = new unsigned[NbrElements];
  this->MapOfGroups = new unsigned[NbrGroups];
  this->InverseMapOfGroups = new unsigned[NbrGroups];
  this->CountOfGroups = new unsigned[NbrGroups];
  for (int i=0; i<NbrGroups; ++i)
    {
      this->MapOfGroups[i]=i;
      this->InverseMapOfGroups[i]=i;
    }

  this->CentralRecursion(this->GetInitialString(),SmallIntegerArray(this->NbrBitsPerGroup), 1);

  this->NbrPermutations = PermutationList.GetNbrElement();
  this->Permutations = new SmallIntegerArray[NbrPermutations];
  this->Multiplicities = new unsigned long[NbrPermutations];

  ListElement<PermutationElement> *ElementPointer = PermutationList.FirstElement;
  int i = 0;
  while (i++<NbrPermutations)
    {
      this->Permutations[i] = ElementPointer->Element.Value;
      this->Multiplicities[i] = ElementPointer->Element.Multiplicity;      
      ElementPointer = ElementPointer->NextPointer;
    }
  this->PermutationList.DeleteList();
  
}


// destructor
GroupedPermutations::~GroupedPermutations()
{
  delete [] MyArray;
  delete [] MapOfGroups;
  delete [] CountOfGroups;
}


// central recursive function that generates all different permutations
void GroupedPermutations::CentralRecursion(SmallIntegerArray remainingElements, SmallIntegerArray permutation, unsigned long multiplicity)
{
  int NbrRemainingElements = remainingElements.GetNbrElements();
  if (NbrRemainingElements>1)
    {      
      for (int i=0; i<NbrRemainingElements; ++i)
	{
	  for (int j=0; j<i; ++j)
	    MyArray[j]=remainingElements.GetElement(j);
	  for (int j=i+1; j<NbrRemainingElements; ++j)
	    MyArray[j-1]=remainingElements.GetElement(j);	  
	  CentralRecursion(SmallIntegerArray(NbrRemainingElements-1, this->NbrBitsPerGroup, MyArray),
			   SmallIntegerArray(permutation, remainingElements.GetElement(i)), 1);
	}
    }
  else // NbrRemainingElements==1
    {
      SmallIntegerArray FinalPermutation(permutation, remainingElements.GetElement(0));
      PermutationElement PE(this->GetPermutationString(FinalPermutation),1);
      int Pos;
      PermutationElement *Duplicate;
      this->PermutationList.Insert(PE, Pos, Duplicate);
      if (Duplicate!=NULL)
	{
	  Duplicate->Multiplicity+=1;
	}
    }
  
}

// translate internal form of permutations to a canonic expression
SmallIntegerArray GroupedPermutations::GetPermutationString(SmallIntegerArray &permutation)
{
  for (int i=0; i<NbrGroups; ++i)
    CountOfGroups[i]=0;
  if (!this->OrderedGroups)
    {
      int PresentGroup=0;
      int PresentElement=1;
      unsigned PresentElementValue;
      this->MapOfGroups[PresentGroup++]=permutation.GetElement(0);
      while (PresentGroup < NbrGroups)
	{
	  PresentElementValue=permutation.GetElement(PresentElement);
	  bool GroupKnown=false;
	  for (int i=0; (i<PresentGroup)&&(!GroupKnown); ++i)
	    if (PresentElementValue==this->MapOfGroups[i])
	      GroupKnown=true;
	  if (!GroupKnown)
	    {
	      this->MapOfGroups[PresentGroup++]=PresentElementValue;
	    }
	}
      for (int i=0; i<NbrGroups; ++i)
	this->InverseMapOfGroups[MapOfGroups[i]]=i;
    }
  for (int i=0; i<NbrElements; ++i)
    {            
      unsigned PresentElementValue=permutation.GetElement(i);
      unsigned PresentGroup=InverseMapOfGroups[PresentElementValue];
      MyArray[i]=PresentGroup*ElementsPerGroup+CountOfGroups[PresentGroup];
      ++CountOfGroups[PresentGroup];
    }
  return SmallIntegerArray(this->NbrElements, this->NbrBitsForElements, MyArray);     
}



// get an initial string without permutations
SmallIntegerArray GroupedPermutations::GetInitialString()
{  
  int count=0;
  for (int g=0; g<NbrGroups; ++g)
    for (int e=0; e<ElementsPerGroup; ++e)
      MyArray[count++]=g;
  return SmallIntegerArray(this->NbrElements, this->NbrBitsPerGroup, MyArray);
}
