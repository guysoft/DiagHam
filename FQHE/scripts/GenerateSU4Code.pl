#!/usr/bin/perl -w

use strict 'vars';

&GenerateHamiltonian();


sub GenerateHamiltonian()
  {
    my @ListPairs = ("up|up", "up|um", "up|dp", "up|dm", "um|um", "um|dp", "um|dm", "dp|dp", "dp|dm", "dm|dm");
    
    my $TmpPair;
    my $SpinIsospin1;
    my $SpinIsospin2;

    my $Source = "";

    foreach $TmpPair (@ListPairs)
      {
	($SpinIsospin1,  $SpinIsospin2) = split (/\|/, $TmpPair);
	$Source .= "  	          Coefficient3 = TmpParticles->A".$SpinIsospin1."A".$SpinIsospin2."(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = this->InteractionFactors".$SpinIsospin1.$SpinIsospin2."[j][(i1 * Lim) >> 2];
		      Coefficient3 *= vSource[i];
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad".$SpinIsospin1."Ad".$SpinIsospin2."(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
			  ++TmpInteractionFactor;
			}
		    }\n";

      }

    print $Source;

  }

sub GenerateHilbertSpace()
  {
    my @ListPairs = ("up|up", "up|um", "up|dp", "up|dm", "um|um", "um|dp", "um|dm", "dp|dp", "dp|dm", "dm|dm");
    
    my $Header = "";
    
    my $TmpPair;
    my $SpinIsospin1;
    my $SpinIsospin2;
    
    foreach $TmpPair (@ListPairs)
      {
	($SpinIsospin1,  $SpinIsospin2) = split (/\|/, $TmpPair);
	$Header .= "  // apply a_n1_".$SpinIsospin1." a_n2_".$SpinIsospin2." operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double A".$SpinIsospin1."A".$SpinIsospin2." (int index, int n1, int n2);

";
  }

    foreach $TmpPair (@ListPairs)
      {
	($SpinIsospin1,  $SpinIsospin2) = split (/\|/, $TmpPair);
	$Header .= "  // apply a^+_m1_".$SpinIsospin1." a^+_m2_".$SpinIsospin2." operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int Ad".$SpinIsospin1."Ad".$SpinIsospin2." (int m1, int m2, double& coefficient);

";
      }
    
    
    print $Header;
    
    my $Source = "";
    my %Shifts = ("up", 3, "um", 2, "dp", 1, "dm" ,0);
    foreach $TmpPair (@ListPairs)
      {
	($SpinIsospin1,  $SpinIsospin2) = split (/\|/, $TmpPair);
	$Source .= "// apply a_n1_".$SpinIsospin1." a_n2_".$SpinIsospin2." operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4Spin::A".$SpinIsospin1."A".$SpinIsospin2." (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 2;\n";
  if ($Shifts{$SpinIsospin1} > 1)
    {
$Source .="  n1 += ".$Shifts{$SpinIsospin1}.";\n";
    }
  else
{
  if ($Shifts{$SpinIsospin1} == 1)
    {
       $Source .="  ++n1;\n";
    }
}
  $Source .="  n2 <<= 2;\n";
  if ($Shifts{$SpinIsospin2} > 1)
    {
$Source .="  n2 += ".$Shifts{$SpinIsospin2}.";\n";
    }
  else
{
  if ($Shifts{$SpinIsospin2} == 1)
    {
       $Source .="  ++n2;\n";
    }
}
$Source .=" if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateLargestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

";
}

foreach $TmpPair (@ListPairs)
  {
    ($SpinIsospin1,  $SpinIsospin2) = split (/\|/, $TmpPair);
    $Source .= "// apply a^+_m1_".$SpinIsospin1." a^+_m2_".$SpinIsospin2." operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4Spin::Ad".$SpinIsospin1."Ad".$SpinIsospin2." (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 2;\n";
  if ($Shifts{$SpinIsospin1} > 1)
    {
$Source .="  m1 += ".$Shifts{$SpinIsospin1}.";\n";
    }
  else
{
  if ($Shifts{$SpinIsospin1} == 1)
    {
       $Source .="  ++m1;\n";
    }
}
  $Source .="  m2 <<= 2;\n";
  if ($Shifts{$SpinIsospin2} > 1)
    {
$Source .="  m2 += ".$Shifts{$SpinIsospin2}.";\n";
    }
  else
{
  if ($Shifts{$SpinIsospin2} == 1)
    {
       $Source .="  ++m2;\n";
    }
}
$Source .="  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

";
}


print $Source;
}
