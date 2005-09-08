////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of matrix main task operation                   //
//                                                                            //
//                        last modification : 10/06/2004                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "MainTask/AbstractMainTask.h"


// constructor 
//
// task = pointer to the main task

MainTaskOperation::MainTaskOperation(AbstractMainTask* task)
{
  this->Task = task;
}


// copy constructor 
//
// operation = reference on operation to copy

MainTaskOperation::MainTaskOperation(const MainTaskOperation& operation)
{
  this->Task = operation.Task;
}

// destructor
//
MainTaskOperation::~MainTaskOperation()
{
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* MainTaskOperation::Clone()
{
  return new MainTaskOperation(*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool MainTaskOperation::RawApplyOperation()
{
  if (this->Task->ExecuteMainTask() != 0)
    return false;
  else
    return true;
}
 
