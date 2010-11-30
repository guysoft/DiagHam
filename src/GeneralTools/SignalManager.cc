#include "SignalManager.h"
#include <signal.h>


SignalManager::SignalManager()
{
  if (NbrSignalManagers>0)
    {
      cout << "Currently, only a single SignalManager is supported in each process"<<endl;
    }    
  NbrSignalManagers++;
}

SignalManager::~SignalManager() 
{
  --NbrSignalManagers;
  if (TermHandlerSet)
    {
      signal(SIGTERM,SIG_DFL);
    }
}


// declaration of data fields for signal checker

/* Flag indicating whether error handler was set */
bool SignalManager::TermHandlerSet = false;

/* If this flag is nonzero, don't handle the signal right away. */
volatile sig_atomic_t SignalManager::TermSignalPending;

/* This is nonzero if a signal arrived and was not handled. */
volatile sig_atomic_t SignalManager::TermSignalDeferred=0;

int SignalManager::NbrSignalManagers=0;

// signal handling routines for deferring exit in  with disk storage
// actual signal handler
void SignalManager::TermSignalHandler(int signum)
{
  if (SignalManager::TermSignalDeferred)
    {
      cout << "Deferring TERM signal in SignalManager"<<endl;
      SignalManager::TermSignalPending = signum;
    }
  else
    {
      cout << "Executing deferred signal"<<endl;
      // execute original term-signal
      signal(SIGTERM, SIG_DFL); // register default SIGTERM handler
      raise(SIGTERM);
    }
}

// routine to call to block Term Signal
void SignalManager::StartToDeferTermSignal()
{
  cout << "Activating deferred TERM signal, value"<<SignalManager::TermHandlerSet<<endl;
  /* set error handler if not done already */
  if (SignalManager::TermHandlerSet==false)
    {
      cout << "Really Activating deferred TERM signal"<<endl;
      signal(SIGTERM, &SignalManager::TermSignalHandler);
      SignalManager::TermHandlerSet=true;
    }
  /* Prevent term signals from having immediate effect. */
  cout << "activated before: "<<SignalManager::TermSignalDeferred<<endl;

  SignalManager::TermSignalDeferred++;
  cout << "activated: "<<SignalManager::TermSignalDeferred<<endl;
}

// routine to call to process pending Term Signal
void SignalManager::ProcessDeferredTermSignal()
{
  cout << "Clearing signals:"<<SignalManager::TermSignalDeferred<<", "<<SignalManager::TermSignalPending<<endl;
  if (SignalManager::TermHandlerSet)
    {
      SignalManager::TermSignalDeferred--;
      if (SignalManager::TermSignalDeferred == 0 && SignalManager::TermSignalPending != 0)
	raise (SignalManager::TermSignalPending);
      SignalManager::TermHandlerSet = false;
      signal(SIGTERM, SIG_DFL); // register default SIGTERM handler
    }
}
