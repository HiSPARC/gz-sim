/***************************************************************************
                          signalhandler.h  -  description
                             -------------------
    begin                : Wed Mar 10 2004
    copyright            : (C) 2004 by Tim Huege
    email                : tim.huege@ik.fzk.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "eventhandler.h"
#include <signal.h>

#ifndef SIGNALHANDLER_H
#define SIGNALHANDLER_H

// see http://www.cs.wustl.edu/~schmidt/signal-patterns.html and http://fara.cs.uni-potsdam.de/~kaufmann/?page=GenCppFaqs&faq=Singleton
// the SignalHandler is implemented as a singleton

class SignalHandler {
public:
	// This method is the only way to get access to the instance.
	// It creates the single instance if it has not yet been created.
	static SignalHandler& get_instance();

	// Register an event handler <eh> for <signum>
	// and return a pointer to any existing <Event_Handler>
	// that was previously registered to handle <signum>.
	EventHandler* register_handler (int signum, EventHandler *eh);

	// Remove the <Event_Handler> for <signum>
	// by setting the slot in the <signal_handlers>
	// table to NULL.
	void remove_handler (int signum);

private:
	// These declarations limit the class to be a true singleton.
	SignalHandler() { }																// Standard constructor private
	SignalHandler(const SignalHandler&);							// Copy constructor private
	SignalHandler& operator= (const SignalHandler&);		// Assignment operator private

	// Entry point adapter installed into <sigaction> (must be a static method or a stand-alone extern "C" function).
	// Can't be a member function of some class because a pointer to a member function cannot be passed to the signal handler!
	// This method looks up the corresponding handler for the given signal and calls it
	static void dispatcher (int signum);

	// Table of pointers to concrete <Event_Handler>s
	// registered by applications.  NSIG is the number of
	// signals defined in </usr/include/sys/signal.h>.
	static EventHandler* signal_handlers[NSIG];
};


#endif
