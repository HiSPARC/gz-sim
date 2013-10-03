/***************************************************************************
                          signalhandler.cpp  -  description
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

#include "/user/bdegier/hisparc/lafebre2/GZ/GZ/signalhandler.h"

// initialization of static data members --- see http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.9
EventHandler* SignalHandler::signal_handlers[NSIG];


SignalHandler& SignalHandler::get_instance()
{
	// The instance is created upon first call and destroyed at the end of the program.
	// This implementation is called "Myers-Singleton".
	static SignalHandler the_instance;
	return the_instance;
}


EventHandler* SignalHandler::register_handler (int signum, EventHandler *eh)
{
	// Copy the <old_eh> from the <signum> slot in
	// the <signal_handlers_> table.
	EventHandler* old_eh = signal_handlers[signum];

	// Store <eh> into the <signum> slot in the
	// <signal_handlers_> table.
	signal_handlers[signum] = eh;

	// Register the <dispatcher> to handle this
	// <signum>.
	struct sigaction sa;
	sa.sa_handler = dispatcher;
	sigemptyset (&sa.sa_mask);
	sa.sa_flags = 0;
	sigaction (signum, &sa, 0);
	return old_eh;
}


void SignalHandler::remove_handler (int signum)
{
	signal_handlers[signum] = 0;
}


void SignalHandler::dispatcher (int signum)
{
	// Perform a sanity check...
	if (signal_handlers[signum] != 0)
	{
		// Dispatch the handler's hook method.
		signal_handlers[signum]->handle_signal (signum);
	}
}
