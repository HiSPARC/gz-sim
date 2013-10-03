/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef SIGINTHANDLER_H
#define SIGINTHANDLER_H

#include "eventhandler.h"
#include <signal.h>

// see http://www.cs.wustl.edu/~schmidt/signal-patterns.html and http://fara.cs.uni-potsdam.de/~kaufmann/?page=GenCppFaqs&faq=Singleton

class SigintHandler : public EventHandler {
public:
	// Initialization.
	SigintHandler() : graceful_quit(0) { }
	virtual ~SigintHandler() { }

	// Hook method. This method is being registered with the SignalHandler and will be called in case of a signal through its dispatcher.
	virtual void handle_signal(int signum) { if (signum == SIGINT) { graceful_quit = 1; } }

	// Accessor.
	sig_atomic_t SigintReceived() { return graceful_quit; }

private:
	sig_atomic_t graceful_quit;
};

#endif
