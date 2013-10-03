/***************************************************************************
                          eventhandler.h  -  description
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

#ifndef EVENTHANDLER_H
#define EVENTHANDLER_H


// see http://www.cs.wustl.edu/~schmidt/signal-patterns.html and http://fara.cs.uni-potsdam.de/~kaufmann/?page=GenCppFaqs&faq=Singleton

class EventHandler	
{
public:
	EventHandler() { }
	virtual ~EventHandler() { }
	// Hook method for the signal hook method.
	// Concrete EventHandlers inherit this from EventHandler, implement it and register it with SignalHandler.
	virtual void handle_signal (int signum) = 0;
};

#endif
