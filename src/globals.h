//      globals.h
//      
//      Copyright 2011 Alexandros Sigalas <alxarch@gmail.com>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#ifndef GLOBALS
#define GLOBALS

#define IN 1
#define FR 1
#define ON 2
#define OUT 4
#define BK 4
#define ONIN 3 	//ON | IN
#define ONFR 3 	//ON | FR
#define ONOUT 6 //ON | OUT
#define ONBK 6 	//ON | BK
#define CROSS 5 //IN | OUT
#define ONCRS 7 //IN | ON | OUT

#define E 8E-1
#define ZEROLENGTH 4E-1
#define ZEROAREA 5E-1

typedef enum axis_pair{ XY, XZ, YZ} AxisPair;
enum R_Mode{ORTHO, PERSP};

#endif
