/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include "nmod_extra.h"

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*              A few extra functionalities for Fp                        */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* returns the smallest i such that 2^i >= x                              */
/*------------------------------------------------------------------------*/
int next_power_of_two(ulong x)
{
    int i = 0, j = 1;
    while (j < (int) x)
    {
        i++;
        j <<= 1;
    }
    return i;
}

/*------------------------------------------------------------------------*/
/* returns 1/p mod 2^k, assuming p is odd                                 */
/* ill-defined when p is even                                             */
/*------------------------------------------------------------------------*/
ulong inverse_mod_power_of_two(ulong p, int k)
{
  ulong ip = 1L;
  ulong old;
  do
  {
      old = ip;
      ip = (2*ip-p*ip*ip);
  }
  while (ip != old);
  return ip & ((1L << k)-1);
}
