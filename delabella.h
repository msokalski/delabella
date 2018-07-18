#ifndef DELABELLA_H
#define DELABELLA_H

// returns: positive value: number triangle indices, negative: number of line segment indices (degenerated input)
//          triangle indices in abc array are always returned in CCW order

int DelaBella(int points, const double* xy/*[points][2]*/, int* abc/*[2*points-5][3]*/);

#endif
