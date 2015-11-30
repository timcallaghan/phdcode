#include "nr.h"
#include "physics.h"
#include "parameters.h"

// Functions that evaluate the equations of motion
// at a given point.

DP mass(Vec_I_DP &vars, const DP phi, const DP wavespeed)
{
	DP value = 0.0;
	value = (vars[0]-Sr*wavespeed*cos(phi))*vars[5]+cos(phi)*vars[1]*vars[8]+(gamma-1.0)/gamma*vars[2]*(vars[3]+cos(phi)*vars[7]-sin(phi)*vars[1]);
	return value;
}

DP lammomen(Vec_I_DP &vars, const DP phi, const DP wavespeed)
{
	DP value = 0.0;
	value = (vars[0]-Sr*wavespeed*cos(phi))*vars[3]+cos(phi)*vars[1]*vars[6]-(cos(phi)/Ro+vars[0])*sin(phi)*vars[1]+vars[5]/(Fr*Fr);
	return value;
}

DP phimomen(Vec_I_DP &vars, const DP phi, const DP wavespeed)
{
	DP value = 0.0;
	value = (vars[0]-Sr*wavespeed*cos(phi))*vars[4]+cos(phi)*vars[1]*vars[7]+(cos(phi)/Ro+vars[0])*sin(phi)*vars[0]+cos(phi)*vars[8]/(Fr*Fr);
	return value;
}
