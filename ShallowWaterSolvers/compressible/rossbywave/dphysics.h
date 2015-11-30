// Function prototypes that evaluate the derivatives
// of the equations of motion
// at a given point.

void dmass(const DP wavespeed, Vec_I_DP &vars,const DP phival,const int i,const int j,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP & Sphi2,Mat_IO_DP &jac);

void dlammomen(const DP wavespeed, Vec_I_DP &vars,const DP phival,const int i,const int j,Mat_I_DP &Ceta,Mat_I_DP SCeta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP & Sphi2,Mat_IO_DP &jac);

void dphimomen(const DP wavespeed, Vec_I_DP &vars,const DP phival,const int i,const int j,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP & Sphi2,Mat_IO_DP &jac);

