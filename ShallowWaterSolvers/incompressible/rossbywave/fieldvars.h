// Function prototypes for evaluating the field
// variables at a given point in the flow grid.

DP ulameval(Vec_I_DP &x, const int i, const int j, const DP phi, Mat_I_DP &Ceta, Mat_I_DP &Cphi1);

DP uphieval(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Seta, Mat_I_DP &Sphi2);

DP heval(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Ceta, Mat_I_DP &Cphi2);


// Function prototypes for evaluating the eta
// derivaties of the field variables at a given 
// point in the flow grid.

DP ulamevaleta(Vec_I_DP &x, const int i, const int j, Mat_I_DP &Seta, Mat_I_DP &Cphi1);

DP uphievaleta(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Ceta, Mat_I_DP &Sphi2);

DP hevaleta(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Seta, Mat_I_DP &Cphi2);

// Function prototypes for evaluating the phi
// derivaties of the field variables at a given 
// point in the flow grid.

DP ulamevalphi(Vec_I_DP &x, const int i, const int j, const DP phi, Mat_I_DP &Ceta, Mat_I_DP &Sphi1);

DP uphievalphi(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Seta, Mat_I_DP &Cphi2);

DP hevalphi(Vec_I_DP &x, const int i, const int j,Mat_I_DP &Ceta, Mat_I_DP &Sphi2);

// Function used in evaluating the volume integral of the free surface
DP surfint(const DP eta, const DP phi, Vec_I_DP &coeffs);

// Functions used to evaluate the jacobian volume elements
DP I10(const int j);

DP I20(const int j, Vec_I_DP &coeffs); 

DP I3(const DP eta, const DP phi, Vec_I_DP &coeffs, const int mval, const int nval);

DP I2i(const int i, const int j, Vec_I_DP &coeffs);