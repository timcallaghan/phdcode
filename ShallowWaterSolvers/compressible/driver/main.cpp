// Program to drive the rossbywave.exe program...
// Basically a script to run it many times...
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

int main(void)
{
	// Set aside memory for start and ending values of the iteration
	int itstart,itend;
	// Now read in how many iterations we are going to perform
	ifstream inNumgoes("numgoes.txt", ios::in);
	if (!inNumgoes)
	{
		cerr << "File numgoes.txt could not be opened" << endl ;
		exit(1);
	}
	while (inNumgoes >> itstart >> itend)
	inNumgoes.close();

	// Now read in what sort of bootstrapping we are doing
	// with this run of the code.
	int bootstrap;
	ifstream inNewparams("newparams.txt", ios::in);
	if (!inNewparams)
	{
		cerr << "File newparams.txt could not be opened" << endl ;
		exit(1);
	}
	while (inNewparams >> bootstrap)
	inNewparams.close();

	// If bootstrap==0 then we are trying to find an
	// appropriate truncation level so we execute the
	// following piece of code which will set up our
	// directories for us...
	if (bootstrap==0)
	{
		for (int i=itstart;i<=itend;i++)
		{
			// Variables to hold the current values of M and N
			// to use in creating a directory...
			int M,N;
			// Read in the old values of M and N
			ifstream inParams("oldparams.txt", ios::in);
			if (!inParams)
			{
				cerr << "File oldparams.txt could not be opened" << endl ;
				exit(1);
			}
			while (inParams >> M >> N)
			inParams.close();
			// Increase their value by 1 to reflect the current truncation
			M++;
			N++;
			// Declare strings for later use
			string MString, NString;
			// Declare string stream for converting from int to string.
			ostringstream StrStream;
			// Copy Integer to String Stream.
			StrStream << M;
			// Assign characters in stream to std::string
			MString = StrStream.str();
			// Do the same for N but first clear StrStream to be the empty string
			StrStream.str("");
			StrStream << N;
			NString = StrStream.str();
			// Make the base directory name based on the values M and N
			string BaseDirectory = "M"+MString+"N"+NString;
			// Make the "make directory" command line string
			string MakeDir = "mkdir "+BaseDirectory;
			// Execute the make directory command
			system(MakeDir.c_str());

			// Execute the bootstraping method with the current values
			// of M and N
			system("rossbywave.exe");

			// Make the copy command to copy the coeffs file to the new directory
			string CopyFiles = "copy coeffs.txt "+BaseDirectory+"\\coeffs.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
			// Make the copy command to copy the current file to the new directory
			CopyFiles = "copy current.txt "+BaseDirectory+"\\current.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
			// Make the copy command to copy the linparams file to the new directory
			CopyFiles = "copy linparams.txt "+BaseDirectory+"\\linparams.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
			// Make the copy command to copy the oldparams file to the new directory
			CopyFiles = "copy oldparams.txt "+BaseDirectory+"\\oldparams.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
			// Make the copy command to copy the residual file to the new directory
			CopyFiles = "copy residual.txt "+BaseDirectory+"\\residual.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
		}
	}
	// If bootstrap==2 then we are trying to find a new
	// solution with a slightly larger amplitude but with
	// the same truncation level so we execute the
	// following piece of code which will set up our
	// directories for us...
	else if (bootstrap==2)
	{
		for (int i=itstart;i<=itend;i++)
		{
			// We will use the current value of "i" to make
			// a new directory for storing our results in.
			// Declare strings for later use
			string iString;
			// Declare string stream for converting from int to string.
			ostringstream StrStream;
			// Copy Integer to String Stream.
			StrStream << i;
			// Assign characters in stream to std::string
			iString = StrStream.str();
			// If i<10 (and >=0) we need to prepend a "0" to the front
			// of the directory structure for easier listing later on...
			string BaseDirectory;
			if (i<10)
			{
				BaseDirectory = "step0"+iString;
			}
			else
			{
				BaseDirectory = "step"+iString;
			}
			// Make the "make directory" command line string
			string MakeDir = "mkdir "+BaseDirectory;
			// Execute the make directory command
			system(MakeDir.c_str());

			// Execute the bootstraping method with the current values
			// of M and N and the new amplitude multiplier that is stored
			// in newparams.txt. NOTE: TO DO: If we need to take smaller steps later on
			// we could modify "scale" here by using "i" somehow...
			system("rossbywave.exe");

			// Make the copy command to copy the coeffs file to the new directory
			string CopyFiles = "copy coeffs.txt "+BaseDirectory+"\\coeffs.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
			// Make the copy command to copy the current file to the new directory
			CopyFiles = "copy current.txt "+BaseDirectory+"\\current.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
			// Make the copy command to copy the linparams file to the new directory
			CopyFiles = "copy linparams.txt "+BaseDirectory+"\\linparams.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
			// Make the copy command to copy the oldparams file to the new directory
			CopyFiles = "copy oldparams.txt "+BaseDirectory+"\\oldparams.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
			// Make the copy command to copy the residual file to the new directory
			CopyFiles = "copy residual.txt "+BaseDirectory+"\\residual.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
			// Make the copy command to copy the newparams file to the new directory
			CopyFiles = "copy newparams.txt "+BaseDirectory+"\\newparams.txt";
			// Execute the copy command
			system(CopyFiles.c_str());
		}
	}
	else
	{
		// We don't have a valid bootstrap value
		cerr << "Invalid bootstrap parameter specified in newparams.txt...must be 0 or 2" << endl ;
		exit(1);
	}
	return 0;
}
