/**************************************
*                                     *
*   Jeff Molofee's Basecode Example   *
*          nehe.gamedev.net           *
*                2001                 *
*                                     *
**************************************/

#include <windows.h>												// Header File For Windows
#include <gl\gl.h>													// Header File For The OpenGL32 Library
#include <gl\glu.h>													// Header File For The GLu32 Library
#include <stdio.h>													// Header File For Standard Input / Output
#include "NeHeGL.h"													// Header File For NeHeGL
#include "glpng.h"		// Header file for accessing the glpng library functions

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>


#include "math.h"												    // NEW: Needed For Sqrtf
#include "ArcBall.h"												// NEW: ArcBall Header

#pragma comment( lib, "opengl32.lib" )								// Search For OpenGL32.lib While Linking
#pragma comment( lib, "glu32.lib" )									// Search For GLu32.lib While Linking

using namespace std;

#ifndef CDS_FULLSCREEN												// CDS_FULLSCREEN Is Not Defined By Some
#define CDS_FULLSCREEN 4											// Compilers. By Defining It This Way,
#endif																// We Can Avoid Errors

#define SCREEN_WIDTH 1024		// We want our screen width 1024 pixels
#define SCREEN_HEIGHT 768		// We want our screen height 768 pixels

GL_Window*	g_window;
Keys*		g_keys;

//////////////////////////////////////////////////////////////
// Need to provide aliases to the pointers and global objects
// so that the font code can work (this code comes from a different source then NeHe)
//HWND g_hWnd = g_window->hWnd;			// This is the handle for the window
//HDC g_hDC = g_window->hDC;				// General HDC - (handle to device context)
//HGLRC g_hRC = g_window->hRC;			// General OpenGL_DC - Our Rendering Context for OpenGL
////////////////////////////////////////


// User Defined Variables
GLUquadricObj *quadratic;											// Used For Our Quadric

///////////////////////
// START: Stuff I added
///////////////////////
const GLuint		numtex = 1;
GLuint		texture[numtex];
GLfloat		xrot		=  -90.0f;						// X Rotation
GLfloat		yrot		=  0.0f;						// Y Rotation
GLfloat		xrotspeed	=  0.0f;						// X Rotation Speed
GLfloat		yrotspeed	=  0.1f;						// Y Rotation Speed
GLfloat		xrosrot		=  0.0f;	// Rossby wave x rotation
GLfloat		yrosrot		=  0.0f;	// Roosby wave y rotation
GLfloat		xrosrotspd	=  0.0f;	// Rossby wave x rotation speed
GLfloat		yrosrotspd	=  0.02f;	// Rossby wave y rotation speed
GLfloat		zoom		= -7.5f;						// Depth Into The Screen
GLfloat		zoomstep	=  0.2f;	// Size of each zoom step
GLfloat		height		=  0.0f;						// Height Of Ball From Floor
bool		isRotate	= true;		// Check to see if we want to rotate the current view
bool		isRotHeld	= false;	// CHeck to see if the right mouse button is being held
bool		isPlusHeld	= false;	// Check to see if the plus key is being held down
bool		isSubHeld	= false;	// Check to see if the minus key is being held.
bool		isfirstload = true;

GLfloat		spheresize	= 2.5f;

// This id is used to know which value of the current texture we are at...
// Initialise to point to the first texture
GLuint texture_id = 1;

// We need to put an upper limit on the number of possible textures. It can
// be any number we choose. Might be better to read this from file
const GLuint max_num_tex = 139;

// Set aside the storage for the wavedata
GLfloat		wavedata[5*max_num_tex];

// This is the maximum number of loaded textures at any one
// instance. If we didn't do things this way we would run
// out of physical RAM very quickly! For now we just
// have 1 loaded texture..possible to speed up if
// we have 3 textures loaded...see how it performs first.
const GLuint max_load_tex = 1;

// Define the structure TextureImage
typedef struct													// Create A Structure
{
	GLubyte	*imageData;											// Image Data (Up To 32 Bits)
	GLuint	bpp;												// Image Color Depth In Bits Per Pixel.
	GLuint	width;												// Image Width
	GLuint	height;												// Image Height
	GLuint	texID;												// Texture ID Used To Select A Texture
} TextureImage;													// Structure Name

// Set aside room for the textures...
TextureImage textures[max_load_tex];		// Storage For max_load_tex textures

double wavespeed = 0.985346363;
double amplitude = 6.764;
double period = 28.25;
int PNGwidth;
int PNGheight;
int PNGalpha;
int PNGdepth;

short zDelta;	// Mouse wheel rotation storage



/////////////////////
// END: Stuff I added
/////////////////////


const float PI2 = 2.0*3.1415926535f;								// PI Squared

Matrix4fT   Transform   = {  1.0f,  0.0f,  0.0f,  0.0f,				// NEW: Final Transform
                             0.0f,  1.0f,  0.0f,  0.0f,
                             0.0f,  0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  0.0f,  1.0f };

Matrix3fT   LastRot     = {  1.0f,  0.0f,  0.0f,					// NEW: Last Rotation
                             0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  1.0f };

Matrix3fT   ThisRot     = {  1.0f,  0.0f,  0.0f,					// NEW: This Rotation
                             0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  1.0f };

ArcBallT    ArcBall(1024.0f, 768.0f);				                // NEW: ArcBall Instance
Point2fT    MousePt;												// NEW: Current Mouse Point
bool        isClicked  = false;										// NEW: Clicking The Mouse?
bool        isRClicked = false;										// NEW: Clicking The Right Mouse Button?
bool        isDragging = false;					                    // NEW: Dragging The Mouse?


// This code takes care of the on screen font stuff
/////// * /////////// * /////////// * NEW * /////// * /////////// * /////////// *
// This is an unsigned int (since we can't have a negative base pointer) that will hold
// the ID for our display list.  It isn't so much of an ID, but it's easier to think of it
// that way.  We we create a display list, we start at 1.  If we create another one we
// then move the list pointer to 2.  g_FontListID holds the BASE number for the display list.
// We are going to create 256 lists for all the characters we need, so we will start at one,
// and the list pointer will then be at 257 if you create another list.  When using display
// lists, we need the base pointer to tell OpenGL which list number we are starting at.
// This will make more sense in the function below.
UINT g_FontListID = 0;

// This will save our old font and select it back in at the end of our program.
// We need this because we use SelectObject() to select in the new font.
// We don't want any memory leaks :)
HFONT hOldFont;

// This define is for the amount of lists we want to create.  There will need to be 1
// for every character.  Since there are 256 ascii characters, we will use 256.
// If we only wanted certain characters used, like alpha numerics, we would give less.
#define MAX_CHARS	256									

// This defines how high we want our font to be
#define FONT_HEIGHT	20


///////////////////////////////// CREATE OPENGL FONT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This function creates a windows font like "Arial" and returns a display list ID
/////
///////////////////////////////// CREATE OPENGL FONT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*

UINT CreateOpenGLFont(LPSTR strFontName, int height)	// Build Our Bitmap Font
{
	UINT	fontListID = 0;								// This will hold the base ID for our display list
	HFONT	hFont;										// This will store the handle to our font

	// Here we generate the lists for each character we want to use.
	// This function then returns the base pointer, which will be 1 because
	// we haven't created any other lists.  If we generated another list after
	// this, the base pointer would be at 257 since the last one used was 256 (which is MAX_CHARS)
	fontListID = glGenLists(MAX_CHARS);					// Generate the list for the font

	// Now we actually need to create the font.  We use a windows function called:
	// CreateFont() that returns a handle to a font (HFONT).  Our CreateOpenGLFont()
	// function allows us to pass in a name and height.  For simplistic reasons, I left
	// other options out, but feel free to add them to your function (like bold, italic, width..)

	hFont = CreateFont(	height,							// Our desired HEIGHT of the font
						0,								// The WIDTH (If we leave this zero it will pick the best width depending on the height)
						0,								// The angle of escapement
						0,								// The angle of orientation
						FW_BOLD,						// The font's weight (We want it bold)
						FALSE,							// Italic - We don't want italic
						FALSE,							// Underline - We don't want it underlined
						FALSE,							// Strikeout - We don't want it strikethrough
						ANSI_CHARSET,					// This is the type of character set
						OUT_TT_PRECIS,					// The Output Precision
						CLIP_DEFAULT_PRECIS,			// The Clipping Precision
						ANTIALIASED_QUALITY,			// The quality of the font - We want anitaliased fonts
						FF_DONTCARE|DEFAULT_PITCH,		// The family and pitch of the font.  We don't care.
						strFontName);					// The font name (Like "Arial", "Courier", etc...)

	// Now that we have created a new font, we need to select that font into our global HDC.
	// We store the old font so we can select it back in when we are done to avoid memory leaks.
	//hOldFont = (HFONT)SelectObject(g_hDC, hFont);
	hOldFont = (HFONT)SelectObject(g_window->hDC, hFont);

	// This function does the magic.  It takes the current font selected in
	// the hdc and makes bitmaps out of each character.  These are called glyphs.
	// The first parameter is the HDC that holds the font to be used.
	// The second parameters is the ASCII value to start from, which is zero in our case.
	// The third parameters is the ASCII value to end on (255 is the last of the ASCII values so we minus 1 from MAX_CHARS)
	// The last parameter is the base pointer for the display lists being used.  This should be 1.

	//wglUseFontBitmaps(g_hDC, 0, MAX_CHARS - 1, fontListID);	// Builds 255 bitmap characters
	wglUseFontBitmaps(g_window->hDC, 0, MAX_CHARS - 1, fontListID);	// Builds 255 bitmap characters

	return fontListID;									// Return the ID to the display list to use later
}

///////////////////////////////// POSITION TEXT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This function sets the drawing position for the text
/////
///////////////////////////////// POSITION TEXT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*

void PositionText( int x, int y )
{
	// If you are to use this font code for your applications,
	// you must be aware that you cannot position the font in 3D,
	// which means you can't rotate and scale it.  That will be covered in
	// the next font tutorial.  BUT, though that might be a drag, this code
	// is useful because when you display the text, it will always be on top
	// of everything else.  This is good if the camera is moving around, and you
	// don't want the text to move.  If the text was positioned in 3D you would have
	// to come up with a tricky way of making it always render in front of the camera.
	// To do this, we need to set the Raster Position.  That is the position that OpenGL
	// starts drawing at.  Since it's in floating point, it's not very intuitive, so what
	// we do is create a new view port, and then always draw the text at (0, 0, 0) in that
	// view port.  The weird part is that the Y is flipped, so (0, 0) is the bottom left corner.
	// Below we do some simple math to flip it back to normal.

	// Before we create a new view port, we need to save the current one we have.
	// This saves our transform (matrix) information and our current viewport information.
	// At the end of this function we POP it back.
	glPushAttrib( GL_TRANSFORM_BIT | GL_VIEWPORT_BIT );

	// Here we use a new projection and modelview matrix to work with.
	glMatrixMode( GL_PROJECTION );						// Set our matrix to our projection matrix
	glPushMatrix();										// Push on a new matrix to work with
	glLoadIdentity();									// reset the matrix
	glMatrixMode( GL_MODELVIEW );						// Set our matrix to our model view matrix
	glPushMatrix();										// Push on a new matrix to work with
	glLoadIdentity();									// Reset that matrix

	// Because the Y is flipped, we want 0 to be at the top, not bottom.
	// If we subtract the font height from the screen height, that should display the
	// font at the top of the screen (if they passed in 0 for Y), but then we subtract
	// the Y from that to get the desired position.  Since the font's drawing point is
	// at the base line of the font, we needed to subtract the font height to make sure
	// if they passed in (0, 0) it wouldn't be off screen.  If you view this in window mode,
	// the top of the window will cut off part of the font, but in full screen it works fine.
	// You just need to add about 25 to the Y to fix that for window mode.

	y = SCREEN_HEIGHT - FONT_HEIGHT - y;				// Calculate the weird screen position

	// Now we create another view port (that is why we saved the old one above).
	// Since glViewPort takes the lower LEFT corner, we needed to change the Y
	// to make it more intuitive when using PositionText().  We minus 1 from the X and Y
	// because 0 is taken into account with the position.  The next 2 parameters are set
	// to 0 for the width and height so it will always draw in the middle of that position.
	// glRasterPos4f() takes (0, 0, 0) as the middle of the viewport, so if we give it a small
	// width/height it will draw at the X and Y given.  Sounds strange, to test this, try
	// using glRasterPos4f(0, 0, 0, 1) instead of PositionText() and you will see, everything
	// will be drawn from the middle.

	glViewport( x - 1, y - 1, 0, 0 );					// Create a new viewport to draw into

	// This is the most important function in here.  This actually positions the text.
	// The parameters are (x, y, z, w).  w should always be 1 , it's a clip coordinate.
	// don't worry about that though.  Because we set the projection and modelview matrix
	// back to the beginning (through LoadIdentity()), the view port is looking at (0, 0, 0).
	// This is the middle, so if we set the drawing position to the middle, it will draw at our
	// X and Y because the width/height of the viewport is 0, starting at X and Y.
	// You can actually call this function (or glRasterPos2f(0, 0)) instead of PositionText(),
	// but it is in floating point and doesn't work as nicely.  You will see why if you try.

	glRasterPos4f( 0, 0, 0, 1 );						// Set the drawing position

	// Now that we positioned the raster position, any text we draw afterwards will start
	// from that position.  Now we just have to put everything else back to normal.

	glPopMatrix();										// Pop the current modelview matrix off the stack
	glMatrixMode( GL_PROJECTION );						// Go back into projection mode
	glPopMatrix();										// Pop the projection matrix off the stack

	glPopAttrib();										// This restores our TRANSFORM and VIEWPORT attributes
}


///////////////////////////////// GL DRAW TEXT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This draws 2D text onto the screen using OpenGL, given an X and Y position 
/////
///////////////////////////////// GL DRAW TEXT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*

void glDrawText(int x, int y, const char *strString, ...)
{
	char		strText[256];							// This will hold our text to display
	va_list		argumentPtr;							// This will hold the pointer to the argument list

	// If you have never used a va_list, listen up.  Remember printf()?
	// or sprintf()?  Well, you can add unlimited arguments into the text like:
	// printf("My name is %s and I am %d years old!", strName, age);
	// Well, that is what va_list's do.  

	// First we need to check if there was even a string given
	if (strString == NULL)								// Check if a string was given
		return;											// Don't render anything then

	// First we need to parse the string for arguments given
	// To do this we pass in a va_list variable that is a pointer to the list of arguments.
	// Then we pass in the string that holds all of those arguments.
	va_start(argumentPtr, strString);					// Parse the arguments out of the string

	// Then we use a special version of sprintf() that takes a pointer to the argument list.
	// This then does the normal sprintf() functionality.
	vsprintf(strText, strString, argumentPtr);			// Now add the arguments into the full string

	va_end(argumentPtr);								// This resets and frees the pointer to the argument list.

	// Before we draw the text, we need to position it with our own function.
	PositionText(x, y);									// Call our own function to position the text on screen

	// Now, before we set the list base, we need to save off the current one.
	glPushAttrib(GL_LIST_BIT);							// This saves the list base information

	// Then we want to set the list base to the font's list base, which should be 1 in our case.
	// That way when we call our display list it will start from the font's lists'.
	glListBase(g_FontListID);							// This sets the lists base

	// Now comes the actually rendering.  We pass in the length of the string,
	// then the data types (which are characters so its a UINT), then the actually char array.
	// This will then take the ASCII value of each character and associate it with a bitmap.
	glCallLists(strlen(strText), GL_UNSIGNED_BYTE, strText);

	glPopAttrib();										// Return the display list back to it's previous state
}


///////////////////////////////// DESTROY FONT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This function frees our display lists
/////
///////////////////////////////// DESTROY FONT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*

void DestroyFont()										
{
	glDeleteLists(g_FontListID, MAX_CHARS);				// Free the display list
	//SelectObject(g_hDC, hOldFont);						// Select the old font back in so we don't have memory leaks
	SelectObject(g_window->hDC, hOldFont);						// Select the old font back in so we don't have memory leaks
}


/////// * /////////// * /////////// * NEW * /////// * /////////// * /////////// *







//////////////////////////////////////////////////////
// START: Stuff that I added to get texture mapping working...
//////////////////////////////////////////////////////

// This is a function for loading a PNG file to a tetxure...works
// very much like LoadTGA except for PNG files
bool LoadPNG(TextureImage *texture, const char *filename)	// Loads a PNG File Into Memory
{
	// Info about the PNG file and if it loaded...
	pngInfo info;

	// Build A Texture id from the Data...only do it if we are in the
	// first texture load call since actually the following function
	// only creates an unique texture ID...not the texture itself etc.
	if (isfirstload == true)
	{
		glGenTextures(1, &texture[0].texID);			// Generate OpenGL texture IDs
		isfirstload = false;		// Set the boolean flag to false now that we've got our first load done.

	}

	//pngSetStandardOrientation(1);

	// Bind Our Texture
	glBindTexture(GL_TEXTURE_2D, texture[0].texID);
	// Set the filter parameters and other texture specific stuff
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);	// Linear Filtered
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);	// Linear Filtered
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	// Now do the loading and texture generation
	if (pngLoad(filename, PNG_NOMIPMAP, PNG_ALPHA, &info)) 
	{
		texture->width = info.Width;
		texture->height = info.Height;
		texture->bpp = info.Depth;
		PNGwidth = texture->width;
		PNGheight = texture->height;
		PNGalpha = info.Alpha;
		PNGdepth = info.Depth;
		return true;
	}
	else 
	{
		// Loading failed so we return false
		return FALSE;
	}
}




//////// Code for loading TGA files
bool LoadTGA(TextureImage *texture, const char *filename)				// Loads A TGA File Into Memory
{    
	GLubyte		TGAheader[12]={0,0,2,0,0,0,0,0,0,0,0,0};		// Uncompressed TGA Header
	GLubyte		TGAcompare[12];									// Used To Compare TGA Header
	GLubyte		header[6];										// First 6 Useful Bytes From The Header
	GLuint		bytesPerPixel;									// Holds Number Of Bytes Per Pixel Used In The TGA File
	GLuint		imageSize;										// Used To Store The Image Size When Setting Aside Ram
	GLuint		temp;											// Temporary Variable
	GLuint		type=GL_RGBA;									// Set The Default GL Mode To RBGA (32 BPP)

	FILE *file = fopen(filename, "rb");							// Open The TGA File

	if(	file==NULL ||											// Does File Even Exist?
		fread(TGAcompare,1,sizeof(TGAcompare),file)!=sizeof(TGAcompare) ||	// Are There 12 Bytes To Read?
		memcmp(TGAheader,TGAcompare,sizeof(TGAheader))!=0				||	// Does The Header Match What We Want?
		fread(header,1,sizeof(header),file)!=sizeof(header))				// If So Read Next 6 Header Bytes
	{
		if (file == NULL)										// Did The File Even Exist? *Added Jim Strong*
			return FALSE;										// Return False
		else													// Otherwise
		{
			fclose(file);										// If Anything Failed, Close The File
			return FALSE;										// Return False
		}
	}

	texture->width  = header[1] * 256 + header[0];				// Determine The TGA Width	(highbyte*256+lowbyte)
	texture->height = header[3] * 256 + header[2];				// Determine The TGA Height	(highbyte*256+lowbyte)
    
 	if(	texture->width	<=0	||									// Is The Width Less Than Or Equal To Zero
		texture->height	<=0	||									// Is The Height Less Than Or Equal To Zero
		(header[4]!=24 && header[4]!=32))						// Is The TGA 24 or 32 Bit?
	{
		fclose(file);											// If Anything Failed, Close The File
		return FALSE;											// Return False
	}

	texture->bpp	= header[4];								// Grab The TGA's Bits Per Pixel (24 or 32)
	bytesPerPixel	= texture->bpp/8;							// Divide By 8 To Get The Bytes Per Pixel
	imageSize		= texture->width*texture->height*bytesPerPixel;	// Calculate The Memory Required For The TGA Data

	// See if we need to reserve memory...only do this if it
	// is the first time we are loading one of the textures
	if (isfirstload == true)
	{
		texture->imageData=(GLubyte *)malloc(imageSize);			// Reserve Memory To Hold The TGA Data
	}

	if(	texture->imageData==NULL ||								// Does The Storage Memory Exist?
		fread(texture->imageData, 1, imageSize, file)!=imageSize)	// Does The Image Size Match The Memory Reserved?
	{
		if(texture->imageData!=NULL)							// Was Image Data Loaded
			free(texture->imageData);							// If So, Release The Image Data

		fclose(file);											// Close The File
		return FALSE;											// Return False
	}

	for(GLuint i=0; i<int(imageSize); i+=bytesPerPixel)			// Loop Through The Image Data
	{															// Swaps The 1st And 3rd Bytes ('R'ed and 'B'lue)
		temp=texture->imageData[i];								// Temporarily Store The Value At Image Data 'i'
		texture->imageData[i] = texture->imageData[i + 2];		// Set The 1st Byte To The Value Of The 3rd Byte
		texture->imageData[i + 2] = temp;						// Set The 3rd Byte To The Value In 'temp' (1st Byte Value)
	}

	fclose (file);												// Close The File

	// Build A Texture From The Data...only do it if we are in the
	// first  texture load call since actually the following function
	// only creates an unique texture ID...not the texture itself etc.
	if (isfirstload == true)
	{
		glGenTextures(1, &texture[0].texID);			// Generate OpenGL texture IDs

	}
	
	glBindTexture(GL_TEXTURE_2D, texture[0].texID);				// Bind Our Texture
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);	// Linear Filtered
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);	// Linear Filtered
	
	if (texture[0].bpp==24)										// Was The TGA 24 Bits
	{
		type=GL_RGB;											// If So Set The 'type' To GL_RGB
	}

	glTexImage2D(GL_TEXTURE_2D, 0, type, texture[0].width, texture[0].height, 0, type, GL_UNSIGNED_BYTE, texture[0].imageData);

	// Free the image data from memory...
	// This is very important and was missed out in the code I got!
	//free(texture->imageData);

	// Set the boolean check to false if this is the first load
	if (isfirstload == true)
	{
		isfirstload = false;
	}

	return true;												// Texture Building Went Ok, Return True
}
////////////////////////////////////


//AUX_RGBImageRec *LoadBMP(char *Filename)				// Loads A Bitmap Image
//{
//	FILE *File=NULL;									// File Handle
//
//	if (!Filename)										// Make Sure A Filename Was Given
//	{
//		return NULL;									// If Not Return NULL
//	}
//
//	File=fopen(Filename,"r");							// Check To See If The File Exists
//
//	if (File)											// Does The File Exist?
//	{
//		fclose(File);									// Close The Handle
//		return auxDIBImageLoad(Filename);				// Load The Bitmap And Return A Pointer
//	}
//
//	return NULL;										// If Load Failed Return NULL
//}

//int LoadGLTextures()                                    // Load Bitmaps And Convert To Textures
//{
//    int Status=FALSE;									// Status Indicator
//    AUX_RGBImageRec *TextureImage[numtex];					// Create Storage Space For The Textures
//    memset(TextureImage,0,sizeof(void *)*numtex);			// Set The Pointer To NULL
//    if ((TextureImage[0]=LoadBMP("data/earthtruecolor.bmp")))// Load the main Earth texture
//   	{   
//		Status=TRUE;									// Set The Status To TRUE
//		glGenTextures(numtex, &texture[0]);					// Create The Texture
//		for (int loop=0; loop<numtex; loop++)				// Loop Through 5 Textures
//		{
//			glBindTexture(GL_TEXTURE_2D, texture[loop]);
//			glTexImage2D(GL_TEXTURE_2D, 0, 3, TextureImage[loop]->sizeX, TextureImage[loop]->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, TextureImage[loop]->data);
//			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
//		}
//		for (loop=0; loop<numtex; loop++)					// Loop Through the Textures
//		{
//			if (TextureImage[loop])						// If Texture Exists
//			{
//				if (TextureImage[loop]->data)			// If Texture Image Exists
//				{
//					free(TextureImage[loop]->data);		// Free The Texture Image Memory
//				}
//				free(TextureImage[loop]);				// Free The Image Structure 
//			}
//		}
//	}
//	return Status;										// Return The Status
//}

//////////////////////////////////////////////////////
// END: Stuff that I added.
//////////////////////////////////////////////////////


BOOL Initialize (GL_Window* window, Keys* keys)						// Any GL Init Code & User Initialiazation Goes Here
{
	/*
	// Load the first texture into
	// positions 0 
	if (!LoadTGA(&textures[0],"Data/1.tga")) 
	{
		return FALSE;											// If Loading Failed, Return False
	}
	//*/

	// Load the PNG texture
	if (!LoadPNG(&textures[0],"Data/1.png"))
	{
		return FALSE;
	}

	// Initialise texture pointer to point to the first texture
	GLuint texture_id = 1;

	// Load all the wavedata from the file wavedata.txt
	// It is stored in the order c, Ae, Ap, Aave, T
	ifstream inWaveData("Data/wavedata.txt", ios::in);
	if (!inWaveData)
	{
		// Couldn't open file
		return FALSE;
	}
	int zz;
	for (zz=0; zz<5*max_num_tex; zz++)
	{
		inWaveData >> wavedata[zz];
	}
	inWaveData.close();
	// Now wavedata contains all the data necessary for display


	// Load the Earth texture from a PNG file
	pngSetStandardOrientation(1);
	pngInfo info;
	texture[0] = pngBind("data/earthtruecolor.png", PNG_NOMIPMAP, PNG_SOLID, &info, GL_CLAMP, GL_LINEAR, GL_LINEAR);
	if (texture[0] == 0)
	{
		return FALSE;
	}

	/*
	// First check to see if we can load the textures
	if (!LoadGLTextures())								// If Loading The Textures Failed
	{
		return FALSE;									// Return False
	}
	*/

	g_window	= window;
	g_keys		= keys;

	// Start Of User Initialization
    isClicked   = false;								            // NEW: Clicking The Mouse?
    isDragging  = false;							                // NEW: Dragging The Mouse?

	// For TGA alpha testures we use the following GL setup
//	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);						// Black Background
//	glClearDepth(1.0f);											// Depth Buffer Setup
//	glDepthFunc(GL_LEQUAL);										// Type Of Depth Testing
//	glEnable(GL_DEPTH_TEST);									// Enable Depth Testing
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);			// Enable Alpha Blending (disable alpha testing)
//	glEnable(GL_BLEND);											// Enable Blending       (disable alpha testing)
//	glAlphaFunc(GL_GREATER,0.1f);								// Set Alpha Testing     (disable blending)
//	glEnable(GL_ALPHA_TEST);									// Enable Alpha Testing  (disable blending)
//	glEnable(GL_TEXTURE_2D);									// Enable Texture Mapping
//	glEnable(GL_CULL_FACE);										// Remove Back Face



	glClearColor (0.0f, 0.0f, 0.0f, 0.0f);							// Black Background
	glClearDepth (1.0f);											// Depth Buffer Setup
	glDepthFunc (GL_LEQUAL);										// The Type Of Depth Testing (Less Or Equal)
	glEnable (GL_DEPTH_TEST);										// Enable Depth Testing
	glShadeModel(GL_SMOOTH);							// Enable Smooth Shading
	glEnable(GL_TEXTURE_2D);							// Enable 2D Texture Mapping
	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);				// Set Perspective Calculations To Most Accurate

	quadratic=gluNewQuadric();										// Create A Pointer To The Quadric Object
	gluQuadricNormals(quadratic, GLU_SMOOTH);						// Create Smooth Normals
	gluQuadricTexture(quadratic, GL_TRUE);							// Create Texture Coords

	glEnable(GL_LIGHT0);											// Enable Default Light
	glEnable(GL_LIGHTING);											// Enable Lighting

	glEnable(GL_COLOR_MATERIAL);									// Enable Color Material
	glEnable(GL_CULL_FACE);							// Remove Back Face


	// After we initialize OpenGL, we want to create our font.  I chose "Arial".
	// We also can pass in a font height.  The width will be chosen according to the height.
	// This then returns the base pointer to the display list for our font.  (Should be 1)
	// We need to be sure and use this to free our display list in DestroyFont().

	g_FontListID = CreateOpenGLFont("Arial", FONT_HEIGHT);

	return TRUE;													// Return TRUE (Initialization Successful)
}

void Deinitialize (void)											// Any User DeInitialization Goes Here
{
	gluDeleteQuadric(quadratic);
	DestroyFont();										// This frees up our font display list
}

void Update (DWORD milliseconds)									// Perform Motion Updates Here
{
	if (g_keys->keyDown [VK_ESCAPE] == TRUE)						// Is ESC Being Pressed?
		TerminateApplication (g_window);							// Terminate The Program

	if (g_keys->keyDown [VK_F1] == TRUE)							// Is F1 Being Pressed?
		ToggleFullscreen (g_window);								// Toggle Fullscreen Mode

	if (g_keys->keyDown ['W'] == TRUE)		// Is 'W' being pressed...if so...zoom in
		zoom += 0.01f;
	if (g_keys->keyDown ['S'] == TRUE)		// Is 'S' being pressed...if so...zoom out
		zoom -= 0.01f;

	if (zDelta > 0)
	{
		// Zoom a bit
		zoom += zoomstep;//0.2f;
		// reset the mouse zoom parameter
		zDelta = 0;
	}
	if (zDelta < 0)
	{
		// Zoom a bit
		zoom -= zoomstep;//0.2f;
		// Reset the mouse zoom parameter
		zDelta = 0;
	}

	// See if the '+' key has been pressed...
	if ((g_keys->keyDown [VK_ADD] == TRUE) && !isPlusHeld)
	{
		isPlusHeld = true;
		// See if it's possible to advance the texture counter
		if (texture_id < max_num_tex)
		{
			// We can advance it so we do it
			texture_id++;
			// Now we need to convert this texture_id into
			// a string so we can load the new texture dynamically.
			// Declare string for later use
			string TexString;
			// Declare string stream for converting from int to string.
			ostringstream StrStream;
			// Copy Integer to String Stream.
			StrStream << texture_id;
			// Assign characters in stream to std::string
			TexString = StrStream.str();
			// Make the base directory  and file name 
			// based on the value of texture_id for TGA
			//string BaseDirectory = "Data\\"+TexString+".tga";
			// Make the base directory  and file name 
			// based on the value of texture_id for PNG
			string BaseDirectory = "Data\\"+TexString+".png";
			// Now load the corresponding texture into textures structure
			LoadPNG(&textures[0],BaseDirectory.c_str());
			//LoadTGA(&textures[0],BaseDirectory.c_str());
		}
		else
		{
			// We are at the last texture so we can't advance any further.
			// Just leave as is and don't do anything...
		}
	}
	// Reset the behaviour of isPlusHeld
	if (g_keys->keyDown [VK_ADD] == FALSE)
	{
		isPlusHeld = false;
	}

	// See if the '-' key has been pressed...
	if ((g_keys->keyDown [VK_SUBTRACT] == TRUE) && !isSubHeld)
	{
		isSubHeld = true;
				// See if it's possible to advance the texture counter
		if (texture_id > 1)
		{
			// We can decrease it so we do it
			texture_id--;
			// Now we need to convert this texture_id into
			// a string so we can load the new texture dynamically.
			// Declare string for later use
			string TexString;
			// Declare string stream for converting from int to string.
			ostringstream StrStream;
			// Copy Integer to String Stream.
			StrStream << texture_id;
			// Assign characters in stream to std::string
			TexString = StrStream.str();
			// Make the base directory  and file name 
			// based on the value of texture_id
			//string BaseDirectory = "Data\\"+TexString+".tga";
			// Make the base directory  and file name 
			// based on the value of texture_id for PNG
			string BaseDirectory = "Data\\"+TexString+".png";
			// Now load the corresponding texture into textures structure
			LoadPNG(&textures[0],BaseDirectory.c_str());
			//LoadTGA(&textures[0],BaseDirectory.c_str());
		}
		else
		{
			// We are at the last texture so we can't advance any further.
			// Just leave as is and don't do anything...
		}
	}
	// Reset the behaviour of isSubHeld
	if (g_keys->keyDown [VK_SUBTRACT] == FALSE)
	{
		isSubHeld = false;
	}

    if (isRClicked && !isRotHeld)	// If Right Mouse Clicked we start/stop Earth rotation based on isRotHeld value
    {
		isRotHeld = true;
		isRotate = !isRotate;		// Reverse the logic of isRotate
    }
	// This is needed to fix the value of holding the image with the right mouse button
	// ...otherwise it will seemingly swap randomly...really just too fast...
	if (!isRClicked)
	{
		// We must have released the right mouse button so reset out boolean flag
		isRotHeld = false;
	}

    if (!isDragging)												// Not Dragging
    {
        if (isClicked)												// First Click
        {
			isDragging = true;										// Prepare For Dragging
			LastRot = ThisRot;										// Set Last Static Rotation To Last Dynamic One
			ArcBall.click(&MousePt);								// Update Start Vector And Prepare For Dragging
        }
    }
    else
    {
        if (isClicked)												// Still Clicked, So Still Dragging
        {
            Quat4fT     ThisQuat;

            ArcBall.drag(&MousePt, &ThisQuat);						// Update End Vector And Get Rotation As Quaternion
            Matrix3fSetRotationFromQuat4f(&ThisRot, &ThisQuat);		// Convert Quaternion Into Matrix3fT
            Matrix3fMulMatrix3f(&ThisRot, &LastRot);				// Accumulate Last Rotation Into This One
            Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);	// Set Our Final Transform's Rotation From This One
        }
        else														// No Longer Dragging
            isDragging = false;
    }
}


void Draw (void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				// Clear Screen And Depth Buffer
	glLoadIdentity();												// Reset The Current Modelview Matrix

	// Turn off lighting before we draw the font to screen
	glDisable(GL_LIGHTING);											// Disble Lighting
	// Draw the text to the screen
	glColor3f(1.0f, 1.0f, 1.0f);	
	// Here we draw our text at (0, 25).  This text will be green.
	glDrawText(20,20, "Wave speed = %.6e ", wavedata[texture_id-1]);
	glDrawText(20,FONT_HEIGHT+20, "Ap = %.2f degrees", wavedata[2*max_num_tex+texture_id-1] );
	glDrawText(20,2*FONT_HEIGHT+20, "Aave = %.2f degrees", wavedata[3*max_num_tex+texture_id-1] );
	glDrawText(20,3*FONT_HEIGHT+20, "Ae = %.2f degrees", wavedata[max_num_tex+texture_id-1] );
	glDrawText(20,4*FONT_HEIGHT+20, "Period = %.2f days", wavedata[4*max_num_tex+texture_id-1] );
	glDrawText(20,5*FONT_HEIGHT+20, "Current tex_id =  %d", texture_id );
	//glDrawText(20,4*FONT_HEIGHT+20, "Width=%d, Height=%d", PNGwidth, PNGheight);
	//glDrawText(20,5*FONT_HEIGHT+20, "Alpha=%d, Depth=%d", PNGalpha, PNGdepth);

	
		
	// Now enable the lighting again
	glEnable(GL_LIGHTING);											// Enable Lighting

	//glTranslatef(0.0f,0.0f,-2.0f);									// Move into The Screen 2.0
	glTranslatef(0.0f, 0.0f, zoom);					// Zoom Camera 	

    glPushMatrix();													// NEW: Prepare Dynamic Transform
		glMultMatrixf(Transform.M);										// NEW: Apply Dynamic Transform
		
		
				
		glRotatef(xrot, 1.0f, 0.0f, 0.0f);					// Rotate On The X Axis
		glRotatef(yrot, 0.0f, 0.0f, 1.0f);					// Rotate On The Y Axis
		// Draw the Earth's surface texture
		glColor3f(1.0f,1.0f,1.0f);
		glBindTexture(GL_TEXTURE_2D, texture[0]);			// Select Texture 1
		gluSphere(quadratic,spheresize,60,60);

		// Here we do the rotation of the rossby wave pattern
		// First we store the current matrix
		glPushMatrix();
			//glRotatef(xrot, 1.0f, 0.0f, 0.0f);					// Rotate On The X Axis
			glRotatef(yrosrot, 0.0f, 0.0f, 1.0f);					// Rotate On The Y Axis
			// Used for TGA alpha channel blending
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);		// Enable Alpha Blending (disable alpha testing)
			
			// Now we enable blending...need to disable depth testing too
			glEnable(GL_BLEND);	// Enable blending
			glDepthMask(GL_FALSE);	// Disable writing to the depth buffer
			// Now we set up masking
	//		glBlendFunc(GL_DST_COLOR,GL_ZERO);	// Blend screen color with Zero (black)
			// Now select the rossby wave mask and bind it and draw it
			glBindTexture(GL_TEXTURE_2D, textures[0].texID);
			gluSphere(quadratic,spheresize*1.01f,60,60);
			// Now apply our rossby wave over the mask
	//		glBlendFunc(GL_ONE,GL_ONE);	// Copy image color to the screen
	//		glBindTexture(GL_TEXTURE_2D, texture[1]);
	//		gluSphere(quadratic,3.05f,60,60);
			glDepthMask(GL_TRUE);	// Enable writing to the depth buffer again
			glDisable(GL_BLEND);	// Disable blending...
		glPopMatrix();
    glPopMatrix();													// NEW: Unapply Dynamic Transform
	glFlush ();	
	
	// Flush The GL Rendering Pipeline

	if (isRotate == true)
	{
		// We have passed the rotate test so we now update all rotation variables.

		// Here we rescale the rotations back to [0.0f,360.0f]...do
		// this so the numbers don't grow too large with time since
		// if left long enough it could eventually overflow etc...
		if (xrot > 360.0f)
		{
			// Reset xrot
			xrot = 0.0f;
		}
		if (yrot > 360.0f)
		{
			// Reset yrot
			yrot = 0.0f;
		}
		if (yrosrot > 360.0f)
		{
			// Reset yrosrot
			yrosrot = 0.0f;
		}

		// Now perform the updates...
		xrot += xrotspeed;				// Update X Rotation Angle By xrotspeed
		yrot += yrotspeed;				// Update Y Rotation Angle By yrotspeed
		yrosrot += yrosrotspd;			// Update rossby wave rotation
	}
	else
	{
		// Don't rotate the earth...just the Rossby wave
		yrosrot += yrosrotspd;
	}
}

