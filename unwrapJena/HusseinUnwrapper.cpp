// Based on the source code by Hussein Abdul-Rahman and Munther Gdeisat on 22nd August 2007
// Based on the paper entitled "Fast three-dimensional phase-unwrapping algorithm based on sorting by reliability following a noncontinuous path"
// by Hussein Abdul-Rahman, Munther A. Gdeisat, David R. Burton, and Michael J. Lalor, published in the Applied Optics, Vol. 46, No. 26, pp. 6623-6635, 2007.
// Enhanced by B.Lehr (Friedrich-Schiller-University Jena, Germany; p2lebe@uni-jena.de) 2008-2009
// Further enhanced by B.Lehr (Unicersity Hospital Jena, Germany; p2lebe@uni-jena.de) 2010

// Changes to previouse version:
// - debugging, debugging, debugging
// - replaced manually implemented q-sort with c++'s qsort
// - added support for analyze (.hdr/.img)-file format (reading header)
// - runtime is meassured
// - now supports 0..4096-16bit-int & -π..π-32bit-float datatype
// - writes out the quality map (capriciousity) of every voxel
// - removed waste variables & code-blocks
// - changed variable names into consistend naming sheme (first letter hints type)
// - has version and SVN-revision number

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "analyze_db.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <wchar.h>
#include <errno.h>

// Information mostly updated by svn
#define VERSION  "1.01"
#define REVISION "$Revision: 763 $"
#define DATE     "$Date$"
#define AUTHOR   "$Author$"


static float    PI = 3.141592654;
static float TWOPI = 6.283185307;
struct dsr hdr;
char bPrint = 1;

/* in the voxel struct two types of information are stored:   /
/  1. the value, increment and capriciousity of the voxel, a  /
/     information directly connected to the very voxel.       /
/  2. the number_of_voxels in_group and a pointer to the      /
/     first, next, and last voxel in the group. Hereby the    /
/     first voxel in a group of voxels stores also the infor- /
/     mation about the whole group.                          */

struct VOXEL
{
    int increment;                  // Number of 2*pi to add to the voxel to unwrap it with respect to other voxels in the same group
    int number_of_voxels_in_group;  // Number of voxels in group (only applies if voxel is head voxel of group)
    float value;                    // Value of the voxel
    float capriciousity;            // Quality meassure, inverse grade of reliability
    struct VOXEL *head;             // Pointer to the first/head voxel in the group
    struct VOXEL *last;             // Pointer to the last voxel in the group
    struct VOXEL *next;             // Pointer to the next (following) voxel in the group
};

struct EDGE                         // Line connecting two voxels
{
    float capriciousity;            // Capriciousity of the edge, equal to the sum of the capriciousity of pointer_1 and pointer_2
    VOXEL *pointer_1, *pointer_2;   // Pointer to both voxels the edge is connecting.
    int increment;                  // Number of 2*pi to add to the 2nd voxel to unwrap it with respect to the first one
};

struct DIMS
{
    int  w, h, d, size;             // volume_width, _height, _depth and number of voxels
};


// Handles the status string initializing a work step
void _beginFunction(const char cInfo[])
{
  if (bPrint)
  {
    printf("%s ", cInfo);
    for (int i = strlen(cInfo); i < 55; i++)
    {
      printf(".");
    }
    printf(" ");
  }
}

// Message displaed when a step has been completed
void _endFunction()
{
    if (bPrint)
    {
        printf("Done\n");
    }
}

void ErrorMessage(const char cErrormessage[])
{
    #define CTRL_SCREEN 0x1B
    #define BRIGHT 1
    #define LOW 0
    printf("%c[%dm%s%c[%dm\n", CTRL_SCREEN, BRIGHT, cErrormessage, CTRL_SCREEN, LOW);
}

void fileError(FILE *fHandler, const char cFilename[])
{
    int iFilenamelength = 21;

    if (fHandler == NULL)
    {
        iFilenamelength += strlen(cFilename);
		char *cErrormessage = new char[iFilenamelength];
        sprintf(cErrormessage, "Error opening file %s", cFilename);
        ErrorMessage(cErrormessage);
        exit(1);
    }
}

// Read the *.hdr file and extract the dataType and Dimension
void read_header(const char cFilename[], DIMS* pdDimensions, int &iDatatype)
{
    _beginFunction("Reading the header file");
    char cHeaderfile[255];
    FILE *fHeaderfile;
    struct image_dimension sDimensions;
    const int iMaxDim = 3;

    strcpy(cHeaderfile, cFilename);
    strcpy(cHeaderfile + strlen(cHeaderfile) - 3, "hdr");
    fHeaderfile = fopen(cHeaderfile, "rb");
    fileError(fHeaderfile, cHeaderfile);
    fread(&hdr, sizeof(struct dsr), 1, fHeaderfile);
    fclose(fHeaderfile);
    sDimensions = hdr.dime;
    _endFunction();

    iDatatype = sDimensions.bitpix;
    if ((sDimensions.dim[0] > iMaxDim) && (sDimensions.dim[4] > 1))
    {
        char *cErrormessage = new char[67];
        sprintf(cErrormessage, "More then %i dimensions are not supported in this version (has %i dimensions)", iMaxDim, sDimensions.dim[0]);
        ErrorMessage(cErrormessage);
        for (int i = iMaxDim; i < sDimensions.dim[0] ; i++) {
            sprintf(cErrormessage, "Dimension %i has a length of %i", i, sDimensions.dim[i+1]);
            ErrorMessage(cErrormessage);
        }
        exit(1);
    }
    if (bPrint)
    {
        printf(" ├─ Datasize : %i x %i x %i (%i byte)\n", sDimensions.dim[1],    sDimensions.dim[2],    sDimensions.dim[3],    iDatatype/8);
        printf(" ├─ Voxelsize: %.2f x %.2f x %.2f mm³\n", sDimensions.pixdim[1], sDimensions.pixdim[2], sDimensions.pixdim[3]);
        printf(" └─ ROI-Size : %.0f x %.0f x %.0f mm³\n", sDimensions.dim[1]*sDimensions.pixdim[1],     sDimensions.dim[2]*sDimensions.pixdim[2],
                                                          sDimensions.dim[3]*sDimensions.pixdim[3]);
    }
    if ((iDatatype != 32) && (iDatatype != 16))
    {
        ErrorMessage("Error: wrong dataformat, should be float (32 bit or 16 bit per voxel)\n");
        exit(1);
    }
    (*pdDimensions).w = sDimensions.dim[1];
    (*pdDimensions).h = sDimensions.dim[2];
    (*pdDimensions).d = sDimensions.dim[3];
    (*pdDimensions).size = sDimensions.dim[1] * sDimensions.dim[2] * sDimensions.dim[3];
}

// Read all data from the *.img file and convert it into float
void read_data(char *cFilename, float *dData, int iLength, int iDatatype)
{

    _beginFunction("Reading wrapped values from binary file");
    FILE *fDatafile;

    strcpy(cFilename + strlen(cFilename) - 3, "img");
    fDatafile = fopen(cFilename,"rb");
    fileError(fDatafile, cFilename);
    switch (iDatatype)
    {
        case 32:
        {
            fread(dData, sizeof(float), iLength, fDatafile);
            _endFunction();
            _beginFunction("Leaving Values in float: (-π..π) ......................");
            break;
        }
        case 16:
        {
		    ushort *uData = new ushort[iLength];
            fread(uData, sizeof(ushort), iLength, fDatafile);
            _endFunction();
            _beginFunction("Converting ushort: (0..4096) to float: (-π..π) ........");
            for (int i = 0; i < iLength; i++)
            {
                dData[i] = ((float)uData[i]/2048 - 1) * PI;
            }
            break;
        }
        default:
            ErrorMessage("Format not supported");
            exit(1);
    }
    fclose(fDatafile);
    _endFunction();
}

// Write all data to *.img file
void write_data(char cFilename[], float *dData, int iLength, const char cDatatype[])
{
    char statusMessage[22+strlen(cDatatype)];
    sprintf(statusMessage, "Writing %s to binary file", cDatatype);
    _beginFunction(statusMessage);
    FILE *fDatafile;

    strcpy(cFilename + strlen(cFilename) - 3, "img");
    fDatafile = fopen(cFilename,"wb");
    fileError(fDatafile, cFilename);
    fwrite(dData, sizeof(float), iLength, fDatafile);
    fclose(fDatafile);
    _endFunction();
}

// Copy the *.hdr-input file and define datatype to be float
void write_header(char *cFilename)
{
    _beginFunction("Writing the header file");
    FILE *fDatafile;

    strcpy(cFilename + strlen(cFilename) - 3, "hdr");
    fDatafile = fopen(cFilename,"wb");
    fileError(fDatafile, cFilename);
    hdr.dime.bitpix = 32;
    hdr.dime.datatype = 16; //DT_FLOAT
    fwrite(&hdr,sizeof(struct dsr), 1, fDatafile);
    fclose(fDatafile);
    _endFunction();
}

// compare function for c++'s quick-sort
int compare(const void *a, const void *b)
{
    float fDiff = ((EDGE*)a)->capriciousity - ((EDGE*)b)->capriciousity;
    if (fDiff < 0)      return -1;
    else if (fDiff > 0) return  1;
    else                return  0;
}
void quicker_sort(EDGE *edge, int iNumEdges)
{
    _beginFunction("Sorting edges with c++'s qsort");
    qsort(edge, iNumEdges, sizeof(EDGE), compare);
    _endFunction();
}

// voxels are initialized as alone-in group with very high capriciousity
void initialiseVOXELs(float *pdWrappedVolume, VOXEL *voxel, DIMS pdDimensions)
{
    _beginFunction("Initializing all voxels");
    for (int i = 0; i < pdDimensions.size; i++)
    {
        voxel->increment = 0;                   // voxels is not wrapped with respect to itself
        voxel->number_of_voxels_in_group = 1;   // group only contains voxel itself
        voxel->value = *(pdWrappedVolume+i);    // voxels value
        voxel->capriciousity = 9999999+rand();  // very high capriciousity (usual values are arround 1)
        voxel->head = voxel;
        voxel->last = voxel;
        voxel->next = NULL;
        voxel++;
    }
    _endFunction();
}

int find_wrap(float dDelta)
{
    if      (dDelta >  PI)  return -1;
    else if (dDelta < -PI)  return  1;
    else                    return  0;
}

// This function was replaced with a macro for performance reasons
#define wrap(dValue) ((dValue) + TWOPI * find_wrap(dValue))
/*float wrap(float dValue)
{
  //return dValue + find_wrap(dValue) * TWOPI;
  if      (dValue >  PI)  return dValue - TWOPI;
  else if (dValue < -PI)  return dValue + TWOPI;
  else                    return dValue;
}*/

// this function-macro-combination was replaced with one macro for performance reasons
#define getValue(w, h, d) (*(pfVoxelBase + ((d) * pdDimensions.h + (h)) * pdDimensions.w + (w)))
/*define getValue(w, h, d) getValue_(w, h, d, pdDimensions, pfVoxelBase)
float getValue_(int w, int h, int d, DIMS pdDimensions, float *pfVoxelBase)
{
    return *(pfVoxelBase + (d * pdDimensions.h + h) * pdDimensions.w + w);
}*/

#define getVoxel(x,y,z) (voxel_base + ((z) * pdDimensions.h + (y)) * pdDimensions.w + (x))

// is calculated by the sum over squares of all difference of the gradients in in the same direction at two opposite sides of the voxel
// see attached animation capriciousity.gif
void calculate_capriciousity(float *pfVoxelBase, VOXEL *voxel_base, DIMS pdDimensions)
{
    _beginFunction("Calculating the capriciousity");
    float dTmp, dCapriciousity, dCenter;
    int w, h, d;

    for (d = 1; d < pdDimensions.d - 1; d++)
    {
        for (h = 1; h < pdDimensions.h - 1; h++)
        {
            for (w = 1; w < pdDimensions.w - 1; w++)
            {
                dCenter = getValue(w, h, d);
                dCapriciousity = 0;

                dTmp = wrap(getValue(w-1, h  , d  ) - dCenter) - wrap(dCenter - getValue(w+1, h  , d  )); dCapriciousity += dTmp * dTmp;

                dTmp = wrap(getValue(w+1, h-1, d  ) - dCenter) - wrap(dCenter - getValue(w-1, h+1, d  )); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w  , h-1, d  ) - dCenter) - wrap(dCenter - getValue(w  , h+1, d  )); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w-1, h-1, d  ) - dCenter) - wrap(dCenter - getValue(w+1, h+1, d  )); dCapriciousity += dTmp * dTmp;

                dTmp = wrap(getValue(w+1, h+1, d-1) - dCenter) - wrap(dCenter - getValue(w-1, h-1, d+1)); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w  , h+1, d-1) - dCenter) - wrap(dCenter - getValue(w  , h-1, d+1)); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w-1, h+1, d-1) - dCenter) - wrap(dCenter - getValue(w+1, h-1, d+1)); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w+1, h  , d-1) - dCenter) - wrap(dCenter - getValue(w-1, h  , d+1)); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w  , h  , d-1) - dCenter) - wrap(dCenter - getValue(w  , h  , d+1)); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w-1, h  , d-1) - dCenter) - wrap(dCenter - getValue(w+1, h  , d+1)); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w+1, h-1, d-1) - dCenter) - wrap(dCenter - getValue(w-1, h+1, d+1)); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w  , h-1, d-1) - dCenter) - wrap(dCenter - getValue(w  , h+1, d+1)); dCapriciousity += dTmp * dTmp;
                dTmp = wrap(getValue(w-1, h-1, d-1) - dCenter) - wrap(dCenter - getValue(w+1, h+1, d+1)); dCapriciousity += dTmp * dTmp;

                getVoxel(w, h, d)->capriciousity = dCapriciousity;
            }
        }
    }
    _endFunction();
}

// remember: #define getVoxel(x,y,z) (voxel_base + ((z) * pdDimensions.h + (y)) * pdDimensions.w + (x))
// here a connection between directly adjacent voxels is created and the capriciousity of this edge is
// set to the sum of both connecting voxel's capriciousity values
void creatingEDGEs(VOXEL *voxel_base, EDGE *edge, DIMS pdDimensions, int &iNumEdges)
{
    _beginFunction("Creating all edges");
    int d, h, w, i;
    VOXEL *activeVoxel, *nextVoxel;

    for (d = 0; d < pdDimensions.d; d++)
    {
        for (h = 0; h < pdDimensions.h; h++)
        {
            for (w = 0; w < pdDimensions.w; w++)
            {
                activeVoxel = getVoxel(w, h, d);
                for (i = 0; i < 3; i++)
                {
                    switch (i)
                    {
                        case 0: { if (w != pdDimensions.w - 1) nextVoxel = getVoxel(w + 1, h, d); else continue; break; }
                        case 1: { if (h != pdDimensions.h - 1) nextVoxel = getVoxel(w, h + 1, d); else continue; break; }
                        case 2: { if (d != pdDimensions.d - 1) nextVoxel = getVoxel(w, h, d + 1); else continue; break; }
                    }
                    edge->pointer_1 = activeVoxel;
                    edge->pointer_2 = nextVoxel;
                    edge->capriciousity = activeVoxel->capriciousity + nextVoxel->capriciousity;
                    edge->increment = find_wrap(activeVoxel->value - nextVoxel->value);
                    edge++;
                    iNumEdges++;
                }
            }
        }
    }
    _endFunction();
}

// in ascending order every edge is processed and connecting voxels are grouped together if not already done
// if a voxel is still in its initial only self-containing group it joins the others voxel group
// if both voxels are in a real group the smaller group is attached to the bigger one
// The number of wraps between a voxel and its new group is modified by the wrap between the wrap between the voxels
// connected by the active edge and the thereby connected groups of voxels
void gatherVOXELs(EDGE *edge, int iNumEdges)
{
    _beginFunction("Gathering voxels");
    VOXEL *voxel1, *voxel2;
    VOXEL *group1, *group2;
    int i, iIncrement;

    for (i = 0; i < iNumEdges; i++)
    {
        voxel1 = edge->pointer_1;
        voxel2 = edge->pointer_2;

        // If both voxels are in the same group
        if (voxel2->head != voxel1->head)
        {
            if (voxel2->head->number_of_voxels_in_group == 1) // If voxel2 is (in) a single group
            {
                group1 = voxel1->head;
                group1->last->next = voxel2;
                group1->last = voxel2;
                (group1->number_of_voxels_in_group)++;
                voxel2->head = group1;
                voxel2->increment = voxel1->increment - edge->increment;
            }
            else if (voxel1->head->number_of_voxels_in_group == 1) // If voxel1 is (in) a single group
            {
                group2 = voxel2->head;
                group2->last->next = voxel1;
                group2->last = voxel1;
                (group2->number_of_voxels_in_group)++;
                voxel1->head = group2;
                voxel1->increment = voxel2->increment + edge->increment;
            }
            else // Else merge both groups
            {
                if (voxel1->head->number_of_voxels_in_group > voxel2->head->number_of_voxels_in_group) // choose smaller group
                {
                    group1 = voxel1->head;
                    group2 = voxel2->head;
                    iIncrement = voxel1->increment - voxel2->increment - edge->increment;
                }
                else
                {
                    group1 = voxel2->head;
                    group2 = voxel1->head;
                    iIncrement = voxel2->increment - voxel1->increment + edge->increment;
                }
                group1->last->next = group2;
                group1->last = group2->last;
                group1->number_of_voxels_in_group = group1->number_of_voxels_in_group + group2->number_of_voxels_in_group;
                while (group2 != NULL)
                {
                    group2->head = group1;
                    group2->increment += iIncrement;
                    group2 = group2->next;
                }
            }
        }
        edge++;
    }
    _endFunction();
}

// every voxel is unwrapped by the number of wraps it has with respect to the one all voxels containing group
void unwrapVolume(VOXEL *voxel, int iNumberElements)
{
    _beginFunction("Unwrapping the volumes voxels");
    for (int i = 0; i < iNumberElements; i++)
    {
        (voxel+i)->value += (voxel+i)->increment * TWOPI;
    }
    _endFunction();
}

// returns the quality map build from the capriciousity value
void returnCapriciousity(VOXEL *voxel, float *pdQualityMap, int iNumberElements)
{
  for (int i=0; i < iNumberElements; i++)
  {
    *(pdQualityMap+i) = (voxel+i)->capriciousity;
  }
}

void returnVolume(VOXEL *voxel, float *pdUnwrappedVolume, int iNumberElements)
{
    for (int i=0; i < iNumberElements; i++)
    {
      *(pdUnwrappedVolume+i) = (voxel+i)->value;
    }
}

int memory_check(long int elements)
{
  int pointersize = 8 * sizeof(void*);
  char prefix = 'G'; // Gigabyte
  if (pointersize == 64)
	{
	  prefix = 'E'; // Exabyte
	}
  long unsigned int adresses = pow(2, pointersize-1);
  printf("This %dbit executable's memory size is limited to %d%cB, ", pointersize, (int)pow(2, pointersize % 10), prefix);
  long unsigned int memory   = (long int)elements*(4+28+3*16+28+4+4/2);
  printf("the selected volume requires %.1LfGB.\n", (long double)(memory/pow(1024,3)*2));
  long int remaining = adresses - memory;
  if (remaining < 0)
	{
	  printf("Unwrapping cannot be performed by this program due to insufficient memory. ");
	  printf("Please use another version%s or unwrapp a smaller (sub)volume.\n", pointersize == 32 ? " (e.g. 64bit)" : "");
	  return -1;
	}
  if (remaining < 100 * (int)pow(1024, 2))
	{
	  printf("Warning: remaining memory may not be sufficient");
	}
  return 0;
}

int main(int argc, char* argv[])
{
  time_t tStart, tEnd;
  int iDatatype, iNumEdges = 0;
  DIMS dimensions;
  char fWrapped[255] = "Phase.img";
  char fUnwrapped[255] = "uw";
  char fQualityMap[255] = "qm";
  
  time (&tStart);
  
  if (argc > 1)
    {
	  if (strcmp(argv[1], "-v") == 0) {
		printf("This is version %s of the HusseinUnwrapper (%s)\n", VERSION, REVISION);
		exit(0);
	  }
	  strcpy(fWrapped, argv[1]);
	  if (argc > 2)
        {
		  bPrint = strcmp("silent", argv[2]);
        }
    }
  strcat(fUnwrapped, fWrapped);
  strcat(fQualityMap, fWrapped);
  read_header(fWrapped, &dimensions, iDatatype);
  
  if (memory_check((long int)dimensions.size) != 0)
	{
	  return -1;
	}
  
  float *pdWrappedVolume   = new float[dimensions.size];
  VOXEL *voxel             = new VOXEL[dimensions.size];
  EDGE  *edge              = new  EDGE[dimensions.size*3];
  float *pdQualityMap      = new float[dimensions.size];
  float *pdUnwrappedVolume = new float[dimensions.size];
  
  read_data(fWrapped, pdWrappedVolume, dimensions.size, iDatatype);
  
  initialiseVOXELs(pdWrappedVolume, voxel, dimensions);
  calculate_capriciousity(pdWrappedVolume, voxel, dimensions);
  returnCapriciousity(voxel, pdQualityMap, dimensions.size);
  write_data(fQualityMap, pdQualityMap, dimensions.size, "quality map");
  write_header(fQualityMap);
  
  creatingEDGEs(voxel, edge, dimensions, iNumEdges);
  quicker_sort(edge, iNumEdges);
  gatherVOXELs(edge, iNumEdges);
  
  unwrapVolume(voxel, dimensions.size);
  returnVolume(voxel, pdUnwrappedVolume, dimensions.size);
  
  write_data(fUnwrapped, pdUnwrappedVolume, dimensions.size, "unwrapped values");
  write_header(fUnwrapped);
  
  free(pdWrappedVolume);
  free(voxel);
  free(edge);
  free(pdUnwrappedVolume);
  
  time (&tEnd);
  if (bPrint)
    {
      printf ("Unwrapped \"%s\" into \"%s\" in  %.1f seconds\n", fWrapped, fUnwrapped, difftime(tEnd,tStart));
    }
  return 1;
}
