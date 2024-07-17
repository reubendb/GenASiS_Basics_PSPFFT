#include "SiloWrappers.h"

void SiloWrappers_MakeDirectory(const char *Pathname, int *Error)
  {
  int error;
  
  if(mkdir(Pathname, S_IRWXU) == -1)
    {
    if(errno /= EEXIST) *Error = errno;
    }
  }


void SiloWrappers_DBGetToc(
       const int *SiloHandle, const int *Rank, int *FileKey)
  {
  int 
    iCount;
  char 
    FileNameOut[PATH_MAX];
  FILE 
    *Handle;
  DBtoc
    *SiloToc;
  DBfile 
     *SiloFile;
  
  srand(time(0));
  *FileKey = rand();
  // Use random number between 1 and 999999 
  *FileKey = 1 + (int) (999999.0 * rand()/ ( RAND_MAX + 1.0 ));
  
  sprintf(FileNameOut, "._%07d_%07d.silotoc", *FileKey, *Rank);
  
  SiloFile = (DBfile *) DBFortranAccessPointer(*SiloHandle);
  SiloToc  = (DBtoc *) DBGetToc(SiloFile);
  
  Handle = fopen(FileNameOut, "w");
  if(Handle == NULL)
    {
    fprintf(stderr, "ERROR: Opening %s for write", FileNameOut);
    perror("Silo_C_DBGetToc:");
    return;
    }

  // Write the TOC
  fprintf(Handle, "ncurve=%d\n", SiloToc->ncurve);
  fprintf(Handle, "curve_names=");
  for(iCount=0; iCount < SiloToc->ncurve-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->curve_names[iCount]);
  if(SiloToc->ncurve > iCount) 
    fprintf(Handle, "%s", SiloToc->curve_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nmultimesh=%d\n", SiloToc->nmultimesh);
  fprintf(Handle, "multimesh_names=");
  for(iCount=0; iCount < SiloToc->nmultimesh-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->multimesh_names[iCount]);
  if(SiloToc->nmultimesh > iCount)
    fprintf(Handle, "%s", SiloToc->multimesh_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nmultimeshadj=%d\n", SiloToc->nmultimeshadj);
  fprintf(Handle, "multimeshadj_names=");
  for(iCount=0; iCount < SiloToc->nmultimeshadj-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->multimeshadj_names[iCount]);
  if(SiloToc->nmultimeshadj > iCount)
    fprintf(Handle, "%s", SiloToc->multimeshadj_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nmultivar=%d\n", SiloToc->nmultivar);
  fprintf(Handle, "multivar_names=");
  for(iCount=0; iCount < SiloToc->nmultivar-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->multivar_names[iCount]);
  if(SiloToc->nmultivar > iCount)
    fprintf(Handle, "%s", SiloToc->multivar_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nmultimat=%d\n", SiloToc->nmultimat);
  fprintf(Handle, "multimat_names=");
  for(iCount=0; iCount < SiloToc->nmultimat-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->multimat_names[iCount]);
  if(SiloToc->nmultimat > iCount)
    fprintf(Handle, "%s", SiloToc->multimat_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nmultimatspecies=%d\n", SiloToc->nmultimatspecies);
  fprintf(Handle, "multimatspecies_names=");
  for(iCount=0; iCount < SiloToc->nmultimatspecies-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->multimatspecies_names[iCount]);
  if(SiloToc->nmultimatspecies > iCount)
    fprintf(Handle, "%s", SiloToc->multimatspecies_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "ncsgmesh=%d\n", SiloToc->ncsgmesh);
  fprintf(Handle, "csgmesh_names=");
  for(iCount=0; iCount < SiloToc->ncsgmesh-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->csgmesh_names[iCount]);
  if(SiloToc->ncsgmesh > iCount)
    fprintf(Handle, "%s", SiloToc->csgmesh_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "ncsgvar=%d\n", SiloToc->ncsgvar);
  fprintf(Handle, "csgvar_names=");
  for(iCount=0; iCount < SiloToc->ncsgvar-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->csgvar_names[iCount]);
  if(SiloToc->ncsgvar > iCount)
    fprintf(Handle, "%s", SiloToc->csgvar_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "ndefvars=%d\n", SiloToc->ndefvars);
  fprintf(Handle, "defvars_names=");
  for(iCount=0; iCount < SiloToc->ndefvars-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->defvars_names[iCount]);
  if(SiloToc->ndefvars > iCount)
    fprintf(Handle, "%s", SiloToc->defvars_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nqmesh=%d\n", SiloToc->nqmesh);
  fprintf(Handle, "qmesh_names=");
  for(iCount=0; iCount < SiloToc->nqmesh-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->qmesh_names[iCount]);
  if(SiloToc->nqmesh > iCount)
    fprintf(Handle, "%s", SiloToc->qmesh_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nqvar=%d\n", SiloToc->nqvar);
  fprintf(Handle, "qvar_names=");
  for(iCount=0; iCount < SiloToc->nqvar-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->qvar_names[iCount]);
  if(SiloToc->nqvar > iCount)
    fprintf(Handle, "%s", SiloToc->qvar_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nucdmesh=%d\n", SiloToc->nucdmesh);
  fprintf(Handle, "ucdmesh_names=");
  for(iCount=0; iCount < SiloToc->nucdmesh-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->ucdmesh_names[iCount]);
  if(SiloToc->nucdmesh > iCount)
    fprintf(Handle, "%s", SiloToc->ucdmesh_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nucdvar=%d\n", SiloToc->nucdvar);
  fprintf(Handle, "ucdvar_names=");
  for(iCount=0; iCount < SiloToc->nucdvar-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->ucdvar_names[iCount]);
  if(SiloToc->nucdvar > iCount)
    fprintf(Handle, "%s", SiloToc->ucdvar_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nptmesh=%d\n", SiloToc->nptmesh);
  fprintf(Handle, "ptmesh_names=");
  for(iCount=0; iCount < SiloToc->nptmesh-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->ptmesh_names[iCount]);
  if(SiloToc->nptmesh > iCount)
    fprintf(Handle, "%s", SiloToc->ptmesh_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nptvar=%d\n", SiloToc->nptvar);
  fprintf(Handle, "ptvar_names=");
  for(iCount=0; iCount < SiloToc->nptvar-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->ptvar_names[iCount]);
  if(SiloToc->nptvar > iCount)
    fprintf(Handle, "%s", SiloToc->ptvar_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nmat=%d\n", SiloToc->nmat);
  fprintf(Handle, "mat_names=");
  for(iCount=0; iCount < SiloToc->nmat-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->mat_names[iCount]);
  if(SiloToc->nmat > iCount)
    fprintf(Handle, "%s", SiloToc->mat_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nmatspecies=%d\n", SiloToc->nmatspecies);
  fprintf(Handle, "matspecies_names=");
  for(iCount=0; iCount < SiloToc->nmatspecies-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->matspecies_names[iCount]);
  if(SiloToc->nmatspecies > iCount)
    fprintf(Handle, "%s", SiloToc->matspecies_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nvar=%d\n", SiloToc->nvar);
  fprintf(Handle, "var_names=");
  for(iCount=0; iCount < SiloToc->nvar-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->var_names[iCount]);
  if(SiloToc->nvar > iCount)
    fprintf(Handle, "%s", SiloToc->var_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "nobj=%d\n", SiloToc->nobj);
  fprintf(Handle, "obj_names=");
  for(iCount=0; iCount < SiloToc->nobj-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->obj_names[iCount]);
  if(SiloToc->nobj > iCount)
    fprintf(Handle, "%s", SiloToc->obj_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "ndir=%d\n", SiloToc->ndir);
  fprintf(Handle, "dir_names=");
  for(iCount=0; iCount < SiloToc->ndir-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->dir_names[iCount]);
  if(SiloToc->ndir > iCount)
    fprintf(Handle, "%s", SiloToc->dir_names[iCount]);
  fprintf(Handle, "\n");
  
  fprintf(Handle, "narrays=%d\n", SiloToc->narray);
  fprintf(Handle, "array_names=");
  for(iCount=0; iCount < SiloToc->narray-1; iCount++)
    fprintf(Handle, "%s,", SiloToc->array_names[iCount]);
  if(SiloToc->narray > iCount)
    fprintf(Handle, "%s", SiloToc->array_names[iCount]);
  fprintf(Handle, "\n");
  
  fclose(Handle);
  }
  

void SiloWrappers_DBFreeToc(const int *FileKey, const int *Rank)
  {
  char FileName[PATH_MAX];
  
  sprintf(FileName, "._%07d_%07d.silotoc", *FileKey, *Rank);
  remove(FileName);
  }


void SiloWrappers_DBGetHeader(const int *SiloHandle, double *Time, int *Cycle)
  {
  DBfile
    *Handle;
  DBquadmesh 
    *Header;
  
  Handle = (DBfile *) DBFortranAccessPointer(*SiloHandle);
  Header = (DBquadmesh *) DBGetQuadmesh(Handle, "Header");
  *Time   = Header->dtime;
  *Cycle  = Header->cycle;
  
  DBFreeQuadmesh(Header);
  }
  

void SiloWrappers_DBGetQuadMeshMtdt(const int *SiloHandle, int *nNodes)
  {
  DBfile
    *Handle;
  DBquadmesh 
    *Mesh;
  
  Handle    = (DBfile *) DBFortranAccessPointer(*SiloHandle);
  Mesh      = (DBquadmesh *) DBGetQuadmesh(Handle, "Mesh");
  nNodes[0] = Mesh->dims[0];
  nNodes[1] = Mesh->dims[1];
  nNodes[2] = Mesh->dims[2];
  
  DBFreeQuadmesh(Mesh);
  }
  

void SiloWrappers_DBGetQuadMeshCrdnt(
       const int *SiloHandle, double *NodeCoordinate_1, 
       double *NodeCoordinate_2, double *NodeCoordinate_3)
      
  {
  int
    iNode;
  double 
    *NodeCoordinateScratch_1, 
    *NodeCoordinateScratch_2, 
    *NodeCoordinateScratch_3; 
  DBfile
    *Handle;
  DBquadmesh 
    *Mesh;
  
  Handle = (DBfile *) DBFortranAccessPointer(*SiloHandle);
  Mesh   = (DBquadmesh *) DBGetQuadmesh(Handle, "Mesh");
  
  NodeCoordinateScratch_1 = (double *) Mesh->coords[0];
  for(iNode = 0; iNode < Mesh->dims[0]; iNode++)
    NodeCoordinate_1[iNode] = NodeCoordinateScratch_1[iNode];
  if(Mesh->ndims > 1)
    {
    NodeCoordinateScratch_2 = (double *) Mesh->coords[1];
    for(iNode = 0; iNode < Mesh->dims[1]; iNode++)
      NodeCoordinate_2[iNode] = NodeCoordinateScratch_2[iNode];
    }
  if(Mesh->ndims > 2)
    {
    NodeCoordinateScratch_3 = (double *) Mesh->coords[2];
    for(iNode = 0; iNode < Mesh->dims[2]; iNode++)
      NodeCoordinate_3[iNode] = NodeCoordinateScratch_3[iNode];
    }
  
  DBFreeQuadmesh(Mesh);
  }
  
void SiloWrappers_DBGetQuadMeshVrbl(
       const int *SiloHandle, const char *VariableName, double *Data)
  {
  int
    iData,
    nData;
  double 
    *DataScratch;
  DBfile
    *Handle;
  DBquadvar 
    *Variable;
    
  Handle   = (DBfile *) DBFortranAccessPointer(*SiloHandle);
  Variable = (DBquadvar *) DBGetQuadvar(Handle, VariableName);
  
  DataScratch = (double *) Variable->vals[0];
  nData = Variable->nels;
  for(iData = 0; iData < nData; iData++)
    Data[iData] = DataScratch[iData];
    
  DBFreeQuadvar(Variable);
  }
       

void SiloWrappers_DBGetUcdMeshMtdt(
       const int *SiloHandle, int *nDimensions, int *nNodes, int *nTotalCells, 
       int *nGhostCells)
  {
  DBfile
    *Handle;
  DBucdmesh 
    *Mesh;
  
  Handle = (DBfile *) DBFortranAccessPointer(*SiloHandle);
  Mesh = (DBucdmesh *) DBGetUcdmesh(Handle, "Mesh");
  
  if(Mesh == NULL)
    {
    *nDimensions = 1;
    return;
    }
  
  *nDimensions 	= Mesh->ndims;
  *nNodes       = Mesh->nnodes;
  *nTotalCells  = Mesh->zones->nzones;
  *nGhostCells
    = *nTotalCells - (Mesh->zones->max_index - Mesh->zones->min_index + 1);
    
  DBFreeUcdmesh(Mesh);
  }


void SiloWrappers_DBGetUcdMeshCrdnt(
       const int *SiloHandle, const int *nDimensions, 
       double *NodeCoordinate_1, double *NodeCoordinate_2, 
       double *NodeCoordinate_3)
  {
  int
    iNode;
  double 
    *NodeCoordinateScratch_1, 
    *NodeCoordinateScratch_2, 
    *NodeCoordinateScratch_3; 
  DBfile
    *Handle;
  DBucdmesh 
    *Mesh;
  
  if(*nDimensions == 1) return;
  
  Handle = (DBfile *) DBFortranAccessPointer(*SiloHandle);
  Mesh = (DBucdmesh *) DBGetUcdmesh(Handle, "Mesh");
  
  if(Mesh == NULL)
    return;
  
  // FIXME: Unit is not yet supported for reading
  if(*nDimensions == 2)
    {
    NodeCoordinateScratch_1 = (double *) Mesh->coords[0];
    NodeCoordinateScratch_2 = (double *) Mesh->coords[1];
    for(iNode = 0; iNode < Mesh->nnodes; iNode++)
      {
      NodeCoordinate_1[iNode] = NodeCoordinateScratch_1[iNode];
      NodeCoordinate_2[iNode] = NodeCoordinateScratch_2[iNode];
      } 
    }
  else if(*nDimensions == 3)
    {
    NodeCoordinateScratch_1 = (double *) Mesh->coords[0];
    NodeCoordinateScratch_2 = (double *) Mesh->coords[1];
    NodeCoordinateScratch_3 = (double *) Mesh->coords[2];
    for(iNode = 0; iNode < Mesh->nnodes; iNode++)
      {
      NodeCoordinate_1[iNode] = NodeCoordinateScratch_1[iNode];
      NodeCoordinate_2[iNode] = NodeCoordinateScratch_2[iNode];
      NodeCoordinate_3[iNode] = NodeCoordinateScratch_3[iNode];
      } 
    }
  
  DBFreeUcdmesh(Mesh);
  }


void SiloWrappers_DBGetUcdMeshVrbl(
       const int *SiloHandle, const int *nProperCells, 
       const char *VariableName, double *Data, int *Error)
  {
  int
    iData, 
    iZone;
  double 
    *Scratch;
  DBfile
    *Handle;
  DBucdmesh 
    *Mesh;
  DBucdvar
    *Variable;
  
  *Error = 0;
  
  Handle = (DBfile *) DBFortranAccessPointer(*SiloHandle);
  Variable = (DBucdvar *) DBGetUcdvar(Handle, VariableName);
  Mesh = (DBucdmesh *) DBGetUcdmesh(Handle, Variable->meshname);
  assert(Variable->nvals == 1);
  Scratch = (double *) Variable->vals[0];
  iData = 0;
  for(iZone=Mesh->zones->min_index; iZone <= Mesh->zones->max_index; iZone++)
    {
    Data[iData] = Scratch[iZone];
    iData++;
    }
  if(iData != *nProperCells) *Error = -1;
  
  DBFreeUcdmesh(Mesh);
  DBFreeUcdvar(Variable);
  }
  
  
void SiloWrappers_DBGetNumberOfBlcks(const int *SiloHandle, int *nBlocks)
  {
  DBfile
    *Handle;
  DBmultimesh
    *Multimesh;
  
  Handle = (DBfile *) DBFortranAccessPointer(*SiloHandle);
  
  // Assume "Level_01/Interior" always exists
  Multimesh = (DBmultimesh *) DBGetMultimesh(Handle, "Level_01/Interior/Mesh");
  *nBlocks = Multimesh->nblocks;
  
  DBFreeMultimesh(Multimesh);
  }
