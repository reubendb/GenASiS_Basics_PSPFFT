
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <libgen.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include "silo.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#ifndef _SILO_WRAPPERS_H
#define _SILO_WRAPPERS_H

#define BUFFER_SIZE 1024

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#ifdef __cplusplus   
extern "C" {
#endif
  
  /* MakeDirectory */
  void SiloWrappers_MakeDirectory(const char *Pathname, int *Error);
  void SiloWrappers_MakeDirectory_(const char *Pathname, int *Error);
  void silowrappers_makedirectory(const char *Pathname, int *Error);
  void silowrappers_makedirectory_(const char *Pathname, int *Error);
  
  void SiloWrappers_MakeDirectory_(const char *Pathname, int *Error)
    {
    SiloWrappers_MakeDirectory(Pathname, Error);
    }
  void silowrappers_makedirectory(const char *Pathname, int *Error)
    {
    SiloWrappers_MakeDirectory(Pathname, Error);
    }
  void silowrappers_makedirectory_(const char *Pathname, int *Error)
    {
    SiloWrappers_MakeDirectory(Pathname, Error);
    }

    
  /* DBGetToc */
  void SiloWrappers_DBGetToc( 
         const int *SiloHandle, const int *Rank, int *FileKey);
  void SiloWrappers_DBGetToc_(
         const int *SiloHandle, const int *Rank, int *FileKey);
  void silowrappers_dbgettoc(
         const int *SiloHandle, const int *Rank, int *FileKey);
  void silowrappers_dbgettoc_(
         const int *SiloHandle, const int *Rank, int *FileKey);
  
  void SiloWrappers_DBGetToc_(
         const int *SiloHandle, const int *Rank, int *FileKey)
    {
    SiloWrappers_DBGetToc(SiloHandle, Rank, FileKey);
    }

  void silowrappers_dbgettoc(
         const int *SiloHandle, const int *Rank, int *FileKey)
    {
    SiloWrappers_DBGetToc(SiloHandle, Rank, FileKey);
    }

  void silowrappers_dbgettoc_(
         const int *SiloHandle, const int *Rank, int *FileKey)
    {
    SiloWrappers_DBGetToc(SiloHandle, Rank, FileKey);
    }
  
  
  /* DBFreeToc */
  void SiloWrappers_DBFreeToc(const int *FileKey, const int *Rank);
  void SiloWrappers_DBFreeToc_(const int *FileKey, const int *Rank);
  void silowrappers_dbfreetoc(const int *FileKey, const int *Rank);
  void silowrappers_dbfreetoc_(const int *FileKey, const int *Rank);
  
  void SiloWrappers_DBFreeToc_(const int *FileKey, const int *Rank)
    {
    SiloWrappers_DBFreeToc(FileKey, Rank);
    }
  void silowrappers_dbfreetoc(const int *FileKey, const int *Rank)
    {
    SiloWrappers_DBFreeToc(FileKey, Rank);
    }
  void silowrappers_dbfreetoc_(const int *FileKey, const int *Rank)
    {
    SiloWrappers_DBFreeToc(FileKey, Rank);
    }
    
  
  /* DBGetHeader */
  void SiloWrappers_DBGetHeader(
         const int *SiloHandle, double *Time, int *Cycle);
  void SiloWrappers_DBGetHeader_(
         const int *SiloHandle, double *Time, int *Cycle);
  void silowrappers_dbgetheader(
         const int *SiloHandle, double *Time, int *Cycle);
  void silowrappers_dbgetheader_(
         const int *SiloHandle, double *Time, int *Cycle);
  
  void SiloWrappers_DBGetHeader_(
         const int *SiloHandle, double *Time, int *Cycle)
    {
    SiloWrappers_DBGetHeader(SiloHandle, Time, Cycle);
    }
  void silowrappers_dbgetheader(
         const int *SiloHandle, double *Time, int *Cycle)
    {
    SiloWrappers_DBGetHeader(SiloHandle, Time, Cycle);
    }
  void silowrappers_dbgetheader_(
         const int *SiloHandle, double *Time, int *Cycle)
    {
    SiloWrappers_DBGetHeader(SiloHandle, Time, Cycle);
    }
  
  
  /* DBGetQuadMeshMetadata */
  void SiloWrappers_DBGetQuadMeshMtdt(const int *SiloHandle, int *nNodes);
  void SiloWrappers_DBGetQuadMeshMtdt_(const int *SiloHandle, int *nNodes);
  void silowrappers_dbgetquadmeshmtdt(const int *SiloHandle, int *nNodes);
  void silowrappers_dbgetquadmeshmtdt_(const int *SiloHandle, int *nNodes);
  
  void SiloWrappers_DBGetQuadMeshMtdt_(const int *SiloHandle, int *nNodes)
    {
    SiloWrappers_DBGetQuadMeshMtdt(SiloHandle, nNodes);
    }
  void silowrappers_dbgetquadmeshmtdt(const int *SiloHandle, int *nNodes)
    {
    SiloWrappers_DBGetQuadMeshMtdt(SiloHandle, nNodes);
    }
  void silowrappers_dbgetquadmeshmtdt_(const int *SiloHandle, int *nNodes)
    {
    SiloWrappers_DBGetQuadMeshMtdt(SiloHandle, nNodes);
    }
    
    
  /* DBGetQuadMeshCoordinate */
  void SiloWrappers_DBGetQuadMeshCrdnt(
         const int *SiloHandle, double *NodeCoordinate_1, 
         double *NodeCoordinate_2, double *NodeCoordinate_3);
  void SiloWrappers_DBGetQuadMeshCrdnt_(
         const int *SiloHandle, double *NodeCoordinate_1, 
         double *NodeCoordinate_2, double *NodeCoordinate_3);
  void silowrappers_dbgetquadmeshcrdnt(
         const int *SiloHandle, double *NodeCoordinate_1, 
         double *NodeCoordinate_2, double *NodeCoordinate_3);
  void silowrappers_dbgetquadmeshcrdnt_(
         const int *SiloHandle, double *NodeCoordinate_1, 
         double *NodeCoordinate_2, double *NodeCoordinate_3);
         
  void SiloWrappers_DBGetQuadMeshCrdnt_(
         const int *SiloHandle, double *NodeCoordinate_1, 
         double *NodeCoordinate_2, double *NodeCoordinate_3)
    {
    SiloWrappers_DBGetQuadMeshCrdnt(
      SiloHandle,NodeCoordinate_1, NodeCoordinate_2, NodeCoordinate_3);
    }
  void silowrappers_dbgetquadmeshcrdnt(
         const int *SiloHandle, double *NodeCoordinate_1, 
         double *NodeCoordinate_2, double *NodeCoordinate_3)
    {
    SiloWrappers_DBGetQuadMeshCrdnt(
      SiloHandle,NodeCoordinate_1, NodeCoordinate_2, NodeCoordinate_3);
    }
  void silowrappers_dbgetquadmeshcrdnt_(
         const int *SiloHandle, double *NodeCoordinate_1, 
         double *NodeCoordinate_2, double *NodeCoordinate_3)
    {
    SiloWrappers_DBGetQuadMeshCrdnt(
      SiloHandle,NodeCoordinate_1, NodeCoordinate_2, NodeCoordinate_3);
    }
  
  
  /* DBGetQuadMeshVariable */
  void SiloWrappers_DBGetQuadMeshVrbl(
         const int *SiloHandle, const char *VariableName, double *Variable);
  void SiloWrappers_DBGetQuadMeshVrbl_(
         const int *SiloHandle, const char *VariableName, double *Variable);
  void silowrappers_dbgetquadmeshvrbl(
         const int *SiloHandle, const char *VariableName, double *Variable);
  void silowrappers_dbgetquadmeshvrbl_(
         const int *SiloHandle, const char *VariableName, double *Variable);
         
  void SiloWrappers_DBGetQuadMeshVrbl_(
         const int *SiloHandle, const char *VariableName, double *Variable)
    {
    SiloWrappers_DBGetQuadMeshVrbl(SiloHandle, VariableName, Variable);
    }
  void silowrappers_dbgetquadmeshvrbl(
         const int *SiloHandle, const char *VariableName, double *Variable)
    {
    SiloWrappers_DBGetQuadMeshVrbl(SiloHandle, VariableName, Variable);
    }
  void silowrappers_dbgetquadmeshvrbl_(
         const int *SiloHandle, const char *VariableName, double *Variable)
    {
    SiloWrappers_DBGetQuadMeshVrbl(SiloHandle, VariableName, Variable);
    }
  
  
  /* DBGetUcdMeshMetadata */
  void SiloWrappers_DBGetUcdMeshMtdt(
         const int *SiloHandle, int *nDimensions, int *nNodes, 
         int *nTotalCells, int *nGhostCells);
  void SiloWrappers_DBGetUcdMeshMtdt_(
         const int *SiloHandle, int *nDimensions, int *nNodes, 
         int *nTotalCells, int *nGhostCells);
  void silowrappers_dbgetucdmeshmtdt(
         const int *SiloHandle, int *nDimensions, int *nNodes, 
         int *nTotalCells, int *nGhostCells);
  void silowrappers_dbgetucdmeshmtdt_(
         const int *SiloHandle, int *nDimensions, int *nNodes, 
         int *nTotalCells, int *nGhostCells);
  
  void SiloWrappers_DBGetUcdMeshMtdt_(
         const int *SiloHandle, int *nDimensions, int *nNodes, 
         int *nTotalCells, int *nGhostCells)
    {
    SiloWrappers_DBGetUcdMeshMtdt(
      SiloHandle, nDimensions, nNodes, nTotalCells, nGhostCells);
    }
  void silowrappers_dbgetucdmeshmtdt(
         const int *SiloHandle, int *nDimensions, int *nNodes, int 
         *nTotalCells, int *nGhostCells)
    {
    SiloWrappers_DBGetUcdMeshMtdt(
      SiloHandle, nDimensions, nNodes, nTotalCells, nGhostCells);
    }
  void silowrappers_dbgetucdmeshmtdt_(
         const int *SiloHandle, int *nDimensions, int *nNodes, 
         int *nTotalCells, int *nGhostCells)
    {
    SiloWrappers_DBGetUcdMeshMtdt(
      SiloHandle, nDimensions, nNodes, nTotalCells, nGhostCells);
    }


  /* DBGetUcdMeshCoordinate */
  void SiloWrappers_DBGetUcdMeshCrdnt(
         const int *SiloHandle, const int *nDimensions, 
         double *NodeCoordinate_1, double *NodeCoordinate_2, 
         double *NodeCoordinate_3);
  void SiloWrappers_DBGetUcdMeshCrdnt_(
         const int *SiloHandle, const int *nDimensions, 
         double *NodeCoordinate_1, double *NodeCoordinate_2, 
         double *NodeCoordinate_3);
  void silowrappers_dbgetucdmeshcrdnt(
         const int *SiloHandle, const int *nDimensions, 
         double *NodeCoordinate_1, double *NodeCoordinate_2, 
         double *NodeCoordinate_3);
  void silowrappers_dbgetucdmeshcrdnt_(
         const int *SiloHandle, const int *nDimensions, 
         double *NodeCoordinate_1, double *NodeCoordinate_2, 
         double *NodeCoordinate_3);
  
  void SiloWrappers_DBGetUcdMeshCrdnt_(
         const int *SiloHandle, const int *nDimensions, 
         double *NodeCoordinate_1, double *NodeCoordinate_2, 
         double *NodeCoordinate_3)
    {
    SiloWrappers_DBGetUcdMeshCrdnt(
      SiloHandle, nDimensions, NodeCoordinate_1, NodeCoordinate_2, 
      NodeCoordinate_3);
    }
  void silowrappers_dbgetucdmeshcrdnt(
         const int *SiloHandle, const int *nDimensions, 
         double *NodeCoordinate_1, double *NodeCoordinate_2, 
         double *NodeCoordinate_3)
    {
    SiloWrappers_DBGetUcdMeshCrdnt(
      SiloHandle, nDimensions, NodeCoordinate_1, NodeCoordinate_2, 
      NodeCoordinate_3);
    }
  void silowrappers_dbgetucdmeshcrdnt_(
         const int *SiloHandle, const int *nDimensions, 
         double *NodeCoordinate_1, double *NodeCoordinate_2, 
         double *NodeCoordinate_3)
    {
    SiloWrappers_DBGetUcdMeshCrdnt(
      SiloHandle, nDimensions, NodeCoordinate_1, NodeCoordinate_2, 
      NodeCoordinate_3);
    }
    
  
  /* DBGetUcdMeshVrbl */
  void SiloWrappers_DBGetUcdMeshVrbl(
         const int *SiloHandle, const int *nProperCells, 
         const char *VariableName, double *Data, int *Error);
  void SiloWrappers_DBGetUcdMeshVrbl_(
         const int *SiloHandle, const int *nProperCells, 
         const char *VariableName, double *Data, int *Error);
  void silowrappers_dbgetucdmeshvrbl(
         const int *SiloHandle, const int *nProperCells, 
         const char *VariableName, double *Data, int *Error);
  void silowrappers_dbgetucdmeshvrbl_(
         const int *SiloHandle, const int *nProperCells, 
         const char *VariableName, double *Data, int *Error);
  
  void SiloWrappers_DBGetUcdMeshVrbl_(
         const int *SiloHandle, const int *nProperCells, 
         const char *VariableName, double *Data, int *Error)
    {
    SiloWrappers_DBGetUcdMeshVrbl(
      SiloHandle, nProperCells, VariableName, Data, Error);
    }
  void silowrappers_dbgetucdmeshvrbl(
         const int *SiloHandle, const int *nProperCells, 
         const char *VariableName, double *Data, int *Error)
    {
    SiloWrappers_DBGetUcdMeshVrbl(
      SiloHandle, nProperCells, VariableName, Data, Error);
    }
  void silowrappers_dbgetucdmeshvrbl_(
         const int *SiloHandle, const int *nProperCells, 
         const char *VariableName, double *Data, int *Error)
    {
    SiloWrappers_DBGetUcdMeshVrbl(
      SiloHandle, nProperCells, VariableName, Data, Error);
    }
  
  
  /* DBGetNumberOfBlocks */
  void SiloWrappers_DBGetNumberOfBlcks(const int *SiloHandle, int *nBlocks);
  void SiloWrappers_DBGetNumberOfBlcks_(const int *SiloHandle, int *nBlocks);
  void silowrappers_dbgetnumberofblcks(const int *SiloHandle, int *nBlocks);
  void silowrappers_dbgetnumberofblcks_(const int *SiloHandle, int *nBlocks);
  
  void SiloWrappers_DBGetNumberOfBlcks_(const int *SiloHandle, int *nBlocks)
    {
    SiloWrappers_DBGetNumberOfBlcks(SiloHandle, nBlocks);
    }
  void silowrappers_dbgetnumberofblcks(const int *SiloHandle, int *nBlocks)
    {
    SiloWrappers_DBGetNumberOfBlcks(SiloHandle, nBlocks);
    }
  void silowrappers_dbgetnumberofblcks_(const int *SiloHandle, int *nBlocks)
    {
    SiloWrappers_DBGetNumberOfBlcks(SiloHandle, nBlocks);
    }
       
  
#ifdef __cplusplus
}
#endif    /* __cplusplus */

#endif   /* _SILO_WRAPPERS_H */
