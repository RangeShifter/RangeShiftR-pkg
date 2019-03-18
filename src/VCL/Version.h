//---------------------------------------------------------------------------

#ifndef VersionH
#define VersionH

#define RSDEBUG 0

#define CLUSTER	0
#define RSWIN64 0
#if RSWIN64
#define RSRANDOM 1  // must be set to 1 if RSWIN64 is 1
#else
#define RSRANDOM 0  // optional
#endif

#define BATCH 0
#define VCL 1

//---------------------------------------------------------------------------
#endif
