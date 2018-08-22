//---------------------------------------------------------------------------

#ifndef VersionH
#define VersionH

#define RSDEBUG 1

#define CLUSTER	0
#define RSWIN64 0
#if RSWIN64
#define RSRANDOM 1  // must be set to 1 if RSWIN64 is 1
#else
#define RSRANDOM 0  // optional
#endif
#define RANDOMCHECK 0
#define BATCH 1
#define VCL 0

//---------------------------------------------------------------------------
#endif
