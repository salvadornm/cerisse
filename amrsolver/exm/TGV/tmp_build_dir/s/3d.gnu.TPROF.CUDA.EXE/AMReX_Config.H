#ifndef AMREX_HAVE_NO_CONFIG_H
#define AMREX_HAVE_NO_CONFIG_H
#ifndef AMREX_EXPORT_DYNAMIC
#define AMREX_EXPORT_DYNAMIC 1
#endif
#ifndef AMREX_TINY_PROFILING
#define AMREX_TINY_PROFILING 1
#endif
#ifndef BL_LAZY
#define BL_LAZY 1
#endif
#ifndef AMREX_LAZY
#define AMREX_LAZY 1
#endif
#ifndef AMREX_USE_CUDA
#define AMREX_USE_CUDA 1
#endif
#ifndef AMREX_USE_GPU
#define AMREX_USE_GPU 1
#endif
#ifndef BL_COALESCE_FABS
#define BL_COALESCE_FABS 1
#endif
#ifndef AMREX_GPU_MAX_THREADS
#define AMREX_GPU_MAX_THREADS 256
#endif
#ifndef AMREX_USE_GPU_RDC
#define AMREX_USE_GPU_RDC 1
#endif
#ifndef BL_SPACEDIM
#define BL_SPACEDIM 3
#endif
#ifndef AMREX_SPACEDIM
#define AMREX_SPACEDIM 3
#endif
#ifndef BL_FORT_USE_UNDERSCORE
#define BL_FORT_USE_UNDERSCORE 1
#endif
#ifndef AMREX_FORT_USE_UNDERSCORE
#define AMREX_FORT_USE_UNDERSCORE 1
#endif
#ifndef BL_Darwin
#define BL_Darwin 1
#endif
#ifndef AMREX_Darwin
#define AMREX_Darwin 1
#endif
#ifndef NDEBUG
#define NDEBUG 1
#endif
#ifndef AMREX_NO_PROBINIT
#define AMREX_NO_PROBINIT 1
#endif
#ifndef AMREX_CRSEGRNDOMP
#define AMREX_CRSEGRNDOMP 1
#endif
#ifdef __cplusplus
#endif
#if defined(AMREX_USE_OMP) && !defined(_OPENMP) && !defined(AMREX_USE_HIP)
#error libamrex was built with OpenMP
#endif
#endif
