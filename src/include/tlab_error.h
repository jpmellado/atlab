
#ifndef DNS_ERROR_H_INCLUDED
#define DNS_ERROR_H_INCLUDED

#define DNS_ERROR_VERSION     10
#define DNS_ERROR_SIMTYPE     11
#define DNS_ERROR_SIMFLOW     12
#define DNS_ERROR_SIMCOOR     13
#define DNS_ERROR_CALCSCALAR  14
#define DNS_ERROR_RKORDER     15
#define DNS_ERROR_PADEORDER   16
#define DNS_ERROR_UNIFORMX    17
#define DNS_ERROR_UNIFORMY    18
#define DNS_ERROR_UNIFORMZ    19
#define DNS_ERROR_IBC         20
#define DNS_ERROR_JBC         21
#define DNS_ERROR_KBC         22
#define DNS_ERROR_PRESSBC     23
#define DNS_ERROR_IVSICBC     24
#define DNS_ERROR_JVSICBC     25
#define DNS_ERROR_INFTYPE     26
#define DNS_ERROR_INFDISCR    27
#define DNS_ERROR_KMAXMISS    28
#define DNS_ERROR_CHECKUNIFX  29
#define DNS_ERROR_CHECKUNIFY  30
#define DNS_ERROR_CHECKUNIFZ  31
#define DNS_ERROR_CHECKFILTX  32
#define DNS_ERROR_CHECKFILTY  33
#define DNS_ERROR_CHECKFILTZ  34
#define DNS_ERROR_CHECKSCALAR 35
#define DNS_ERROR_ZERODELTAU  36
#define DNS_ERROR_CGORDER     37
#define DNS_ERROR_NEGDENS     38
#define DNS_ERROR_NEGPRESS    39
#define DNS_ERROR_SMALLTSTEP  40
#define DNS_ERROR_MAXPROC     41
#define DNS_ERROR_MINPROC     42
#define DNS_ERROR_RECLEN      43
#define DNS_ERROR_PARPARTITION 45
#define DNS_ERROR_MPITYPECHECK 46
#define DNS_ERROR_INFLOWDOMAIN 47
#define DNS_ERROR_DIMGRID      48
#define DNS_ERROR_MPIBUFFER    49
#define DNS_ERROR_MPICOUNT     50
#define DNS_ERROR_MPITYPE      51
#define DNS_ERROR_MPIRANK      52
#define DNS_ERROR_MPITAG       53
#define DNS_ERROR_MPICOMM      54
#define DNS_ERROR_IMAXP        55
#define DNS_ERROR_JMAXP        56
#define DNS_ERROR_KMAXP        57
#define DNS_ERROR_INFIMAXP     58
#define DNS_ERROR_INFJMAXP     59
#define DNS_ERROR_INFKMAXP     60
#define DNS_ERROR_INVALOPT     61
#define DNS_ERROR_POLARRMAX    62
#define DNS_ERROR_POLARTMAX    63
#define DNS_ERROR_POLAREVEN    64
#define DNS_ERROR_MAXRADIUS    65
#define DNS_ERROR_REGRID       66
#define DNS_ERROR_BISPEV       67
#define DNS_ERROR_VORTSUBDOM   68
#define DNS_ERROR_MINWRKAVG    69
#define DNS_ERROR_CURFIT       70
#define DNS_ERROR_LESEVEN      71
#define DNS_ERROR_LESN         72
#define DNS_ERROR_MAXVORT      73
#define DNS_ERROR_VTYPEUNK     74
#define DNS_ERROR_FTYPEUNK     75
#define DNS_ERROR_SUBVORT      76
#define DNS_ERROR_CALCFLOW     77
#define DNS_ERROR_POISSON      78
#define DNS_ERROR_TMPSIZE      79
#define DNS_ERROR_ALLOC        80
#define DNS_ERROR_KVSICBC      81
#define DNS_ERROR_SPPR         82
#define DNS_ERROR_WRKOVERFLW   83
#define DNS_ERROR_FLTDIM       84
#define DNS_ERROR_OPTION       85
#define DNS_ERROR_PARAMETER    86
#define DNS_ERROR_WRKSIZE      87
#define DNS_ERROR_NOTIMPL      88
#define DNS_ERROR_INFFLTDOM    89
#define DNS_ERROR_WRONGSCALNB  90
#define DNS_ERROR_THERMOFORMAT 91
#define DNS_ERROR_THERMOCONT   92
#define DNS_ERROR_MECHFORMAT   93
#define DNS_ERROR_TRSPFORMAT   94
#define DNS_ERROR_SCALFORMAT   95
#define DNS_ERROR_TOTALVARS    96
#define DNS_ERROR_UNITS        97
#define DNS_ERROR_NBREAC       98
#define DNS_ERROR_IJMISMATCH   99
#define DNS_ERROR_FLTWIDTH     100
#define DNS_ERROR_NORMALIZE    101
#define DNS_ERROR_STATZERO     102
#define DNS_ERROR_GRID_SCALE   103
#define DNS_ERROR_UNDEVELOP    104
#define DNS_ERROR_READMREAL    105
#define DNS_ERROR_READMINT     106
#define DNS_ERROR_FAVFLUCT     107
#define DNS_ERROR_STFILE       108
#define DNS_ERROR_CRAY         109
#define DNS_ERROR_PSFFT        110
#define DNS_ERROR_DILATATION   111
#define DNS_ERROR_AVGTMP       112
#define DNS_ERROR_THOMAS       113
#define DNS_ERROR_BUFFER       114

#define DNS_ERROR_PARTICLE          200
#define DNS_ERROR_PARTICLE_UNMATCH  201
#define DNS_ERROR_CALCPARTICLE      202
#define DNS_ERROR_CALCTRAJECTORIES  203
#define DNS_ERROR_RESIDENCERESET    204

#define DNS_ERROR_CUBIC_SPLINE      300

#define DNS_ERROR_IBM_SPLINE        400
#define DNS_ERROR_IBM_GEOMETRY      401
#define DNS_ERROR_IBM_MISS_GEO      402
#define DNS_ERROR_IBM_GAMMA         403
#define DNS_ERROR_IBM_SHEAR         404
#define DNS_ERROR_IBM_INITIALIZE    405

#define DNS_ERROR_PRESSURE_DECOMPOSITION    500

#define DNS_ERROR_AVG_PHASE         600

#endif
