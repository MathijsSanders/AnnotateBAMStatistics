#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Algorithm.o \
	${OBJECTDIR}/AnnovarFile.o \
	${OBJECTDIR}/BamFile.o \
	${OBJECTDIR}/Mutex.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -pthread
CXXFLAGS=-m64 -pthread

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L./samtools -lbam -lz

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/AnnotateBAMStatistics

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/AnnotateBAMStatistics: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/AnnotateBAMStatistics ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/Algorithm.o: Algorithm.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -O -Werror -I./samtools -std=c++0x -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Algorithm.o Algorithm.cpp

${OBJECTDIR}/AnnovarFile.o: AnnovarFile.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -O -Werror -I./samtools -std=c++0x -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/AnnovarFile.o AnnovarFile.cpp

${OBJECTDIR}/BamFile.o: BamFile.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -O -Werror -I./samtools -std=c++0x -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BamFile.o BamFile.cpp

${OBJECTDIR}/Mutex.o: Mutex.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -O -Werror -I./samtools -std=c++0x -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Mutex.o Mutex.cpp

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -O -Werror -I./samtools -std=c++0x -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
