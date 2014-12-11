################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/sambamio/Alignment.cpp \
../src/sambamio/BamParser.cpp \
../src/sambamio/FastaParser.cpp \
../src/sambamio/Parser.cpp \
../src/sambamio/SamParser.cpp 

OBJS += \
./src/sambamio/Alignment.o \
./src/sambamio/BamParser.o \
./src/sambamio/FastaParser.o \
./src/sambamio/Parser.o \
./src/sambamio/SamParser.o 

CPP_DEPS += \
./src/sambamio/Alignment.d \
./src/sambamio/BamParser.d \
./src/sambamio/FastaParser.d \
./src/sambamio/Parser.d \
./src/sambamio/SamParser.d 


# Each subdirectory must supply rules for building sources it contributes
src/sambamio/%.o: ../src/sambamio/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/software/ngm/src/bamtools/include -I/software/ngm/lib/tclap/ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


