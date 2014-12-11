################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BasisClass.cpp \
../src/Ident_hist.cpp \
../src/ParseSNP.cpp \
../src/VShape.cpp 

OBJS += \
./src/BasisClass.o \
./src/Ident_hist.o \
./src/ParseSNP.o \
./src/VShape.o 

CPP_DEPS += \
./src/BasisClass.d \
./src/Ident_hist.d \
./src/ParseSNP.d \
./src/VShape.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/software/ngm/src/bamtools/include -I/software/ngm/lib/tclap/ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


