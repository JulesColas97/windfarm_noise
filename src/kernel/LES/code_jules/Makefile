FC = h5pfc -fpp
#FC += -O0 -g -traceback -check bounds -debug all -warn all
#FC += -O3 -shared-intel -heap-arrays -mcmodel=medium # -xAVX -axCORE-AVX2
FC += -O3# -mkl #-fma -ip -finline-functions -align array64byte -xAVX -axCORE-AVX2
#FC+= -save-temps #Generates version of the actually compiled code, removing pre-processor flags (*.i90 files)

# Use binary output for output files
windbreaks = no
turbines = no
scalar   = no
coriolis = no
baroclinic = no
lvlset   = no
# Data for higher order statistics
statske= yes
# Stats for txx,txy,txz,tyy,tyz,tzz
statstau = yes
# Check tavg_pdf before use (mpi problems when restarting)
tavg_pdf = no
# tavg_spectrum has to be used with care (see spec.f90)
tavg_spectrum = no
#tavg_sp dumps out time averaged files in single precision 
tavg_sp = no

# Concurrent precursor simulation
# Requires configuration of two builds: upstream and downstream
cps = yes
cps_y = no

# Actuator line/section model
atm = no
# Actuator disk method with rotation
# if admr=yes, there must be atm=yes
admr = no

# Recommended settings for production simulations
# USE_DEBUG=no # Yes: FFTW library in debug mode
# USE_MEMORY=yes #No: Use less memory, which may be faster: yes: use less memory, which may be slower
debug  = no
memory = yes

ifeq ($(statske), yes)
FC += -DSTATSKE
endif

ifeq ($(statstau), yes)
FC += -DSTATSTAU
endif

ifeq ($(tavg_pdf), yes)
FC += -DTAVG_PDF
endif

ifeq ($(tavg_spectrum), yes)
FC += -DTAVG_SPECTRUM
endif

ifeq ($(tavg_sp), yes)
FC += -DTAVG_SP
endif

ifeq ($(scalar), yes)
FC += -DSCALAR
endif

ifeq ($(debug), yes)
FC += -DDEBUG
endif

ifeq ($(coriolis), yes)
FC += -DCORIOLIS
endif

ifeq ($(baroclinic), yes)
FC += -DBAROCLINIC
endif

ifeq ($(windbreaks), yes)
FC += -DWINDBREAKS
endif

ifeq ($(turbines), yes)
FC += -DTURBINES
endif

ifeq ($(admr), yes)
FC += -DADMR
endif

ifeq ($(atm), yes)
FC += -DATM
endif

ifeq ($(lvlset), yes)
FC += -DLVLSET
endif

ifeq ($(cps), yes)
FC += -DCPS
endif

ifeq ($(cps_y), yes)
FC += -DCPS_Y
endif

ifeq ($(memory), yes)
FC += -DMEMORY
endif

#======================================================================
# Library
#======================================================================
LDFLAGS = $(FFTW3_LDFLAGS)
FC += $(FFTW3_FFLAGS) -module $(OBJDIR)
#The above option would not work on a gfortran compiler
#FC += -fintrinsic-modules-path $(OBJDIR)

#====================================================================== 
# Make program
#======================================================================
PROGRAM = LESsolver

ifeq ($(lvlset), yes)
PROGRAM = LESsolver_lvst
endif
ifeq ($(cps), yes)
PROGRAM = LESsolver_cps
endif
ifeq ($(turbines), yes)
PROGRAM = LESsolver_turbines
endif
ifeq ($(turbines), yes)
ifeq ($(lvlset), yes)
PROGRAM = LESsolver_turbines_lvst
endif
endif
ifeq ($(turbines), yes)
ifeq ($(cps),yes)
PROGRAM = LESsolver_cps_turbines
endif
endif
ifeq ($(lvlset), yes)
ifeq ($(cps),yes)
PROGRAM = LESsolver_cps_lvst
endif
endif
ifeq ($(turbines), yes)
ifeq ($(lvlset), yes)
ifeq ($(cps),yes)
PROGRAM = LESsolver_cps_turbines_lvst
endif
endif
endif

ifeq ($(atm), yes)
PROGRAM = LESsolver_atm
ifeq ($(admr), yes)
PROGRAM = LESsolver_admr
endif
ifeq ($(cps),yes)
PROGRAM = LESsolver_cps_atm
ifeq ($(admr), yes)
PROGRAM = LESsolver_cps_admr
endif
ifeq ($(lvlset), yes)
PROGRAM = LESsolver_cps_atm_lvst
endif
endif
endif

ifeq ($(windbreaks), yes)
PROGRAM = LESsolver_windbreaks
endif
ifeq ($(turbines), yes)
ifeq ($(windbreaks), yes)
PROGRAM = LESsolver_turbines_windbreaks
endif
endif
ifeq ($(cps), yes)
ifeq ($(windbreaks), yes)
PROGRAM = LESsolver_cps_windbreaks
endif
endif
ifeq ($(cps), yes)
ifeq ($(turbines), yes)
ifeq ($(windbreaks), yes)
PROGRAM = LESsolver_cps_turbines_windbreaks
endif
endif
endif

#Fortran source files (excluding the module files) 
FFILES = amd_model.f90 cfl.f90 convec.f90 derivatives.f90 divstress_uv.f90 divstress_w.f90 \
         emul_complex.f90 fft.f90 filters.f90 forcing.f90 fringe.f90 grid.f90 \
         ic.f90 initial.f90 initialize.f90 input_util.f90 interpolag_Sdep.f90 \
         io.f90 lagrange_Sdep.f90 main.f90 monitor.f90 mpi_defs.f90 padd.f90 \
         pdf.f90 press_stag.f90 ran3.f90 sgs_stag.f90 smag_model_init.f90 spec.f90 \
	 tridag.f90 trilinear_interp.f90 unpadd.f90 wallstress.f90 

#Files for reading and writing
FFILES +=  read_2Dplane.f90 read_3Dfield.f90 write_2Dplane.f90 write_3Dfield.f90 write_3Dfield_SinglePrecision.f90 

#Object and module directory 

#Module files 
MFILES = param.f90 sim_param.f90 sgs_param.f90

#Files for printing only 1D time averaged data
ifneq ($(scalar),yes)
FFILES += read_1Dfield.f90 write_1Dfield.f90 
endif

ifeq ($(windbreaks), yes)
MFILES += windbreak_mod.f90
FFILES += windbreak.f90
endif

ifeq ($(turbines), yes)
MFILES += turbine_mod.f90 turbine_indicator.f90
FFILES += turbine.f90
endif

ifeq ($(atm), yes)
MFILES += atm_base.f90 atm_input_util.f90 actuator_turbine_model.f90 atm_lesgo_interface.f90
endif

ifeq ($(cps),yes)
FFILES += concurrent_precursor.f90
endif

ifeq ($(coriolis),yes)
FFILES += wind_angle_controller.f90
endif

ifeq ($(scalar),yes)
MFILES += scalars_param.f90 
FFILES += damping_layer.f90 ic_scal.f90 interpolag_Sdep_Scalar.f90 \
          scalars_module.f90 wallstress_scalar.f90
FFILES += read_1Dfield.f90 write_1Dfield.f90
endif

ifeq ($(tavg_pdf),yes)
FFILES += read_1Dfield.f90 read_PDF.f90 write_1Dfield.f90 write_PDF.f90
endif

ifeq ($(tavg_spectrum), yes)
FFILES += write_streamw_spec.f90 write_spanw_spec.f90 read_streamw_spec.f90 read_spanw_spec.f90
endif

ifeq ($(lvlset),yes)
FFILES += level_set.f90
endif

OBJDIR =obj
OBJS  := $(FFILES:%.f90=$(OBJDIR)/%.o)
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

#=======================================================================
# Linking
#=======================================================================
all: directories $(PROGRAM)

$(PROGRAM): $(MOBJS) $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)  

#=======================================================================
# Dependencies
#=======================================================================

$(OBJDIR)/param.o: param.f90
	$(FC) -c -o $@ $< $(LDFLAGS)  
$(OBJDIR)/sim_param.o: sim_param.f90 $(OBJDIR)/param.o 
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/sgs_param.o: sgs_param.f90
	$(FC) -c -o $@ $< $(LDFLAGS)

ifeq ($(scalar),yes) 
$(OBJDIR)/scalars_param.o: scalars_param.f90 
	$(FC) -c -o $@ $< $(LDFLAGS)
endif 

ifeq ($(turbines), yes)
$(OBJDIR)/turbine_indicator.o: turbine_indicator.f90 
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/turbine_mod.o: turbine_mod.f90 $(OBJDIR)/turbine_indicator.o
	$(FC) -c -o $@ $< $(LDFLAGS)
endif

ifeq ($(atm), yes)
$(OBJDIR)/atm_base.o: atm_base.f90 
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/atm_input_util.o: atm_input_util.f90
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/actuator_turbine_model.o: actuator_turbine_model.f90 $(OBJDIR)/atm_input_util.o
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/atm_lesgo_interface.o: atm_lesgo_interface.f90 $(OBJDIR)/actuator_turbine_model.o
	$(FC) -c -o $@ $< $(LDFLAGS)
endif

ifeq ($(windbreaks),yes)
$(OBJDIR)/windbreak_mod.o: windbreak_mod.f90 
	$(FC) -c -o $@ $< $(LDFLAGS)
endif

$(OBJDIR)/%.o:   %.f90 $(MOBJS)
	$(FC) -c -o $@ $< $(LDFLAGS)


#=======================================================================
# Clean up
#=======================================================================
clean: 
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*genmod* $(OBJDIR)/*.o a.out *genmod* *.mod *.o obj\

.PHONY: directories
directories: $(OBJDIR)
$(OBJDIR): 
	mkdir -p ${OBJDIR}
