
Simulation += Simulation_data.o sim_readProfile.o sim_scaleProfile.o compute_tcool.o compute_Lambda.o compute_index.o

Simulation_initBlock.o : Simulation_data.o \
			 Driver_interface.o Grid_interface.o Eos_interface.o
Simulation_init.o : Simulation_data.o

Cool += Cool_exact_integrate.o compute_tcool.o compute_Lambda.o compute_index.o Cool_data.o
Cool += Heat_exact_integrate.o
Cool += Cool_compute_eint.o
#Cool_init.o : Cool_data.o


Driver += Driver_s.o

