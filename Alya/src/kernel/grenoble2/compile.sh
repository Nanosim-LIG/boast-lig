rm -rf *.o *.mod
gfortran -c -J./ -I./ -ffree-line-length-none -x f95-cpp-input -O3 ./def_parameters.f90
gfortran -c -J./ -I./ -ffree-line-length-none -x f95-cpp-input -O3 ./mod_nsi_element_assembly.f90
gfortran -c -J./ -I./ -ffree-line-length-none -x f95-cpp-input -O3 ./mod_nsi_assembly.f90
gfortran -c -J./ -I./ -ffree-line-length-none -x f95-cpp-input -O3 ./grenoble.f90
gfortran -o grenoble.x *.o 
