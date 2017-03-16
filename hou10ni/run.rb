require 'BOAST'
include BOAST
require './KOssRef.rb'
require 'narray_ffi'

class Experiment

 def self.run()
  k_ref_params = {:kernel => :ref, :preprocessor => false, :LDFLAGS => "-lgfortran"}  

  set_lang(FORTRAN)
  kernels={}

  # Creating ref kernel
  kernels[k_ref_params] = KOssRef::new(k_ref_params)
  kernels[k_ref_params].generate
  puts kernels[k_ref_params].kernel
  kernels[k_ref_params].kernel.build(:LDFLAGS => k_ref_params[:LDFLAGS])
 
  # Setting the values of arguments
  seed = 10
  ANArray.srand(seed) unless seed.nil?
  ## Ces valeurs different selon le processus MPI
  ## P0 AfluSize=16952560  Size=879732  Nflu_inner=210025  idx_vec_flu(12,214811) idx_mat_flu(214812)
  ## P1 AfluSize=16955792  Size=881612  Nflu_inner=209629  idx_vec_flu(12,214805) idx_mat_flu(214806)
  ## P2 AfluSize=16950864  Size=880460  Nflu_inner=209830  idx_vec_flu(12,214820) idx_mat_flu(214821)
  ## P3 AfluSize=16953312  Size=880240  Nflu_inner=209986  idx_vec_flu(12,214811) idx_mat_flu(214812)
  #@@Nflu_inner=210025
  @@Nflu_inner=2
  @@Nflusol_inner=0
  @@nb_rhs=1
  @@afluSize=16952560
  @@pSize=879732
  @@vecSize_2=214811
  @@vecSize_1=12
  @@matSize=214812
  @@idx_vec_flu=ANArray.int(64,@@vecSize_1,@@vecSize_2).random!
  @@idx_mat_flu=ANArray.int(64,@@matSize).random!
  @@P_new=ANArray.float(64,@@pSize,@@nb_rhs).random!
  @@P_old=ANArray.float(64,@@pSize,@@nb_rhs).random!
  @@P_inter=ANArray.float(64,@@pSize,@@nb_rhs).random!
  @@A_flu=ANArray.float(64,@@afluSize).random!

  puts kernels[k_ref_params].kernel.run(@@Nflu_inner, @@Nflusol_inner, @@nb_rhs, @@idx_vec_flu, @@idx_mat_flu, @@P_new, @@P_inter, @@P_old, @@A_flu, @@afluSize, @@pSize,  @@vecSize_1, @@vecSize_2, @@matSize)

  return kernels
 end
end

Experiment.run
