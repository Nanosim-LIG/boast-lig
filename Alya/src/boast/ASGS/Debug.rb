require '../Common/CommonArgs.rb'
require './KASGSRef.rb'
require './KASGSBoast.rb'

class Debug
  include CommonArgs
  def self.run
    vector_length = 1
    seed = 1
    nests = (1..23).to_a
    nb_instances = 1

    pgaus = 6
    pnode = 6
    dimension = 2

    kernels = []

    CommonArgs.init(pgaus,pnode,dimension,vector_length,nb_instances,seed)
    k_orig_params = {:vector_length => vector_length, :preprocessor => false, :nests => nests}
    
    # set_lang(FORTRAN)
    # set_fortran_line_length(300)
    # kernels[0] = KASGSRef::new(k_orig_params)
    # kernels[0].generate
    # kernels[0].kernel.build(:FCFLAGS => "-ffree-line-length-none -x f95-cpp-input -O3")

    set_lang(FORTRAN)
    kernels[1] = KASGSBoast::new(k_orig_params)
    kernels[1].generate
    puts kernels[1].kernel
    kernels[1].kernel.build(:FCFLAGS => "-ffree-line-length-none -x f95-cpp-input -O3")
    
  end
end

Debug.run
