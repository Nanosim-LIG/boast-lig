require 'BOAST'
require './KBoast.rb'
require 'narray_ffi'

include BOAST

set_default_real_size(4)

#k_boast_params = {:kernel => :boast, :LDFLAGS => "-lgfortran -L/usr/lib/ -lblas", :FCFLAGS => "-fimplicit-none -fexternal-blas -O3"}  
k_boast_params = {:kernel => :boast, :LDFLAGS => "-lgfortran", :FCFLAGS => "-fimplicit-none -O3"}

set_lang(FORTRAN)
stats = []
repeat = 5

# Creating boast kernel
k = KBoast::new(k_boast_params)
k.generate
puts k.kernel
k.kernel.build(:LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS] )

inputs = k.kernel.load_ref_inputs()
outputs = k.kernel.load_ref_outputs()


l=[]

inputs.each_key { |key|

  k.kernel.run(*(inputs[key]))
  repeat.times {
    stats.push k.kernel.run(*(inputs[key]))
  }

  l[1]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_L1_DCM"])
  l[2]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_L2_DCM"])
  l[3]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_TOT_CYC"])
  l[4]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_TOT_INS"])
  l[5]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_BR_MSP"])
  l[6]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_L1_LDM"])
  l[7]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_L1_STM"])
  l[8]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_VEC_SP"])
  l[9]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_BR_INS"])
  l[10]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_L2_TCA"])
  l[11]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_LD_INS"])
  l[12]=k.kernel.run(*(inputs[key]), :PAPI =>["PAPI_SR_INS"])

  #puts k.kernel.compare_ref(outputs[key], inputs[key]).inspect
}


stats.sort_by! { |a| a[:duration] }
p stats
stats = stats.first

puts " "
puts "-------------------------- PAPI -------------------------"
puts "-> Data cache misses: L1=#{l[1][:PAPI]["PAPI_L1_DCM"]}  L2=#{l[2][:PAPI]["PAPI_L2_DCM"]}"
puts "-> Total cache accesses: L2=#{l[10][:PAPI]["PAPI_L2_TCA"]}"
puts "-> Level 1 load misses = #{l[6][:PAPI]["PAPI_L1_LDM"]}"
puts "-> Level 1 store misses = #{l[7][:PAPI]["PAPI_L1_STM"]}"
puts "-> Load instructions = #{l[11][:PAPI]["PAPI_LD_INS"]}"
puts "-> Store instructions = #{l[12][:PAPI]["PAPI_SR_INS"]}"
puts "-> Total cycles = #{l[3][:PAPI]["PAPI_TOT_CYC"]}"
puts "-> Instructions completed = #{l[4][:PAPI]["PAPI_TOT_INS"]}"
puts "=> Instruction/cycle = #{l[4][:PAPI]["PAPI_TOT_INS"]/l[3][:PAPI]["PAPI_TOT_CYC"]}"
puts "-> Conditional branch instructions mispredicted = #{l[5][:PAPI]["PAPI_BR_MSP"]}"
puts "-> Single precision vector/SIMD instructions = #{l[8][:PAPI]["PAPI_VEC_SP"]}"
puts "-> Branch instructions = #{l[9][:PAPI]["PAPI_BR_INS"]}"
puts " "
puts "Floating point operations = 68914896"
puts "#{k.kernel.procedure.name}: #{stats[:duration]*1.0e3} ms    #{ 68914896 / (stats[:duration]*1.0e9)} GFlops"
puts " "
