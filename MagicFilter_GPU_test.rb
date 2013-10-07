require './MagicFilter_GPU.rb'
require './MagicFilter.rb'

FILTER = [ "8.4334247333529341094733325815816e-7",
       "-0.1290557201342060969516786758559028e-4",
       "0.8762984476210559564689161894116397e-4",
       "-0.30158038132690463167163703826169879e-3",
       "0.174723713672993903449447812749852942e-2",
       "-0.942047030201080385922711540948195075e-2",
       "0.2373821463724942397566389712597274535e-1",
       "0.612625895831207982195380597e-1",
       "0.9940415697834003993178616713",
       "-0.604895289196983516002834636e-1",
       "-0.2103025160930381434955489412839065067e-1",
       "0.1337263414854794752733423467013220997e-1",
       "-0.344128144493493857280881509686821861e-2",
       "0.49443227688689919192282259476750972e-3",
       "-0.5185986881173432922848639136911487e-4",
       "2.72734492911979659657715313017228e-6" ]

require 'opencl'
require 'statsample'
platform = OpenCL::Platform::get_platforms.pop
device = OpenCL::Device::get_devices(platform, OpenCL::Device::TYPE_ALL).pop

def bench (k, name,n1,n2,n3,input,output,ref=nil)
#  puts "#{n1} #{n2} #{n3} #{[rndup(n3,16), rndup(n1*n2,16),1]}"
  epsilon = 10e-15
  stats = []
  k.build(:CLFLAGS => "-cl-mad-enable")
  prop = k.procedure.properties
  wgs = prop[:reqd_work_group_size]
  if not ref then
    20.times { stats.push( k.run(n3, n1*n2, input, output, :global_work_size => [rndup(n3,wgs[0]), rndup(n1*n2,wgs[1]),wgs[2]], :local_work_size => wgs )) }
  else
    20.times { stats.push( k.run(n3, n1*n2, input, output, :global_work_size => [wgs[0], rndup(n1*n2,wgs[1]),wgs[2]], :local_work_size => wgs )) }
  end
  durations = stats.collect{ |elem| elem[:duration] }
  durations.sort!
  puts "#{n3} #{k.procedure.name}#{name}: #{durations.first*1.0e3} #{32*n1*n2*n3 / (durations.first*1.0e9)} GFlops"
  if ref then
    diff = (ref - output).abs
    diff.each { |elem|
      raise "Error benching #{n3} for #{k.procedure.name}#{name}" if elem > epsilon
    }
  end
end
def rndup( val, div)
  return (val%div) == 0 ? val : val + div - (val%div)
end
(16..64).each {|n3|

n1 = 62
n2 = 66
input = NArray.float(n1,n2,n3).random
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1,n2,n3)
tmp = NArray.float(n1*n2*n3)

#k = ConvolutionGenerator::magicfilter_per_ref
#stats = k.run(n1, n2*n3, input, output_ref)
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#stats = k.run(n2, n3*n1, output_ref, tmp)
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#stats = k.run(n3, n1*n2, tmp, output_ref)
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#

bench(ConvolutionGenerator::magicfilter_GPU_per_ref,"",n1,n2,n3,input,output_ref)
#stats = k.run(n3, n1*n2, input, output, :global_work_size => [rndup(n3,16), rndup(n1*n2,16),1], :local_work_size => [16,16,1] )
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#stats = k.run(n2, n3*n1, output, tmp, :global_work_size => [rndup(n2,16), rndup(n3*n1,16),1], :local_work_size => [16,16,1] )
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#stats = k.run(n1, n2*n3, tmp, output, :global_work_size => [rndup(n1,16), rndup(n2*n3,16),1], :local_work_size => [16,16,1] )
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#
#diff = (output_ref - output).abs
#diff.each { |elem|
#  puts "Warning: residue too big: #{elem}" if elem > epsilon
#}
#
#
#stats = []
#stats.push( k.run(n3, n1*n2, input, output_ref, :global_work_size => [rndup(n3,16), rndup(n1*n2,16),1], :local_work_size => [16,16,1] ))
#durations = stats.collect{ |elem| elem[:duration] }
#durations.sort!
#puts "#{k.procedure.name}: #{durations.first*1.0e3} #{32*n1*n2*n3 / (durations.first*1.0e9)} GFlops"

bench(ConvolutionGenerator::magicfilter_GPU(FILTER,8,n3,256, device.local_mem_size ),"-1",n1,n2,n3,input,output,output_ref)

#k = ConvolutionGenerator::magicfilter_GPU(FILTER,8,n3,256, device.local_mem_size )
#stats = []
#stats.push(k.run(n3, n1*n2, input, output, :global_work_size => [8, rndup(n1*n2,8),1], :local_work_size => [8,8,1] ))
#puts "#{k.procedure.name}-1: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#
#
#diff = (output_ref - output).abs
#diff.each { |elem|
#  puts "Warning: residue too big: #{elem}" if elem > epsilon
#}
#
output.fill!(0.0)

bench(ConvolutionGenerator::magicfilter_GPU_next(FILTER,8,n3,256, device.local_mem_size ),"-2",n1,n2,n3,input,output,output_ref)
output.fill!(0.0)
#puts k.print
#k.build(:CLFLAGS => "-cl-mad-enable")
#stats = []
#stats.push( k.run(n3, n1*n2, input, output, :global_work_size => [16, rndup(n1*n2,16),1], :local_work_size => [16,16,1] ))
#puts "#{k.procedure.name}-2: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#
#diff = (output_ref - output).abs
#diff.each { |elem|
##  puts "Warning: residue too big: #{elem}" if elem > epsilon
#  raise "Error benching #{n3}" if elem > epsilon
#}
bench(ConvolutionGenerator::magicfilter_GPU_next_next(FILTER,8,n3,256, device.local_mem_size ),"-3",n1,n2,n3,input,output,output_ref)
output.fill!(0.0)
$stdout.flush
#puts k.print
#k.build(:CLFLAGS => "-cl-mad-enable")
#stats = []
#stats.push( k.run(n3, n1*n2, input, output, :global_work_size => [16, rndup(n1*n2,16),1], :local_work_size => [16,16,1] ))
#puts "#{k.procedure.name}-3: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} G,n1,n2,n3,input,output,output_refFlops"

#diff = (output_ref - output).abs
#diff.each { |elem|
#  puts "Warning: residue too big: #{elem}" if elem > epsilon
#  raise "Error benching #{n3}" if elem > epsilon
#}
}
