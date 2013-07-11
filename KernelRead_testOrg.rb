require'./KernelRead.rb'
#puts k.print

m_start = 0
m_cycle = 1024*4
m_stride = 1
buffer_size = 1024

element = 4

NArray.srand(42)
output = NArray.int(1024*12).random(1000)

(1..2).each { |m|
element *= m
(1..10).each { |loop|
# WRITING CODE AND COMPILING:
puts "* Unrolled (#{loop}) times, with element size (#{element*8}b):"
k = ConvolutionGenerator::kernel_read_ref(loop, element)
puts "** Code:"
puts k.print
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
# WARMUP:
(1..10).each { |page|
  k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
}
# EXECUTION:
puts "** Results:"
puts "NAME LOOP ELEMENT PAGE TIME BANDWIDTH ID" 
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
      puts "#{k.procedure.name}: #{loop} #{element*8} #{page} #{stats[:duration]*1.0e3} #{8*buffer_size*page*m_cycle/(2*stats[:duration]*1.0e9)} #{stats[:return]}"
}
puts "#################"
}
}
