require'./KernelRead.rb'
#puts k.print

m_start = 0
m_cycle = 1024*5
m_stride = 1
buffer_size = 1024

output = NArray.int(1024*12).random


k = ConvolutionGenerator::kernel_read_ref
puts k.print
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page, output)
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{4*buffer_size*page*m_cycle/(stats[:duration]*1.0e9)}"
}

k = ConvolutionGenerator::kernel_read_ref(8)
puts k.print
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page, output)
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{4*buffer_size*page*m_cycle/(stats[:duration]*1.0e9)}"
}

k = ConvolutionGenerator::kernel_read_ref(4, 8)
puts k.print
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{8*buffer_size*page*m_cycle/(2*stats[:duration]*1.0e9)}"
}

k = ConvolutionGenerator::kernel_read_ref(6, 8)
puts k.print
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{8*buffer_size*page*m_cycle/(2*stats[:duration]*1.0e9)}"
}

k = ConvolutionGenerator::kernel_read_ref(8, 8)
puts k.print
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{8*buffer_size*page*m_cycle/(2*stats[:duration]*1.0e9)}"
}

k = ConvolutionGenerator::kernel_read_ref(16, 8)
puts k.print
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{8*buffer_size*page*m_cycle/(2*stats[:duration]*1.0e9)}"
}
