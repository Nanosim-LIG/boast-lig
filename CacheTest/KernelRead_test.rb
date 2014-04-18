require'./KernelRead.rb'
require'./KernelRead2.rb'
#puts k.print

m_start = 0
m_cycle = 1024*5
m_stride = 1
page_size = 4096

output = NArray.int(1024*12).random(1000)

k = BOAST::kernel_read_vectorized2(8)
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, page_size*page, output)
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{page_size*page*m_cycle/(stats[:duration]*1.0e9)}"
}

#k = BOAST::kernel_read_vectorized(4,8,4)
#k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
#(1..10).each { |page|
#  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page, output)
#  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{4*buffer_size*page*m_cycle/(stats[:duration]*1.0e9)}"
#}
k = BOAST::kernel_read_vectorized2(8,4,4)
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, page_size*page, output)
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{page_size*page*m_cycle/(stats[:duration]*1.0e9)}"
}

#k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
#(1..10).each { |page|
#  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page, output)
#  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{4*buffer_size*page*m_cycle/(stats[:duration]*1.0e9)}"
#}
#
#k = BOAST::kernel_read_ref(8)
#puts k.print
#k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
#(1..10).each { |page|
#  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page, output)
#  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{4*buffer_size*page*m_cycle/(stats[:duration]*1.0e9)}"
#}
#
#k = BOAST::kernel_read_ref(4, 8)
#puts k.print
#k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
#(1..10).each { |page|
#  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
#  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{8*buffer_size*page*m_cycle/(2*stats[:duration]*1.0e9)}"
#}
#
#k = BOAST::kernel_read_ref(6, 8)
#puts k.print
#k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
#(1..10).each { |page|
#  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
#  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{8*buffer_size*page*m_cycle/(2*stats[:duration]*1.0e9)}"
#}
#
#k = BOAST::kernel_read_ref(8, 8)
#puts k.print
#k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
#(1..10).each { |page|
#  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
#  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{8*buffer_size*page*m_cycle/(2*stats[:duration]*1.0e9)}"
#}
#
#k = BOAST::kernel_read_ref(16, 8)
#puts k.print
#k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
#(1..10).each { |page|
#  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/2, output)
#  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{8*buffer_size*page*m_cycle/(2*stats[:duration]*1.0e9)}"
#}
