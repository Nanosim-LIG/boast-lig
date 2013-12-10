require'./KernelRead.rb'
#puts k.print

machine="sse3"

vect=nil
flag=nil
if machine == "intel" then
vect="sse3"
flag="-msse3"
elsif machine == "arm" then
vect="neon"
flag="-mfpu=neon"
elsif machine == "avx" then
vect="avx"
flag="-mavx"
end

m_start = 0
m_cycle = 1024*8
m_stride = 1
buffer_size = 1024*4

element = 4

NArray.srand(42)
output = NArray.int(1024*12).random(1000)

(1..2).each { |m|
element *= m
(1..10).each { |loop|
# WRITING CODE AND COMPILING:
puts "* Unrolled (#{loop}) times, with element size (#{element*8}b):"
k = BOAST::kernel_read_ref(loop, element)
puts "** Code:"
puts k.print
k.build({:CC => 'gcc',:FCFLAGS => "-O3"})
# WARMUP:
(1..10).each { |page|
  k.run(m_start, m_cycle, m_stride, buffer_size*page/element, output)
}
# EXECUTION:
puts "** Results:"
puts "NAME LOOP ELEMENT PAGE TIME BANDWIDTH ID" 
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, buffer_size*page/element, output)
      puts "#{k.procedure.name}: #{loop} #{element*8} #{page} #{stats[:duration]*1.0e3} #{buffer_size*page*m_cycle/(stats[:duration]*1.0e9)} #{stats[:return]}"
}
puts "#################"
}
}

# VECTORIZED INSTRUCTIONS
element = 8
length = 2
if machine == "avx" then
length = 4
end

(1..10).each { |loop|
# WRITING CODE AND COMPILING:
puts "* Vectorized (#{vect}), unrolled (#{loop}) times, with element size (#{element*length*8}b):"
k = BOAST::kernel_read_vectorized(loop, element, length, vect)
puts "** Code:"
puts k.print
    k.build({:CC => 'gcc',:CFLAGS => "-O3 #{flag}"})
# WARMUP:
(1..10).each { |page|
  k.run(m_start, m_cycle, m_stride, (buffer_size*page) / (element*length), output)
}
# EXECUTION:
puts "** Results:"
puts "NAME LOOP ELEMENT PAGE TIME BANDWIDTH ID" 
(1..10).each { |page|
  stats = k.run(m_start, m_cycle, m_stride, (buffer_size*page) / (element*length), output)
      puts "#{k.procedure.name}: #{loop} #{element*length*8} #{page} #{stats[:duration]*1.0e3} #{buffer_size*page*m_cycle/(stats[:duration]*1.0e9)} #{stats[:return]}"
}
puts "#################"
}

