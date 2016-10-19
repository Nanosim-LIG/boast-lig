require './WaveletFilters.rb'
require 'bigdecimal'

sym8_lp = SYM8_LP.collect { |e| BigDecimal.new(e) }
sym8_hp = sym8_lp.each_with_index.collect { |e,i| i % 2 == 0 ? e : -e }.reverse
sym8_mf = SYM8_MF.collect { |e| BigDecimal.new(e) }

sym8_lp_center = sym8_lp.length / 2
sym8_hp_center = sym8_hp.length / 2
sym8_mf_center = 7
sym8_mf_lr = [ -sym8_mf_center, sym8_mf.length - sym8_mf_center - 1 ]
sym8_lp_lr = [ -sym8_lp_center, sym8_lp.length - sym8_lp_center - 1 ]
sym8_hp_lr = [ -sym8_hp_center, sym8_hp.length - sym8_hp_center - 1 ]
sym8_mf_lp_lr = [ sym8_mf_lr[0] - sym8_lp_lr[1], sym8_mf_lr[1] - sym8_lp_lr[0] ]
sym8_mf_hp_lr = [ sym8_mf_lr[0] - sym8_hp_lr[1], sym8_mf_lr[1] - sym8_hp_lr[0] ]

sym8_mf_lp = (sym8_mf_lp_lr[0]..sym8_mf_lp_lr[1]).collect { |i|
  sum = 0
  (sym8_lp_lr[0]..sym8_lp_lr[1]).each { |j|
    sum += sym8_lp[j+sym8_lp_center] * sym8_mf[i-j-1+sym8_mf_center] if (0...sym8_mf.length).include?(i-j-1+sym8_mf_center)
  }
  sum
}
sym8_mf_lp.each { | e| puts e }

puts "------------------------------------------------"

sym8_mf_hp = (sym8_mf_hp_lr[0]..sym8_mf_hp_lr[1]).collect { |i|
  sum = 0
  (sym8_hp_lr[0]..sym8_hp_lr[1]).each { |j|
    sum += sym8_hp[j+sym8_hp_center] * sym8_mf[i-j-1+sym8_mf_center] if (0...sym8_mf.length).include?(i-j-1+sym8_mf_center)
  }
  sum
}
sym8_mf_hp.each { | e| puts e }

