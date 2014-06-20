require 'narray'

L =[-0.12940952255092145,
    0.22414386804185735,
    0.836516303737469,
    0.48296291314469025]


#L= [0.0014009155259146807,
#    0.0006197808889855868,
#    -0.013271967781817119,
#    -0.01152821020767923,
#    0.03022487885827568,
#    0.0005834627461258068,
#    -0.05456895843083407,
#    0.238760914607303,
#    0.717897082764412,
#    0.6173384491409358,
#    0.035272488035271894,
#    -0.19155083129728512,
#    -0.018233770779395985,
#    0.06207778930288603,
#    0.008859267493400484,
#    -0.010264064027633142,
#    -0.0004731544986800831,
#    0.0010694900329086053 ]

#L = [-0.0033824159510050025955,
#     -0.00054213233180001068935,
#      0.031695087811525991431,
#      0.0076074873249766081919,
#     -0.14329423835127266284,
#     -0.061273359067811077843,
#      0.48135965125905339159,
#      0.77718575169962802862,
#      0.36444189483617893676,
#     -0.051945838107881800736,
#     -0.027219029917103486322,
#      0.049137179673730286787,
#      0.0038087520138944894631,
#     -0.014952258337062199118,
#     -0.00030292051472413308126,
#      0.0018899503327676891843]

h = []
L.each_with_index { |e,i|
  if i % 2 == 0 then
    h.push(e)
  else
    h.push(-e)
  end
}

H = h.reverse

LRE = L.reverse.values_at(*(0..(L.length-1)).step(2).collect)
LRO = L.reverse.values_at(*(1..(L.length-1)).step(2).collect)
HRE = H.reverse.values_at(*(0..(L.length-1)).step(2).collect)
HRO = H.reverse.values_at(*(1..(L.length-1)).step(2).collect)

rng=Random::new(1234)
a=[]
34.times { |i|
  a[i] = 0 #rng.rand
}

a[0]=1
a[-1]=-1
#32.times { |i|
#  if i == 0 then
#    a[i] = 1.0
#  elsif i == 31 then
#    a[i] = 2.0
#  else
#    a[i] = 0.0
#  end
#}
def dwt(data, center)
  l = []
  h = []

  (data.length/2).times { |i|
    li = 0.0
    hi = 0.0
    L.length.times { |indx|
      e_indx = (indx - center + 1 + 2*i) % data.length
      li += L[indx]*data[e_indx]
      hi += H[indx]*data[e_indx]
    }
    l[i] = li
    h[i] = hi
  }
  return [l, h]
end

def idwt(l,h,c)
  d = []
  center = (L.length - c)/2
  l.length.times { |i|
    h1ei = 0.0
    h1oi = 0.0
    h2ei = 0.0
    h2oi = 0.0
    
    LRE.length.times { |indx|
      e_indx = (indx - center + i) % l.length
     
      h1ei += l[e_indx]*LRE[indx]
      h1oi += l[e_indx]*LRO[indx]
      h2ei += h[e_indx]*HRE[indx]
      h2oi += h[e_indx]*HRO[indx]
    }
    if c % 2 == 0 then
      ind = (2*i-1) % (l.length*2)
      d[2*i] = h1ei + h2ei
      d[ind] = h1oi + h2oi
    else
      d[2*i] = h1oi + h2oi
      d[2*i+1] = h1ei + h2ei
    end
  }
  return d
end


EPSILON = 1e-12
d = nil
H.length.times { |j|
puts "center: #{j}"
  l, h = dwt(a, j)
  d = idwt(l, h, j)
  puts l.inspect
  puts h.inspect
  d.each_index { |i|
    raise "Error!#{d[i]-a[i]}" if (d[i]-a[i]).abs > EPSILON
  }
}


(8..8).each{ | center|
puts "test#{center}"
sum_a=0
a.each{ |e| sum_a+=e*e}
l, h = dwt(a, center)
sum_l=0
l.each{ |e| sum_l+=e*e}
sum_h=0
h.each{ |e| sum_h+=e*e}
#puts 'sum_a,sum_l,sum_h,sum_l+h',sum_a,sum_l,sum_h,sum_l+sum_h
#puts 'diff',sum_a - (sum_l+sum_h)
#puts l.inspect
#puts h.inspect
d = idwt(l, h, center)
}


l=[]
h=[]
17.times { |i|
  l[i] = 0#rng.rand
  h[i] = 0 #rng.rand
}
l[0]=1
l[-1]=-1
h[0]=1
h[-1]=-1


H.length.times { |j|
  puts "center: #{j}"
  d = idwt(l, h, j)
  puts d.inspect
  l2, h2 = dwt(d, j)
  l.each_index { |i|
    raise "Error!" if (l[i]-l2[i]).abs > EPSILON
    raise "Error!" if (h[i]-h2[i]).abs > EPSILON
  }
}
