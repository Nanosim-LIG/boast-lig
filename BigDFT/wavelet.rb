require 'narray'

L = [-0.0033824159510050025955,
     -0.00054213233180001068935,
      0.031695087811525991431,
      0.0076074873249766081919,
     -0.14329423835127266284,
     -0.061273359067811077843,
      0.48135965125905339159,
      0.77718575169962802862,
      0.36444189483617893676,
     -0.051945838107881800736,
     -0.027219029917103486322,
      0.049137179673730286787,
      0.0038087520138944894631,
     -0.014952258337062199118,
     -0.00030292051472413308126,
      0.0018899503327676891843]
h = []
L.each_with_index { |e,i|
  if i % 2 == 0 then
    h.push(e)
  else
    h.push(-e)
  end
}

H = h.reverse
#H = [-0.0018899503327594609,
#     -0.0003029205147213668,
#      0.01495225833704823,
#      0.003808752013890615,
#     -0.049137179673607506,
#     -0.027219029917056003,
#      0.05194583810770904,
#      0.3644418948353314,
#     -0.7771857517005235,
#      0.4813596512583722,
#      0.061273359067658524,
#     -0.1432942383508097,
#     -0.007607487324917605,
#      0.03169508781149298,
#      0.0005421323317911481,
#     -0.0033824159510061256]

R = L.reverse
RE = R.values_at(*(0..(R.length-1)).step(2).collect)
RU = R.values_at(*(1..(R.length-1)).step(2).collect)
SE = L.values_at(*(0..(R.length-1)).step(2).collect)
SU = L.values_at(*(1..(R.length-1)).step(2).collect).collect { |e| -e }

LE = L.values_at(*(0..(L.length-1)).step(2).collect).reverse
LO = L.values_at(*(1..(L.length-1)).step(2).collect).reverse
HE = H.values_at(*(0..(H.length-1)).step(2).collect).reverse
HO = H.values_at(*(1..(H.length-1)).step(2).collect).reverse

a=[]
32.times { |i|
  if i == 0 then
    a[i] = 1
  elsif i == 31 then
    a[i] = 2
  else
    a[i] = 0
  end
}
puts a.inspect
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
  l.length.times { |i|
    h1ei = 0.0
    h1oi = 0.0
    h2ei = 0.0
    h2oi = 0.0
    center = c/2
    center -= c - (c/2)*2
    
    RE.length.times { |indx|
      e_indx = (indx - center + i) % l.length
#      puts RE[indx]
#      puts RU[indx]
#      puts SE[indx]
#      puts SU[indx]
#      puts
     
      h1ei += l[e_indx]*RE[indx]
      h1oi += l[e_indx]*RU[indx]
      h2ei += h[e_indx]*SE[indx]
      h2oi += h[e_indx]*SU[indx]
    }
    if center % 2 == 0 then
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
#puts a.inspect
l, h = dwt(a,8)
puts l.inspect
puts h.inspect
puts idwt(l, h, 8).inspect

l, h = dwt(a,9)
puts l.inspect
puts h.inspect
puts idwt(l, h, 9).inspect
