arr = nil
File::open("arm_neon.txt", "r") { |f|
  text = f.read.squeeze("\n")
  arr = text.scan( /.*?ARMv.*?$/m ).collect{ |e| e.split("\n").reject{|x| x.strip == ""}.values_at(0,-1) }
}
#p arr
arr.each { |e|
  splt = e[0].split(" ")
  #puts e[0] if splt.length != 2 or not splt[1].match("_")
  puts "  #{splt[1].inspect} => #{e[1].gsub("ARMv7","neon").gsub("ARMv8\(AArch64\)","asimd").gsub("ARMv8","asimd").split(", ")}"
}
