require 'yaml'
require 'active_support/core_ext/hash/conversions'
f = File::open('data-3.3.14.xml')
hash = Hash.from_xml(f.read)
hash["intrinsics_list"]["intrinsic"].each { |intrinsic|
  list = []
  if intrinsic["CPUID"] then
    cpuid_list = [intrinsic["CPUID"]].flatten
    if cpuid_list.length > 1 then
      list += [cpuid_list]
    elsif cpuid_list.length == 1 then
      list += cpuid_list.first.split("/")
    end
    intrinsic["CPUID"]=list
  end
  list = []
  intrinsic["parameter"] = [intrinsic["parameter"]].flatten if intrinsic["parameter"]
  intrinsic["type"] = [intrinsic["type"]].flatten if intrinsic["type"]
  intrinsic["category"] = [intrinsic["category"]].flatten if intrinsic["category"]
  intrinsic["perfdata"] = [intrinsic["perfdata"]].flatten if intrinsic["perfdata"]
  intrinsic["instruction"] = [intrinsic["instruction"]].flatten if intrinsic["instruction"]
}
puts YAML::dump hash
