require 'yaml'

module LogInfo
  @doc = {}
  @path = ""
  def self.init(path)
    @path = path 
    @doc[:header] = {}
    @doc[:body] = {}
    @doc[:body][:environment] = {}
    @doc[:body][:environment][:hardware] = {}
    @doc[:body][:environment][:software] = {} 
    @doc[:body][:kernel_info] = {}
  end
  def self.get_info
    @doc[:header][:title] = "Experiment information"
    @doc[:header][:time] = Time.now.utc
    @doc[:header][:machine] = `hostname`

    @doc[:body][:environment][:hardware][:cpu] = `less /proc/cpuinfo`
    @doc[:body][:environment][:hardware][:cpu_frequency] = `cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq`
    begin
      @doc[:body][:environment][:hardware][:gpu] = `nvidia-smi -q`
    rescue Errno::ENOENT
      @doc[:body][:environment][:hardware][:gpu] = `lshw -numeric -C display`
    end
    @doc[:body][:environment][:software][:linux] = `cat /proc/version` if File.exists?("/proc/version")
    @doc[:body][:environment][:software][:environment_variables] = `env`
    @doc[:body][:environment][:software][:cpu_governor] = `cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor`
    @doc[:body][:environment][:software][:running_softwares] = `ps -le`
    @doc[:body][:environment][:software][:users] = `who`
  end
  def self.register_kernel_info(key,infos)
    @doc[:body][:kernel_info][key] = infos
  end
  def self.dump_info
    File::open( @path, "w") { |f|
      f.print YAML::dump(@doc)
    }
  end
end
