module CacheHelper
  CACHE_DIR = File.join(File.dirname(File.expand_path(__FILE__)), '..', 'cache')
  FileUtils.mkdir_p(CACHE_DIR) unless File.exist?(CACHE_DIR)
 
  class CacheLocked < Exception; end
  def self.reset
    FileUtils.rm Dir.glob(CACHE_DIR + '*')
  end
 
  def self.reset_locks
    FileUtils.rm Dir.glob(CACHE_DIR + '*.lock')
  end
 
 
  def self.build_filename(name, key)
    File.join(CACHE_DIR, name + ": " + Digest::MD5.hexdigest(key.to_s))
  end
 
  def self.do(filename, block)
    FileUtils.touch(filename + '.lock')
    t = Time.now
    data = block.call
    STDERR.puts "#{ filename } time: #{Time.now - t}"
    File.open(filename, 'w'){|f| f.write data}
    FileUtils.rm(filename + '.lock')
    return data
  end
 
  def clean(name)
    FileUtils.rm Dir.glob(CACHE_DIR + "#{ name }*")
  end
 
  def cache_ready?(name, key)
    filename = CacheHelper.build_filename(name, key)
    File.exist?(filename)
  end
 
  def cache(name, key = [], wait = nil, &block)
    filename = CacheHelper.build_filename(name, key)
    begin
      case
      when File.exist?(filename)
        return File.open(filename){|f| f.read}
      when File.exist?(filename + '.lock')
        raise CacheLocked
      else
        if wait.nil?
          CacheHelper.do(filename, block)
        else
          Thread.new{CacheHelper.do(filename, block)}
          return wait
        end
 
      end
    rescue CacheLocked
      if wait.nil?
        sleep 30
        retry
      else
        return wait
      end
    rescue Exception
      FileUtils.rm(filename + '.lock') if File.exist?(filename + '.lock')
      raise $!
    end
  end
 
  def marshal_cache(name, key = [])
    Marshal::load( cache(name, key) do
      Marshal::dump(yield)
    end)
  end
end
