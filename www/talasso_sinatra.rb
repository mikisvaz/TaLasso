require 'rubygems'
require 'sinatra'
require 'matrix_format'
require 'haml'
require 'sass'
require 'simplews'
require 'base64'
require 'yaml'
require 'paginator'
require 'soap/wsdlDriver'

RESULTS_DIR = File.join(File.dirname(File.expand_path(__FILE__)), 'public', 'results')
FileUtils.mkdir_p RESULTS_DIR unless File.exists? RESULTS_DIR

WSDL_FILE = File.join(File.dirname(File.dirname(File.expand_path(__FILE__))), 'webservice', 'wsdl', 'TaLassoWS.wsdl')
$driver = SOAP::WSDLDriverFactory.new(WSDL_FILE).create_rpc_driver

get '/favicon.ico' do
  ""
end

get '/aplication.css' do
  headers 'Content-Type' => 'text/css'
  sass :screen  
end

get '/wsdl' do
  send_file(WSDL_FILE, :filename => 'xmippWS.wsdl')
end

get '/example' do
  send_file(EXAMPLE_FILE, :filename => 'PolAB_msk4.spi', :type => 'application/xmipp')
end


get '/documentation' do
  haml :documentation
end

get '/' do
  @title ="Home"
  haml :index
end

post '/' do
  name  = params[:name] unless params[:name].nil?

  (gefile = params[:GE][:tempfile]) unless params[:GE].nil?
  (mefile = params[:ME][:tempfile]) unless params[:ME].nil?
  (genesfile = params[:genes][:tempfile]) unless params[:genes].nil?
  (mirnasfile = params[:mirnas][:tempfile]) unless params[:mirnas].nil?
  (samplesfile = params[:samples][:tempfile]) unless params[:samples].nil?

  data_id = $driver.upload(gefile.read, mefile.read, genesfile.read, mirnasfile.read, samplesfile.read)
  
  putatives = Array.new
  %w(microrna mirbase mirecords tarbase hoctar).each{|p|
    if params[p] == "1"
     putatives << p
    end
   }
   %w(DIANAmicroT  I_mirandaXL_pictar4way_targetscans  I_pictar4way_targetscans  Union  microrna  mirbase  pictar4way  pictar5way  targetscans).each{|p|
     if params["mirgen_#{p}"] == "1"
       putatives << "mirgen/#{p}"
     end
   }

  if name.nil? || name.empty?
    filename = params[:GE][:filename]
    name = File.basename(filename).sub(/\.[^\.]*$/,'')
  else
    name = name.gsub(/\s+/,"_")
  end

  puts data_id
  putatives.each{|p|
	puts p
   }


  job = case params[:algorithm]
    when "TaLasso" then $driver.lasso(data_id,putatives,name)
    when "Correlation" then $driver.correlation(data_id,putatives,name)
    when "GenMir" then $driver.gen_mir(data_id,putatives,name)
    when "SimpleGenMir" then $driver.simple_gen_mir(data_id,putatives,name)
  end

  # Change this information to match you actual web serice
  redirect "/" + job
end

get '/help' do
  haml :help
end


get '/:job' do
  @job   = params[:job]
  @title = @job

  case 
  when $driver.error(@job)
    @error = $driver.messages(@job).last
    @title += " [Error]"
    haml :error

  when ! $driver.done(@job)
    @status   = $driver.status(@job)
    @messages = $driver.messages(@job)
    @title += " [#{@status}]"
    haml :wait

  else
    @results = $driver.results(@job)
    File.open(File.join(RESULTS_DIR,"#{@job}_targets.txt"), 'w') do |f| f.write Base64.decode64 $driver.result(@results[0]) end unless File.exists? File.join(RESULTS_DIR,"#{@job}_targets.txt")
    File.open(File.join(RESULTS_DIR,"#{@job}_gene.txt"), 'w') do |f| f.write Base64.decode64 $driver.result(@results[1]) end    unless File.exists? File.join(RESULTS_DIR,"#{@job}_gene.txt")
    File.open(File.join(RESULTS_DIR,"#{@job}_mirna.txt"), 'w') do |f| f.write Base64.decode64 $driver.result(@results[2]) end   unless File.exists? File.join(RESULTS_DIR,"#{@job}_mirna.txt")

    @info = $driver.info(@job)
    @title += " [Done]"
    
	#def marshal_cache(name, key = [])

    targets = MatrixFormat::matrix2list( File.join(RESULTS_DIR,"#{@job}_targets.txt"), File.join(RESULTS_DIR,"#{@job}_gene.txt"), File.join(RESULTS_DIR,"#{@job}_mirna.txt"))
    sorted_targets = targets.collect{|gen,info| info.collect{|mirna,score| {:gen => gen, :mirna => mirna, :score => score.to_f}}}.flatten.sort_by{|p| p[:score]}.reverse
    @pager = Paginator.new(sorted_targets.size, 50) do |offset, per_page|
    	sorted_targets[offset, per_page]
    end
    @page = @pager.page(params[:page])
    haml :results
  end

end

get '/:job/m/:mirna' do
  @job   = params[:job]
  @title = @job
  @mirna = params[:mirna]

  targets = MatrixFormat::matrix2list( File.join(RESULTS_DIR,"#{@job}_targets.txt"), File.join(RESULTS_DIR,"#{@job}_gene.txt"), File.join(RESULTS_DIR,"#{@job}_mirna.txt"))
  sorted_targets = targets.collect{|gen,info| info.collect{|mirna,score| {:gen => gen, :mirna => mirna, :score => score.to_f}}}.flatten.sort_by{|p| p[:score]}.reverse
  @mirna_targets = sorted_targets.select{|p| p[:mirna] == @mirna}

  haml :mirna_results
end

get '/:job/g/:gen' do
  @job   = params[:job]
  @title = @job
  @gen = params[:gen]
  targets = MatrixFormat::matrix2list( File.join(RESULTS_DIR,"#{@job}_targets.txt"), File.join(RESULTS_DIR,"#{@job}_gene.txt"), File.join(RESULTS_DIR,"#{@job}_mirna.txt"))
  sorted_targets = targets.collect{|gen,info| info.collect{|mirna,score| {:gen => gen, :mirna => mirna, :score => score.to_f}}}.flatten.sort_by{|p| p[:score]}.reverse
  @gene_targets = sorted_targets.select{|p| p[:gen] == @gen}

  haml :gene_results
end
